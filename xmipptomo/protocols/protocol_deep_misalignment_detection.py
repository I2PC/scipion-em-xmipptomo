# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os

import numpy as np

from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pyworkflow.object import Set, Float
from pyworkflow.protocol import FileParam, PointerParam, EnumParam
from tomo.objects import SetOfCoordinates3D, SetOfTomograms
from tomo.protocols import ProtTomoBase
from tensorflow.keras.models import load_model
from xmipptomo import utils

COORDINATES_FILE_NAME = '3dCoordinates.xmd'
BOX_SIZE = 32


class XmippProtDeepDetectMisalignment(EMProtocol, ProtTomoBase):
    """
    Wrapper protocol to Xmipp image peak high contrast applied to any volume
    """

    _label = 'detect misalignment from fiducials'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.alignedTomograms = None
        self.misalignedTomograms = None
        self.outputSubtomos = None

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSetOfCoordinates',
                      PointerParam,
                      pointerClass=SetOfCoordinates3D,
                      label='Fiducial 3D coordinates',
                      help='3D coordinates in dicating the location of the fiducials (gold beads) in the tomogram. '
                           'These fiducails will be the ones used to study misalignment artifacts over them. '
                           'The coordinate denotes the center of the subtomogram')

        form.addParam('tomoSource', EnumParam,
                      choices=['same as picking', 'other'],
                      default=0,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Tomogram source',
                      help='By default the subtomograms will be extracted from the tomogram used in the picking step '
                           '( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide a different tomogram to extract from. \n'
                           '*Note*: In the _other_ case, ensure that provided tomogram and coordinates are related ')

        form.addParam('inputSetOfTomograms',
                      PointerParam,
                      pointerClass=SetOfTomograms,
                      condition='tomoSource==1',
                      allowsNull=True,
                      label='Input tomograms',
                      help='Tomograms from which extract the fiducials (gold beads) at the specified coordinates '
                           'locations.')

        # form.addParam('inputSetOfSubTomograms',
        #               PointerParam,
        #               pointerClass='SetOfSubTomograms',
        #               important=True,
        #               label='Input set of subtomos',
        #               help='Set of 3D coordinates of the location of the fiducials used to predict if the tomogram '
        #                    'presents misalignment.')

        form.addParam('firstModelPath', FileParam,
                      label='Input model for first split',
                      help='Input model for first split prediction')

        form.addParam('secondModelPath', FileParam,
                      label='Input model for second split',
                      help='Input model for second split prediction')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self.tomoDict = self.getTomoDict()

        for key in self.tomoDict.keys():
            tomo = self.tomoDict[key]
            coordFilePath = self.getExtraPath(os.path.join(tomo.getTsId()), COORDINATES_FILE_NAME)

            self._insertFunctionStep(utils.writeMdCoordinates,
                                     self.inputSetOfCoordinates.get(),
                                     tomo,
                                     coordFilePath)
            self._insertFunctionStep(self.extractSubtomos,
                                     key,
                                     coordFilePath)
            self._insertFunctionStep(self.subtomoPrediction)
            self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEP functions --------------------------------
    def extractSubtomos(self, key, coordFilePath):
        tomo = self.tomoDict[key]

        outputPath = self.getExtraPath(os.path.join(tomo.getTsId()))
        tomoFn = tomo.getFileName()

        paramsExtractSubtomos = {
            'tomogram': tomoFn,
            'coordinates': coordFilePath,
            'boxsize': BOX_SIZE,
            'threads': 1,  # ***
            'output': outputPath
        }

        argsExtractSubtomos = "--tomogram %(tomogram)s " \
                              "--coordinates %(coordinates)s " \
                              "--boxsize %(boxsize)d " \
                              "--threads %(threads)d " \
                              "--o %(outputPath)s " \
                              "--subtomo "

        self.runJob('xmipp_tomo_extract_subtomograms', argsExtractSubtomos % paramsExtractSubtomos)

    def subtomoPrediction(self):
        totalNumberOfSubtomos = len(self.inputSetOfSubTomograms.get())

        print("total number of subtomos: " + str(totalNumberOfSubtomos))

        self.loadModels()

        subtomoList = []
        currentVolId = None

        self.getOutputSetOfSubtomos()

        for index, subtomo in enumerate(self.inputSetOfSubTomograms.get().iterSubtomos(orderBy="_volId")):

            if (subtomo.getVolId() != currentVolId and currentVolId is not None) or (
                    (index + 1) == totalNumberOfSubtomos):
                # Make prediction
                overallPrediction, predictionArray = self.makePrediction(subtomoList)

                print("For volume id " + str(currentVolId) + " obtained prediction from " +
                      str(len(subtomoList)) + " subtomos is " + str(overallPrediction))

                self.addTomoToOutput(volId=currentVolId, overallPrediction=overallPrediction)

                # Generate output set of subtomograms with a prediction score
                for i, s in enumerate(subtomoList):
                    s._misaliScore = Float(predictionArray[i])

                    self.outputSubtomos.append(s)
                    self.outputSubtomos.write()
                    self._store()

                currentVolId = subtomo.getVolId()
                subtomoList = []

            if currentVolId is None:
                currentVolId = subtomo.getVolId()

            subtomoList.append(subtomo.clone())

    def closeOutputSetsStep(self):
        if self.alignedTomograms:
            self.alignedTomograms.setStreamState(Set.STREAM_CLOSED)
            self.alignedTomograms.write()

        if self.misalignedTomograms:
            self.misalignedTomograms.setStreamState(Set.STREAM_CLOSED)
            self.misalignedTomograms.write()

        self.outputSubtomos.setStreamState(Set.STREAM_CLOSED)
        self.outputSubtomos.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getTomoDict(self):
        if self.tomoSource.get() == 0:
            self.inputSetOfCoordinates.get().getPrecedentsInvolved()
        else:
            tomoDict = {}
            for tomo in self.inputSetOfTomograms.get():
                tsId = tomo.getTsId()
                tomoDict[tsId] = tomo

        return tomoDict

    def makePrediction(self, subtomoList):
        ih = ImageHandler()

        numberOfSubtomos = len(subtomoList)

        subtomoArray = np.zeros((numberOfSubtomos, 32, 32, 32), dtype=np.float64)

        for index, subtomo in enumerate(subtomoList):
            subtomoDataTmp = ih.read(subtomo.getFileName()).getData()

            subtomoArray[index, :, :, :] = subtomoDataTmp[:, :, :]

        std = subtomoArray.std()
        mean = subtomoArray.mean()

        subtomoArray = (subtomoArray - mean) / std

        predictionArray = self.firstModel.predict(subtomoArray)

        overallPrediction = self.determineOverallPrediction(predictionArray)

        if overallPrediction:
            predictionArray = self.secondModel.predict(subtomoArray)

            overallPrediction = self.determineOverallPrediction(predictionArray)

        return overallPrediction, predictionArray

    def addTomoToOutput(self, volId, overallPrediction):
        if overallPrediction:  # Ali
            self.getOutputSetOfAlignedTomograms()
            newTomogram = self.inputSetOfSubTomograms.get().getTomograms()[volId].clone()

            self.alignedTomograms.append(newTomogram)
            self.alignedTomograms.write()
            self._store()

        else:  # Misali
            self.getOutputSetOfMisalignedTomograms()
            newTomogram = self.inputSetOfSubTomograms.get().getTomograms()[volId].clone()

            self.misalignedTomograms.append(newTomogram)
            self.misalignedTomograms.write()
            self._store()

    def loadModels(self):
        self.firstModel = load_model(self.firstModelPath.get())
        print(self.firstModel.summary())

        self.secondModel = load_model(self.secondModelPath.get())
        print(self.firstModel.summary())

    @staticmethod
    def determineOverallPrediction(predictionList):
        predictionClasses = np.round(predictionList)

        overallPrediction = 0

        for i in predictionClasses:
            overallPrediction += i

        print("Subtomo analysis: " + str(overallPrediction) + " aligned vs " +
              str(predictionList.size - overallPrediction) + "misaligned")

        overallPrediction = overallPrediction / predictionList.size

        return True if overallPrediction > 0.5 else False  # aligned (1) or misaligned (0)

    def getOutputSetOfAlignedTomograms(self):
        if self.alignedTomograms:
            self.alignedTomograms.enableAppend()

        else:
            alignedTomograms = self._createSetOfTomograms(suffix='Ali')

            alignedTomograms.copyInfo(self.inputSetOfSubTomograms.get())

            alignedTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(alignedTomograms=alignedTomograms)
            self._defineSourceRelation(self.inputSetOfSubTomograms, alignedTomograms)

        return self.alignedTomograms

    def getOutputSetOfMisalignedTomograms(self):
        if self.misalignedTomograms:
            self.misalignedTomograms.enableAppend()

        else:
            misalignedTomograms = self._createSetOfTomograms(suffix='Misali')

            misalignedTomograms.copyInfo(self.inputSetOfSubTomograms.get())

            misalignedTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(misalignedTomograms=misalignedTomograms)
            self._defineSourceRelation(self.inputSetOfSubTomograms, misalignedTomograms)

        return self.misalignedTomograms

    def getOutputSetOfSubtomos(self):
        if self.outputSubtomos:
            self.outputSubtomos.enableAppend()

        else:
            outputSubtomos = self._createSetOfSubTomograms()

            outputSubtomos.copyInfo(self.inputSetOfSubTomograms.get())

            outputSubtomos.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSubtomos=outputSubtomos)
            self._defineSourceRelation(self.inputSetOfSubTomograms, outputSubtomos)

        return self.outputSubtomos

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        methods = []

        return methods
