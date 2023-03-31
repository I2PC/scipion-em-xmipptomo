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
from pyworkflow import BETA
from pyworkflow.object import Set, Float
from pyworkflow.protocol import FileParam, PointerParam, EnumParam
from tomo import constants
from tomo.objects import SetOfCoordinates3D, SetOfTomograms, SubTomogram, Coordinate3D
from tomo.protocols import ProtTomoBase
from tensorflow.keras.models import load_model
from xmipptomo import utils

COORDINATES_FILE_NAME = 'subtomo_coords.xmd'
BOX_SIZE = 32
TARGET_SAMPLING_RATE = 6.25


class XmippProtDeepDetectMisalignment(EMProtocol, ProtTomoBase):
    """
    Wrapper protocol to Xmipp image peak high contrast applied to any volume
    """

    _label = 'detect misalignment from fiducials'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.alignedTomograms = None
        self.misalignedTomograms = None
        self.outputSubtomos = None
        self.isot = None

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
            coordFilePath = self._getExtraPath(os.path.join(tomo.getTsId()), COORDINATES_FILE_NAME)

            self._insertFunctionStep(utils.writeMdCoordinates,
                                     self.inputSetOfCoordinates.get(),
                                     tomo,
                                     coordFilePath)
            self._insertFunctionStep(self.extractSubtomos,
                                     key,
                                     coordFilePath)
            self._insertFunctionStep(self.subtomoPrediction,
                                     key,
                                     coordFilePath)
            self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEP functions --------------------------------
    def extractSubtomos(self, key, coordFilePath):
        tomo = self.tomoDict[key]

        outputPath = self._getExtraPath(os.path.join(tomo.getTsId()))
        tomoFn = tomo.getFileName()

        paramsExtractSubtomos = {
            'tomogram': tomoFn,
            'coordinates': coordFilePath,
            'boxsize': BOX_SIZE,
            'threads': 1,  # ***
            'outputPath': outputPath,
            'downsample': TARGET_SAMPLING_RATE / tomo.getSamplingRate(),
        }

        argsExtractSubtomos = "--tomogram %(tomogram)s " \
                              "--coordinates %(coordinates)s " \
                              "--boxsize %(boxsize)d " \
                              "--threads %(threads)d " \
                              "-o %(outputPath)s " \
                              "--downsample %(downsample)f" \
                              "--subtomo "

        self.runJob('xmipp_tomo_extract_subtomograms', argsExtractSubtomos % paramsExtractSubtomos)

    def subtomoPrediction(self, key, coordFilePath):
        tomo = self.tomoDict[key]
        subtomoPathList = self.getSubtomoPathList(coordFilePath)

        subtomoCoordList = utils.retrieveXmipp3dCoordinatesIntoList(coordFilePath, xmdFormat=1)

        totalNumberOfSubtomos = len(subtomoPathList)
        tsId = tomo.getTsId()

        print("Analyzing tomogram " + tsId)
        print("Total number of subtomos: " + str(totalNumberOfSubtomos))

        self.loadModels()

        overallPrediction, predictionArray = self.makePrediction(subtomoPathList)

        print("For volume id " + str(tsId) + " obtained prediction from " +
              str(len(subtomoPathList)) + " subtomos is " + str(overallPrediction))

        self.addTomoToOutput(tomo=tomo, overallPrediction=overallPrediction)

        # Generate output set of subtomograms with a prediction score
        self.getOutputSetOfSubtomos()

        for i, subtomoPath in enumerate(subtomoPathList):
            newCoord3D = Coordinate3D()
            newCoord3D.setVolume(tomo)
            newCoord3D.setVolId(i)
            newCoord3D.setX(subtomoCoordList[i][0], constants.BOTTOM_LEFT_CORNER)
            newCoord3D.setY(subtomoCoordList[i][1], constants.BOTTOM_LEFT_CORNER)
            newCoord3D.setZ(subtomoCoordList[i][2], constants.BOTTOM_LEFT_CORNER)

            subtomogram = SubTomogram()
            subtomogram.setLocation(subtomoPath)
            subtomogram.setCoordinate3D(newCoord3D)
            subtomogram.setSamplingRate(TARGET_SAMPLING_RATE)
            subtomogram.setVolName(tomo.getTsId())
            subtomogram._misaliScore = Float(predictionArray[i])

            self.outputSubtomos.append(subtomogram)
            self.outputSubtomos.write()
            self._store()

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
            self.isot = self.inputSetOfCoordinates.get().getPrecedents()
            tomoDict = self.inputSetOfCoordinates.get().getPrecedentsInvolved()
        else:
            tomoDict = {}
            self.isot = self.inputSetOfTomograms.get()

            for tomo in self.inputSetOfTomograms.get():
                tsId = tomo.getTsId()
                tomoDict[tsId] = tomo

        return tomoDict

    def makePrediction(self, subtomoPathList):
        ih = ImageHandler()

        numberOfSubtomos = len(subtomoPathList)

        subtomoArray = np.zeros((numberOfSubtomos, 32, 32, 32), dtype=np.float64)

        for index, subtomo in enumerate(subtomoPathList):
            subtomoDataTmp = ih.read(subtomo)
            subtomoDataTmp = subtomoDataTmp.getData()

            print("subtomo " + str(index) + " mean " + str(subtomoDataTmp.mean()))
            print("subtomo " + str(index) + " std " + str(subtomoDataTmp.std()))

            subtomoArray[index, :, :, :] = subtomoDataTmp[:, :, :]

        std = subtomoArray.std()
        mean = subtomoArray.mean()

        subtomoArray = (subtomoArray - mean) / std

        predictionArray = self.firstModel.predict(subtomoArray)

        overallPrediction = self.determineOverallPrediction(predictionArray)

        print("first prediction array")
        print(predictionArray)
        print("first overall prediction " + str(overallPrediction))

        if overallPrediction:
            predictionArray = self.secondModel.predict(subtomoArray)

            overallPrediction = self.determineOverallPrediction(predictionArray)

            print("second prediction array")
            print(predictionArray)
            print("second overall prediction " + str(overallPrediction))

        return overallPrediction, predictionArray

    def addTomoToOutput(self, tomo, overallPrediction):
        if overallPrediction:  # Ali
            self.getOutputSetOfAlignedTomograms()
            self.alignedTomograms.append(tomo)
            self.alignedTomograms.write()
            self._store()

        else:  # Misali
            self.getOutputSetOfMisalignedTomograms()
            self.misalignedTomograms.append(tomo)
            self.misalignedTomograms.write()
            self._store()

    def loadModels(self):
        self.firstModel = load_model(self.firstModelPath.get())
        print(self.firstModel.summary())

        self.secondModel = load_model(self.secondModelPath.get())
        print(self.firstModel.summary())

    @staticmethod
    def determineOverallPrediction(predictionList):
        # predictionClasses = np.round(predictionList)
        #
        # overallPrediction = 0
        #
        # for i in predictionClasses:
        #     overallPrediction += i
        #
        # print("Subtomo analysis: " + str(overallPrediction) + " aligned vs " +
        #       str(predictionList.size - overallPrediction) + "misaligned")
        #
        # overallPrediction = overallPrediction / predictionList.size
        #
        # return True if overallPrediction > 0.5 else False  # aligned (1) or misaligned (0)

        predictionClasses = predictionList

        overallPrediction = 0

        for i in predictionClasses:
            overallPrediction += i

        overallPrediction = overallPrediction / predictionList.size

        print("Subtomo analysis preditcion: " + overallPrediction)

        return True if overallPrediction > 0.5 else False  # aligned (1) or misaligned (0)
    def getOutputSetOfAlignedTomograms(self):
        if self.alignedTomograms:
            self.alignedTomograms.enableAppend()

        else:
            alignedTomograms = self._createSetOfTomograms(suffix='Ali')

            alignedTomograms.copyInfo(self.isot)

            alignedTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(alignedTomograms=alignedTomograms)
            self._defineSourceRelation(self.isot, alignedTomograms)

        return self.alignedTomograms

    def getOutputSetOfMisalignedTomograms(self):
        if self.misalignedTomograms:
            self.misalignedTomograms.enableAppend()

        else:
            misalignedTomograms = self._createSetOfTomograms(suffix='Misali')

            misalignedTomograms.copyInfo(self.isot)

            misalignedTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(misalignedTomograms=misalignedTomograms)
            self._defineSourceRelation(self.isot, misalignedTomograms)

        return self.misalignedTomograms

    def getOutputSetOfSubtomos(self):
        if self.outputSubtomos:
            self.outputSubtomos.enableAppend()

        else:
            outputSubtomos = self._createSetOfSubTomograms()

            outputSubtomos.copyInfo(self.isot)
            outputSubtomos.setSamplingRate(TARGET_SAMPLING_RATE)

            outputSubtomos.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSubtomos=outputSubtomos)
            self._defineSourceRelation(self.isot, outputSubtomos)

        return self.outputSubtomos

    @staticmethod
    def getSubtomoPathList(coordFilePath):
        coordFilePath_noExt = os.path.splitext(coordFilePath)[0]
        counter = 1

        subtomoPathList = []

        while True:
            subtomoPath = coordFilePath_noExt + '-' + str(counter) + '.mrc'

            if not os.path.exists(subtomoPath):
                break

            subtomoPathList.append(subtomoPath)
            counter += 1

        return subtomoPathList

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        methods = []

        return methods
