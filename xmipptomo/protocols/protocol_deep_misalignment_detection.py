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

import numpy as np
from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pyworkflow.object import Set, Float
from pyworkflow.protocol import FileParam, PointerParam
from tomo.protocols import ProtTomoBase
from tensorflow.keras.models import load_model


class XmippProtDeepDetectMisalignment(EMProtocol, ProtTomoBase):
    """
    Wrapper protocol to Xmipp image peak high contrast applied to any volume
    """

    _label = 'detect misalignment from fiducials'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.alignedTomograms = None
        self.misalignedTomograms = None
        self.outputSetOfSubtomograms = None

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSetOfSubTomograms',
                      PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      label='Input set of subtomos',
                      help='Set of 3D coordinates of the location of the fiducials used to predict if the tomogram '
                           'presents misalignment.')

        form.addParam('firstModelPath', FileParam,
                      label='Input model for first split',
                      help='Input model for first split prediction')

        form.addParam('secondModelPath', FileParam,
                      label='Input model for second split',
                      help='Input model for second split prediction')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.subtomoPrediction)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEP functions --------------------------------
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
                for i, subtomo in enumerate(subtomoList):
                    subtomo._misaliScore = Float(predictionArray[i])

                    self.outputSetOfSubtomograms.append(subtomo)
                    self.outputSetOfSubtomograms.write()
                    self._store()

                currentVolId = subtomo.getVolId()
                subtomoList = []

            if currentVolId is None:
                currentVolId = subtomo.getVolId()

            subtomoList.append(subtomo)

    def closeOutputSetsStep(self):
        if self.alignedTomograms:
            self.alignedTomograms.setStreamState(Set.STREAM_CLOSED)
            self.alignedTomograms.write()

        if self.misalignedTomograms:
            self.misalignedTomograms.setStreamState(Set.STREAM_CLOSED)
            self.misalignedTomograms.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
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
              str(predictionList.size-overallPrediction) + "misaligned")
        
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
        if self.outputSetOfSubtomograms:
            self.outputSetOfSubtomograms.enableAppend()

        else:
            outputSetOfSubtomograms = self._createSetOfSubTomograms()

            outputSetOfSubtomograms.copyInfo(self.inputSetOfSubTomograms.get())

            outputSetOfSubtomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfSubtomograms=outputSetOfSubtomograms)
            self._defineSourceRelation(self.inputSetOfSubTomograms, outputSetOfSubtomograms)

        return self.outputSetOfSubtomograms

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        methods = []

        return methods
