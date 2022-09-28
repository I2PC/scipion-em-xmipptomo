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
from pyworkflow.object import Set
from pyworkflow.protocol import FileParam, PointerParam
from tomo.protocols import ProtTomoBase
from tensorflow.keras.models import load_model


class XmippProtDeepDetectMisalignment(EMProtocol, ProtTomoBase):
    """
    Wrapper protocol to Xmipp image peak high contrast applied to any volume
    """

    _label = 'detect misalignment from fiducials'

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

        form.addParam('model', FileParam,
                      label='Input model',
                      help='Input model for prediction')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.subtomoPrediction)

    # --------------------------- STEP functions --------------------------------
    def subtomoPrediction(self):
        totalNumberOfSubtomos = len(self.inputSetOfSubTomograms.get())

        print("total number of subtomos: " + str(totalNumberOfSubtomos))

        self.loadModel()

        subtomoFilePathList = []
        currentVolId = None

        for index, subtomo in enumerate(self.inputSetOfSubTomograms.get().iterSubtomos(orderBy="_volId")):

            if (subtomo.getVolId() != currentVolId and currentVolId is not None) or (
                    (index + 1) == totalNumberOfSubtomos):
                # Make prediction
                overallPrediction = self.makePrediction(subtomoFilePathList)

                print("For volume id " + str(currentVolId) + " obtained prediction from " +
                      str(len(subtomoFilePathList)) + " subtomos is " + str(overallPrediction))

                self.addTomoToOutput(volId=currentVolId, overallPrediction=overallPrediction)

                currentVolId = subtomo.getVolId()
                subtomoFilePathList = []

            if currentVolId is None:
                currentVolId = subtomo.getVolId()

            subtomoFilePathList.append(subtomo.getFileName())

    # --------------------------- UTILS functions ----------------------------
    def makePrediction(self, subtomoFilePathList):
        ih = ImageHandler()

        numberOfSubtomos = len(subtomoFilePathList)

        subtomoArray = np.zeros((numberOfSubtomos, 32, 32, 32), dtype=np.float64)

        for index, subtomoFilePath in enumerate(subtomoFilePathList):
            subtomoDataTmp = ih.read(subtomoFilePath).getData()

            subtomoArray[index, :, :, :] = subtomoDataTmp[:, :, :]

        predictionArray = self.model.predict(subtomoArray)

        overallPrediction = self.determineOverallPrediction(predictionArray)

        return overallPrediction

    def addTomoToOutput(self, volId, overallPrediction):
        if overallPrediction:  # Ali
            self.getOutputSetOfAlignedTomograms()
            newTomogram = self.inputSetOfSubTomograms.get().getTomograms()[volId]

            self.outputSetOfAlignedTomograms.append(newTomogram)
            self.outputSetOfAlignedTomograms.update(newTomogram)
            self.outputSetOfAlignedTomograms.write()
            self._store()

        else:  # Misali
            self.getOutputSetOfMisalignedTomograms()
            newTomogram = self.inputSetOfSubTomograms.get().getTomograms()[volId]

            self.outputSetOfMisalignedTomograms.append(newTomogram)
            self.outputSetOfMisalignedTomograms.update(newTomogram)
            self.outputSetOfMisalignedTomograms.write()
            self._store()

    def loadModel(self):
        self.model = load_model(self.model.get())
        print(self.model.summary())

    @staticmethod
    def determineOverallPrediction(predictionList):
        predictionClasses = np.round(predictionList)

        overallPrediction = 0

        for i in predictionClasses:
            overallPrediction += i

        overallPrediction = overallPrediction / predictionList.size

        return True if overallPrediction > 0.5 else False  # aligned (1) or misaligned (0)

    def getOutputSetOfAlignedTomograms(self):
        outputSetOfAlignedTomograms = self._createSetOfTomograms()

        outputSetOfAlignedTomograms.setAcquisition(self.inputSetOfSubTomograms.get().getAcquisition())
        outputSetOfAlignedTomograms.setSamplingRate(self.inputSetOfSubTomograms.get().getSamplingRate())

        outputSetOfAlignedTomograms.setStreamState(Set.STREAM_OPEN)

        self._defineOutputs(outputSetOfAlignedTomograms=outputSetOfAlignedTomograms)
        self._defineSourceRelation(self.inputSetOfSubTomograms, outputSetOfAlignedTomograms)

        return self.outputSetOfAlignedTomograms

    def getOutputSetOfMisalignedTomograms(self):
        outputSetOfMisalignedTomograms = self._createSetOfTomograms()

        outputSetOfMisalignedTomograms.setAcquisition(self.inputSetOfSubTomograms.get().getAcquisition())
        outputSetOfMisalignedTomograms.setSamplingRate(self.inputSetOfSubTomograms.get().getSamplingRate())

        outputSetOfMisalignedTomograms.setStreamState(Set.STREAM_OPEN)

        self._defineOutputs(outputSetOfMisalignedTomograms=outputSetOfMisalignedTomograms)
        self._defineSourceRelation(self.inputSetOfSubTomograms, outputSetOfMisalignedTomograms)

        return self.outputSetOfMisalignedTomograms

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        methods = []

        return methods
