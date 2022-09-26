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
        ih = ImageHandler()
        totalNumberOfSubtomos = len(self.inputSetOfSubTomograms.get()) # *** esto no va a funcionar en cuanto haya coordenadas de mas de 1 tomo

        print("total number of subtomos: " + str(totalNumberOfSubtomos))

        self.loadModel()
        # self.generateVolumeIdList()

        subtomoList = []
        predictionList = np.array([0])
        # subtomoList = self.loadSubtomoSubset()

        currentVolId = None

        for index, subtomo in enumerate(self.inputSetOfSubTomograms.get().iterSubtomos(orderBy="_volId")):

            if (subtomo.getVolId() != currentVolId and currentVolId is not None) or ((index+1) == totalNumberOfSubtomos):
                # Make prediction
                print("Making prediction for subtomos with volId=" + str(currentVolId))
                print(predictionList.size)
                for s in subtomoList:
                    prediction = self.model.predict(s)
                    np.append(predictionList, prediction)

                overallPrediction = self.determineOverallPrediction(predictionList)

                print("for volume id " + currentVolId + "obtained prediction from " + subtomoList + " subtomos is "
                      + overallPrediction)

                print(predictionList)

                currentVolId = subtomo.getVolId()
                subtomoList = []

            if currentVolId is None:
                currentVolId = subtomo.getVolId()

            subtomoList.append(ih.read(subtomo.getFileName()).getData())

            print("subtomoList len " + str(len(subtomoList)))

    # --------------------------- UTILS functions ----------------------------
    def loadModel(self):
        self.model = load_model(self.model.get())
        print(self.model.summary())

    def loadSubtomoSubset(self):
        subtomoSubset = []

        return subtomoSubset

    def determineOverallPrediction(self, predictionList):
        predictionClasses = np.round(predictionList)

        overallPrediction = 0

        for i in predictionClasses:
            overallPrediction += i

        overallPrediction = overallPrediction / predictionList.size

        return True if overallPrediction > 0.5 else False

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        methods = []

        return methods
