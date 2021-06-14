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
from xmipptomo import utils
from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase


class XmippProtRotateAstigmatism(EMProtocol, ProtTomoBase):
    """
    Rotate the astigmatism of a set of ctf tilt-series estimation given a set of transformation matrices comming from
    a set of tilt series.
    """

    _label = 'Astigmatism rotation'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('getTMSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      help='Set of tilt-series from which transformation matrices will be obtained.',
                      label='Set of tilt-series from which get transform')

        form.addParam('inputSetOfCtfTomoSeries',
                      params.PointerParam,
                      label="input tilt-series CTF estimation",
                      pointerClass='SetOfCTFTomoSeries',
                      help='Select the CTF estimation from the set of tilt-series.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.getTMSetOfTiltSeries.get():
            self._insertFunctionStep(self.rotateAstimatism, ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def rotateAstimatism(self, tsObjId):
        outputAssignedTransformSetOfTiltSeries = self.getOutputAssignedTransformSetOfTiltSeries()

        getTMTS = self.getTMSetOfTiltSeries.get()[tsObjId]
        inputCtfTomoSeries = self.inputSetOfCtfTomoSeries.get()[tsObjId]
        tsId = getTMTS.getTsId()

        self.getOutputSetOfCTFTomoSeries()

        newCTFTomoSeries = tomoObj.CTFTomoSeries()
        newCTFTomoSeries.copyInfo(inputCtfTomoSeries)
        newCTFTomoSeries.setTiltSeries(getTMTS)
        newCTFTomoSeries.setTsId(tsId)
        newCTFTomoSeries.setObjId(tsObjId)
        newCTFTomoSeries.setIMODDefocusFileFlag(inputCtfTomoSeries.getIMODDefocusFileFlag())
        newCTFTomoSeries.setNumberOfEstimationsInRange(inputCtfTomoSeries.getNumberOfEstimationsInRange())

        self.outputSetOfCTFTomoSeries.append(newCTFTomoSeries)

        for index, tiltImageGetTM, inputCtfTomo in enumerate(zip(getTMTS, inputCtfTomoSeries)):
            newCTFTomo = tomoObj.CTFTomo()
            newCTFTomo.copyInfo(inputCtfTomo)

            rotationAngle = utils.calculateRotationAngleFromTM(tiltImageGetTM)

            print(rotationAngle)

        newCTFTomoSeries.append(newCTFTomo)

        newCTFTomoSeries.write(properties=False)

        self.outputSetOfCTFTomoSeries.update(newCTFTomoSeries)
        self.outputSetOfCTFTomoSeries.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfCTFTomoSeries(self):
        if hasattr(self, "outputSetOfCTFTomoSeries"):
            self.outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = tomoObj.SetOfCTFTomoSeries.create(self._getPath(),
                                                                         template='CTFmodels%s.sqlite')
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(self.getTMSetOfTiltSeries.get())
            outputSetOfCTFTomoSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfCTFTomoSeries=outputSetOfCTFTomoSeries)
        return self.outputSetOfCTFTomoSeries

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        for tsGetTM, tsSetTM in zip(self.getTMSetOfTiltSeries.get(), self.setTMSetOfTiltSeries.get()):
            if not tsGetTM.getFirstItem().hasTransform():
                validateMsgs.append("Some tilt-series from the input set of tilt-series does not have a "
                                    "transformation matrix assigned.")

            if tsGetTM.getSize() != tsSetTM.getSize():
                validateMsgs.append("Some tilt-series from the input set of tilt-series and its target in the assign "
                                    "transformation set of tilt-series size's do not match. Every input tilt-series "
                                    "and its target must have the same number of elements")

        if self.getTMSetOfTiltSeries.get().getSize() != self.setTMSetOfTiltSeries.get().getSize():
            validateMsgs.append("Both input sets of tilt-series size's do not match. Both sets must have the same "
                                "number of elements.")

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputAssignedTransformSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices assigned: %d.\n"
                           % (self.getTMSetOfTiltSeries.get().getSize(),
                              self.outputAssignedTransformSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputAssignedTransformSetOfTiltSeries'):
            methods.append("The transformation matrix has been assigned to %d Tilt-series from the input set.\n"
                           % (self.outputAssignedTransformSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
