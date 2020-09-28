# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
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
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from pwem.emlib.image import ImageHandler


class XmippProtApplyTransformationMatrixTS(EMProtocol, ProtTomoBase):
    """
    Compute the interpolated tilt-series from its transform matrix.
    """

    _label = 'apply alignment tilt-series'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('applyAlignmentStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def applyAlignmentStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        path.makePath(extraPrefix)
        outputTsFileName = os.path.join(extraPrefix, "%s.mrc" % tsId)

        ts.applyTransform(outputTsFileName)

        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputInterpolatedSetOfTiltSeries.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1, outputTsFileName)
            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))
        newTs.write()

        outputInterpolatedSetOfTiltSeries.update(newTs)
        outputInterpolatedSetOfTiltSeries.updateDim()
        outputInterpolatedSetOfTiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputInterpolatedSetOfTiltSeries(self):
        if not hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)
        return self.outputInterpolatedSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        for ts in self.inputSetOfTiltSeries.get():
            if not ts.getFirstItem().hasTransform():
                validateMsgs.append("Some tilt-series from the input set of tilt-series is missing from a "
                                    "transformation matrix.")

        return validateMsgs

    def _summary(self):
        summary = []
        if hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nInterpolations applied: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The interpolation has been computed for %d "
                           "Tilt-series using the IMOD newstack program.\n"
                           % (self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
