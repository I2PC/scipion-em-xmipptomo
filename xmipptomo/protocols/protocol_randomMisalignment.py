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
import pwem.objects as data
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
from imod import Plugin
from pwem.emlib.image import ImageHandler


class XmippProtRandomMisalignment(EMProtocol, ProtTomoBase):
    """
    Introduce a random misalignment in the transformation matrix of a tilt-series
    """

    _label = 'random misalignment'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('shiftMisalignment',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Shift misalignment',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Apply misalignment in the shifting terms of the transformation matrix.')

        form.addParam('shiftSigma',
                      params.FloatParam,
                      default=0.1,
                      label='Shift sigma',
                      condition='shiftMisalignment==0',
                      help='Sigma value for the shift misalignment')

        form.addParam('angleMisalignment', params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Angle misalignment',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Apply misalignment in the angle terms of the transformation matrix.')

        form.addParam('angleSigma',
                      params.FloatParam,
                      default=0.1,
                      label='Angle sigma',
                      condition='angleMisalignment==0',
                      help='Sigma value for the angle misalignment (radians).')

        form.addParam('computeAlignment', params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Generate interpolated tilt-series',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Generate and save the interpolated tilt-series applying the'
                           'obtained transformation matrix.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('introduceRandomMisalignment', ts.getObjId())
            # if self.computeAlignment.get() == 0:
            #     self._insertFunctionStep('computeInterpolation', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def introduceRandomMisalignment(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        outputMisalignedSetOfTiltSeries = self.getOutputMisalignedSetOfTiltSeries()
        missAliTs = tomoObj.TiltSeries(tsId=tsId)
        missAliTs.copyInfo(ts)
        outputMisalignedSetOfTiltSeries.append(missAliTs)

        for index, ti in enumerate(ts):
            missAliTi = tomoObj.TiltImage()
            missAliTi.copyInfo(ti, copyId=True)
            missAliTi.setLocation(ti.getLocation())

            if ti.hasTransform():
                transformMat = ti.getTransform().getMatrix()
            else:
                transformMat = np.identity(3)
            newTransformMat = self.modifyTransformMatrix(transformMat)

            newTransform = data.Transform()
            newTransform.setMatrix(newTransformMat)
            missAliTi.setTransform(newTransform)

            missAliTs.append(missAliTi)

        missAliTs.write()

        outputMisalignedSetOfTiltSeries.update(missAliTs)
        outputMisalignedSetOfTiltSeries.write()

        self._store()

        if self.computeAlignment.get() == 0:
            outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

            extraPrefix = self._getExtraPath(tsId)
            path.makePath(extraPrefix)
            outputTsFileName = os.path.join(extraPrefix, "%s_missAli.st" % tsId)

            """Apply the transformation form the input tilt-series"""
            missAliTs.applyTransform(outputTsFileName)

            missAliInterTs = tomoObj.TiltSeries(tsId=tsId)
            missAliInterTs.copyInfo(ts)
            outputInterpolatedSetOfTiltSeries.append(missAliInterTs)

            for index, tiltImage in enumerate(ts):
                missAliInterTi = tomoObj.TiltImage()
                missAliInterTi.copyInfo(tiltImage, copyId=True)
                missAliInterTi.setLocation(index + 1, outputTsFileName)
                missAliInterTs.append(missAliInterTi)

            missAliInterTs.write()

            outputInterpolatedSetOfTiltSeries.update(missAliInterTs)
            outputInterpolatedSetOfTiltSeries.write()

            self._store()

    def computeInterpolation(self, tsObjId):
        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        path.makePath(extraPrefix)
        outputTsFileName = os.path.join(extraPrefix, "%s_missAli.st" % tsId)

        """Apply the transformation form the input tilt-series"""
        ts.applyTransform(outputTsFileName)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputInterpolatedSetOfTiltSeries.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1, (os.path.join(extraPrefix, '%s_missAli.st' % tsId)))
            newTs.append(newTi)

        newTs.write()

        outputInterpolatedSetOfTiltSeries.update(newTs)
        outputInterpolatedSetOfTiltSeries.updateDim()
        outputInterpolatedSetOfTiltSeries.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def modifyTransformMatrix(self, transformMatrix):
        if self.shiftMisalignment.get() == 0:
            transformMatrix[0, 2] = np.random.normal(transformMatrix[0, 2], self.shiftSigma.get())
            transformMatrix[1, 2] = np.random.normal(transformMatrix[1, 2], self.shiftSigma.get())

        if self.angleMisalignment.get() == 0:
            angle = np.arccos(transformMatrix[0, 0])
            newAngle = angle + np.random.normal(angle, self.angleSigma.get())
            transformMatrix[0, 0] = np.cos(newAngle)
            transformMatrix[0, 1] = - np.sin(newAngle)
            transformMatrix[1, 0] = np.sin(newAngle)
            transformMatrix[1, 1] = np.cos(newAngle)

        return transformMatrix

    def getOutputMisalignedSetOfTiltSeries(self):
        if not hasattr(self, "outputMisalignedSetOfTiltSeries"):
            outputMisalignedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Misaligned')
            outputMisalignedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputMisalignedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputMisalignedSetOfTiltSeries=outputMisalignedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputMisalignedSetOfTiltSeries)
        return self.outputMisalignedSetOfTiltSeries

    def getOutputInterpolatedSetOfTiltSeries(self):
        if not hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)
        return self.outputInterpolatedSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputMisalignedSetOfTiltSeries.getSize()))

        elif hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           "Interpolated Tilt-Series: %d.\n"
                           % (self.outputMisalignedSetOfTiltSeries.getSize(),
                              self.outputMisalignedSetOfTiltSeries.getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("New transformation matrices has been calculated for %d Tilt-series.\n"
                           % (self.outputMisalignedSetOfTiltSeries.getSize()))

        elif hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("New transformation matrices has been calculated for %d Tilt-series.\n"
                           "Also, interpolation has been completed for %d Tilt-series.\n"
                           % (self.outputMisalignedSetOfTiltSeries.getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
