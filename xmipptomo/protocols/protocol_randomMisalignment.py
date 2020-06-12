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


class XmippProtRandomMisalignment(EMProtocol, ProtTomoBase):
    """
    Introduce a random misalignment in the transformation matrix of a tilt-series
    """

    _label = 'random tilt-series misalignment'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series.')

        """ Options to introduce misalignment in the X axis shift"""
        form.addParam('shiftXNoiseType',
                      params.EnumParam,
                      choices=['None', 'Constant', 'Incremental', 'Sine lobe', 'Sine cycle', 'Random'],
                      default=0,
                      label='Shift misalignment type',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Introduce an specific type of noise in the shift alignment value in the X axis:\n'
                           '- None: no noise is introduced in any image of the tilt-series.\n'
                           '- Constant: the same error value is introduced in every image of the tilt-series.\n'
                           '- Incremental: an incremental error is introduced in each image given an initial and final '
                           'error for the first and last image of the tilt-series.\n'
                           '- Sine lobe: the introduced error presents a "half sine" shape, indicating the amplitude '
                           'of the maximum error and a phase to displace the error function a given number of image.\n'
                           '- Sine cycle: the error introduced presents a "full sine cycle", indicating the amplitude '
                           'of the maximum error and a phase to displace the error function a given number of images.\n'
                           '- Random: a random error is introduced in every image of the tilt-series given a sigma '
                           'value.')

        form.addParam('shiftXConstantError',
                      params.FloatParam,
                      default=0.0,
                      label='Error value',
                      condition='shiftXNoiseType==1',
                      help='Constant shift to add in the X axis for every image of the tilt-series.')

        form.addParam('shiftXIncrementalErrorInitial',
                      params.FloatParam,
                      default=0.0,
                      label='Initial error value',
                      condition='shiftXNoiseType==2',
                      help='Initial shift value in the X axis for the first image (lowest angle) of the tilt-series.')

        form.addParam('shiftXIncrementalErrorFinal',
                      params.FloatParam,
                      default=0.0,
                      label='Final error value',
                      condition='shiftXNoiseType==2',
                      help='Final shift value in the X axis for the last image (highest angle) of the tilt-series.')

        form.addParam('shiftXSineErrorAmplitude',
                      params.FloatParam,
                      default=0.0,
                      label='Error amplitude',
                      condition='shiftXNoiseType==3 or shiftXNoiseType==4',
                      help='Maximum shift value in the X axis for the error function.')

        form.addParam('shiftXSineErrorPhase',
                      params.IntParam,
                      default=0,
                      label='Error phase',
                      condition='shiftXNoiseType==3 or shiftXNoiseType==4',
                      help='Phase (displacement) of the error function. The number introduced corresponds to the '
                           'number of images from the tilt-series that the origin of the error function is going to be '
                           'displaced.')

        form.addParam('shiftXRandomErrorSigma',
                      params.FloatParam,
                      default=2.0,
                      label='Shift sigma',
                      condition='shiftXNoiseType==5',
                      help='Sigma value for random error introduced in the shift X.')

        """ Options to introduce misalignment in the Y axis shift"""
        form.addParam('shiftYNoiseType',
                      params.EnumParam,
                      choices=['None', 'Constant', 'Incremental', 'Sine lobe', 'Sine cycle', 'Random'],
                      default=0,
                      label='Shift misalignment type',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Introduce an specific type of noise in the shift alignment value in the Y axis:\n'
                           '- None: no noise is introduced in any image of the tilt-series.\n'
                           '- Constant: the same error value is introduced in every image of the tilt-series\n'
                           '- Incremental: an incremental error is introduced in each image given an initial and final '
                           'error for the first and last image of the tilt-series.\n'
                           '- Sine lobe: the introduced error presents a "half sine" shape, indicating the amplitude '
                           'of the maximum error and a phase to displace the error function a given number of image.\n'
                           '- Sine cycle: the error introduced presents a "full sine cycle", indicating the amplitude '
                           'of the maximum error and a phase to displace the error function a given number of image.\n'
                           '- Random: a random error is introduced in every image of the tilt-series given a sigma '
                           'value.')

        form.addParam('shiftYConstantError',
                      params.FloatParam,
                      default=0.0,
                      label='Error value',
                      condition='shiftYNoiseType==1',
                      help='Constant shift to add in the Y axis for every image of the tilt-series.')

        form.addParam('shiftYIncrementalErrorInitial',
                      params.FloatParam,
                      default=0.0,
                      label='Initial error value',
                      condition='shiftYNoiseType==2',
                      help='Initial shift value in the Y axis for the first image (lowest angle) of the tilt-series.')

        form.addParam('shiftYIncrementalErrorFinal',
                      params.FloatParam,
                      default=0.0,
                      label='Final error value',
                      condition='shiftYNoiseType==2',
                      help='Final shift value in the Y axis for the last image (highest angle) of the tilt-series.')

        form.addParam('shiftYSineErrorAmplitude',
                      params.FloatParam,
                      default=0.0,
                      label='Error amplitude',
                      condition='shiftYNoiseType==3 or shiftYNoiseType==4',
                      help='Maximum shift value in the Y axis for the error function.')

        form.addParam('shiftYSineErrorPhase',
                      params.IntParam,
                      default=0,
                      label='Error phase',
                      condition='shiftYNoiseType==3 or shiftYNoiseType==4',
                      help='Phase (displacement) of the error function. The number introduced corresponds to the '
                           'number of images from the tilt-series that the origin of the error function is going to be '
                           'displaced.')

        form.addParam('shiftYRandomErrorSigma',
                      params.FloatParam,
                      default=2.0,
                      label='Shift sigma',
                      condition='shiftYNoiseType==5',
                      help='Sigma value for random error introduced in the shift Y.')

        """ Options to introduce misalignment in the angle"""
        form.addParam('angleNoiseType',
                      params.EnumParam,
                      choices=['None', 'Constant', 'Incremental', 'Sine lobe', 'Sine cycle', 'Random'],
                      default=0,
                      label='Angle misalignment type',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Introduce an specific type of noise in the angle alignment value in the Y axis:\n'
                           '- None: no noise is introduced in any image of the tilt-series.\n'
                           '- Constant: the same error value is introduced in every image of the tilt-series\n'
                           '- Incremental: an incremental error is introduced in each image given an initial and final '
                           'error for the first and last image of the tilt-series.\n'
                           '- Sine lobe: the introduced error presents a "half sine" shape, indicating the amplitude '
                           'of the maximum error and a phase to displace the error function a given number of image.\n'
                           '- Sine cycle: the error introduced presents a "full sine cycle", indicating the amplitude '
                           'of the maximum error and a phase to displace the error function a given number of image.\n'
                           '- Random: a random error is introduced in every image of the tilt-series given a sigma '
                           'value.')

        form.addParam('angleConstantError',
                      params.FloatParam,
                      default=0.0,
                      label='Error value',
                      condition='angleNoiseType==1',
                      help='Constant angle to add in the Y axis for every image of the tilt-series. Angles are '
                           'measured in radians.')

        form.addParam('angleIncrementalErrorInitial',
                      params.FloatParam,
                      default=0.0,
                      label='Initial error value',
                      condition='angleNoiseType==2',
                      help='Initial angle error value for the first image (lowest angle) of the tilt-series. Angles '
                           'are measured in radians.')

        form.addParam('angleIncrementalErrorFinal',
                      params.FloatParam,
                      default=0.0,
                      label='Final error value',
                      condition='angleNoiseType==2',
                      help='Final angle error value for the last image (highest angle) of the tilt-series. Angles '
                           'are measured in radians.')

        form.addParam('angleSineErrorAmplitude',
                      params.FloatParam,
                      default=0.0,
                      label='Error amplitude',
                      condition='angleNoiseType==3 or angleNoiseType==4',
                      help='Maximum angle error value for the error function. Angles are measured in radians.')

        form.addParam('angleSineErrorPhase',
                      params.IntParam,
                      default=0,
                      label='Error phase',
                      condition='angleNoiseType==3 or angleNoiseType==4',
                      help='Phase (displacement) of the error function. The number introduced corresponds to the '
                           'number of images from the tilt-series that the origin of the error function is going to be '
                           'displaced.')

        form.addParam('angleRandomErrorSigma',
                      params.FloatParam,
                      default=2.0,
                      label='Angle sigma',
                      condition='angleNoiseType==5',
                      help='Sigma value for random error introduced in the angle. Angles are measured in radians')

        """ Options for misalignment interpolation"""
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
            newTransformMat = self.modifyTransformMatrix(transformMat, index, ts.getSize())

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

    # --------------------------- UTILS functions ----------------------------
    def modifyTransformMatrix(self, transformMatrix, index, size):
        """Shift in X axis modifications"""
        if self.shiftXNoiseType.get() == 1:  # Constant
            transformMatrix[0, 2] = transformMatrix[0, 2] + self.shiftXConstantError.get()
        if self.shiftXNoiseType.get() == 2:  # Incremental
            transformMatrix[0, 2] = transformMatrix[0, 2] + \
                                    self.getIncrementalNoiseValue(index,
                                                                  size,
                                                                  self.shiftXIncrementalErrorInitial.get(),
                                                                  self.shiftXIncrementalErrorFinal.get())
        if self.shiftXNoiseType.get() == 5:  # Random
            transformMatrix[0, 2] = np.random.normal(transformMatrix[0, 2], self.shiftXRandomErrorSigma.get())

        """Shift in Y axis modifications"""
        if self.shiftYNoiseType.get() == 1:  # Constant
            transformMatrix[0, 2] = transformMatrix[0, 2] + self.shiftYConstantError.get()
        if self.shiftYNoiseType.get() == 2:  # Incremental
            transformMatrix[0, 2] = transformMatrix[0, 2] + \
                                    self.getIncrementalNoiseValue(index,
                                                                  size,
                                                                  self.shiftYIncrementalErrorInitial.get(),
                                                                  self.shiftYIncrementalErrorFinal.get())
        if self.shiftYNoiseType.get() == 5:  # Random
            transformMatrix[1, 2] = np.random.normal(transformMatrix[1, 2], self.shiftYRandomErrorSigma.get())

        """Angle modifications"""
        angleModified = False
        if not self.angleMisalignment.get() == 0:
            oldAngle = np.arccos(transformMatrix[0, 0])
        if self.angleNoiseType.get() == 1:  # Constant
            newAngle = oldAngle + self.angleConstantError.get()
            angleModified = True
        if self.angleNoiseType.get() == 2:  # Incremental
            transformMatrix[0, 2] = transformMatrix[0, 2] + \
                                    self.getIncrementalNoiseValue(index,
                                                                  size,
                                                                  self.angleIncrementalErrorInitial.get(),
                                                                  self.angleIncrementalErrorFinal.get())
            angleModified = True
        if self.angleMisalignment.get() == 5:  # Random
            newAngle = oldAngle + np.random.normal(oldAngle, self.angleRandomErrorSigma.get())
            angleModified = True

        if angleModified:
            transformMatrix[0, 0] = np.cos(newAngle)
            transformMatrix[0, 1] = - np.sin(newAngle)
            transformMatrix[1, 0] = np.sin(newAngle)
            transformMatrix[1, 1] = np.cos(newAngle)

        return transformMatrix

    def getIncrementalNoiseValue(self, index, size, low, high):
        slope = (high-low)/size
        increment = low + (slope * index)
        return increment


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
