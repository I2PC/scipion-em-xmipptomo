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


class XmippProtMisalignTiltSeries(EMProtocol, ProtTomoBase):
    """
    Introduce misalignment in the transformation matrix of a tilt-series
    """

    _label = 'misalign tilt-series '

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

        """ Options to introduce misalignment in the X axis shift"""
        form.addParam('shiftXNoiseToggle',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Introduce misalignment in shift X?',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Introduce noise in the shift alignment value in the X axis. Characterize the noise '
                           'behaviour through the parameters in the following formula:\n'
                           '\n'
                           'dx = a0 + a1 * i + a2 * sin((i + a3) / S * pi) + a4 * sin((i + a5) / S * 2 * pi) '
                           '+ N(0,a6)\n'
                           '\n'
                           'Being i the index position of the image inside the tilt-series, S the size of it and N a '
                           'normal distribution.'
                           'These parameters characterize the following behaviours:\n'
                           '- Constant (a0): an offset error (a0) is introduced in every image of the tilt-series.\n'
                           '- Incremental (a1): a constant incremental error (a1) is propagated through the '
                           'tilt-series.\n'
                           '- Sine lobe (a2, a3): the introduced error presents a half sine shape, characterized by '
                           'the error amplitude (a2) and the phase to displace the error function a given number of '
                           'images inside the tilt-series (a3).\n'
                           '- Sine cycle (a4, a5): the introduced error presents a full sine cycle shape, '
                           'characterized by the error amplitude (a4) and the phase to displace the error function a '
                           'given number of images inside the tilt-series (a5).\n'
                           '- Random (a6): a random error is introduced in every image of the tilt-series given a '
                           'sigma value (a6).\n')

        groupShiftX = form.addGroup('Misalignment parameters in shift X',
                                    condition='shiftXNoiseToggle==0')

        groupShiftX.addParam('a0param',
                             params.FloatParam,
                             default=0.0,
                             label='Offset error (a0)',
                             help='Offset shift error introduced in the X axis for every image of the tilt-series.')

        groupShiftX.addParam('a1param',
                             params.FloatParam,
                             default=0.0,
                             label='Incremental error (a1)',
                             help='Incremental shift error introduced in the X axis for every image of the '
                                  'tilt-series.')

        groupShiftX.addParam('a2param',
                             params.FloatParam,
                             default=0.0,
                             label='Sine lobe error amplitude (a2)',
                             help='Maximum amplitude of the sine lobe error function introduced in the X axis.')

        groupShiftX.addParam('a3param',
                             params.IntParam,
                             default=0,
                             label='Sine lobe error phase (a3)',
                             help='Phase (displacement) of the sine lobe error function. The introduced number '
                                  'corresponds to the number of images from the tilt-series that the origin of the '
                                  'error function will be displaced.')

        groupShiftX.addParam('a4param',
                             params.FloatParam,
                             default=0.0,
                             label='Sine error amplitude (a4)',
                             help='Maximum amplitude of the sine error function introduced in the X axis.')

        groupShiftX.addParam('a5param',
                             params.IntParam,
                             default=0,
                             label='Sine error phase (a5)',
                             help='Phase (displacement) of the sine error function. The introduced number corresponds '
                                  'to the number of images from the tilt-series that the origin of the error function '
                                  'will be displaced.')

        groupShiftX.addParam('a6param',
                             params.FloatParam,
                             default=0.0,
                             label='Random error sigma (a6)',
                             help='Sigma value for the random error introduced in the shift X.')

        """ Options to introduce misalignment in the Y axis shift"""
        form.addParam('shiftYNoiseToggle',
                      params.EnumParam,
                      choices=['Yes', 'No'],
                      default=1,
                      label='Introduce misalignment in shift Y?',
                      important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Introduce noise in the shift alignment value in the Y axis. Characterize the noise '
                           'behaviour through the parameters in the following formula:\n'
                           '\n'
                           'dY = b0 + b1 * i + b2 * sin((i + b3) / S * pi) + b4 * sin((i + b5) / S * 2 * pi) '
                           '+ N(0,b6)\n'
                           '\n'
                           'Being i the index position of the image inside the tilt-series, S the size of it and N a '
                           'normal distribution.'
                           'These parameters characterize the following behaviours:\n'
                           '- Constant (b0): an offset error (b0) is introduced in every image of the tilt-series.\n'
                           '- Incremental (b1): a constant incremental error (b1) is propagated through the '
                           'tilt-series.\n'
                           '- Sine lobe (b2, b3): the introduced error presents a half sine shape, characterized by '
                           'the error amplitude (b2) and the phase to displace the error function a given number of '
                           'images inside the tilt-series (b3).\n'
                           '- Sine cycle (b4, b5): the introduced error presents a full sine cycle shape, '
                           'characterized by the error amplitude (b4) and the phase to displace the error function a '
                           'given number of images inside the tilt-series (b5).\n'
                           '- Random (b6): a random error is introduced in every image of the tilt-series given a '
                           'sigma value (b6).\n')

        groupShiftY = form.addGroup('Misalignment parameters in shift Y',
                                    condition='shiftYNoiseToggle==0')

        groupShiftY.addParam('b0param',
                             params.FloatParam,
                             default=0.0,
                             label='Offset error (b0)',
                             help='Offset shift error introduced in the Y axis for every image of the tilt-series.')

        groupShiftY.addParam('b1param',
                             params.FloatParam,
                             default=0.0,
                             label='Incremental error (b1)',
                             help='Incremental shift error introduced in the Y axis for every image of the '
                                  'tilt-series.')

        groupShiftY.addParam('b2param',
                             params.FloatParam,
                             default=0.0,
                             label='Sine lobe error amplitude (b2)',
                             help='Maximum amplitude of the sine lobe error function introduced in the Y axis.')

        groupShiftY.addParam('b3param',
                             params.IntParam,
                             default=0,
                             label='Sine lobe error phase (b3)',
                             help='Phase (displacement) of the sine lobe error function. The introduced number '
                                  'corresponds to the number of images from the tilt-series that the origin of the '
                                  'error function will be displaced.')

        groupShiftY.addParam('b4param',
                             params.FloatParam,
                             default=0.0,
                             label='Sine error amplitude (b4)',
                             help='Maximum amplitude of the sine error function introduced in the Y axis.')

        groupShiftY.addParam('b5param',
                             params.IntParam,
                             default=0,
                             label='Sine error phase (b5)',
                             help='Phase (displacement) of the sine error function. The introduced number corresponds '
                                  'to the number of images from the tilt-series that the origin of the error function '
                                  'will be displaced.')

        groupShiftY.addParam('b6param',
                             params.FloatParam,
                             default=0.0,
                             label='Random error sigma (b6)',
                             help='Sigma value for the random error introduced in the shift Y.')

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
                      help='Constant angle error to add for every image of the tilt-series. Angles are measured in '
                           'radians.')

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

        form.addParam('angleSineErrorOffset',
                      params.FloatParam,
                      default=0.0,
                      label='Error offset',
                      condition='angleNoiseType==3 or angleNoiseType==4',
                      help='Offset of the error function.')

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

        if self.shiftXNoiseType.get() == 3:  # Sine lobe
            transformMatrix[0, 2] = transformMatrix[0, 2] + \
                                    self.getSineLobeNoiseValue(index,
                                                               size,
                                                               self.shiftXSineErrorAmplitude.get(),
                                                               self.shiftXSineErrorPhase.get(),
                                                               self.shiftXSineErrorOffset.get())

        if self.shiftXNoiseType.get() == 4:  # Sine
            transformMatrix[0, 2] = transformMatrix[0, 2] + \
                                    self.getSineNoiseValue(index,
                                                           size,
                                                           self.shiftXSineErrorAmplitude.get(),
                                                           self.shiftXSineErrorPhase.get(),
                                                           self.shiftXSineErrorOffset.get())

        if self.shiftXNoiseType.get() == 5:  # Random
            transformMatrix[0, 2] = np.random.normal(transformMatrix[0, 2], self.shiftXRandomErrorSigma.get())

        """Shift in Y axis modifications"""
        if self.shiftYNoiseType.get() == 1:  # Constant
            transformMatrix[1, 2] = transformMatrix[1, 2] + self.shiftYConstantError.get()

        if self.shiftYNoiseType.get() == 2:  # Incremental
            transformMatrix[1, 2] = transformMatrix[1, 2] + \
                                    self.getIncrementalNoiseValue(index,
                                                                  size,
                                                                  self.shiftYIncrementalErrorInitial.get(),
                                                                  self.shiftYIncrementalErrorFinal.get())

        if self.shiftYNoiseType.get() == 3:  # Sine lobe
            transformMatrix[1, 2] = transformMatrix[1, 2] + \
                                    self.getSineLobeNoiseValue(index,
                                                               size,
                                                               self.shiftYSineErrorAmplitude.get(),
                                                               self.shiftYSineErrorPhase.get(),
                                                               self.shiftYSineErrorOffset.get())

        if self.shiftYNoiseType.get() == 4:  # Sine
            transformMatrix[1, 2] = transformMatrix[1, 2] + \
                                    self.getSineNoiseValue(index,
                                                           size,
                                                           self.shiftYSineErrorAmplitude.get(),
                                                           self.shiftYSineErrorPhase.get(),
                                                           self.shiftYSineErrorOffset.get())

        if self.shiftYNoiseType.get() == 5:  # Random
            transformMatrix[1, 2] = np.random.normal(transformMatrix[1, 2], self.shiftYRandomErrorSigma.get())

        """Angle modifications"""
        angleModified = False
        oldAngle = np.arccos(transformMatrix[0, 0])

        if self.angleNoiseType.get() == 1:  # Constant
            newAngle = oldAngle + self.angleConstantError.get()
            angleModified = True

        if self.angleNoiseType.get() == 2:  # Incremental
            newAngle = oldAngle + self.getIncrementalNoiseValue(index,
                                                                size,
                                                                self.angleIncrementalErrorInitial.get(),
                                                                self.angleIncrementalErrorFinal.get())
            angleModified = True

        if self.angleNoiseType.get() == 3:  # Sine lobe
            newAngle = oldAngle + self.getSineLobeNoiseValue(index,
                                                             size,
                                                             self.angleSineErrorAmplitude.get(),
                                                             self.angleSineErrorPhase.get(),
                                                             self.angleSineErrorOffset.get())
            angleModified = True

        if self.angleNoiseType.get() == 4:  # Sine
            newAngle = oldAngle + self.getSineNoiseValue(index,
                                                         size,
                                                         self.angleSineErrorAmplitude.get(),
                                                         self.angleSineErrorPhase.get(),
                                                         self.angleSineErrorOffset.get())
            angleModified = True

        if self.angleNoiseType.get() == 5:  # Random
            newAngle = oldAngle + np.random.normal(oldAngle, self.angleRandomErrorSigma.get())
            angleModified = True

        if angleModified:
            transformMatrix[0, 0] = np.cos(newAngle)
            transformMatrix[0, 1] = - np.sin(newAngle)
            transformMatrix[1, 0] = np.sin(newAngle)
            transformMatrix[1, 1] = np.cos(newAngle)

        return transformMatrix

    @staticmethod
    def getIncrementalNoiseValue(index, size, low, high):
        slope = (high - low) / size

        increment = low + (slope * index)

        return increment

    @staticmethod
    def getSineLobeNoiseValue(index, size, amplitude, phase, offset):
        abscissa = (index + phase) / size
        if abscissa > 1:
            abscissa = abscissa - 1

        increment = offset + amplitude * np.sin(abscissa * np.pi)

        return increment

    @staticmethod
    def getSineNoiseValue(index, size, amplitude, phase, offset):
        abscissa = (index + phase) / size
        if abscissa > 1:
            abscissa = abscissa - 1

        increment = offset + amplitude * np.sin(abscissa * 2 * np.pi)

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
