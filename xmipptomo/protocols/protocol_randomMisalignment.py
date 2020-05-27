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
import imod.utils as utils
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
                      default=1.0,
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
                      default=1.0,
                      label='Angle sigma',
                      condition='angleMisalignment==0',
                      help='Sigma value for the angle misalignment')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('convertInputStep', ts.getObjId())
            self._insertFunctionStep('computeXcorrStep', ts.getObjId())
            if self.computeAlignment.get() == 0:
                self._insertFunctionStep('computeInterpolatedStackStep', ts.getObjId())

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)
        outputTsFileName = os.path.join(tmpPrefix, "%s.st" % tsId)

        """Apply the transformation form the input tilt-series"""
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, "%s.rawtlt" % tsId)
        ts.generateTltFile(angleFilePath)

    def computeXcorrStep(self, tsObjId):
        """Compute transformation matrix for each tilt series"""
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsXcorr = {
            'input': os.path.join(tmpPrefix, '%s.st' % tsId),
            'output': os.path.join(extraPrefix, '%s.prexf' % tsId),
            'tiltfile': os.path.join(tmpPrefix, '%s.rawtlt' % tsId),
            'RotationAngle': self.rotationAngle.get(),
            'FilterSigma1': 0.03,
            'FilterSigma2': 0.05,
            'FilterRadius2': 0.25
        }
        argsXcorr = "-input %(input)s " \
                    "-output %(output)s " \
                    "-tiltfile %(tiltfile)s " \
                    "-RotationAngle %(RotationAngle)f " \
                    "-FilterSigma1 %(FilterSigma1)f " \
                    "-FilterSigma2 %(FilterSigma2)f " \
                    "-FilterRadius2 %(FilterRadius2)f"
        Plugin.runImod(self, 'tiltxcorr', argsXcorr % paramsXcorr)

        paramsXftoxg = {
            'input': os.path.join(extraPrefix, '%s.prexf' % tsId),
            'goutput': os.path.join(extraPrefix, '%s.prexg' % tsId),
        }
        argsXftoxg = "-input %(input)s " \
                     "-goutput %(goutput)s"
        Plugin.runImod(self, 'xftoxg', argsXftoxg % paramsXftoxg)

        """Generate output tilt series"""
        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries()
        tsId = ts.getTsId()
        alignmentMatrix = utils.formatTransformationMatrix(self._getExtraPath('%s/%s.prexg' % (tsId, tsId)))
        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputSetOfTiltSeries.append(newTs)
        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(tiltImage.getLocation())
            transform = data.Transform()
            transform.setMatrix(alignmentMatrix[:, :, index])
            newTi.setTransform(transform)
            newTs.append(newTi)
        newTs.write()
        outputSetOfTiltSeries.update(newTs)
        outputSetOfTiltSeries.write()
        self._store()

    def computeInterpolatedStackStep(self, tsObjId):
        outputInterpolatedSetOfTiltSeries = self.getOutputInterpolatedSetOfTiltSeries()

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        paramsAlignment = {
            'input': os.path.join(tmpPrefix, '%s.st' % tsId),
            'output': os.path.join(extraPrefix, '%s_preali.st' % tsId),
            'xform': os.path.join(extraPrefix, "%s.prexg" % tsId),
            'bin': int(self.binning.get()),
            'imagebinned': 1.0
        }
        argsAlignment = "-input %(input)s " \
                        "-output %(output)s " \
                        "-xform %(xform)s " \
                        "-bin %(bin)d " \
                        "-imagebinned %(imagebinned)s"
        Plugin.runImod(self, 'newstack', argsAlignment % paramsAlignment)

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputInterpolatedSetOfTiltSeries.append(newTs)

        if self.binning > 1:
            newTs.setSamplingRate(ts.getSamplingRate() * int(self.binning.get()))

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True)
            newTi.setLocation(index + 1, (os.path.join(extraPrefix, '%s_preali.st' % tsId)))
            if self.binning > 1:
                newTi.setSamplingRate(tiltImage.getSamplingRate() * int(self.binning.get()))
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
    def getOutputSetOfTiltSeries(self):
        if not hasattr(self, "outputSetOfTiltSeries"):
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries

    def getOutputInterpolatedSetOfTiltSeries(self):
        if not hasattr(self, "outputInterpolatedSetOfTiltSeries"):
            outputInterpolatedSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Interpolated')
            outputInterpolatedSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputInterpolatedSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            if self.binning > 1:
                samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()
                samplingRate *= self.binning.get()
                outputInterpolatedSetOfTiltSeries.setSamplingRate(samplingRate)
            self._defineOutputs(outputInterpolatedSetOfTiltSeries=outputInterpolatedSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputInterpolatedSetOfTiltSeries)
        return self.outputInterpolatedSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        if not self.inputSetOfTiltSeries.get().getFirstItem().getFirstItem().hasTransform():
            validateMsgs = "The input set of tilt-series must have an associated transform matrix"

        return validateMsgs

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputSetOfTiltSeries.getSize()))
        elif hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           "Interpolated Tilt-Series: %d.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputSetOfTiltSeries.getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The transformation matrix has been calculated for %d "
                           "Tilt-series using the IMOD procedure.\n"
                           % (self.outputSetOfTiltSeries.getSize()))
        elif hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The transformation matrix has been calculated for %d "
                           "Tilt-series using the IMOD procedure.\n"
                           "Also, interpolation has been completed for %d Tilt-series.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
