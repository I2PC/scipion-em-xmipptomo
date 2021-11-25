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

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol.params import EnumParam, FloatParam, IntParam, BooleanParam, PointerParam
import pyworkflow.protocol.constants as const
from tomo.protocols import ProtTomoBase
import os
import pyworkflow.utils.path as path
from pyworkflow.object import Set
import tomo.objects as tomoObj
from pwem.emlib.image import ImageHandler
from xmipptomo.protocols.protocol_resize_base import XmippProtResizeBase


class XmippProtResizeTiltSeries(XmippProtResizeBase):
    """
    Wrapper protocol to Xmipp image resize applied on tilt-series
    """

    _label = 'resize tilt-series'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSetOfTiltSeries',
                      PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series',
                      help='Select a set of tilt-series to be resized.')
        self._defineParamsReSize(form)
    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):

        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('resizeTiltSeries', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())
        self._insertFunctionStep('closeStreamStep')

    # --------------------------- STEP functions --------------------------------
    def resizeTiltSeries(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        print(ts)
        tsId = ts.getTsId()
        print(tsId)
        extraPrefix = self._getExtraPath(tsId)
        print(extraPrefix)
        path.makePath(extraPrefix)

        samplingRate = self.inputSetOfTiltSeries.get().getSamplingRate()

        args = self.resizeCommonArgsResize(self, samplingRate)

        for ti in ts:
            self.runJob("xmipp_image_resize", self._ioArgs(ti)+args)

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]

        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries()

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputSetOfTiltSeries.append(newTs)

        newTs.setSamplingRate(self.samplingRate)

        fileName = ts.getFirstItem().parseFileName()

        for index, ti in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(ti, copyId=True)
            newTi.setLocation(index + 1, os.path.join(extraPrefix, fileName))

            if ti.hasTransform():
                newTi.setTransform(ti.getTransform())

            newTi.setSamplingRate(self.samplingRate)

            newTs.append(newTi)

        ih = ImageHandler()
        x, y, z, _ = ih.getDimensions(newTs.getFirstItem().getFileName())
        newTs.setDim((x, y, z))
        newTs.write(properties=False)

        outputSetOfTiltSeries.update(newTs)
        outputSetOfTiltSeries.updateDim()
        outputSetOfTiltSeries.write()

        self._store()

    def closeStreamStep(self):
        self.getOutputSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions -------------------------------
    def getTsSize(self):
        """ Get the X dimension of the tilt-series from the set """
        return self.inputSetOfTiltSeries.get().getDim()[0]

    def _ioArgs(self, ti):
        tsId = ti.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        tiIndex = ti.getLocation()[0]
        tsPath = os.path.join(extraPrefix, ti.parseFileName())

        inputTs = str(tiIndex) + ":mrcs@" + ti.getFileName()
        outputTs = str(tiIndex) + "@" + tsPath

        return "-i %s -o %s " % (inputTs, outputTs)


    def getOutputSetOfTiltSeries(self):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputSetOfTiltSeries.setSamplingRate(self.samplingRate)
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\n"
                           "Interpolations applied: %d.\n"
                           "Target sampling rate: %.2f A/px.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputSetOfTiltSeries.getSize(),
                              self.outputSetOfTiltSeries.getSamplingRate()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfTiltSeries'):
            methods.append("%d tilt-series have been interpolated using the xmipp_image_resize program to a %.2f A/px "
                           "target sampling rate.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputSetOfTiltSeries.getSamplingRate()))
        else:
            methods.append("Output classes not ready yet.")
        return methods

