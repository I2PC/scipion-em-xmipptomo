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
import pyworkflow.utils.path as path
import pyworkflow.utils as pwutils
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
import tomo.objects as tomoObj



class XmippProtSplitTiltSeries(EMProtocol, ProtTomoBase):
    """
    Wrapper protocol to Xmipp split Odd Even on tilt-series
    """
    _label = 'split tilt-series'

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series',
                      help='Select a set of tilt-series to be split into two sets (odd and even).'
                           'It means, the set of tilt-series is split in two subsets.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('splitTiltSeries', ts.getObjId())
            self._insertFunctionStep('convertXmdToStackStep', ts.getObjId())
            self._insertFunctionStep('createOutputStep', ts.getObjId())

    # --------------------------- STEPS functions -------------------------------
    def splitTiltSeries(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        tsFileName = ts.getFirstItem().getFileName()
        path.makePath(self._getExtraPath(tsId))

        tsFileNameOdd = pwutils.removeExt(os.path.basename(tsFileName)) + "_odd.xmd"
        tsFileNameEven = pwutils.removeExt(os.path.basename(tsFileName)) + "_even.xmd"

        paramsOddEven = {
            'inputImg': tsFileName,
            'outputOdd': self._getExtraPath(os.path.join(tsId, tsFileNameOdd)),
            'outputEven': self._getExtraPath(os.path.join(tsId, tsFileNameEven)),
            'type': "frames",
        }

        argsOddEven = "--img %(inputImg)s " \
                      "-o %(outputOdd)s " \
                      "-e %(outputEven)s " \
                      "--type %(type)s "

        self.runJob('xmipp_image_odd_even', argsOddEven % paramsOddEven)

    def convertXmdToStackStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        tsFileName = ts.getFirstItem().getFileName()

        tsFileNameOdd = pwutils.removeExt(os.path.basename(tsFileName)) + "_odd.xmd"
        tsFileNameEven = pwutils.removeExt(os.path.basename(tsFileName)) + "_even.xmd"

        tsFileNameOddMrc = pwutils.removeExt(os.path.basename(tsFileNameOdd)) + ".mrc"
        tsFileNameEvenMrc = pwutils.removeExt(os.path.basename(tsFileNameEven)) + ".mrc"

        paramsConvertOdd = {
            'inputXmdOdd': self._getExtraPath(os.path.join(tsId, tsFileNameOdd)),
            'outputMrcOdd': self._getExtraPath(os.path.join(tsId, tsFileNameOddMrc)),
        }

        argsConvertOdd = "-i %(inputXmdOdd)s " \
                         "-o %(outputMrcOdd)s "

        self.runJob('xmipp_image_convert', argsConvertOdd % paramsConvertOdd)

        paramsConvertEven = {
            'inputXmdEven': self._getExtraPath(os.path.join(tsId, tsFileNameEven)),
            'outputMrcEven': self._getExtraPath(os.path.join(tsId, tsFileNameEvenMrc)),
        }

        argsConvertEven = "-i %(inputXmdEven)s " \
                          "-o %(outputMrcEven)s "

        self.runJob('xmipp_image_convert', argsConvertEven % paramsConvertEven)

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()
        tsFileName = ts.getFirstItem().getFileName()

        """Output even set"""
        outputOddSetOfTiltSeries = self.getOutputEvenSetOfTiltSeries()
        tsFileNameEvenMrc = pwutils.removeExt(os.path.basename(tsFileName)) + "_even.mrc"

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputOddSetOfTiltSeries.append(newTs)

        dimCounter = 0
        for index, tiltImage in enumerate(ts):
            if index % 2 == 0:
                dimCounter += 1
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(dimCounter, self._getExtraPath(os.path.join(tsId, tsFileNameEvenMrc)))
                newTs.append(newTi)
        newTs.write()
        outputOddSetOfTiltSeries.update(newTs)
        outputOddSetOfTiltSeries.write()
        self._store()

        """Output odd set"""
        outputOddSetOfTiltSeries = self.getOutputOddSetOfTiltSeries()
        tsFileNameOddMrc = pwutils.removeExt(os.path.basename(tsFileName)) + "_odd.mrc"

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)
        outputOddSetOfTiltSeries.append(newTs)

        dimCounter = 0
        for index, tiltImage in enumerate(ts):
            if index % 2 == 1:
                dimCounter += 1
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setLocation(dimCounter, self._getExtraPath(os.path.join(tsId, tsFileNameOddMrc)))
                newTs.append(newTi)
        newTs.write()
        outputOddSetOfTiltSeries.update(newTs)
        outputOddSetOfTiltSeries.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputEvenSetOfTiltSeries(self):
        if not hasattr(self, "outputEvenSetOfTiltSeries"):
            outputEvenSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Even')
            outputEvenSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputEvenSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputEvenSetOfTiltSeries=outputEvenSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputEvenSetOfTiltSeries)
        return self.outputEvenSetOfTiltSeries

    def getOutputOddSetOfTiltSeries(self):
        if not hasattr(self, "outputOddSetOfTiltSeries"):
            outputOddSetOfTiltSeries = self._createSetOfTiltSeries(suffix='Odd')
            outputOddSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputOddSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            self._defineOutputs(outputOddSetOfTiltSeries=outputOddSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputOddSetOfTiltSeries)
        return self.outputOddSetOfTiltSeries

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        messages = []
        if hasattr(self, 'resolution_Volume'):
            messages.append(
                'Information about the method/article in ')
        return messages

    def _summary(self):
        summary = []

        return summary

    def _citations(self):
        return ['']
