# **************************************************************************
# *
# * Authors:       Jose Luis Vilas Prieto (jlvilas@cnb.csic.es) [1]
# *                Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
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
# **************************************************************************ç

import os

import numpy as np

from pwem.emlib.image import ImageHandler
from pwem.objects import Micrograph
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, StringParam, IntParam
import pyworkflow.utils.path as path
from pyworkflow.object import String, Float
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
import xmipptomo.utils as utils


class XmippProtAverageViewTiltSeries(EMProtocol, ProtTomoBase):
    """
    Scipion protocol to average a subset of tilt-images for gaining SNR to posteriorly use a SPA 2D picker.
    """

    _label = 'Average tilt-series views'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtTomoBase.__init__(self)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('avgAngleList',
                      StringParam,
                      important=True,
                      label='Angles of average',
                      help='List of angles (split by commas) indicating the angles at which to perform the average.'
                           'For example: -25, 0, 25')

        form.addParam('numberViewsAverage',
                      IntParam,
                      important=True,
                      label='Number of views to average',
                      help='Number of views to average centered in each angle given in the "Angles of average" list')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):

        for ts in self.inputSetOfTiltSeries.get():
            tsObjId = ts.getObjId()
            self._insertFunctionStep(self.convertInputStep,
                                     tsObjId)

            self._insertFunctionStep(self.averageViews,
                                     tsObjId)

            self._insertFunctionStep(self.createOutputStep,
                                     tsObjId)

    # --------------------------- STEPS functions ----------------------------

    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        firstItem = ts.getFirstItem()

        """Apply the transformation form the input tilt-series"""
        # Use Xmipp interpolation via Scipion

        if firstItem.hasTransform():
            avgRotAngle = utils.calculateRotationAngleFromTM(ts)
            swap = True if (avgRotAngle > 45 or avgRotAngle < -45) else False

            outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
            ts.applyTransform(outputTsFileName, swapXY=swap)

        else:
            outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
            ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))
        utils.writeXmippMetadataTiltAngleList(ts, angleFilePath)

    def averageViews(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()

        tiltAngleList = self.getTiltAngleList(ts)

        tmpTiltImage = os.path.join(tmpPrefix, firstItem.parseFileName(suffix="_tmp", extension=".mrc"))
        avgAngleList = self.avgAngleList.get().split(',')

        ih = ImageHandler()

        ih.createEmptyImage(fnOut=tmpTiltImage,
                            xDim=firstItem.getXDim(),
                            yDim=firstItem.getYDim(),
                            nDim=1)

        for a, avgAngle in enumerate(avgAngleList):
            difference = 999

            for i, angle in enumerate(tiltAngleList):
                if abs(float(avgAngle) - angle) < difference:
                    difference = abs(float(avgAngle) - angle)
                    index = i

            outputFilePath = os.path.join(extraPrefix, firstItem.parseFileName(suffix="_" + avgAngle.strip(), extension=".mrc"))

            ih.createEmptyImage(fnOut=outputFilePath,
                                xDim=firstItem.getXDim(),
                                yDim=firstItem.getYDim(),
                                nDim=1)

            for i in range(index - int(self.numberViewsAverage.get() / 2),
                           index + int(self.numberViewsAverage.get() / 2) + 1):

                t = np.array([[np.cos(np.radians(tiltAngleList[index]) - np.radians(tiltAngleList[i])), 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]])

                ih.applyTransform(inputFile=str(i) + "@" + os.path.join(tmpPrefix, firstItem.parseFileName()),
                                  outputFile=str(1) + "@" + tmpTiltImage,
                                  transformMatrix=t.flatten(),
                                  shape=(firstItem.getYDim(), firstItem.getXDim()))

                paramsImageOperate = {
                    'i1': str(1) + "@" + tmpTiltImage,
                    'i2': str(1) + "@" + outputFilePath,
                    'out': str(1) + "@" + outputFilePath,
                }

                argsImageOperate = "-i %(i1)s " \
                                   "--plus %(i2)s " \
                                   "-o %(out)s "

                self.runJob('xmipp_image_operate', argsImageOperate % paramsImageOperate)

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        avgAngleList = self.avgAngleList.get().split(',')

        setOfMicrographs = self._createSetOfMicrographs(suffix='_ts_average')

        for a, angle in enumerate(avgAngleList):
            tsAvg = Micrograph()

            outputFilePath = os.path.join(extraPrefix, firstItem.parseFileName(suffix="_"+angle.strip(), extension=".mrc"))

            tsAvg.setFileName(outputFilePath)
            tsAvg.setSamplingRate(firstItem.getSamplingRate())
            tsAvg._tsId = String(tsId)
            tsAvg._avgAngle = Float(angle)

            setOfMicrographs.append(tsAvg)

        setOfMicrographs.copyInfo(ts)
        setOfMicrographs.setSamplingRate(ts.getSamplingRate())

        self._defineOutputs(tsAverage=setOfMicrographs)
        self._defineSourceRelation(self.inputSetOfTiltSeries, setOfMicrographs)

    # --------------------------- UTILS functions ----------------------------
    @staticmethod
    def getTiltAngleList(ts):
        angleList = []

        for ti in ts:
            angleList.append(ti.getTiltAngle())

        return angleList

    # --------------------------- INFO functions ----------------------------
