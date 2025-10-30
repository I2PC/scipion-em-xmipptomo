# **************************************************************************
# *
# * Authors:       Jose Luis Vilas Prieto (jlvilas@cnb.csic.es) [1]
# *                Federico P. de Isidro-Gomez (fp.deisidro@cnb.csi.es) [1]
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
import shutil

import numpy as np

from pwem.emlib.image import ImageHandler
from pwem.objects import Micrograph
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, FloatParam, IntParam, BooleanParam
import pyworkflow.utils.path as path
from pyworkflow.object import String, Float, Integer
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
import xmipptomo.utils as utils
from tomo.objects import TiltImage


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

        line = form.addLine('Angle range',
                            important=True,
                            help="Angle range over which the tilt series images will be averaged.")

        line.addParam('minAngle',
                      FloatParam,
                      default=-60,
                      label='Min')

        line.addParam('maxAngle',
                      FloatParam,
                      default=60,
                      label='Max')

        form.addParam('numberViewsAverage',
                      IntParam,
                      important=True,
                      label='Number of images to average',
                      help='Number of tilt-images to be averaged for each calculated mean. The averaging is '
                           'symmetrical so same number of images left and right to the center image wil be used for '
                           'averaging.')

        form.addParam('gaussFilter',
                      BooleanParam,
                      important=True,
                      default=False,
                      label="Apply Gaussian filter",
                      help='Filter calculated averages using a Gaussian kernel.')

        form.addParam('gaussStd',
                      FloatParam,
                      important=True,
                      condition='gaussFilter',
                      default=3,
                      label="Standard deviation",
                      help='Specify the standard deviation of the Gaussian filter.')

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
        outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())

        if firstItem.hasTransform():
            avgRotAngle = utils.calculateRotationAngleFromTM(ts)
            swap = True if (avgRotAngle > 45 or avgRotAngle < -45) else False

            ts.applyTransform(outputTsFileName, swapXY=swap)

        else:
            ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))
        utils.writeXmippMetadataTiltAngleList(ts, angleFilePath)

    def averageViews(self, tsObjId):
        ih = ImageHandler()

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()
        interpolatedTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
        cosStretchTiltImage = os.path.join(tmpPrefix, firstItem.parseFileName(suffix="_CS_tmp", extension=".mrc"))
        sliceStretchTiltImage = os.path.join(tmpPrefix, firstItem.parseFileName(suffix="_extractSlice", extension=".mrc"))

        ih.createEmptyImage(fnOut=cosStretchTiltImage,
                            xDim=firstItem.getXDim(),
                            yDim=firstItem.getYDim(),
                            zDim=1,
                            nDim=1)

        tiltAngleList = self.getTiltAngleList(ts)
        avgIndexList = [i for i, x in enumerate(tiltAngleList) if self.minAngle.get() <= x <= self.maxAngle.get()]
        sideImagesForAvg = int(float(self.numberViewsAverage.get()) / 2)
        maxIdx = len(tiltAngleList)

        for index in avgIndexList:
            print("----------- Processing image " + str(index) + " at angle " + str(tiltAngleList[index]))

            outputFilePathTmp = os.path.join(tmpPrefix,
                                             firstItem.parseFileName(suffix="_" + str(index + 1), extension=".mrc"))
            outputFilePathExtra = os.path.join(extraPrefix,
                                               firstItem.parseFileName(suffix="_" + str(index + 1), extension=".mrc"))

            ih.createEmptyImage(fnOut=outputFilePathTmp,
                                xDim=firstItem.getXDim(),
                                yDim=firstItem.getYDim(),
                                nDim=1)

            for i in range(index - sideImagesForAvg, index + sideImagesForAvg + 1):

                if i < 0 or i >= maxIdx:
                    continue

                centralAngle = tiltAngleList[index]
                projectedAngle = tiltAngleList[i]

                cosineStretchingFactor = np.cos(np.radians(centralAngle)) / np.cos(np.radians(projectedAngle))

                t = np.array([[cosineStretchingFactor, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]])

                # Extract image
                paramsImageOperateSlice = {
                    'i1': interpolatedTsFileName + ":mrc",
                    'slice': i,
                    'out': sliceStretchTiltImage,
                }

                argsImageOperateSlice = "-i %(i1)s " \
                                        "--slice %(slice)d " \
                                        "-o %(out)s "

                self.runJob('xmipp_image_operate', argsImageOperateSlice % paramsImageOperateSlice)

                # Apply cosine stretching
                print("\033[92m applyTransform(inputFile=" + sliceStretchTiltImage +
                      ", outputFile=" + cosStretchTiltImage + ", transformMatrix=transformMatrix, shape=(" +
                      str(firstItem.getYDim()) + "," + str(firstItem.getXDim()) + ")\033[0m")

                ih.applyTransform(inputFile=sliceStretchTiltImage,
                                  outputFile=cosStretchTiltImage,
                                  transformMatrix=t.flatten(),
                                  shape=(firstItem.getYDim(), firstItem.getXDim()))

                # Add to average
                paramsImageOperate = {
                    'i1': str(1) + "@" + cosStretchTiltImage,
                    'i2': str(1) + "@" + outputFilePathTmp,
                    'out': str(1) + "@" + outputFilePathTmp,
                }

                argsImageOperate = "-i %(i1)s " \
                                   "--plus %(i2)s " \
                                   "-o %(out)s "

                self.runJob('xmipp_image_operate', argsImageOperate % paramsImageOperate)

            if self.gaussFilter.get():
                paramsTransformFilter = {
                    'i': outputFilePathTmp,
                    'out': outputFilePathExtra,
                    'std': self.gaussStd.get(),
                }

                argsTransformFilter = "-i %(i)s " \
                                      "-o %(out)s " \
                                      "--fourier real_gaussian %(std)d"

                self.runJob('xmipp_transform_filter', argsTransformFilter % paramsTransformFilter)

            else:
                shutil.move(outputFilePathTmp, outputFilePathExtra)

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        acqOrderList = self.getAcqOrderList(ts)
        tiltAngleList = self.getTiltAngleList(ts)
        avgIndexList = [i for i, x in enumerate(tiltAngleList) if self.minAngle.get() <= x <= self.maxAngle.get()]

        setOfMicrographs = self._createSetOfMicrographs(suffix='_ts_average')

        for idx in avgIndexList:
            acqOrder = acqOrderList[idx]
            angle = tiltAngleList[idx]
            tsAvg = Micrograph()

            outputFilePath = os.path.join(extraPrefix,
                                          firstItem.parseFileName(suffix="_" + str(idx+1),
                                                                  extension=".mrc"))

            tsAvg.setFileName(outputFilePath)
            tsAvg.setSamplingRate(firstItem.getSamplingRate())

            setattr(tsAvg, TiltImage.TS_ID_FIELD, String(tsId))
            setattr(tsAvg, TiltImage.ACQ_ORDER_FIELD , Integer(acqOrder))
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
    
    @staticmethod
    def getAcqOrderList(ts):
        acqOrderList = []

        for ti in ts:
            acqOrderList.append(ti.getAcquisitionOrder())

        return acqOrderList

    # --------------------------- INFO functions ----------------------------
