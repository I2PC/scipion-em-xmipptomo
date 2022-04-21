# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csi.es) [1]
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
from pyworkflow.object import Set
import tomo.objects as tomoObj


class XmippProtResizeBase(EMProtocol, ProtTomoBase):
    """
    Base class to resize 2D and 3D objects in tomography
    """

    RESIZE_SAMPLINGRATE = 0
    RESIZE_FACTOR = 1
    RESIZE_PYRAMID = 2

    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    @classmethod
    def _defineParamsReSize(self, form):
        form.addParam('resizeOption',
                      EnumParam,
                      choices=['Sampling Rate', 'Factor', 'Pyramid'],
                      default=self.RESIZE_SAMPLINGRATE,
                      label="Resize option",
                      isplay=EnumParam.DISPLAY_COMBO,
                      help='Select an option to resize the images: \n '
                           '_Sampling Rate_: Set the desire sampling rate to resize. \n'
                           '_Factor_: Set a resize factor to resize. \n '
                           '_Pyramid_: Use positive level value to expand and negative to reduce. \n'
                           'Pyramid uses spline pyramids for the interpolation. All the rest uses normally \n'
                           'interpolation (cubic B-spline or bilinear interpolation).')

        form.addParam('resizeSamplingRate',
                      FloatParam,
                      default=1.0,
                      condition='resizeOption==%d' % self.RESIZE_SAMPLINGRATE,
                      label='Resize sampling rate (Å/px)',
                      help='Set the new output sampling rate.')

        form.addParam('resizeFactor',
                      FloatParam,
                      default=0.5,
                      condition='resizeOption==%d' % self.RESIZE_FACTOR,
                      label='Resize factor',
                      help='New size is the old one x resize factor.')

        form.addParam('resizeLevel',
                      IntParam,
                      default=0,
                      condition='resizeOption==%d' % self.RESIZE_PYRAMID,
                      label='Pyramid level',
                      help='Use positive value to expand and negative to reduce.')

        form.addParam('hugeFile', BooleanParam, default=False, expertLevel=const.LEVEL_ADVANCED,
                      label='Huge file',
                      help='If the file is huge, very likely you may have problems doing the antialiasing filter '
                           '(because there is no memory for the input and its Fourier transform). This option '
                           'removes the antialiasing filter (meaning you will get aliased results), and performs '
                           'a bilinear interpolation (to avoid having to produce the B-spline coefficients).')

    # --------------------------- UTILS functions -------------------------------
    def resizeCommonArgsResize(cls, protocol, samplingRate):
        if protocol.resizeOption == cls.RESIZE_SAMPLINGRATE:
            newSamplingRate = protocol.resizeSamplingRate.get()
            factor = samplingRate / newSamplingRate
            args = " --factor %f" % factor
        elif protocol.resizeOption == cls.RESIZE_FACTOR:
            factor = protocol.resizeFactor.get()
            newSamplingRate = samplingRate / factor
            args = " --factor %f" % factor
        elif protocol.resizeOption == cls.RESIZE_PYRAMID:
            level = protocol.resizeLevel.get()
            factor = 2 ** level
            newSamplingRate = samplingRate / factor
            args = " --pyramid %d" % level
        if protocol.hugeFile:
            args += " --interp linear"

        protocol.samplingRate = newSamplingRate
        protocol.samplingRateOld = samplingRate
        protocol.factor = factor

        return args

class XmippProtCropBase2D(EMProtocol, ProtTomoBase):
    """
    Base class to crop 2D objects in tomography
    """
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    @classmethod
    def _defineParamsReSize2D(self, form):
        linex = form.addLine('Keep pixels along X',
                             help='Select the number of pixels to be cropped from the borders'
                                  ' along the x axis')
        linex.addParam('xcrop0', FloatParam, default=100, label='from')
        linex.addParam('xcropF', FloatParam, default=100, label='to')

        liney = form.addLine('Keep pixels along Y',
                             help='Select the number of pixels to be cropped from the borders'
                                  ' along the y axis')
        liney.addParam('ycrop0', FloatParam, default=100, label='from')
        liney.addParam('ycropF', FloatParam, default=100, label='to')


class XmippProtCropBase3D(XmippProtCropBase2D):
    """
    Base class to crop 3D objects in tomography
    """
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    @classmethod
    def _defineParamsReSize3D(self, form):

        linez = form.addLine('Keep pixels along Z',
                             help='Select the number of pixels to be cropped from the borders'
                                  ' along the z axis')
        linez.addParam('zcrop0', FloatParam, default=100, label='from')
        linez.addParam('zcropF', FloatParam, default=100, label='to')
