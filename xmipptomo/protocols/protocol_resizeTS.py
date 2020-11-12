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
from pyworkflow.protocol.params import EnumParam, FloatParam, IntParam, BooleanParam
import pyworkflow.protocol.constants as const
from tomo.protocols import ProtTomoBase


class XmippProtResizeTiltSeries(EMProtocol, ProtTomoBase):
    """
    Wrapper protocol to Xmipp image resize applied on tilt-series
    """

    RESIZE_SAMPLINGRATE = 0
    RESIZE_DIMENSIONS = 1
    RESIZE_FACTOR = 2
    RESIZE_PYRAMID = 3

    _label = 'resize tilt-series'

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('resizeOption',
                      EnumParam,
                      choices=['Sampling Rate', 'Dimensions', 'Factor', 'Pyramid'],
                      condition='doResize',
                      default=self.RESIZE_SAMPLINGRATE,
                      label="Resize option",
                      isplay=EnumParam.DISPLAY_COMBO,
                      help='Select an option to resize the images: \n '
                           '_Sampling Rate_: Set the desire sampling rate to resize. \n'
                           '_Dimensions_: Set the output dimensions. Resize operation can be done in Fourier space.\n'
                           '_Factor_: Set a resize factor to resize. \n '
                           '_Pyramid_: Use positive level value to expand and negative to reduce. \n'
                           'Pyramid uses spline pyramids for the interpolation. All the rest uses normally \n'
                           'interpolation (cubic B-spline or bilinear interpolation). If you set the method to \n'
                           'dimensions, you may choose between interpolation and Fourier cropping.')

        form.addParam('resizeSamplingRate',
                      FloatParam,
                      default=1.0,
                      condition='doResize and resizeOption==%d' % self.RESIZE_SAMPLINGRATE,
                      label='Resize sampling rate (Å/px)',
                      help='Set the new output sampling rate.')

        form.addParam('doFourier',
                      BooleanParam,
                      default=False,
                      condition='doResize and resizeOption==%d' % self.RESIZE_DIMENSIONS,
                      label='Use fourier method to resize?',
                      help='If you set to *True*, the final dimensions must be lower than the original ones.')

        form.addParam('resizeDim',
                      IntParam,
                      default=0,
                      condition='doResize and resizeOption==%d' % self.RESIZE_DIMENSIONS,
                      label='New image size (px)',
                      help='Size in pixels of the particle images <x> <y=x> <z=x>.')

        form.addParam('resizeFactor',
                      FloatParam,
                      default=0.5,
                      condition='doResize and resizeOption==%d' % self.RESIZE_FACTOR,
                      label='Resize factor',
                      help='New size is the old one x resize factor.')

        form.addParam('resizeLevel',
                      IntParam,
                      default=0,
                      condition='doResize and resizeOption==%d' % self.RESIZE_PYRAMID,
                      label='Pyramid level',
                      help='Use positive value to expand and negative to reduce.')

        form.addParam('hugeFile', BooleanParam, default=False, expertLevel=const.LEVEL_ADVANCED,
                      label='Huge file',
                      help='If the file is huge, very likely you may have problems doing the antialiasing filter '
                           '(because there is no memory for the input and its Fourier tranform). This option '
                           'removes the antialiasing filter (meaning you will get aliased results), and performs '
                           'a bilinear interpolation (to avoid having to produce the B-spline coefficients).')
