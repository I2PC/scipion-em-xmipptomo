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
from pyworkflow.object import String
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
import xmipptomo.utils as utils


class XmippProtAverageViewTiltSeries(EMProtocol, ProtTomoBase):
    """
    Scipion protocol to calculate the 3D coordinate based on a set of 2D coordinates obtained from the tilt-series
    using a SPA 2D picker.
    """

    _label = 'Tilt-series calculate coordinates 3D'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtTomoBase.__init__(self)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfCoordinates',
                      PointerParam,
                      pointerClass='SetOfCoordinate',
                      important=True,
                      label='Input set of 2D coordinates')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):


    # --------------------------- STEPS functions ----------------------------

