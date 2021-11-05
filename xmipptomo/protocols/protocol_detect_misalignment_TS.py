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
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
import xmipptomo.utils as utils

class XmippProtDetectMisalignmentTiltSeries(EMProtocol, ProtTomoBase):
    """
    Scipion protocol for xmipp_tomo_detect_misalignment_trajectory. Detect misalignment in a tilt series.
    """

    _label = 'detect misaligned tilt-series '
    _devStatus = BETA

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

        form.addParam('fiducialSize',
                      params.FloatParam,
                      important=True,
                      label='Fiducial size (A)',
                      help='Fiducial size in Angstroms (A).')

        form.addParam('sdThreshold',
                      params.FloatParam,
                      advanced=True,
                      label='Coordinate value SD threshold',
                      help='Number of SD a coordinate value must be over the mean to consider that it belongs to a '
                           'high contrast feature.')

        form.addParam('numberOfCoordinatesThr',
                      params.IntParam,
                      advanced=True,
                      label='Number of coordinates threshold',
                      help='Minimum number of coordinates attracted to a center of mass to consider it as a high '
                           'contrast feature.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            tsObjId = ts.getObjId()
            self._insertFunctionStep(self.convertInputStep, tsObjId)
            self._insertFunctionStep(self.detectMisalignment, ts.getObjId())

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
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))
        utils.writeXmippTiltAngleList(ts, angleFilePath)

    def detectMisalignment(self, tsObjId):
        pass
