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
from pyworkflow.object import Set, List, String, Boolean
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
import xmipptomo.utils as utils
import emtable

METADATA_INPUT_COORDINATES = "fiducialCoordinates.xmd"
VMC_FILE_NAME = "vCM.xmd"


class XmippProtdetectLandmarkTS(EMProtocol, ProtTomoBase):
    """
    Scipion protocol for xmipp_tomo_detect_landmarks. Detect landmarks in a tilt series.
    """

    _label = 'detect landmarks TS'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtTomoBase.__init__(self)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('landmarkSize',
                      params.FloatParam,
                      important=True,
                      label='Landmark size (nm)',
                      help='Landmark size in nanometers (nm).')

        form.addParam('thrSD',
                      params.FloatParam,
                      advanced=True,
                      default=3,
                      label='Coordinate value SD threshold',
                      help='Number of SD a coordinate value must be over the mean to consider that it belongs to a '
                           'high contrast feature.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):


    # --------------------------- STEPS functions ----------------------------


    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfLandmarkModels(self):
        if hasattr(self, "outputSetOfLandmarkModels"):
            self.outputSetOfLandmarkModels.setStreamState(Set.STREAM_OPEN)
            self.outputSetOfLandmarkModels.enableAppend()

        else:
            outputSetOfLandmarkModels = self._createSetOfLandmarkModels("_ali")

            outputSetOfLandmarkModels.copyInfo(self.inputSetOfTiltSeries.get())

            outputSetOfLandmarkModels.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfLandmarkModels=outputSetOfLandmarkModels)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfLandmarkModels)

        return self.outputSetOfLandmarkModels

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []

        return methods
