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


class XmippProtCalculate3dCoordinatesFromTS(EMProtocol, ProtTomoBase):
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
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label='Input set of 2D coordinates form tilt-series')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)

    # --------------------------- STEPS functions ----------------------------

    def convertInputStep(self):
        print(self.inputSetOfCoordinates.get().getMicrographs())
        for c in self.inputSetOfCoordinates.get():
            print(c)
            print(c.getMicId())
        # """ Read the input metadatata. """
        # # Get the converted input micrographs in Xmipp format
        # makePath(self._getExtraPath("untilted"))
        # makePath(self._getExtraPath("tilted"))
        #
        # uSet = self.untiltedSet.get()
        # tSet = self.tiltedSet.get()
        #
        # # Get the untilted and tilted coordinates, depending on the input type
        # if isinstance(uSet, SetOfParticles):
        #     uCoords = uSet.getCoordinates()
        #     tCoords = tSet.getCoordinates()
        #
        #     # If there are not Coordinates associated to particles
        #     # we need to create and fill the set of coordinates
        #     if uCoords is None or tCoords is None:
        #         micTiltedPairs = self.inputMicrographsTiltedPair.get()
        #         uCoords = self._coordsFromParts(micTiltedPairs.getUntilted(),
        #                                         uSet, '_untilted')
        #         tCoords = self._coordsFromParts(micTiltedPairs.getTilted(),
        #                                         tSet, '_tilted')
        # else:
        #     uCoords = uSet
        #     tCoords = tSet
        #
        # writeSetOfCoordinates(self._getExtraPath("untilted"), uCoords)
        # writeSetOfCoordinates(self._getExtraPath("tilted"), tCoords)


    def assignmentStep(self, fnuntilt, fntilt, fnmicsize, fnposUntilt, fnposTilt):
        params = ' --untiltcoor %s' % fnuntilt
        params += ' --tiltcoor %s' % fntilt
        params += ' --tiltmicsize %s' % fnmicsize
        params += ' --maxshift %f' % self.maxShift
        params += ' --particlesize %d' % self._getBoxSize()
        params += ' --threshold %f' % self.threshold
        params += ' --odir %s' % self._getExtraPath()
        self.runJob('xmipp_image_assignment_tilt_pair', params)

        # Estimate the tilt axis
        params = ' --untilted %s' % fnposUntilt
        params += ' --tilted %s' % fnposTilt
        params += ' -o %s' % self._getPath('input_micrographs.xmd')
        self.runJob('xmipp_angular_estimate_tilt_axis', params)