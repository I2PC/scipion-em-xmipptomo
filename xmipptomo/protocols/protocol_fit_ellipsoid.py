# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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
import numpy as np
from pyworkflow.protocol.params import PointerParam
import pyworkflow.utils as pwutlis
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.objects import Mesh


class XmippProtFitEllipsoid(EMProtocol, ProtTomoBase):
    """ This protocol adjust a SetOfCoordinates to an ellipsoid (for example, a vesicle), defining a region of interest
    (ROI) as output."""

    _label = 'fit vesicles'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input coordinates')
        form.addParam('inputCoordinates', PointerParam, label="Coordinates",
                      pointerClass='SetOfCoordinates3D', help='Select the coordinates.')
        form.addParam('inputTomos', PointerParam, label="Tomograms",
                      pointerClass='SetOfTomograms', help='Select tomograms.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fitEllipsoidStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def fitEllipsoidStep(self):
        pass

    def createOutputStep(self):
        outSet = self._createSetOfMeshes()
        for file in os.listdir(self._getExtraPath()):
            if file.endswith(".txt"):
                data = np.loadtxt(self._getExtraPath(file), delimiter=',')
                groups = np.unique(data[:, 3]).astype(int)
                for group in groups:
                    mesh = Mesh(group=group, path=self._getExtraPath(file))
                    for tomo in self.inputTomos.get().iterItems():
                        if file[:-4] == pwutlis.removeBaseExt(tomo.getFileName()):
                            mesh.setVolume(tomo.clone())
                    outSet.append(mesh)
        outSet.setVolumes(self.inputTomos.get())
        self._defineOutputs(outputMeshes=outSet)
        self._defineSourceRelation(self.inputTomos.get(), outSet)

    # --------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        if not os.listdir(self._getExtraPath()):
            summary.append("Output vesicles not ready yet.")
        else:
            summary.append("")
        return summary

    def _methods(self):
        return ["Fit an ellipsoid (vesicle) into a 3D set of points."]
