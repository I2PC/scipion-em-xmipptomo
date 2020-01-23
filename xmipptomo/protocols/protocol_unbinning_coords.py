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

from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import PointerParam, FloatParam
from tomo.protocols import ProtTomoBase

class XmippProtUnbinningCoord(EMProtocol, ProtTomoBase):
    """ This protocol takes a set of coordinates and multiplies them by a binning factor to get the coordinates of the
     unbinning tomogram."""

    _label = 'unbinning coordinates'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input coordinates')
        form.addParam('inputCoordinates', PointerParam, label="Input Coordinates",
                      pointerClass='SetOfCoordinates3D', help='Select the SetOfCoordinates3D.')
        form.addParam('factor', FloatParam, label='Binning factor', help='binning factor')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('unbinningCoords')

    # --------------------------- STEPS functions -------------------------------
    def unbinningCoords(self):
        inputSet = self.inputCoordinates.get()
        outputSet = self._createSetOfCoordinates3D(inputSet.getPrecedents())
        outputSet.copyInfo(inputSet)
        outputSet.setBoxSize(inputSet.getBoxSize())
        for coord in inputSet:
            coord.setX(coord.getX() * self.factor.get())
            coord.setY(coord.getY() * self.factor.get())
            coord.setZ(coord.getZ() * self.factor.get())
            outputSet.append(coord)
        self._defineOutputs(outputCoordinates=outputSet)
        self._defineSourceRelation(inputSet, outputSet)

    # --------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        summary.append("")
        return summary

    def _methods(self):
        methods = []
        methods.append("")
        return methods
