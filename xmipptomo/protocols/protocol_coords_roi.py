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
from pyworkflow.protocol.params import PointerParam, EnumParam, IntParam
from tomo.protocols import ProtTomoBase


class XmippProtCCroi(EMProtocol, ProtTomoBase):
    """ This protocol adjust a SetOfCoordinates (which usually will come from a connected componnent) to a ROI (region
    of interest) previously defined"""

    _label = 'coordinates to roi'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input coordinates')
        form.addParam('inputCoordinates', PointerParam, label="Input Coordinates",
                      pointerClass='SetOfCoordinates3D', help='Select the SetOfCoordinates3D.')
        form.addParam('inputMesh', PointerParam, label="Input mesh",
                      pointerClass='Mesh, SetOfCoordinates', help='Select the mesh')  # REMOVE SETOFCOORD!!!! (just for the test to work)
        form.addParam('selection', EnumParam, choices=['Whole cc', 'Points in roi'], default=0, label='Selection',
                      display=EnumParam.DISPLAY_HLIST,
                      help='Selection options:\n*Whole cc*: It takes the whole connected componnent (cc) if all the '
                           'points in the cc belongs to the ROI. If a "Number of points" is introduced in the following'
                           ' field, the whole cc will be taken if that number of points from the cc belongs to the ROI.'
                           '\n*Points in roi*: It takes just the points of the cc which belongs to the roi')
        form.addParam('points', IntParam, label="Number of points", condition='selection == 0', allowsNull=True,
                      help='see "Selection" help')
        form.addParam('distance', IntParam, label='Distance', default=0,
                      help='Maximum radial distance (in pixels) between mesh vertex and a coordinate to consider that '
                           'it belongs to the ROI. Wizard returns three times the box size of the input coordinates.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeDistances')

    # --------------------------- STEPS functions -------------------------------
    def computeDistances(self):
        inputCoor = self.inputCoordinates.get()
        inputMesh = self.inputMesh.get()  # for REAL MESHES, inputMesh.getMesh() + .iterItems()???
        distance = self.distance.get()
        outputSet = self._createSetOfCoordinates3D(inputCoor.getPrecedents())
        outputSet.copyInfo(inputCoor)
        outputSet.setBoxSize(inputCoor.getBoxSize())
        for coorm in inputMesh.iterItems():  # TODO: This works if mesh is a SetOfCoord
            for coorc in inputCoor.iterItems():
                if abs(coorm.getX() - coorc.getX()) <= distance and abs(coorm.getY() - coorc.getY()) <= distance \
                        and abs(coorm.getZ() - coorc.getZ()) <= distance:
                    outputSet.append(coorc)
        # selection options??
        self._defineOutputs(outputCoordinates=outputSet)
        self._defineSourceRelation(inputCoor, outputSet)

    # --------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        summary.append("")
        return summary

    def _methods(self):
        methods = []
        methods.append("")
        return methods

