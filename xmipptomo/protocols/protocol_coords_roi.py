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

import numpy as np
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import MultiPointerParam, IntParam, PointerParam
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
        form.addParam('inputCoordinates', MultiPointerParam, label="Input connected components",
                      pointerClass='SetOfCoordinates3D', help='Select the Connected components (SetOfCoordinates3D).')
        form.addParam('inputMeshes', PointerParam, label="Input ROIs",
                      pointerClass='SetOfMeshes', help='Select the ROIs (Regions Of Interest)')
        # form.addParam('selection', EnumParam, choices=['Whole cc', 'Points in roi'], default=0, label='Selection',
        #               display=EnumParam.DISPLAY_HLIST,
        #               help='Selection options:\n*Whole cc*: It takes the whole connected component (cc) if all the '
        #                  'points in the cc belongs to the ROI. If a "Number of points" is introduced in the following'
        #                  ' field, the whole cc will be taken if that number of points from the cc belongs to the ROI.'
        #                    '\n*Points in roi*: It takes just the points of the cc which belongs to the roi')
        form.addParam('points', IntParam, label="Percentage of coordinates in ROI",  default=80,
                      # condition='selection == 0',
                      allowsNull=True, help='Percentage of coordinates from a connected component that should be inside'
                                            ' the ROI to consider that connected component.')
        form.addParam('distance', IntParam, label='Distance', default=128,
                      help='Maximum euclidean distance (in pixels) between ROI vertex and a coordinate to consider that'
                           ' it belongs to the ROI. Wizard returns three times the box size of the input coordinates.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeDistances')

    # --------------------------- STEPS functions -------------------------------
    def computeDistances(self):
        for ix, inputSetCoor in enumerate(self.inputCoordinates):
            i = 0
            perc = self._percentage(inputSetCoor)
            for mesh in self.inputMeshes.get().iterItems():
                if inputSetCoor.get().getPrecedents().getFirstItem().getFileName() == mesh.getVolume().getFileName():
                    for coorcc in inputSetCoor.get():
                        for coormesh in mesh.getMesh():
                            if self._euclideanDistance(coorcc, coormesh) <= self.distance.get():
                                i += 1
                                break
            print("---------", i)
            if i >= perc and i<= inputSetCoor.get().getSize():
                outputSet = self._createSetOfCoordinates3D(inputSetCoor.get().getPrecedents(), ix+1)
                outputSet.copyInfo(inputSetCoor)
                outputSet.copyItems(inputSetCoor.get())
                outputSet.setBoxSize(inputSetCoor.get().getBoxSize())
                outputSet.setSamplingRate(inputSetCoor.get().getSamplingRate())
                name = 'output3DCoordinates%s' % str(ix+1)
                args = {}
                args[name] = outputSet
                self._defineOutputs(**args)
                self._defineSourceRelation(inputSetCoor, outputSet)

    # --------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        summary.append("")
        return summary

    def _methods(self):
        methods = []
        methods.append("")
        return methods

    # --------------------------- UTILS functions --------------------------------------------
    def _percentage(self, inputSetCoor):
        return (self.points.get()*inputSetCoor.get().getSize())/100

    def _euclideanDistance(self, coorcc, cm):
        return np.sqrt((coorcc.getX() - int(cm[0])) + (coorcc.getY() - int(cm[1])) + (coorcc.getZ() - int(cm[2])))
