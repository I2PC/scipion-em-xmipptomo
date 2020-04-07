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

from pyworkflow.object import Set
from pyworkflow.protocol.params import MultiPointerParam, IntParam, PointerParam, EnumParam
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase


class XmippProtCCroi(EMProtocol, ProtTomoBase):
    """ This protocol adjust a SetOfCoordinates (which usually will come from a
    connected componnent) to a ROI (region of interest) previously defined"""

    _label = 'connected components to ROIs'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input coordinates')
        form.addParam('inputCoordinates', MultiPointerParam, label="Input connected components",
                      pointerClass='SetOfCoordinates3D', help='Select the Connected components.')
        form.addParam('inputMeshes', PointerParam, label="Input ROIs",
                      pointerClass='SetOfMeshes', help='Select the ROIs (Regions Of Interest)')
        form.addParam('selection', EnumParam, choices=['Connected component', 'Points'], default=0, label='Selection',
                      display=EnumParam.DISPLAY_HLIST,
                      help='Selection options:\n*Connected component*: It takes the whole connected component (cc) if '
                           'a percentage of the points (introduced in the next field) in the cc belongs to the ROI. '
                           '\n*Points*: It takes just the points of the cc which belongs to the roi')
        form.addParam('points', IntParam, label="Percentage of coordinates in ROI",  default=80,
                      condition='selection == 0',
                      allowsNull=True, help='Percentage of coordinates from a connected component that should be inside'
                                            ' the ROI to consider that connected component.')
        form.addParam('distance', IntParam, label='Distance', default=50,
                      help='Maximum euclidean distance (in pixels) between ROI vertex and a coordinate to consider that'
                           ' it belongs to the ROI.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeDistances')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def computeDistances(self):
        sel = self.selection.get()
        for ix, inputSetCoor in enumerate(self.inputCoordinates):
            perc = self._percentage(inputSetCoor)
            for mesh in self.inputMeshes.get().iterItems():
                if os.path.basename(inputSetCoor.get().getFirstItem().getVolName()) == os.path.basename(mesh.getVolume().getFileName()):
                    if sel == 0:
                        i = 0
                    else:
                        outputSetList = []
                    for coorcc in inputSetCoor.get():
                        for coormesh in mesh.getMesh():
                            if self._euclideanDistance(coorcc, coormesh) <= self.distance.get():
                                if sel == 0:
                                    i += 1
                                else:
                                    outputSetList.append(coorcc.getObjId())
                                break
                    if sel == 0:
                        if i >= perc:
                            outputSet = self._createSetOfCoordinates3D(inputSetCoor.get().getPrecedents(), ix + 1)
                            outputSet.copyInfo(inputSetCoor)
                            outputSet.copyItems(inputSetCoor.get())
                            outputSet.setBoxSize(inputSetCoor.get().getBoxSize())
                            outputSet.setSamplingRate(inputSetCoor.get().getSamplingRate())
                            name = 'output3DCoordinates%s' % str(ix+1)
                            args = {}
                            args[name] = outputSet
                            outputSet.setStreamState(Set.STREAM_OPEN)
                            self._defineOutputs(**args)
                            self._defineSourceRelation(inputSetCoor, outputSet)

                    else:
                        if len(outputSetList) != 0:
                            outputSet = self._createSetOfCoordinates3D(inputSetCoor.get().getPrecedents(), ix + 1)
                            outputSet.copyInfo(inputSetCoor)
                            outputSet.setBoxSize(inputSetCoor.get().getBoxSize())
                            outputSet.setSamplingRate(inputSetCoor.get().getSamplingRate())
                            for coor3D in inputSetCoor.get().iterItems():
                                if coor3D.getObjId() in outputSetList:
                                    outputSet.append(coor3D)
                            name = 'output3DCoordinates%s' % str(ix + 1)
                            args = {}
                            args[name] = outputSet
                            outputSet.setStreamState(Set.STREAM_OPEN)
                            self._defineOutputs(**args)
                            self._defineSourceRelation(inputSetCoor, outputSet)

    def createOutputStep(self):
        for outputset in self._iterOutputsNew():
            outputset[1].setStreamState(Set.STREAM_CLOSED)
        self._store()

    # --------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        if self.selection.get() == 0:
            summary.append("Percentage of coordinates in ROI: %d" % self.points.get())
        summary.append("Max distance to ROI: %d\nConnected components in ROIs: %d"
                       % (self.distance.get(), len(self._outputs)))
        return summary

    def _methods(self):
        methods = []
        methods.append("%d connected components detected" % len(self._outputs))
        if self.selection.get() == 0:
            methods.append("with at least %d percent of points" % self.points.get())
        methods.append("at a maximun distance of %d pixels of a ROI." % self.distance.get())
        return methods

    # --------------------------- UTILS functions --------------------------------------------
    def _percentage(self, inputSetCoor):
        return (self.points.get()*inputSetCoor.get().getSize())/100

    def _euclideanDistance(self, coorcc, cm):
        return np.sqrt((coorcc.getX() - int(cm[0])) + (coorcc.getY() - int(cm[1])) + (coorcc.getZ() - int(cm[2])))
