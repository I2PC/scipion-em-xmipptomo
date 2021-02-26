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
from pyworkflow.protocol.params import IntParam, PointerParam, EnumParam
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
import tomo.constants as const


class XmippProtCCroi(EMProtocol, ProtTomoBase):
    """ This protocol adjust a SetOfCoordinates (which usually will come from a
    connected componnent) to a ROI (region of interest) previously defined"""

    _label = 'connected components to ROIs'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', PointerParam, label="Connected components",
                      pointerClass='SetOfCoordinates3D', help='Select the Connected components (SetOfCoordinates3D).')
        form.addParam('inputMeshes', PointerParam, label="ROIs",
                      pointerClass='SetOfMeshes', help='Select the ROIs (Regions Of Interest) they are SetOfMeshes')
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
        """Compare all connected components with all meshes; a connected component can be just in one mesh
        (first one found)"""

        inputSetCoor = self.inputCoordinates.get()
        inputSetMeshes = self.inputMeshes.get()

        # group connected components by groupId
        listOfCCs = []
        listOfgroupIdCCs = []
        for coor in inputSetCoor.iterCoordinates():
            tomoNameC = coor.getVolName()
            groupIdC = coor.getGroupId()
            tupleC = (groupIdC, tomoNameC)
            if tupleC not in listOfgroupIdCCs:
                listOfgroupIdCCs.append(tupleC)
                listOfCCs.append(self._createSetOfCoordinates3D(inputSetCoor.getPrecedents(), groupIdC))
            coorLoc = listOfCCs[listOfgroupIdCCs.index(tupleC)]
            if not coorLoc:
                coorLoc.append(coor)
            else:
                if coorLoc.getFirstItem().getVolName() == tomoNameC:
                    coorLoc.append(coor)

        # group meshes by groupId and tomoName
        listOfMeshes = []
        listOfgroupIdMeshes = []
        for meshPoint in inputSetMeshes.iterCoordinates():
            tomoNameM = meshPoint.getVolName()
            tomoId = meshPoint.getVolume().getObjId()
            groupIdM = meshPoint.getGroupId()
            tupleM = (groupIdM, tomoNameM)
            if tupleM not in listOfgroupIdMeshes:
                listOfgroupIdMeshes.append(tupleM)
                listOfMeshes.append(self._createSetOfMeshes(inputSetMeshes.getPrecedents(), '%d_%d' % (tomoId,groupIdM)))
            meshPointLoc = listOfMeshes[listOfgroupIdMeshes.index(tupleM)]
            if not meshPointLoc:
                meshPointLoc.append(meshPoint)
            else:
                if meshPointLoc.getFirstItem().getVolName() == tomoNameM:
                    meshPointLoc.append(meshPoint)

        sel = self.selection.get()
        self.outputSet = self._createSetOfCoordinates3D(inputSetCoor.getPrecedents())
        self.outputSet.copyInfo(inputSetCoor)
        self.outputSet.setBoxSize(inputSetCoor.getBoxSize())
        self.outputSet.setSamplingRate(inputSetCoor.getSamplingRate())

        for cc in listOfCCs:
            for mesh in listOfMeshes:
                if os.path.basename(cc.getFirstItem().getVolName()) == os.path.basename(mesh.getFirstItem().getVolName()):
                    if sel == 0:
                        i = 0
                    else:
                        outputSetList = []
                    for coorcc in cc.iterCoordinates():
                        for meshPoint in mesh.iterCoordinates():
                            if self._euclideanDistance(coorcc, meshPoint) <= self.distance.get():
                                if sel == 0:
                                    i += 1
                                else:
                                    outputSetList.append(coorcc.getObjId())
                                break

                    if sel == 0:
                        perc = self._percentage(cc)
                        if i >= perc:
                            self.groupIdCoorCC = coorcc.getGroupId()
                            self.outputSet.copyItems(inputSetCoor, updateItemCallback=self._updateItem)

                    else:
                        if len(outputSetList) != 0:
                            for coor3D in inputSetCoor.iterItems():
                                if coor3D.getObjId() in outputSetList:
                                    self.outputSet.append(coor3D)
                break

    def createOutputStep(self):
        self._defineOutputs(outputSet=self.outputSet)
        self._defineSourceRelation(self.inputCoordinates, self.outputSet)

    # --------------------------- INFO functions --------------------------------
    def _summary(self):
        summary = []
        if self.selection.get() == 0:
            summary.append("Percentage of coordinates in ROI: %d" % self.points.get())
        if self.selection.get() == 0:
            sel = 'complete connected components'
        else:
            sel = 'points in connected components'
        summary.append("Max distance to ROI: %d\nSelect: %s"
                       % (self.distance.get(), sel))
        return summary

    def _methods(self):
        methods = "Connected components detected"
        if self.selection.get() == 0:
            methods.append("with at least %d percent of points" % self.points.get())
        methods.append("at a maximun distance of %d pixels of a ROI." % self.distance.get())
        return methods

    # --------------------------- UTILS functions --------------------------------------------
    def _percentage(self, inputSetCoor):
        return (self.points.get()*inputSetCoor.getSize())/100

    def _euclideanDistance(self, coorcc, cm):
        return np.sqrt((coorcc.getX(const.BOTTOM_LEFT_CORNER) - int(cm.getX(const.BOTTOM_LEFT_CORNER))) +
                       (coorcc.getY(const.BOTTOM_LEFT_CORNER) - int(cm.getY(const.BOTTOM_LEFT_CORNER))) +
                       (coorcc.getZ(const.BOTTOM_LEFT_CORNER) - int(cm.getZ(const.BOTTOM_LEFT_CORNER))))

    def _updateItem(self, item, row):
        if item.getGroupId() != self.groupIdCoorCC:
            setattr(item, "_appendItem", False)
