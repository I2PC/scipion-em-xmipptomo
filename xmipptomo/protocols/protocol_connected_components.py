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

from os.path import basename
import numpy as np
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, FloatParam
from tomo.protocols import ProtTomoBase


class XmippProtConnectedComponents(EMProtocol, ProtTomoBase):
    """ This protocol takes a set of coordinates and identifies connected
    components among the picked particles."""

    _label = 'connected components'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input coordinates')
        form.addParam('inputCoordinates', PointerParam, label="Input Coordinates",
                      pointerClass='SetOfCoordinates3D', help='Select the SetOfCoordinates3D.')
        form.addParam('distance', FloatParam, label='Distance', help='Maximum radial distance (in voxels) between '
                                                                      'particles to consider that they are in the same '
                                                                      'connected component. Wizard returns three times '
                                                                      'the box size of the input coordinates.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeConnectedComponentsStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def computeConnectedComponentsStep(self):
        inputCoors = self.inputCoordinates.get()
        coorSetList = []
        tomoList = []

        # Create a separate setOfCoordinates for each tomogram
        for coor in inputCoors.iterCoordinates():
            tomoName = coor.getVolName()
            if tomoName not in tomoList:
                tomoList.append(tomoName)
                tomocoorset = self._createSetOfCoordinates3D(inputCoors.getPrecedents(), '_' + basename(tomoName))
                coorSetList.append(tomocoorset)
            idx = tomoList.index(tomoName)
            coorSetList[idx].append(coor)

        # For each tomogram coordinates, perform "connected components" logic
        outputsetIndex = 0
        self.outputSet = self._createSetOfCoordinates3D(inputCoors.getPrecedents())
        self.outputSet.copyInfo(inputCoors)
        self.outputSet.setBoxSize(inputCoors.getBoxSize())
        self.outputSet.setSamplingRate(inputCoors.getSamplingRate())
        for coorSet in coorSetList:
            coorSet.write()
            minDist = self.distance.get()
            coorlist = []
            for i, coor in enumerate(coorSet.iterItems()):
                coorlist.append([coor.getX(), coor.getY(), coor.getZ()])
            A = np.zeros([len(coorlist), len(coorlist)])
            for j, coor1 in enumerate(coorlist):
                for k, _ in enumerate(coorlist, start=j+1):
                    if k == len(coorlist):
                        break
                    else:
                        coor2 = coorlist[k]
                        if abs(coor1[0]-coor2[0]) <= minDist and abs(coor1[1]-coor2[1]) <= minDist \
                                and abs(coor1[2]-coor2[2]) <= minDist:
                            A[j, k] = 1
                            A[k, j] = 1
            tomoNameshort = basename(coorSet.getFirstItem().getVolName())
            np.savetxt(self._getExtraPath('adjacency_matrix_%s' % tomoNameshort), A)
            D = np.diag(A.sum(axis=1))
            np.savetxt(self._getExtraPath('degree_matrix_%s' % tomoNameshort), D)
            L = D - A
            np.savetxt(self._getExtraPath('laplacian_matrix_%s' % tomoNameshort), L)
            vals, vecs = np.linalg.eig(L)
            vals = vals.real
            vecs = vecs.real
            np.savetxt(self._getExtraPath('eigenvecs_matrix_%s' % tomoNameshort), vecs)
            np.savetxt(self._getExtraPath('eigenvalues_matrix_%s' % tomoNameshort), vals)
            vals0list = [i for i, x in enumerate(vals.real) if abs(x) < (1/(np.sqrt(len(coorSet)))*1e-3)]
            listOfSets = [[] for x in range(len(vals0list))]
            for k in range(len(coorSet)):
                row = []
                for j in vals0list:
                    row.append(vecs[k, j])
                row = [abs(number)for number in row]
                ixMax = np.argmax(row)
                listOfSets[ixMax].append(k)

            for ix, coorInd in enumerate(listOfSets):
                if len(coorInd) != 0:
                    outputsetIndex += 1
                    for id, coor3D in enumerate(coorSet.iterItems()):
                        if id in coorInd:
                            coor3D.setGroupId(outputsetIndex)
                            self.outputSet.append(coor3D)

    def createOutputStep(self):
        self._defineOutputs(outputSetOfCoordinates3D=self.outputSet)
        self._defineSourceRelation(self.inputCoordinates, self.outputSet)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        return validateMsgs

    def _summary(self):
        summary = []
        summary.append("Maximum radial distance between particles in the same connected component: %d voxels"
                       % (self.distance.get()))
        return summary

    def _methods(self):
        methods = []
        methods.append("Coordinates grouped in connected component with a maximum radial distance of %d voxels."
                       % self.distance.get())
        return methods

