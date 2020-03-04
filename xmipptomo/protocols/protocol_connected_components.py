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
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import PointerParam, FloatParam
from tomo.protocols import ProtTomoBase


class XmippProtConnectedComponents(EMProtocol, ProtTomoBase):
    """ This protocol takes a set of coordinates and identifies connected components among the picked particles."""

    _label = 'connected components'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input coordinates')
        form.addParam('inputCoordinates', PointerParam, label="Input Coordinates",
                      pointerClass='SetOfCoordinates3D', help='Select the SetOfCoordinates3D.')
        form.addParam('distance', FloatParam, label='Distance1', help='Maximum radial distance (in voxels) between '
                                                                      'particles to consider that they are in the same '
                                                                      'connected component. Wizard returns three times '
                                                                      'the box size of the input coordinates.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeConnectedComponentsStep')

    # --------------------------- STEPS functions -------------------------------
    def computeConnectedComponentsStep(self):
        self.tomoList = []
        for coor in self.inputCoordinates.get().iterItems():
            coorName = coor.getVolName()
            if coorName not in self.tomoList:
                self.tomoList.append(coorName)
        fi = 0
        for tomoname in self.tomoList:
            tomoNameshort = os.path.basename(tomoname)
            for prec in self.inputCoordinates.get().getPrecedents().iterItems():
                if prec.getFileName() == tomoname:
                    precSet = self._createSetOfTomograms('_'+tomoNameshort)
                    precSet.copyInfo(self.inputCoordinates.get().getPrecedents())
                    precSet.append(prec)
                    tomocoorset = self._createSetOfCoordinates3D(precSet, '_'+tomoNameshort)
                    precSet.write()
            for coorinset in self.inputCoordinates.get().iterItems():
                if coorinset.getVolName() == tomoname:
                    tomocoorset.append(coorinset)
            tomocoorset.write()
            inputCoor = tomocoorset
            minDist = self.distance.get()
            coorlist = []
            for i, coor in enumerate(inputCoor.iterItems()):
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
            vals0list = [i for i, x in enumerate(vals.real) if abs(x) < (1/(np.sqrt(len(inputCoor)))*1e-3)]
            self.listOfSets = [[] for x in range(len(vals0list))]
            for k in range(len(inputCoor)):
                row = []
                for j in vals0list:
                    row.append(vecs[k, j])
                row = [abs(number)for number in row]
                ixMax = np.argmax(row)
                self.listOfSets[ixMax].append(k)
            for ix, coorInd in enumerate(self.listOfSets):
                fi += 1
                outputSet = self._createSetOfCoordinates3D(self.inputCoordinates.get().getPrecedents(), fi)
                outputSet.copyInfo(self.inputCoordinates.get())
                outputSet.setBoxSize(self.inputCoordinates.get().getBoxSize())
                outputSet.setSamplingRate(self.inputCoordinates.get().getSamplingRate())
                for id, coor3D in enumerate(inputCoor.iterItems()):
                    if id in coorInd:
                        outputSet.append(coor3D)
                name = 'output3DCoordinates%s' % str(fi)
                args = {}
                args[name] = outputSet
                self._defineOutputs(**args)
                self._defineSourceRelation(self.inputCoordinates.get(), outputSet)

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
        methods.append("%d connected component identified, with a maximum radial distance of %d voxels."
                       % (len(self._outputs), self.distance.get()))
        return methods

