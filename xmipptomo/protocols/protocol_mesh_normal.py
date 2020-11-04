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

from os import path
import numpy as np
from pyworkflow.protocol.params import PointerParam, FloatParam, LEVEL_ADVANCED, BooleanParam, IntParam, EnumParam
import pyworkflow.utils as pwutlis
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.utils import delaunayTriangulation, computeNormals, normalFromMatrix


class XmippProtFilterbyNormal(EMProtocol, ProtTomoBase):
    """ This protocol takes surfaces or ROIs (SetOfMeshes) and a SetOfSubtomograms and filters them by different
    criteria related with the normal direction."""

    _label = 'filter subtomos by normal'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubtomos', PointerParam, pointerClass="SetOfSubTomograms",
                      label='Subtomograms', help='SetOfSubtomograms to filter.')
        form.addParam('inputMeshes', PointerParam, label="Vesicles", pointerClass='SetOfMeshes',
                      help='Select the vesicles in which the subtomograms are.')
        form.addParam('tol', FloatParam, default=0.1, label='Tolerance', expertLevel=LEVEL_ADVANCED,
                      help='Tolerance for normal comparisson')
        # form.addParam('topBottom', BooleanParam, default=True,
        #               label='Remove particles in the top and bottom of the vesicle',
        #               help='Remove the particles that have been picked from the top and bottom parts because they '
        #                    'had a different view.')
        # form.addParam('tilt', IntParam, default=60, label='Maximun allowed tilt',
        #               help='Remove the particles that have been picked from the top and bottom with a tilt bigger
        #               than the one specified in here.')
        # form.addParam('mwDir', BooleanParam, default=True,
        #               label='Remove particles in the missing wedge direction',
        #               help='Remove the particles that are in the missing wedge direction because they are highly '
        #                    'affected by the missing wedge.')
        # form.addParam('mwDir', EnumParam, default=True, label='Missing wedge direction',
        #               help='Missing wedge direction of the tomograms.')
        # form.addParam('normalDir', BooleanParam, default=True,
        #               label='Remove particles if they are not perpendicular to the membrane of the vesicle',
        #               help='Remove the particles that have a normal direction not equal to the normal direction of '
        #                    'the vesicle in the coordinate of the particle.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeNormalStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def computeNormalStep(self):
        self.outSet = self._createSetOfSubTomograms()
        self.outSet.copyInfo(self.inputSubtomos.get())
        tol = self.tol.get()
        for subtomo in self.inputSubtomos.get():
            for mesh in self.inputMeshes.get().iterItems():
                pathV = pwutlis.removeBaseExt(path.basename(mesh.getPath())).split('_vesicle_')
                if pwutlis.removeBaseExt(path.basename(subtomo.getVolName())) == pathV[0]:
                    if self._getVesicleId(subtomo) == pathV[1]:
                        normalsList = self._getNormalVesicleList(mesh)
                        normSubtomo, normVesicle = self._getNormalVesicle(normalsList, subtomo)
                        if abs(normSubtomo[0]-normVesicle[0]) < tol and abs(normSubtomo[1]-normVesicle[1]) < tol \
                                and abs(normSubtomo[2]-normVesicle[2]) < tol:
                            self.outSet.append(subtomo)

    def createOutputStep(self):
        self._defineOutputs(outputset=self.outSet)
        self._defineSourceRelation(self.inputSubtomos.get(), self.outSet)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        if not self.inputMeshes.get().getFirstItem().hasDescription():
            validateMsgs.append('Vesicles not adjusted, please use protocol "fit vesicles" previously')
        return validateMsgs

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Output vesicles not ready yet.")
        else:
            # if self.topBottom:
            #     topBottom = 'Yes'
            # else:
            #     topBottom = 'No'
            # if self.mwDir:
            #     mwDir = 'Yes'
            # else:
            #     mwDir = 'No'
            # if self.normalDir:
            #     normalDir = 'Yes'
            # else:
            #     normalDir = 'No'
            # summary.append("Filter criteria:\nTop/Bottom %s\nMW direction %s\nNormal direction %s" %
            #                (topBottom, mwDir, normalDir))
            summary.append("Remove particles in normal direction")
            summary.append("Tolerance: %0.2f" % self.tol.get())
        return summary

    def _methods(self):
        methods = []
        if not self.isFinished():
            methods.append("Output vesicles not ready yet.")
        else:
            methods.append("%d subtomograms filtered from %d input subtomograms." %
                           ((len(self.inputSubtomos.get()) - len(self.outputset.get())), len(self.inputSubtomos.get())))
            # if self.topBottom:
            #     methods.append("Particles in the top and bottom parts of the vesicles have been removed.")
            # if self.mwDir:
            #     methods.append("Particles in the missing wedge direction have been removed.")
            if self.normalDir:
                methods.append("Particles that are not perpendicular to the membrane have been removed.")
        return methods

    # --------------------------- UTILS functions --------------------------------------------
    def _getVesicleId(self, subtomo):
        if subtomo.getCoordinate3D().hasAttribute('_vesicleId'):
            vesicleId = subtomo.getFileName().split('tid_')[1]
            vesicleId = vesicleId.split('.')[0]
        else:  # For now it works with several vesicles in the same tomo just for Pyseg subtomos
            vesicleId = 1
        return vesicleId

    def _getNormalVesicleList(self, mesh):
        triangulation = delaunayTriangulation(mesh.getMesh())
        normalsList = computeNormals(triangulation, associateCoords=True)
        return normalsList

    def _getNormalVesicle(self, normalsList, subtomo):
        normSubtomo = normalFromMatrix(subtomo.getTransform().getMatrix())
        coord = subtomo.getCoordinate3D()
        coors = np.asarray([coord.getX(), coord.getY(), coord.getZ()])
        points, normals = zip(*normalsList)
        points = np.asarray(points)
        idx = np.argmin(np.sum((points - coors) ** 2, axis=1))
        print('---coord---', coors)
        print('---point---', points[idx])
        print('---normalS---', normSubtomo)
        print('---normalV---', normals[idx])
        return normSubtomo, normals[idx]
