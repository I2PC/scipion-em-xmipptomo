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

from os import path, listdir
import numpy as np
import random
import math
from pyworkflow.protocol.params import PointerParam
import pyworkflow.utils as pwutlis
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.objects import Mesh, Ellipsoid
from xmipp3.utils import fit_ellipsoid


class XmippProtFitEllipsoid(EMProtocol, ProtTomoBase):
    """ This protocol adjust a SetOfSubtomograms (with coordinates), to a vesicle (ellipsoid), defining regions of
    interest (SetOfMeshes) for each vesicle as output."""

    _label = 'fit vesicles'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubtomos', PointerParam, pointerClass="SetOfSubTomograms",
                      label='Subtomograms', help='Subtomograms in vesicles, they should come from Pyseg.')
        form.addParam('inputTomos', PointerParam, label="Tomograms", pointerClass='SetOfTomograms',
                      help='Select tomograms from which subtomograms come.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fitEllipsoidStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def fitEllipsoidStep(self):
        inputSubtomos = self.inputSubtomos.get()
        vesicleList = []
        self.tomoList = []
        vesicleIdList = []

        # Split input particles by tomograms and by vesicles in each tomogram
        for subtomo in inputSubtomos.iterItems():
            tomoName = subtomo.getVolName()
            subtomoName = subtomo.getFileName()
            vesicleId = self._getVesicleId(subtomoName)
            if tomoName not in self.tomoList:
                self.tomoList.append(tomoName)
            if vesicleId not in vesicleIdList:
                vesicleIdList.append(vesicleId)
                vesicle = self._createSetOfSubTomograms('_' + pwutlis.removeBaseExt(path.basename(tomoName)) +
                                                        '_vesicle_' + vesicleId)
                vesicleList.append(vesicle)
            idx = vesicleIdList.index(vesicleId)
            vesicleList[idx].append(subtomo)

        self.outSet = self._createSetOfMeshes(suffix='_allVesicles')

        # For each vesicle:
        for vesicle in vesicleList:
            x = []
            y = []
            z = []
            for subtomo in vesicle.iterItems():
                coord = subtomo.getCoordinate3D()
                x.append(coord.getX())
                y.append(coord.getY())
                z.append(coord.getZ())

            [center, radii, evecs, v, chi2] = fit_ellipsoid(np.array(x), np.array(y), np.array(z))  # TODO: Review func
            algDesc = '%f*x*x + %f*y*y + %f*z*z + 2*%f*x*y + 2*%f*x*z + 2*%f*y*z + 2*%f*x + 2*%f*y + 2*%f*z + %f = 0' \
                      % (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9])
            adjEllipsoid = Ellipsoid()
            adjEllipsoid.setAlgebraicDesc(algDesc)
            adjEllipsoid.setCenter(str(center))
            adjEllipsoid.setRadii(str(radii))

            print("Algebraic description of the ellipsoid: %s" % algDesc)

            fnVesicle = self._getExtraPath(path.basename(subtomo.getVolName()).split('.')[0] + '_vesicle_' +
                                           self._getVesicleId(subtomo.getFileName()) + '.txt')
            fhVesicle = open(fnVesicle, 'w')
            radii = np.diag(radii)

            print('--------floorx--------', center[0]-radii[0])
            print('--------ceilx--------', center[0]+radii[0])
            print('--------center--------', center)
            print('--------radii-------MUST BE > 0!!!!----', radii)
            print('--------floory--------', center[1]-radii[1])
            print('--------ceily--------', center[1]+radii[1])

            for xcoor in np.linspace((center[0]-radii[0]), (center[0]+radii[0])):  # Default num of samples = 50
                ycoor = random.uniform((center[1]-radii[1]), (center[1]+radii[1]))
                zcoor = np.sqrt((1 - (xcoor**2/radii[0]**2) - (ycoor**2/radii[1]**2)) * radii[2])
                zcoorn = zcoor * (-1)  # write coors for +-z => num samples linspace * 2 = 100 now for each vesicle
                fhVesicle.write('%f,%f,%f,%d\n' % (xcoor, ycoor, zcoor, vesicleList.index(vesicle)))
                fhVesicle.write('%f,%f,%f,%d\n' % (xcoor, ycoor, zcoorn, vesicleList.index(vesicle)))

                print('--------x--------', xcoor)
                print('--------y--------', ycoor)
                print('--------z--------', zcoor)
                print('-------- -z -------', zcoorn)

            fhVesicle.close()
            data = np.loadtxt(fnVesicle, delimiter=',')
            groups = np.unique(data[:, 3]).astype(int)
            for group in groups:
                mesh = Mesh(group=group, path=fnVesicle)  # Group = vesicle in this case
                mesh.setDescription(adjEllipsoid)
                for tomo in self.inputTomos.get().iterItems():
                    if pwutlis.removeBaseExt(fnVesicle[:-4].split('_vesicle')[0]) == \
                            pwutlis.removeBaseExt(tomo.getFileName()):
                        mesh.setVolume(tomo.clone())
                self.outSet.append(mesh)

    def createOutputStep(self):
        self.outSet.setVolumes(self.inputTomos.get())
        self._defineOutputs(outputMeshes=self.outSet)
        self._defineSourceRelation(self.inputTomos.get(), self.outSet)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        if self.inputSubtomos.get().getSize() < 9:
            validateMsgs.append('At least 9 subtomograms are required to fit a unique ellipsoid')
        if not self.inputSubtomos.get().getFirstItem().hasCoordinate3D():
            validateMsgs.append('Subtomograms should have coordinates')
        return validateMsgs

    def _summary(self):
        summary = []
        if not listdir(self._getExtraPath()):
            summary.append("Output vesicles not ready yet.")
        else:
            summary.append("%d subtomograms from %d tomograms\n%d vesicles adjusted" %
                           (self.inputSubtomos.get().getSize(), self.inputTomos.get().getSize(),
                            self.outputMeshes.getSize()))
        return summary

    def _methods(self):
        return ["Fit an ellipsoid and compute a 3D set of points (mesh) for each of the %d vesicles adjusted." %
                self.outputMeshes.getSize()]

    # --------------------------- UTILS functions --------------------------------------------
    def _getVesicleId(self, subtomoName):
        vesicleId = subtomoName.split('tid_')[1]
        vesicleId = vesicleId.split('.')[0]
        return vesicleId
