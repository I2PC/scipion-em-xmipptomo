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

from os.path import basename, dirname, join
import numpy as np
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
import pyworkflow.utils as pwutlis
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.objects import MeshPoint, Ellipsoid, SubTomogram
from tomo.utils import fit_ellipsoid, generatePointCloud
import tomo.constants as const


class XmippProtFitEllipsoid(EMProtocol, ProtTomoBase):
    """ This protocol adjust a SetOfSubtomograms with coordinates assigned or a SetOfCoordinates3D, to a vesicle
    (ellipsoid), defining regions of interest (SetOfMeshes) for each vesicle as output."""

    _label = 'fit vesicles'
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input', PointerParam, pointerClass="SetOfSubTomograms, SetOfCoordinates3D",
                      label='Subtomograms/Coordinates3D',
                      help='Subtomograms or coordinates3D picked in vesicles. If there are more than one vesicle per '
                           'tomogram, input subtomograms or coordinates should have assigned groupId.')
        form.addParam('inputTomos', PointerParam, label="Tomograms", pointerClass='SetOfTomograms',
                      help='Select tomograms from which subtomograms come.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fitEllipsoidStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def fitEllipsoidStep(self):
        input = self.input.get()
        inputTomos = self.inputTomos.get()
        self.outSet = self._createSetOfMeshes(inputTomos)
        totalMeshes = 0
        for tomo in inputTomos.iterItems(iterate=False):
            vesicleList = []
            vesicleIdList = []
            tomoName = basename(tomo.getFileName())
            tomoDim = [float(d) for d in tomo.getDim()]
            firstItem = input.getFirstItem()
            if self._getInputisSubtomo(firstItem):
                dirName = dirname(firstItem.getVolName())
                for item in input.iterItems(where="_volName='%s'" % join(dirName, tomoName)):
                    vesicleId = self._getVesicleId(item)
                    if vesicleId not in vesicleIdList:
                        vesicleIdList.append(vesicleId)
                        vesicle = self._createSetOfSubTomograms('_' + pwutlis.removeBaseExt(tomoName) +
                                                                '_vesicle_' + str(vesicleId))
                        vesicleList.append(vesicle)
                    idx = vesicleIdList.index(vesicleId)
                    vesicleList[idx].append(item)

            else:
                for item in input.iterCoordinates(volume=tomo):
                    vesicleId = self._getVesicleId(item)
                    if vesicleId not in vesicleIdList:
                        vesicleIdList.append(vesicleId)
                        vesicle = self._createSetOfCoordinates3D(inputTomos, '_' + pwutlis.removeBaseExt(tomoName)
                                                                 + '_vesicle_' + str(vesicleId))
                        vesicleList.append(vesicle)
                    idx = vesicleIdList.index(vesicleId)
                    vesicleList[idx].append(item)

            totalMeshes += len(vesicleList)

            for vesicle in vesicleList:
                x = []
                y = []
                z = []

                if self._getInputisSubtomo(input.getFirstItem()):
                    for item in vesicle.iterItems():
                        coord = self._getCoor(item)
                        x.append(float(coord.getX(const.BOTTOM_LEFT_CORNER)) / tomoDim[0])
                        y.append(float(coord.getY(const.BOTTOM_LEFT_CORNER)) / tomoDim[1])
                        z.append(float(coord.getZ(const.BOTTOM_LEFT_CORNER)) / tomoDim[2])
                else:
                    for item in vesicle.iterCoordinates(volume=tomo):
                        coord = self._getCoor(item)
                        x.append(float(coord.getX(const.BOTTOM_LEFT_CORNER)) / tomoDim[0])
                        y.append(float(coord.getY(const.BOTTOM_LEFT_CORNER)) / tomoDim[1])
                        z.append(float(coord.getZ(const.BOTTOM_LEFT_CORNER)) / tomoDim[2])

                [center, radii, v, _, chi2] = fit_ellipsoid(np.array(x), np.array(y), np.array(z))
                algDesc = '%f*x*x + %f*y*y + %f*z*z + 2*%f*x*y + 2*%f*x*z + 2*%f*y*z + 2*%f*x + 2*%f*y + 2*%f*z + %f ' \
                          '= 0' % (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9])
                adjEllipsoid = Ellipsoid()
                adjEllipsoid.setAlgebraicDesc(algDesc)
                adjEllipsoid.setCenter(str(center))
                adjEllipsoid.setRadii(str(radii))
                print(algDesc)
                print('Chi2: ', chi2)
                pointCloud = generatePointCloud(v, tomo.getDim())
                if not pointCloud:
                    raise Exception("It does not seem like any output is produced!")

                for point in pointCloud:
                    meshPoint = MeshPoint()
                    meshPoint.setVolume(tomo.clone())
                    meshPoint.setX(point[0], const.BOTTOM_LEFT_CORNER)
                    meshPoint.setY(point[1], const.BOTTOM_LEFT_CORNER)
                    meshPoint.setZ(point[2], const.BOTTOM_LEFT_CORNER)
                    meshPoint.setGroupId(self._getVesicleId(item))
                    meshPoint.setDescription(adjEllipsoid)
                    meshPoint.setVolumeName(basename(tomo.getFileName()))
                    self.outSet.append(meshPoint)
        self.outSet.setPrecedents(inputTomos)
        self.outSet.setNumberOfMeshes(totalMeshes)

    def createOutputStep(self):
        self._defineOutputs(outputMeshes=self.outSet)
        self._defineSourceRelation(self.inputTomos.get(), self.outSet)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        input = self.input.get()
        if input.getSize() < 9:
            validateMsgs.append('At least 9 subtomograms/coordinates are required to fit a unique ellipsoid')
        if self._getInputisSubtomo(self.input.get().getFirstItem()):
            if not input.getFirstItem().hasCoordinate3D():
                validateMsgs.append('Subtomograms should have coordinates')
        return validateMsgs

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Output vesicles not ready yet.")
        else:
            summary.append("%d subtomograms/coordinates3D from %d tomograms\n%d vesicles adjusted" %
                           (self.input.get().getSize(), self.inputTomos.get().getSize(),
                            self.outputMeshes.getNumberOfMeshes()))
        return summary

    def _methods(self):
        if not self.isFinished():
            return ["Protocol not finished yet"]
        else:
            return ["Fit an ellipsoid and compute a 3D set of points (mesh) for each of the %d vesicles adjusted." %
                self.outputMeshes.getNumberOfMeshes()]

    # --------------------------- UTILS functions --------------------------------------------
    def _getInputisSubtomo(self, item):
        if isinstance(item, SubTomogram):
            return True
        else:
            return False

    def _getCoor(self, item):
        if self._getInputisSubtomo(item):
            coor = item.getCoordinate3D()
        else:
            coor = item
        return coor

    def _getVesicleId(self, item):
        coor = self._getCoor(item)
        if coor.hasGroupId():
            vesicleId = coor.getGroupId()
        else:
            vesicleId = '1'  # For now it works with several vesicles in the same tomo just for subtomos with groupId
        return vesicleId

    def _evaluateQuadric(self, v, x, y, z):
        return v[0]*x*x + v[1]*y*y + v[2]*z*z + 2*v[3]*x*y + 2*v[4]*x*z + 2*v[5]*y*z + 2*v[6]*x + 2*v[7]*y + 2*v[8]*z +\
               v[9]
