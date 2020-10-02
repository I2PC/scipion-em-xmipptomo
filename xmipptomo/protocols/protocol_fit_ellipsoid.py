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
from pyworkflow.protocol.params import PointerParam
import pyworkflow.utils as pwutlis
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
from tomo.objects import Mesh
from xmipp3.utils import fit_ellipsoid


class XmippProtFitEllipsoid(EMProtocol, ProtTomoBase):
    """ This protocol adjust a SetOfCoordinates to an ellipsoid (for example, a vesicle), defining a region of interest
    (ROI) as output."""

    _label = 'fit vesicles'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input coordinates')
        form.addParam('inputSubtomos', PointerParam, pointerClass="SetOfSubTomograms",
                      label='Subtomograms', help='Subtomograms in vesicles')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fitEllipsoidStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def fitEllipsoidStep(self):
        inputSubtomos = self.inputSubtomos.get()
        coorSetList = []
        tomoList = []

        # Split input particles by tomograms
        # TEST IT WITH PARTICLES FROM DIFF TOMOS
        for subtomo in inputSubtomos.iterItems():
            tomoName = subtomo.getVolName()
            if tomoName not in tomoList:
                tomoList.append(tomoName)
                tomocoorset = self._createSetOfSubTomograms('_' + os.path.basename(tomoName))
                coorSetList.append(tomocoorset)
            idx = tomoList.index(tomoName)
            coorSetList[idx].append(subtomo)

        # Split input particles in each tomo by vesicles
        for coorSet in coorSetList:
            x = []
            y = []
            z = []
            for subtomo in coorSet.iterItems():
                coord = subtomo.getCoordinate3D()
                x.append(coord.getX())
                y.append(coord.getY())
                z.append(coord.getZ())

            [center, radii, evecs, v, chi2] = fit_ellipsoid(np.array(x), np.array(y), np.array(z))
            algDesc = '%f*x*x + %f*y*y + %f*z*z + 2*%f*x*y + 2*%f*x*z + 2*%f*y*z + 2*%f*x + 2*%f*y + 2*%f*z + %f = 0" ' \
                      % (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9])
            print("Center: ", center)
            print("Radii: ", radii)
            print("Evecs: ", evecs)
            print("Chi2: ", chi2)
            # print("----------v-----", v)
            print("Algebraic description of the ellipsoid: %s" % algDesc)

    def createOutputStep(self):
        pass
        # outSet = self._createSetOfMeshes()
        # for file in os.listdir(self._getExtraPath()):
        #     if file.endswith(".txt"):
        #         data = np.loadtxt(self._getExtraPath(file), delimiter=',')
        #         groups = np.unique(data[:, 3]).astype(int)
        #         for group in groups:
        #             mesh = Mesh(group=group, path=self._getExtraPath(file))
        #             for tomo in self.inputTomos.get().iterItems():
        #                 if file[:-4] == pwutlis.removeBaseExt(tomo.getFileName()):
        #                     mesh.setVolume(tomo.clone())
        #             outSet.append(mesh)
        # outSet.setVolumes(self.inputTomos.get())
        # self._defineOutputs(outputMeshes=outSet)
        # self._defineSourceRelation(self.inputTomos.get(), outSet)

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
        if not os.listdir(self._getExtraPath()):
            summary.append("Output vesicles not ready yet.")
        else:
            summary.append("")
        return summary

    def _methods(self):
        return ["Fit an ellipsoid (vesicle) into a 3D set of points."]
