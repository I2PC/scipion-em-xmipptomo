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


class XmippProtFitEllipsoid(EMProtocol, ProtTomoBase):
    """ This protocol adjust a SetOfCoordinates to an ellipsoid (for example, a vesicle), defining a region of interest
    (ROI) as output."""

    _label = 'fit vesicles'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input coordinates')
        form.addParam('inputCoordinates', PointerParam, label="Coordinates",
                      pointerClass='SetOfCoordinates3D', help='Select the coordinates.')
        form.addParam('inputTomos', PointerParam, label="Tomograms",
                      pointerClass='SetOfTomograms', help='Select tomograms.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fitEllipsoidStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def fitEllipsoidStep(self):
        """ Fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
        and A + B + C = 3 constraint removing one extra parameter. """

        # From matlab function "Fit Ellipsoid"
        for coord in self.inputCoordinates.get():
            x = coord.getX()
            y = coord.getY()
            z = coord.getZ()

            D = np.array([[x*x + y*y - 2*z*z, x*x + z*z - 2*y*y, 2*x*y], [2*x*z, 2*y*z, 2*x], [2*y, 2*z, 1 + 0*x]])  # efg: write as 3d matrix

            # Solve the normal system of equations
            d2 = x*x + y*y + z*z  # The RHS of the llsq problem (y's)
            cD = D.conj().transpose()
            a = cD*D
            b = cD*d2
            u = np.linalg.lstsq(a, b, rcond=None)[0]  # Solution to the normal equations
            u = u.flatten('F')  # Flats matrix u as in Matlab to use same indexes

            # Find the ellipsoid parameters
            # Convert back to the conventional algebraic form
            v = np.empty(10)
            v[0] = u[0] + u[1] - 1  # unique indexing for matrix in matlab
            v[1] = u[0] - 2 * u[1] - 1
            v[2] = u[1] - 2 * u[0] - 1
            v[3:9] = u[2:8]
            v = v.conj().transpose()

            # Form the algebraic form of the ellipsoid
            A = np.array([[v[0], v[3], v[4], v[6]],
                          [v[3], v[1], v[5], v[7]],
                          [v[4], v[5], v[2], v[8]],
                          [v[6], v[7], v[8], v[9]]])

            # Find the center of the ellipsoid
            center = np.array(np.linalg.lstsq(-A[0:2, 0:2], v[6:8], rcond=None))[0]

            # Form the corresponding translation matrix
            T = np.eye(4)
            T[3, 0:2] = center.conj().transpose()
            cT = T.conj().transpose()

            # Translate to the center
            R = T * A * cT

            # Solve the eigenproblem
            [evals, evecs] = np.linalg.eig(R[0:2, 0:2] / -R[3, 3])
            radii = np.sqrt(1/np.diag(np.abs(evals)))
            sgns = np.sign(np.diag(evals))
            radii = radii*sgns

            # Calculate difference of the fitted points from the actual data normalized by the conic radii
            d = np.array([x-center[0], y-center[1], z-center[2]])  # shift data to origin
            d = d * evecs  # Rotate to cardinal axes of the conic
            d = [d[:, 0] / radii[0], d[:, 1] / radii[1], d[:, 2] / radii[2]]  # normalize to the conic radii
            chi2 = np.sum(np.abs(1 - np.sum(d**2 * np.tile(sgns.conj().transpose(), (np.shape(d), 2)))))

            if np.abs(v[-1]) > 1e-6:
                v = -v / v[-1]  # Normalize to the more conventional form with constant term = -1
            else:
                v = -np.sign(v[-1]) * v

    def createOutputStep(self):
        outSet = self._createSetOfMeshes()
        for file in os.listdir(self._getExtraPath()):
            if file.endswith(".txt"):
                data = np.loadtxt(self._getExtraPath(file), delimiter=',')
                groups = np.unique(data[:, 3]).astype(int)
                for group in groups:
                    mesh = Mesh(group=group, path=self._getExtraPath(file))
                    for tomo in self.inputTomos.get().iterItems():
                        if file[:-4] == pwutlis.removeBaseExt(tomo.getFileName()):
                            mesh.setVolume(tomo.clone())
                    outSet.append(mesh)
        outSet.setVolumes(self.inputTomos.get())
        self._defineOutputs(outputMeshes=outSet)
        self._defineSourceRelation(self.inputTomos.get(), outSet)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        if self.inputCoordinates.get().getSize() < 9:
            validateMsgs.append('Set of coordinates must have at least 9 points to fit a unique ellipsoid')
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
