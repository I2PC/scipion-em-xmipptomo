# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Pablo Conesa
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
import enum
import numpy as np
from pwem.convert.transformations import euler_matrix
from pwem.objects.data import Integer
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol.params import IntParam, FloatParam, StringParam
from tomo.protocols import ProtTomoBase
from tomo.objects import TomoAcquisition, Coordinate3D, SetOfCoordinates3D, SetOfTomograms,Tomogram
import tomo.constants as const
from pwem.convert.headers import setMRCSamplingRate

class OutputPhantomTomos(enum.Enum):
    tomograms = SetOfTomograms
    coordinates3D = SetOfCoordinates3D

class XmippProtPhantomTomo(EMProtocol, ProtTomoBase):
    """ Create phantom tomograms with phantom particles and its coordinates with the right Scipion transformation matrix """

    _label = 'phantom tomograms'
    _devStatus = BETA
    _possibleOutputs = OutputPhantomTomos

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('dimensions', StringParam, label='Tomogram dimensions',
                      default='200 200 100',
                      help="Tomogram dimensions: X, Y, Z")
        form.addParam('sampling', FloatParam, label='Sampling rate', default=4)

        form.addParam('nparticles', IntParam, label='Number of particles', default=10,
                      help="How many particles in each tomogram")

        form.addParam('ntomos', IntParam, label='Number of tomograms', default=1,
                      help="How many tomograms")

        # form.addParam('mwfilter', BooleanParam, label='Apply missing wedge?', default=False,
        #               help='Apply a filter to simulate the missing wedge along Y axis.')
        form.addParam('mwangle', IntParam, label='Missing wedge angle', default=60, #condition='mwfilter==True',
                      help='Missing wedge (along y) for data between +- this angle.')

        # Angles
        form.addSection(label='Rotation')
        form.addParam('rotmin', IntParam, label='Min rot angle', default=0,
                      help='Minimum and maximum range for each Euler angle in degrees')
        form.addParam('rotmax', IntParam, label='Max rot angle', default=60)
        form.addParam('tiltmin', IntParam, label='Min tilt angle', default=0)
        form.addParam('tiltmax', IntParam, label='Max tilt angle', default=60)
        form.addParam('psimin', IntParam, label='Min psi angle', default=0)
        form.addParam('psimax', IntParam, label='Max psi angle', default=60)


    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createPhantomsStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def createPhantomsStep(self):

        # Create the set of tomograms
        self.tomoSet = self._createSetOfTomograms(self._getOutputSuffix(SetOfTomograms))
        self.tomoSet.setSamplingRate(self.sampling.get())

        # Create the set of coordinates
        self.coords = self._createSetOfCoordinates3D(self.tomoSet)

        # Create acquisition
        mwangle= self.mwangle.get()
        acq = TomoAcquisition()
        acq.setAngleMax(mwangle)
        acq.setAngleMin(mwangle * -1)

        # Description string to generate the phantom
        desc = self.dimensions.get() + " 0 \n"
        dims = desc.split()
        xT, yT, zT = int(dims[0]), int(dims[1]), int(dims[2])
        minDim = min(xT, yT, zT)
        self.info("Minimum dimension: %s" % minDim)

        radius = max(4, int(minDim * 0.05))
        height = max(10, int(minDim * 0.1))
        self.info("Cone radius, height :%s, %s" % (radius, height))

        # Reduce dimensions form center to avoid particles in the border
        xT = xT/2 - height
        yT = yT/2 - height
        zT = zT/2 - height
        self.info("Valid offset form center, x, y,z: %s, %s, %s" % (xT, yT, zT))
        # Create as many tomograms as indicated
        for i in range(self.ntomos.get()):

            self._addTomogram(i, acq, xT, yT, zT, radius, height, desc)


    def _addTomogram(self, index, acq, xT, yT, zT, radius, height, desc):
        """ Creates creates and adds a phantom tomogram as specified"""

        # Temporary list to keep Coordinated3D
        coords = []

        # Write the description file
        fnDescr = self._getExtraPath("phantom.descr")

        with open(fnDescr, 'w') as fhDescr:

            # For each particle
            for i in range(self.nparticles.get()):

                rot, tilt, psi = self._getRandomAngles()
                x = np.random.randint(-xT, xT)
                y = np.random.randint(-yT, yT)
                z = np.random.randint(-zT, zT)

                # Want to generate a cone --> con + 3 0 0 0 8 30 0 0 0
                desc += "con = 1 %s %s %s %s %s %s %s %s\n" % (x, y, z, radius, height, rot, tilt, psi)

                coords.append(self._createCoordinate(x,y,z,height, rot, tilt, psi))

            # Write the description
            fhDescr.write(desc)

        # Create the phantom based on the description file
        fnVol = self._getExtraPath("phantom_tomo%d.mrc" % index)
        self.runJob("xmipp_phantom_create", " -i %s -o %s" % (fnDescr, fnVol))
        setMRCSamplingRate(fnVol, self.sampling.get())

        # Instantiate the Tomogram object
        tomo = Tomogram()
        tomo.setAcquisition(acq)
        tomo.setLocation(fnVol)

        self.tomoSet.append(tomo)

        # Now that we have the tomogram persisted, we persist the coordinates
        self.coords.setBoxSize(height)
        for coord in coords:
            coord.setVolume(tomo)
            self.coords.append(coord)

    def _getRandomAngles(self):
        """ Returns random rot, tilt, psi in range"""

        rot = np.random.randint(self.rotmin.get(), self.rotmax.get())
        tilt = np.random.randint(self.tiltmin.get(), self.tiltmax.get())
        psi = np.random.randint(self.psimin.get(), self.psimax.get())

        # Get the random values
        return rot, tilt, psi

    def _createCoordinate(self, x, y, z, height, rot, tilt, psi):
        """ Creates a Coordinate3D with the right transfomation matrix"""

        coord = Coordinate3D()
        coord.setX(x, const.SCIPION)
        coord.setY(y, const.SCIPION)
        coord.setZ(z, const.SCIPION)

        coord.setObjComment(
            "Angle Rot, tilt, psi: %d, %d, %d" % (rot, tilt, psi))
        coord.phantom_rot = Integer(rot)
        coord.phantom_tilt = Integer(tilt)
        coord.phantom_psi = Integer(psi)

        # Scipion alignment matrix
        A = euler_matrix(np.deg2rad(psi), np.deg2rad(tilt), np.deg2rad(rot), 'szyz')
        shifts = np.transpose(np.array([0, 0, 0]))

        # Translation matrix
        T = np.eye(4)
        T[:3, 3] = shifts
        M = A@T

        coord.setMatrix(M)

        return coord

    def createOutputStep(self):

        self._defineOutputs(**{self._possibleOutputs.tomograms.name:self.tomoSet})
        self._defineOutputs(**{self._possibleOutputs.coordinates3D.name:self.coords})
        self._defineSourceRelation(self.coords, self.tomoSet)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if self.rotmin.get() >= self.rotmax.get():
            errors.append("rot max must be bigger than rot min")
        if self.tiltmin.get() >= self.tiltmax.get():
            errors.append("tilt max must be bigger than tilt min")
        if self.psimin.get() >= self.psimax.get():
            errors.append("psi max must be bigger than psi min")
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, self._possibleOutputs.tomograms.name):
            summary.append("Output phantom not ready yet.")
        else:
            summary.append("%s phantoms created with random orientations" % self.nsubtomos.get())
            if self.mwfilter.get():
                summary.append("Missing wedge applied between +-%d along Y axis" % self.mwangle.get())
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, self._possibleOutputs.tomograms.name):
            methods.append("Output phantoms not ready yet.")
            return methods
        else:
            methods.append("%s phantoms created with random orientations." % self.nsubtomos.get())
            # if self.mwfilter.get():
            #     methods.append("Missing wedge applied between +-%d along Y axis." % self.mwangle.get())
            return methods
