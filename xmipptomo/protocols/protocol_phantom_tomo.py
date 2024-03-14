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
import math

import numpy as np
from pwem.convert.transformations import euler_matrix
from pwem.objects.data import Integer
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol.params import IntParam, FloatParam, StringParam, BooleanParam
from pyworkflow.utils import removeBaseExt
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

        form.addParam('addNoise', BooleanParam, label="Add noise to the tomogram.", default=True,
                      help="Add noise using xmipp_transform.")

        form.addParam('heterogeneous', BooleanParam, label="2 particles?",
                      default=False, help="Add 2 different particles to allow for 3d classification")

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

        # Hard coded Acquisition
        acq = TomoAcquisition(angleMin=-60, angleMax=60, step=3,
                              accumDose=0, tiltAxisAngle=90,
                              voltage=300, amplitudeContrast= 0.1,
                              sphericalAberration=2.0, magnification=20000,
                              doseInitial=0, dosePerFrame=1
                              )
        self.tomoSet.setAcquisition(acq)

        # Create the set of coordinates
        self.coords = self._createSetOfCoordinates3D(self.tomoSet)
        self.coords.setSamplingRate(self.sampling.get())

        # Create acquisition
        mwangle= self.mwangle.get()
        acq = TomoAcquisition()
        acq.setAngleMax(mwangle)
        acq.setAngleMin(mwangle * -1)

        # Description string to generate the phantom
        desc = self.dimensions.get() + " 1 \n"
        dims = desc.split()
        xT, yT, zT = int(dims[0]), int(dims[1]), int(dims[2])
        minDim = min(xT, yT, zT)
        self.info("Minimum dimension: %s" % minDim)

        boxSize = max(32, int(minDim * 0.1))
        self.info("BoxSize :%s" % boxSize)

        # Create as many tomograms as indicated
        for i in range(self.ntomos.get()):

            self._addTomogram(i, acq, xT, yT, zT, boxSize, desc)


    def _addTomogram(self, index, acq, xT, yT, zT, boxSize, desc):
        """ Creates creates and adds a phantom tomogram as specified"""

        # Temporary list to keep Coordinated3D
        coords = []

        # Write the description file
        fnDescr = self._getExtraPath("phantom.descr")

        # Reduce dimensions from center to avoid particles in the border
        halfHeight = boxSize/2
        xT = xT - halfHeight
        yT = yT - halfHeight
        zT = zT - halfHeight
        self.info("Valid offset from center, x, y,z: %s, %s, %s" % (xT, yT, zT))

        # Create a list with possible positions to warranty no overlapping
        xLen = math.floor(xT/boxSize)
        yLen = math.floor(yT/boxSize)
        zLen = math.floor(zT/boxSize)
        self.info("Possible position matrix (x, y, z): %s, %s, %s." % (xLen, yLen, zLen))

        positions = []

        # Calculate halves of the dimensions
        xHalf = int(xT / 2)
        yHalf = int(yT / 2)
        zHalf = int(zT / 2)

        for z in range(0,zLen):
            for y in range(0, yLen):
                for x in range(0, xLen):

                    # Convert positions to actual coordinates with the origin in the center
                    xPos = int((x * boxSize) + halfHeight) - xHalf
                    yPos = int((y * boxSize) + halfHeight) - yHalf
                    zPos = int((z * boxSize) + halfHeight) - zHalf
                    positions.append((xPos,yPos,zPos))

        with open(fnDescr, 'w') as fhDescr:

            # For each particle
            for i in range(self.nparticles.get()):

                if len(positions) == 0:
                    self.warning("There is no more space to add subtomograms without overlapping")
                    break

                rot, tilt, psi = self._getRandomAngles()

                rng = np.random.default_rng()
                posIndex = rng.integers(len(positions))
                pos = positions[posIndex]
                del positions[posIndex]

                # Heterogeneity, if active every even particle.
                if i % 2 == 1 and self.heterogeneous.get():
                    classId = 2
                    altShape= True
                else:
                    classId = 1
                    altShape = False

                # Want to generate a cone --> con + 3 0 0 0 8 30 0 0 0
                desc += self.getParticleShape(boxSize, i, pos, rot, tilt, psi, alternative=altShape)

                coords.append(self._createCoordinate(pos[0],pos[1],pos[2], rot, tilt, psi, classId))


            # Write the description
            fhDescr.write(desc)

        # Create the phantom based on the description file
        fnVol = self._getExtraPath("phantom_tomo%d.mrc" % index)
        self.runJob("xmipp_phantom_create", " -i %s -o %s" % (fnDescr, fnVol))

        # Add noise
        if self.addNoise.get():
            self.runJob("xmipp_transform_add_noise",  "-i %s -o %s --type gaussian 60 0" % (fnVol, fnVol))

        setMRCSamplingRate(fnVol, self.sampling.get())

        # Instantiate the Tomogram object
        tomo = Tomogram()
        tomo.setAcquisition(acq)
        tomo.setLocation(fnVol)
        tomo.setTsId(removeBaseExt(fnVol))
        self.tomoSet.append(tomo)

        # Now that we have the tomogram persisted, we persist the coordinates
        self.coords.setBoxSize(boxSize)
        for coord in coords:
            coord.setVolume(tomo)
            self.coords.append(coord)

    def getParticleShape(self, boxSize, i, pos, rot, tilt, psi, alternative=False):
        """ Returns the description for the phantom create. See specs here:
        https://web.archive.org/web/20180813105422/http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/FileFormats#Phantom_metadata_file


            param: alternative (False) return an alternative shape to have a second particle
        """

        value = -100 - i
        # Do not use the whole boxSize. Provide some padding
        maxDim = int(boxSize*0.75)

        # Try a more sophisticated one
        # X axis: a bar

        shift = 4 if not alternative else -4

        shape = "cub = %s %s %s %s %s 5 5 %s %s %s\n" % (value, pos[0]-shift, pos[1], pos[2], maxDim, rot, tilt, psi)
        # # Y axis: an ellipsoid
        shape = shape + "ell = %s %s %s %s 5 %s 5 %s %s %s\n" % (value, pos[0], pos[1]+2, pos[2], maxDim/3, rot, tilt, psi)
        # # Z axis: a cone
        shape = shape + "con = %s %s %s %s %s %s %s %s %s\n" % (value, pos[0], pos[1], pos[2], maxDim/4, maxDim, rot, tilt, psi)

        return shape

    def _getRandomAngles(self):
        """ Returns random rot, tilt, psi in range"""

        rng = np.random.default_rng()
        rot = rng.integers(self.rotmin.get(), self.rotmax.get())
        tilt = rng.integers(self.tiltmin.get(), self.tiltmax.get())
        psi = rng.integers(self.psimin.get(), self.psimax.get()

        # Get the random values
        return rot, tilt, psi

    def _createCoordinate(self, x, y, z, rot, tilt, psi, groupId):
        """ Creates a Coordinate3D with the right transfomation matrix"""

        coord = Coordinate3D()
        coord.setX(x, const.SCIPION)
        coord.setY(y, const.SCIPION)
        coord.setZ(z, const.SCIPION)
        coord.setGroupId(groupId)
        coord.setObjComment(
            "Angle Rot, tilt, psi: %d, %d, %d | group: %s" % (rot, tilt, psi, groupId))
        coord.phantom_rot = Integer(rot)
        coord.phantom_tilt = Integer(tilt)
        coord.phantom_psi = Integer(psi)

        # Scipion alignment matrix
        A = euler_matrix(np.deg2rad(psi), np.deg2rad(tilt), np.deg2rad(rot), 'szyz')
        shifts = np.transpose(np.array([0, 0, 0]))

        # Translation matrix
        T = np.eye(4)
        T[:3, 3] = shifts

        # Phantom create assumes that rot, tilt, psi moves the particle to the reference
        # In Scipion convention rot, tilt, psi move the reference to the particle position
        # therefore we need the inverse rotation matrix
        A=np.linalg.inv(A)

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
        if not self._hasOutputs():
            summary.append("Output phantom not ready yet.")
        else:
            # In case not all requested particles fit in the tomogram we calculate the final count
            summary.append("%s phantom tomograms created with %s asymmetric "
                           "anvil-like particles with random orientations" % (self.ntomos ,self._getParticlesPerTomogram()))
        return summary

    def _hasOutputs(self):
        return hasattr(self, self._possibleOutputs.tomograms.name)

    def _getParticlesPerTomogram(self):
        if self._hasOutputs():
            return int(self.coordinates3D.getSize()/self.tomograms.getSize())
        else:
            return self.nparticles.get()

    def _methods(self):
        methods = []

        methods.append("%s synthetic tomograms were created with %s asymmetric anvil-like particles each." % (self.ntomos, self._getParticlesPerTomogram()))
        methods.append("Particle's angles were randomly assigned following the criteria:")
        methods.append("Rot : %s --> %s" % (self.rotmin, self.rotmax))
        methods.append("Tilt: %s --> %s" % (self.tiltmin, self.tiltmax))
        methods.append("Psi : %s --> %s" % (self.psimin, self.psimax))
        methods.append("The corresponding set of 3D coordinates was created with %s elements." % (self.ntomos.get() * self._getParticlesPerTomogram()))

        return methods
