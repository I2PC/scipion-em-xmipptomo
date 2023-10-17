# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.L. Vilas
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
import os
from pwem.convert.transformations import euler_matrix
from pwem.objects.data import Integer
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pwem.emlib import lib
from pwem.emlib.image import ImageHandler
import pwem.emlib.metadata as md
from pyworkflow.protocol.params import IntParam, FloatParam, StringParam, BooleanParam, PointerParam
from pyworkflow.utils import removeBaseExt
from tomo.protocols import ProtTomoBase
from pyworkflow.object import Set
from tomo.objects import TomoAcquisition, Coordinate3D, SetOfCoordinates3D, TiltSeries, SetOfTiltSeries, TiltImage
import tomo.constants as const
from pwem.convert.headers import setMRCSamplingRate

COORDINATES_FN = 'generatedCoordinates.xmd'


class XmippProtPhantomTiltSeries(EMProtocol, ProtTomoBase):
    """ Create phantom tilt series with particles """

    _label = 'simulate tilt series'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', PointerParam,
                      pointerClass="Volume",
                      label='Volume',
                      help="Volume in the sample (tomogram) that will appear as projections in the tilt series")

        form.addParam('numberOfTiltSeries', IntParam, label='Number of tilt series', default=3,
                      help="How many tilt series will be created")

        form.addParam('nparticles', IntParam, label='Particles per tilt series', default=10,
                      help="How many particles in each tomogram")

        lineDim = form.addLine('Tilt series dimensions (px)',
                               help="")
        lineDim.addParam('xdim', IntParam, label='X-size', default=1000)
        lineDim.addParam('ydim', IntParam, label='Y-size', default=1000)
        form.addParam('thickness', IntParam, label='Sample Thickness (px)', default=300)

        form.addParam('sampling', FloatParam, label='Sampling rate', default=1.0)

        # Acquisition
        form.addSection(label='Acquisition')
        line = form.addLine('Tilting Range (degrees)',
                            help="Simulates the set of angles used to tilt the sample and then acquire "
                                 "tilt images in the microscope. The minumum and maximum tilt values limit"
                                 "the tilting range of the sample and the tilt step defines the angular"
                                 "separatation between two consequtive tilt image.")
        line.addParam('minTilt', FloatParam, label='Min Tilt', default=-60.0)
        line.addParam('maxTilt', FloatParam, label='Max Tilt', default=60.0)
        line.addParam('tiltStep', FloatParam, label='step', default=3.0)
        # TODO: Add tilt axis
        # TODO: Add magnification

        form.addParam('voltage', FloatParam, label="Voltage (keV)", default=300.0,
                      help="Voltage of the electron gun.")
        form.addParam('sphericalAberration', FloatParam, label="Spherical Aberration (mm)", default=2.7,
                      help="Spherical aberration of the microscope")
        form.addParam('amplcontrst', FloatParam, label="Amplitude Contrast", default=0.1,
                      help="Amplitude contrast.")
        form.addParam('magnification', FloatParam, label="Magnification", default=50000.0,
                      help="Magnification of the acquired images.")

        form.addSection(label='Alignment')
        form.addParam('addNoise', BooleanParam, label="Add noise to the tilt series.", default=True,
                      help="Add noise using xmipp_transform.")

        form.addParam('sigmaNoise', FloatParam, label="Standard deviation", default=20, condition='addNoise',
                      help="The noise will be Gaussian. This parameter determines the standard deviation of the "
                           "noise distribution.")

        form.addParam('heterogeneous', BooleanParam, label="2 particles?",
                      default=False, help="Add 2 different particles to allow for 3d classification")

        form.addSection(label='Alignment')
        form.addParam('addFiducials', BooleanParam, label="Add fiducials?",
                      default=True, help="Fiducials are useful for aligning the sample. They are gold beads"
                                         "that due to their high atomic number can be identified as high contrasts"
                                         "points in the images.")
        form.addParam('fiducialSize', FloatParam, label='Fiducial diameter (nm)', default=10.0,
                      help='Defines the fiducial diameter')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        for i in range(0, self.numberOfTiltSeries.get()):
            tsId = 'TS_%i' % i
            self._insertFunctionStep(self.createTiltSeriesStep, tsId)
            self._insertFunctionStep(self.createOutputStep, tsId)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions --------------------------------------------
    def createTiltSeriesStep(self, tsId):

        pathTs = self._getExtraPath(tsId)

        if not os.path.isdir(pathTs):
            os.mkdir(pathTs)

        tsName = os.path.join(pathTs, tsId + '.mrcs')

        coordinatesFn = self.createCoordinateFile()

        params = ' --coordinates %s ' % coordinatesFn
        params += ' --vol %s ' % self.inputVolume.get().getFileName()
        params += ' --xdim %i ' % self.xdim.get()
        params += ' --ydim %i ' % self.ydim.get()
        params += ' --minTilt %i ' % self.minTilt.get()
        params += ' --maxTilt %i ' % self.maxTilt.get()
        params += ' --tiltStep %i ' % self.tiltStep.get()
        params += ' --thickness %i ' % self.thickness.get()
        params += ' --sampling %f ' % self.inputVolume.get().getSamplingRate()
        params += ' --sigmaNoise %f ' % self.sigmaNoise.get()

        if self.addFiducials.get():
            params += ' --fiducialSize %f ' % self.fiducialSize.get()
        params += ' -o %s ' % tsName

        self.runJob('xmipp_tomo_simulate_tilt_series', params)

    def closeOutputSetsStep(self):
        self.outputSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.outputSetOfTiltSeries.write()
        self._store()

    def createCoordinateFile(self):
        boxsize = self.inputVolume.get().getXDim()
        coordinatesFn = self._getExtraPath(COORDINATES_FN)

        import random

        listOfCoords = []
        # This is a reasonable number of attempts to set particles
        attemps = 1000
        attemptnumber = 0
        while len(listOfCoords) < self.nparticles.get() and attemptnumber < attemps:
            xpos = random.randint(0, self.xdim.get())
            ypos = random.randint(0, self.ydim.get())
            zpos = random.randint(0, self.thickness.get())

            pos = [xpos, ypos, zpos]
            if self.checkdistance(pos, boxsize, listOfCoords):
                listOfCoords.append(pos)
            else:
                attemptnumber += 1

        mdCoords = lib.MetaData()
        for coord in listOfCoords:
            nRow = md.Row()

            nRow.setValue(lib.MDL_XCOOR, coord[0])
            nRow.setValue(lib.MDL_YCOOR, coord[1])
            nRow.setValue(lib.MDL_ZCOOR, coord[2])
            nRow.addToMd(mdCoords)

        mdCoords.write(coordinatesFn)

        return coordinatesFn

    def checkdistance(self, pos, boxsize, listOfCoords):

        addCoor = True
        for coord in listOfCoords:
            if not np.linalg.norm(np.array(pos) - np.array(coord)) > (np.sqrt(2) * boxsize):
                addCoor = False
                break
        return addCoor

    def createOutputStep(self, tsId):
        """Generate output filtered tilt series"""

        sampling = self.inputVolume.get().getSamplingRate()
        fnTs = os.path.join(self._getExtraPath(tsId), tsId+'.mrcs')
        outputSetOfTiltSeries = self.getOutputSetOfTiltSeries(sampling, fnTs)
        acquisitionParams = TomoAcquisition(angleMin=self.minTilt.get(), angleMax=self.maxTilt.get(), step=self.tiltStep.get(),
                 accumDose=None, tiltAxisAngle=0.0)
        newTs = TiltSeries(tsId=tsId)
        outputSetOfTiltSeries.append(newTs)

        acquisitionParams.setVoltage(self.voltage.get())
        acquisitionParams.setMagnification(self.magnification.get())
        acquisitionParams.setAmplitudeContrast(self.amplcontrst.get())
        acquisitionParams.setSphericalAberration(self.sphericalAberration.get())

        outputSetOfTiltSeries.setAcquisition(acquisitionParams)

        newTs.setSamplingRate(sampling)
        newTs.setAcquisition(acquisitionParams)

        fnts = os.path.join(self._getExtraPath(tsId), "%s.xmd" % tsId)

        if os.path.isfile(fnts):
            for row in md.iterRows(fnts):
                fn = row.getValue(lib.MDL_IMAGE)
                tilt = row.getValue(lib.MDL_ANGLE_TILT)

                newTi = TiltImage()
                index = int(fn.split('@')[0])
                fnTi = fn.split('@')[1]
                newTi.setObjId(index)
                newTi.setTiltAngle(tilt)
                newTi.setLocation(index + 1, fnTi)
                newTi.setSamplingRate(sampling)

                newTs.append(newTi)

        newTs.write(properties=False)

        outputSetOfTiltSeries.update(newTs)
        outputSetOfTiltSeries.updateDim()
        outputSetOfTiltSeries.write()

        self._store()

    def getOutputSetOfTiltSeries(self, sampling, fnTs):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            ih = ImageHandler()
            x, y, z, n = ih.getDimensions(fnTs)
            print('n=', z)
            dims = (x, y, z, n)
            #outputSetOfTiltSeries.setDim(dims)
            outputSetOfTiltSeries.setSamplingRate(sampling)
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputVolume, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries

    def closeOutputSetsStep(self):
        self.outputSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.outputSetOfTiltSeries.write()
        self._store()

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        if not self._hasOutputs():
            summary.append("Output phantom not ready yet.")
        else:
            # In case not all requested particles fit in the tomogram we calculate the final count
            summary.append("%s phantom tomograms created with %s asymmetric "
                           "anvil-like particles with random orientations" % (
                           self.ntomos, self._getParticlesPerTomogram()))
        return summary

    def _hasOutputs(self):
        return hasattr(self, self._possibleOutputs.tomograms.name)

    def _getParticlesPerTomogram(self):
        if self._hasOutputs():
            return int(self.coordinates3D.getSize() / self.tomograms.getSize())
        else:
            return self.nparticles.get()

    def _methods(self):
        methods = []

        methods.append("%s synthetic tomograms were created with %s asymmetric anvil-like particles each." % (
        self.ntomos, self._getParticlesPerTomogram()))
        methods.append("Particle's angles were randomly assigned following the criteria:")
        methods.append("Rot : %s --> %s" % (self.rotmin, self.rotmax))
        methods.append("Tilt: %s --> %s" % (self.tiltmin, self.tiltmax))
        methods.append("Psi : %s --> %s" % (self.psimin, self.psimax))
        methods.append("The corresponding set of 3D coordinates was created with %s elements." % (
                    self.ntomos.get() * self._getParticlesPerTomogram()))

        return methods
