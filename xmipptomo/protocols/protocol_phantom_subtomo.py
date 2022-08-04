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
import enum
import numpy as np
from pwem.convert.transformations import euler_matrix
from pwem.objects.data import Transform, Integer
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol.params import IntParam, FloatParam, EnumParam, PointerParam, TextParam, BooleanParam
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms, SubTomogram, TomoAcquisition, Coordinate3D, SetOfCoordinates3D, \
    SetOfTomograms
import tomo.constants as const
from pwem.convert.headers import setMRCSamplingRate

FN_PARAMS = 'projection.params'
FN_PHANTOM_DESCR = 'phantom.descr'
FN_PHANTOM = 'phantom_'
MRC_EXT = '.mrc'


class OutputPhantomSubtomos(enum.Enum):
    outputSubtomograms = SetOfSubTomograms
    outputCoord = SetOfCoordinates3D


class XmippProtPhantomSubtomo(EMProtocol, ProtTomoBase):
    """ Create subtomogram phantoms """

    _label = 'phantom create subtomo'
    _devStatus = BETA
    _possibleOutputs = OutputPhantomSubtomos

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('option', EnumParam,
                      choices=['Import volume', 'creating a phantom'], default=0,
                      display=EnumParam.DISPLAY_HLIST, label=' ',
                      help="Import a volume or create 'base' phantom manually")
        form.addParam('inputVolume', PointerParam, pointerClass="Volume", label='Input volume',
                      condition="option==0", help="Volume used as 'base' phantom", allowsNull=True)
        form.addParam('create', TextParam, label='Phantom description', condition="option==1",
                      default='40 40 40 0\ncyl + 1 0 0 0 15 15 2 0 0 0\nsph + 1 0 0 5 2\ncyl + 1 0 0 -5 2 2 10 0 90 0\n'
                              'sph + 1 0 -5 5 2',
                      help="create a phantom description: x y z backgroundValue geometry(cyl, sph...) +(superimpose) "
                           "density value origin radius height rot tilt psi. More info at https://web.archive.org/web/20180813105422/http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/FileFormats#Phantom_metadata_file")

        form.addParam('simulateTiltSeries', BooleanParam, label='Simulating tilt series', default=False)
        #lineAngSamp = form.addLine('Tilt Sampling (degrees)', condition='simulateTiltSeries',
        #                           help='The tilt series is acquires from -angle to + angle, in steps of d degrees. '
        #                                'Angular Sampling in degrees. The tilt series was acquired by taking images in intervals of X angles. X is the angualr sampling')
        #lineAngSamp.addParam('mintilt', IntParam, label='min', default=-60)
        #lineAngSamp.addParam('maxtilt', IntParam, label='max', default=60)
        #lineAngSamp.addParam('angularSampling', IntParam, label='step', default=3)
        form.addParam('sampling', FloatParam, label='Sampling rate (A/px)', default=1)
        form.addParam('nsubtomos', IntParam, label='Number of subtomograms', default=50,
                      help="How many phantom subtomograms")

        form.addParam('mwfilter', BooleanParam, label='Apply missing wedge?', default=False,
                      condition='not simulateTiltSeries',
                      help='Apply a filter to simulate the missing wedge along Y axis.')
        form.addParam('mwangle', IntParam, label='Missing wedge angle', default=60,
                      condition='mwfilter==True and (not simulateTiltSeries)',
                      help='Missing wedge (along y) for data between +- this angle.')
        # Angles
        form.addSection(label='Rotation')
        form.addParam('rotate', BooleanParam, label='Apply rotation?', default=False,
                      help='Apply a random rotation to the generated subtomograms. The subtomograms will present a random'
                           'orientation.')
        form.addParam('uniformAngularDistribution', BooleanParam, label='Randomly distributed?', default=True, condition='rotate',
                      help='Apply a random rotation to the generated subtomograms. The subtomograms will present a random'
                           'orientation.')

        form.addParam('stdError', BooleanParam, label='Introduce random error with std', default=False, condition='rotate and uniformAngularDistribution',
                      help='It introduces angular assignment errors with standard deviation given by the introduced value'
                           'It is assumed that the errors are Gaussian.')
        form.addParam('sigma', FloatParam, label='std', default=0,
                      condition='rotate and uniformAngularDistribution and stdError',
                      help='It introduces angular assignment errors with standard deviation given by the introduced value'
                           'It is assumed that the errors are Gaussian.')

        lineRot = form.addLine('Rot Angle (degrees)',
                               help="This is the rot angle (in-plane rotation). Minimum and maximum range for each Euler angle in degrees",
                               condition='rotate and (not uniformAngularDistribution)')
        lineRot.addParam('rotmin', IntParam, label='Min', default=0, condition='rotate and (not uniformAngularDistribution)')
        lineRot.addParam('rotmax', IntParam, label='Max', default=60, condition='rotate and (not uniformAngularDistribution)')
        lineTilt = form.addLine('Tilt Angle (degrees)',
                                help="This is the Tilt angle (latitude). Minimum and maximum range for each Euler angle in degrees",
                                condition='rotate and (not uniformAngularDistribution)')
        lineTilt.addParam('tiltmin', IntParam, label='Min', default=0, condition='rotate and (not uniformAngularDistribution)')
        lineTilt.addParam('tiltmax', IntParam, label='Max', default=60, condition='rotate and (not uniformAngularDistribution)')
        linePsi = form.addLine('Psi Angle (degrees)',
                               help="This is the Psi angle. Minimum and maximum range for each Euler angle in degrees",
                               condition='rotate and (not uniformAngularDistribution)')
        linePsi.addParam('psimin', IntParam, label='Min', default=0, condition='rotate and (not uniformAngularDistribution)')
        linePsi.addParam('psimax', IntParam, label='Max', default=60, condition='rotate and (not uniformAngularDistribution)')


        # Shifts
        form.addSection(label='Shifts')
        form.addParam('applyShift', BooleanParam, label='Apply random shift?', default=False,
                      help='Apply a random shit to the generated subtomograms. The subtomograms will present a random'
                           'displacement from the center of the box.')
        lineX = form.addLine('Shift along X-axis:',
                             help='It considers a random shift between the minumum and the maximum values, '
                                  'following a uniform distribution.', condition='applyShift')
        lineX.addParam('xmin', IntParam, label='Min (px)', default=0, condition='applyShift')
        lineX.addParam('xmax', IntParam, label='Max (px)', default=5, condition='applyShift')
        lineY = form.addLine('Shift along Y-axis:',
                             help='It considers a random shift between the minumum and the maximum values, '
                                  'following a uniform distribution.', condition='applyShift')
        lineY.addParam('ymin', IntParam, label='Min (px)', default=0, condition='applyShift')
        lineY.addParam('ymax', IntParam, label='Max (px)', default=5, condition='applyShift')
        lineZ = form.addLine('Shift along Z-axis:',
                             help='It considers a random shift between the minumum and the maximum values, '
                                  'following a uniform distribution.', condition='applyShift')
        lineZ.addParam('zmin', IntParam, label='Min (px)', default=0, condition='applyShift')
        lineZ.addParam('zmax', IntParam, label='Max (px)', default=5, condition='applyShift')

        form.addSection(label='Coordinates')
        form.addParam('coords', BooleanParam, label='Assign random coordinates?', default=False,
                      help='Create random x, y, z coordinates for each subtomogram.')
        form.addParam('tomos', PointerParam, pointerClass="SetOfTomograms", label='Tomograms',
                      condition="coords==True", help="Tomograms to get dimension for random creation of coordinates")

        form.addSection(label='Noise')
        form.addParam('addNoise', BooleanParam, label='Add gaussian noise to subtomograms?', default=False,
                      help='Select true to generate noisy subtomograms, and False to obtain clean subtomograms.')
        form.addParam('differentStatistics', BooleanParam, condition='addNoise',
                      label='Add variable gaussian noise to subtomograms?', default=True,
                      help='(False) Each subtomogram will follow a different noise distribution. All distribution are Gaussian, the noise of two'
                           'different subtomograms will follow two different Gaussian distributions. It means, given two subtomograms'
                           'A and B, the noise of the subtomogram A will follow a Gaussian distribution with mean mu_A, and std s_A,'
                           'in contrast, the subtomograms B will follow also a Gaussian distribution but with different mean, mu_B, '
                           'and different standard deviation s_B. If True, the noise of both subtomograms will follow the same Gaussian'
                           'distribution')

        lineStatsStd = form.addLine('range of standard deviation between',
                                    help='The guassian standard deviation will be a random number between the provided range (min,max)',
                                    condition='addNoise and differentStatistics')
        lineStatsStd.addParam('minstd', IntParam, label='Min', default=0, condition='differentStatistics')
        lineStatsStd.addParam('maxstd', IntParam, label='Max', default=50, condition='differentStatistics')
        lineStat = form.addLine('Gaussian Noise mean and std',
                                help='Defines the statistics of noise by providing the mean and standard deviation'
                                     'of the Gaussian distribution.',
                                condition='addNoise and not differentStatistics')
        lineStat.addParam('meanNoise', IntParam, label='mean', default=0, condition='not differentStatistics')
        lineStat.addParam('stdNoise', IntParam, label='std', default=40, condition='not differentStatistics')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):

        #NOTE: This protocol was discussed with the ScipionTeam about if the subtomograms have or do not have tomogramId.
        # The agreement was that the tomogramId should not appear if subtomograms are imported or phantoms created
        self._insertFunctionStep(self.createSubtomogramsStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def createSubtomogramsStep(self):
        self.createPhantomSubtomograms()

    def createPhantomSubtomograms(self):
        if self.mwfilter.get():
            mwangle = self.mwangle.get()
        else:
            mwangle = 90

        fnInVol = None
        dim = None
        if self.option == 0:
            inputVol = self.inputVolume.get()
            fnInVol = inputVol.getFileName()
            dim = inputVol.getDim()
        if self.option == 1:
            dim, fnVol = self.createGeometricalPhantom()

        self.definingOrientationsAndRegisteringInformation(dim, mwangle, fnInVol)


    def definingOrientationsAndRegisteringInformation(self, dim, mwangle, fnVol):
        self.createOutputSet(dim)
        tomo = None
        coordsBool = self.generateCoordinates()
        if coordsBool:
            tomos = self.tomos.get()
            tomo = tomos.getFirstItem()
            self.coordsSet = self._createSetOfCoordinates3D(tomos)
            self._store(self.coordsSet)

        # Create acquisition
        acq = TomoAcquisition()
        acq.setAngleMax(mwangle)
        acq.setAngleMin(mwangle * -1)

        for i in range(int(self.nsubtomos.get())):
            fnPhantomi = self._getExtraPath(FN_PHANTOM + str(int(i+1)) + MRC_EXT)

            self.addNoiseToPhantom(fnVol, fnPhantomi)

            rot, tilt, psi, shiftX, shiftY, shiftZ = self.applyRandomOrientation(fnPhantomi, fnPhantomi)

            if self.mwfilter.get():
                self.applyMissingWedge(mwangle, fnPhantomi, fnPhantomi)

            # Add the subtomogram and the coordinate if applies
            self._addSubtomogram(tomo, acq, fnPhantomi, rot, tilt, psi, shiftX, shiftY, shiftZ)

        if coordsBool:
            self.outputSet.setCoordinates3D(self.coordsSet)

    def createGeometricalPhantom(self):
        fnVol = self._getExtraPath(FN_PHANTOM+MRC_EXT)
        desc = self.create.get()
        fnDescr = self._getExtraPath(FN_PHANTOM_DESCR)
        fhDescr = open(fnDescr, 'w')
        fhDescr.write(desc)
        fhDescr.close()
        dim = [desc.split()[0], desc.split()[1], desc.split()[2]]

        self.runJob("xmipp_phantom_create", " -i %s -o %s" % (fnDescr, fnVol))
        setMRCSamplingRate(fnVol, self.sampling.get())
        return dim, fnVol


    def applyRandomOrientation(self, fnIn, fnOut):
        rot = 0
        tilt = 0
        psi = 0
        shiftX = 0
        shiftY = 0
        shiftZ = 0

        rotErr = 0
        tiltErr = 0

        if self.uniformAngularDistribution:
            rot = 2*np.pi * np.random.uniform(0, 1)*180/np.pi
            tilt = np.arccos(2*np.random.uniform(0, 1) - 1)*180/np.pi

            #It is neccesary to create a new variable because of the random errors
            rotErr = rot
            tiltErr = tilt

            if self.stdError:
                rotErr = rot + np.random.normal(0, self.sigma.get())
                tiltErr = tilt + np.random.normal(0, self.sigma.get())
        else:
            if self.rotate:
                rot = np.random.randint(self.rotmin.get(), self.rotmax.get())
                tilt = np.random.randint(self.tiltmin.get(), self.tiltmax.get())
                psi = np.random.randint(self.psimin.get(), self.psimax.get())
                rotErr = rot
                tiltErr = rot

        if self.applyShift:
            # Shifts
            shiftX = np.random.randint(self.xmin.get(), self.xmax.get())
            shiftY = np.random.randint(self.ymin.get(), self.ymax.get())
            shiftZ = np.random.randint(self.zmin.get(), self.zmax.get())

        if self.rotate or self.applyShift:
            self.runJob("xmipp_transform_geometry",
                        " -i %s -o %s --rotate_volume euler %d %d %d --shift %d %d %d --dont_wrap"
                        % (fnIn, fnOut, rotErr, tiltErr, psi, shiftX, shiftY, shiftZ))

        return rot, tilt, psi, shiftX, shiftY, shiftZ


    def addNoiseToPhantom(self, fnIn, fnOut):
        if self.addNoise.get():
            params_noise = ' -i %s ' % fnIn
            if self.option == 0:
                params_noise += ' --save_metadata_stack'
            import random
            if self.differentStatistics.get():
                sigmaNoise = random.uniform(self.minstd.get(), self.maxstd.get())
                params_noise += ' --type gaussian %f ' % sigmaNoise
            else:
                meanNoise = self.meanNoise.get()
                sigmaNoise = self.stdNoise.get()
                params_noise += ' --type gaussian %f %f ' % (sigmaNoise, meanNoise)

            params_noise += ' -o %s ' % fnOut
            self.runJob('xmipp_transform_add_noise', params_noise)


    def creatingprojections(self, mystack):

        params_phantom = ' -i %s ' % self.inputVolume.get().getFileName()
        params_phantom += ' --method real_space '
        params_phantom += ' --params %s ' % (self._getExtraPath(FN_PARAMS))
        params_phantom += ' --sampling_rate %f ' % (self.sampling.get())
        params_phantom += ' -o %s ' % mystack
        self.runJob('xmipp_phantom_project', params_phantom)

    def reconstructSubtomo(self, mystack, subtomoRecosntruct):

        params_fourier = ' -i %s ' % mystack
        params_fourier += ' -o %s ' % subtomoRecosntruct
        params_fourier += ' -thr 4'
        self.runJob('xmipp_reconstruct_fourier', params_fourier)

    def createSubtomogramsSimulatingTiltSeries(self):
        # This function defines the output to be updated as soon as the subtomos are created
        self.createOutputSet(self.inputVolume.get().getDim())

        self.createParamsFile()

        for idx in range(self.nsubtomos.get()):
            # The inputs and outputs of the steps are defined
            mystack = self._getExtraPath('projectionstack_' + str(idx) + '.xmd')
            subtomoRecosntruct = self._getExtraPath('subtomo_' + str(idx) + MRC_EXT)

            # The different steps to create the noisy phantoms are launched here
            self.creatingprojections(mystack)
            #self.addNoiseToPhantom(mystack, mystack2)
            self.reconstructSubtomo(mystack, subtomoRecosntruct)

            # The subtomogram is stored as Scipion object
            self._addSubtomogram(None, None, subtomoRecosntruct, 0, 0, 0, 0, 0, 0)

    def applyMissingWedge(self, mwangle, fnIn, fnOut):
        self.runJob("xmipp_transform_filter", " --fourier wedge -%d %d 0 0 0 -i %s -o %s"
                    % (mwangle, mwangle, fnIn, fnOut))

    def createOutputSet(self, dim):
        self.outputSet = self._createSetOfSubTomograms(self._getOutputSuffix(SetOfSubTomograms))
        self.outputSet.setDim(dim)
        self.outputSet.setSamplingRate(self.sampling.get())

    def _addSubtomogram(self, tomo, acq, phantomfn, rot, tilt, psi, shiftX, shiftY, shiftZ):
        """ Creates and adds a the phantom subtomogram to the set. It creates the coordinate as well if active"""
        subtomo = SubTomogram()
        subtomo.setAcquisition(acq)
        subtomo.setLocation(phantomfn)
        subtomo.setSamplingRate(self.sampling.get())
        subtomo.setObjComment(
            "Angle Rot, tilt, psi: %d, %d, %d \nShifts X,Y,Z: %d, %d, %d" % (rot, tilt, psi, shiftX, shiftY, shiftZ))

        subtomo.phantom_rot = Integer(rot)
        subtomo.phantom_tilt = Integer(tilt)
        subtomo.phantom_psi = Integer(psi)
        subtomo.phantom_shiftX = Integer(shiftX)
        subtomo.phantom_shiftY = Integer(shiftY)
        subtomo.phantom_shiftZ = Integer(shiftZ)

        # Scipion alignment matrix
        A = euler_matrix(np.deg2rad(psi), np.deg2rad(tilt), np.deg2rad(rot), 'szyz')
        shifts = np.transpose(np.array([-shiftX, -shiftY, -shiftZ]))

        # Translation matrix
        T = np.eye(4)
        T[:3, 3] = shifts
        M = A @ T

        transform = Transform()
        transform.setMatrix(M)
        subtomo.setTransform(transform)

        self._addCoordinate(subtomo, tomo)

        self.outputSet.append(subtomo)

    def _addCoordinate(self, subtomo, tomo):
        """ Adds a Coordinate3D (if apply) to the coordinate set and fills the subtomogram with the coordinate"""

        if self.generateCoordinates():
            self.debug("Adding a coordinate to subtomo %s" % subtomo)
            coor = Coordinate3D()
            coor.setVolume(tomo)
            tomoDim = tomo.getDim()
            coor.setX(np.random.randint(0, tomoDim[0]), const.BOTTOM_LEFT_CORNER)
            coor.setY(np.random.randint(0, tomoDim[1]), const.BOTTOM_LEFT_CORNER)
            coor.setZ(np.random.randint(0, tomoDim[2]), const.BOTTOM_LEFT_CORNER)

            self.coordsSet.append(coor)
            self.coordsSet.setBoxSize(subtomo.getDim()[0])
            subtomo.setCoordinate3D(coor)
            subtomo.setVolName(tomo.getFileName())

    def generateCoordinates(self):
        return self.coords.get()

    def createOutputStep(self):
        self._defineOutputs(outputSubtomograms=self.outputSet)
        if self.option.get() == 0:
            self._defineSourceRelation(self.inputVolume.get(), self.outputSet)
        if self.generateCoordinates():
            self._defineOutputs(outputCoord=self.coordsSet)
            self._defineSourceRelation(self.tomos.get(), self.outputSet)


    def createParamsFile(self):
        fn_params = self._getExtraPath(FN_PARAMS)
        dim = self.inputVolume.get().getDim()

        f = open(fn_params, 'w')
        f.write('# XMIPP_STAR_1 *\n')
        f.write('# Projection Parameters\n')
        f.write('data_block1\n')
        f.write('# X and Y projection dimensions [Xdim Ydim]\n')
        f.write('_dimensions2D   \'%d %d\' \n' % (dim[0], dim[1]))
        f.write('# Rotation range and number of samples [Start Finish NSamples]\n')
        f.write('_projRotRange    \'%d %d %d\' \n' % (self.rotmin.get(), self.rotmax.get(), 1))
        f.write('# Rotation angle added noise  [noise (bias)]\n')
        f.write('_projRotNoise   \'%d\'\n' % self.rotStd.get())
        f.write('# Tilt range and number of samples for Tilt\n')
        f.write('_projTiltRange    \'%d %d %d\' \n' % (self.tiltmin.get(), self.tiltmax.get(),
                                                       round(
                                                           abs(self.tiltmax.get() - self.tiltmin.get()) / self.angularSampling.get())))
        f.write('# Tilt angle added noise\n')
        f.write('_projTiltNoise   \'%d\' \n' % self.tiltStd.get())
        f.write('# Psi range and number of samples\n')
        f.write('_projPsiRange    \'0 0 0\'\n')
        f.write('# Psi added noise\n')
        f.write('_projPsiNoise   \'0\'\n')
        f.write('# Noise applied to pixels [noise (bias)]\n')

        import random
        sigmaNoise  = 0
        if self.addNoise.get():
            if self.differentStatistics.get():
                sigmaNoise = random.uniform(self.minstd.get(), self.maxstd.get())
            else:
                sigmaNoise = self.stdNoise.get()

        f.write('_noisePixelLevel   \'%d\'\n' % sigmaNoise)
        f.write('# Noise applied to particle center coordenates [noise (bias)]\n')
        f.write('_noiseCoord   \'0\'\n')
        f.close()

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if self.rotmin.get() > self.rotmax.get():
            errors.append("rot max must be bigger than rot min")
        if self.tiltmin.get() > self.tiltmax.get():
            errors.append("tilt max must be bigger than tilt min")
        if self.psimin.get() > self.psimax.get():
            errors.append("psi max must be bigger than psi min")
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputSubtomograms'):
            summary.append("Output phantom not ready yet.")
        else:
            summary.append("%s phantoms created with random orientations" % self.nsubtomos.get())
            if self.mwfilter.get():
                summary.append("Missing wedge applied between +-%d along Y axis" % self.mwangle.get())
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputSubtomograms'):
            methods.append("Output phantoms not ready yet.")
            return methods
        else:
            methods.append("%s phantoms created with random orientations." % self.nsubtomos.get())
            if self.mwfilter.get():
                methods.append("Missing wedge applied between +-%d along Y axis." % self.mwangle.get())
            return methods

