# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
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
import glob
from pwem.emlib import lib
import numpy as np
from pwem.convert.transformations import euler_matrix
from pwem.objects.data import Transform, Integer
from pwem.protocols import EMProtocol
from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from pyworkflow import BETA
from pyworkflow.protocol.params import IntParam, FloatParam, LEVEL_ADVANCED, PointerParam, BooleanParam, StringParam, \
    GT, GE, Range, EnumParam, USE_GPU, GPU_LIST
from tomo.protocols import ProtTomoBase
import pwem.emlib.metadata as md

import pwem
from xmipp3.convert import writeSetOfVolumes
import xmipp3
from xmipp3.convert import writeSetOfVolumes
from ..utils import calculateRotationAngleAndShiftsFromTM, setGeometricalParametersToRow
from ..objects import SetOfTiltSeriesParticle, TiltSeriesParticle, TiltParticle

from pwem.objects import Volume

FN_INPUTPARTICLES = 'ts_'
XMD_EXT = '.xmd'
MRCS_EXT = '.mrcs'
SUFFIX_CTF_CORR = '_ctf_corrected'
SUFFIX_ZEROPARTICLES = 'zeroParticles'


class XmippProtReconstructInitVol(EMProtocol, ProtTomoBase, xmipp3.XmippProtocol):
    """ This protocol performs initial volumes for subtomogram averaging """

    _label = 'initial volume'
    _devStatus = BETA
    _conda_env = 'xmipp_swiftalign'
    tsIdList = []
    GRADIENT_DESCENT_OPTIONS = [
        'None',
        'RMSProp']

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStacks', PointerParam,
                      pointerClass="SetOfTiltSeriesParticle",
                      label='Particles',
                      help="Input Tilt series particles")

        form.addParam('provideInitialVolume', BooleanParam, default=False,
                      label='Provide Initial Volume?',
                      help="blablablbala")

        form.addParam('initVol', PointerParam, condition='provideInitialVolume',
                      pointerClass="Volume",
                      label='Initial Volume',
                      help="Input Tilt series particles")

        form.addParam('correctCTF', BooleanParam, default=True,
                      label='Correct CTF?',
                      help="The Set of tilt series particles have a CTF associated. "
                           "Set Yes for correcting th CTF. No will not correct the CTF.")

        form.addParam('padding_factor', FloatParam, default=2,
                      label='Padding factor',
                      help="Paddign factor for CTF wiener correction.")

        form.addParam('wiener_constant', FloatParam, default=-1, expertLevel=LEVEL_ADVANCED,
                      label="Wiener constant",
                      help=' Wiener-filter constant (if < 0: use FREALIGN default)')
        form.addParam('correctEnvelope', BooleanParam, default='False', expertLevel=LEVEL_ADVANCED,
                      label="Correct for CTF envelope",
                      help=' Only in cases where the envelope is well estimated correct for it')

        form.addSection(label='CTF')
        form.addParam('considerInputCtf', BooleanParam, label='Consider CTF',
                      default=True,
                      help='Consider the CTF of the particles')

        form.addSection(label='Refinement')
        form.addParam('reconstructLast', BooleanParam, label='Reconstruct last volume', default=True)
        form.addParam('numberOfIterations', IntParam, label='Number of iterations', default=3, validators=[GT(0)])
        form.addParam('numberOfLocalIterations', IntParam, label='Number of local iterations', default=1,
                      validators=[GT(0)])
        form.addParam('numberOfAlignmentRepetitions', IntParam, label='Number of repetitions', default=2,
                      validators=[GT(0)])
        form.addParam('maximumResolution', FloatParam, label="Maximum alignment resolution (A)", default=8.0,
                      validators=[GT(0)],
                      help='Image comparison resolution limit of the refinement')
        form.addParam('nextResolutionCriterion', FloatParam, label="FSC criterion", default=0.5,
                      validators=[Range(0, 1)],
                      expertLevel=LEVEL_ADVANCED,
                      help='The resolution of the reconstruction is defined as the inverse of the frequency at which ' \
                           'the FSC drops below this value. Typical values are 0.143 and 0.5')
        form.addParam('initialMaxPsi', FloatParam, label='Maximum psi (deg)', default=180.0, validators=[Range(0, 180)],
                      expertLevel=LEVEL_ADVANCED,
                      help='Maximum psi parameter of the particles')
        form.addParam('initialMaxShift', FloatParam, label='Maximum shift (px)', default=16.0, validators=[GT(0)],
                      help='Maximum shift of the particle in pixels')
        form.addParam('useAutomaticStep', BooleanParam, label='Use automatic step', default=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='Automatically determine the step used when exploring the projection landscape')
        stepGroup = form.addGroup('Steps', condition='not useAutomaticStep', expertLevel=LEVEL_ADVANCED)
        stepGroup.addParam('angleStep', FloatParam, label='Angle step (deg)', default=5.0, validators=[GT(0)])
        stepGroup.addParam('shiftStep', FloatParam, label='Shift step (px)', default=2.0, validators=[GT(0)])

        form.addParam('gradientDescent', EnumParam, label='Gradient descent', choices=self.GRADIENT_DESCENT_OPTIONS,
                      default=0,
                      help='If selected, each iteration will be subdivided in multiple iterations of a gradient descent')
        form.addParam('minibatchSize', IntParam, label='Mini-batch size', default=16384,
                      condition='gradientDescent > 0', validators=[GT(0)],
                      help='Size of the minibatch when performing a gradient descent')
        rmspropGroup = form.addGroup('RMSProp', condition='gradientDescent == 1', expertLevel=LEVEL_ADVANCED)
        rmspropGroup.addParam('rmspropGamma', FloatParam, label='Gradient learning rate', default=0.9,
                              validators=[Range(0, 1)])
        rmspropGroup.addParam('rmspropNu', FloatParam, label='Learning rate', default=0.001, validators=[GE(0)])
        rmspropGroup.addParam('rmspropEps', FloatParam, label='Stability factor', default=1e-8, validators=[GE(0)])

        form.addParam('reconstructPercentage', FloatParam, label='Reconstruct percentage (%)', default=50,
                      validators=[Range(0, 100)],
                      help='Percentage of best particles used for reconstruction')
        form.addParam('numberOfMatches', IntParam, label='Number of matches', default=1, validators=[GE(0)],
                      help='Number of reference matches for each particle')

        form.addSection(label='Compute')
        form.addParam('databaseRecipe', StringParam, label='Database recipe',
                      default='OPQ48_192,IVF32768,PQ48', expertLevel=LEVEL_ADVANCED,
                      help='FAISS database structure. Please refer to '
                           'https://github.com/facebookresearch/faiss/wiki/The-index-factory')
        form.addParam('useFloat16', BooleanParam, label='Use float16', default=False, expertLevel=LEVEL_ADVANCED,
                      help='When enabled, FAISS will be prompted to use half precision floating point '
                           'numbers. This may improve performance and/or memory footprint at some '
                           'accuracy cost. Only supported for GPUs')
        form.addParam('usePrecomputed', BooleanParam, label='Precompute centroid distances', default=False,
                      help='When using PQ encoding precompute pairwise distances between centroids')
        form.addParam('databaseTrainingSetSize', IntParam, label='Database training set size',
                      default=int(2e6),
                      help='Number of data-augmented particles to used when training the database',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('databaseMaximumSize', IntParam, label='Database size limit',
                      default=int(10e6),
                      help='Maximum number of elements that can be stored in the database '
                           'before performing an alignment and flush')
        form.addParam('batchSize', IntParam, label='Batch size',
                      default=8192,
                      help='It is recommended to use powers of 2. Using numbers around 8192 works well')
        form.addParam('copyParticles', BooleanParam, label='Copy particles to scratch', default=False,
                      help='Copy input particles to scratch directory. Note that if input file format is '
                           'incompatible the particles will be converted into scratch anyway')

        form.addHidden(USE_GPU, BooleanParam, default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation.\
                       Select the one you want to use.")
        form.addHidden(GPU_LIST, StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):

        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.wienerCorrectionZeroDegStep)
        self._insertFunctionStep(self.wienerCorrectionStep)

        for iteration in range(0, self.numberOfIterations.get()):
            iterFolder = self._getExtraPath('iteration_%04d' % iteration)
            if not os.path.isdir(iterFolder):
                os.mkdir(iterFolder)
            self._insertFunctionStep(self.createGalleryStep, iteration)
            self._insertFunctionStep(self.trainDatabaseStep, iteration)
            self._insertFunctionStep(self.aligningFaissStep, iteration)
            # self._insertFunctionStep(self.validateAligningStep)
            # self._insertFunctionStep(self.reconstructionStep)
        # self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        # To obtain a list with the tsIds and write the particles with such tsId as a metadata (.xmd) file
        # These .xmd will be used as input in the wienerCorrectionStep to correct the CTF
        sotsp = self.inputStacks.get()

        # List of tsIds
        self.tsIdList = []

        for tp in sotsp.iterItems():
            tsId = str(tp.getOriginalTs())
            if tsId not in self.tsIdList:
                # If not, add it to the list
                os.mkdir(self._getExtraPath(tsId))
                self.tsIdList.append(tsId)

        psi = 0.0
        mdtpZeroTilt = lib.MetaData()
        fnParticlesZero = self._getExtraPath(SUFFIX_ZEROPARTICLES) + XMD_EXT
        for tsId in self.tsIdList:
            mdtp = lib.MetaData()

            fnParticles = os.path.join(self._getExtraPath(tsId), tsId + XMD_EXT)

            for tsparticle in sotsp:
                tsIdOrig = str(tsparticle.getOriginalTs())
                if tsId == tsIdOrig:
                    # A huge an non sense tilt number
                    mintilt = 1e38
                    for tp in tsparticle:
                        tilt = tp.getTiltAngle()
                        rot, sx, sy = calculateRotationAngleAndShiftsFromTM(tp)
                        nRow = md.Row()
                        nRowZeroTilt = md.Row()
                        fn = tp.getFileName()
                        ctf = tp.getCTF()
                        defU = ctf.getDefocusU()
                        defV = ctf.getDefocusV()
                        defAng = ctf.getDefocusAngle()
                        nRow = setGeometricalParametersToRow(nRow, fn, rot, tilt, psi, sx, sy, defU, defV, defAng)
                        if np.abs(tilt) < mintilt:
                            mintilt = np.abs(tilt)
                            fnZero = fn
                            rotZero, tiltZero, psiZero, sxZero, syZero, defUZero, defVZero, defAngZero = rot, tilt, psi, sx, sy, defU, defV, defAng
                        nRow.addToMd(mdtp)
                    nRowZeroTilt = setGeometricalParametersToRow(nRowZeroTilt, fnZero, rotZero, tiltZero, psiZero,
                                                                 sxZero, syZero, defUZero, defVZero, defAngZero)

                    nRowZeroTilt.addToMd(mdtpZeroTilt)
            mdtp.write(fnParticles)
        mdtpZeroTilt.write(fnParticlesZero)

    def mergeIndividualXmdFiles(self):
        pass
        # for tsId in self.tsIdList:
        #    and

    def wienerJob(self, inputdata, outputfolder, svmetadata):
        # TODO: Check if is phase flipped
        sampling = self.inputStacks.get().getSamplingRate()
        paddingFactor = self.padding_factor.get()
        wc = self.wiener_constant.get()
        isphaseflipped = self.inputStacks.get().isPhaseFlipped()
        corrEnvelope = self.correctEnvelope
        nMpi = self.numberOfMpi.get()
        nThr = self.numberOfThreads.get()

        params = '  -i %s' % inputdata
        params += '  -o %s' % outputfolder
        params += '  --save_metadata_stack %s' % svmetadata
        params += '  --pad %s' % paddingFactor
        params += '  --wc %s' % wc
        params += '  --sampling_rate %s' % sampling
        if isphaseflipped:
            params += '  --phase_flipped '
        if corrEnvelope:
            params += '  --correct_envelope '
        self.runJob('xmipp_ctf_correct_wiener2d', params, numberOfMpi=nMpi, numberOfThreads=nThr)

    def wienerCorrectionZeroDegStep(self):
        inputdata = self._getExtraPath(SUFFIX_ZEROPARTICLES + XMD_EXT)
        outputfolder = self._getExtraPath(SUFFIX_ZEROPARTICLES + SUFFIX_CTF_CORR + MRCS_EXT)
        svmetadata = self._getExtraPath(SUFFIX_ZEROPARTICLES + SUFFIX_CTF_CORR + XMD_EXT)

        self.wienerJob(inputdata, outputfolder, svmetadata)

    def wienerCorrectionStep(self):
        for tsId in self.tsIdList:
            inputdata = os.path.join(self._getExtraPath(tsId), tsId + XMD_EXT)
            outputfolder = os.path.join(self._getExtraPath(tsId), tsId + SUFFIX_CTF_CORR + MRCS_EXT)
            svmetadata = os.path.join(self._getExtraPath(tsId), tsId + SUFFIX_CTF_CORR + XMD_EXT)

            self.wienerJob(inputdata, outputfolder, svmetadata)

    def createGalleryStep(self, iteration):
        # TODO: Check if an initial volume is not provided
        # TODO: Adjust resolution in projection
        if self.provideInitialVolume.get():
            fnVol = self.initVol.get().getFileName()
        else:
            fnVol = ''
        boxsize = self.inputStacks.get().getFirstItem().getBoxSize()
        sampling = self.inputStacks.get().getSamplingRate()
        angular_sampling = np.arctan2(2.0 * self.maximumResolution.get(), float(boxsize) * sampling) * 180 / 3.14159

        params = ' -i %s' % fnVol
        params += ' -o %s' % os.path.join(self._getExtraPath('iteration_%04d' %iteration), 'gallery_%s.mrcs' % iteration)
        params += ' --sampling_rate %f' % angular_sampling

        self.runJob('xmipp_angular_project_library', params)

    def writeParticleStackToMd(self, sotsp, fn):

        md = lib.MetaData()
        psi = 0.0
        for tsparticle in sotsp:
            for tp in tsparticle:
                tilt = tp.getTiltAngle()
                rot, sx, sy = calculateRotationAngleAndShiftsFromTM(tp)
                nRow = md.Row()
                nRow.setValue(lib.MDL_IMAGE, fn)
                defU = tp.getDefocusU()
                defV = tp.getDefocusV()
                defAng = tp.getDefocusAngle()
                nRow.setValue(lib.MDL_CTF_DEFOCUSU, defU)
                nRow.setValue(lib.MDL_CTF_DEFOCUSV, defV)
                nRow.setValue(lib.MDL_CTF_DEFOCUS_ANGLE, defAng)

                nRow.setValue(lib.MDL_ANGLE_TILT, tilt)
                nRow.setValue(lib.MDL_ANGLE_ROT, rot)
                nRow.setValue(lib.MDL_ANGLE_PSI, psi)
                nRow.setValue(lib.MDL_SHIFT_X, sx)
                nRow.setValue(lib.MDL_SHIFT_Y, sy)
                nRow.addToMd(md)
        md.write(fn)

    def setInfoFirstIteration(self, iteration, fnSetInfo):
        sampling = self.inputStacks.get().getSamplingRate()
        maxfreqDig = sampling/self.maximumResolution.get()
        self.setInfoIteration(iteration, self.maximumResolution.get(), maxfreqDig,
                               self.initialMaxPsi.get(), self.initialMaxShift.get(),
                               self.shiftStep.get(), self.angleStep.get(), fnSetInfo)

    def setInfoIteration(self, iteration, maxResolution, maxfreqDig, maxPsi, maxShift, shiftStep, angleStep, fnSetInfo):
        # Write to metadata
        md = lib.MetaData()
        id = md.addObject()
        md.setValue(lib.MDL_RESOLUTION_FREQ, maxResolution, id)
        md.setValue(lib.MDL_RESOLUTION_FREQREAL, maxfreqDig, id)
        md.setValue(lib.MDL_ANGLE_PSI, maxPsi, id)
        md.setValue(lib.MDL_SHIFT_X, maxShift, id)
        md.setValue(lib.MDL_SHIFT_Y, maxShift, id)
        md.setValue(lib.MDL_SHIFT_DIFF, shiftStep, id)
        md.setValue(lib.MDL_ANGLE_DIFF, angleStep, id)
        md.write(fnSetInfo)

    def getParticleSize(self):
        return self.inputStacks.get().getFirstItem().getBoxSize()

    def trainDatabaseStep(self, iteration: int):
        #TODO: Change resolution and freqDig when will be updated
        trainingSize = int(self.databaseTrainingSetSize)
        recipe = self.databaseRecipe

        imageSize = self.getParticleSize()

        fnSetInfo = self._getIterationParametersFilename(iteration)

        if iteration == 0:
            self.setInfoFirstIteration(iteration, fnSetInfo)
            maxShiftPx = self.initialMaxShift.get()
            maxPsi = self.initialMaxPsi.get()
            sampling = self.inputStacks.get().getSamplingRate()
            maxFrequencyDig = sampling / self.maximumResolution.get()
        else:
            self.setInfoIteration(iteration, fnSetInfo)
            md = lib.MetaData(fnSetInfo)
            maxFrequencyDig = md.getValue(lib.MDL_RESOLUTION_FREQREAL, 1)
            maxPsi = md.getValue(lib.MDL_ANGLE_PSI, 1)
            maxShiftPx = md.getValue(lib.MDL_SHIFT_X, 1)

        maxShift = maxShiftPx / float(imageSize)

        fnGallery = os.path.join(self._getExtraPath('iteration_%04d' %iteration), 'gallery_%s.doc' % iteration)

        args = []
        args += ['-i', fnGallery]
        args += ['-o', self._getTrainingIndexFilename(iteration)]
        args += ['--recipe', recipe]
        # args += ['--weights', self._getWeightsFilename(iteration)]
        args += ['--max_shift', maxShift]  # FIXME fails with != 180
        args += ['--max_psi', maxPsi]
        args += ['--max_frequency', maxFrequencyDig]
        args += ['--training', trainingSize]
        args += ['--batch', self.batchSize]
        args += ['--scratch', self._getTrainingScratchFilename()]
        if self.useGpu:
            args += ['--device'] + self._getDeviceList()
        if self.useFloat16:
            args += ['--fp16']
        if self.usePrecomputed:
            args += ['--use_precomputed']

        env = self.getCondaEnv()
        env['LD_LIBRARY_PATH'] = ''  # Torch does not like it
        self.runJob('xmipp_swiftalign_train', args, numberOfMpi=1, env=env)

    def _getIterationPath(self, iteration: int, *paths):
        return self._getExtraPath('iteration_%04d' % iteration, *paths)

    def _getAlignmentRepetitionMdFilename(self, iteration: int):
        return self._getIterationPath(iteration, 'aligned.xmd')

    def _getIterationParametersFilename(self, iteration: int):
        return self._getIterationPath(iteration, 'params.xmd')

    def _getTrainingScratchFilename(self):
        return self._getTmpPath('scratch.bin')

    def _getTrainingMdFilename(self, iteration: int):
        return self._getIterationPath(iteration, 'training.xmd')

    def _getTrainingIndexFilename(self, iteration: int):
        return self._getIterationPath(iteration, 'database.idx')

    def _getGalleryMdFilename(self, iteration: int, repetition: int):
        return self._getIterationPath(iteration, 'gallery%04d.xmd' % repetition)

    def _getDeviceList(self):
        gpus = self.getGpuList()
        return list(map('cuda:{:d}'.format, gpus))

    def getAlignmentInfo(self, iteration):
        md = lib.MetaData(self._getIterationParametersFilename(iteration))
        imageSize = self.getParticleSize()
        maxFrequency = md.getValue(lib.MDL_RESOLUTION_FREQREAL, 1)
        maxPsi = md.getValue(lib.MDL_ANGLE_PSI, 1)
        maxShiftPx = md.getValue(lib.MDL_SHIFT_X, 1)
        maxShift = maxShiftPx / float(imageSize)
        nShift = round((2 * maxShiftPx) / md.getValue(lib.MDL_SHIFT_DIFF, 1)) + 1
        nRotations = round(360 / md.getValue(lib.MDL_ANGLE_DIFF, 1))

        return maxFrequency, maxPsi, maxShift, nShift, nRotations

    def aligningFaissStep(self, iteration):

        maxFrequency, maxPsi, maxShift, nShift, nRotations = self.getAlignmentInfo(iteration)

        if iteration == 0:
            inputMdFilename = self._getExtraPath(SUFFIX_ZEROPARTICLES) + XMD_EXT
        else:
            inputMdFilename = self._getIterationInputParticleMdFilename(iteration)
        fnGallery = os.path.join(self._getExtraPath('iteration_%04d' % iteration), 'gallery_%s.doc' % iteration)

        local = 0
        import math
        localFactor = math.pow(2, -local)
        maxPsi *= localFactor
        maxShift *= localFactor

        # Perform the alignment
        args = []
        args += ['-i', inputMdFilename]
        args += ['-o', self._getAlignmentRepetitionMdFilename(iteration)]
        args += ['-r', fnGallery]
        args += ['--index', self._getTrainingIndexFilename(iteration)]
        # args += ['--weights', self._getWeightsFilename(iteration)]
        args += ['--max_shift', maxShift]
        args += ['--max_psi', maxPsi]
        args += ['--rotations', nRotations]
        args += ['--shifts', nShift]
        args += ['--max_frequency', maxFrequency]
        args += ['--dropna']
        args += ['--batch', self.batchSize]
        args += ['--max_size', self.databaseMaximumSize]
        args += ['-k', self.numberOfMatches]
        args += ['--reference_labels', 'angleRot', 'angleTilt', 'ref3d', 'imageRef']
        if self.useGpu:
            args += ['--device'] + self._getDeviceList()
        if local > 0:
            args += ['--local']
        if self.usePrecomputed:
            args += ['--use_precomputed']

        env = self.getCondaEnv()
        env['LD_LIBRARY_PATH'] = ''  # Torch does not like it
        self.runJob('xmipp_swiftalign_query', args, numberOfMpi=1, env=env)


    def validateAligningStep(self, objId):

        stack = self.inputStacks.get()[objId]
        fn = self._getExtraPath('caca.xmd')
        self.writeParticleStackToMd(stack, fn)
        params = ' -i %s ' % self._getExtraPath(fn)
        if self.correctCTF:
            params += ' --useCTF '
        params += ' --sampling %f ' % (self.subtomos.get().getSamplingRate())
        params += ' -o %s ' % self._getExtraPath()
        params += ' --thr %d ' % self.numberOfThreads.get()

        self.runJob('xmipp_tomo_align_subtomo_stack', params)

    def reconstructionStep(self, objId):

        stack = self.inputStacks.get()[objId]
        fn = self._getExtraPath('caca.xmd')
        self.writeParticleStackToMd(stack, fn)
        params = ' -i %s ' % self._getExtraPath(fn)
        if self.correctCTF:
            params += ' --useCTF '
        params += ' --sampling %f ' % (self.subtomos.get().getSamplingRate())
        params += ' -o %s ' % self._getExtraPath()
        params += ' --thr %d ' % self.numberOfThreads.get()

        self.runJob('xmipp_reconstruct_fourier', params)

    def createOutputStep(self):
        pass

        # --------------------------- INFO functions --------------------------------------------

    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        methods = []
        return methods
