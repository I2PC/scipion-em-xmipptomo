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
from pyworkflow.protocol.params import IntParam, FloatParam, LEVEL_ADVANCED, PointerParam, BooleanParam
from tomo.protocols import ProtTomoBase
import pwem.emlib.metadata as md

import pwem
from xmipp3.convert import writeSetOfVolumes
from ..utils import calculateRotationAngleAndShiftsFromTM, setGeometricalParametersToRow
from ..objects import SetOfTiltSeriesParticle, TiltSeriesParticle, TiltParticle

from pwem.objects import Volume

FN_INPUTPARTICLES = 'ts_'
XMD_EXT = '.xmd'
MRCS_EXT = '.mrcs'
SUFFIX_CTF_CORR = '_ctf_corrected'


class XmippProtReconstructInitVol(EMProtocol, ProtTomoBase):
    """ This protocol performs initial volumes for subtomogram averaging """

    _label = 'initial volume'
    _devStatus = BETA
    tsIdList = []

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

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):

        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.wienerCorrectionStep)
        self._insertFunctionStep(self.createGalleryStep)
        self._insertFunctionStep(self.aligningStep)
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
        for tsId in self.tsIdList:
            mdtp = lib.MetaData()
            mdtpZeroTilt = lib.MetaData()
            fnParticles = os.path.join(self._getExtraPath(tsId), tsId + XMD_EXT)
            fnParticlesZero = os.path.join(self._getExtraPath(tsId), tsId +'_zeroTilt'+ XMD_EXT)
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
                        nRow.setValue(lib.MDL_IMAGE, fn)
                        ctf = tp.getCTF()
                        defU = ctf.getDefocusU()
                        defV = ctf.getDefocusV()
                        defAng = ctf.getDefocusAngle()
                        nRow = setGeometricalParametersToRow(nRow, rot, tilt, psi, sx, sy, defU, defV, defAng)
                        if np.abs(tilt)<mintilt:
                            mintilt = tilt
                            rotZero, tiltZero, psiZero, sxZero, syZero, defUZero, defVZero, defAngZero = rot, tilt, psi, sx, sy, defU, defV, defAng
                            nRowZeroTilt = setGeometricalParametersToRow(nRowZeroTilt, rotZero, tiltZero, psiZero,
                                                                         sxZero, syZero, defUZero, defVZero, defAngZero)
                            nRowZeroTilt.addToMd(mdtpZeroTilt)
                        nRow.addToMd(mdtp)

            mdtp.write(fnParticles)
            mdtpZeroTilt.write(fnParticlesZero)

    def mergeIndividualXmdFiles(self):
        for tsId in self.tsIdList:
            and


    def wienerCorrectionStep(self):
        # TODO: Check if is phase flipped
        sampling = self.inputStacks.get().getSamplingRate()

        for tsId in self.tsIdList:
            params = '  -i %s' % os.path.join(self._getExtraPath(tsId), tsId + XMD_EXT)
            params += '  -o %s' % os.path.join(self._getExtraPath(tsId), tsId + SUFFIX_CTF_CORR + MRCS_EXT)
            params += '  --save_metadata_stack %s' % os.path.join(self._getExtraPath(tsId),
                                                                  tsId + SUFFIX_CTF_CORR + XMD_EXT)
            params += '  --pad %s' % self.padding_factor.get()
            params += '  --wc %s' % self.wiener_constant.get()
            params += '  --sampling_rate %s' % sampling

            if self.inputStacks.get().isPhaseFlipped():
                params += '  --phase_flipped '

            if self.correctEnvelope:
                params += '  --correct_envelope '

            nproc = self.numberOfMpi.get()
            nT = self.numberOfThreads.get()

            self.runJob('xmipp_ctf_correct_wiener2d', params, numberOfMpi=nproc, numberOfThreads=nT)




    def createGalleryStep(self):
        #TODO: Check if an initial volume is not provided
        sampling = self.inputStacks.get().getSamplingRate()

        if self.provideInitialVolume.get():
            fnVol = self.initVol.get().getFileName()
        else:
            fnVol = ''

        params = ' -i %s' % fnVol
        params += ' -o %s' % self._getExtraPath('gallery.mrcs')
        params += ' --sampling_rate %f' % sampling
        params += ' --psi_sampling %f' % 5

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

        '''
        for ti in sotsp:
            tilt = ti.getTiltAngle()
            rot, sx, sy = calculateRotationAngleAndShiftsFromTM(ti)
            nRow = md.Row()
            nRow.setValue(lib.MDL_IMAGE, fn)
            defU = ti.getDefocusU()
            defV = ti.getDefocusV()
            defAng = ti.getDefocusAngle()
            nRow.setValue(lib.MDL_CTF_DEFOCUSU, defU)
            nRow.setValue(lib.MDL_CTF_DEFOCUSV, defV)
            nRow.setValue(lib.MDL_CTF_DEFOCUS_ANGLE, defAng)

            nRow.setValue(lib.MDL_ANGLE_TILT, tilt)
            nRow.setValue(lib.MDL_ANGLE_ROT, rot)
            nRow.setValue(lib.MDL_ANGLE_PSI, psi)
            nRow.setValue(lib.MDL_SHIFT_X, sx)
            nRow.setValue(lib.MDL_SHIFT_Y, sy)
            nRow.addToMd(md)
        '''
        md.write(fn)

    def aligningStep(self, objId):

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
