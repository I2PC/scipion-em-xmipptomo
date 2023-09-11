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
import pwem
from xmipp3.convert import writeSetOfVolumes
from ..utils import calculateRotationAngleAndShiftsFromTM
from ..objects import SetOfTiltSeriesParticle, TiltSeriesParticle, TiltParticle

from pwem.objects import Volume


FN_INPUTPARTICLES = 'ts_'
XMD_EXT = '.xmd'
MRC_EXT = '.mrc'
SUFFIX_CTF_CORR = '_ctf_corrected'



class XmippProtReconstructInitVol(EMProtocol, ProtTomoBase):
    """ This protocol performs initial volumes for subtomogram averaging """

    _label = 'initial volume'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStacks', PointerParam,
                      pointerClass="SetOfTiltSeriesParticle, SetOfParticles",
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

        form.addParallelSection(threads = 1, mpi = 0)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):

        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.wienerCorrectionStep)
        self._insertFunctionStep(self.createGalleryStep)
        #self._insertFunctionStep(self.reconstructStep)
        #self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def wienerCorrectionStep(self, tsIdList):
        #TODO: Check if is phase flipped
        sampling = 1

        for tsId in tsIdList:
            params = '  -i %s' % os.path.join(self._getExtraPath(tsId), tsId+XMD_EXT)
            params += '  -o %s' % os.path.join(self._getExtraPath(tsId), tsId+SUFFIX_CTF_CORR+XMD_EXT)
            params += '  --save_metadata_stack %s' % os.path.join(self._getExtraPath(tsId), tsId+SUFFIX_CTF_CORR+MRC_EXT)
            params += '  --pad %s' % self.padding_factor.get()
            params += '  --wc %s' % self.wiener_constant.get()
            params += '  --sampling_rate %s' % sampling

            if self.inputParticles.get().isPhaseFlipped():
                params += '  --phase_flipped '

            if self.correctEnvelope:
                params += '  --correct_envelope '

            nproc = self.numberOfMpi.get()
            nT = self.numberOfThreads.get()

            self.runJob('xmipp_ctf_correct_wiener2d', params, numberOfMpi=nproc, numberOfThreads=nT)


    def convertInputStep(self):
        # To obtain a list with the tsIds
        sotp = self.inputStacks.get()
        for tp in sotp.iterItems():
            tsId = tp.getTsId()
            fnParticles = os.path.join(self._getExtraPath(tsId), tsId+XMD_EXT)
            self.writeParticleStackToMd(sotp, fnParticles)


    def createGalleryStep(self):

        if self.provideInitialVolume.get():
            fnVol = self.initVol.get()
        else:
            fnVol = ''

        params =  ' -i %s' % fnVol
        params += ' -o %s' % self._getExtraPath('gallery')
        params += ' --sampling_rate %f' % 5
        params += ' ----psi_sampling %f' % 5

        self.runJob('xmipp_angular_project_library', params)


    def writeParticleStackToMd(self, stack, fn):

        md = lib.MetaData()
        psi = 0.0
        for ti in stack:
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

        md.write(fn)

    def reconstructStep(self, objId):

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


