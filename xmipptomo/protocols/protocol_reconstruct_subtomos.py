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
import pwem.emlib.metadata as md
from pwem.convert.transformations import euler_matrix
from pwem.objects.data import Transform, Integer
from pwem.protocols import EMProtocol
from pwem.convert.transformations import translation_from_matrix, euler_from_matrix
from pyworkflow import BETA
from pyworkflow.protocol.params import IntParam, FloatParam, EnumParam, PointerParam, BooleanParam
from tomo.protocols import ProtTomoBase
import pwem
from xmipp3.convert import writeSetOfVolumes
from ..utils import calculateRotationAngleAndShiftsFromTM
from ..objects import SetOfTiltSeriesParticle

from pwem.objects import Volume

FN_INPUTPARTICLES = 'ts_'
XMD_EXT = '.xmd'


class XmippProtReconstructSubtomos(EMProtocol, ProtTomoBase):
    """ This protocol performs a smart STA """

    _label = 'reconstruct subtomos'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStacks', PointerParam,
                      pointerClass="SetOfTiltSeriesParticle",
                      label='Tilt series Particle',
                      help="Input Tilt series particles")
        form.addParam('correctCTF', BooleanParam, default=True,
                      label='Correct CTF?',
                      help="The Set of tilt series particles have a CTF associated. "
                           "Set Yes for correcting th CTF. No will not correct the CTF.")

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions -------------------------
    def _insertAllSteps(self):

        self._insertFunctionStep(self.reconstructStep)
        # self._insertFunctionStep(self.reconstructStep)
        # self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def reconstructStep(self):

        stack = self.inputStacks.get()

        s_id = 0
        for s in stack.iterItems():
            print(s)
            idx = s.getObjId()
            tsId = s.getOriginalTs()
            print(tsId)
            folderTsId = self._getExtraPath(str(tsId))
            print(tsId)
            if not os.path.exists(folderTsId):
                os.makedirs(folderTsId)
            fn = os.path.join(folderTsId, 'subtomo_%i.xmd' % idx)
            print(fn)
            self.writeParticleStackToMd(s, fn)
            s_id += 1

            self.reconstructCmd(fn, folderTsId, s_id)

    def writeParticleStackToMd(self, particlestack, fn):

        mtd = lib.MetaData()
        for ti in particlestack.iterItems():
            tilt = ti.getTiltAngle()
            rot, sx, sy = calculateRotationAngleAndShiftsFromTM(ti)

            nRow = md.Row()
            print(ti.getFileName())
            print(ti.getLocation())
            fntp = str(ti.getLocation()[0]+1)+'@'+ti.getLocation()[1]
            nRow.setValue(lib.MDL_IMAGE, fntp)
            if self.correctCTF:
                ctf = ti.getCTF()
                defU = ctf.getDefocusU()
                defV = ctf.getDefocusV()
                defAng = ctf.getDefocusAngle()
                nRow.setValue(lib.MDL_CTF_DEFOCUSU, defU)
                nRow.setValue(lib.MDL_CTF_DEFOCUSV, defV)
                nRow.setValue(lib.MDL_CTF_DEFOCUS_ANGLE, defAng)
            nRow.setValue(lib.MDL_ANGLE_TILT, tilt)
            nRow.setValue(lib.MDL_ANGLE_ROT, rot)
            nRow.setValue(lib.MDL_ANGLE_PSI, 0.0)
            nRow.setValue(lib.MDL_SHIFT_X, 0.0)
            nRow.setValue(lib.MDL_SHIFT_Y, 0.0)
            nRow.addToMd(mtd)

        mtd.write(fn)

    def reconstructCmd(self, fn, outputfolder, idx):

        params = ' -i %s ' % fn
        if self.correctCTF:
            params += ' --useCTF '
        params += ' --sampling %f ' % (self.inputStacks.get().getSamplingRate())
        params += ' -o %s ' % os.path.join(outputfolder,'subtomo_%i.mrc' % idx)
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
