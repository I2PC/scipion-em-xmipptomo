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
from tomo.objects import SetOfSubTomograms, SubTomogram, TomoAcquisition, SetOfTomograms
import tomo.constants as const
from pwem.convert.headers import setMRCSamplingRate
import pwem
from xmipp3.convert import writeSetOfVolumes

from pwem.objects import Volume


FN_INPUTSUBTOMOS = 'input_subtomos.xmd'
STANDARD_STA = 'sta.mrc'
SMART_STA = 'STA_smart.mrc'


class XmippProtsmartSTA(EMProtocol, ProtTomoBase):
    """ This protocol performs a smart STA """

    _label = 'smart STA'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('subtomos', PointerParam,
                      pointerClass="SetOfSubTomograms",
                      label='Input subtomograms',
                      help="Input subtomograms aligned")
        form.addParam('hasReference', BooleanParam, default=False,
                      label='Reference Map',
                      help="Reference Map")
        form.addParam('reference', PointerParam, condition='hasReference',
                      pointerClass="Volume",
                      label='Input subtomograms',
                      help="Input subtomograms aligned")
        form.addParam('fscResolution', FloatParam,
                      label='fscResolution',
                      help="FSC resolution")

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.smartSTAStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):

        self.subtomoIn = self._createSetOfSubTomograms('input')
        self.subtomoIn.setSamplingRate(self.subtomos.get().getSamplingRate())
        subtomosFn = self._getExtraPath(FN_INPUTSUBTOMOS)
        for subtomo in self.subtomos.get().iterItems():
            s = subtomo.clone()
            if subtomo.getFileName().endswith('.mrc') or subtomo.getFileName().endswith('.map'):
                s.setFileName(subtomo.getFileName() + ':mrc')
            self.subtomoIn.append(s)
        writeSetOfVolumes(self.subtomoIn, subtomosFn, alignType=pwem.ALIGN_3D)

    def smartSTAStep(self):

        nthreads = 1

        params_phantom = ' --subtomos %s ' % self._getExtraPath(FN_INPUTSUBTOMOS)
        if self.hasReference.get():
            params_phantom += ' --reference %s ' % self.reference.get().getFileName()
        params_phantom += ' --fscResolution %s ' % (self.fscResolution.get())
        params_phantom += ' --sampling %f ' % (self.subtomos.get().getSamplingRate())
        params_phantom += ' -o %s ' % self._getExtraPath()
        params_phantom += ' --threads %d ' % nthreads
        self.runJob('xmipp_tomo_sta_deblurring', params_phantom)

    def createOutputStep(self):
        volume = Volume()
        volume.setFileName(self._getExtraPath(STANDARD_STA))
        volume.setSamplingRate(self.subtomos.get().getSamplingRate())
        self._defineOutputs(standard_sta=volume)
        self._defineSourceRelation(self.subtomos, volume)

        volume.setFileName(self._getExtraPath(SMART_STA))
        volume.setSamplingRate(self.subtomos.get().getSamplingRate())
        self._defineOutputs(smart_sta=volume)
        self._defineSourceRelation(self.subtomos, volume)


        # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputSubtomograms'):
            summary.append("Output phantom not ready yet.")
        else:
            summary.append("%s phantoms created with random orientations" % 50)
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

