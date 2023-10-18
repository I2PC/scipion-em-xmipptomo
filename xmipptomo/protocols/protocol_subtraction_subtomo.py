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

from pyworkflow import BETA, UPDATED, NEW, PROD
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, FloatParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pwem.protocols import EMProtocol
from pwem.constants import ALIGN_3D
from tomo.protocols import ProtTomoBase
from xmipptomo.convert import writeSetOfSubtomograms, readSetOfSubtomograms

OUTPUT = "output_subtomograms.xmd"
INPUT = "input_subtomograms.xmd"

class XmippProtSubtractionSubtomo(EMProtocol, ProtTomoBase):
    """ This protocol subtracts a subtomogram average to a SetOfSubtomograms, which are internally aligned and
    numerically adjusted in order to obtain reliable results. The adjustment and subtraction is perfomed by
    xmipp_volume_subtraction program. A mask can be provided if the user wants to perform the subtraction in a
    determined region."""

    _label = 'subtomo subtraction'
    _devStatus = UPDATED

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubtomos', PointerParam, pointerClass='SetOfSubTomograms', label="Subtomograms ",
                      help='Select the SetOfSubTomograms with transform matrix which will be subtracted.')
        form.addParam('average', PointerParam, pointerClass='SubTomogram', label="Average subtomogram ",
                      help='Select an average subtomogram to be subtracted.')
        form.addParam('maskBool', BooleanParam, label='Mask subtomograms?', default=True,
                      help='The mask are not mandatory but highly recommendable.')
        form.addParam('mask', PointerParam, pointerClass='VolumeMask', label="Average mask",
                      condition='maskBool', help='Specify a mask for the average.')
        form.addParam('maskSub', PointerParam, pointerClass='VolumeMask', label="Subtraction mask", allowsNull=True,
                      condition='maskBool', help='Optional, specify a mask for the region of subtraction')
        form.addParam('resol', FloatParam, label="Filter at resolution: ", default=3, allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      help='Resolution (A) at which subtraction will be performed, filtering the input volumes.'
                           'Value 0 implies no filtering.')
        form.addParam('sigma', FloatParam, label="Decay of the filter (sigma): ", default=3, condition='resol',
                      help='Decay of the filter (sigma parameter) to smooth the mask transition',
                      expertLevel=LEVEL_ADVANCED)
        form.addParam('iter', IntParam, label="Number of iterations: ", default=5, expertLevel=LEVEL_ADVANCED)
        form.addParam('rfactor', FloatParam, label="Relaxation factor (lambda): ", default=1,
                      expertLevel=LEVEL_ADVANCED,
                      help='Relaxation factor for Fourier amplitude projector (POCS), it should be between 0 and 1, '
                           'being 1 no relaxation and 0 no modification of volume 2 amplitudes')
        form.addParam('saveFiles', BooleanParam, label='Save intermediate files?', default=False,
                      expertLevel=LEVEL_ADVANCED, help='Save input volume 1 (first subtomogram of the set) filtered '
                                                       'and input volume 2 (average) adjusted, which are the volumes '
                                                       'that are really subtracted.')
        form.addParallelSection(threads=0, mpi=4)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertStep')
        self._insertFunctionStep('subtractionStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def convertStep(self):
        writeSetOfSubtomograms(self.inputSubtomos.get(), self._getExtraPath(INPUT), alignType=ALIGN_3D)

    def subtractionStep(self):
        """Subtract reference to each of the subtomogram in the input Set"""
        average = self.average.get().getFileName()
        if average.endswith('.mrc'):
            average += ':mrc'
        resol = self.resol.get()
        iter = self.iter.get()
        program = "xmipp_subtomo_subtraction"
        args = '-i %s -o % s --ref %s --oroot %s --iter %s --lambda %s --sub --subtomos' % \
               (self._getExtraPath(INPUT), self._getExtraPath(OUTPUT), average,
                self._getExtraPath("subtracted_subtomo"), iter, self.rfactor.get())
        if resol:
            fc = self.inputSubtomos.get().getSamplingRate() / resol
            args += ' --cutFreq %f --sigma %d' % (fc, self.sigma.get())
        if self.maskBool:
            args += ' --mask1 %s' % (self.mask.get().getFileName())
            maskSub = self.maskSub.get()
            if maskSub:
                args += ' --maskSub %s' % maskSub.getFileName()
        if self.saveFiles:
            args += ' --saveV1 %s --saveV2 %s' % (self._getExtraPath('v1.mrc'), self._getExtraPath('v2.mrc'))
        self.runJob(program, args)

    def createOutputStep(self):
        inputSubtomos = self.inputSubtomos.get()
        outputSet = self._createSetOfSubTomograms()
        outputSet.copyInfo(inputSubtomos)
        readSetOfSubtomograms(self._getExtraPath(OUTPUT), outputSet, alignType=ALIGN_3D)
        self._defineOutputs(outputSubtomograms=outputSet)
        self._defineSourceRelation(inputSubtomos, outputSet)

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputSubtomograms'):
            summary.append("Output subtomograms not ready yet.")
        else:
            summary.append("Subtraction performed to %s subtomograms." % self.inputSubtomos.get().getSize())
            summary.append("Average subtomogram subtracted: %s" % self.average.get().getFileName())
            if self.maskBool:
                summary.append("Mask: %s" % self.mask.get().getFileName())
            if self.resol.get() != 0:
                summary.append("Subtraction at resolution %0.2f A" % self.resol.get())
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputSubtomograms'):
            methods.append("Output subtomograms not ready yet.")
        else:
            methods.append("Subtraction of average %s performed to %s subtomograms" %
                           (self.average.get().getFileName(), self.inputSubtomos.get().getSize()))
            if self.resol.get() != 0:
                methods.append(" at resolution %0.2f A" % self.resol.get())
        return methods

    def _validate(self):
        validateMsgs = []
        rfactor = self.rfactor.get()
        if rfactor < 0 or rfactor > 1:
            validateMsgs.append('Relaxation factor (lambda) must be between 0 and 1.')
        for subtomo in self.inputSubtomos.get().iterItems():
            if not subtomo.hasTransform():
                validateMsgs.append(
                    'Please provide subtomograms which have transformation matrix.')
        return validateMsgs
