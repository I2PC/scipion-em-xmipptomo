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

from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, FloatParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils import replaceBaseExt
from pwem.convert import headers
from pwem.objects import Volume, Transform
from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase


class XmippProtSubtractionSubtomo(EMProtocol, ProtTomoBase):
    """ This protocol subtracts a subtomogram average to a SetOfSubtomograms, which are internally aligned and
    numerically adjusted in order to obtain reliable results. The adjustment and subtraction is perfomed by
    xmipp_volume_subtraction program. A mask can be provided if the user wants to perform the subtraction in a
    determined region."""

    _label = 'subtomo subtraction'

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubtomos', PointerParam, pointerClass='SubTomogram', label="Subtomograms ",
                      help='Select the SetOfSubTomograms with transform matrix which will be subtracted.')
        form.addParam('average', PointerParam, pointerClass='SubTomogram', label="Average subtomogram ",
                      help='Select an average subtomogram to be subtracted.')
        form.addParam('maskBool', BooleanParam, label='Mask subtomograms?', default=True,
                      help='The mask are not mandatory but highly recommendable.')
        form.addParam('mask', PointerParam, pointerClass='VolumeMask', label="Mask",
                      condition='maskBool', help='Specify a mask for the region of subtraction.')
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
                      expertLevel=LEVEL_ADVANCED, help='Save input volume 1 filtered and input volume 2 adjusted, which'
                                                       'are the volumes that are really subtracted.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('subtractionStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def subtractionStep(self):
        inputSubtomos = self.inputSubtomos.get().clone()
        fileName1 = inputSubtomos.getFileName()
        if fileName1.endswith('.mrc'):
            fileName1 += ':mrc'

        average = self.average.get().getFileName()
        if average.endswith('.mrc'):
            average += ':mrc'

        resol = self.resol.get()
        iter = self.iter.get()
        program = "xmipp_volume_subtraction"
        args = '--i1 %s --i2 %s -o %s --iter %s --lambda %s --sub' % \
               (inputSubtomos.getFileName(), average, self._getExtraPath("output_volume.mrc"), iter, self.rfactor.get())
        if resol:
            fc = inputSubtomos.getSamplingRate() / resol
            args += ' --cutFreq %f --sigma %d' % (fc, self.sigma.get())

        if self.maskBool:
            args += ' --mask1 %s' % (self.mask.get().getFileName())

        if self.saveFiles:
            args += ' --saveV1 %s --saveV2 %s' % (self._getExtraPath('vol1_filtered.mrc'),
                                                  self._getExtraPath('average_adjusted.mrc'))
        self.runJob(program, args)

    def createOutputStep(self):
        inputSubtomos = self.inputSubtomos.get()
        volume = Volume()
        volume.setSamplingRate(inputSubtomos.getSamplingRate())
        if inputSubtomos.getFileName().endswith('mrc'):
            origin = Transform()
            ccp4header = headers.Ccp4Header(inputSubtomos.getFileName(), readHeader=True)
            shifts = ccp4header.getOrigin()
            origin.setShiftsTuple(shifts)
            volume.setOrigin(origin)
        volume.setFileName(self._getExtraPath("output_volume.mrc"))
        filename = volume.getFileName()
        if filename.endswith('.mrc') or filename.endswith('.map'):
            volume.setFileName(filename + ':mrc')
        self._defineOutputs(outputVolume=volume)

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = ["Volume 1: %s" % self.inputSubtomos.get().getFileName()]

        summary.append("Volume 2: %s" % self.average.get().getFileName())
        if self.maskBool:
            summary.append("Input mask 1: %s" % self.mask.get().getFileName())
        if self.resol.get() != 0:
            summary.append("Subtraction at resolution %f A" % self.resol.get())
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputVolume'):
            methods.append("Output volume not ready yet.")
        else:
            methods.append("Volume %s subtracted from volume %s" % (self.average.get().getFileName(),
                                                                    self.inputSubtomos.get().getFileName()))
            if self.resol.get() != 0:
                methods.append(" at resolution %f A" % self.resol.get())
        return methods

    def _validate(self):
        errors = []
        rfactor = self.rfactor.get()
        if rfactor < 0 or rfactor > 1:
            errors.append('Relaxation factor (lambda) must be between 0 and 1')
        return errors

    # --------------------------- UTLIS functions --------------------------------------------
    def _getPdbFileName(self):
        if self.inputPdbData == self.IMPORT_OBJ:
            return self.pdbObj.get().getFileName()
        else:
            return self.pdbFile.get()

    def _getVolName(self):
        return self._getExtraPath(replaceBaseExt(self._getPdbFileName(), "vol"))
