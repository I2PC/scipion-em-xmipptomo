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

from pwem import ALIGN_3D
from pwem.convert import headers
from pwem.emlib import MDL_IMAGE
import pwem.emlib.metadata as md
from pwem.objects import Volume, Transform
from pwem.protocols import EMProtocol

from tomo.objects import SubTomogram
from tomo.protocols import ProtTomoBase
from xmipp3.convert import xmippToLocation, writeSetOfVolumes, alignmentToRow


class XmippProtSubtractionSubtomo(EMProtocol, ProtTomoBase):
    """ This protocol subtracts a subtomogram average to a SetOfSubtomograms, which are internally aligned and
    numerically adjusted in order to obtain reliable results. The adjustment and subtraction is perfomed by
    xmipp_volume_subtraction program. A mask can be provided if the user wants to perform the subtraction in a
    determined region."""

    _label = 'subtomo subtraction'

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSubtomos', PointerParam, pointerClass='SetOfSubTomograms', label="Subtomograms ",
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
                      expertLevel=LEVEL_ADVANCED, help='Save input volume 1 (first subtomogram of the set) filtered '
                                                       'and input volume 2 (average) adjusted, which are the volumes '
                                                       'that are really subtracted.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('applyAlignStep')
        self._insertFunctionStep('subtractionStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def applyAlignStep(self):
        inputSet = self.inputSubtomos.get()
        outputFn = self._getExtraPath('input_subtomos.xmd')
        if inputSet.getFirstItem().getFileName().endswith('.mrc') or \
                inputSet.getFirstItem().getFileName().endswith('.map'):
            S = self._createSetOfSubTomograms()
            S.setSamplingRate(inputSet.getSamplingRate())
            for subtomo in inputSet:
                s = subtomo.clone()
                s.setFileName(subtomo.getFileName() + ':mrc')
                S.append(s)
            writeSetOfVolumes(S, outputFn, alignType=ALIGN_3D)
        else:
            writeSetOfVolumes(inputSet, outputFn, alignType=ALIGN_3D)

        # Window subtomograms twice their size
        windowedStk = self._getExtraPath('windowed_subtomograms.stk')
        self.runJob('xmipp_transform_window', '-i %s -o %s --size %d --save_metadata_stack' %
                    (outputFn, windowedStk, 2 * inputSet.getFirstItem().getDim()[0]), numberOfMpi=1)

        # Add input transform matrix to md generated by xmipp_transform_window
        mdWindow = md.MetaData(self._getExtraPath('windowed_subtomograms.xmd'))
        mdWindowTransform = md.MetaData()
        idList = list(inputSet.getIdSet())
        for row in md.iterRows(mdWindow):
            rowOut = md.Row()
            rowOut.copyFromRow(row)
            id = row.getValue(MDL_IMAGE)
            id = id.split('@')[0]
            id = id.strip('0')
            alignmentToRow(inputSet[(idList[int(id) - 1])].getTransform(), rowOut, ALIGN_3D)
            rowOut.addToMd(mdWindowTransform)
        mdWindowTransform.write(self._getExtraPath("window_with_original_geometry.xmd"))
        # Align subtomograms
        self.runJob('xmipp_transform_geometry', '-i %s -o %s --apply_transform' %
                    (self._getExtraPath("window_with_original_geometry.xmd"),
                     self._getExtraPath('aligned_subtomograms.stk')))
        # Window subtomograms to their original size
        alignStk = self._getExtraPath('aligned_subtomograms.stk')
        outputStk = self._getExtraPath('output_subtomograms.stk')
        self.runJob('xmipp_transform_window', '-i %s -o %s --size %d ' %
                    (alignStk, outputStk, inputSet.getFirstItem().getDim()[0]),
                    numberOfMpi=1)

        self.alignedSet = self._createSetOfSubTomograms()
        self.alignedSet.copyInfo(inputSet)
        inputMd = self._getExtraPath('output_subtomograms.stk')
        self.alignedSet.copyItems(inputSet,
                             updateItemCallback=self._updateItemAlign,
                             itemDataIterator=md.iterRows(inputMd, sortByLabel=md.MDL_ITEM_ID))


    def subtractionStep(self):
        average = self.average.get().getFileName()
        if average.endswith('.mrc'):
            average += ':mrc'
        resol = self.resol.get()
        iter = self.iter.get()
        program = "xmipp_volume_subtraction"
        args = '--i2 % s --iter %s --lambda %s --sub' % (average, iter, self.rfactor.get())
        if resol:
            fc = self.inputSubtomos.get().getSamplingRate() / resol
            args += ' --cutFreq %f --sigma %d' % (fc, self.sigma.get())
        if self.maskBool:
            args += ' --mask1 %s' % (self.mask.get().getFileName())

        for subtomo in self.alignedSet.iterItems():
            ix = subtomo.getObjId()
            argsSubtomo = ' --i1 %06d@%s -o % s' % (ix, subtomo.getFileName(),
                                                    self._getExtraPath("output_subtomo%d.mrc" % ix))
            if self.saveFiles and ix == 1:
                argsSubtomo += ' --saveV1 %s --saveV2 %s' % (self._getExtraPath('vol1_filtered.mrc'),
                                                      self._getExtraPath('average_adjusted.mrc'))
            print('\n-----Subtomogram %d-----' % ix)
            self.runJob(program, args + argsSubtomo)
            # filename = subtomo.getFileName()
            # subtomo.setFileName(filename + ':mrc')
            # origin = Transform()
            # ccp4header = headers.Ccp4Header(subtomo.getFileName(), readHeader=True)
            # shifts = ccp4header.getOrigin()
            # origin.setShiftsTuple(shifts)
            # subtomo.setOrigin(origin)

            # APPLY INVERSE TRANSF


    def createOutputStep(self):
        inputSubtomos = self.inputSubtomos.get()
        outputSet = self._createSetOfSubTomograms()
        outputSet.copyItems(inputSubtomos, updateItemCallback=self._updateItemOutput)
        # volume = Volume()
        # volume.setSamplingRate(inputSubtomos.getSamplingRate())
        # if inputSubtomos.getFileName().endswith('mrc'):
        #     origin = Transform()
        #     ccp4header = headers.Ccp4Header(inputSubtomos.getFileName(), readHeader=True)
        #     shifts = ccp4header.getOrigin()
        #     origin.setShiftsTuple(shifts)
        #     volume.setOrigin(origin)
        # volume.setFileName(self._getExtraPath("output_volume.mrc"))
        # filename = volume.getFileName()
        # if filename.endswith('.mrc') or filename.endswith('.map'):
        #     volume.setFileName(filename + ':mrc')
        self._defineOutputs(outputSet=outputSet)
        self._defineSourceRelation(inputSubtomos, outputSet)

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
        validateMsgs = []
        rfactor = self.rfactor.get()
        if rfactor < 0 or rfactor > 1:
            validateMsgs.append('Relaxation factor (lambda) must be between 0 and 1')
        for subtomo in self.inputSubtomos.get().iterItems():
            if not subtomo.hasTransform():
                validateMsgs.append(
                    'Please provide subtomograms which have transformation matrix in "inputAlignment".')
        return validateMsgs

    # --------------------------- UTLIS functions --------------------------------------------
    def _getPdbFileName(self):
        if self.inputPdbData == self.IMPORT_OBJ:
            return self.pdbObj.get().getFileName()
        else:
            return self.pdbFile.get()

    def _getVolName(self):
        return self._getExtraPath(replaceBaseExt(self._getPdbFileName(), "vol"))

    def _updateItemAlign(self, item, row):
        newFn = row.getValue(md.MDL_IMAGE)
        newLoc = xmippToLocation(newFn)
        item.setLocation(newLoc)
        item.setTransform(None)

    def _updateItemOutput(self, item, row):
        item.setFileName(self._getExtraPath("output_subtomo%d.mrc" % item.getObjId()))

