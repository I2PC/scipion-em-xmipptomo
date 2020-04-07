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

import numpy as np
import pwem
from pwem.emlib import lib
from pwem.objects import Transform
import pwem.emlib.metadata as md
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam
from xmipp3.convert import alignmentToRow, xmippToLocation, writeSetOfVolumes


class XmippProtUndoAlignSubtomo(EMProtocol):
    """ Apply inverse alignment matrix from one set of subtomograms to an equivalent one, which have been previously
    aligned, in order to undo the aligment. Note that the ids of the subtomograms should match."""

    _label = 'undo alignment subtomo'

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('alignedSubtomograms', PointerParam, pointerClass="SetOfSubTomograms",
                      label='Aligned subtomograms', help="Set of aligned subtomograms")
        form.addParam('matrixSubtomograms', PointerParam, pointerClass="SetOfSubTomograms",
                      label='Subtomograms with transformation', help="Set of subtomograms with transformation matrix")
        form.addParallelSection(threads=0, mpi=8)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        subtomosFn = self._getExtraPath('aligned_input_subtomos.xmd')
        self._insertFunctionStep('convertInputStep', subtomosFn)
        self._insertFunctionStep('applyAlignmentStep', subtomosFn)
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self, outputFn):
        writeSetOfVolumes(self.alignedSubtomograms.get(), outputFn, alignType=pwem.ALIGN_3D)
        return [outputFn]

    def applyAlignmentStep(self, inputFn):
        inputSt = self.alignedSubtomograms.get()
        # Window subtomograms twice their size
        windowedStk = self._getExtraPath('windowed_subtomograms.stk')
        self.runJob('xmipp_transform_window', '-i %s -o %s --size %d --save_metadata_stack' %
                    (inputFn, windowedStk, 2*inputSt.getFirstItem().getDim()[0]), numberOfMpi=1)
        # Add input transform matrix to md generated by xmipp_transform_window
        mdWindow = md.MetaData(self._getExtraPath('windowed_subtomograms.xmd'))
        mdWindowTransform = md.MetaData()
        for row in md.iterRows(mdWindow):
            rowOut = md.Row()
            rowOut.copyFromRow(row)
            id = row.getValue(lib.MDL_IMAGE)
            id = id.split('@')[0]
            id = id.strip('0')
            inverseMatrix = np.linalg.inv(self.matrixSubtomograms.get().__getitem__(id).getTransform().getMatrix())
            transform = Transform()
            transform.setMatrix(inverseMatrix)
            alignmentToRow(transform, rowOut, pwem.ALIGN_3D)
            rowOut.addToMd(mdWindowTransform)
        mdWindowTransform.write(self._getExtraPath("window_with_original_geometry.xmd"))
        # Align subtomograms
        unalignStk = self._getExtraPath('unaligned_subtomograms.stk')
        self.runJob('xmipp_transform_geometry', '-i %s -o %s --apply_transform' %
                    (self._getExtraPath("window_with_original_geometry.xmd"), unalignStk))
        # Window subtomograms to their original size
        outputStk = self._getExtraPath('output_subtomograms.stk')
        self.runJob('xmipp_transform_window', '-i %s -o %s --size %d  --save_metadata_stack' %
                    (unalignStk, outputStk, self.alignedSubtomograms.get().getFirstItem().getDim()[0]),
                    numberOfMpi=1)
        return [outputStk]

    def createOutputStep(self):
        matrixSubtomos = self.matrixSubtomograms.get()
        unAlignedSet = self._createSetOfSubTomograms()
        unAlignedSet.copyInfo(matrixSubtomos)
        outputMd = self._getExtraPath('output_subtomograms.xmd')
        unAlignedSet.copyItems(matrixSubtomos,
                             updateItemCallback=self._updateItem,
                             itemDataIterator=md.iterRows(outputMd, sortByLabel=md.MDL_ITEM_ID))
        self._defineOutputs(outputSubtomograms=unAlignedSet)
        self._defineSourceRelation(self.alignedSubtomograms, unAlignedSet)
        self._defineOutputs(outputSubtomograms=unAlignedSet)
        self._defineSourceRelation(self.matrixSubtomograms, unAlignedSet)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        validateMsgs = []
        for subtomo in self.matrixSubtomograms.get().iterItems():
            if not subtomo.hasTransform():
                validateMsgs.append('Please provide subtomograms which have transformation matrix as "subtomograms with'
                                    ' transformation".')
                break
        for subtomo in self.alignedSubtomograms.get().iterItems():
            if subtomo.hasTransform():
                validateMsgs.append('Please provide subtomograms which have been previously aligned as '
                                    '"aligned subtomograms"')
                break
        return validateMsgs

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputSubtomograms'):
            summary.append("Output subtomograms not ready yet.")
        else:
            summary.append("Alignment undone to %s subtomograms." % self.alignedSubtomograms.get().getSize())
        return summary

    def _methods(self):
        if not hasattr(self, 'outputSubtomograms'):
            return ["Output subtomograms not ready yet."]
        else:
            return ["We undo alignment to %s subtomograms %s and the output produced is %s."
                    % (self.alignedSubtomograms.get().getSize(), self.getObjectTag('alignedSubtomograms'),
                       self.getObjectTag('outputSubtomograms'))]

    # --------------------------- UTILS functions --------------------------------------------
    def _updateItem(self, item, row):
        newFn = row.getValue(md.MDL_IMAGE)
        newLoc = xmippToLocation(newFn)
        item.setLocation(newLoc)
