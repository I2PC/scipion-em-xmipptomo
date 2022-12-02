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

import pwem
from pwem.emlib import MDL_IMAGE
from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
import pwem.emlib.metadata as md

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, BooleanParam

from tomo.objects import AverageSubTomogram, SetOfSubTomograms
from tomo.protocols import ProtTomoBase
from xmipp3.convert import xmippToLocation, writeSetOfVolumes, alignmentToRow

MRC_EXT = '.mrc'
XMD_EXT = '.xmd'
PADDED_SUBTOMOS_BASE_FN = 'padded_subtomograms'
INPUT_SUBTOMOS_BASE_FN = 'input_subtomos'
FN_BASE_STACK = 'input_subtomos_'

FN_OUTPUTSUBTOMOS = 'output_subtomograms.stk'
STA_FN = 'sta.mrc'

class OutputApplyTransform(enum.Enum):
    outputAverage = AverageSubTomogram
    outputSubtomograms = SetOfSubTomograms

class XmippProtApplyTransformSubtomo(EMProtocol, ProtTomoBase):
    """ Apply alignment matrix and produce a new setOfSubtomograms, with each subtomogram aligned to its reference. """

    _label = 'apply alignment subtomo'
    _devStatus = BETA
    _possibleOutputs = OutputApplyTransform

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('inputSubtomograms', PointerParam, pointerClass="SetOfSubTomograms",
                      label='Set of subtomograms', help="Set of subtomograms to be alignment")

        form.addParam('applyPadding', BooleanParam, default=False,
                      label='pad to interpolate?', help="The rotations that take place when the alignment is carried"
                                                        "out can result in a loose of quality. To alleviate it a padding"
                                                        "can be performed, but the price to pay is a high computational burden.")

        form.addParallelSection(threads=1, mpi=4)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        subtomosFn = self._getPath(INPUT_SUBTOMOS_BASE_FN+XMD_EXT)
        self._insertFunctionStep(self.convertInputStep, subtomosFn)
        self._insertFunctionStep(self.applyAlignmentStep, subtomosFn)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self, outputFn):
        inputSet = self.inputSubtomograms.get()
        if inputSet.getFirstItem().getFileName().endswith(MRC_EXT) or \
                inputSet.getFirstItem().getFileName().endswith('.map'):
            S = self._createSetOfSubTomograms()
            S.setSamplingRate(inputSet.getSamplingRate())
            for subtomo in self.inputSubtomograms.get():
                s = subtomo.clone()
                s.setFileName(subtomo.getFileName() + ':mrc')
                S.append(s)
            writeSetOfVolumes(S, outputFn, alignType=pwem.ALIGN_3D)
        else:
            writeSetOfVolumes(inputSet, outputFn, alignType=pwem.ALIGN_3D)
        return [outputFn]

    def estimateNumberOfStacks(self, inputFn):

        inputSt = self.inputSubtomograms.get()

        numberSubtomos = inputSt.getSize()
        bytesOfFloat = 8
        xdim, ydim, zdim = inputSt.getDim()
        sizeOfsubtomo = xdim * ydim * zdim * bytesOfFloat

        # The maximum size of a stack will be 64Gb
        maxStackSize = 64000
        numberSubtomosStack = maxStackSize / sizeOfsubtomo
        numberStacks = round(numberSubtomos/numberSubtomosStack)

        for stk in range(0, numberStacks):

            fn = self._getExtraPath(FN_BASE_STACK + str(stk) + XMD_EXT)
            mdWindow = md.MetaData(fn_mdTomos)
            mdWindowTransform = md.MetaData()
            idList = list(inputSt.getIdSet())
            for row in md.iterRows(mdWindow):
                rowOut = md.Row()
                rowOut.copyFromRow(row)
                id = row.getValue(MDL_IMAGE)
                id = id.split('@')[0]
                id = id.lstrip('0')
                objId = (idList[int(id) - 1])
                alignmentToRow(inputSt[objId].getTransform(), rowOut, pwem.ALIGN_3D)
                rowOut.addToMd(mdWindowTransform)
                fn_mdTomos = fn_mdTomos_orig

            mdWindowTransform.write(fn_mdTomos)


    def applyAlignmentStep(self, inputFn):
        inputSt = self.inputSubtomograms.get()

        fn_mdTomos = inputFn

        # Window subtomograms twice their size
        if self.applyPadding:
            fn_mdTomos_orig = self._getExtraPath("window_with_original_geometry.xmd")
            windowedStk = self._getExtraPath(PADDED_SUBTOMOS_BASE_FN+MRC_EXT)
            self.runJob('xmipp_transform_window', '-i %s -o %s --size %d --save_metadata_stack' %
                        (fn_mdTomos, windowedStk, 2 * inputSt.getFirstItem().getDim()[0]), numberOfMpi=self.numberOfMpi.get())
            fn_mdTomos = self._getExtraPath(PADDED_SUBTOMOS_BASE_FN+XMD_EXT)
            # Add input transform matrix to md generated by xmipp_transform_window
            mdWindow = md.MetaData(fn_mdTomos)
            mdWindowTransform = md.MetaData()
            idList = list(inputSt.getIdSet())
            for row in md.iterRows(mdWindow):
                rowOut = md.Row()
                rowOut.copyFromRow(row)
                id = row.getValue(MDL_IMAGE)
                id = id.split('@')[0]
                id = id.lstrip('0')
                objId = (idList[int(id) - 1])
                alignmentToRow(inputSt[objId].getTransform(), rowOut, pwem.ALIGN_3D)
                rowOut.addToMd(mdWindowTransform)
                fn_mdTomos = fn_mdTomos_orig

            mdWindowTransform.write(fn_mdTomos)

        fn_alignSubtomos = self._getExtraPath('aligned_subtomograms.stk')
        # Align subtomograms
        self.runJob('xmipp_transform_geometry', '-i %s -o %s --apply_transform' %
                    (fn_mdTomos, fn_alignSubtomos), numberOfMpi=self.numberOfMpi.get())

        # Window subtomograms to their original size
        if self.applyPadding:
            outputStk = self._getPath(FN_OUTPUTSUBTOMOS)
            self.runJob('xmipp_transform_window', '-i %s -o %s --size %d ' %
                        (fn_alignSubtomos, outputStk, self.inputSubtomograms.get().getFirstItem().getDim()[0]),
                        numberOfMpi=self.numberOfMpi.get())
            return [outputStk]



    def createOutputStep(self):
        subtomograms = self.inputSubtomograms.get()
        alignedSet = self._createSetOfSubTomograms()
        alignedSet.copyInfo(subtomograms)
        inputMd = self._getPath(FN_OUTPUTSUBTOMOS)
        alignedSet.copyItems(subtomograms,
                             updateItemCallback=self._updateItem,
                             itemDataIterator=md.iterRows(inputMd, sortByLabel=md.MDL_ITEM_ID))
        alignedSet.setAlignment(pwem.ALIGN_NONE)
        avgFile = self._getExtraPath(STA_FN)
        imgh = ImageHandler()
        avgImage = imgh.computeAverage(alignedSet)
        avgImage.write(avgFile)
        avg = AverageSubTomogram()
        avg.setLocation(avgFile)
        avg.copyInfo(alignedSet)

        outputDic = {self._possibleOutputs.outputAverage.name:avg,
                     self._possibleOutputs.outputSubtomograms.name:alignedSet}
        self._defineOutputs(** outputDic)
        self._defineSourceRelation(self.inputSubtomograms, avg)
        self._defineSourceRelation(self.inputSubtomograms, alignedSet)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputSubtomograms'):
            summary.append("Output subtomograms not ready yet.")
        else:
            summary.append("Alignment applied to %s subtomograms." % self.inputSubtomograms.get().getSize())
        return summary

    def _methods(self):
        if not hasattr(self, 'outputSubtomograms'):
            return ["Output subtomograms not ready yet."]
        else:
            return ["We applied alignment to %s subtomograms %s and the output produced is %s."
                    % (self.inputSubtomograms.get().getSize(), self.getObjectTag('inputSubtomograms'),
                       self.getObjectTag('outputSubtomograms'))]

    # --------------------------- UTILS functions --------------------------------------------
    def _updateItem(self, item, row):
        # By default update the item location (index, filename) with the new binary data location
        newFn = row.getValue(md.MDL_IMAGE)
        newLoc = xmippToLocation(newFn)
        item.setLocation(newLoc)
        item.setTransform(None)
