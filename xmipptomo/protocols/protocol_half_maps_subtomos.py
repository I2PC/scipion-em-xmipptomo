# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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

import pwem
from pwem.emlib import MDL_IMAGE
from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
import pwem.emlib.metadata as md

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
import pyworkflow.utils as pwutils

from tomo.objects import AverageSubTomogram, SetOfAverageSubTomograms
from tomo.protocols import ProtTomoBase
from xmipp3.convert import xmippToLocation, writeSetOfVolumes, alignmentToRow


class XmippProtHalfMapsSubtomo(EMProtocol, ProtTomoBase):
    """ Create half maps from a SetOfSubtomograms and its alignment """

    _label = 'half maps'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('inputSubtomograms', PointerParam, pointerClass="SetOfSubTomograms",
                      label='Set of subtomograms', help="Set of subtomograms to be alignment")

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        subtomosEvenFn = self._getPath('input_subtomos_even.xmd')
        subtomosOddFn = self._getPath('input_subtomos_odd.xmd')
        self._insertFunctionStep('convertInputStep', subtomosEvenFn, subtomosOddFn)
        self._insertFunctionStep('applyAlignmentStep', subtomosEvenFn, subtomosOddFn)
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self, evenFn, oddFn):
        inputSet = self.inputSubtomograms.get()
        self.S_even, self.S_odd = self._createSetOfSubTomograms('_even'), self._createSetOfSubTomograms('_odd')
        self.S_even.setSamplingRate(inputSet.getSamplingRate())
        self.S_odd.setSamplingRate(inputSet.getSamplingRate())
        for subtomo in inputSet.iterItems():
            s = subtomo.clone()
            if subtomo.getFileName().endswith('.mrc') or \
                    subtomo.getFileName().endswith('.map'):
                s.setFileName(subtomo.getFileName() + ':mrc')
            self.S_even.append(s) if subtomo.getObjId() % 2 == 0 else self.S_odd.append(s)
        writeSetOfVolumes(self.S_even, evenFn, alignType=pwem.ALIGN_3D)
        writeSetOfVolumes(self.S_odd, oddFn, alignType=pwem.ALIGN_3D)
        return evenFn, oddFn

    def applyAlignmentStep(self, evenFn, oddFn):
        inputSt = self.inputSubtomograms.get()
        # Window subtomograms twice their size
        windowedStk_even = self._getExtraPath('windowed_subtomograms_even.stk')
        windowedStk_odd = self._getExtraPath('windowed_subtomograms_odd.stk')
        self.runJob('xmipp_transform_window', '-i %s -o %s --size %d --save_metadata_stack' %
                    (evenFn, windowedStk_even, 2*inputSt.getFirstItem().getDim()[0]), numberOfMpi=1)
        self.runJob('xmipp_transform_window', '-i %s -o %s --size %d --save_metadata_stack' %
                    (oddFn, windowedStk_odd, 2*inputSt.getFirstItem().getDim()[0]), numberOfMpi=1)

        # Add input transform matrix to md generated by xmipp_transform_window
        outputStks = []
        for path, stack in [(windowedStk_even, 'even'), (windowedStk_odd, 'odd')]:
            mdWindow = md.MetaData(pwutils.replaceExt(path, 'xmd'))
            mdWindowTransform = md.MetaData()
            idList = list(inputSt.getIdSet())
            print(idList)
            for row in md.iterRows(mdWindow):
                rowOut = md.Row()
                rowOut.copyFromRow(row)
                id = row.getValue(MDL_IMAGE)
                id = id.split('@')[0]
                id = id.strip('0')
                if stack == 'even':
                    id = 2 * (int(id) - 1) + 1
                elif stack == 'odd':
                    id = 2 * (int(id) - 1)
                alignmentToRow(inputSt[(idList[id])].getTransform(), rowOut, pwem.ALIGN_3D)
                rowOut.addToMd(mdWindowTransform)
            if "even" in path:
                ori_geo_file = self._getExtraPath("window_with_original_geometry_even.xmd")
                aligned_subtomos_fn = self._getExtraPath("aligned_subtomograms_even.xmd")
                alignStk = self._getExtraPath('aligned_subtomograms_even.stk')
                outputStk = self._getPath('output_subtomograms_even.stk')
            else:
                ori_geo_file = self._getExtraPath("window_with_original_geometry_odd.xmd")
                aligned_subtomos_fn = self._getExtraPath("aligned_subtomograms_odd.xmd")
                alignStk = self._getExtraPath('aligned_subtomograms_odd.stk')
                outputStk = self._getPath('output_subtomograms_odd.stk')
            mdWindowTransform.write(ori_geo_file)
            # Align subtomograms
            self.runJob('xmipp_transform_geometry', '-i %s -o %s --apply_transform' %
                        (ori_geo_file, aligned_subtomos_fn))
            # Window subtomograms to their original size
            self.runJob('xmipp_transform_window', '-i %s -o %s --size %d ' %
                        (alignStk, outputStk, inputSt.getFirstItem().getDim()[0]),
                        numberOfMpi=1)
            outputStks.append(outputStk)
        return outputStks

    def createOutputStep(self):
        setOfAverageSubTomograms = self._createSet(SetOfAverageSubTomograms, 'halfmaps%s.sqlite', "")
        setOfAverageSubTomograms.copyInfo(self.inputSubtomograms.get())
        for path in ['output_subtomograms_even.stk', 'output_subtomograms_odd.stk']:
            subtomograms = self.S_even if "even" in path else self.S_odd
            alignedSet = self._createSetOfSubTomograms()
            alignedSet.copyInfo(subtomograms)
            inputMd = self._getPath(path)
            alignedSet.copyItems(subtomograms,
                                 updateItemCallback=self._updateItem,
                                 itemDataIterator=md.iterRows(inputMd, sortByLabel=md.MDL_ITEM_ID))
            alignedSet.setAlignment(pwem.ALIGN_NONE)
            if "even" in path:
                avgFile = self._getExtraPath("average_even.mrc")
            else:
                avgFile = self._getExtraPath("average_odd.mrc")
            imgh = ImageHandler()
            avgImage = imgh.computeAverage(alignedSet)
            avgImage.write(avgFile)
            avg = AverageSubTomogram()
            avg.setLocation(1, avgFile)
            avg.copyInfo(alignedSet)
            setOfAverageSubTomograms.append(avg)
        self._defineOutputs(halfMaps=setOfAverageSubTomograms)
        self._defineSourceRelation(self.inputSubtomograms, setOfAverageSubTomograms)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'halfMaps'):
            summary.append("Output subtomograms not ready yet.")
        else:
            summary.append("Alignment applied to %s subtomograms." % self.inputSubtomograms.get().getSize())
        return summary

    def _methods(self):
        if not hasattr(self, 'halfMaps'):
            return ["Output subtomograms not ready yet."]
        else:
            return ["We applied alignment to %s subtomograms %s and the output produced is %s."
                    % (self.inputSubtomograms.get().getSize(), self.getObjectTag('inputSubtomograms'),
                       self.getObjectTag('halfMaps'))]

    # --------------------------- UTILS functions --------------------------------------------
    def _updateItem(self, item, row):
        # By default update the item location (index, filename) with the new binary data location
        newFn = row.getValue(md.MDL_IMAGE)
        newLoc = xmippToLocation(newFn)
        item.setLocation(newLoc)
        item.setTransform(None)
        item.setObjId(item.getObjId() - 1 if item.getObjId() == 1 else item.getObjId())
