# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from pyworkflow.object import Set
from pyworkflow.protocol.params import (PointerParam, FloatParam)
from pwem.protocols import EMProtocol


from xmippBase.protocols.morphology import Morphology

from tomo.objects import TomoMask, SetOfTomograms, SetOfTomoMasks
from pyworkflow import BETA
import pyworkflow.utils.path as path

MRCEXT = '.mrc'
XMDEXT = '.xmd'
SUFFIX_SEGMENTATION_MRC = '_segmentation.mrc'

OUTPUT_SEGMENTATION_NAME = "TomoMasks"
class XmippProtMorphology(EMProtocol):
    """
    This protocol segments a tomogram by means of thresholding, different kind of filters and morphological operations.
    """
    _label = 'tomogram morphology'
    _devStatus = BETA
    _possibleOutputs = {OUTPUT_SEGMENTATION_NAME: SetOfTomoMasks}

    def __init__(self, **args):
        self.morphology = Morphology(self)
        EMProtocol.__init__(self, **args)

        self.TomoMasks  = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSetOfTomograms', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms',
                      important=True,
                      help='Set of tomograms to be segmented or cleaned.')

        #TODO: add wizard
        form.addParam('threshold', FloatParam, default=0.0,
                      label='Threshold',
                      help="Select the threshold. Gray values lesser than the threshold" \
                           "will be set to zero, otherwise will be one (mask area).")
        # Postprocessing
        self.morphology.addPostprocessingSection(form)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):

        for tom in self.inputSetOfTomograms.get():
            tomId = tom.getObjId()
            self._insertFunctionStep(self.morphologyStep, tomId)
            self._insertFunctionStep(self.createOutputStep, tomId)
        self._insertFunctionStep(self.closeOutputSetsStep)

    def morphologyStep(self, tomId):
        inTomograms = self.inputSetOfTomograms.get()

        ts = inTomograms[tomId]
        fnTomo = ts.getFileName()
        tsId = ts.getTsId()

        # Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        path.makePath(tomoPath)

        fnSegmentation = os.path.join(tomoPath, tsId+SUFFIX_SEGMENTATION_MRC)

        self.runJob("xmipp_transform_threshold", "-i %s -o %s --select below %f --substitute binarize"
                    % (fnTomo, fnSegmentation, self.threshold.get()))

        if self.doSmall:
            self.morphology.removeSmallObjects(fnSegmentation, self.smallSize.get())

        if self.doBig:
            self.morphology.keepBiggest(fnSegmentation)

        if self.doMorphological:
            self.morphology.doMorphological(fnSegmentation, self.elementSize.get(), self.getEnumText('morphologicalOperation'))

        if self.doInvert:
            Morphology.doInvert(fnSegmentation)

        if self.doSmooth:
            self.morphology.doSmooth(fnSegmentation, self.sigmaConvolution.get())

    def getOutputSetOfTomoMasks(self, inputSet, setPath, binning=1) -> SetOfTomoMasks:

        if self.TomoMasks:
            getattr(self, OUTPUT_SEGMENTATION_NAME).enableAppend()
        else:
            outputSetOfTomoMasks = SetOfTomoMasks.create(setPath, suffix='tomoMasks')

            if isinstance(inputSet, SetOfTomograms):
                outputSetOfTomoMasks.copyInfo(inputSet)

            if binning > 1:
                samplingRate = inputSet.getSamplingRate()
                samplingRate *= self.binning.get()
                outputSetOfTomoMasks.setSamplingRate(samplingRate)

            outputSetOfTomoMasks.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_SEGMENTATION_NAME: outputSetOfTomoMasks})
            self._defineSourceRelation(inputSet, outputSetOfTomoMasks)

        return self.TomoMasks


    def createOutputStep(self, tomId):
        output = self.getOutputSetOfTomoMasks(self.inputSetOfTomograms.get(), self.getPath())

        tomo = self.inputSetOfTomograms.get()[tomId]
        newTomoMask = TomoMask()
        tsId = tomo.getTsId()

        tomoPath = self._getExtraPath(tsId)
        fnSegmentation = os.path.join(tomoPath, tsId+SUFFIX_SEGMENTATION_MRC)

        newTomoMask.copyInfo(tomo)
        newTomoMask.setTsId(tsId)

        newTomoMask.setLocation(fnSegmentation)
        newTomoMask.setVolName(fnSegmentation)

        newTomoMask.copyAttributes(tomo, '_origin')
        output.append(newTomoMask)
        output.updateDim()
        output.update(newTomoMask)

        output.write()
        self._store()

    def closeOutputSetsStep(self):
        for _, output in self.iterOutputAttributes():
            output.setStreamState(Set.STREAM_CLOSED)
            output.write()
        self._store()

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        messages = []
        if self.Tomograms:
            messages.append("%d tomograms have been segmented .\n"
                           % (self.Tomograms.getSize()))
        return messages

    def _summary(self):
        summary = []
        summary.append("A set of segmentations was created")
        return summary
