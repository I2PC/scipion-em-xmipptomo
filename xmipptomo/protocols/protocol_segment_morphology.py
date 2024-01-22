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

from pwem.protocols import ProtAnalysis3D
from ..xTomoIO import xTomoIO

from tomo.objects import Tomogram
from pyworkflow import BETA
from xmipp3.scipionSharedUtils.morphology import Morphology

MRCEXT = '.mrc'
XMDEXT = '.xmd'


class XTomoProtMorphology(xTomoIO):
    """
    This protocol segments a tomogram by means of thresholding, different kind of filters and morphological operations.
    """
    _label = 'tomogram morphology'
    _devStatus = BETA

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inTomos', PointerParam,
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
        Morphology.addPostprocessingSection(form)

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):

        for tom in self.inTomos.get():
            tomId = tom.getObjId()
            self._insertFunctionStep(self.morphologyStep, tomId)
        self._insertFunctionStep(self.createOutputStep(tomId))

    def morphologyStep(self, tomId):

        inTomograms = self.inTomos.get()

        ts = inTomograms[tomId]
        fnTomo = ts.getFileName()
        tsId = ts.getTsId()

        # Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        if not os.path.exists(tomoPath):
            os.mkdir(tomoPath)
        fnSegmentation = os.path.join(tomoPath, tsId+'segmented.mrc')

        self.runJob("xmipp_transform_threshold", "-i %s -o %s --select below %f --substitute binarize"
                    % (fnTomo, fnSegmentation, self.threshold.get()))

        if self.doSmall:
            Morphology.removeSmallObjects(fnSegmentation, self.smallSize.get())

        if self.doBig:
            Morphology.keepBiggest(fnSegmentation)

        if self.doMorphological:
            Morphology.doMorphological(fnSegmentation, self.elementSize.get(), self.getEnumText('morphologicalOperation'))

        if self.doInvert:
            Morphology.doInvert(fnSegmentation)

        if self.doSmooth:
            Morphology.doSmooth(fnSegmentation, self.sigmaConvolution.get())

        outputSegmentation = self.getOutputSetOfTomograms(inTomograms, inTomograms.getSamplingRate())

        newSegmentation = Tomogram()
        newSegmentation.copyInfo(ts)
        newSegmentation.copyAttributes(ts, '_origin')

        newSegmentation.setLocation(fnSegmentation)

        newSegmentation.setSamplingRate(ts.getSamplingRate())
        outputSegmentation.append(newSegmentation)
        outputSegmentation.update(newSegmentation)
        outputSegmentation.write()
        self._store()

    def createOutputStep(self):
        '''
        This function closes the generated output of the protocol
        '''
        inTomograms = self.inTomos.get()
        self.getOutputSetOfTomograms(inTomograms, inTomograms.getSamplingRate()).setStreamState(Set.STREAM_CLOSED)
        self._store()

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        messages = []
        if hasattr(self, 'LocalResolution_Tomogram'):
            messages.append('A set of segmentation was created')
        return messages

    def _summary(self):
        summary = []
        summary.append("A set of segmentation was created")
        return summary

