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
from pyworkflow.protocol.params import PointerParam
import pyworkflow.utils.path as path

from pwem.protocols import EMProtocol
from tomo.objects import Tomogram, SetOfTomograms
from pyworkflow import BETA


MRCEXT = '.mrc'
XMDEXT = '.xmd'
OUTPUT_TOMOGRAMS_NAME = "Tomograms"

class XmippProtApplySegmentation(EMProtocol):
    """
    This protocol applies a segmentation to a set of tomograms.
    """
    _label = 'apply segmentation'
    _devStatus = BETA
    _OUTPUT_NAME = OUTPUT_TOMOGRAMS_NAME
    _possibleOutputs = {OUTPUT_TOMOGRAMS_NAME: SetOfTomograms}

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSetOfTomograms', PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Tomograms',
                      important=True,
                      help='The segmentations or masks will be applied to this set of tomograms.')

        form.addParam('inputSetOfTomoMasks', PointerParam,
                      pointerClass='SetOfTomoMasks',
                      label='Segmentations',
                      important=True,
                      help='Set of segmentation to be applied to the set of tomograms.')

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):

        for tom in self.inputSetOfTomograms.get():
            tomId = tom.getObjId()
            self._insertFunctionStep(self.applySegmentationStep, tomId)
            self._insertFunctionStep(self.createOutputStep, tomId)
        self._insertFunctionStep(self.closeOutputSetsStep)

    def applySegmentationStep(self, tomId):

        inTomograms = self.inputSetOfTomograms.get()
        ts = inTomograms[tomId]
        fnTomo = ts.getFileName()
        tsId = ts.getTsId()

        segmentations = self.inputSetOfTomoMasks.get()

        for seg in segmentations.iterItems():
            tsidSeg = seg.getTsId()

            if tsidSeg == tsId:
                # Defining the output folder
                tomoPath = self._getExtraPath(tsId)
                path.makePath(tomoPath)
                fnSegmentation = seg.getFileName()
                fnMasked = self.createOutputExtraPath(tsId, '_masked'+MRCEXT)
                self.runJob("xmipp_image_operate", "-i %s --mult %s -o %s" % (fnTomo, fnSegmentation, fnMasked))

    def createOutputStep(self, tsObjId):
        ts = self.inputSetOfTomograms.get()[tsObjId]
        tsId = ts.getTsId()

        fullTomogramName = self.createOutputExtraPath(tsId, '_masked'+MRCEXT)

        if os.path.exists(fullTomogramName):
            output = self.getOutputSetOfTomograms(self.inputSetOfTomograms.get())

            newTomogram = Tomogram()
            newTomogram.setLocation(fullTomogramName)

            newTomogram.setTsId(tsId)
            newTomogram.setSamplingRate(ts.getSamplingRate())

            # Set default tomogram origin
            newTomogram.setOrigin(newOrigin=None)
            newTomogram.setAcquisition(ts.getAcquisition())

            output.append(newTomogram)
            output.update(newTomogram)
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
        messages.append('A set of segmentations were applied to the tomograms')
        return messages

    def _summary(self):
        summary = []
        summary.append("Segmentations were applied to the tomograms")
        return summary
