# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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
import numpy as np

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
import pyworkflow.utils as pwutils
from pwem.viewers.showj import *
from pwem.protocols import ProtAnalysis2D

from tomo.objects import MeshPoint
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider, TomogramsDialog
from tomo.protocols.protocol_base import ProtTomoBase
import tomo.constants as const


class XmippProtRoiIJ(ProtAnalysis2D, ProtTomoBase):
    """ Tomogram ROI selection in IJ """

    _label = 'imagej roi'
    _devStatus = BETA

    def __init__(self, **kwargs):
        ProtAnalysis2D.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='General parameters')
        form.addParam('inputTomos', PointerParam, label="Input Tomograms",
                      pointerClass='SetOfTomograms',
                      help='Select tomograms.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # Launch Boxing GUI
        self._insertFunctionStep('launchIJGUIStep', interactive=True)

    def _createOutput(self):
        inputTomos = self.inputTomos.get()
        outSet = self._createSetOfMeshes(inputTomos)
        for tomo in inputTomos.iterItems():
            tomoName = pwutils.removeBaseExt(tomo.getFileName())
            outFile = self._getExtraPath(tomoName + '.txt')
            if os.path.isfile(outFile):
                data = np.loadtxt(outFile, delimiter=',')
                for coord in data:
                    mesh = MeshPoint()
                    mesh.setVolume(tomo)
                    mesh.setPosition(coord[0], coord[1], coord[2], const.BOTTOM_LEFT_CORNER)
                    mesh.setGroupId(coord[3])
                    outSet.append(mesh)
        outSet.setPrecedents(inputTomos)
        self._defineOutputs(outputMeshes=outSet)
        self._defineSourceRelation(inputTomos, outSet)

    # --------------------------- STEPS functions -----------------------------
    def launchIJGUIStep(self):

        tomoList = []
        for tomo in self.inputTomos.get().iterItems():
            tomogram = tomo.clone()
            tomoName = pwutils.removeBaseExt(tomo.getFileName())
            outFile = self._getExtraPath(tomoName + '.txt')
            if os.path.isfile(outFile):
                tomogram.count = np.loadtxt(outFile, delimiter=',').shape[0]
            else:
                tomogram.count = 0
            tomoList.append(tomogram)

        tomoProvider = TomogramsTreeProvider(tomoList, self._getExtraPath(), 'txt')

        path = self._getExtraPath()
        self.dlg = TomogramsDialog(None, False, provider=tomoProvider, path=path)

        self._createOutput()


    def _summary(self):
        summary = []
        if not os.listdir(self._getExtraPath()):
            summary.append("Output ROIs not ready yet.")
        else:
            count = 0
            for file in os.listdir(self._getExtraPath()):
                if file.endswith(".txt"):
                    count += 1
            summary.append("Rois defined for %d/%d files have been saved in Scipion (%s)." % (
                count, self.inputTomos.get().getSize(), self._getExtraPath()))
        return summary

    def _methods(self):
        tomos = self.inputTomos.get()
        return [
            "ROI selection and extraction using ImageJ",
            "A total of %d tomograms of dimensions %s were used"
            % (tomos.getSize(), tomos.getDimensions()),
        ]
