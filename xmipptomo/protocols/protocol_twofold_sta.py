# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
# *              Oier Lauzirika Zarrabeitia
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************

from typing import Dict, Tuple, Sequence
import numpy as np
import scipy.linalg

from pyworkflow import BETA
from pyworkflow.protocol import params

from pwem import emlib
from pwem.convert import euler_matrix

from tomo.protocols import ProtTomoSubtomogramAveraging
from tomo.objects import SubTomogram, SetOfSubTomograms


class XmippProtTwofoldSta(ProtTomoSubtomogramAveraging):
    '''Protocol to rotate align a set of subtomograms'''

    _label = 'two-fold sta'
    _devStatus = BETA

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolumes', params.PointerParam,
                      pointerClass=SetOfSubTomograms,
                      label='Input subtomograms', important=True)


    # --------------------------- INSERT steps functions ---------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.alignSubtomogramPairsStep)
        self._insertFunctionStep(self.conciliateAlignmentsStep)
        self._insertFunctionStep(self.averageSubtomogramsStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ---------------------------
    def convertInputStep(self):
        pass
    
    def alignSubtomogramPairsStep(self):
        pass
    
    def conciliateAlignmentsStep(self):
        images = list(map(SubTomogram.getLocation  self.inputVolumes.get())) # TODO
        
        pairwiseAlignment = self._readPairwiseAlignment()
        pairwiseAlignmentMatrix = self._buildPairwiseAlignmentMatrix(pairwiseAlignment, images)
        alignment = self._conciliateAlignmentMatrix(pairwiseAlignmentMatrix)
        
        # TODO write
        
    def averageSubtomogramsStep(self):
        pass
    
    def createOutputStep(self):
        pass


    # --------------------------- UTILS functions ---------------------------
    def _getPairwiseAlignmentMdFilename(self) -> str:
        self._getExtraPath('pairwise.xmd')
    
    def _readPairwiseAlignment(self) -> Dict[Tuple[str, str], np.ndarray]:
        result = dict()
        md = emlib.MetaData(self._getPairwiseAlignmentMdFilename())
        for objId in md:
            # Retrieve from metadata
            image1 = md.getValue(emlib.MDL_IMAGE1, objId)
            image2 = md.getValue(emlib.MDL_IMAGE2, objId)
            rot = md.getValue(emlib.MDL_ANGLE_ROT, objId)
            tilt = md.getValue(emlib.MDL_ANGLE_TILT, objId)
            psi = md.getValue(emlib.MDL_ANGLE_PSI, objId)
            
            # Convert to 3D rotation matrix
            matrix = euler_matrix(
                np.deg2rad(psi), np.deg2rad(tilt), np.deg2rad(rot), 
                'szyz'
            )
            
            # Commit item
            result[(image1, image2)] = matrix[:3, :3] # Only rotation

        return result
    
    def _buildPairwiseAlignmentMatrix(alignments: Dict[Tuple[str, str], np.ndarray], 
                                      images: Sequence[str] ) -> np.ndarray:
        alignmentMatrix = np.eye(3*len(images))
        
        for (image1, image2), alignment in alignments.items():
            # Find the indices where to write the alignment matrix
            index1 = images.index(image1)
            index2 = images.index(image2)
            start1 = 3*index1
            start2 = 3*index2
            end1 = start1 + 3
            end2 = start2 + 3
            
            # Write on symmetric positions. TODO maybe the other way around
            alignmentMatrix[start1:end1,start2:end2] = alignment
            alignmentMatrix[start2:end2,start1:end1] = alignment.T

        return alignmentMatrix

    def _conciliateAlignmentMatrix(matrix: np.ndarray) -> np.ndarray:
        n = len(matrix) // 3
        w, v = scipy.linalg.eigh(matrix, subset_by_index=[n-3, n-1])
        v *= np.sqrt(w)
        return v.reshape((n, 3, 3))

    # --------------------------- INFO functions ---------------------------
    def _summary(self):
        pass