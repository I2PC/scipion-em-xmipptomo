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

from typing import Dict, Tuple, Sequence, List
import numpy as np
import scipy.linalg
import scipy.spatial

from pyworkflow import BETA
from pyworkflow.protocol import params, LEVEL_ADVANCED

from pwem import emlib
from pwem.convert.transformations import euler_matrix, euler_from_matrix
from pwem.emlib.image import ImageHandler

from tomo.protocols import ProtTomoSubtomogramAveraging
import pwem.objects as objects
import pwem.constants


class XmippProtTwofoldSta(ProtTomoSubtomogramAveraging):
    '''Protocol to rotate align a set of subtomograms'''

    _label = 'two-fold sta'
    _devStatus = BETA

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputVolumes', params.PointerParam,
                      pointerClass=objects.SetOfVolumes,
                      label='Input subtomograms', important=True)
        form.addParam('maxTilt', params.FloatParam, label='Maximum tilt',
                      default=60.0,
                      help='Maximum tilt angle in degrees')
        form.addParam('maxRes', params.FloatParam, label='Alignment resolution',
                      default=10.0, expertLevel=LEVEL_ADVANCED,
                      help='Maximum tilt angle in angstroms')
        form.addParam('angularSampling', params.FloatParam, label='Angular sampling rate',
                      default=5.0, expertLevel=LEVEL_ADVANCED,
                      help='Angular sampling rate in degrees')


    # --------------------------- INSERT steps functions ---------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.alignSubtomogramPairsStep)
        self._insertFunctionStep(self.conciliateAlignmentsStep)
        self._insertFunctionStep(self.averageSubtomogramsStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ---------------------------
    def convertInputStep(self):
        from xmipp3.convert import writeSetOfVolumes
        writeSetOfVolumes(self.inputVolumes.get(), self._getInputVolumesMdFilename())
    
    def alignSubtomogramPairsStep(self):
        maxFreq = self.inputVolumes.get().getSamplingRate() / self.maxRes.get()
        
        cmd = 'xmipp_tomo_volume_align_twofold'
        args = []
        args += ['-i', self._getInputVolumesMdFilename()]
        args += ['-o', self._getPairwiseAlignmentMdFilename()]
        args += ['--maxTilt', self.maxTilt.get()]
        args += ['--maxFreq', maxFreq]
        args += ['--padding', 2.0]
        args += ['--angularSampling', self.angularSampling.get()]
        args += ['--interp', 1]
    
        self.runJob(cmd, args, numberOfMpi=1)
    
    def conciliateAlignmentsStep(self):       
        images = self._getVolumeFilenames()
        
        pairwiseAlignment = self._readPairwiseAlignment()
        pairwiseAlignmentMatrix = self._buildPairwiseAlignmentMatrix(pairwiseAlignment, images)
        
        alignment = self._conciliateAlignmentMatrix(pairwiseAlignmentMatrix, len(images))
        md = self._createAlignmentMetadata(images, alignment)
        md.write(self._getAlignedMdFilename())
        
    def averageSubtomogramsStep(self):
        pass
    
    def createOutputStep(self):
        from xmipp3.convert import readSetOfVolumes
        outputSetOfVolumes = self._createSetOfVolumes()
        outputSetOfVolumes.copyInfo(self.inputVolumes.get())
        readSetOfVolumes(
            self._getAlignedMdFilename(), 
            outputSetOfVolumes, 
            alignType=pwem.constants.ALIGN_3D
        )
        
        self._defineOutputs(volumes = outputSetOfVolumes)
        self._defineSourceRelation(self.inputVolumes, outputSetOfVolumes)

    # --------------------------- UTILS functions ---------------------------
    def _getVolumeFilenames(self) -> List[str]:
        result = []
        for image in self.inputVolumes.get():
            result.append(ImageHandler.locationToXmipp(image.getLocation()))
        return result
    
    def _getInputVolumesMdFilename(self) -> str:
        return self._getExtraPath('volumes.xmd')
        
    def _getPairwiseAlignmentMdFilename(self) -> str:
        return self._getExtraPath('pairwise.xmd')
    
    def _getAlignedMdFilename(self) -> str:
        return self._getExtraPath('aligned.xmd')

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
    
    def _buildPairwiseAlignmentMatrix(self,
                                      alignments: Dict[Tuple[str, str], np.ndarray], 
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

    def _conciliateAlignmentMatrix(self, matrix: np.ndarray, n: int) -> np.ndarray:
        w, v = scipy.linalg.eigh(matrix, subset_by_index=[3*n-3, 3*n-1])
        v *= np.sqrt(w)
        v[:,-1] *= -1 # TODO determine if conditional
        return v.reshape((n, 3, 3))

    def _createAlignmentMetadata(self, images: Sequence[str], alignment: np.ndarray):
        result = emlib.MetaData()

        row = emlib.metadata.Row()
        for i, image in enumerate(images):
            psi, tilt, rot = np.rad2deg(euler_from_matrix(alignment[i], axes='szyz'))
            row.setValue(emlib.MDL_IMAGE, image)
            row.setValue(emlib.MDL_ANGLE_ROT, rot)
            row.setValue(emlib.MDL_ANGLE_TILT, tilt)
            row.setValue(emlib.MDL_ANGLE_PSI, psi)
            
            row.addToMd(result)
            
        return result

    # --------------------------- INFO functions ---------------------------
    def _summary(self):
        pass