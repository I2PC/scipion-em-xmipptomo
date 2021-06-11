# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
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


import numpy as np

from pyworkflow import BETA
from pyworkflow.object import Float
from pyworkflow.protocol import params

from tomo.protocols import ProtTomoPicking
from tomo.objects import SubTomogram, TomoAcquisition

from pwem.convert.transformations import listQuaternions, quaternion_distance


class XmippProtScoreTransform(ProtTomoPicking):
    '''Protocol to score a series of alignments stored in a SetOfSubtomograms by
    quaternion distance analysis'''

    _label = 'score transformations'
    _devStatus = BETA

    def __init__(self, **args):
        ProtTomoPicking.__init__(self, **args)
        # self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('firstSubtomos', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      label="First Subtomograms to compare", important=True)
        form.addParam('secondSubtomos', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      label="Second Subtomograms to compare", important=True)

    # --------------------------- INSERT steps functions ---------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('scoreTransformStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ---------------------------
    def scoreTransformStep(self):
        self.first_subtomos = self.firstSubtomos.get()
        self.second_subtomos = self.secondSubtomos.get()

        # Extract Transformation Matrices from input SubTomograms
        first_matrices = self.queryMatrices(self.first_subtomos)
        second_matrices = self.queryMatrices(self.second_subtomos)

        # Convert Trasnformation Matrices to Quaternions
        aux = list(zip(*first_matrices))
        first_quaternions = list(zip(aux[0], listQuaternions(aux[1])))
        aux = list(zip(*second_matrices))
        second_quaternions = list(zip(aux[0], listQuaternions(aux[1])))

        # Compute distance matrix from quaternions
        self.dist = [(t1[0], quaternion_distance(t1[1], t2[1]))
                     for t1, t2 in zip(first_quaternions, second_quaternions)
                     if t1[0] == t2[0]]

        # Save summary to use it in the protocol Info
        only_distances = np.asarray(list(zip(*self.dist))[1])
        mean_dist = np.mean(only_distances)
        std_dist = np.std(only_distances)
        percentage_outliers = np.sum(only_distances > mean_dist + 3 * std_dist) \
                              + np.sum(only_distances < mean_dist - 3 * std_dist)
        percentage_outliers = 100 * percentage_outliers / len(only_distances)
        summary = self._getExtraPath('Summary.txt')
        with open(summary, 'w') as file:
            file.write('Mean distance between the two sets: %.2f\n' % mean_dist)
            file.write('Estimated percentage of outliers: %.2f%%\n' % percentage_outliers)

    def createOutputStep(self):
        outSubtomos = self._createSetOfSubTomograms()
        outSubtomos.setSamplingRate(self.second_subtomos.getSamplingRate())
        outSubtomos.setCoordinates3D(self.second_subtomos.getCoordinates3D())
        # acquisition = TomoAcquisition()
        # acquisition.setAngleMin(self.second_subtomos.getFirstItem().getAcquisition().getAngleMin())
        # acquisition.setAngleMax(self.second_subtomos.getFirstItem().getAcquisition().getAngleMax())
        # acquisition.setStep(self.second_subtomos.getFirstItem().getAcquisition().getStep())
        # outSubtomos.setAcquisition(acquisition)
        for ids, inSubtomo in enumerate(self.first_subtomos.iterItems()):
            subtomogram = SubTomogram()
            subtomogram.setObjId(self.dist[ids][0])
            subtomogram.setLocation(inSubtomo.getLocation())
            subtomogram.setCoordinate3D(inSubtomo.getCoordinate3D())
            subtomogram.setTransform(inSubtomo.getTransform())
            subtomogram.setVolName(inSubtomo.getVolName())
            subtomogram.distanceScore = Float(self.dist[ids][1])
            outSubtomos.append(subtomogram)
        self._defineOutputs(outputSetOfSubtomogram=outSubtomos)
        self._defineSourceRelation(self.firstSubtomos, outSubtomos)
        self._defineSourceRelation(self.secondSubtomos, outSubtomos)

    # --------------------------- UTILS functions ---------------------------
    def queryMatrices(self, subtomos):
        matrices = [(subtomo.getObjId(), subtomo.getTransform().getMatrix())
                    for subtomo in subtomos.iterItems()]
        return matrices

    # --------------------------- INFO functions ---------------------------
    def _summary(self):
        summary = []

        if self.getOutputsSize() >= 1:
            summary_file = self._getExtraPath('Summary.txt')
            with open(summary_file, 'r') as file:
                summary.append(file.read())
        else:
            summary.append('Output not ready yet')

        return summary

    def _validate(self):
        validateMsgs = []
        firstTransform = self.firstSubtomos.get().getFirstItem().getTransform()
        secondTransform = self.secondSubtomos.get().getFirstItem().getTransform()

        if firstTransform is None or secondTransform is None:
            validateMsgs.append('Please provide subtomograms which have transformation matrix in "inputAlignment".')

        return validateMsgs
