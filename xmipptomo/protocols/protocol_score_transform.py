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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import enum
import math

import numpy as np

from pyworkflow import BETA
from pyworkflow.object import Float
from pyworkflow.protocol import params

from tomo.protocols import ProtTomoPicking
from tomo.objects import SetOfSubTomograms

from pwem.convert.transformations import listQuaternions, quaternion_distance


class ScoreTransformOutputs(enum.Enum):
    Subtomograms = SetOfSubTomograms

class XmippProtScoreTransform(ProtTomoPicking):
    """Protocol to score a series of alignments stored in a SetOfSubtomograms by
    quaternion distance analysis.

    xmipp_alignmentDistance ranges from 0º to 180º. Therefore, a 0º distance is the best and means alignment is the same.
    The lower the score the more similar is the alignment.
    """

    _label = 'subtomo alignment consensus'
    _devStatus = BETA
    _possibleOutputs = ScoreTransformOutputs
    SCORE_ATTR = "xmipp_alignmentDistance"

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
        # For convenience
        self.first_subtomos = self.firstSubtomos.get()
        self.second_subtomos = self.secondSubtomos.get()

        self._insertFunctionStep(self.scoreTransformStep)

    # --------------------------- STEPS functions ---------------------------
    def scoreTransformStep(self):

        # Extract Transformation Matrices from input SubTomograms
        first_matrices = self.queryMatrices(self.first_subtomos)
        second_matrices = self.queryMatrices(self.second_subtomos)

        # Convert Trasnformation Matrices to Quaternions
        aux = list(zip(*first_matrices))
        first_quaternions = list(zip(aux[0], listQuaternions(aux[1])))
        aux = list(zip(*second_matrices))
        second_quaternions = list(zip(aux[0], listQuaternions(aux[1])))

        # Compute distance matrix from quaternions
        dist = [(t1[0], quaternion_distance(t1[1], t2[1]))
                for t1, t2 in zip(first_quaternions, second_quaternions)
                if t1[0] == t2[0]]

        # Crete the output here. Continuation is not possible since "dist" is needed.
        self.createOutput(dist)

        # Save summary to use it in the protocol Info
        only_distances = np.asarray(list(zip(*dist))[1])
        mean_dist = np.mean(only_distances)
        std_dist = np.std(only_distances)
        percentage_outliers = np.sum(only_distances > mean_dist + 3 * std_dist) \
                              + np.sum(only_distances < mean_dist - 3 * std_dist)
        percentage_outliers = 100 * percentage_outliers / len(only_distances)

        self._mean_dist = Float(mean_dist)
        self._percentage_outliers = Float(percentage_outliers)
        self._store()


    def createOutput(self, distanceScores):

        outSubtomos = self.firstSubtomos.get().create(self._getPath())
        outSubtomos.copyInfo(self.second_subtomos)

        def addScoreToSubtomogram(subtomo, row):

            score = distanceScores[subtomo.getObjId()-1][1]
            distance =math.degrees(score)
            setattr(subtomo, self.SCORE_ATTR, Float(distance))

        outSubtomos.copyItems(self.first_subtomos, updateItemCallback=addScoreToSubtomogram)

        self._defineOutputs(**{ScoreTransformOutputs.Subtomograms.name:outSubtomos})
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

            summary.append('Mean distance between the two sets: %.2f\n' % self._mean_dist)
            summary.append('Estimated percentage of outliers: %.2f%%\n' % self._percentage_outliers)

        else:
            summary.append('Output not ready yet')

        return summary

    def _validate(self):
        validateMsgs = []
        firstTransform = self.firstSubtomos.get().getFirstItem().getTransform()
        secondTransform = self.secondSubtomos.get().getFirstItem().getTransform()

        if firstTransform is None or secondTransform is None:
            validateMsgs.append('Please provide subtomograms which have transformation matrices".')

        return validateMsgs
