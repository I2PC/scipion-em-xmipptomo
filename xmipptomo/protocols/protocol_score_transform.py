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
from scipy.spatial.distance import cdist

from pyworkflow import BETA
from pyworkflow.object import Float
from pyworkflow.protocol import params
from tomo.constants import SCIPION

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
        self._percentage_outliers = None
        self.scoresIndex = None
        self._mean_dist = None
        self.first_subtomos = None
        self.second_subtomos = None

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
        # first_matrices, second_matrices = self.getMatricesFromCommonItems()
        commonSet1, commonSet2 = self.getMatchingCoordsFromSubtomos()
        first_matrices, second_matrices = self.getTransformMatrices(commonSet1, commonSet2)

        # Convert Trasnformation Matrices to Quaternions
        first_quaternions = listQuaternions(first_matrices)
        second_quaternions = listQuaternions(second_matrices)

        # Compute distance matrix from quaternions --> dict {subtomo: min_quaternion_distance}
        distDict = {subtomo: min(quaternion_distance(t1, t2), quaternion_distance(t1, -t2))
                    for subtomo, t1, t2 in zip(commonSet1, first_quaternions, second_quaternions)}

        # Crete the output here. Continuation is not possible since "dist" is needed.
        self.createOutput(distDict)

        # Save summary to use it in the protocol Info
        only_distances = list(distDict.values())
        mean_dist = np.mean(only_distances)
        std_dist = np.std(only_distances)
        percentage_outliers = np.sum(only_distances > mean_dist + 3 * std_dist) \
                              + np.sum(only_distances < mean_dist - 3 * std_dist)
        percentage_outliers = 100 * percentage_outliers / len(only_distances)

        self._mean_dist = Float(mean_dist)
        self._percentage_outliers = Float(percentage_outliers)
        self._store()

    def createOutput(self, distanceScoresDict):

        outSubtomos = SetOfSubTomograms.create(self._getPath(), template='submograms%s.sqlite')
        outSubtomos.copyInfo(self.second_subtomos)
        for subtomo, distance in distanceScoresDict.items():
            setattr(subtomo, self.SCORE_ATTR, Float(math.degrees(distance)))
            outSubtomos.append(subtomo)


        # self.scoresIndex = 0
        #
        # def addScoreToSubtomogram(subtomo, row):
        #
        #     try:
        #
        #         score = distanceScores[self.scoresIndex][1]
        #         distance = math.degrees(score)
        #
        #     except Exception as e:
        #         self.info("Can't find score for %s. Adding -181." % subtomo.getObjId())
        #         distance = -181
        #
        #     setattr(subtomo, self.SCORE_ATTR, Float(distance))
        #
        #     self.scoresIndex += 1
        #
        # outSubtomos.copyItems(self.first_subtomos, updateItemCallback=addScoreToSubtomogram)

        self._defineOutputs(**{ScoreTransformOutputs.Subtomograms.name: outSubtomos})
        self._defineSourceRelation(self.firstSubtomos, outSubtomos)
        self._defineSourceRelation(self.secondSubtomos, outSubtomos)

    # --------------------------- UTILS functions ---------------------------
    # def queryMatrices(self, subtomos):
    #     matrices = [(subtomo.getObjId(), subtomo.getTransform().getMatrix())
    #                 for subtomo in subtomos.iterItems()]
    #     return matrices

    def getMatricesFromCommonItems(self):
        objIds1 = [str(item.getObjId()) for item in self.first_subtomos]  # Id list present in the first set
        objIds2 = [str(item.getObjId()) for item in self.second_subtomos]  # The same for the second set
        commonInds = set(objIds1) & set(objIds2)  # Get the Ids present in the intersection of both lists
        whereCond = f'id IN ({",".join(commonInds)})'
        # Get the transformation lists from both sets considering only the common ids and iterated following the same
        # order
        tr1 = [item.getTransform().getMatrix() for item in self.first_subtomos._getMapper().selectAll(where=whereCond)]
        tr2 = [item.getTransform().getMatrix() for item in self.second_subtomos._getMapper().selectAll(where=whereCond)]
        return tr1, tr2

    def getMatchingCoordsFromSubtomos(self, minDistance=1):
        coordDictList1 = {part.clone(): coord.getPosition(SCIPION) for part, coord in
                          zip(self.first_subtomos, self.first_subtomos.getCoordinates3D())}
        coordDictList2 = {part.clone(): coord.getPosition(SCIPION) for part, coord in
                          zip(self.second_subtomos, self.second_subtomos.getCoordinates3D())}
        coordList1 = np.array(list(coordDictList1.values()))
        coordList2 = np.array(list(coordDictList2.values()))
        numel1 = len(coordList1)
        numel2 = len(coordList2)
        if numel1 == numel2 or numel1 < numel2:
            distances = cdist(coordList1, coordList2)
        else:
            distances = cdist(coordList2, coordList1)

        indices = np.argmin(distances, axis=1)  # Indices of min distance
        # m x 2 array, being the min distance indices the first column and the min distance value the second
        minDistances = np.min(distances, axis=1)
        indsAndDistances1 = np.column_stack((range(len(indices)), minDistances))
        indsAndDistances2 = np.column_stack((indices, minDistances))

        # Filter elements using the min distance value provided
        finalIndices1 = np.where(indsAndDistances1[:, 1] < minDistance)
        finalIndices2 = np.where(indsAndDistances2[:, 1] < minDistance)

        # Generate the final lists of coordinates objects
        subtomos1 = list(coordDictList1.keys())
        subtomos2 = list(coordDictList2.keys())
        # [0] in the comprehensions below is because the previous calls to np.where() returns a 1-element tuple
        # containing the desired list of indices
        finalList1 = [subtomos1[ind] for ind in finalIndices1[0]]
        finalList2 = [subtomos2[ind] for ind in finalIndices2[0]]
        return finalList1, finalList2

    def getTransformMatrices(self, subtomos1, subtomos2):
        transforms1 = []
        transforms2 = []
        for part1, part2 in zip(subtomos1, subtomos2):
            transforms1.append(part1.getTransform().getMatrix())
            transforms2.append(part2.getTransform().getMatrix())

        return transforms1, transforms2

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
