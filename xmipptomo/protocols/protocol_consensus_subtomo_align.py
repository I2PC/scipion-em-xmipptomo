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
from pyworkflow.protocol import params, LEVEL_ADVANCED
from tomo.constants import SCIPION

from tomo.protocols import ProtTomoPicking
from tomo.objects import SetOfSubTomograms

from pwem.convert.transformations import listQuaternions, quaternion_distance


class ScoreTransformOutputs(enum.Enum):
    Subtomograms = SetOfSubTomograms


class XmippProtSubtomoAlignConsensus(ProtTomoPicking):
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
        form.addParam('minDistance', params.IntParam,
                      default=0,
                      label='Min distance [Å]',
                      expertLevel=LEVEL_ADVANCED,
                      help='Minimum distance to be considered the same particle [Å]')
        form.addParam('maxDistance', params.IntParam,
                      default=10,
                      label='Max distance [Å]',
                      expertLevel=LEVEL_ADVANCED,
                      help='Maximum distance to be considered different particles [Å]')

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
        first_matrices, second_matrices = self._getTransformMatrices(commonSet1, commonSet2)

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

        self._defineOutputs(**{ScoreTransformOutputs.Subtomograms.name: outSubtomos})
        self._defineSourceRelation(self.firstSubtomos, outSubtomos)
        self._defineSourceRelation(self.secondSubtomos, outSubtomos)

    # --------------------------- UTILS functions ---------------------------
    def getMatchingCoordsFromSubtomos(self):
        coordDictList1 = {part.clone(): coord.getPosition(SCIPION) for part, coord in
                          zip(self.first_subtomos, self.first_subtomos.getCoordinates3D())}
        coordDictList2 = {part.clone(): coord.getPosition(SCIPION) for part, coord in
                          zip(self.second_subtomos, self.second_subtomos.getCoordinates3D())}
        subtomos1 = list(coordDictList1.keys())
        subtomos2 = list(coordDictList2.keys())
        coordList1 = np.array(list(coordDictList1.values()))
        coordList2 = np.array(list(coordDictList2.values()))
        # Coords are in pixels, so the threshold must be, too
        sRate = self.first_subtomos.getSamplingRate()
        minDistThresholdPix = self.minDistance.get() * sRate
        maxDistThresholdPix = self.maxDistance.get() * sRate
        numel1 = len(coordList1)
        numel2 = len(coordList2)
        if numel1 == numel2 or numel1 < numel2:
            distances = cdist(coordList1, coordList2)
        else:
            distances = cdist(coordList2, coordList1)
        indices = np.argmin(distances, axis=1)  # Indices of min distance
        minDist = np.min(distances, axis=1)  # Min distance values
        # Apply threshold to the results
        filterInd = np.logical_and(minDistThresholdPix <= minDist, minDist <= maxDistThresholdPix)
        finalIndices = indices[filterInd]
        # Index the introduced sets properly
        finalList1 = np.array(subtomos1)[finalIndices]
        finalList2 = np.array(subtomos2)[finalIndices]
        return finalList1, finalList2

    @staticmethod
    def _getTransformMatrices(subtomos1, subtomos2):
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
        firstSubtomos = self.firstSubtomos.get()
        secondSubtomos = self.secondSubtomos.get()
        firstSRate = firstSubtomos.getSamplingRate()
        secondSRate = secondSubtomos.getSamplingRate()
        sRateTol = 1e-4
        if not firstSubtomos.hasAlignment3D():
            validateMsgs.append('The first set of subtomograms provided does not contain aligned particles.')
        if not secondSubtomos.hasAlignment3D():
            validateMsgs.append('The second set of subtomograms provided does not contain aligned particles.')
        if abs(firstSRate - secondSRate) >= sRateTol:
            validateMsgs.append('The sampling rate of both sets of subtomograms is expected to be the same within '
                                'tolerance %.4f: %.4f != %.4f' % (sRateTol, firstSRate, secondSRate))
        return validateMsgs
