# **************************************************************************
# *
# * Authors:    Federico P. de Isidro-Gomez
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
import numpy as np
import pyworkflow.utils as pwutils
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSet
from tomo.protocols import (ProtImportCoordinates3D,
                            ProtImportTomograms)

from xmipptomo.protocols import XmippProtScoreCoordinates, XmippProtResizeTomograms


class XmippTomoScoreCoordinates(BaseTest):
    """This class check if the protocol project top works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('pyseg')
        cls.tomogram = cls.dataset.getFile('tomo')

    def _runPreviousProtocols(self):
        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=6.87,
                                              objLabel='Import Tomogram')
        self.launchProtocol(protImportTomogram)
        self.assertSetSize(protImportTomogram.Tomograms, size=1, msg=
        "There was a problem with import tomograms output")

        protResizeTomogram = self.newProtocol(XmippProtResizeTomograms,
                                              inputSet=protImportTomogram.Tomograms,
                                              resizeOption=1,
                                              resizeFactor=0.5)
        self.launchProtocol(protResizeTomogram)
        self.assertSetSize(protResizeTomogram.outputSetOfTomograms, size=1, msg=
        "There was a problem with tomogram resizing output")

        # Create phantom coords file
        tomogram_file = protResizeTomogram.outputSetOfTomograms.getFirstItem().getFileName()
        coords_path = protResizeTomogram._getExtraPath(pwutils.removeBaseExt(tomogram_file) + ".txt")
        coords = np.asarray([[172.74201474, -126.69410319, 6.],
                             [116.21130221, -13.63267813, 6.],
                             [64.98034398, 44.66461916, 6.],
                             [24.34889435, -165.55896806, 6.],
                             [3.14987715, 113.56142506, 6.],
                             [-26.88206388, -119.62776413, 6.],
                             [-41.01474201, -174.39189189, 6.],
                             [-55.14742015, 168.32555283, 6.],
                             [-56.91400491, -303.35257985, 6.],
                             [-99.31203931, -57.7972973, 6.],
                             [-125.81081081, 216.02334152, 6.],
                             [-152.30958231, -25.9987715, 6.],
                             [-169.97542998, 87.06265356, 6.],
                             [-171.74201474, 34.06511057, 6.],
                             [-184.10810811, -135.52702703, 6.],
                             [-185.87469287, -234.45577396, 6.]])
        origin = np.asarray([256., 356, 75.])[None, ...]
        np.savetxt(coords_path, coords + origin)

        protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                   objLabel='Import Coordinates - TXT',
                                                   auto=0,
                                                   filesPath=coords_path,
                                                   importTomograms=protResizeTomogram.outputSetOfTomograms,
                                                   filesPattern='', boxSize=64,
                                                   samplingRate=2.0 * 6.87)
        self.launchProtocol(protImportCoordinates3d)
        output = getattr(protImportCoordinates3d, 'outputCoordinates', None)
        self.assertIsNotNone(output,
                             "There was a problem with import coordinates output")

        return protImportCoordinates3d

    def _createProtScoreCarbon(self, protCoords):
        protScoreCoordinates = self.newProtocol(XmippProtScoreCoordinates,
                                                objLabel='Filter Carbon - Threshold 0.8',
                                                inputCoordinates=protCoords.outputCoordinates,
                                                filter=1,
                                                outliers=False,
                                                carbonThreshold=0.5)
        self.launchProtocol(protScoreCoordinates)
        return getattr(protScoreCoordinates, 'outputCoordinates')

    def _createProtScoreOutliers(self, protCoords):
        protScoreCoordinates = self.newProtocol(XmippProtScoreCoordinates,
                                                objLabel='Filter Outliers - Threshold 1',
                                                inputCoordinates=protCoords.outputCoordinates,
                                                filter=1,
                                                carbon=False)
        self.launchProtocol(protScoreCoordinates)
        return getattr(protScoreCoordinates, 'outputCoordinates')

    def test_score_coordinates(self):
        # Run imports
        protCoords = self._runPreviousProtocols()

        # Test Carbon based filtering
        filteredCoords = self._createProtScoreCarbon(protCoords)
        # self.assertTrue(filteredCoords,
        #                 "There was a problem with score coordinates output")
        self.assertTrue(filteredCoords.getSize() < 16)
        self.assertTrue(filteredCoords.getBoxSize() == 64)
        self.assertTrue(filteredCoords.getSamplingRate() == 2.0 * 6.87)

        # Test Outlier based filtering
        filteredCoords = self._createProtScoreOutliers(protCoords)
        # self.assertTrue(filteredCoords,
        #                 "There was a problem with score coordinates output")
        self.assertTrue(filteredCoords.getSize() == 11)
        self.assertTrue(filteredCoords.getBoxSize() == 64)
        self.assertTrue(filteredCoords.getSamplingRate() == 2.0 * 6.87)

# class TestXmippTomoScoreTransform(BaseTest):
#     """This class check if the protocol score_transform works properly."""
#
#     @classmethod
#     def setUpClass(cls):
#         setupTestProject(cls)
#         cls.dataset = DataSet.getDataSet('tomo-em')
#         cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')
#
#     def _runImportSubtomos(self, label):
#         protImport = self.newProtocol(ProtImportSubTomograms,
#                                       filesPath=self.setOfSubtomograms,
#                                       samplingRate=5,
#                                       label=label)
#         self.launchProtocol(protImport)
#         self.assertIsNotNone(protImport.outputSubTomograms,
#                              "There was a problem with import subtomograms output")
#         return protImport.outputSubTomograms
#
#     def _scoreTransformations(self, firstSubtomos, secondSubtomos):
#         scoredTr = self.newProtocol(XmippProtScoreTransform,
#                                    firstSubtomos=firstSubtomos,
#                                    secondSubtomos=secondSubtomos,
#                                    label='Score Subtomogram Transformations')
#         self.launchProtocol(scoredTr)
#         self.assertIsNotNone(scoredTr.outputSetOfSubtomogram,
#                              "There was a problem with score subtomograms output")
#         return scoredTr.outputSetOfSubtomogram
#
#     def test_scoreTransformations(self):
#         firstSubtomos = self._runImportSubtomos('First Subtomograms')
#         secondSubtomos = self._runImportSubtomos('Second Subtomograms')
#         scoredSubtomos = self._scoreTransformations(firstSubtomos, secondSubtomos)
#         self.assertTrue(scoredSubtomos.getSize() == 4)
#         dist = scoredSubtomos.aggregate(["MAX"], "distanceScore", ["distanceScore"])
#         dist = np.asarray([d["distanceScore"] for d in dist])
#         meanDist = np.mean(dist)
#         self.assertTrue(meanDist == 0)
#         return scoredSubtomos
