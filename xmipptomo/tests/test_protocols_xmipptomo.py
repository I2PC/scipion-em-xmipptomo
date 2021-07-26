# **************************************************************************
# *
# * Authors:    Estrella Fernandez Gimenez [1]
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

# import numpy as np

from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSet
from tomo.protocols import (ProtImportCoordinates3D,
                            ProtImportTomograms,
                            ProtImportSubTomograms)

from xmipptomo.protocols import XmippProtSubtomoProject, XmippProtConnectedComponents, XmippProtApplyTransformSubtomo, \
    XmippProtSubtomoMapBack, XmippProtPhantomSubtomo, XmippProtScoreCoordinates, XmippProtScoreTransform


class TestXmipptomoProtCC(BaseTest):
    """ This class check if the protocol to compute connected components works
    properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')

    def _runPreviousProtocols(self):
        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)
        self.launchProtocol(protImportTomogram)
        protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   boxSize=32,
                                                   samplingRate=5)
        self.launchProtocol(protImportCoordinates3d)
        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")
        self.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                             "There was a problem with coordinates 3d output")
        return protImportCoordinates3d

    def test_cc(self):
        protImport = self._runPreviousProtocols()
        protConnectedComponents = self.newProtocol(XmippProtConnectedComponents,
                                                   inputCoordinates=protImport.outputCoordinates,
                                                   distance=120)
        self.launchProtocol(protConnectedComponents)
        self.assertTrue(protConnectedComponents.outputSetOfCoordinates3D)
        self.assertEqual(protConnectedComponents.outputSetOfCoordinates3D.getSize(), 5)
        return protConnectedComponents


class TestXmipptomoProtProjectZ(BaseTest):
    """This class check if the protocol project top works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    def _runPreviousProtocols(self):
        protImport = self.newProtocol(ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5)
        self.launchProtocol(protImport)
        return protImport

    def _createProjZ(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms)
        self.launchProtocol(proj)
        self.assertIsNotNone(proj.outputParticles,
                             "There was a problem with particles output")
        return proj

    def _createProjY(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms,
                                dirParam=1)
        self.launchProtocol(proj)
        self.assertIsNotNone(proj.outputParticles,
                             "There was a problem with particles output")
        return proj

    def _createProjZ_range(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms,
                                rangeParam=1,
                                cropParam=20)
        self.launchProtocol(proj)
        self.assertIsNotNone(proj.outputParticles,
                             "There was a problem with particles output")
        return proj

    def _createRadAvgZ(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms,
                                radAvg=True)
        self.launchProtocol(proj)
        self.assertIsNotNone(proj.outputParticles,
                             "There was a problem with particles output")
        return proj

    def test_top(self):
        projection = self._createProjZ()
        outputParticles = getattr(projection, 'outputParticles')
        self.assertTrue(outputParticles)
        return projection

    def test_y(self):
        projection = self._createProjY()
        outputParticles = getattr(projection, 'outputParticles')
        self.assertTrue(outputParticles)
        return projection

    def test_top_range(self):
        projection = self._createProjZ_range()
        outputParticles = getattr(projection, 'outputParticles')
        self.assertTrue(outputParticles)
        return projection

    def test_top_radAvg(self):
        projection = self._createRadAvgZ()
        outputParticles = getattr(projection, 'outputParticles')
        self.assertTrue(outputParticles)
        return projection


class TestXmipptomoApplyTransf(BaseTest):
    """This class check if the protocol apply_alignment_subtomo works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')

    def _runPreviousProtocols(self):
        protPhantom = self.newProtocol(XmippProtPhantomSubtomo, option=1, nsubtomos=5)
        self.launchProtocol(protPhantom)
        self.assertIsNotNone(protPhantom.outputSubtomograms,
                             "There was a problem with subtomograms output")
        return protPhantom

    def _applyAlignment(self):
        protPhantom = self._runPreviousProtocols()
        apply = self.newProtocol(XmippProtApplyTransformSubtomo,
                                 inputSubtomograms=protPhantom.outputSubtomograms)
        self.launchProtocol(apply)
        self.assertIsNotNone(apply.outputSubtomograms,
                             "There was a problem with subtomograms output")
        self.assertIsNotNone(apply.outputAverage,
                             "There was a problem with average output")
        return apply

    def test_applyAlignment(self):
        align = self._applyAlignment()
        self.assertTrue(getattr(align, 'outputSubtomograms'))
        self.assertTrue(getattr(align, 'outputAverage'))
        return align


class TestXmipptomoMapback(BaseTest):
    """This class check if the protocol mapback works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')

    def _runPreviousProtocols(self):
        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)
        self.launchProtocol(protImportTomogram)
        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")

        protPhantom = self.newProtocol(XmippProtPhantomSubtomo,
                                       option=1,
                                       sampling=4,
                                       nsubtomos=5,
                                       coords=True,
                                       tomos=protImportTomogram.outputTomograms)
        self.launchProtocol(protPhantom)
        self.assertIsNotNone(protPhantom.outputSubtomograms,
                             "There was a problem with subtomograms output")

        return protImportTomogram, protPhantom

    def _mapback(self):
        protImportTomogram, protPhantom = self._runPreviousProtocols()
        mapback = self.newProtocol(XmippProtSubtomoMapBack,
                                   selection=1,
                                   inputSubtomos=protPhantom.outputSubtomograms,
                                   inputRef=protPhantom,
                                   inputTomograms=protImportTomogram.outputTomograms,
                                   removeBackground=True)
        mapback.inputRef.setExtended("outputSubtomograms.1")
        self.launchProtocol(mapback)
        self.assertIsNotNone(mapback.outputTomograms,
                             "There was a problem with tomograms output")
        return mapback

    def test_mapback(self):
        mapback = self._mapback()
        self.assertTrue(getattr(mapback, 'outputTomograms'))
        return mapback


class TestXmipptomoPhantom(BaseTest):
    """This class check if the protocol create phantom subtomo works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _phantom(self):
        phantom = self.newProtocol(XmippProtPhantomSubtomo, option=1)
        self.launchProtocol(phantom)
        self.assertIsNotNone(phantom.outputSubtomograms,
                             "There was a problem with subtomograms output")
        return phantom

    def test_phantom(self):
        phantom = self._phantom()
        self.assertTrue(getattr(phantom, 'outputSubtomograms'))
        self.assertEqual(phantom.outputSubtomograms.getFirstItem().getAcquisition().getAngleMax(), 90)
        return phantom

    def _phantomMW(self):
        phantomMW = self.newProtocol(XmippProtPhantomSubtomo,
                                     option=1,
                                     mwfilter=True)
        self.launchProtocol(phantomMW)
        self.assertIsNotNone(phantomMW.outputSubtomograms,
                             "There was a problem with subtomograms output")
        return phantomMW

    def test_phantomMW(self):
        phantomMW = self._phantomMW()
        self.assertTrue(getattr(phantomMW, 'outputSubtomograms'))
        self.assertEqual(phantomMW.outputSubtomograms.getFirstItem().getAcquisition().getAngleMax(), 60)
        self.assertEqual(phantomMW.outputSubtomograms.getFirstItem().getAcquisition().getAngleMin(), -60)
        return phantomMW


class XmippTomoScoreCoordinates(BaseTest):
    """This class check if the protocol project top works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('overview_wbp.em')

    def _runPreviousProtocols(self):
        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                      filesPath=self.tomogram,
                                      samplingRate=2,
                                      objLabel='Import Tomogram')
        self.launchProtocol(protImportTomogram)
        output = getattr(protImportTomogram, 'outputTomograms', None)
        self.assertIsNotNone(output,
                             "There was a problem with import tomograms output")

        protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                   objLabel='Import Coordinates - TXT',
                                                   auto=0,
                                                   filesPath=self.dataset.getPath(),
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   filesPattern='*.txt', boxSize=32,
                                                   samplingRate=5)
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
                                                outliers=False)
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
        self.assertTrue(filteredCoords.getSize() == 5)
        self.assertTrue(filteredCoords.getBoxSize() == 32)
        self.assertTrue(filteredCoords.getSamplingRate() == 5)

        # Test Outlier based filtering
        filteredCoords = self._createProtScoreOutliers(protCoords)
        # self.assertTrue(filteredCoords,
        #                 "There was a problem with score coordinates output")
        self.assertTrue(filteredCoords.getSize() == 0)
        self.assertTrue(filteredCoords.getBoxSize() == 32)
        self.assertTrue(filteredCoords.getSamplingRate() == 5)


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
