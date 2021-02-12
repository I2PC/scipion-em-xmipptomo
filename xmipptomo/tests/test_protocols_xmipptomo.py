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
    XmippProtSubtomoMapBack, XmippProtPhantomSubtomo, XmippProtScoreTransform


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
        self.assertTrue(protConnectedComponents.output3DCoordinates1)
        self.assertEqual(protConnectedComponents.output3DCoordinates1.getSize(), 2)
        self.assertTrue(protConnectedComponents.output3DCoordinates2)
        self.assertEqual(protConnectedComponents.output3DCoordinates2.getSize(), 2)
        self.assertTrue(protConnectedComponents.output3DCoordinates3)
        self.assertEqual(protConnectedComponents.output3DCoordinates3.getSize(), 1)
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
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    def _runPreviousProtocols(self):
        protImport = self.newProtocol(ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5)
        self.launchProtocol(protImport)
        from xmipp2.protocols import Xmipp2ProtMLTomo
        protMltomo = self.newProtocol(Xmipp2ProtMLTomo,
                                      inputVolumes=protImport.outputSubTomograms,
                                      randomInitialization=True,
                                      numberOfReferences=1,
                                      numberOfIters=3,
                                      angularSampling=30)
        self.launchProtocol(protMltomo)
        self.assertIsNotNone(protMltomo.outputSubtomograms,
                         "There was a problem with SetOfSubtomogram output")
        self.assertIsNotNone(protMltomo.outputClassesSubtomo,
                         "There was a problem with SetOfSubtomogram output")
        return protMltomo

    def _applyAlignment(self):
        protMltomo = self._runPreviousProtocols()
        apply = self.newProtocol(XmippProtApplyTransformSubtomo,
                                 inputSubtomograms=protMltomo.outputSubtomograms)
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
        protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   boxSize=32,
                                                   samplingRate=5)
        self.launchProtocol(protImportCoordinates3d)
        self.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                             "There was a problem with coordinates 3d output")
        from emantomo.protocols import EmanProtTomoExtraction
        protTomoExtraction = self.newProtocol(EmanProtTomoExtraction,
                                              inputTomograms=protImportTomogram.outputTomograms,
                                              inputCoordinates=protImportCoordinates3d.outputCoordinates,
                                              boxSize=32)
        self.launchProtocol(protTomoExtraction)
        self.assertIsNotNone(protTomoExtraction.outputSetOfSubtomogram,
                         "There was a problem with SetOfSubtomogram output")
        from xmipp2.protocols import Xmipp2ProtMLTomo
        protMltomo = self.newProtocol(Xmipp2ProtMLTomo,
                                      inputVolumes=protTomoExtraction.outputSetOfSubtomogram,
                                      randomInitialization=True,
                                      numberOfReferences=1,
                                      numberOfIters=3,
                                      angularSampling=30)
        self.launchProtocol(protMltomo)
        self.assertIsNotNone(protMltomo.outputSubtomograms,
                         "There was a problem with SetOfSubtomogram output")
        self.assertIsNotNone(protMltomo.outputClassesSubtomo,
                         "There was a problem with SetOfSubtomogram output")
        return protImportTomogram, protMltomo

    def _mapback(self):
        protImportTomogram, protMltomo = self._runPreviousProtocols()
        mapback = self.newProtocol(XmippProtSubtomoMapBack,
                                   inputClasses=protMltomo.outputClassesSubtomo,
                                   inputTomograms=protImportTomogram.outputTomograms)
        self.launchProtocol(mapback)
        self.assertIsNotNone(mapback.outputTomograms,
                             "There was a problem with tomograms output")

        mapback2 = self.newProtocol(XmippProtSubtomoMapBack,
                                   selection=1,
                                   inputSubtomos=protMltomo.outputSubtomograms,
                                   inputRef=protMltomo,
                                   inputTomograms=protImportTomogram.outputTomograms)
        mapback2.inputRef.setExtended("outputSubtomograms.1")
        self.launchProtocol(mapback2)
        self.assertIsNotNone(mapback2.outputTomograms,
                             "There was a problem with tomograms output")
        return mapback, mapback2

    def test_mapback(self):
        mapback, mapback2 = self._mapback()
        self.assertTrue(getattr(mapback, 'outputTomograms'))
        self.assertTrue(getattr(mapback2, 'outputTomograms'))
        return mapback, mapback2

    # def _mapback_subtomos(self):
        # protImportTomogram, protMltomo = self._runPreviousProtocols()
        # mapback = self.newProtocol(XmippProtSubtomoMapBack,
        #                            selection=1,
        #                            inputSubtomos=protMltomo.outputSubtomograms,
        #                            inputTomograms=protImportTomogram.outputTomograms)
        # mapback.inputRef.setExtended("outputSubtomograms.1")
        # self.launchProtocol(mapback)
        # self.assertIsNotNone(mapback.outputTomograms,
        #                      "There was a problem with tomograms output")
        # return mapback

    # def test_mapback_subtomos(self):
    #     mapback = self._mapback_subtomos()
    #     self.assertTrue(getattr(mapback, 'outputTomograms'))
    #     return mapback


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
        return phantom

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
