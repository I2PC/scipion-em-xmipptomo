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
import os

import numpy as np


from pwem.protocols import ProtUnionSet
from pwem.emlib.image import ImageHandler
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSet
from tomo.protocols import (ProtImportCoordinates3D,
                            ProtImportTomograms,
                            ProtImportSubTomograms,
                            ProtImportTs)


from xmipptomo.protocols import XmippProtSubtomoProject, XmippProtConnectedComponents, XmippProtApplyTransformSubtomo, \
    XmippProtSubtomoMapBack, XmippProtPhantomSubtomo, XmippProtScoreCoordinates, XmippProtScoreTransform, \
    XmippProtHalfMapsSubtomo, XmippProtPhantomTomo, XmippProtSubtractionSubtomo, XmippProtCCroi, XmippProtFitEllipsoid, XmippProtSplitTiltSeries

from xmipptomo.protocols.protocol_phantom_tomo import OutputPhantomTomos
from xmipptomo.protocols.protocol_project_top import SubtomoProjectOutput


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

    importProt = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    @classmethod
    def _runPreviousProtocols(cls):
        if cls.importProt is None:
            protImport = cls.newProtocol(ProtImportSubTomograms,
                                         filesPath=cls.setOfSubtomograms,
                                         samplingRate=5)
            cls.launchProtocol(protImport)
            cls.importProt = protImport
        return cls.importProt

    def _createProjZ(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms)
        self.launchProtocol(proj)
        output = getattr(proj, SubtomoProjectOutput.particles.name)
        self.assertSetSize(output, 4, "There was a problem with particles output")
        return proj

    def _createProjY(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms,
                                dirParam=1)
        self.launchProtocol(proj)

        output = getattr(proj, SubtomoProjectOutput.particles.name)
        self.assertSetSize(output, 4, "There was a problem with particles output")
        return proj

    def _createProjZ_range(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms,
                                rangeParam=1,
                                cropParam=20)
        self.launchProtocol(proj)
        output = getattr(proj, SubtomoProjectOutput.particles.name)
        self.assertSetSize(output, 4, "There was a problem with particles output")
        return proj

    def _createRadAvgZ(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms,
                                radAvg=True)
        self.launchProtocol(proj)
        output = getattr(proj, SubtomoProjectOutput.particles.name)
        self.assertSetSize(output, 4, "There was a problem with particles output")
        return proj

    def test_top(self):
        self._createProjZ()

    def test_y(self):
        self._createProjY()

    def test_top_range(self):
        self._createProjZ_range()

    def test_top_radAvg(self):
        self._createRadAvgZ()


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

    def test_mapback(self):
        protImportTomogram, protPhantom = self._runPreviousProtocols()
        mapback = self.newProtocol(XmippProtSubtomoMapBack,
                                   selection=1,
                                   inputSubtomos=protPhantom.outputSubtomograms,
                                   inputRef=protPhantom,
                                   removeBackground=True)

        mapback.inputRef.setExtended("outputSubtomograms.1")
        self.launchProtocol(mapback)

        self.assertSetSize(getattr(mapback, XmippProtSubtomoMapBack._possibleOutputs.tomograms.name), 1,
                           "There was a problem with tomograms output")


class TestXmippTomoPhantom(BaseTest):
    """This class check if the protocol create phantom tomograms works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def test_phantom(self):
        phantom = self.newProtocol(XmippProtPhantomTomo,
                                   ntomos=2,
                                   nparticles=10)
        self.launchProtocol(phantom)
        tomograms = getattr(phantom, OutputPhantomTomos.tomograms.name)
        self.assertSetSize(tomograms, 2, "There was a problem with subtomograms output")

        self.assertIsNotNone(tomograms.getAcquisition().getMagnification(),
                             "Acquisition magnification not populated properly")

        coordinates = getattr(phantom, OutputPhantomTomos.coordinates3D.name)
        self.assertSetSize(coordinates, 20, "There was a problem with subtomograms output")

        return phantom


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


class TestXmipptomoHalfMaps(BaseTest):
    """This class check if the protocol to create half maps from a SetOfSubtomograms
    works properly"""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _phantom(self):
        phantom = self.newProtocol(XmippProtPhantomSubtomo, option=1)
        self.launchProtocol(phantom)
        self.assertIsNotNone(phantom.outputSubtomograms,
                             "There was a problem with subtomograms output")
        return phantom

    def _half_maps(self, subtomos):
        half_maps = self.newProtocol(XmippProtHalfMapsSubtomo, inputSubtomograms=subtomos)
        self.launchProtocol(half_maps)
        self.assertIsNotNone(half_maps.halfMaps,
                             "There was a problem with half maps output")
        return half_maps

    def test_half_maps(self):
        phantom = self._phantom()
        half_maps = self._half_maps(phantom.outputSubtomograms)
        map_even = np.squeeze(ImageHandler().read(half_maps.halfMaps[1].getFileName()).getData())
        map_odd = np.squeeze(ImageHandler().read(half_maps.halfMaps[2].getFileName()).getData())
        error = np.sqrt(np.sum((map_even - map_odd) ** 2) / map_even.size)
        self.assertAlmostEqual(error, 0.0, delta=0.1, msg="Unexpected half maps")
        return half_maps


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

class TestXmipptomoSubtractionSubtomo(BaseTest):
    """This class check if the protocol subtraction_subtomo works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _runPreviousProtocols(self):
        createAverage = self.newProtocol(XmippProtPhantomSubtomo,
                                         option=1,
                                         nsubtomos=1)
        self.launchProtocol(createAverage)
        self.assertIsNotNone(createAverage.outputSubtomograms,
                             "There was a problem with subtomogram average phantom output")
        createSubtomos = self.newProtocol(XmippProtPhantomSubtomo,
                                          option=1,
                                          nsubtomos=10,
                                          rotate=True,
                                          applyShift=True)
        self.launchProtocol(createSubtomos)
        self.assertIsNotNone(createSubtomos.outputSubtomograms,
                             "There was a problem with subtomogram phantoms output")
        return createAverage, createSubtomos

    def test_subtraction(self):
        createAverage, createSubtomos = self._runPreviousProtocols()
        subtraction = self.newProtocol(XmippProtSubtractionSubtomo,
                                       inputSubtomos=createSubtomos.outputSubtomograms,
                                       average=createAverage.outputSubtomograms,
                                       maskBool=False)
        self.launchProtocol(subtraction)
        self.assertIsNotNone(subtraction.outputSubtomograms,
                             "There was a problem with subtracted subtomograms output")
        self.assertSetSize(subtraction.outputSubtomograms, 10, "There was a problem with particles output")
        self.assertTrue(subtraction.outputSubtomograms.getSamplingRate() == 1.0)
        self.assertTrue(subtraction.outputSubtomograms.getFirstItem().getDim() == (40, 40, 40))


class TestXmipptomoProtFitVesicles(BaseTest):
    """ This class check if the protocol fit vesicles works properly."""
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

        protImportCoordinates3d_2 = self.newProtocol(ProtImportCoordinates3D,
                                                     filesPath=self.coords3D,
                                                     importTomograms=protImportTomogram.outputTomograms,
                                                     boxSize=32,
                                                     samplingRate=5)
        self.launchProtocol(protImportCoordinates3d_2)
        self.assertIsNotNone(protImportCoordinates3d_2.outputCoordinates,
                             "There was a problem with coordinates 3d output 2")

        protJoinCoordinates = self.newProtocol(ProtUnionSet,
                                               inputSets=[protImportCoordinates3d.outputCoordinates,
                                                          protImportCoordinates3d_2.outputCoordinates])
        self.launchProtocol(protJoinCoordinates)
        self.assertIsNotNone(protJoinCoordinates.outputSet,
                             "There was a problem with join coordinates output")

        return protJoinCoordinates, protImportTomogram

    def test_fitEllipsoid(self):
        protJoinCoordinates, protImportTomogram = self._runPreviousProtocols()
        protFitVesicles = self.newProtocol(XmippProtFitEllipsoid,
                                           input=protJoinCoordinates.outputSet,
                                           inputTomos=protImportTomogram.outputTomograms)
        self.launchProtocol(protFitVesicles)
        self.assertIsNotNone(protFitVesicles.outputMeshes, "There was a problem with output vesicles (SetOfMeshes)")
        self.assertEqual(protFitVesicles.outputMeshes.getSize(), 11184)
        return protFitVesicles


class TestXmipptomoProtCCtoROI(BaseTest):
    """ This class check if the protocol connected componnents to ROIs works properly."""

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

        protImportCoordinates3d_2 = self.newProtocol(ProtImportCoordinates3D,
                                                     filesPath=self.coords3D,
                                                     importTomograms=protImportTomogram.outputTomograms,
                                                     boxSize=32,
                                                     samplingRate=5)
        self.launchProtocol(protImportCoordinates3d_2)
        self.assertIsNotNone(protImportCoordinates3d_2.outputCoordinates,
                             "There was a problem with coordinates 3d output 2")

        protJoinCoordinates = self.newProtocol(ProtUnionSet,
                                               inputSets=[protImportCoordinates3d.outputCoordinates,
                                                          protImportCoordinates3d_2.outputCoordinates])
        self.launchProtocol(protJoinCoordinates)
        self.assertIsNotNone(protJoinCoordinates.outputSet,
                             "There was a problem with join coordinates output")

        protFitVesicles = self.newProtocol(XmippProtFitEllipsoid,
                                           input=protJoinCoordinates.outputSet,
                                           inputTomos=protImportTomogram.outputTomograms)
        self.launchProtocol(protFitVesicles)
        self.assertIsNotNone(protFitVesicles.outputMeshes, "There was a problem with output vesicles (SetOfMeshes)")
        return protFitVesicles, protJoinCoordinates

    def test_CCtoROIs(self):
        protFitVesicles, protJoinCoordinates = self._runPreviousProtocols()
        protCCtoROI = self.newProtocol(XmippProtCCroi,
                                       inputCoordinates=protJoinCoordinates.outputSet,
                                       inputMeshes=protFitVesicles.outputMeshes)
        self.launchProtocol(protCCtoROI)
        self.assertIsNotNone(protCCtoROI.outputSet, "There was a problem with CCtoROIs output")
        self.assertEqual(protCCtoROI.outputSet.getSize(), 10)
        self.assertEqual(protCCtoROI.outputSet.getBoxSize(), 32)
        return protCCtoROI

    def test_CCtoROIsPoints(self):
        protFitVesicles, protJoinCoordinates = self._runPreviousProtocols()
        protCCtoROI = self.newProtocol(XmippProtCCroi,
                                       inputCoordinates=protJoinCoordinates.outputSet,
                                       inputMeshes=protFitVesicles.outputMeshes,
                                       selection=1)
        self.launchProtocol(protCCtoROI)
        self.assertIsNotNone(protCCtoROI.outputSet, "There was a problem with CCtoROIs points output")
        self.assertEqual(protCCtoROI.outputSet.getSize(), 10)
        self.assertEqual(protCCtoROI.outputSet.getBoxSize(), 32)


class TestXmipptomoSplitTS(BaseTest):
    """This class check if the protocol split tilt-series works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('tomo-em')
        cls.inputSoTS = cls.inputDataSet.getFile('ts1_imod')

    def runWorkflow(self):
        protImportTS = self.newProtocol(ProtImportTs,
                                        filesPath=os.path.split(self.inputSoTS)[0],
                                        filesPattern="BB{TS}.st",
                                        voltage=300,
                                        anglesFrom=0,
                                        magnification=105000,
                                        sphericalAberration=2.7,
                                        amplitudeContrast=0.1,
                                        samplingRate=20.2,
                                        doseInitial=0,
                                        dosePerFrame=3.0,
                                        minAngle=-55,
                                        maxAngle=65.0,
                                        stepAngle=2.0,
                                        tiltAxisAngle=-12.5)

        self.launchProtocol(protImportTS)

        protSplitTS = self.newProtocol(XmippProtSplitTiltSeries,
                                       inputSetOfTiltSeries=protImportTS.outputTiltSeries)

        self.launchProtocol(protSplitTS)

        return protSplitTS

    def testSplitTS(self):
        protSplitTS = self.runWorkflow()

        self.assertIsNotNone(protSplitTS.outputEvenSetOfTiltSeries,
                             "Even set of tilt series has not been generated")
        self.assertIsNotNone(protSplitTS.outputOddSetOfTiltSeries,
                             "Odd set of tilt series has not been generated")

        self.assertTrue(protSplitTS.outputEvenSetOfTiltSeries.getSamplingRate() == 20.20)
        self.assertTrue(protSplitTS.outputOddSetOfTiltSeries.getSamplingRate() == 20.20)

        self.assertSetSize(protSplitTS.outputEvenSetOfTiltSeries, 2, "Missing TS in even set")
        self.assertSetSize(protSplitTS.outputOddSetOfTiltSeries, 2, "Missing TS in odd set")

        self.assertTrue(protSplitTS.outputEvenSetOfTiltSeries.getFirstItem().getDim() == (512, 512, 30))
        self.assertTrue(protSplitTS.outputOddSetOfTiltSeries.getFirstItem().getDim() == (512, 512, 31))
