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

import pwem
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSet
from tomo.protocols import (ProtImportCoordinates3D,
                            ProtImportTomograms,
                            ProtImportSubTomograms)

from xmipptomo.protocols import XmippProtSubtomoProject, XmippProtConnectedComponents, XmippProtApplyTransformSubtomo


class TestXmippProtCC(BaseTest):
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


class TestXmippProtProjectZ(BaseTest):
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

class TestXmippApplyTransf(BaseTest):
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
        mltomo = pwem.Domain.importFromPlugin('xmipp2.protocols', 'Xmipp2ProtMLTomo')
        protMltomo = self.newProtocol(mltomo,
                                      inputVolumes=protImport.outputSubTomograms,
                                      randomInitialization=True,
                                      numberOfReferences=1,
                                      numberOfIters=3,
                                      angularSampling=30)
        self.launchProtocol(protMltomo)
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
