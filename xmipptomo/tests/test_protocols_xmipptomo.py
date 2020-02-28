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

from pyworkflow.utils import importFromPlugin
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSet
from xmipptomo.protocols import XmippProtSubtomoProject, XmippProtConnectedComponents
ProtImportCoordinates3D = importFromPlugin("tomo.protocols", "ProtImportCoordinates3D")
ProtImportTomograms = importFromPlugin("tomo.protocols", "ProtImportTomograms")
ProtImportSubTomograms = importFromPlugin("tomo.protocols", "ProtImportSubTomograms")

class TestXmippProtCC(BaseTest):
    """ This class check if the protocol to compute connected components works properly."""

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
    """This class check if the protocol project top in Xmipp works properly."""

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


