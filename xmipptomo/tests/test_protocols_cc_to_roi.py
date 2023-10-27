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


from pwem.protocols import ProtUnionSet
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSet
from tomo.protocols import (ProtImportCoordinates3D,
                            ProtImportTomograms,
                            TomoProtFitEllipsoid)

from xmipptomo.protocols import XmippProtCCroi


class TestXmipptomoProtCCtoROI(BaseTest):
    """ This class check if the protocol connected components to ROIs works properly."""

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
        self.assertIsNotNone(protImportTomogram.Tomograms,
                             "There was a problem with tomogram output")

        protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.Tomograms,
                                                   boxSize=32,
                                                   samplingRate=5)
        self.launchProtocol(protImportCoordinates3d)
        self.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                             "There was a problem with coordinates 3d output")

        protImportCoordinates3dBis = self.newProtocol(ProtImportCoordinates3D,
                                                      filesPath=self.coords3D,
                                                      importTomograms=protImportTomogram.Tomograms,
                                                      boxSize=32,
                                                      samplingRate=5)
        self.launchProtocol(protImportCoordinates3dBis)
        self.assertIsNotNone(protImportCoordinates3dBis.outputCoordinates,
                             "There was a problem with coordinates 3d output 2")

        protJoinCoordinates = self.newProtocol(ProtUnionSet,
                                               inputSets=[protImportCoordinates3d.outputCoordinates,
                                                          protImportCoordinates3dBis.outputCoordinates])
        self.launchProtocol(protJoinCoordinates)
        self.assertIsNotNone(protJoinCoordinates.outputSet,
                             "There was a problem with join coordinates output")

        protFitVesicles = self.newProtocol(TomoProtFitEllipsoid,
                                           input=protJoinCoordinates.outputSet,
                                           inputTomos=protImportTomogram.Tomograms)
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
