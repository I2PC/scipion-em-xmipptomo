# **************************************************************************
# *
# * Authors:    J.L. Vilas (jlvilas@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.protocols.protocol_import_coordinates import IMPORT_FROM_EMAN
from ..protocols import *
import tomo.protocols

from xmipptomo.protocols.protocol_extract_subtomos import OUTPUTATTRIBUTE
from xmipptomo.protocols import XmippProtExtractSubtomos


class TestXmippProtExtractSubtomosBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def setData(cls, projectData='tomo-em'):
        from tomo.tests import DataSet
        cls.dataset = DataSet.getDataSet(projectData)
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')
        cls.coords3D_Large = cls.dataset.getFile('overview_wbp_large.txt')
        cls.inputSetOfSubTomogram = cls.dataset.getFile('subtomo')
        cls.smallTomogram = cls.dataset.getFile('coremask_normcorona.mrc')

class TestXmippProtExtractSubtomosProts(TestXmippProtExtractSubtomosBase):
    """
    This prepares the protocols to perform the necessary tests.
    """
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippProtExtractSubtomosProts.setData()
        
    def _runImportCoordinatesAndTomograms(self):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                   auto=IMPORT_FROM_EMAN,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   filesPattern='', boxSize=32,
                                                   samplingRate=5)

        self.launchProtocol(protImportCoordinates3d)
        self.assertSetSize(protImportTomogram.outputTomograms, 1,
                           "There was a problem with tomogram output")
        self.assertSetSize(protImportCoordinates3d.outputCoordinates, 5,
                           "There was a problem with coordinates 3d output")
        return protImportCoordinates3d, protImportTomogram

    def _runXmippTomoExtraction(self, doInvert=False, boxSize=32, differenttomogram = False):
        protImportCoordinates3d, protImportTomogram = self._runImportCoordinatesAndTomograms()
        if differenttomogram:
            protTomoExtraction = self.newProtocol(XmippProtExtractSubtomos,
                                                  tomograms=protImportTomogram.outputTomograms,
                                                  coords=protImportCoordinates3d.outputCoordinates,
                                                  invertContrast=doInvert,
                                                  boxSize=boxSize)
        else:
            protTomoExtraction = self.newProtocol(XmippProtExtractSubtomos,
                                                  coords=protImportCoordinates3d.outputCoordinates,
                                                  invertContrast=doInvert,
                                                  boxSize=boxSize)
        self.launchProtocol(protTomoExtraction)
        self.assertSetSize(getattr(protTomoExtraction, OUTPUTATTRIBUTE), 5,
                           "There was a problem with SetOfSubtomogram output")
        return protTomoExtraction

class TestXmippProtExtractSubtomos(TestXmippProtExtractSubtomosProts):
    """
    This class check if the protocol to extract subtomograms in Xmipptomo works properly.
    """
    def test_extractParticlesWithDoInvert(self):
        protTomoExtraction = self._runXmippTomoExtraction(doInvert=True)
        output = getattr(protTomoExtraction, OUTPUTATTRIBUTE)
        self.assessOutput(output)

    def test_extractParticlesModifiedBoxSize(self):
        protTomoExtraction = self._runXmippTomoExtraction(boxSize=64)
        output =getattr(protTomoExtraction, OUTPUTATTRIBUTE)
        self.assessOutput(output)

    def test_extractParticlesFromDifferentTomogram(self):
        protTomoExtraction = self._runXmippTomoExtraction(differenttomogram=True)
        output =getattr(protTomoExtraction, OUTPUTATTRIBUTE)
        self.assessOutput(output)

    def assessOutput(self, outputSet, size=5):
        self.assertSetSize(outputSet, size)
        self.assertTrue(outputSet.hasCoordinates3D())