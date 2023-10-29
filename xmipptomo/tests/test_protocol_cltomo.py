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
import tomo.protocols

from xmipptomo.protocols.protocol_cltomo import XmippProtCLTomo
from xmipptomo.protocols.protocol_extract_subtomos import XmippProtExtractSubtomos

## Tomogram type constants for particle extraction
#OUTPUTATTRIBUTE = 'Subtomograms'

class TestXmippProtCLTomoBase(BaseTest):
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


class TestXmippProtCLtomo(TestXmippProtCLTomoBase):
    """
    This prepares the protocols to perform the necessary tests.
    """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestXmippProtCLtomo.setData()

    def _runImportCoordinatesAndTomograms(self):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                   auto=IMPORT_FROM_EMAN,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.Tomograms,
                                                   filesPattern='', boxSize=32,
                                                   samplingRate=5)

        self.launchProtocol(protImportCoordinates3d)
        self.assertSetSize(protImportTomogram.Tomograms, 1,
                           "There was a problem with tomogram output")
        self.assertSetSize(protImportCoordinates3d.outputCoordinates, 5,
                           "There was a problem with coordinates 3d output")

        protTomoExtraction = self.newProtocol(XmippProtExtractSubtomos,
                                                  tomograms=protImportTomogram.Tomograms,
                                                  coords=protImportCoordinates3d.outputCoordinates,
                                                  invertContrast=True,
                                                  boxSize=32)

        self.launchProtocol(protTomoExtraction)
        self.assertSetSize(getattr(protTomoExtraction, 'Subtomograms'), 5,
                           "There was a problem with SetOfSubtomogram output")

        return protTomoExtraction

    def _runXmippCLtomo(self):
        protTomoExtraction = self._runImportCoordinatesAndTomograms()
        protCLTomo = self.newProtocol(XmippProtCLTomo, inputVolumes=protTomoExtraction.Subtomograms)

        self.launchProtocol(protCLTomo)


        return protCLTomo


class TestXmippProtCLtomo(TestXmippProtCLtomo):
    """
    This class check if CLtomos in Xmipptomo works properly.
    """

    def test_CLtomo(self):
        protCLTomo = self._runXmippCLtomo()
        msg = "There was a problem with CLtomo"
        self.assertIsNotNone(protCLTomo.outputClasses, msg)

