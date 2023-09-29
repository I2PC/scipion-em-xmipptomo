# **************************************************************************
# *
# * Authors:    Martín Salinas (ssalinasmartin@gmail.com)
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

# General imports
import os

# Scipion em imports
from pyworkflow.tests import setupTestProject

# External plugin imports
from tomo.protocols import ProtImportSubTomograms

# Plugin imports
from xmipptomo.utils import removeTmpElements
from xmipptomo.protocols.protocol_extract_subtomos import OUTPUTATTRIBUTE as EXTRACT_SUBTOMOS_OUTPUTATTRIBUTE
from xmipptomo.protocols.protocol_project_subtomograms import OUTPUTATTRIBUTE as PROJECT_SUBTOMOGRAMS_OUTPUTATTRIBUTE
from xmipptomo.protocols import XmippProtProjectSubtomograms
from xmipptomo.tests.test_protocol_extract_subtomos import TestXmippProtExtractSubtomosProts, TestXmippProtExtractSubtomos

class TestXmippProtProjectSubtomograms(TestXmippProtExtractSubtomosProts):
    """This class check if the protocol to generate projections from subtomograms in Xmipptomo works properly."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.tmpElements = []

        # Getting XmippExtractSubtomos test so it can produce necessary inputs for this test
        extractSubtomos = TestXmippProtExtractSubtomos()
        extractSubtomosPath = cls.getOutputPath().replace(cls.__name__, type(extractSubtomos).__name__)

        # Check if temp project already exists
        if not os.path.isdir(extractSubtomosPath):
            # If not, add to delete after end of this test so no extra files and directories are left behind
            cls.tmpElements.append(extractSubtomosPath)

        # Launching XmippProtExtractSubtomos and obtaining output
        extractSubtomos.setUpClass()
        protocolOutput = extractSubtomos._runXmippTomoExtraction(doInvert=True)
        subtomograms = getattr(protocolOutput, EXTRACT_SUBTOMOS_OUTPUTATTRIBUTE, None)
        cls.assertIsNotNone(subtomograms, "There was an error obtaining the input subtomograms.")

        # Getting output details to input for import subtomograms protocol
        samplingRate = subtomograms.getSamplingRate()
        subtomogramsPath = os.path.abspath(protocolOutput._getExtraPath('overview_wbp'))

        # Setting up project again to overwrite temp project variables
        setupTestProject(cls)

        # Importing subtomograms
        cls.subtomograms = cls._runImportSetOfSubTomograms(cls, subtomogramsPath, samplingRate)

    def _runImportSetOfSubTomograms(self, subtomogramsPath, samplingRate):
        """This function imports the files produced by XmippExtractSubtomo."""
        protImportSetOfSubTomograms = self.newProtocol(
            ProtImportSubTomograms,
            filesPath=subtomogramsPath,
            filesPattern='*.mrc',
            copyFiles=True,
            samplingRate=samplingRate
        )
        self.launchProtocol(protImportSetOfSubTomograms)
        outputSubtomograms = getattr(protImportSetOfSubTomograms, "outputSubTomograms", None)
        self.assertIsNotNone(outputSubtomograms, "There was a problem importing subtomograms.")
        return outputSubtomograms

    def _runXmippProjectSubtomograms(self):
        """This function creates and runs a XmippProtProjectSubtomograms protocol with controlled params."""
        protXmippProjectSubtomograms = self.newProtocol(
            XmippProtProjectSubtomograms,
            inputSubtomograms=self.subtomograms,
            tiltRangeNSamples=40,
            numberOfMpi=1,
            numberOfThreads=5
        )
        self.launchProtocol(protXmippProjectSubtomograms)
        self.assertIsNotNone(getattr(protXmippProjectSubtomograms, PROJECT_SUBTOMOGRAMS_OUTPUTATTRIBUTE, None), "Projections were not properly generated.")
    
    def test_projectSubtomos(self):
        """This function runs XmippProtProjectSubtomograms using the output of XmippExtractSubtomos as input."""
        self._runXmippProjectSubtomograms()
        # Last test calls cleaning function so it does not count as a separate test
        removeTmpElements(self.tmpElements)