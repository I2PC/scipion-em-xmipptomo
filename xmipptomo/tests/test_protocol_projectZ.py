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


from pyworkflow.tests import BaseTest, setupTestProject

from tomo.tests import DataSet
from tomo.protocols import ProtImportSubTomograms

from xmipptomo.protocols import XmippProtSubtomoProject

from xmipptomo.protocols.protocol_project_top import SubtomoProjectOutput


class TestXmipptomoProtProjectZ(BaseTest):
    """This class check if the protocol project top works properly."""

    importProt = None
    errorMsg = "There was a problem with particles output"



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
        self.assertSetSize(output, 4, self.errorMsg)
        return proj

    def _createProjY(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms,
                                dirParam=1)
        self.launchProtocol(proj)

        output = getattr(proj, SubtomoProjectOutput.particles.name)
        self.assertSetSize(output, 4, self.errorMsg)
        return proj

    def _createProjZ_range(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms,
                                rangeParam=1,
                                cropParam=20)
        self.launchProtocol(proj)
        output = getattr(proj, SubtomoProjectOutput.particles.name)
        self.assertSetSize(output, 4, self.errorMsg)
        return proj

    def _createRadAvgZ(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms,
                                radAvg=True)
        self.launchProtocol(proj)
        output = getattr(proj, SubtomoProjectOutput.particles.name)
        self.assertSetSize(output, 4, self.errorMsg)
        return proj

    def test_top(self):
        self._createProjZ()

    def test_y(self):
        self._createProjY()

    def test_top_range(self):
        self._createProjZ_range()

    def test_top_radAvg(self):
        self._createRadAvgZ()
