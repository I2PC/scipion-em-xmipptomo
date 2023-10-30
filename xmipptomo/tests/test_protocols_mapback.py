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
from tomo.protocols import ProtImportTomograms

from xmipptomo.protocols import XmippProtSubtomoMapBack, XmippProtPhantomSubtomo


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
        self.assertIsNotNone(protImportTomogram.Tomograms,
                             "There was a problem with tomogram output")

        protPhantom = self.newProtocol(XmippProtPhantomSubtomo,
                                       option=1,
                                       sampling=4,
                                       nsubtomos=5,
                                       coords=True,
                                       tomos=protImportTomogram.Tomograms)
        self.launchProtocol(protPhantom)
        self.assertIsNotNone(protPhantom.outputSubtomograms,
                             "There was a problem with subtomograms output")

        return protImportTomogram, protPhantom

    def test_mapback(self):
        _, protPhantom = self._runPreviousProtocols()
        mapback = self.newProtocol(XmippProtSubtomoMapBack,
                                   selection=1,
                                   inputSubtomos=protPhantom.outputSubtomograms,
                                   inputRef=protPhantom,
                                   removeBackground=True)

        mapback.inputRef.setExtended("outputSubtomograms.1")
        self.launchProtocol(mapback)

        self.assertSetSize(getattr(mapback, XmippProtSubtomoMapBack._possibleOutputs.tomograms.name), 1,
                           "There was a problem with tomograms output")
