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

from xmipptomo.protocols import XmippProtPhantomSubtomo, XmippProtSubtractionSubtomo


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
        self.assertIsNotNone(createAverage.outputVolume,
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
                                       average=createAverage.outputVolume,
                                       maskBool=False)
        self.launchProtocol(subtraction)
        self.assertIsNotNone(subtraction.outputSubtomograms,
                             "There was a problem with subtracted subtomograms output")
        self.assertSetSize(subtraction.outputSubtomograms, 10, "There was a problem with particles output")

        srTol = 0.01
        srDiff = abs(subtraction.outputSubtomograms.getSamplingRate() - 1.0)
        self.assertTrue(srDiff < srTol)
        self.assertTrue(subtraction.outputSubtomograms.getFirstItem().getDim() == (40, 40, 40))
