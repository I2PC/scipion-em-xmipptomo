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

from xmipptomo.protocols import XmippProtApplyTransformSubtomo, XmippProtPhantomSubtomo


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
