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
import numpy as np
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSet

from xmipptomo.protocols import XmippProtScoreTransform
from tomo.protocols import ProtImportSubTomograms

class TestXmippTomoScoreTransform(BaseTest):
    """This class check if the protocol score_transform works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    def _runImportSubtomos(self, label):
        protImport = self.newProtocol(ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5,
                                      label=label)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputSubTomograms,
                             "There was a problem with import subtomograms output")
        return protImport.outputSubTomograms

    def _scoreTransformations(self, firstSubtomos, secondSubtomos):
        scoredTr = self.newProtocol(XmippProtScoreTransform,
                                   firstSubtomos=firstSubtomos,
                                   secondSubtomos=secondSubtomos,
                                   label='Score Subtomogram Transformations')
        self.launchProtocol(scoredTr)
        self.assertIsNotNone(scoredTr.outputSetOfSubtomogram,
                             "There was a problem with score subtomograms output")
        return scoredTr.outputSetOfSubtomogram

    def test_scoreTransformations(self):
        firstSubtomos = self._runImportSubtomos('First Subtomograms')
        secondSubtomos = self._runImportSubtomos('Second Subtomograms')
        scoredSubtomos = self._scoreTransformations(firstSubtomos, secondSubtomos)
        self.assertTrue(scoredSubtomos.getSize() == 4)
        dist = scoredSubtomos.aggregate(["MAX"], "distanceScore", ["distanceScore"])
        dist = np.asarray([d["distanceScore"] for d in dist])
        meanDist = np.mean(dist)
        self.assertTrue(meanDist == 0)
        return scoredSubtomos
