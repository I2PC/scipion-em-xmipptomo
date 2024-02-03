# **************************************************************************
# *
# * Authors:    Jose Luis Vilas Prieto (jlvilas@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from tomo.protocols import ProtImportTomograms
from xmipptomo.protocols import XmippProtResizeTomograms


class TestReSizeBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='monotomo'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.tomos = cls.dataset.getFile('even_tomogram_rx*.mrc')

    @classmethod
    def runImportTomograms(cls, pattern, samplingRate):
        """ Run an Import tomograms protocol. """
        cls.protImport = cls.newProtocol(ProtImportTomograms,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                         )
        cls.launchProtocol(cls.protImport)
        return cls.protImport


class TestReSizeTomograms(TestReSizeBase):
    _objLabel = 'Resize tomos'

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestReSizeBase.setData()
        cls.protImportTomos = cls.runImportTomograms(cls.tomos, 16.14)

    def testReSizeTomogramsSamplingRate(self):
        protResizeTomos = XmippProtResizeTomograms()
        reSize = self.newProtocol(XmippProtResizeTomograms,
                                  objLabel=self._objLabel,
                                  inputSet=self.protImportTomos.Tomograms,
                                  resizeOption=protResizeTomos.RESIZE_SAMPLINGRATE,
                                  resizeSamplingRate=32.28)
        self.launchProtocol(reSize)
        self.assertTrue(reSize)
        self.assertSetSize(reSize.outputSetOfTomograms, 2,
                           "resize has failed in the samplingrate option probably related with the use "
                           "of a SetOfTomograms (processing the second tomogram)")

    def testReSizeTomogramsFactor(self):
        protResizeTomos = XmippProtResizeTomograms()
        reSize = self.newProtocol(XmippProtResizeTomograms,
                                  objLabel=self._objLabel,
                                  inputSet=self.protImportTomos.Tomograms,
                                  resizeOption=protResizeTomos.RESIZE_FACTOR,
                                  resizeFactor=0.5)
        self.launchProtocol(reSize)
        self.assertTrue(reSize)
        self.assertSetSize(reSize.outputSetOfTomograms, 2,
                           "resize has failed in the Factor option probably related with the use "
                           "of a SetOfTomograms (processing the second tomogram)")

    def testReSizeTomogramsPiramid(self):
        protResizeTomos = XmippProtResizeTomograms()
        reSize = self.newProtocol(XmippProtResizeTomograms,
                                  objLabel=self._objLabel,
                                  inputSet=self.protImportTomos.Tomograms,
                                  resizeOption=protResizeTomos.RESIZE_PYRAMID,
                                  resizeLevel=0)
        self.launchProtocol(reSize)
        self.assertTrue(reSize)
        self.assertSetSize(reSize.outputSetOfTomograms, 2,
                           "Resize has failed in the pyramid option probably related with the use "
                           "of a SetOfTomograms (processing the second tomogram)")
