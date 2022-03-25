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
from os.path import exists, join, split
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from tomo.protocols import ProtImportTomograms
from xmipptomo.protocols import XmippProtCropTomograms
from xmipptomo.protocols.protocol_crop_tomograms import SUFIXCROPPED


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


class TestCropTomograms(TestReSizeBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestReSizeBase.setData()
        cls.protImportTomos = cls.runImportTomograms(cls.tomos, 16.14)

    def testCropTomogramsSamplingRate(self):
        Rrb = XmippProtCropTomograms()
        crop = self.newProtocol(XmippProtCropTomograms,
                                    objLabel='crop tomos',
                                    inputSet=self.protImportTomos.outputTomograms,
                                    xcrop0=100,
                                    xcropF=900,
                                    ycrop0=100,
                                    ycropF=900,
                                    zcrop0=10,
                                    zcropF=100)
        self.launchProtocol(crop)

        self.assertTrue(crop)
        self.assertSetSize(crop.outputSetOfTomograms, 2,
                           "There was a problem with the output (SetOfSubtomogram)")
