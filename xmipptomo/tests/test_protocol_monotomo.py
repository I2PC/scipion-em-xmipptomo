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

from os.path import exists

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwem.protocols import ProtImportVolumes

from tomo.protocols import ProtImportTomograms

from xmipptomo.protocols import XmippProtMonoTomo


class TestMonoTomoBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.map3D = cls.dataset.getFile('betagal')
        cls.half1 = cls.dataset.getFile('even_tomogram_rx*.mrc')
        cls.half2 = cls.dataset.getFile('odd_tomogram_rx*.mrc')

    @classmethod
    def runImportTomograms(cls, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportTomograms,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                         )
        cls.launchProtocol(cls.protImport)
        return cls.protImport


class TestMonoTomo(TestMonoTomoBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestMonoTomoBase.setData()
        cls.protImportHalf1 = cls.runImportTomograms(cls.half1, 16.14)
        cls.protImportHalf2 = cls.runImportTomograms(cls.half2, 16.14)

    def testMonoTomo(self):
        MonoTomo = self.newProtocol(XmippProtMonoTomo,
                                    objLabel='two halves monotomo',
                                    oddTomograms=self.protImportHalf1.outputTomograms,
                                    evenTomograms=self.protImportHalf2.outputTomograms,
                                    useMask=False,
                                    minRes=1,
                                    maxRes=150,
                                    )
        self.launchProtocol(MonoTomo)
        self.assertTrue(exists(MonoTomo._getExtraPath('tomo_1/fullTomogram_1.mrc')), "MonoTomo has failed")
