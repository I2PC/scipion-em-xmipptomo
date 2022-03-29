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
from xmipptomo.protocols import XmippProtResizeTomograms, XmippProtResizeTiltSeries
from xmipptomo.protocols.protocol_resize_base import XmippProtResizeBase
import os


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
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestReSizeBase.setData()
        cls.protImportTomos = cls.runImportTomograms(cls.tomos, 16.14)

    def testReSizeTomogramsSamplingRate(self):
        Rrb = XmippProtResizeTomograms()
        reSize = self.newProtocol(XmippProtResizeTomograms,
                                    objLabel='Resize tomos',
                                    inputSet=self.protImportTomos.outputTomograms,
                                    resizeOption = Rrb.RESIZE_SAMPLINGRATE,
                                    resizeSamplingRate = 32.28)
        self.launchProtocol(reSize)
        head1, tail1 = split(Rrb.outputTomoFileName(self.protImportTomos.outputTomograms[1]))
        head2, tail2 = split(Rrb.outputTomoFileName(self.protImportTomos.outputTomograms[2]))
        cond1 = exists(reSize._getExtraPath(join('tomo_1',tail1)))
        cond2 = exists(reSize._getExtraPath(join('tomo_2',tail2)))

        errstr = ''
        if not cond1:
            errstr = 'resize has failed in the samplingrate option'
        if not cond2:
            errstr = errstr + '\n resize has failed in the samplingrate option probably related with the use ' \
                     'of a SetOfTomograms (processing the second tomogram)'

        self.assertTrue(cond1 and cond2, errstr)

    def testReSizeTomogramsFactor(self):
        Rrb = XmippProtResizeTomograms()
        reSize = self.newProtocol(XmippProtResizeTomograms,
                                    objLabel='Resize tomos',
                                    inputSet=self.protImportTomos.outputTomograms,
                                    resizeOption = Rrb.RESIZE_FACTOR,
                                    resizeFactor = 0.5)
        self.launchProtocol(reSize)
        head1, tail1 = split(Rrb.outputTomoFileName(self.protImportTomos.outputTomograms[1]))
        head2, tail2 = split(Rrb.outputTomoFileName(self.protImportTomos.outputTomograms[2]))
        cond1 = exists(reSize._getExtraPath(join('tomo_1',tail1)))
        cond2 = exists(reSize._getExtraPath(join('tomo_2',tail2)))

        errstr = ''
        if not cond1:
            errstr = 'resize has failed in the Factor option'
        if not cond2:
            errstr = errstr + '\n resize has failed in the Factor option probably related with the use ' \
                     'of a SetOfTomograms (processing the second tomogram)'

        self.assertTrue(cond1 and cond2, errstr)

    def testReSizeTomogramsPiramid(self):
        Rrb = XmippProtResizeTomograms()
        reSize = self.newProtocol(XmippProtResizeTomograms,
                                    objLabel='Resize tomos',
                                    inputSet=self.protImportTomos.outputTomograms,
                                    resizeOption = Rrb.RESIZE_PYRAMID,
                                    resizeLevel = 0)
        self.launchProtocol(reSize)
        head1, tail1 = split(Rrb.outputTomoFileName(self.protImportTomos.outputTomograms[1]))
        head2, tail2 = split(Rrb.outputTomoFileName(self.protImportTomos.outputTomograms[2]))
        cond1 = exists(reSize._getExtraPath(join('tomo_1', tail1)))
        cond2 = exists(reSize._getExtraPath(join('tomo_2', tail2)))

        errstr = ''
        if not cond1:
            errstr = 'resize has failed in the pyramid option'
        if not cond2:
            errstr = errstr + '\n resize has failed in the pyramid option probably related with the use ' \
                              'of a SetOfTomograms (processing the second tomogram)'

        self.assertTrue(cond1 and cond2, errstr)


