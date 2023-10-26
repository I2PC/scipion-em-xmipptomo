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


import os

from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSet
from tomo.protocols import ProtImportTs

from xmipptomo.protocols import XmippProtSplitTiltSeries


class TestXmipptomoSplitTS(BaseTest):
    """This class check if the protocol split tilt-series works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('tomo-em')
        cls.inputSoTS = cls.inputDataSet.getFile('tutorialData/BB*.st')


    def runWorkflow(self):
        protImportTS = self.newProtocol(ProtImportTs,
                                        filesPath=os.path.split(self.inputSoTS)[0],
                                        filesPattern="BB{TS}.st",
                                        voltage=300,
                                        anglesFrom=0,
                                        magnification=105000,
                                        sphericalAberration=2.7,
                                        amplitudeContrast=0.1,
                                        samplingRate=20.2,
                                        doseInitial=0,
                                        dosePerFrame=3.0,
                                        minAngle=-55,
                                        maxAngle=65.0,
                                        stepAngle=2.0,
                                        tiltAxisAngle=-12.5)

        self.launchProtocol(protImportTS)

        protSplitTS = self.newProtocol(XmippProtSplitTiltSeries,
                                       inputSetOfTiltSeries=protImportTS.outputTiltSeries)

        self.launchProtocol(protSplitTS)

        return protSplitTS

    def testSplitTS(self):
        protSplitTS = self.runWorkflow()

        self.assertIsNotNone(protSplitTS.outputEvenSetOfTiltSeries,
                             "Even set of tilt series has not been generated")
        self.assertIsNotNone(protSplitTS.outputOddSetOfTiltSeries,
                             "Odd set of tilt series has not been generated")

        self.assertTrue(protSplitTS.outputEvenSetOfTiltSeries.getSamplingRate() == 20.20)
        self.assertTrue(protSplitTS.outputOddSetOfTiltSeries.getSamplingRate() == 20.20)

        self.assertSetSize(protSplitTS.outputEvenSetOfTiltSeries, 2, "Missing TS in even set")
        self.assertSetSize(protSplitTS.outputOddSetOfTiltSeries, 2, "Missing TS in odd set")

        self.assertTrue(protSplitTS.outputEvenSetOfTiltSeries.getFirstItem().getDim() == (512, 512, 30))
        self.assertTrue(protSplitTS.outputOddSetOfTiltSeries.getFirstItem().getDim() == (512, 512, 31))