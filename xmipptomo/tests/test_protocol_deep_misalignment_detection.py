# *****************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
# *****************************************************************************

import os

from pyworkflow.tests import setupTestProject, DataSet, BaseTest
from tomo.protocols.protocol_import_tomograms import ProtImportTomograms
from xmipptomo.protocols.protocol_peak_high_contrast import XmippProtPeakHighContrast
from xmipptomo.protocols.protocol_deep_misalignment_detection import XmippProtDeepDetectMisalignment


class TestDeepMisaligmentDetectionBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def _runImportTomograms(cls, filesPath, pattern, samplingRate, objLabel):
        cls.protImportTomo = cls.newProtocol(ProtImportTomograms,
                                             filesPath=filesPath,
                                             filesPattern=pattern,
                                             samplingRate=samplingRate,
                                             objLabel=objLabel)

        cls.launchProtocol(cls.protImportTomo)

        return cls.protImportTomo

    @classmethod
    def _runPHC(cls, inputSoT, fiducialSize, boxSize, relaxedModeBool, relaxedModeThr, sampSlices, sdThr, numCoordThr,
                minCorrThr, mahalaThr):
        cls.protPHC = cls.newProtocol(XmippProtPeakHighContrast,
                                      inputSetOfTomograms=inputSoT,
                                      fiducialSize=fiducialSize,
                                      boxSize=boxSize,
                                      relaxedMode=relaxedModeBool,
                                      relaxedModeThr=relaxedModeThr,
                                      numberSampSlices=sampSlices,
                                      sdThr=sdThr,
                                      numberOfCoordinatesThr=numCoordThr,
                                      mirrorCorrelationThr=minCorrThr,
                                      mahalanobisDistanceThr=mahalaThr)

        cls.launchProtocol(cls.protPHC)

        return cls.protPHC

    @classmethod
    def _runDeepMisaliDetection(cls, inputSoC, tomoSource, inputSetOfT, misaliThrBool, misaliThr, modelPick,
                                misalignmentCriteria):
        cls.protDMD = cls.newProtocol(XmippProtDeepDetectMisalignment,
                                      inputSetOfCoordinates=inputSoC,
                                      tomoSource=tomoSource,
                                      inputSetOfTomograms=inputSetOfT,
                                      misaliThrBool=misaliThrBool,
                                      misaliThr=misaliThr,
                                      modelPick=modelPick,
                                      misalignmentCriteria=misalignmentCriteria)

        cls.launchProtocol(cls.protDMD)

        return cls.protDMD


class TestDeepMisaligmentDetection(TestDeepMisaligmentDetectionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.inputDataSet = DataSet.getDataSet('deepMisaliTomo')
        cls.inputTomo = cls.inputDataSet.getFile('tomo1')
        cls.inputTomoSR = 18.92

        cls.protImportTomo = cls._runImportTomograms(filesPath=os.path.split(cls.inputTomo)[0],
                                                     pattern=os.path.split(cls.inputTomo)[1],
                                                     samplingRate=cls.inputTomoSR,
                                                     objLabel="Import Tomograms")

        cls.fidSize = 8.0
        cls.boxSize = 32
        cls.relaxedModeBool = 1
        cls.relaxedModeThr = 3
        cls.sampSlices = 400
        cls.sdThr = 2.0
        cls.numCoordThr = 10
        cls.minCorrThr = 0.2
        cls.mahalaThr = 2.0

        cls.protPHC = cls._runPHC(inputSoT=cls.protImportTomo.Tomograms,
                                  fiducialSize=cls.fidSize,
                                  boxSize=cls.boxSize,
                                  relaxedModeBool=cls.relaxedModeBool,
                                  relaxedModeThr=cls.relaxedModeThr,
                                  sampSlices=cls.sampSlices,
                                  sdThr=cls.sdThr,
                                  numCoordThr=cls.numCoordThr,
                                  minCorrThr=cls.minCorrThr,
                                  mahalaThr=cls.mahalaThr)

        cls.tomoSource = 0
        cls.misaliThrBool = 1
        cls.misaliThr = 0.45
        cls.modelPick = 0
        cls.misalignmentCriteria = 0

        cls.protDMD = cls._runDeepMisaliDetection(inputSoC=cls.protPHC.outputSetOfCoordinates3D,
                                                  tomoSource=cls.tomoSource,
                                                  inputSetOfT=None,
                                                  misaliThrBool=cls.misaliThrBool,
                                                  misaliThr=cls.misaliThr,
                                                  modelPick=cls.modelPick,
                                                  misalignmentCriteria=cls.misalignmentCriteria)


    def test_importTomo(self):
        tomos = self.protImportTomo.Tomograms

        self.assertSetSize(tomos, size=1)

    def test_PHC(self):
        coords = self.protPHC.outputSetOfCoordinates3D

        self.assertSetSize(coords, size=55)
        self.assertEqual(coords.getSamplingRate(), self.inputTomoSR)

    def test_DMD(self):
        subtomos = self.protDMD.outputSubtomos
        tomosAli = self.protDMD.alignedTomograms

        self.assertSetSize(subtomos, size=55)
        self.assertSetSize(tomosAli, size=1)
        self.assertIsNone(self.protDMD.strongMisalignedTomograms)
        self.assertIsNone(self.protDMD.weakMisalignedTomograms)
