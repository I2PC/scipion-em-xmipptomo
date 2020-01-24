# **************************************************************************
# *
# * Authors:    Estrella Fernandez Gimenez [1]
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

from pyworkflow.utils import importFromPlugin
from pyworkflow.tests import BaseTest, setupTestProject
from tomo.tests import DataSet
from xmipptomo.protocols import XmippProtUnbinningCoord, XmippProtCCroi, XmippProtSubtomoProject
ProtImportCoordinates3D = importFromPlugin("tomo.protocols", "ProtImportCoordinates3D")
ProtImportTomograms = importFromPlugin("tomo.protocols", "ProtImportTomograms")
ProtImportSubTomograms = importFromPlugin("tomo.protocols", "ProtImportSubTomograms")


class TestXmippProtUnbinningCoord(BaseTest):
    """ This class check if the protocol to unbinning coordinates works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')

    def _runPreviousProtocols(self):
        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=self.tomogram,  # TODO: import binning 4 tomogram
                                              samplingRate=5)
        self.launchProtocol(protImportTomogram)
        protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                   auto=ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   filesPattern='', boxSize=32,
                                                   samplingRate=5)
        self.launchProtocol(protImportCoordinates3d)
        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")
        self.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                             "There was a problem with coordinates 3d output")
        return protImportCoordinates3d

    def _runUnbinningCoords(self):
        protImport = self._runPreviousProtocols()
        unbinning = self.newProtocol(XmippProtUnbinningCoord,
                                     inputCoordinates=protImport.outputCoordinates,
                                     factor=4)
        self.launchProtocol(unbinning)
        self.assertIsNotNone(unbinning.outputCoordinates,
                             "There was a problem with SetOfCoordinates output")
        return unbinning

    def test_basicUnbinning(self):
        xmipptomoUnbinning = self._runUnbinningCoords()
        outputCoordinates = getattr(xmipptomoUnbinning, 'outputCoordinates')
        self.assertTrue(outputCoordinates)
        self.assertTrue(outputCoordinates.getFirstItem().getX() == 1256)
        self.assertTrue(outputCoordinates.getFirstItem().getY() == 1400)
        self.assertTrue(outputCoordinates.getFirstItem().getZ() == 1024)
        return xmipptomoUnbinning


class TestXmippProtCCroi(BaseTest):
    """ This class check if the protocol to adjust coordinates to a roi works properly."""

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
        protImportCoordinates3d = self.newProtocol(ProtImportCoordinates3D,
                                                   auto=ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.outputTomograms,
                                                   filesPattern='', boxSize=32,
                                                   samplingRate=5)
        self.launchProtocol(protImportCoordinates3d)
        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")
        self.assertIsNotNone(protImportCoordinates3d.outputCoordinates,
                             "There was a problem with coordinates 3d output")
        return protImportCoordinates3d

    def _runCoordsRoi(self):
        protImport = self._runPreviousProtocols()
        coordsRoi = self.newProtocol(XmippProtCCroi,
                                     inputCoordinates=protImport.outputCoordinates,
                                     inputMesh=protImport.outputCoordinates,
                                     selection=0)
        self.launchProtocol(coordsRoi)
        self.assertIsNotNone(coordsRoi.outputCoordinates,
                             "There was a problem with SetOfCoordinates output")
        return coordsRoi

    def test_basicCoordsRoi(self):
        xmipptomoCoordsRoi = self._runCoordsRoi()
        outputCoordinates = getattr(xmipptomoCoordsRoi, 'outputCoordinates')
        self.assertTrue(outputCoordinates)
        return xmipptomoCoordsRoi


class TestXmippProtProjectZ(BaseTest):
    """This class check if the protocol project top in Xmipp works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    def _runPreviousProtocols(self):
        protImport = self.newProtocol(ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5)
        self.launchProtocol(protImport)
        return protImport

    def _createProjZ(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms)
        self.launchProtocol(proj)
        self.assertIsNotNone(proj.outputParticles,
                             "There was a problem with particles output")
        return proj

    def _createProjY(self):
        protImport = self._runPreviousProtocols()
        proj = self.newProtocol(XmippProtSubtomoProject,
                                input=protImport.outputSubTomograms,
                                dirParam=1)
        self.launchProtocol(proj)
        self.assertIsNotNone(proj.outputParticles,
                             "There was a problem with particles output")
        return proj

    def test_top(self):
        projection = self._createProjZ()
        outputParticles = getattr(projection, 'outputParticles')
        self.assertTrue(outputParticles)
        return projection

    def test_y(self):
        projection = self._createProjY()
        outputParticles = getattr(projection, 'outputParticles')
        self.assertTrue(outputParticles)
        return projection


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

import unittest, sys

from pyworkflow.em import exists
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.protocol import ProtImportVolumes

from xmipp3.protocols import XmippProtCropResizeVolumes
from xmipptomo.protocols import XmippProtMonoTomo

class TestMonoTomoBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='resmap'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.map3D = cls.dataset.getFile('betagal')
        cls.half1 = cls.dataset.getFile('betagal_half1')
        cls.half2 = cls.dataset.getFile('betagal_half2')

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                         )
        cls.launchProtocol(cls.protImport)
        cls.protResize = cls.newProtocol(XmippProtCropResizeVolumes,
                                         inputVolumes=cls.protImport.outputVolume,
                                         doResize=True,
                                         resizeOption=2,
                                         resizeFactor=5)
        cls.launchProtocol(cls.protResize)
        return cls.protResize


class TestMonoTomo(TestMonoTomoBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestMonoTomoBase.setData()
        cls.protImportHalf1 = cls.runImportVolumes(cls.half1, 3.54)
        cls.protImportHalf2 = cls.runImportVolumes(cls.half2, 3.54)

    def testMonoTomo(self):
        MonoTomo = self.newProtocol(XmippProtMonoTomo,
                                    objLabel='two halves monores',
                                    inputVolume=self.protImportHalf1.outputVol,
                                    inputVolume2=self.protImportHalf2.outputVol,
                                    provideMaskInHalves=True,
                                    useMask=False,
                                    minRes=1,
                                    maxRes=25,
                                    )
        self.launchProtocol(MonoTomo)
        self.assertTrue(exists(MonoTomo._getExtraPath('mgresolution.mrc')),
                        "MonoTomo has failed")

