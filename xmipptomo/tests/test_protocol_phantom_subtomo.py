# **************************************************************************
# *
# * Authors:    Jose Luis Vilas [1]
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from xmipp3.protocols import XmippProtConvertPdb

from xmipptomo.protocols import XmippProtPhantomSubtomo


class TestXmippSubtomoPhantom(BaseTest):
    """This class check if the protocol create phantom subtomo works properly."""

    @classmethod
    def setUpClass(cls):
        cls.dataset = DataSet.getDataSet('nma')
        cls.pdb = cls.dataset.getFile('pdb')
        setupTestProject(cls)

        setSize = True
        sizeX = 128
        sizeY = 128
        sizeZ = 128
        pdbid = "3j3i"
        sampling = 1
        cls.protConvert = cls._runVolumeFromPDB(XmippProtConvertPdb, pdbid, sampling, setSize, sizeZ, sizeY, sizeX)


        inputVolume = cls.protConvert.outputVolume
        option = 0
        nsubtomos  = 20
        mwfilter = False
        rotate = False
        applyShift = False
        coords = False
        addNoise = False
        cls.phantom = cls._runPhantomSubtomo_NoMW(XmippProtPhantomSubtomo, option, inputVolume, sampling,
                                  nsubtomos, mwfilter, rotate, applyShift, coords, addNoise)

        addNoise = True
        differentStatistics = True
        minstd = 2
        maxStd = 10
        cls.phantom_NoMW_noisy = cls._runPhantomSubtomo_NoMW_noisy(XmippProtPhantomSubtomo, option, inputVolume, sampling,
                                nsubtomos, mwfilter, rotate, applyShift, coords, addNoise, differentStatistics,
                                minstd, maxStd)

        mwfilter = True
        mwangle = 60
        cls.phantom_MW_noisy = cls._runPhantomSubtomo_MW_noisy(XmippProtPhantomSubtomo, option, inputVolume, sampling,
                                nsubtomos, mwfilter, mwangle, rotate, applyShift, coords, addNoise, differentStatistics,
                                minstd, maxStd)

        rotate = True
        uniformAngularDistribution = True
        applyShift = True
        xmin = 0
        xmax = 5
        ymin = 0
        ymax = 5
        zmin = 0
        zmax = 5
        cls.phantom_MW_noisy_Randomrotation_shift = cls._runPhantomSubtomo_MW_noisy_Randomrotation_shift(XmippProtPhantomSubtomo,
                                option, inputVolume, sampling, nsubtomos, mwfilter, mwangle, rotate, uniformAngularDistribution,
                                applyShift, xmin, xmax, ymin, ymax, zmin, zmax, coords, addNoise, differentStatistics, minstd, maxStd)

    @classmethod
    def _runVolumeFromPDB(cls, protocolName, pdbid, sampling, setSize, size_z, size_y, size_x):

        print("Run convert a pdb from database")
        cls.protConvert = cls.newProtocol(protocolName, pdbId=pdbid, sampling=sampling, setSize=setSize,
                                       size_z=size_z, size_y=size_y, size_x=size_x)
        cls.launchProtocol(cls.protConvert)
        return cls.protConvert

    @classmethod
    def _runPhantomSubtomo_NoMW(cls, protocolName, option, inputVolume, sampling, nsubtomos, mwfilter,
                                  rotate, applyShift, coords, addNoise):

        cls.phantom_NoMW = cls.newProtocol(protocolName, option=option, inputVolume=inputVolume, sampling=sampling,
                                nsubtomos= nsubtomos, mwfilter= mwfilter, rotate =rotate,
                                      applyShift=applyShift, coords=coords, addNoise=addNoise)
        cls.launchProtocol(cls.phantom_NoMW)

        return cls.phantom_NoMW

    @classmethod
    def _runPhantomSubtomo_NoMW_noisy(cls, protocolName, option, inputVolume, sampling, nsubtomos, mwfilter,
                                  rotate, applyShift, coords, addNoise, differentStatistics, minstd, maxStd):

        cls.phantom_NoMW_noisy = cls.newProtocol(protocolName, option=option, inputVolume=inputVolume, sampling=sampling,
                                nsubtomos= nsubtomos, mwfilter= mwfilter, rotate =rotate,
                                applyShift=applyShift, coords=coords, addNoise=addNoise,
                                differentStatistics=differentStatistics, minstd=minstd, maxStd=maxStd)
        cls.launchProtocol(cls.phantom_NoMW_noisy)

        return cls.phantom_NoMW_noisy

    @classmethod
    def _runPhantomSubtomo_MW_noisy(cls, protocolName, option, inputVolume, sampling, nsubtomos, mwfilter, mwangle,
                                  rotate, applyShift, coords, addNoise, differentStatistics, minstd, maxStd):

        cls.phantom_MW_noisy = cls.newProtocol(protocolName, option=option, inputVolume=inputVolume, sampling=sampling,
                                nsubtomos= nsubtomos, mwfilter= mwfilter, mwangle=mwangle, rotate =rotate,
                                applyShift=applyShift, coords=coords, addNoise=addNoise,
                                differentStatistics=differentStatistics, minstd=minstd, maxStd=maxStd)
        cls.launchProtocol(cls.phantom_MW_noisy)

        return cls.phantom_MW_noisy



    @classmethod
    def _runPhantomSubtomo_MW_noisy_Randomrotation_shift(cls, protocolName, option, inputVolume, sampling, nsubtomos, mwfilter, mwangle,
                                rotate, uniformAngularDistribution, applyShift, xmin, xmax, ymin, ymax, zmin, zmax,
                                coords, addNoise, differentStatistics, minstd, maxStd):

        cls.phantom_MW_noisy_Randomrotation_shift = cls.newProtocol(protocolName, option=option, inputVolume=inputVolume,
                                sampling=sampling, nsubtomos= nsubtomos, mwfilter= mwfilter, mwangle=mwangle, rotate =rotate,
                                uniformAngularDistribution = uniformAngularDistribution, applyShift=applyShift,
                                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax, coords=coords,
                                addNoise=addNoise, differentStatistics=differentStatistics, minstd=minstd, maxStd=maxStd)
        cls.launchProtocol(cls.phantom_MW_noisy_Randomrotation_shift)

        return cls.phantom_MW_noisy_Randomrotation_shift


    def checkResults(self, outputObject):
        self.assertIsNotNone(outputObject.outputSubtomograms,
                            "There was a problem with the phantom subtomograms created from a given volume")
        self.assertAlmostEqual(outputObject.outputSubtomograms.getSamplingRate(), outputObject.sampling.get(), places=1, msg="Problem with the sampling")
        self.assertEqual(outputObject.outputSubtomograms.getSize(), outputObject.nsubtomos.get(), msg=("Problem with the number of subtomograms"))

    def _geometricalphantom(self):
        phantom = self.newProtocol(XmippProtPhantomSubtomo,
                                   option=1)
        self.launchProtocol(phantom)
        self.assertIsNotNone(phantom.outputSubtomograms,
                             "There was a problem with the phantom subtomograms created from a given volume")
        return phantom

    def _geometricalphantomMW(self):
        geometricalphantomMW = self.newProtocol(XmippProtPhantomSubtomo,
                                                option=1,
                                                rotate=True)
        self.launchProtocol(geometricalphantomMW)
        self.assertIsNotNone(geometricalphantomMW.outputSubtomograms,
                             "There was a problem with subtomograms output")
        return geometricalphantomMW

    def test_VolumeFromPDB(cls):
        cls.assertIsNotNone(cls.protConvert.outputVolume.getFileName(), "There was a problem with the conversion")
        cls.assertAlmostEqual(cls.protConvert.outputVolume.getSamplingRate(), cls.protConvert.sampling.get(), places=1,
                               msg=("Problem with the sampling", "volume"))
        cls.assertEqual(cls.protConvert.outputVolume.getDim()[0], cls.protConvert.size_z.get(),
                               msg=("Problem with the dimensions", "volume"))

    def test_PhantomSubtomos_noMW(cls):
        cls.checkResults(cls.phantom_NoMW)

    def test_PhantomSubtomos_NoMW_noisy(cls):
        cls.checkResults(cls.phantom_NoMW_noisy)

    def test_PhantomSubtomos_MW_noisy(cls):
        cls.checkResults(cls.phantom_MW_noisy)

    def test_PhantomSubtomos_MW_noisy(cls):
        cls.checkResults(cls.phantom_MW_noisy_Randomrotation_shift)


    def test_geometricalphantomMW(self):
        geometricalphantomMW = self._geometricalphantom()
        self.assertTrue(getattr(geometricalphantomMW, 'outputSubtomograms'))
        #self.assertEqual(geometricalphantomMW.outputSubtomograms.getFirstItem().getAcquisition().getAngleMax(), 60)
        #self.assertEqual(geometricalphantomMW.outputSubtomograms.getFirstItem().getAcquisition().getAngleMin(), -60)
        return geometricalphantomMW

    def test_geometricalphantom(self):
        phantom = self._geometricalphantom()
        self.assertTrue(getattr(phantom, 'outputSubtomograms'))
        self.assertEqual(phantom.outputSubtomograms.getFirstItem().getAcquisition().getAngleMax(), 90)
        return phantom



