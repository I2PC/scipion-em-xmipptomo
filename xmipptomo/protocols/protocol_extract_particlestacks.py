# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
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

import os
import glob
from pwem.emlib import lib
from pwem.objects import Particle, Volume, Transform, String, SetOfVolumes, SetOfParticles, ImageDim
from pwem.protocols import EMProtocol
from pwem import ALIGN_PROJ
import pwem.emlib.metadata as md
import pwem.emlib.image.image_handler as ih

from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, IntParam, BooleanParam

from tomo.objects import SetOfCTFTomoSeries, SetOfTiltSeries, SetOfCoordinates3D, TomoAcquisition, MATRIX_CONVERSION
from ..objects import SetOfTiltSeriesParticle, TiltSeriesParticle
from ..utils import writeMdTiltSeries, calculateRotationAngleAndShiftsFromTM
import tomo.constants as const

from tomo.protocols import ProtTomoBase
from xmipp3.convert import alignmentToRow, readSetOfParticles

COORD_BASE_FN = 'coords'
EXT_XMD = '.xmd'

# Tomogram type constants for particle extraction
OUTPUTATTRIBUTE = 'TiltSeriesParticle'


class XmippProtExtractParticleStacks(EMProtocol, ProtTomoBase):
    """
    Extract a set of particle stacks from a set of tilt series given a set of coordinates in the tomogram.
    """
    _label = 'extract particle stacks'
    _devStatus = BETA
    _possibleOutputs = {OUTPUTATTRIBUTE: SetOfTiltSeriesParticle}
    lines = []
    tomoFiles = []
    listofTsId = []

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Parameters')

        form.addParam('coords', PointerParam, pointerClass=SetOfCoordinates3D,
                      label='Coordinates',
                      help='3D coordinates to use in the extraction process.'
                           'The coordinate denotes the center of the particle stack '
                           '(subtomogram position)',
                      important=True)

        form.addParam('tiltseries', PointerParam, pointerClass=SetOfTiltSeries,
                      allowsNull=True,
                      label='Tilt Series', help='The particle stacks will be extracted from this set.',
                      important=True)

        form.addParam('setCTFinfo', BooleanParam,
                      default=True,
                      label='Associate CTF?', help='Set of CTF associated to the used tilt series. Please'
                                                   'tae into accourt the alignment information XXXXX',
                      important=True)

        form.addParam('ctfTomoSeries', PointerParam, pointerClass=SetOfCTFTomoSeries,
                      condition='setCTFinfo',
                      label='Set Of CTFs ', help='Set of CTF associated to the used tilt series. Please'
                                                 'tae into accourt the alignment information XXXXX',
                      important=True)

        form.addParam('boxSize', IntParam,
                      label='Box size',
                      help='The particle stack are extracted as squares. The box size defines the edge of the square',
                      important=True)

        form.addParam('invertContrast', BooleanParam,
                      label='Invert Contrast', default=True,
                      help='Normally, tilt series has the contrast inverted with respect to standard for tilt series particles. '
                           'It means, in the tomograms/recosntructions the structure is black and the noise is white. Generally, the'
                           'tilt series particles are white with a black background. This means that the tilt series particles have the '
                           'contrast inverted with respect to the tilt series. Put this flag as True if the contrast'
                           'need to be inverted.')
        '''
        form.addParam('swapXY', BooleanParam,
                      label='swap X - Y coordinates?', default=True,
                      help='The dimensions along the the axis X and Y must be equal to the dimenions along the same axis of the tilt series.'
                           'Sometimes the tilt series or the tomograms present a rotation of 90º. This means that the X axis of the tomogram'
                           'corresponds to the Y axis of the tilt series. If this is the case, this option must be set as true. \n '
                           '\n Example 1:\n'
                           'Tomogram dimensions (1024x1000x300) \n'
                           'Tilt series dimensions (1000x1024) \n'
                           'In this example the swap option must be true, because the dimensions of Tilt series and the tomogram are different.\n'
                           '\n Example 2:\n'
                           'Tomogram dimensions (1024x1000x300)'
                           'Tilt series dimensions (1024x1000)'
                           'In this example the swap option must be False, because the dimensions of Tilt series and the tomogram are identical '
                      )
        '''

        form.addParam('normalize', BooleanParam, default=True,
                      label='Normalize particles',
                      help='Generally, the noise presents zero mean and standard deviation equal. If true, the particles will be normalize'
                           'to fulfill this property. If false, no action will be carried out')

        form.addParam('asSPAparticles', BooleanParam, default=False,
                      label='Extract as SPA particles',
                      help='False for a pure tilt series particles averaging treatment of the problem. True if the extractted particles will be used in'
                           'a single Particle Analysis workflow.')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self.checkSwapXY()
        for ts in self.tiltseries.get():
            self._insertFunctionStep(self.extractStackStep, ts.getObjId())
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -------------------------------

    def writeMdCoordinates(self, ts, tomoPath):
        """
            Returns the filename of a metadata with the coordinates.
        """
        mdCoor = lib.MetaData()

        tsid = ts.getTsId()
        coordDict = []

        for item in self.coords.get().iterCoordinates():
            coord = item
            transform = Transform(matrix=item.getMatrix(convention=MATRIX_CONVERSION.XMIPP))

            if coord.getTomoId() == tsid:
                nRow = md.Row()
                nRow.setValue(lib.MDL_ITEM_ID, int(coord.getObjId()))
                nRow.setValue(lib.MDL_XCOOR, int(coord.getX(const.SCIPION)))
                nRow.setValue(lib.MDL_YCOOR, int(coord.getY(const.SCIPION)))
                nRow.setValue(lib.MDL_ZCOOR, int(coord.getZ(const.SCIPION)))

                alignmentToRow(transform, nRow, ALIGN_PROJ)
                nRow.addToMd(mdCoor)

                newCoord = item.clone()
                newCoord.setVolume(coord.getVolume())
                coordDict.append(newCoord)
                self.lines.append(coordDict)
                self.listofTsId.append(tsid)

        fnCoor = os.path.join(tomoPath, "coor_%s.xmd" % tsid)
        mdCoor.write(fnCoor)

        return fnCoor

    def checkSwapXY(self):
        firstTomo = self.coords.get().getPrecedents().getFirstItem()
        tomoDims = firstTomo.getDimensions()
        TSdims = self.tiltseries.get().getDimensions()

        swap = False
        if (tomoDims[0] == TSdims[1]) and (tomoDims[1] == TSdims[0]):
            swap = True

        return swap


    def writeMdTiltSeriesWithCTF(self, ts, tomoPath):
        """
        Returns a metadata with the tilt series information, TsID, filename and tilt angle.
        """
        doseValueInTS = []
        tsList = []
        tsId = ts.getTsId()
        for ti in ts:
            doseValueInTS.append(ti.getAcquisition().getAccumDose())
            tsList.append(ti.clone())

        sorted_Idx_doseValueInTS = sorted(range(len(doseValueInTS)), key=lambda k: doseValueInTS[k])

        doseValueInTSCTF = []
        ctfTomo = self.getCTFfromId(self.ctfTomoSeries.get(), tsId)
        ctfList = []
        tiltSeriesFromCTF = ctfTomo.getTiltSeries()

        for ctf in ctfTomo:
            ctfList.append(ctf.clone())
            index = ctf.getIndex().get()
            ti = tiltSeriesFromCTF[index]
            doseValueInTSCTF.append(ti.getAcquisition().getAccumDose())

        sorted_doseValueInCTF = sorted(range(len(doseValueInTS)), key=lambda k: doseValueInTSCTF[k])

        warningStr = ''
        mdts = lib.MetaData()
        for idxTS, idxCTF in zip(sorted_Idx_doseValueInTS, sorted_doseValueInCTF):
            if doseValueInTS[idxTS] != doseValueInTSCTF[idxCTF]:
                warningStr += '%s \n' % tsId
                break
            else:
                ti = tsList[idxTS]
                ctfti = ctfList[idxCTF]
                transform = ti.getTransform()
                if transform is None:
                    rot = 0
                    Sx = 0
                    Sy = 0
                else:
                    rot, Sx, Sy = calculateRotationAngleAndShiftsFromTM(ti)

                tiIndex = ti.getLocation()[0]
                fn = str(tiIndex) + "@" + ti.getFileName()
                nRow = md.Row()
                nRow.setValue(lib.MDL_IMAGE, fn)

                defU = ctfti.getDefocusU()
                defV = ctfti.getDefocusV()
                defAng = ctfti.getDefocusAngle()
                nRow.setValue(lib.MDL_CTF_DEFOCUSU, defU)
                nRow.setValue(lib.MDL_CTF_DEFOCUSV, defV)
                nRow.setValue(lib.MDL_CTF_DEFOCUS_ANGLE, defAng)

                nRow.setValue(lib.MDL_TSID, tsId)
                tilt = ti.getTiltAngle()
                nRow.setValue(lib.MDL_ANGLE_TILT, tilt)
                nRow.setValue(lib.MDL_ANGLE_ROT, rot)
                nRow.setValue(lib.MDL_SHIFT_X, Sx)
                nRow.setValue(lib.MDL_SHIFT_Y, Sy)
                nRow.addToMd(mdts)

        fnts = os.path.join(tomoPath, "%s_ts.xmd" % tsId)

        mdts.write(fnts)
        return fnts

    def getCTFfromId(self, setofctftomoseries, targetTsId):
        for ctfts in setofctftomoseries:
            ctfId = ctfts.getTsId()
            if targetTsId == ctfId:
                return ctfts

    def extractStackStep(self, objId):
        """ This function executes xmipp_tomo_extract_particlestacks. To do that
            1) It defines the set of tiltseries to be used
            2) A folder where the particlestacks will be stored is created. The name of this folder is the tsId
            3) Launches the xmipp_tomo_extract_particlestacks
        """
        ts = self.tiltseries.get()[objId]

        tsId = ts.getTsId()
        tomoPath = self._getExtraPath(tsId)
        os.mkdir(tomoPath)

        tomoFn = ts.getFileName()

        fnCoords = self.writeMdCoordinates(ts, tomoPath)
        if self.setCTFinfo:
            fnTs = self.writeMdTiltSeriesWithCTF(ts, tomoPath)
        else:
            fnTs = writeMdTiltSeries(ts, tomoPath, fnXmd=tsId + '_ts.xmd')

        params = ' --tiltseries %s' % fnTs
        params += ' --coordinates %s' % fnCoords
        params += ' --boxsize %i' % self.boxSize.get()
        if self.invertContrast.get():
            params += ' --invertContrast'
        if self.normalize.get():
            params += ' --normalize'
        if self.checkSwapXY():
            params += ' --swapXY'
        params += ' --threads %i' % 1
        params += ' -o %s' % tomoPath

        self.runJob('xmipp_tomo_extract_particlestacks', params)

        self.tomoFiles.append(tomoFn)

    def createOutputStep(self):
        """
            This function creates the output of the protocol
        """
        ts = self.tiltseries.get()
        firstItem = ts.getFirstItem()
        acquisitonInfo = firstItem.getAcquisition()
        # TODO: Check the sampling if the tomograms are different than the picked ones
        # TODO: Check the sampling rate if a downsampling option is implemented
        outputSet = None

        self.outputParticleStackSet = SetOfTiltSeriesParticle.create(self._getPath())
        self.outputParticleStackSet.setSamplingRate(ts.getSamplingRate())
        bs = self.boxSize.get()
        self.outputParticleStackSet.setDim(ImageDim(bs, bs, bs))
        self.outputParticleStackSet.setAnglesCount(ts.getAnglesCount())
        if firstItem.getAcquisition():
            acquisition = TomoAcquisition()
            acquisition.setAngleMin(acquisitonInfo.getAngleMin())
            acquisition.setAngleMax(acquisitonInfo.getAngleMax())
            acquisition.setStep(acquisitonInfo.getStep())
            self.outputParticleStackSet.setAcquisition(acquisition)

        counter = 0

        for item in ts.iterItems():
            for ind, tomoFile in enumerate(self.tomoFiles):
                if os.path.basename(tomoFile) == os.path.basename(item.getFileName()):
                    coordSet = self.lines[ind]
                    tsId = item.getTsId()
                    outputSet, counter = self.readSetOfParticleStack(tomoFile,
                                                                     self.outputParticleStackSet,
                                                                     counter, tsId)

        self._defineOutputs(**{OUTPUTATTRIBUTE: outputSet})
        self._defineSourceRelation(self.coords, outputSet)

        if self.asSPAparticles:
            fn = self._getExtraPath('allparticles.xmd')
            outputSet = self._createSetOfParticles()
            self.createMdWithall2DTiltParticles(fn)
            readSetOfParticles(fn, outputSet)
            outputSet.setSamplingRate(self.tiltseries.get().getSamplingRate())

            self._defineOutputs(outputParticles=outputSet)
            self._defineSourceRelation(self.coords, outputSet)
            self._defineSourceRelation(self.tiltseries, outputSet)

    def createMdWithall2DTiltParticles(self, fn):

        mdAllParticles = md.MetaData()
        for ts in self.tiltseries.get():
            tsId = ts.getTsId()

            basemdpath = os.path.join(self._getExtraPath(), tsId)
            fnMd = os.path.abspath(os.path.join(basemdpath, tsId + EXT_XMD))
            if os.path.exists(fnMd):
                auxMd = lib.MetaData(fnMd)
                for row in md.iterRows(auxMd):
                    rowglobal = md.Row()
                    fnImg = row.getValue(md.MDL_IMAGE).split('@')
                    ts_orig = row.getValue(md.MDL_TSID)

                    rowglobal.setValue(md.MDL_IMAGE, fnImg[0] + '@' + os.path.join(basemdpath, fnImg[1]))
                    rowglobal.setValue(md.MDL_TSID, ts_orig)
                    rowglobal.addToMd(mdAllParticles)
        mdAllParticles.write(fn)


    def readSetOfParticleStack(self, tomoFile, outputTsParticlesSet, counter, tsId):
        """
            This function set the corresponing attributes to each tilt series particle. Coordinates and transformation matrix
            The output is the set of tilt series particle
        """
        self.info("Registering tilt series particle for %s" % tomoFile)

        outRegex = os.path.join(self._getExtraPath(tsId), tsId + '-*.mrcs')
        print('.............................')
        print(self._getExtraPath(tsId), tsId + '-*.mrc')
        subtomoFileList = sorted(glob.glob(outRegex))
        print(subtomoFileList)

        for idx, subtomoFile in enumerate(subtomoFileList):
            self.debug("Registering tilt series particle %s - %s" % (counter, subtomoFile))

            tsparticle = TiltSeriesParticle()
            tsparticle.cleanObjId()
            tsparticle.setLocation(subtomoFile)
            tsparticle.setVolName(tsId)
            outputTsParticlesSet.append(tsparticle)
            counter += 1
        return outputTsParticlesSet, counter

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        toms = self.coords.get().getPrecedents()
        return ["A set of %d tilt series particles with dimensions %s was obtained."
                % (toms.getSize(), toms.getDimensions())]

    def _validate(self):
        validateMsgs = []

        if self.setCTFinfo and not self.ctfTomoSeries.get():
            validateMsgs.append("A set of CTF must be provided. \n")

        firstTomo = self.coords.get().getPrecedents().getFirstItem()
        tomoDims = firstTomo.getDimensions()
        TSdims = self.tiltseries.get().getDimensions()
        if (tomoDims[0] != TSdims[0]):
            if (tomoDims[0] != TSdims[1]):
                validateMsgs.append('The tilt series dimensions and the tomogram dimensions do not match.'
                                    'The used tomograms for picking have different dimensions than the '
                                    'tilt series.')
        return validateMsgs

    def _warnings(self):
        warningMsgs = []
        if self.checkSwapXY():
            strBold = '*This message is only informative:*'
            warningMsgs.append(strBold)
            warningMsgs.append('The tilt series and the tomograms used for picking are swapped.'
                               'Scipion will carried out the swap for you. '
                               'The dimensions x and y of the tomogram and the tilt series are swapped'
                               'This means that if the tomogam has dimensions (x, y, z) the tilt '
                               'series has dimensions (y, x). ')
        return warningMsgs

    def _summary(self):
        summary = []
        return summary
