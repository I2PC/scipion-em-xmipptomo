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
from pwem.objects import Particle, Volume, Transform, String, SetOfVolumes, SetOfParticles
from pwem.protocols import EMProtocol
from pwem import ALIGN_PROJ
import pwem.emlib.metadata as md

from pyworkflow import BETA
from pyworkflow import utils as pwutils
from pyworkflow.protocol.params import PointerParam, EnumParam, IntParam, BooleanParam

from tomo.objects import SetOfTomograms, SetOfSubTomograms, SubTomogram, SetOfCoordinates3D, TomoAcquisition, MATRIX_CONVERSION
import tomo.constants as const

from tomo.protocols import ProtTomoBase
from xmipp3.convert import alignmentToRow


COORD_BASE_FN = 'coords'
# Tomogram type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1


class XmippProtExtractSubtomos(EMProtocol, ProtTomoBase):
    """
    Extract a set of subtomograms from a set of tomograms given a set of coordinates.
    """
    _label = 'extract subtomos'
    _devStatus = BETA
    #_possibleOutputs = SetOfSubTomograms
    #OUTPUT_PREFIX = _possibleOutputs.subtomograms.name
    lines = []
    tomoFiles = []

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Parameters')

        form.addParam('tomograms', PointerParam, pointerClass=SetOfTomograms,
                      label='Tomograms', help='The subtomograms will be extracted from this set.')

        form.addParam('coords', PointerParam, pointerClass=SetOfCoordinates3D,
                  label='Coordinates', help='The coordinates that definethe position of the set of subtomograms to be extracted.'
                                            'The coordinate denotes the center of the subtomogram')

        form.addParam('boxSize', IntParam,
                      label='Box size',
                      help='The subtomograms are extracted as a cube. The box size defines the edge of the cube'
                           'The wizard selects same box size as picking')

        form.addParam('invertContrast', BooleanParam,
                      label='Invert Contrast',  default=True,
                      help='Normally, tomograms has the contrast inverted with respect to standard for subtomograms. '
                           'It means, in the tomograms the structure is black and the noise is white. Generally, the'
                           'subtomograms are white with a black background. This means that the subtomograms has the '
                           'contrast inverted with respect to the tomograms. Put this flag as True if the contrast'
                           'need to be inverted.')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):

        for tom in self.tomograms.get():
            tomId = tom.getObjId()
            self._insertFunctionStep(self.extractStep, self.tomograms.get(), tomId)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -------------------------------

    def writeMdCoordinates(self, tomo, tomoPath):
        mdCoor = lib.MetaData()

        tsid = tomo.getTsId()
        coordDict = []

        for item in self.coords.get().iterCoordinates(volume=tomo):
            coord = item
            transform = Transform(matrix=item.getMatrix(convention=MATRIX_CONVERSION.XMIPP))

            if coord.getTomoId() == tsid:

                nRow = md.Row()
                nRow.setValue(lib.MDL_ITEM_ID, int(coord.getObjId()))
                coord.setVolume(tomo)

                nRow.setValue(lib.MDL_XCOOR, int(coord.getX(const.BOTTOM_LEFT_CORNER)))
                nRow.setValue(lib.MDL_YCOOR, int(coord.getY(const.BOTTOM_LEFT_CORNER)))
                nRow.setValue(lib.MDL_ZCOOR, int(coord.getZ(const.BOTTOM_LEFT_CORNER)))

                alignmentToRow(transform, nRow, ALIGN_PROJ)
                nRow.addToMd(mdCoor)

                newCoord = item.clone()
                newCoord.setVolume(coord.getVolume())
                coordDict.append(newCoord)
                self.lines.append(coordDict)

        fnCoor = os.path.join(tomoPath, "%s.xmd" % tsid)
        mdCoor.write(fnCoor)

        return fnCoor


    def extractStep(self, inTomograms, tomId):

        ts = inTomograms[tomId]
        tsId = ts.getTsId()

        # Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        os.mkdir(tomoPath)

        tomoFn = inTomograms[tomId].getFileName()

        fnCoords = self.writeMdCoordinates(inTomograms[tomId], tomoPath)

        params = ' --tomogram %s' % tomoFn
        params += ' --coordinates %s' % fnCoords
        params += ' --boxsize %i' % self.boxSize.get()
        if self.invertContrast.get():
            params += ' --invertContrast'
        params += ' --subtomo'
        params += ' --threads %i' % 1

        params += ' -o %s' % tomoPath
        self.runJob('xmipp_tomo_extract_subtomograms', params)

        self.tomoFiles.append(tomoFn)

    def createOutputStep(self):

        outputSet = None
        self.outputSubTomogramsSet = self._createSetOfSubTomograms(self._getOutputSuffix(SetOfSubTomograms))
        self.outputSubTomogramsSet.setSamplingRate(self.tomograms.get().getSamplingRate())
        self.outputSubTomogramsSet.setCoordinates3D(self.coords)
        if self.tomograms.get().getFirstItem().getAcquisition():
            acquisition = TomoAcquisition()
            acquisition.setAngleMin(self.tomograms.get().getFirstItem().getAcquisition().getAngleMin())
            acquisition.setAngleMax(self.tomograms.get().getFirstItem().getAcquisition().getAngleMax())
            acquisition.setStep(self.tomograms.get().getFirstItem().getAcquisition().getStep())
            self.outputSubTomogramsSet.setAcquisition(acquisition)

        counter = 0

        for item in self.tomograms.get().iterItems():
            for ind, tomoFile in enumerate(self.tomoFiles):
                if os.path.basename(tomoFile) == os.path.basename(item.getFileName()):
                    coordSet = self.lines[ind]
                    tomId = item.getObjId()
                    tsId = self.tomograms.get()[tomId].getTsId()
                    outputSet, counter = self.readSetOfSubTomograms(tomoFile,
                                                                    self.outputSubTomogramsSet,
                                                                    coordSet, 1, counter, tsId)

        self._defineOutputs(outputSetOfSubtomogram=outputSet)
        self._defineSourceRelation(self.coords, outputSet)

    def readSetOfSubTomograms(self, tomoFile, outputSubTomogramsSet, coordSet, factor, counter, tsId):

        self.info("Registering subtomograms for %s" % tomoFile)

        outRegex = os.path.join(self._getExtraPath(tsId), pwutils.removeBaseExt(tomoFile)+'-*.mrc')

        subtomoFileList = sorted(glob.glob(outRegex))

        for idx, subtomoFile in enumerate(subtomoFileList):

            self.debug("Registering subtomogram %s - %s" % (counter, subtomoFile))

            subtomogram = SubTomogram()
            subtomogram.cleanObjId()
            subtomogram.setLocation(subtomoFile)
            subtomogram.setCoordinate3D(coordSet[idx])
            transformation = coordSet[idx]._eulerMatrix
            shift_x, shift_y, shift_z = transformation.getShifts()
            transformation.setShifts(factor * shift_x,
                                     factor * shift_y,
                                     factor * shift_z)
            subtomogram.setTransform(transformation, convention=const.TR_EMAN)
            subtomogram.setVolName(tomoFile)
            outputSubTomogramsSet.append(subtomogram)
            counter += 1
        return outputSubTomogramsSet, counter

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        toms = self.tomograms.get()
        return ["A set of %d subtomograms with dimensions %s was obtained."
                % (toms.getSize(), toms.getDimensions())]

    def _summary(self):
        summary = []
        return summary


