# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
# *              Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es)
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
from pwem.objects import Transform
from pwem.protocols import EMProtocol
from pwem import ALIGN_PROJ
import pwem.emlib.metadata as md

from pyworkflow import BETA
from pyworkflow import utils as pwutils
from pyworkflow.protocol.params import PointerParam, FloatParam, IntParam, BooleanParam, EnumParam

from tomo.objects import SetOfTomograms, SetOfSubTomograms, SubTomogram, SetOfCoordinates3D, TomoAcquisition, \
    MATRIX_CONVERSION
import tomo.constants as const

from tomo.protocols import ProtTomoBase
from xmipp3.convert import alignmentToRow

COORD_BASE_FN = 'coords'

# Tomogram type constants for particle extraction
OUTPUTATTRIBUTE = 'Subtomograms'


class XmippProtExtractSubtomos(EMProtocol, ProtTomoBase):
    """
    Extract a set of subtomograms from a set of tomograms given a set of coordinates.
    """
    _label = 'extract subtomos'
    _devStatus = BETA
    _possibleOutputs = {OUTPUTATTRIBUTE: SetOfSubTomograms}

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Parameters')

        form.addParam('coords',
                      PointerParam,
                      pointerClass=SetOfCoordinates3D,
                      label='Coordinates',
                      help='3D coordinates to be used in the extraction process.'
                           'The coordinate denotes the center of the subtomogram.')

        form.addParam('tomoSource',
                      EnumParam,
                      choices=['same as picking', 'other'],
                      default=0,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Tomogram source',
                      help='By default the subtomograms will be extracted from the tomogram used in the picking step '
                           '( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide a different tomogram to extract from. \n'
                           '*Note*: In the _other_ case, ensure that provided tomogram and coordinates are related ')

        form.addParam('tomograms',
                      PointerParam,
                      pointerClass=SetOfTomograms,
                      condition='tomoSource==1',
                      allowsNull=True,
                      label='Tomograms',
                      help='The subtomograms will be extracted from this set. If this parameter is empty'
                           'the subtomograms will be extracted from the set of tomograms used for picking the '
                           'coordinates')

        form.addParam('boxSize',
                      IntParam,
                      label='Box size (pixels)',
                      help='The subtomograms are extracted as a cube. The box size defines the edge of the cube. '
                           'This is the boxsize measured on the original tomogram, it means without '
                           'downsampling/binning. ')

        form.addParam('dowsamplingFactor',
                      FloatParam,
                      label='Downsampling factor',
                      default=1.0,
                      help='This factor will reduce the size of the extracted subtomograms. The size of the output '
                           'will be size = boxsize/downSamplingFactor. It must be greater than 1.')

        form.addParam('invertContrast',
                      BooleanParam,
                      label='Invert Contrast',
                      default=True,
                      help='Normally, tomograms has the contrast inverted with respect to standard for subtomograms. '
                           'It means, in the tomograms the structure is black and the noise is white. Generally, the'
                           'subtomograms are white with a black background. This means that the subtomograms has the '
                           'contrast inverted with respect to the tomograms. Put this flag as True if the contrast'
                           'need to be inverted.')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):

        tomodict = self.coords.get().getPrecedentsInvolved()
        for key in tomodict.keys():
            tom = tomodict[key]
            tsId = tom.getTsId()
            self._insertFunctionStep(self.extractStep, tsId)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -------------------------------

    def writeMdCoordinates(self, tomo, tomoPath):
        """
            Returns the filename of a metadata with the coordinates.
        """
        mdCoor = lib.MetaData()

        tsid = tomo.getTsId()

        for item in self.coords.get().iterCoordinates(volume=tomo):
            coord = item
            transform = Transform(matrix=item.getMatrix(convention=MATRIX_CONVERSION.XMIPP))

            if coord.getTomoId() == tsid:
                nRow = md.Row()
                coord.setVolume(tomo)
                nRow.setValue(lib.MDL_XCOOR, int(coord.getX(const.BOTTOM_LEFT_CORNER)))
                nRow.setValue(lib.MDL_YCOOR, int(coord.getY(const.BOTTOM_LEFT_CORNER)))
                nRow.setValue(lib.MDL_ZCOOR, int(coord.getZ(const.BOTTOM_LEFT_CORNER)))
                nRow.setValue(lib.MDL_PARTICLE_ID, int(coord.getObjId()))

                alignmentToRow(transform, nRow, ALIGN_PROJ)
                nRow.addToMd(mdCoor)

                newCoord = item.clone()
                newCoord.setVolume(coord.getVolume())

        fnCoor = os.path.join(tomoPath, "%s.xmd" % tsid)
        mdCoor.write(fnCoor)

        return fnCoor

    def getTomograms(self):
        """
            Returns the set of tomograms that will be used for extracting the coordinates.

            If the unique input is the set of coordinates the output will be their corresponding tomograms
            If a set of tomograms is provided in the form, then the output will be such set of tomograms
        """
        if self.tomoSource.get() == 0:
            inTomograms = self.coords.get().getPrecedents()
        else:
            inTomograms = self.tomograms.get()
        return inTomograms

    def extractStep(self, tsId):
        """
            This function executes xmipp_tomo_extract_subtomos. To do that
            1) It defines the set of tomograms to be used (self.getTomograms)
            2) A folder where the subtomograms will be stored is created. The name of this folder is the tsId
            3) Lanches the xmipp_tomo_extract_subtomos

        """

        inTomograms = self.getTomograms()

        tomo = inTomograms[{'_tsId': tsId}]

        if tomo is None:
            self.warning('Tomogram not found for tsId %s' % tsId)
            return

        # Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        os.mkdir(tomoPath)

        tomoFn = tomo.getFileName()

        fnCoords = self.writeMdCoordinates(tomo, tomoPath)

        params = ' --tomogram %s' % tomoFn
        params += ' --coordinates %s' % fnCoords
        params += ' --boxsize %i' % self.boxSize.get()
        if self.invertContrast.get():
            params += ' --invertContrast'
        params += ' --threads %i' % 1
        if self.dowsamplingFactor.get() != 1:
            params += ' --downsample %f' % self.dowsamplingFactor.get()
        params += ' -o %s ' % tomoPath
        self.runJob('xmipp_tomo_extract_subtomograms', params)

    def createOutputStep(self):
        """
            This function creates the output of the protocol
        """

        precedents = self.getTomograms()
        firstItem = precedents.getFirstItem()
        acquisitonInfo = firstItem.getAcquisition()

        scaleFactor = self.dowsamplingFactor.get()
        newSamplingRate = precedents.getSamplingRate() * scaleFactor

        self.outputSubTomogramsSet = self._createSetOfSubTomograms()
        self.outputSubTomogramsSet.setSamplingRate(newSamplingRate)
        self.outputSubTomogramsSet.setCoordinates3D(self.coords)
        if firstItem.getAcquisition():
            acquisition = TomoAcquisition()
            acquisition.setAngleMin(acquisitonInfo.getAngleMin())
            acquisition.setAngleMax(acquisitonInfo.getAngleMax())
            acquisition.setStep(acquisitonInfo.getStep())
            self.outputSubTomogramsSet.setAcquisition(acquisition)

        for tomo in precedents:
            tsId = tomo.getTsId()
            self.writeSetOfSubtomograms(tsId, self.outputSubTomogramsSet, newSamplingRate, scaleFactor)

        self._defineOutputs(**{OUTPUTATTRIBUTE: self.outputSubTomogramsSet})
        self._defineSourceRelation(self.coords, self.outputSubTomogramsSet)

    def writeSetOfSubtomograms(self, tsId, outputSubTomogramsSet, sampling, scaleFactor):
        fnSubtomos = os.path.join(self._getExtraPath(tsId), tsId + '_extracted.xmd')

        mdsubtomos = lib.MetaData(fnSubtomos)
        import numpy as np

        coords = self.coords.get()

        for row in md.iterRows(mdsubtomos):
            subtomo = SubTomogram()
            idx = row.getValue(md.MDL_PARTICLE_ID)
            fn = row.getValue(md.MDL_IMAGE)
            index, filename = fn.split('@')
            subtomo.setIndex(index)
            subtomo.setFileName(os.path.join(self._getExtraPath(tsId), filename))
            subtomo.setSamplingRate(sampling)
            coord = coords[idx]
            subtomo.setCoordinate3D(coord)
            subtomo.setVolName(tsId)
            trMatrix = coord.getMatrix()
            transform = Transform()
            if scaleFactor != 1:
                shifts = np.array([trMatrix[0, 3], trMatrix[1, 3], trMatrix[2, 3]])
                scaledShifts = scaleFactor * shifts
                trMatrix[0, 3] = scaledShifts[0]
                trMatrix[1, 3] = scaledShifts[1]
                trMatrix[2, 3] = scaledShifts[2]
            transform.setMatrix(trMatrix)
            subtomo.setTransform(transform, convention=const.TR_SCIPION)
            outputSubTomogramsSet.append(subtomo)

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        toms = self.coords.get().getPrecedents()
        return ["A set of %d subtomograms with dimensions %s was obtained."
                % (toms.getSize(), toms.getDimensions())]

    def _validate(self):
        errors = []
        inTomograms = self.tomograms.get()
        if inTomograms is not None:
            if abs(self.coords.get().getPrecedents().getSamplingRate() - inTomograms.getSamplingRate()) > 0.01:  # 0.01A/px is an aceptable error
                errors.append("The coordinates has a different sampling rate than the selected tomograms."
                              "Tomograms and coordinates must be at the same scale. Please ensure a matching in the "
                              "sample rate")
        if self.dowsamplingFactor.get() < 1:
            errors.append("Downsampling factor must be greater than 1.")
        return errors

    def _summary(self):
        summary = []
        toms = self.coords.get().getPrecedents()
        summary.append("A set of %d subtomograms was obtained." % (toms.getSize()))
        return summary
