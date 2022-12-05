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
import enum

import numpy as np
from pwem.emlib.image import ImageHandler as ih
from pwem.emlib import lib
from pwem.objects import Particle, Volume, Transform, String, SetOfVolumes, SetOfParticles
from pwem.protocols import ProtAnalysis3D
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, EnumParam, IntParam, BooleanParam
from tomo.objects import SetOfTomograms, SetOfSubTomograms, SetOfCoordinates3D, TomoAcquisition, MATRIX_CONVERSION

import pwem.emlib as emlib
import pwem.emlib.metadata as md
from xmipptomo import utils
from xmipp3.convert import alignmentToRow
import tomo.constants as const


class SubtomoProjectOutput(enum.Enum):
    particles = SetOfParticles
    average = Particle

class XmippProtExtractSubtomos(ProtAnalysis3D):
    """
    Extract a set of subtomograms from a tomogram given their coordinates.
    """
    _label = 'extract subtomos'
    _devStatus = BETA
    _possibleOutputs = SubtomoProjectOutput

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Parameters')
        form.addParam('tomograms', PointerParam, pointerClass=SetOfTomograms,
                      label='Tomograms', help='This protocol can *not* work with .em files *if* the input is a set'
                                                  ' of tomograms or a set of volumes, ')

        form.addParam('coords', PointerParam, pointerClass=SetOfCoordinates3D,
                  label='Coordinates', help='This protocol can *not* work with .em files *if* the input is a set'
                                              ' of tomograms or a set of volumes, ')

        form.addParam('boxSize', IntParam,
                      label='Box size',
                      help='The subtomograms are extracted as a cubic box of this size. '
                           'The wizard selects same box size as picking')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):

        for tom in self.tomograms.get():
            tomId = tom.getObjId()
            self._insertFunctionStep(self.extractStep, self.tomograms.get(), tomId)
            #self._insertFunctionStep(self.createHistrogram, tomId)
        #self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -------------------------------
    '''
    def writeSetOfCoordinates3D(self):
        tomoList = []
        for tomo in self.getInputTomograms():
            tomoList.append(tomo.clone())

        for tomo in tomoList:
            coordDict = []
            self.coordsFileName = self._getExtraPath(pwutils.replaceBaseExt(tomo.getFileName(), 'coords'))

            with open(self.coordsFileName, "w") as out:
                coords = self.inputCoordinates.get()
                for coord3D in coords.iterCoordinates(volume=tomo):
                    if os.path.basename(tomo.getFileName()) == os.path.basename(coord3D.getVolName()):
                        out.write("%d\t%d\t%d\n" % (coord3D.getX(const.BOTTOM_LEFT_CORNER),
                                                    coord3D.getY(const.BOTTOM_LEFT_CORNER),
                                                    coord3D.getZ(const.BOTTOM_LEFT_CORNER)))
                        newCoord = coord3D.clone()
                        newCoord.setVolume(coord3D.getVolume())
                        coordDict.append(newCoord)

            if coordDict:
                self.lines.append(coordDict)
                self.tomoFiles.append(tomo.getFileName())
                self.samplingRateTomo = tomo.getSamplingRate()
    '''

    def writeMdCoordinates(self, tomo):
        mdCoor = lib.MetaData()

        tsid = tomo.getTsId()

        for item in self.coords.get().iterItems():
            coord = item
            transform = Transform(matrix=item.getMatrix(convention=MATRIX_CONVERSION.XMIPP))

            if coord.getTomoId() == tsid:
                nRow = md.Row()
                nRow.setValue(lib.MDL_ITEM_ID, int(coord.getObjId()))
                coord.setVolume(tomo)
                nRow.setValue(lib.MDL_XCOOR, int(coord.getX(const.BOTTOM_LEFT_CORNER)) )
                nRow.setValue(lib.MDL_YCOOR, int(coord.getY(const.BOTTOM_LEFT_CORNER)) )
                nRow.setValue(lib.MDL_ZCOOR, int(coord.getZ(const.BOTTOM_LEFT_CORNER)) )
                # Compute inverse matrix
                # A = subtomo.getTransform().getMatrix()
                # subtomo.getTransform().setMatrix(np.linalg.inv(A))
                # Convert transform matrix to Euler Angles (rot, tilt, psi)
                from pwem import ALIGN_PROJ
                alignmentToRow(transform, nRow, ALIGN_PROJ)
            nRow.addToMd(mdCoor)

        fnCoor = self._getExtraPath("coords%s.xmd" % tsid)
        mdCoor.write(fnCoor)



    def extractStep(self, inTomograms, tomId):


        print('----------------------')

        self.writeMdCoordinates(self, inTomograms[tomId])
        print('----------------------')

        '''
        import os
        tomoFn = inTomograms[tomId].getFileName()

        ts = inTomograms[tomId]
        tsId = ts.getTsId()

        # Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        os.mkdir(tomoPath)

        fnPath = os.path.join(tomoPath, str(tomId) + ext)

        # Defining outfiles

        fullTomogramName = self.createOutputPath(FULL_TOMOGRAM_FILE, tsId, MRCEXT)

        input = self.tomograms.get()
        x, y, z = input.getDim()
        dir = self.dirParam.get()
        if self.rangeParam.get() == 1:
            cropParam = self.cropParam.get()

        fnProj = self._getExtraPath("projections.mrcs")
        lib.createEmptyFile(fnProj, x, y, 1, input.getSize())

        for subtomo in input.iterItems():
            fn = subtomo.getLocation()
            if fn[1].endswith('.mrc'):
                fn = list(fn)
                fn[1] += ':mrc'
                fn = tuple(fn)
                subtomo.setFileName(fn[1])
            vol = Volume()
            vol.setLocation('%d@%s' % fn)
            vol = ih().read(vol.getLocation())
            img = ih().createImage()
            if self.radAvg.get():
                img = vol.radialAverageAxis()
            else:
                volData = vol.getData()
                proj = np.empty([x, y])
                if dir == 0:
                    if self.rangeParam.get() == 1:
                        volData = volData[:, :, int(x/2 - cropParam):int(x/2 + cropParam):1]
                    for zi in range(z):
                        for yi in range(y):
                            proj[zi, yi] = np.sum(volData[zi, yi, :])
                img.setData(proj)

            img.write('%d@%s' % (subtomo.getObjId(), fnProj))

        '''

    def createOutputStep(self):
        self._defineOutputs(outputSubtomograms=self.outputSet)
        self._defineSourceRelation(self.tomos.get(), self.outputSet)

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        vols = self.input.get()
        return ["Projection of %d volumes with dimensions %s obtained."
                % (vols.getSize(), vols.getDimensions())]

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Output views not ready yet.")

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                summary.append("*%s:* \n %s " % (key, output.getObjComment()))
        return summary

    def getSummary(self, imgSetOut):
        summary = []
        summary.append("Number of projections generated: %s" % imgSetOut.getSize())
        return "\n".join(summary)
