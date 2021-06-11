# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
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

from pwem.protocols import EMProtocol
import pyworkflow.protocol.params as params
from pyworkflow.object import Set
from tomo import constants
from tomo.protocols import ProtTomoBase
import tomo.objects as tomoObj
from xmipptomo import utils


class XmippProtPeakHighContrast(EMProtocol, ProtTomoBase):
    """
    Wrapper protocol to Xmipp image peak high contrast applied to any volume
    """

    RESIZE_SAMPLINGRATE = 0
    RESIZE_DIMENSIONS = 1
    RESIZE_FACTOR = 2
    RESIZE_PYRAMID = 3

    _label = 'peak high contrast'

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSetOfVolumes',
                      params.PointerParam,
                      pointerClass='SetOfVolumes',
                      important=True,
                      label='Input set of volumnes',
                      help='Select a set of volumes to peak high contrast regions.')

        form.addParam('boxSize',
                      params.IntParam,
                      label='Box size',
                      default='32',
                      help="Size of the box containing the high contrast feature in pixels.")

        form.addParam('fiducialSize',
                      params.FloatParam,
                      label='Fiducial size (nm)',
                      default='10',
                      help="Size of the fiducial markers (or any other object) to be peaked in nanometers.")

        form.addParam('sdThreshold',
                      params.FloatParam,
                      label='Threshold for initial coordinates (SD)',
                      default='5',
                      help="Number of standard deviations (SD) that a coordinate value must be over the mean in other "
                           "to consider it a member of a high contrast feature.")

        form.addParam('centerFeatures',
                      params.BooleanParam,
                      label='Center features',
                      default=False,
                      help="Center the peaked high contrast features by sliding the symmetric volume and maximizing "
                           "the correlation .")

        form.addParam('numberSampSlices',
                      params.IntParam,
                      label='Number of sampling slices',
                      default='10',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Number of slices used as a sample to calculate the threshold pixel value, for posterior "
                           "high contrast regions detection.")

        form.addParam('numberCenterOfMass',
                      params.IntParam,
                      label='Number of initial center of mass',
                      default='10',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Number of slices used as a sample to calculate the threshold pixel value, for posterior "
                           "high contrast regions detection.")

        form.addParam('distanceThr',
                      params.FloatParam,
                      label='Center of mass distance threshold',
                      default='10',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Threshold distance to consider that a 3D coordinate belong to a center of mass.")

        form.addParam('numberOfCoordinatesThr',
                      params.IntParam,
                      label='Number of coordinates threshold',
                      default='50',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Number of coordinates that must be attracted by a center of mass to consider it a "
                           "plausible high contrast feature.")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        for vol in self.inputSetOfVolumes.get():
            peakId = self._insertFunctionStep('peakHighContrastStep',
                                              vol.getObjId(),
                                              prerequisites=[])

            self._insertFunctionStep('createOutputStep',
                                     vol.getObjId(),
                                     prerequisites=[peakId])

        self._insertFunctionStep('closeOutputSetStep')

    # --------------------------- STEP functions --------------------------------
    def peakHighContrastStep(self, volId):
        vol = self.inputSetOfVolumes.get()[volId]

        inputFilePath = vol.getFileName()
        outputFileName = os.path.splitext(os.path.split(inputFilePath)[1])[0] + ".xmd"
        outputFilePath = os.path.join(self._getExtraPath(), outputFileName)

        paramsPeakHighContrast = {
            'inputVol': inputFilePath + ":mrc",
            'output': outputFilePath,
            'boxSize': self.boxSize.get(),
            'fiducialSize': self.fiducialSize.get() * 10,
            'sdThreshold': self.sdThreshold.get(),
            'numberSampSlices': self.numberSampSlices.get(),
            'numberCenterOfMass': self.numberCenterOfMass.get(),
            'distanceThr': self.distanceThr.get(),
            'numberOfCoordinatesThr': self.numberOfCoordinatesThr.get(),
            'samplingRate': self.inputSetOfVolumes.get().getSamplingRate(),
        }

        argsPeakHighContrast = "--vol %(inputVol)s " \
                               "-o %(output)s " \
                               "--boxSize %(boxSize)d " \
                               "--fiducialSize %(fiducialSize)f " \
                               "--sdThreshold %(sdThreshold)f " \
                               "--numberSampSlices %(numberSampSlices)d " \
                               "--numberCenterOfMass %(numberCenterOfMass)d " \
                               "--distanceThr %(distanceThr)f " \
                               "--numberOfCoordinatesThr %(numberOfCoordinatesThr)s " \
                               "--samplingRate %(samplingRate)f "

        if(self.centerFeatures.get()):
            argsPeakHighContrast += "--centerFeatures "

        self.runJob('xmipp_image_peak_high_contrast', argsPeakHighContrast % paramsPeakHighContrast)

    def createOutputStep(self, volId):
        vol = self.inputSetOfVolumes.get()[volId]
        volObjId = vol.getObjId()

        volFileName = vol.getFileName()
        outputFileName = os.path.splitext(os.path.split(volFileName)[1])[0] + ".xmd"
        outputFilePath = os.path.join(self._getExtraPath(), outputFileName)

        outputSetOfCoordinates3D = self.getOutputSetOfCoordinates3Ds()

        coordList = utils.retrieveXmipp3dCoordinatesIntoList(outputFilePath)

        for element in coordList:
            newCoord3D = tomoObj.Coordinate3D()
            newCoord3D.setVolume(vol)
            newCoord3D.setX(element[0], constants.BOTTOM_LEFT_CORNER)
            newCoord3D.setY(element[1], constants.BOTTOM_LEFT_CORNER)
            newCoord3D.setZ(element[2], constants.BOTTOM_LEFT_CORNER)

            newCoord3D.setVolId(volObjId)
            outputSetOfCoordinates3D.append(newCoord3D)
            outputSetOfCoordinates3D.update(newCoord3D)

        outputSetOfCoordinates3D.write()
        self._store()

    def closeOutputSetStep(self):
        self.getOutputSetOfCoordinates3Ds().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfCoordinates3Ds(self):
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.outputSetOfCoordinates3D.enableAppend()
        else:
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=self.inputSetOfVolumes.get(),
                                                                      suffix='')
            outputSetOfCoordinates3D.setSamplingRate(self.inputSetOfVolumes.get().getSamplingRate())
            outputSetOfCoordinates3D.setPrecedents(self.inputSetOfVolumes)
            outputSetOfCoordinates3D.setBoxSize(self.boxSize.get())
            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
            self._defineSourceRelation(self.inputSetOfVolumes, outputSetOfCoordinates3D)
        return self.outputSetOfCoordinates3D

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfCoordinates3D'):
            summary.append("High contrast features found in %d volumes: %d."
                           % (self.inputSetOfVolumes.get().getSize(),
                              self.outputSetOfCoordinates3D.getSize()))

        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfCoordinates3D'):
            methods.append("%d high contrast features have been found in %d volumes."
                           % (self.outputSetOfCoordinates3D.getSize(),
                              self.inputSetOfVolumes.get().getSize()))

        else:
            methods.append("Output classes not ready yet.")
        return methods
