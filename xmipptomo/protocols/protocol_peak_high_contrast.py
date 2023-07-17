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
from pyworkflow import BETA
from tomo import constants
from tomo.protocols import ProtTomoBase
from tomo.objects import Coordinate3D, SetOfCoordinates3D
from xmipptomo import utils


class XmippProtPeakHighContrast(EMProtocol, ProtTomoBase):
    """
    Wrapper protocol to Xmipp image peak high contrast applied to any volume
    """

    _label = 'peak high contrast'
    _devStatus = BETA
    _possibleOutputs = {"outputSetOfCoordinates3D": SetOfCoordinates3D}

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSetOfTomograms',
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Input set of tomograms',
                      important=True,
                      help='Select a set of volumes to peak high contrast regions.')

        form.addParam('fiducialSize',
                      params.FloatParam,
                      label='Fiducial size (nm)',
                      default='10',
                      important=True,
                      help="Size of the fiducial markers (or any other object) to be peaked in nanometers.")

        form.addParam('boxSize',
                      params.IntParam,
                      label='Box size',
                      default='32',
                      important=True,
                      help="Size of the box containing the high contrast feature in pixels.")

        form.addParam('relaxedMode',
                      params.BooleanParam,
                      default=True,
                      label='Run in relaxed mode?',
                      help="If this option is selected coordinates are kept when none of them pass the mirror "
                           "correlation filter. If not, and empty output is possible. This second case might happen "
                           "if the tomogram does not present any gold bead or if it presents misalignment")

        form.addParam('relaxedModeThr',
                      params.IntParam,
                      label='Relaxed mode threshold',
                      default='3',
                      condition='relaxedMode==True',
                      help="Minimum number of surviving coordinates to enter in relaxed mode.")

        # Advanced parameters
        form.addParam('numberSampSlices',
                      params.IntParam,
                      label='Number of sampling slices',
                      default='10',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Number of slices used as a sample to calculate the threshold pixel value, for posterior "
                           "high contrast regions detection.")

        form.addParam('sdThr',
                      params.FloatParam,
                      label='Threshold for initial coordinates (SD)',
                      default='3',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Number of standard deviations (SD) that a coordinate value must be over the mean in other "
                           "to consider it a member of a high contrast feature.")

        form.addParam('numberOfCoordinatesThr',
                      params.IntParam,
                      label='Number of coordinates threshold',
                      default='50',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Number of coordinates that must be attracted by a center of mass to consider it a "
                           "plausible high contrast feature.")

        form.addParam('mirrorCorrelationThr',
                      params.FloatParam,
                      label='Minimum mirror correlation',
                      default='0.1',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Minimum correlation between a feature and its mirror to consider it a fiducial.")

        form.addParam('mahalanobisDistanceThr',
                      params.FloatParam,
                      label='Mahalanobis distance threshold',
                      default='2',
                      expertLevel=params.LEVEL_ADVANCED,
                      help="Maximum Mahalanobis distance of the radial average of the gold bead between all the "
                           "peaked coordinates.")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        for vol in self.inputSetOfTomograms.get():
            peakId = self._insertFunctionStep('peakHighContrastStep',
                                              vol.getObjId(),
                                              prerequisites=[])

            self._insertFunctionStep('createOutputStep',
                                     vol.getObjId(),
                                     prerequisites=[peakId])

        self._insertFunctionStep('closeOutputSetStep')

    # --------------------------- STEP functions --------------------------------
    def peakHighContrastStep(self, volId):
        vol = self.inputSetOfTomograms.get()[volId]

        inputFilePath = vol.getFileName()
        outputFileName = os.path.splitext(os.path.split(inputFilePath)[1])[0] + ".xmd"
        outputFilePath = os.path.join(self._getExtraPath(), outputFileName)

        paramsPeakHighContrast = {
            'inputVol': inputFilePath + ":mrc",
            'output': outputFilePath,
            'boxSize': self.boxSize.get(),
            'fiducialSize': self.fiducialSize.get() * 10,
            'sdThr': self.sdThr.get(),
            'mirrorCorrelationThr': self.mirrorCorrelationThr.get(),
            'numberSampSlices': self.numberSampSlices.get(),
            'numberOfCoordinatesThr': self.numberOfCoordinatesThr.get(),
            'samplingRate': self.inputSetOfTomograms.get().getSamplingRate(),
            'mahalanobisDistanceThr': self.mahalanobisDistanceThr.get(),
        }

        argsPeakHighContrast = "--vol %(inputVol)s " \
                               "-o %(output)s " \
                               "--boxSize %(boxSize)d " \
                               "--fiducialSize %(fiducialSize)f " \
                               "--sdThr %(sdThr)f " \
                               "--mirrorCorrelationThr %(mirrorCorrelationThr)f " \
                               "--mahalanobisDistanceThr %(mahalanobisDistanceThr)f " \
                               "--numberSampSlices %(numberSampSlices)d " \
                               "--numberOfCoordinatesThr %(numberOfCoordinatesThr)s " \
                               "--samplingRate %(samplingRate)f "

        if self.relaxedMode.get():
            paramsPeakHighContrast['relaxedModeThr'] = self.relaxedModeThr.get()

            argsPeakHighContrast += "--relaxedModeThr %(relaxedModeThr)d "

        self.runJob('xmipp_image_peak_high_contrast', argsPeakHighContrast % paramsPeakHighContrast)

    def createOutputStep(self, volId):
        vol = self.inputSetOfTomograms.get()[volId]
        volObjId = vol.getObjId()

        volFileName = vol.getFileName()
        outputFileName = os.path.splitext(os.path.split(volFileName)[1])[0] + ".xmd"
        outputFilePath = os.path.join(self._getExtraPath(), outputFileName)

        outputSetOfCoordinates3D = self.getOutputSetOfCoordinates3Ds()

        coordList = utils.retrieveXmipp3dCoordinatesIntoList(outputFilePath)

        if not coordList:
            print("WARNING: no coordinates picked in tomogram " + volFileName)

        for element in coordList:
            newCoord3D = Coordinate3D()
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
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=self.inputSetOfTomograms.get(),
                                                                      suffix='')
            outputSetOfCoordinates3D.setSamplingRate(self.inputSetOfTomograms.get().getSamplingRate())
            outputSetOfCoordinates3D.setPrecedents(self.inputSetOfTomograms)
            outputSetOfCoordinates3D.setBoxSize(self.boxSize.get())
            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
            self._defineSourceRelation(self.inputSetOfTomograms, outputSetOfCoordinates3D)
        return self.outputSetOfCoordinates3D

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfCoordinates3D'):
            summary.append("High contrast features found in %d volumes: %d."
                           % (self.inputSetOfTomograms.get().getSize(),
                              self.outputSetOfCoordinates3D.getSize()))

        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfCoordinates3D'):
            methods.append("%d high contrast features have been found in %d volumes."
                           % (self.outputSetOfCoordinates3D.getSize(),
                              self.inputSetOfTomograms.get().getSize()))

        else:
            methods.append("Output classes not ready yet.")
        return methods
