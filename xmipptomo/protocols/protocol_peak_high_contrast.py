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

from pwem.protocols import EMProtocol
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as const
from tomo.protocols import ProtTomoBase
import os


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

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        for vol in self.inputSetOfVolumes.get():
            self._insertFunctionStep('peakHighContrastStep', vol.getObjId())
            self._insertFunctionStep('createOutputStep', vol.getObjId())

    # --------------------------- STEP functions --------------------------------
    def peakHighContrastStep(self, volId):
        vol = self.inputSetOfVolumes.get()[volId]

        volFileName = vol.getFileName()
        outputFileName = os.path.split(os.path.splitext(volFileName)[0])[0] + "xmd"
        outputFilePath = os.path.join(self._getExtraPath(), outputFileName)

        paramsPeakHighContrast = {
            'inputVol': volFileName,
            'output': outputFilePath,
        }

        argsPeakHighContrast = "--vol %(inputVol)s " \
                               "-o %(output)s "

        self.runJob('xmipp_image_peak_high_contrast', argsPeakHighContrast % paramsPeakHighContrast)

    def createOutputStep(self, volId):
        vol = self.inputSetOfVolumes.get()[volId]

        volFileName = vol.getFileName()
        outputFileName = os.path.split(os.path.splitext(volFileName)[0])[0] + "xmd"
        outputFilePath = os.path.join(self._getExtraPath(), outputFileName)

        outputSetOfCoordinates3D = self.getOutputSetOfCoordinates3Ds()

        xVol, yVol, zVol = vol.getDim()

        for element in coordList:
            newCoord3D = tomoObj.Coordinate3D(x=element[0],
                                              y=element[1],
                                              z=element[2])
            newCoord3D.setVolume(ts)
            newCoord3D.setVolId(tsObjId)
            outputSetOfCoordinates3D.append(newCoord3D)
            outputSetOfCoordinates3D.update(newCoord3D)
        outputSetOfCoordinates3D.write()
        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfCoordinates3Ds(self):
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.outputSetOfCoordinates3D.enableAppend()
        else:
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=self.getOutputSetOfTiltSeries(),
                                                                      suffix='LandmarkModel')
            outputSetOfCoordinates3D.setSamplingRate(self.inputSetOfTiltSeries.get().getSamplingRate())
            outputSetOfCoordinates3D.setPrecedents(self.inputSetOfTiltSeries)
            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfCoordinates3D)
        return self.outputSetOfCoordinates3D
