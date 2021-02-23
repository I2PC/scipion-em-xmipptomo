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

    # --------------------------- STEP functions --------------------------------
    def peakHighContrastStep(self, volId):
        vol = self.inputSetOfVolumes.get()[volId]

        volFileName = vol.getFileName()

        outputFile = os.path.split(os.path.splitext(volFileName)[0])[0] + "xmd"

        print(outputFile)

        paramsPeakHighContrast = {
            'inputVol': volFileName,
            'output': outputFile,
        }

        argsPeakHighContrast = "--vol %(inputVol)s " \
                               "-o %(output)s "

        self.runJob('xmipp_image_peak_high_contrast', argsPeakHighContrast % paramsPeakHighContrast)
