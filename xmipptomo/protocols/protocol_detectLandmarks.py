# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
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
from pyworkflow import BETA
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path
from pyworkflow.object import Set, List, String, Boolean
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
import xmipptomo.utils as utils
import emtable

OUTPUT_COORDS_FILENAME = "outputLandmarkCoordinates.xmd"


class XmippProtDetectLandmarkTS(EMProtocol, ProtTomoBase):
    """
    Scipion protocol for xmipp_tomo_detect_landmarks. Detect landmarks in a tilt series.
    """

    _label = 'detect landmarks TS'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtTomoBase.__init__(self)
        self.stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('fiducialSize',
                      params.FloatParam,
                      important=True,
                      label='Landmark size (nm)',
                      help='Landmark size in nanometers (nm).')

        # Advanced params
        form.addParam('thrSD',
                      params.FloatParam,
                      advanced=True,
                      default=5,
                      label='Coordinate value SD threshold',
                      help='Number of SD a coordinate value must be over the mean to consider that it belongs to a '
                           'high contrast feature.')

        form.addParam('targetLMsize',
                      params.FloatParam,
                      advanced=True,
                      default=8,
                      label='Target landmark size (px)',
                      help='Target landmark size to adjust down sampling before filtering. Default value is 8 px.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        allcosId = []

        for ts in self.inputSetOfTiltSeries.get():
            tsObjId = ts.getObjId()
            cisID = self._insertFunctionStep(self.convertInputStep,
                                             tsObjId,
                                             prerequisites=[])
            dlsID = self._insertFunctionStep(self.convertInputStep,
                                             tsObjId,
                                             prerequisites=[cisID])
            cosID = self._insertFunctionStep(self.createOutputStep,
                                             tsObjId,
                                             prerequisites=[dlsID])
            allcosId.append(cosID)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=allcosId)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        firstItem = ts.getFirstItem()

        """Apply the transformation form the input tilt-series"""
        # Use Xmipp interpolation via Scipion
        if firstItem.hasTransform():
            avgRotAngle = utils.calculateAverageRotationAngleFromTM(ts)
            swap = True if (avgRotAngle > 45 or avgRotAngle < -45) else False

            outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
            ts.applyTransform(outputTsFileName, swapXY=swap)

        else:
            outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
            ts.applyTransform(outputTsFileName)

    def detectLandmarksStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()

        params = {
            'i': os.path.join(tmpPrefix, firstItem.parseFileName() + ":mrcs"),
            'o': os.path.join(extraPrefix, OUTPUT_COORDS_FILENAME),
            'thrSD': self.thrSD.get(),
            'samplingRate': self.inputSetOfTiltSeries.get().getSamplingRate(),
            'fiducialSize': self.fiducialSize.get() * 10,
            'targetLMsize': self.targetLMsize.get()
        }

        args = "-i %(i)s " \
               "--tlt %(tlt)s " \
               "--inputCoord %(inputCoord)s " \
               "-o %(o)s " \
               "--thrSDHCC %(thrSDHCC).2f " \
               "--thrNumberCoords %(thrNumberCoords).2f " \
               "--samplingRate %(samplingRate).2f " \
               "--fiducialSize %(fiducialSize).2f " \
               "--thrChainDistanceAng %(thrChainDistanceAng).2f " \
               "--thrFiducialDistance %(thrFiducialDistance).2f " \
               "--avgResidPercentile_LocalAlignment %(avgResidPercentile_LocalAlignment).4f "

        self.runJob('xmipp_tomo_detect_landmarks', args % params)

    def generateOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        self.getOutputSetOfLandmarkModels()

        lmFileName = os.path.join(extraPrefix, firstItem.parseFileName(suffix='_lm', extension='.sfid'))
        lm = tomoObj.LandmarkModel(tsId=tsId,
                                   fileName=lmFileName,
                                   modelName=None,
                                   size=self.fiducialSize.get() * 10,
                                   applyTSTransformation=Boolean(False))
        lm.setTiltSeries(ts)

        lmList = self.parseLandmarkCoordinatesFile(vcmFilePath=os.path.join(extraPrefix, OUTPUT_COORDS_FILENAME))

        for i, lmInfo in enumerate(lmList):
            lm.addLandmark(xCoor=lmInfo[0],
                           yCoor=lmInfo[1],
                           tiltIm=lmInfo[2],
                           chainId=i,
                           xResid=0.0,
                           yResid=0.0)

        self._store()

    def closeOutputSetsStep(self):
        if hasattr(self, "outputSetOfLandmarkModels"):
            self.outputSetOfLandmarkModels.setStreamState(Set.STREAM_CLOSED)
            self.outputSetOfLandmarkModels.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getOutputSetOfLandmarkModels(self):
        if hasattr(self, "outputSetOfLandmarkModels"):
            self.outputSetOfLandmarkModels.setStreamState(Set.STREAM_OPEN)
            self.outputSetOfLandmarkModels.enableAppend()

        else:
            outputSetOfLandmarkModels = self._createSetOfLandmarkModels("_ali")

            outputSetOfLandmarkModels.copyInfo(self.inputSetOfTiltSeries.get())

            outputSetOfLandmarkModels.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfLandmarkModels=outputSetOfLandmarkModels)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfLandmarkModels)

        return self.outputSetOfLandmarkModels

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []

        return methods
