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

METADATA_INPUT_COORDINATES = "fiducialCoordinates.xmd"
VRESMOD_FILE_NAME = "vResMod.xmd"


class XmippProtDetectMisalignmentTiltSeries(EMProtocol, ProtTomoBase):
    """
    Scipion protocol for xmipp_tomo_detect_misalignment_trajectory. Detect misalignment in a tilt series.
    """

    _label = 'detect misaligned TS'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtTomoBase.__init__(self)
        self.stepsExecutionMode = STEPS_PARALLEL

        self.alignmentReport = List([])

        # Global variable to check if coordinates have been input for a specific tilt-series
        self.check = True

        self.inputSetOfTiltSeries = None
        self.inputSetOfLandmarkModels = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSet',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries, SetOfLandmarkModels',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('inputSetOfCoordinates',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeriesCoordinates',
                      important=True,
                      label='Input set of coordinates 3D',
                      help='Set of 3D coordinates indicating the position in space of the fiducials. This set should '
                           'be obtained from the previous alignment step of the tilt-series.')

        form.addParam('fiducialSize',
                      params.FloatParam,
                      important=True,
                      label='Fiducial size (nm)',
                      help='Fiducial size in nanometers (nm).')

        # Advanced parameters
        form.addParam('thrSDHCC',
                      params.FloatParam,
                      advanced=True,
                      default=3,
                      label='Coordinate value SD threshold',
                      help='Number of SD a coordinate value must be over the mean to consider that it belongs to a '
                           'high contrast feature.')

        form.addParam('thrFiducialDistance',
                      params.FloatParam,
                      advanced=True,
                      default=0.5,
                      label='Landmark distance threshold',
                      help='Threshold times of fiducial size as maximum distance to consider a match between the 3d '
                           'coordinate projection and the detected fiducial.')

        # Advanced parameters
        form.addParam('targetLMsize',
                      params.FloatParam,
                      default=8,
                      label='Target fiducial size (px)',
                      help='Target fiducial size in pixels to calculate the downsampling when detecting landmarks'
                           'on the tilt series.',
                      expertLevel=params.LEVEL_ADVANCED)

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self.alignmentReport = List([])
        self.alignmentReport.clear()

        allcossId = []

        if isinstance(self.inputSet.get(), tomoObj.SetOfTiltSeries):
            self.inputSetOfTiltSeries = self.inputSet.get()

            for ts in self.inputSetOfTiltSeries:
                tsObjId = ts.getObjId()
                cisID = self._insertFunctionStep(self.convertInputStep,
                                                 tsObjId,
                                                 prerequisites=[])

                crvID = self._insertFunctionStep(self.calculateResidualVectors,
                                                 tsObjId,
                                                 prerequisites=[cisID])

                dmsID = self._insertFunctionStep(self.detectMisalignmentStep,
                                                 tsObjId,
                                                 prerequisites=[crvID])

                gosID = self._insertFunctionStep(self.generateOutputStep,
                                                 tsObjId,
                                                 prerequisites=[dmsID])

                allcossId.append(gosID)

            self._insertFunctionStep(self.closeOutputSetsStep,
                                     prerequisites=allcossId)

        else:  # SetOfLandmarkModels
            self.inputSetOfLandmarkModels = self.inputSet.get()
            self.inputSetOfTiltSeries = self.inputSet.get().getSetOfTiltSeries(pointer=True)

            for ts in self.inputSetOfTiltSeries:
                tsObjId = ts.getObjId()
                grfID = self._insertFunctionStep(self.generateResidualFileFromLandmarkModel,
                                                 tsObjId,
                                                 prerequisites=[])

                dmsID = self._insertFunctionStep(self.detectMisalignmentStep,
                                                 tsObjId,
                                                 prerequisites=[grfID])

                gosID = self._insertFunctionStep(self.generateOutputStep,
                                                 tsObjId,
                                                 prerequisites=[dmsID])

                allcossId.append(gosID)

            self._insertFunctionStep(self.closeOutputSetsStep,
                                     prerequisites=allcossId)

    # --------------------------- STEPS functions ----------------------------

    def convertInputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        firstItem = ts.getFirstItem()

        """Apply the transformation form the input tilt-series"""
        # Use Xmipp interpolation via Scipion
        swap = False
        if firstItem.hasTransform():
            avgRotAngle = utils.calculateAverageRotationAngleFromTM(ts)
            swap = True if (avgRotAngle > 45 or avgRotAngle < -45) else False

            outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
            ts.applyTransform(outputTsFileName, swapXY=swap)

        else:
            outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
            ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))
        utils.writeXmippMetadataTiltAngleList(ts, angleFilePath)

        """Generate 3D coordinates metadata"""
        sots_soc = self.inputSetOfCoordinates.get().getSetOfTiltSeries()
        ts_soc = sots_soc.getTiltSeriesFromTsId(tsId)

        xDim, yDim, _ = ts_soc.getFirstItem().getDimensions()

        if firstItem.hasTransform():
            if swap:
                xHalf = yDim / 2
                yHalf = xDim / 2
            else:
                xHalf = xDim / 2
                yHalf = yDim / 2
        else:
            xHalf = firstItem.getDimensions()[0] / 2
            yHalf = firstItem.getDimensions()[1] / 2

        self.check = utils.writeOutputTiltSeriesCoordinates3dXmdFile(self.inputSetOfCoordinates.get(),
                                                                     os.path.join(extraPrefix,
                                                                                  METADATA_INPUT_COORDINATES),
                                                                     ts.getSamplingRate(),
                                                                     xHalf,
                                                                     yHalf,
                                                                     tsId)

        if not self.check:
            print("No input coordinates for ts %s. Skipping this tilt-series for analysis." % tsId)

    def calculateResidualVectors(self, tsObjId):
        if self.check:
            ts = self.inputSetOfTiltSeries[tsObjId]
            tsId = ts.getTsId()

            extraPrefix = self._getExtraPath(tsId)
            tmpPrefix = self._getTmpPath(tsId)

            firstItem = ts.getFirstItem()

            angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))

            paramsLandmarkResiduals = {
                'i': os.path.join(tmpPrefix, firstItem.parseFileName() + ":mrcs"),
                'tlt': angleFilePath,
                'inputCoord': os.path.join(extraPrefix, METADATA_INPUT_COORDINATES),
                'o': os.path.join(extraPrefix, VRESMOD_FILE_NAME),
                'samplingRate': self.inputSetOfTiltSeries.getSamplingRate(),
                'fiducialSize': self.fiducialSize.get() * 10,
                'thrSDHCC': self.thrSDHCC.get(),
                'thrFiducialDistance': self.thrFiducialDistance.get(),
                'targetLMsize': self.targetLMsize.get()
            }

            argsLandmarkResiduals = "-i %(i)s " \
                                    "--tlt %(tlt)s " \
                                    "--inputCoord %(inputCoord)s " \
                                    "-o %(o)s " \
                                    "--samplingRate %(samplingRate).2f " \
                                    "--fiducialSize %(fiducialSize).2f " \
                                    "--thrSDHCC %(thrSDHCC).2f " \
                                    "--thrFiducialDistance %(thrFiducialDistance).2f " \
                                    "--targetLMsize %(targetLMsize).2f"

            self.runJob('xmipp_tomo_calculate_landmark_residuals', argsLandmarkResiduals % paramsLandmarkResiduals)

    def detectMisalignmentStep(self, tsObjId):
        if self.check:
            ts = self.inputSetOfTiltSeries[tsObjId]
            tsId = ts.getTsId()

            extraPrefix = self._getExtraPath(tsId)
            tmpPrefix = self._getTmpPath(tsId)

            firstItem = ts.getFirstItem()

            paramsDetectMisali = {
                'i': os.path.join(tmpPrefix, firstItem.parseFileName() + ":mrcs"),
                'inputResInfo': os.path.join(extraPrefix, VRESMOD_FILE_NAME),
                'o': os.path.join(extraPrefix, firstItem.parseFileName(suffix='_alignmentReport', extension='.xmd')),
                'samplingRate': self.inputSetOfTiltSeries.getSamplingRate(),
                'fiducialSize': self.fiducialSize.get() * 10,
                'thrFiducialDistance': self.thrFiducialDistance.get(),
            }

            argsDetectMisali = "-i %(i)s " \
                               "--inputResInfo %(inputResInfo)s " \
                               "-o %(o)s " \
                               "--samplingRate %(samplingRate).2f " \
                               "--fiducialSize %(fiducialSize).2f " \
                               "--thrFiducialDistance %(thrFiducialDistance).2f "

            self.runJob('xmipp_tomo_detect_misalignment_residuals', argsDetectMisali % paramsDetectMisali)

    def generateOutputStep(self, tsObjId):
        if self.check:
            ts = self.inputSetOfTiltSeries[tsObjId]
            tsId = ts.getTsId()

            extraPrefix = self._getExtraPath(tsId)

            firstItem = ts.getFirstItem()

            xmdEnableTiltImages = os.path.join(extraPrefix,
                                               firstItem.parseFileName(suffix='_alignmentReport', extension='.xmd'))

            enableInfoList = utils.readXmippMetadataEnabledTiltImages(xmdEnableTiltImages)

            self.generateAlignmentReportDictionary(enableInfoList, tsId)

            # Generate output sets of aligned and misaligned tilt series
            aligned = True

            # Check if some tilt image presents misalignment
            for line in enableInfoList:
                if float(line[0]) != 1:
                    aligned = False
                    break

            newTs = tomoObj.TiltSeries(tsId=tsId)
            newTs.copyInfo(ts)

            # Generate output landmark model
            lmFileName = os.path.join(extraPrefix, firstItem.parseFileName(suffix='_lm', extension='.sfid'))
            lm = tomoObj.LandmarkModel(tsId=tsId,
                                       fileName=lmFileName,
                                       modelName=None,
                                       size=self.fiducialSize.get() * 10,
                                       applyTSTransformation=Boolean(True))
            lm.setTiltSeries(newTs)

            vcmInfoList = self.parseVCMFile(vcmFilePath=os.path.join(extraPrefix, VRESMOD_FILE_NAME))

            for lmInfo in vcmInfoList:
                lm.addLandmark(xCoor=lmInfo[0],
                               yCoor=lmInfo[1],
                               tiltIm=lmInfo[2],
                               chainId=lmInfo[3],
                               xResid=lmInfo[4],
                               yResid=lmInfo[5])

            if aligned:
                self.getOutputSetOfTiltSeries()
                self.outputSetOfTiltSeries.append(newTs)

                self.getOutputSetOfAlignedLandmarkModels()
                self.outputSetOfAlignedLandmarkModels.append(lm)
            else:
                self.getOutputSetOfMisalignedTiltSeries()
                self.outputSetOfMisalignedTiltSeries.append(newTs)

                self.getOutputSetOfMisalignedLandmarkModels()
                self.outputSetOfMisalignedLandmarkModels.append(lm)

            for index, tiltImage in enumerate(ts):
                newTi = tomoObj.TiltImage()
                newTi.copyInfo(tiltImage, copyId=True)
                newTi.setAcquisition(tiltImage.getAcquisition())
                newTi.setLocation(tiltImage.getLocation())
                newTs.append(newTi)

            newTs.write(properties=False)

            if aligned:
                self.outputSetOfTiltSeries.update(newTs)
                self.outputSetOfTiltSeries.write()
            else:
                self.outputSetOfMisalignedTiltSeries.update(newTs)
                self.outputSetOfMisalignedTiltSeries.write()

            self._store()

    def closeOutputSetsStep(self):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)
            self.outputSetOfTiltSeries.write()

        if hasattr(self, "outputSetOfMisalignedTiltSeries"):
            self.outputSetOfMisalignedTiltSeries.setStreamState(Set.STREAM_CLOSED)
            self.outputSetOfMisalignedTiltSeries.write()

        if hasattr(self, "outputSetOfAlignedLandmarkModels"):
            self.outputSetOfAlignedLandmarkModels.setStreamState(Set.STREAM_CLOSED)
            self.outputSetOfAlignedLandmarkModels.write()

        if hasattr(self, "outputSetOfMisalignedTiltSeries"):
            self.outputSetOfMisalignedLandmarkModels.setStreamState(Set.STREAM_CLOSED)
            self.outputSetOfMisalignedLandmarkModels.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def generateAlignmentReportDictionary(self, enableInfoList, tsId):
        self.alignmentReport.clear()
        self.alignmentReport = List([])

        globalMisalignment = True
        aligned = True

        misalignedTiltImages = String("%s misalignment detected in images: " % tsId)

        for line in enableInfoList:
            if float(line[0]) != 1:
                aligned = False

                previousLine = misalignedTiltImages.get()
                line = " " + str(line[1]) + ","

                misalignedTiltImages = String(previousLine + line)
            else:
                globalMisalignment = False

        if globalMisalignment:
            self.alignmentReport.append(String("%s global misalignment detected" % tsId))
        elif not aligned:
            # Remove final comma
            previousLine = misalignedTiltImages.get()
            previousLine = previousLine[:-1]
            self.alignmentReport.append(String(previousLine))

        # Un comment if it makes sense of report residual statistics. This requires fixes in code.
        # for line in self.convertResidualStatisticInString(tsId):
        #     self.alignmentReport.append(String(line))

        self._store()

    @staticmethod
    def parseVCMFile(vcmFilePath):
        """ Method to retrieve the information contained in the vcm (vector of coordinate models) file in a list of
        lists"""

        vCMinfo = []
        table = emtable.Table(fileName=vcmFilePath)

        for row in table.iterRows(fileName='noname@' + vcmFilePath):
            chainId = row.get('frameId')
            xResid = row.get('shiftX')
            yResid = row.get('shiftY')
            xCoor = row.get('x')
            yCoor = row.get('y')
            tiltIm = row.get('z')

            vCMinfo.append([xCoor, yCoor, tiltIm, chainId, xResid, yResid])

        return vCMinfo

    def getOutputSetOfTiltSeries(self):
        """ Method to generate output classes of set of tilt-series"""

        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self.outputSetOfTiltSeries.enableAppend()

        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries(suffix="_ali")

            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries)
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.getDim())

            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSet, outputSetOfTiltSeries)

        return self.outputSetOfTiltSeries

    def getOutputSetOfMisalignedTiltSeries(self):
        """ Method to generate output classes of set of tilt-series"""

        if hasattr(self, "outputSetOfMisalignedTiltSeries"):
            self.outputSetOfMisalignedTiltSeries.setStreamState(Set.STREAM_OPEN)
            self.outputSetOfMisalignedTiltSeries.enableAppend()

        else:
            outputSetOfMisalignedTiltSeries = self._createSetOfTiltSeries(suffix="_misali")

            outputSetOfMisalignedTiltSeries.copyInfo(self.inputSetOfTiltSeries)
            outputSetOfMisalignedTiltSeries.setDim(self.inputSetOfTiltSeries.getDim())

            outputSetOfMisalignedTiltSeries.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfMisalignedTiltSeries=outputSetOfMisalignedTiltSeries)
            self._defineSourceRelation(self.inputSet, outputSetOfMisalignedTiltSeries)

        return self.outputSetOfMisalignedTiltSeries

    def getOutputSetOfAlignedLandmarkModels(self):
        if hasattr(self, "outputSetOfAlignedLandmarkModels"):
            self.outputSetOfAlignedLandmarkModels.setStreamState(Set.STREAM_OPEN)
            self.outputSetOfAlignedLandmarkModels.enableAppend()

        else:
            outputSetOfAlignedLandmarkModels = self._createSetOfLandmarkModels("_ali")

            outputSetOfAlignedLandmarkModels.copyInfo(self.inputSetOfTiltSeries)
            outputSetOfAlignedLandmarkModels.setSetOfTiltSeries(self.getOutputSetOfTiltSeries())

            outputSetOfAlignedLandmarkModels.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfAlignedLandmarkModels=outputSetOfAlignedLandmarkModels)
            self._defineSourceRelation(self.inputSet, outputSetOfAlignedLandmarkModels)

        return self.outputSetOfAlignedLandmarkModels

    def getOutputSetOfMisalignedLandmarkModels(self):
        if hasattr(self, "outputSetOfMisalignedLandmarkModels"):
            self.outputSetOfMisalignedLandmarkModels.setStreamState(Set.STREAM_OPEN)
            self.outputSetOfMisalignedLandmarkModels.enableAppend()

        else:
            outputSetOfMisalignedLandmarkModels = self._createSetOfLandmarkModels("_misali")
            outputSetOfMisalignedLandmarkModels.setSetOfTiltSeries(self.getOutputSetOfMisalignedTiltSeries())

            outputSetOfMisalignedLandmarkModels.copyInfo(self.inputSetOfTiltSeries)

            outputSetOfMisalignedLandmarkModels.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfMisalignedLandmarkModels=outputSetOfMisalignedLandmarkModels)
            self._defineSourceRelation(self.inputSet, outputSetOfMisalignedLandmarkModels)

        return self.outputSetOfMisalignedLandmarkModels

    def convertResidualStatisticInString(self, tsId):
        extraPrefix = self._getExtraPath(tsId)
        statisticsInfoTable = utils.readResidualStatisticsXmdFile(os.path.join(extraPrefix, "residualStatistics.xmd"))

        outputStrings = []

        for key in statisticsInfoTable.keys():
            outputStrings.append("Info form coordinate " + str(str(statisticsInfoTable[key][4]) +
                                                               ":  convex hull area: " + str(
                statisticsInfoTable[key][0]) +
                                                               ", convex hull perimeter: " + str(
                statisticsInfoTable[key][1]) +
                                                               ", passed tests: " + str(statisticsInfoTable[key][2]) +
                                                               ", failed test: " + str(statisticsInfoTable[key][3])))

        return outputStrings

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = ["MISALIGNMENT REPORT"]

        if not hasattr(self, 'outputSetOfMisalignedTiltSeries'):
            summary.append("No tilt series present misalignment")

            # for ts in self.outputSetOfTiltSeries:
            #     summary.extend(self.convertResidualStatisticInString(ts.getTsId()))
        else:
            summary.append("From the %d analyzed tilt series, %d presents misalignment:" %
                           (self.inputSet.get().getSize(),
                            self.outputSetOfMisalignedTiltSeries.getSize()))

            # for ts in self.outputSetOfMisalignedTiltSeries:
            #     summary.extend(self.convertResidualStatisticInString(ts.getTsId()))

            for line in self.alignmentReport:
                summary.append(line.get())

        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfTiltSeries'):
            methods.append("%d tilt-series have been classified as properly aligned using the Xmipp program "
                           "xmipp_tomo_detect_misalignment.\n"
                           % (self.outputSetOfTiltSeries.getSize()))
        elif hasattr(self, 'outputSetOfMisalignedTiltSeries'):
            methods.append("%d tilt-series have been classified as misaligned using the Xmipp program "
                           "xmipp_tomo_detect_misalignment.\n"
                           % (self.outputSetOfMisalignedTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
