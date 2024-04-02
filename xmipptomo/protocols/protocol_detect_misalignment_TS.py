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

from pwem.emlib import lib
import pwem.emlib.metadata as md
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
VRESMOD_FILE_NAME_EXT = "_vResMod.xmd"


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

        # Global varibale to keep the quality of the tilt-series under study.
        # If the tilt series is not aligned the detection of subtle misalignment is avoided
        self.aligned = True

        self.inputSetOfTiltSeries = None
        self.inputSetOfLandmarkModels = None

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSet',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries, SetOfLandmarkModels',
                      important=True,
                      label='Input set of tilt-series',
                      help='Input set of tilt-series or landmark models:\n'
                           '- Tilt-series: landmarks are detected in the tilt-series to posteriorly calculate the '
                           'residual vectors and finally analyze them for misalignment detection.\n'
                           '- Landmark models: calculated residuals by the alignment algorithm are directly analyze to '
                           'detect misalignment.')

        form.addParam('inputSetOfCoordinates',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeriesCoordinates',
                      condition="inputSet is not None and inputSet.getClassName()=='SetOfTiltSeries'",
                      important=True,
                      label='Input set of coordinates 3D',
                      help='Set of 3D coordinates indicating the position in space of the fiducials. This set should '
                           'be obtained from the previous alignment step of the tilt-series.')

        form.addParam('fiducialSize',
                      params.FloatParam,
                      important=True,
                      label='Fiducial size (nm)',
                      help='Fiducial size in nanometers (nm).')

        form.addParam('maxMisaliImages',
                      params.IntParam,
                      default=3,
                      label='Max. misaligned images',
                      help='Maximum number of tilt-images that might present misalignment to keep series as aligned. '
                           'Default value is 3, meaning that if 3 or less tilt-images present misalignment they are '
                           'annotated but the tilt-series is not classified as misaligned as a whole.')

        form.addParam('subtleMisaliToggle',
                      params.BooleanParam,
                      default=True,
                      label='Subtle misalignment analysis',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Run local alignment to detect subtle misalignment. This analysis detects lower alignment '
                           'errors but it requires some computational extra load.')

        form.addParam('subtleMisalignmentTolerance',
                      params.IntParam,
                      default=3,
                      condition='subtleMisaliToggle',
                      label='Misalignment tolerance (px)',
                      help='Maximum displacement between to consecutive images to consider that misalignment is '
                           'present. This displacement is calculated by correlation.')

        form.addParam('inputTsFromLm',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      condition="inputSet is not None and inputSet.getClassName()=='SetOfLandmarkModels' and "
                                "subtleMisaliToggle",
                      allowsNull=True,
                      label='Tilt-series associated to the LM',
                      help='Input set of tilt-series associated to the provided landmark model. This tilt-series is '
                           'mean to be the tilt-series aligned by the algorithm that has produces the set of landmark'
                           'models. It might be the interpolated tilt-series or a non-interpolated one with its '
                           'associated alignment.')

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

        form.addParam('targetLMsize',
                      params.FloatParam,
                      default=8,
                      label='Target fiducial size (px)',
                      help='Target fiducial size in pixels to calculate the downsampling when detecting landmarks'
                           'on the tilt series.',
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParallelSection(threads=4, mpi=1)

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
                                                 True,
                                                 prerequisites=[])

                crvID = self._insertFunctionStep(self.calculateResidualVectors,
                                                 tsObjId,
                                                 prerequisites=[cisID])

                dmsID = self._insertFunctionStep(self.detectMisalignmentStep,
                                                 tsObjId,
                                                 prerequisites=[crvID])

                if self.subtleMisaliToggle.get():
                    dsmsID = self._insertFunctionStep(self.detectSubtleMisalignment,
                                                      tsObjId,
                                                      prerequisites=[dmsID])

                    gosID = self._insertFunctionStep(self.generateOutputStep,
                                                     tsObjId,
                                                     prerequisites=[dsmsID])

                    allcossId.append(gosID)
                else:
                    gosID = self._insertFunctionStep(self.generateOutputStep,
                                                     tsObjId,
                                                     prerequisites=[dmsID])

                    allcossId.append(gosID)

            self._insertFunctionStep(self.closeOutputSetsStep,
                                     prerequisites=allcossId)

        else:  # SetOfLandmarkModels
            self.inputSetOfLandmarkModels = self.inputSet.get()
            if self.subtleMisaliToggle.get():
                self.inputSetOfTiltSeries = self.inputTsFromLm.get()
            else:
                self.inputSetOfTiltSeries = self.inputSet.get().getSetOfTiltSeries()

            for lm in self.inputSetOfLandmarkModels:
                lmTsId = lm.getTsId()
                ts = self.inputSetOfTiltSeries.getTiltSeriesFromTsId(lmTsId)
                tsObjId = ts.getObjId()

                cisID = self._insertFunctionStep(self.convertInputStep,
                                                 tsObjId,
                                                 False,
                                                 prerequisites=[])

                grfID = self._insertFunctionStep(self.generateResidualFileFromLandmarkModel,
                                                 tsObjId,
                                                 prerequisites=[cisID])

                dmsID = self._insertFunctionStep(self.detectMisalignmentStep,
                                                 tsObjId,
                                                 prerequisites=[grfID])

                if self.subtleMisaliToggle.get():
                    dsmsID = self._insertFunctionStep(self.detectSubtleMisalignment,
                                                      tsObjId,
                                                      prerequisites=[dmsID])

                    gosID = self._insertFunctionStep(self.generateOutputStep,
                                                     tsObjId,
                                                     prerequisites=[dsmsID])

                    allcossId.append(gosID)
                else:
                    gosID = self._insertFunctionStep(self.generateOutputStep,
                                                     tsObjId,
                                                     prerequisites=[dmsID])

                    allcossId.append(gosID)

                allcossId.append(gosID)

            self._insertFunctionStep(self.closeOutputSetsStep,
                                     prerequisites=allcossId)

    # --------------------------- STEPS functions ----------------------------
    def convertInputStep(self, tsObjId, modeTs):
        ts = self.inputSetOfTiltSeries[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        firstItem = ts.getFirstItem()

        if modeTs or self.subtleMisaliToggle.get():
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

        if modeTs:
            """Generate 3D coordinates metadata"""
            xDim, yDim, _ = ts.getFirstItem().getDimensions()

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

    def generateResidualFileFromLandmarkModel(self, tsObjId):
        ts = self.inputSetOfTiltSeries[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        lm = self.inputSetOfLandmarkModels.getLandmarkModelFromTsId(tsId)

        lmInfoTable = lm.retrieveInfoTable()
        resModFilePath = os.path.join(extraPrefix, (tsId + VRESMOD_FILE_NAME_EXT))

        mdlm = lib.MetaData()

        # Check that there is residual information in the input landmark models
        self.check = False

        for infoLine in lmInfoTable:
            nRow = md.Row()

            self.check = True

            # Only consider those landmark models with residual information
            if infoLine[3] == 'nan' or infoLine[4] == 'nan':
                continue

            nRow.setValue(lib.MDL_X, float(infoLine[0]))
            nRow.setValue(lib.MDL_Y, float(infoLine[1]))
            nRow.setValue(lib.MDL_Z, float(infoLine[2]) - 1)
            nRow.setValue(lib.MDL_FRAME_ID, int(infoLine[3]) - 1)
            nRow.setValue(lib.MDL_SHIFT_X, float(infoLine[4]))
            nRow.setValue(lib.MDL_SHIFT_Y, float(infoLine[5]))

            nRow.addToMd(mdlm)

        mdlm.write(resModFilePath)

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
                'o': os.path.join(extraPrefix, (tsId + VRESMOD_FILE_NAME_EXT)),
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
                'inputResInfo': os.path.join(extraPrefix, (tsId + VRESMOD_FILE_NAME_EXT)),
                'o': os.path.join(extraPrefix, firstItem.parseFileName(suffix='_alignmentReport', extension='.xmd')),
                'samplingRate': self.inputSetOfTiltSeries.getSamplingRate(),
                'fiducialSize': self.fiducialSize.get() * 10,
                'thrFiducialDistance': self.thrFiducialDistance.get(),
            }

            # "-i %(i)s " \
            argsDetectMisali = "--inputResInfo %(inputResInfo)s " \
                               "-o %(o)s " \
                               "--samplingRate %(samplingRate).2f " \
                               "--fiducialSize %(fiducialSize).2f " \
                               "--thrFiducialDistance %(thrFiducialDistance).2f "

            self.runJob('xmipp_tomo_detect_misalignment_residuals', argsDetectMisali % paramsDetectMisali)

            # Detect if tilt-series presents misalignment. If not subtle misalignment detection will be executed
            xmdEnableTiltImages = os.path.join(extraPrefix,
                                               firstItem.parseFileName(suffix='_alignmentReport', extension='.xmd'))

            enableInfoList = utils.readXmippMetadataEnabledTiltImages(xmdEnableTiltImages)

            self.generateAlignmentReportDictionary(enableInfoList, tsId)

            # Check number of locally misaligned tilt images
            misaliTi = 0

            for line in enableInfoList:
                if float(line[0]) != 1:
                    misaliTi += 1

            self.aligned = True if misaliTi <= self.maxMisaliImages.get() else False

    def detectSubtleMisalignment(self, tsObjId):
        if self.check and self.aligned:
            ts = self.inputSetOfTiltSeries[tsObjId]
            tsId = ts.getTsId()

            extraPrefix = self._getExtraPath(tsId)
            tmpPrefix = self._getTmpPath(tsId)

            firstItem = ts.getFirstItem()

            angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))

            paramsDetectMisali = {
                'i': os.path.join(tmpPrefix, firstItem.parseFileName() + ":mrcs"),
                'tlt': angleFilePath,
                'o': os.path.join(extraPrefix,
                                  firstItem.parseFileName(suffix='_subtleMisalignmentReport', extension='.xmd')),
                'samplingRate': self.inputSetOfTiltSeries.getSamplingRate(),
                'shiftTol': self.subtleMisalignmentTolerance.get(),
            }

            argsDetectMisali = "-i %(i)s " \
                               "--tlt %(tlt)s " \
                               "-o %(o)s " \
                               "--samplingRate %(samplingRate).2f " \
                               "--shiftTol %(shiftTol).2f "

            self.runJob('xmipp_tomo_tiltseries_detect_misalignment_corr', argsDetectMisali % paramsDetectMisali)

    def generateOutputStep(self, tsObjId):
        if self.check:
            ts = self.inputSetOfTiltSeries[tsObjId]
            tsId = ts.getTsId()

            extraPrefix = self._getExtraPath(tsId)

            firstItem = ts.getFirstItem()

            # If tilt-series is still aligned search for subtle misalignment
            if self.aligned and self.subtleMisaliToggle.get():
                xmdEnableTiltImages = os.path.join(
                    extraPrefix,
                    firstItem.parseFileName(suffix='_subtleMisalignmentReport', extension='.xmd'))

                enableInfoList = utils.readXmippMetadataEnabledTiltImages(xmdEnableTiltImages)

                # Check number of locally misaligned tilt images
                misaliTi = 0

                for line in enableInfoList:
                    if float(line[0]) != 1:
                        misaliTi += 1

                self.aligned = True if misaliTi <= self.maxMisaliImages.get() else False

            # Generate output sets of aligned and misaligned tilt series
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

            vcmInfoList = self.parseVCMFile(vcmFilePath=os.path.join(extraPrefix, (tsId + VRESMOD_FILE_NAME_EXT)))

            for lmInfo in vcmInfoList:
                lm.addLandmark(xCoor=lmInfo[0],
                               yCoor=lmInfo[1],
                               tiltIm=lmInfo[2],
                               chainId=lmInfo[3],
                               xResid=lmInfo[4],
                               yResid=lmInfo[5])

            if self.aligned:
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

            if self.aligned:
                self.outputSetOfTiltSeries.update(newTs)
                self.outputSetOfTiltSeries.write()
            else:
                self.outputSetOfMisalignedTiltSeries.update(newTs)
                self.outputSetOfMisalignedTiltSeries.write()

            self._store()

        # Return check to true for the next iteration
        self.check = True

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
    @staticmethod
    def is_float_zero(value, tolerance=1e-9):
        return abs(value) < tolerance

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
    def _validate(self):
        errors = []

        if isinstance(self.inputSet.get(), tomoObj.SetOfLandmarkModels):
            if not self.inputSet.get().hasResidualInfo():
                errors.append("Input set of landmark models has no residual information. "
                              "Impossible to study misalignment.")
                errors.append("Typically only landmark models obtained after alignment has residual information.")

        return errors

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
