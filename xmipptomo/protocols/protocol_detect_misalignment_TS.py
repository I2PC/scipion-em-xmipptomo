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
from pyworkflow.object import Set, List, String
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pwem.protocols import EMProtocol
import tomo.objects as tomoObj
from tomo.protocols import ProtTomoBase
import xmipptomo.utils as utils

METADATA_INPUT_COORDINATES = "fiducialCoordinates.xmd"


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

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series')

        form.addParam('inputSetOfCoordinates',
                      params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      important=True,
                      label='Input set of coordinates 3D',
                      help='Set of 3D coordinates indicating the position in space of the fiducials. This set should '
                           'be obtained from the previous alignment step of the tilt-series.')

        form.addParam('fiducialSize',
                      params.FloatParam,
                      important=True,
                      label='Fiducial size (nm)',
                      help='Fiducial size in nanometers (nm).')

        form.addParam('thrSDHCC',
                      params.FloatParam,
                      advanced=True,
                      default=5,
                      label='Coordinate value SD threshold',
                      help='Number of SD a coordinate value must be over the mean to consider that it belongs to a '
                           'high contrast feature.')

        form.addParam('thrNumberCoords',
                      params.IntParam,
                      advanced=True,
                      default=10,
                      label='Number of coordinates threshold',
                      help='Minimum number of coordinates attracted to a center of mass to consider it as a high '
                           'contrast feature.')

        form.addParam('thrChainDistanceAng',
                      params.FloatParam,
                      advanced=True,
                      default=10,
                      label='Landmark distance threshold',
                      help='Threshold maximum distance in angstroms of a detected landmark to consider it belongs to '
                           'a chain.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self.alignmentReport = List([])

        allcossId = []

        for ts in self.inputSetOfTiltSeries.get():
            tsObjId = ts.getObjId()
            cisID = self._insertFunctionStep(self.convertInputStep,
                                             tsObjId,
                                             prerequisites=[])
            dmsID = self._insertFunctionStep(self.detectMisalignmentStep,
                                             tsObjId,
                                             prerequisites=[cisID])
            gosID = self._insertFunctionStep(self.generateOutputStep,
                                             tsObjId,
                                             prerequisites=[dmsID])

            allcossId.append(gosID)

        self._insertFunctionStep(self.closeOutputSetsStep,
                                 prerequisites=allcossId)

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
        outputTsFileName = os.path.join(tmpPrefix, firstItem.parseFileName())
        ts.applyTransform(outputTsFileName)

        """Generate angle file"""
        angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))
        utils.writeXmippMetadataTiltAngleList(ts, angleFilePath)

        """Generate 3D coordinates metadata"""
        utils.writeOutputCoordinates3dXmdFile(self.inputSetOfCoordinates.get(),
                                              os.path.join(extraPrefix, METADATA_INPUT_COORDINATES),
                                              tsObjId)

    def detectMisalignmentStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        firstItem = ts.getFirstItem()

        angleFilePath = os.path.join(tmpPrefix, firstItem.parseFileName(extension=".tlt"))

        paramsDetectMisali = {
            'i': os.path.join(tmpPrefix, firstItem.parseFileName() + ":mrcs"),
            'inputCoord': os.path.join(extraPrefix, METADATA_INPUT_COORDINATES),
            'tlt': angleFilePath,
            'o': os.path.join(extraPrefix, firstItem.parseFileName(suffix='_alignmentReport', extension='.xmd')),
            'thrSDHCC': self.thrSDHCC.get(),
            'thrNumberCoords': self.thrNumberCoords.get(),
            'samplingRate': self.inputSetOfTiltSeries.get().getSamplingRate(),
            'fiducialSize': self.fiducialSize.get() * 10,
            'thrChainDistanceAng': self.thrChainDistanceAng.get(),
        }

        argsDetectMisali = "-i %(i)s " \
                           "--tlt %(tlt)s " \
                           "--inputCoord %(inputCoord)s " \
                           "-o %(o)s " \
                           "--thrSDHCC %(thrSDHCC).2f " \
                           "--thrNumberCoords %(thrNumberCoords).2f " \
                           "--samplingRate %(samplingRate).2f " \
                           "--fiducialSize %(fiducialSize).2f " \
                           "--thrChainDistanceAng %(thrChainDistanceAng).2f"

        self.runJob('xmipp_tomo_detect_misalignment_trajectory', argsDetectMisali % paramsDetectMisali)

    def generateOutputStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        firstItem = ts.getFirstItem()

        xmdEnableTiltImages = os.path.join(extraPrefix,
                                           firstItem.parseFileName(suffix='_alignmentReport', extension='.xmd'))

        enableInfoList = utils.readXmippMetadataEnabledTiltImages(xmdEnableTiltImages)

        self.generateAlignmentReportDictionary(enableInfoList, tsId)

        """ Generate output sets of aligned and misaligned tilt series """
        aligned = True

        # Check if some tilt image presents misalignment
        for line in enableInfoList:
            if float(line[0]) != 1:
                aligned = False
                break

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)

        if aligned:
            self.getOutputSetOfTiltSeries()
            self.outputSetOfTiltSeries.append(newTs)
        else:
            self.getOutputSetOfMisalignedTiltSeries()
            self.outputSetOfMisalignedTiltSeries.append(newTs)

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

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def generateAlignmentReportDictionary(self, enableInfoList, tsId):
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

        for line in self.convertResidualStatisticInString(tsId):
            self.alignmentReport.append(String(line))

        self._store()

    def getOutputSetOfTiltSeries(self):
        """ Method to generate output classes of set of tilt-series"""

        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self.outputSetOfTiltSeries.enableAppend()

        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries(suffix="_ali")

            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())

            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)

        return self.outputSetOfTiltSeries

    def getOutputSetOfMisalignedTiltSeries(self):
        """ Method to generate output classes of set of tilt-series"""

        if hasattr(self, "outputSetOfMisalignedTiltSeries"):
            self.outputSetOfMisalignedTiltSeries.setStreamState(Set.STREAM_OPEN)
            self.outputSetOfMisalignedTiltSeries.enableAppend()

        else:
            outputSetOfMisalignedTiltSeries = self._createSetOfTiltSeries(suffix="_misali")

            outputSetOfMisalignedTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfMisalignedTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())

            outputSetOfMisalignedTiltSeries.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfMisalignedTiltSeries=outputSetOfMisalignedTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfMisalignedTiltSeries)

        return self.outputSetOfMisalignedTiltSeries

    def getOutputSetOfLandmarkModels(self):
        if hasattr(self, "outputSetOfLandmarkModels"):
            self.outputSetOfLandmarkModels.enableAppend()

        else:
            outputSetOfLandmarkModels = self._createSetOfLandmarkModels()

            outputSetOfLandmarkModels.copyInfo(self.inputSetOfTiltSeries.get())

            outputSetOfLandmarkModels.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfLandmarkModels=outputSetOfLandmarkModels)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfLandmarkModels)

        return self.outputSetOfLandmarkModels

    def convertResidualStatisticInString(self, tsId):
        extraPrefix = self._getExtraPath(tsId)
        statisticsInfoTable = utils.readResidualStatisticsXmdFile(os.path.join(extraPrefix, "residualStatistics.xmd"))

        outputStrings = []

        for key in statisticsInfoTable.keys():
            outputStrings.append("Info form coordinate " + str(str(statisticsInfoTable[key][4]) +
                                 ":  convex hull area: " + str(statisticsInfoTable[key][0]) +
                                 ", convex hull perimeter: " + str(statisticsInfoTable[key][1]) +
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
                           (self.inputSetOfTiltSeries.get().getSize(),
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
