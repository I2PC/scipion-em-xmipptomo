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

from pwem.emlib import MetaData, MDL_MAX, MDL_MIN
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Set
from pyworkflow.protocol import PointerParam, EnumParam, FloatParam, BooleanParam, LEVEL_ADVANCED
from tomo.objects import SetOfCoordinates3D, SetOfTomograms, Coordinate3D, SubTomogram, SetOfSubTomograms
from tomo.protocols import ProtTomoBase
from tomo import constants
from xmipp3 import XmippProtocol
from xmipptomo import utils

COORDINATES_FILE_NAME = 'subtomo_coords.xmd'
TARGET_BOX_SIZE = 32


class XmippProtDeepDetectMisalignment(EMProtocol, ProtTomoBase, XmippProtocol):
    """
    Wrapper protocol to Xmipp xmipp_deep_misalignment_detection for misalignment detection
    in tomographic reconstructions based on artifacted landmarks
    """

    _label = 'detect misalignment from fiducials'
    _devStatus = BETA
    _conda_env = 'xmipp_DLTK_v1.0'
    _possibleOutputs = {"strongMisalignedTomograms": SetOfTomograms,
                        "weakMisalignedTomograms": SetOfTomograms,
                        "alignedTomograms": SetOfTomograms,
                        "outputSubtomos": SetOfSubTomograms}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.alignedTomograms = None
        self.strongMisalignedTomograms = None
        self.weakMisalignedTomograms = None
        self.outputSubtomos = None
        self.isot = None

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSetOfCoordinates',
                      PointerParam,
                      pointerClass=SetOfCoordinates3D,
                      label='Fiducial 3D coordinates',
                      help='3D coordinates indicating the location of the fiducials (gold beads) in the tomogram. '
                           'These fiducails will be the ones used to study misalignment artifacts over them. '
                           'The coordinate denotes the center of the subtomogram')

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

        form.addParam('inputSetOfTomograms',
                      PointerParam,
                      pointerClass=SetOfTomograms,
                      condition='tomoSource==1',
                      allowsNull=True,
                      label='Input tomograms',
                      help='Tomograms from which extract the fiducials (gold beads) at the specified coordinates '
                           'locations.')

        form.addParam('fiducialSize',
                      FloatParam,
                      default=10,
                      label='Fiducial size (nm)',
                      help='Peaked gold bead size in nanometers.')

        form.addParam('misaliThrBool',
                      BooleanParam,
                      default=True,
                      label='Use misalignment threshold?',
                      help='Threshold to settle if a tomogram presents weak or strong misalignment. If this value is '
                           'not provided two output set of tomograms are generated, those discarded which present '
                           'strong misalignment and those which do not. If this value is provided the second group of '
                           'tomograms is splitted into two, using this threshold to settle if the tomograms present'
                           'or not a weak misalignment.')

        form.addParam('misaliThr',
                      FloatParam,
                      default=0.45,
                      condition='misaliThrBool==True',
                      label='Misalignment threshold',
                      help='Threshold value to settle if a tomogram presents weak or strong misalignment. Value '
                           'ranged between (0, 1).')

        # Advanced parameters
        form.addParam('modelPick',
                      EnumParam,
                      choices=['strict', 'loose'],
                      default=0,
                      display=EnumParam.DISPLAY_HLIST,
                      expertLevel=LEVEL_ADVANCED,
                      label='Model for weak misalignment estimation',
                      help='Choose model for weak misalignment estimation. By default, strict model is picked in '
                           'order to avoid false positives. In case loose model is chosen, less good aligned '
                           'tomograms are lost. As a tradeoff, the number of false positives will increase.')

        form.addParam('misalignmentCriteria',
                      EnumParam,
                      choices=['mean', 'votes'],
                      default=0,
                      display=EnumParam.DISPLAY_HLIST,
                      expertLevel=LEVEL_ADVANCED,
                      label='Misalignment criteria',
                      help='Criteria used for making a decision on the presence of misalignment on the tomogram '
                           'based on the individual scores of each subtomogram. By default the mean of this scores '
                           'is calculated. The other option is to implement a voting system based on if each subtomo '
                           'score is closer to 0 o 1.')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self.tomoDict = self.getTomoDict()

        # Target sampling that fits the fiducial in 16 px (half od the box size to feed the network).
        self.targetSamplingRate = self.fiducialSize.get() / 1.6

        for key in self.tomoDict.keys():
            tomo = self.tomoDict[key]
            coordFilePath = self._getExtraPath(os.path.join(tomo.getTsId()), COORDINATES_FILE_NAME)

            self._insertFunctionStep(utils.writeMdCoordinates,
                                     self.inputSetOfCoordinates.get(),
                                     tomo,
                                     coordFilePath)
            self._insertFunctionStep(self.extractSubtomos,
                                     key,
                                     coordFilePath)
            self._insertFunctionStep(self.subtomoPrediction,
                                     key)
            self._insertFunctionStep(self.createOutputStep,
                                     key,
                                     coordFilePath)

        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEP functions --------------------------------
    def extractSubtomos(self, key, coordFilePath):
        tomo = self.tomoDict[key]

        outputPath = self._getExtraPath(os.path.join(tomo.getTsId()))
        tomoFn = tomo.getFileName()

        dsFactor = self.targetSamplingRate / tomo.getSamplingRate()

        paramsExtractSubtomos = {
            'tomogram': tomoFn,
            'coordinates': coordFilePath,
            'boxsize': TARGET_BOX_SIZE,
            'threads': 1,  # ***
            'outputPath': outputPath,
            'downsample': dsFactor,
        }

        argsExtractSubtomos = "--tomogram %(tomogram)s " \
                              "--coordinates %(coordinates)s " \
                              "--boxsize %(boxsize)d " \
                              "--threads %(threads)d " \
                              "-o %(outputPath)s " \
                              "--downsample %(downsample)f " \
                              "--normalize " \
                              "--fixedBoxSize "

        self.runJob('xmipp_tomo_extract_subtomograms', argsExtractSubtomos % paramsExtractSubtomos)

    def subtomoPrediction(self, key):
        tomo = self.tomoDict[key]

        subtomoFilePath = self._getExtraPath(os.path.join(tomo.getTsId()), COORDINATES_FILE_NAME)

        paramsMisaliPrediction = {
            'modelPick': self.modelPick.get(),
            'subtomoFilePath': subtomoFilePath
        }

        argsMisaliPrediction = "--modelPick %(modelPick)d " \
                               "--subtomoFilePath %(subtomoFilePath)s "

        # Set misalignment threshold
        if self.misaliThrBool.get():
            paramsMisaliPrediction['misaliThr'] = self.misaliThr.get()

            argsMisaliPrediction += "--misaliThr %(misaliThr)f "

        # Set misalignment criteria
        if self.misalignmentCriteria.get() == 1:
            argsMisaliPrediction += "--misalignmentCriteriaVotes "

        self.runJob('xmipp_deep_misalignment_detection',
                    argsMisaliPrediction % paramsMisaliPrediction,
                    env=self.getCondaEnv())

    def createOutputStep(self, key, coordFilePath):
        tomo = self.tomoDict[key]
        tsId = tomo.getTsId()
        subtomoPathList = self.getSubtomoPathList(coordFilePath)

        subtomoFilePath = self._getExtraPath(os.path.join(tsId), COORDINATES_FILE_NAME)
        outputSubtomoXmdFilePath = os.path.join(os.path.dirname(subtomoFilePath), "misalignmentSubtomoStatistics.xmd")
        outputTomoXmdFilePath = os.path.join(os.path.dirname(subtomoFilePath), "misalignmentTomoStatistics.xmd")

        if len(subtomoPathList) != 0:
            self.getOutputSetOfSubtomos()

            subtomoCoordList = utils.retrieveXmipp3dCoordinatesIntoList(coordFilePath, xmdFormat=1)

            firstPredictionArray, secondPredictionArray = self.readPredictionArrays(outputSubtomoXmdFilePath)
            overallPrediction, predictionAverage = self.readTomoScores(outputTomoXmdFilePath)

            print("For volume id " + str(tsId) + " obtained prediction from " + str(len(subtomoPathList)) +
                  " subtomos is " + str(overallPrediction))

            tomo._misaliScore = predictionAverage
            self.addTomoToOutput(tomo=tomo, overallPrediction=overallPrediction)

            for i, subtomoPath in enumerate(subtomoPathList):
                newCoord3D = Coordinate3D()
                newCoord3D.setVolume(tomo)
                newCoord3D.setVolId(i)
                newCoord3D.setX(subtomoCoordList[i][0], constants.BOTTOM_LEFT_CORNER)
                newCoord3D.setY(subtomoCoordList[i][1], constants.BOTTOM_LEFT_CORNER)
                newCoord3D.setZ(subtomoCoordList[i][2], constants.BOTTOM_LEFT_CORNER)

                subtomogram = SubTomogram()
                subtomogram.setLocation(subtomoPath)
                subtomogram.setCoordinate3D(newCoord3D)
                subtomogram.setSamplingRate(self.targetSamplingRate)
                subtomogram.setVolName(tomo.getTsId())
                subtomogram._strongMisaliScore = firstPredictionArray[i]
                subtomogram._weakMisaliScore = secondPredictionArray[i]

                self.outputSubtomos.append(subtomogram)
                self.outputSubtomos.write()
                self._store()

        else:
            print("WARNING: NO SUBTOMOGRAM ESTRACTED FOR TOMOGRAM " + tomo.getTsId() + "IMPOSSIBLE TO STUDY " +
                  "MISALIGNMENT!")

    def closeOutputSetsStep(self):
        if self.alignedTomograms:
            self.alignedTomograms.setStreamState(Set.STREAM_CLOSED)
            self.alignedTomograms.write()

        if self.weakMisalignedTomograms:
            self.weakMisalignedTomograms.setStreamState(Set.STREAM_CLOSED)
            self.weakMisalignedTomograms.write()

        if self.strongMisalignedTomograms:
            self.strongMisalignedTomograms.setStreamState(Set.STREAM_CLOSED)
            self.strongMisalignedTomograms.write()

        self.outputSubtomos.setStreamState(Set.STREAM_CLOSED)
        self.outputSubtomos.write()

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getTomoDict(self):
        if self.tomoSource.get() == 0:
            self.isot = self.inputSetOfCoordinates.get().getPrecedents()
            tomoDict = self.inputSetOfCoordinates.get().getPrecedentsInvolved()
        else:
            tomoDict = {}
            self.isot = self.inputSetOfTomograms.get()

            for tomo in self.inputSetOfTomograms.get():
                tsId = tomo.getTsId()
                tomoDict[tsId] = tomo

        return tomoDict

    def addTomoToOutput(self, tomo, overallPrediction):
        if overallPrediction == 1:  # Strong misali
            self.getOutputSetOfStrongMisalignedTomograms()
            self.strongMisalignedTomograms.append(tomo)
            self.strongMisalignedTomograms.write()
            self._store()

        elif overallPrediction == 2:  # Weak misali
            self.getOutputSetOfWeakMisalignedTomograms()
            self.weakMisalignedTomograms.append(tomo)
            self.weakMisalignedTomograms.write()
            self._store()

        elif overallPrediction == 3:  # Ali
            self.getOutputSetOfAlignedTomograms()
            self.alignedTomograms.append(tomo)
            self.alignedTomograms.write()
            self._store()

    @staticmethod
    def readPredictionArrays(outputSubtomoXmdFilePath):
        mData = MetaData()
        mData.read(outputSubtomoXmdFilePath)

        firstPredictionArray = []
        secondPredictionArray = []

        for objId in mData:
            firstPredictionArray.append(mData.getValue(MDL_MAX, objId))
            secondPredictionArray.append(mData.getValue(MDL_MIN, objId))

        return firstPredictionArray, secondPredictionArray

    @staticmethod
    def readTomoScores(outputTomoXmdFilePath):
        mData = MetaData()
        mData.read(outputTomoXmdFilePath)

        for objId in mData:
            overallPrediction = mData.getValue(MDL_MAX, objId)
            predictionAverage = mData.getValue(MDL_MIN, objId)

        return overallPrediction, predictionAverage

    def getOutputSetOfAlignedTomograms(self):
        if self.alignedTomograms:
            self.alignedTomograms.enableAppend()

        else:
            alignedTomograms = self._createSetOfTomograms(suffix='Ali')

            alignedTomograms.copyInfo(self.isot)

            alignedTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(alignedTomograms=alignedTomograms)
            self._defineSourceRelation(self.isot, alignedTomograms)

        return self.alignedTomograms

    def getOutputSetOfStrongMisalignedTomograms(self):
        if self.strongMisalignedTomograms:
            self.strongMisalignedTomograms.enableAppend()

        else:
            strongMisalignedTomograms = self._createSetOfTomograms(suffix='StrongMisali')

            strongMisalignedTomograms.copyInfo(self.isot)

            strongMisalignedTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(strongMisalignedTomograms=strongMisalignedTomograms)
            self._defineSourceRelation(self.isot, strongMisalignedTomograms)

        return self.strongMisalignedTomograms

    def getOutputSetOfWeakMisalignedTomograms(self):
        if self.weakMisalignedTomograms:
            self.weakMisalignedTomograms.enableAppend()

        else:
            weakMisalignedTomograms = self._createSetOfTomograms(suffix='WeakMisali')

            weakMisalignedTomograms.copyInfo(self.isot)

            weakMisalignedTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(weakMisalignedTomograms=weakMisalignedTomograms)
            self._defineSourceRelation(self.isot, weakMisalignedTomograms)

        return self.weakMisalignedTomograms

    def getOutputSetOfSubtomos(self):
        if self.outputSubtomos:
            self.outputSubtomos.enableAppend()

        else:
            outputSubtomos = self._createSetOfSubTomograms(suffix="FS")

            outputSubtomos.copyInfo(self.isot)
            outputSubtomos.setSamplingRate(self.targetSamplingRate)

            outputSubtomos.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSubtomos=outputSubtomos)
            self._defineSourceRelation(self.isot, outputSubtomos)

        return self.outputSubtomos

    @staticmethod
    def getSubtomoPathList(coordFilePath):
        coordFilePath_noExt = os.path.splitext(coordFilePath)[0]
        counter = 1

        subtomoPathList = []

        while True:
            subtomoPath = coordFilePath_noExt + '-' + str(counter) + '.mrc'

            if not os.path.exists(subtomoPath):
                break

            subtomoPathList.append(subtomoPath)
            counter += 1

        return subtomoPathList

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = ["Misalignment analysis:"]

        if self.alignedTomograms:
            summary.append("Aligned tomograms: %d"
                           % (self.getOutputSetOfAlignedTomograms().getSize()))

        if self.weakMisalignedTomograms:
            summary.append("Weak misaligned tomograms: %d"
                           % (self.getOutputSetOfWeakMisalignedTomograms().getSize()))

        if self.strongMisalignedTomograms:
            summary.append("Strong misaligned tomograms: %d"
                           % (self.getOutputSetOfStrongMisalignedTomograms().getSize()))

        summary.append("From %d subtomos analyzed."
                       % self.outputSubtomos.getSize())

        return summary

    def _methods(self):
        size = 0
        if self.alignedTomograms:
            size += self.getOutputSetOfAlignedTomograms().getSize()

        if self.weakMisalignedTomograms:
            size += self.getOutputSetOfWeakMisalignedTomograms().getSize()

        if self.strongMisalignedTomograms:
            size += self.getOutputSetOfStrongMisalignedTomograms().getSize()

        methods = ["%d tomograms have been analyzed using the deep_misalignment_detection xmipp method."
                   % size]

        return methods
