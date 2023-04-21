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

import numpy as np

from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.object import Set, Float
from pyworkflow.protocol import FileParam, PointerParam, EnumParam, FloatParam, BooleanParam
from tomo import constants
from tomo.objects import SetOfCoordinates3D, SetOfTomograms, SubTomogram, Coordinate3D
from tomo.protocols import ProtTomoBase
from tensorflow.keras.models import load_model
from xmipptomo import utils

COORDINATES_FILE_NAME = 'subtomo_coords.xmd'
BOX_SIZE = 32
TARGET_SAMPLING_RATE = 6.25


class XmippProtDeepDetectMisalignment(EMProtocol, ProtTomoBase):
    """
    Wrapper protocol to Xmipp image peak high contrast applied to any volume
    """

    _label = 'detect misalignment from fiducials'
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.alignedTomograms = None
        self.strongMisalignedTomograms = None
        self.weakMisalignedTomograms = None
        self.outputSubtomos = None
        self.outputSubtomosBis = None
        self.isot = None

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSetOfCoordinates',
                      PointerParam,
                      pointerClass=SetOfCoordinates3D,
                      label='Fiducial 3D coordinates',
                      help='3D coordinates in dicating the location of the fiducials (gold beads) in the tomogram. '
                           'These fiducails will be the ones used to study misalignment artifacts over them. '
                           'The coordinate denotes the center of the subtomogram')

        form.addParam('tomoSource', EnumParam,
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
                      help='Threshold value to settle if a tomogram presents weak or strong misalignment.')

        # Only for devel purposes
        form.addParam('firstModelPath', FileParam,
                      label='Input model for first split',
                      help='Input model for first split prediction')

        form.addParam('secondModelPath', FileParam,
                      label='Input model for second split',
                      help='Input model for second split prediction')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self.tomoDict = self.getTomoDict()

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
                                     key,
                                     coordFilePath)
            self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEP functions --------------------------------
    def extractSubtomos(self, key, coordFilePath):
        tomo = self.tomoDict[key]

        outputPath = self._getExtraPath(os.path.join(tomo.getTsId()))
        tomoFn = tomo.getFileName()

        paramsExtractSubtomos = {
            'tomogram': tomoFn,
            'coordinates': coordFilePath,
            'boxsize': BOX_SIZE,
            'threads': 1,  # ***
            'outputPath': outputPath,
            'downsample': TARGET_SAMPLING_RATE / tomo.getSamplingRate(),
        }

        argsExtractSubtomos = "--tomogram %(tomogram)s " \
                              "--coordinates %(coordinates)s " \
                              "--boxsize %(boxsize)d " \
                              "--threads %(threads)d " \
                              "-o %(outputPath)s " \
                              "--downsample %(downsample)f" \
                              "--subtomo "

        self.runJob('xmipp_tomo_extract_subtomograms', argsExtractSubtomos % paramsExtractSubtomos)

    def subtomoPrediction(self, key, coordFilePath):
        tomo = self.tomoDict[key]
        subtomoPathList = self.getSubtomoPathList(coordFilePath)

        subtomoCoordList = utils.retrieveXmipp3dCoordinatesIntoList(coordFilePath, xmdFormat=1)

        totalNumberOfSubtomos = len(subtomoPathList)
        tsId = tomo.getTsId()

        print("Analyzing tomogram " + tsId)
        print("Total number of subtomos: " + str(totalNumberOfSubtomos))

        self.loadModels()

        overallPrediction, predictionAverage, firstPredictionArray, secondPredictionArray = \
            self.makePrediction(subtomoPathList)

        print("For volume id " + str(tsId) + " obtained prediction from " +
              str(len(subtomoPathList)) + " subtomos is " + str(overallPrediction))

        tomo._misaliScore = Float(predictionAverage)
        self.addTomoToOutput(tomo=tomo, overallPrediction=overallPrediction)

        # -------------------------------------------------------------------------------------FS
        # Generate output set of subtomograms with a prediction score
        self.getOutputSetOfSubtomos()

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
            subtomogram.setSamplingRate(TARGET_SAMPLING_RATE)
            subtomogram.setVolName(tomo.getTsId())
            subtomogram._strongMisaliScore = Float(firstPredictionArray[i])
            subtomogram._weakMisaliScore = Float(secondPredictionArray[i])

            self.outputSubtomos.append(subtomogram)
            self.outputSubtomos.write()
            self._store()

        # -------------------------------------------------------------------------------------SS
        # Generate output set of subtomograms with a prediction score
        print(self.outputSubtomosBis)
        self.getOutputSetOfSubtomosBis()
        print(self.outputSubtomosBis)

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
            subtomogram.setSamplingRate(TARGET_SAMPLING_RATE)
            subtomogram.setVolName(tomo.getTsId())
            subtomogram._misaliScore = Float(secondPredictionArray[i])

            self.outputSubtomosBis.append(subtomogram)
            self.outputSubtomosBis.write()
            self._store()

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

        self.outputSubtomosBis.setStreamState(Set.STREAM_CLOSED)
        self.outputSubtomosBis.write()

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

    def makePrediction(self, subtomoPathList):
        """
        :param subtomoPathList: list to every subtomo extracted to be analyzed
        :return: overallPrediction: alignment statement for the whole tomograms obtained from the estimations of each
        subtomo:
            1: strong misalignment (first split negative)
            2: weak misalignment (second split negative). Implies the existence of an input alignment threshold
            3: alignment (second split positive)
        """
        ih = ImageHandler()

        numberOfSubtomos = len(subtomoPathList)

        subtomoArray = np.zeros((numberOfSubtomos, 32, 32, 32), dtype=np.float64)

        for index, subtomo in enumerate(subtomoPathList):
            subtomoDataTmp = ih.read(subtomo)
            subtomoDataTmp = subtomoDataTmp.getData()

            # print("subtomo " + str(index) + " mean " + str(subtomoDataTmp.mean()))
            # print("subtomo " + str(index) + " std " + str(subtomoDataTmp.std()))

            subtomoArray[index, :, :, :] = subtomoDataTmp[:, :, :]

        std = subtomoArray.std()
        mean = subtomoArray.mean()

        subtomoArray = (subtomoArray - mean) / std

        firstPredictionArray = self.firstModel.predict(subtomoArray)

        overallPrediction, predictionAverage = self.determineOverallPrediction(firstPredictionArray, overallCriteria=1)

        if not overallPrediction:
            overallPrediction = 1  # Strong misalignment

            # Set misalignment score to -1 if subtomos removed by the first network
            secondPredictionArray = np.full(firstPredictionArray.shape, -1)

        # print("first prediction array")
        # print(firstPredictionArray)
        # print("first overall prediction " + str(overallPrediction))

        else:
            secondPredictionArray = self.secondModel.predict(subtomoArray)

            overallPrediction, predictionAverage = self.determineOverallPrediction(secondPredictionArray,
                                                                                   overallCriteria=1)

            if self.misaliThrBool.get():  # Using threshold

                if predictionAverage > self.misaliThr.get():
                    overallPrediction = 3  # Alignment
                else:
                    overallPrediction = 2  # Weak misalignment

            # print("second prediction array")
            # print(secondPredictionArray)
            # print("second overall prediction " + str(overallPrediction))

        return overallPrediction, predictionAverage, firstPredictionArray, secondPredictionArray

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

    def loadModels(self):
        self.firstModel = load_model(self.firstModelPath.get())
        print(self.firstModel.summary())

        self.secondModel = load_model(self.secondModelPath.get())
        print(self.firstModel.summary())

    @staticmethod
    def determineOverallPrediction(predictionList, overallCriteria):
        """
        This method return an overall prediction based on the different singular predictions for each gold bead. This
        can be estimated with a voting system (no considering the actual score value) or as the average of the obtained
        scores for each gold beads.
        :param predictionList: vector with the score values predicted for each gold bead
        :param overallCriteria: criteria to be used to calculate the overall prediction as the most voted option (0) or
        the average of all the scores (1)
        :return: bool indicating if the tomogram present misalignment or not
        :return average of the predicted scores
        """

        predictionAvg = np.average(predictionList)

        if overallCriteria == 0:
            predictionClasses = np.round(predictionList)

            overallPrediction = 0

            for prediction in predictionClasses:
                overallPrediction += prediction

            print("Subtomo analysis: " + str(overallPrediction) + " aligned vs " +
                  str(predictionList.size - overallPrediction) + "misaligned")

            overallPrediction = overallPrediction / predictionList.size

            # aligned (1) or misaligned (0)
            return (True if overallPrediction > 0.5 else False), predictionAvg

        elif overallCriteria == 1:
            print("prediction list:")
            print(predictionList)

            print("Subtomo analysis preditcion: " + str(predictionAvg))

            # aligned (1) or misaligned (0)
            return (True if predictionAvg > 0.5 else False), predictionAvg

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
            outputSubtomos.setSamplingRate(TARGET_SAMPLING_RATE)

            outputSubtomos.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSubtomos=outputSubtomos)
            self._defineSourceRelation(self.isot, outputSubtomos)

        return self.outputSubtomos

    def getOutputSetOfSubtomosBis(self):
        if self.outputSubtomosBis:
            self.outputSubtomosBis.enableAppend()

        else:
            outputSubtomosBis = self._createSetOfSubTomograms(suffix="SS")

            outputSubtomosBis.copyInfo(self.isot)
            outputSubtomosBis.setSamplingRate(TARGET_SAMPLING_RATE)

            outputSubtomosBis.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSubtomosBis=outputSubtomosBis)
            self._defineSourceRelation(self.isot, outputSubtomosBis)

        return self.outputSubtomosBis

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
        summary = []

        return summary

    def _methods(self):
        methods = []

        return methods
