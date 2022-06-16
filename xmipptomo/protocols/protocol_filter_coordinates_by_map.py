# **************************************************************************
# *
# * Authors:    Jose Luis Vilas (jlvilas@cnb.csic.es)
# *             Federico P. de Isidro-Gomez (fp.deisidro@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************

import os

from pyworkflow import BETA
from pyworkflow.protocol.params import FloatParam, BooleanParam, PointerParam, EnumParam, LEVEL_ADVANCED

from tomo import constants
from tomo.protocols import ProtTomoBase
from tomo.objects import Coordinate3D

import pyworkflow.utils.path as path
from pyworkflow.object import Set, Float
from pwem.protocols import EMProtocol
import pwem.emlib as emlib
import pwem.emlib.metadata as md

from xmipptomo import utils

METADATA_COORDINATES_STATS = 'coordinateStats_'
METADATA_INPUT_COORDINATES = "inputCoordinates"
XMD_EXT = '.xmd'
OUTPUT_XMD_COORS = 'outCoors.xmd'


class XmippProtFilterCoordinatesByMap(EMProtocol, ProtTomoBase):
    '''Filter coordinate by map both given a mask or a resolucion map from a tomogram'''

    _label = 'Filter coordinates by map'
    _devStatus = BETA

    NO_FILTER = 0
    FILTER_AVERAGE = 1
    FILTER_STD = 2

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        ProtTomoBase.__init__(self)

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputCoordinates',
                      PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input 3D coordinates",
                      important=True,
                      help='Select the set of 3D coordinates to be filtered')

        form.addParam('inputSetOfTomograms',
                      PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Input Tomogram",
                      important=True,
                      help='Select the Set Of Tomograms to be used. The coordinates'
                           'make references to their corresponding tomograms, then, the'
                           'statistics of the the enviroment of each coordinates will'
                           'be calculated. Thus it is possible to associate a mean, and'
                           'a standard deviation to each coordinate.')

        form.addParam('radius',
                      FloatParam,
                      default=50,
                      label="Radius",
                      help='Radius of the ball with center at the coordinate')

        form.addParam('filterOption',
                      EnumParam,
                      choices=['No Filter', 'Average', 'Standard Deviation'],
                      default=self.NO_FILTER,
                      label="Filter option",
                      isplay=EnumParam.DISPLAY_COMBO,
                      help='Select an option to filter the coordinates: \n '
                           '_Average_: Filter by Average value. \n'
                           '_StandardDeviation_: Filter by Standard deviation value.')

        form.addParam('averageFilter',
                      FloatParam,
                      allowsNull=True,
                      condition='filterOption==%d' % self.FILTER_AVERAGE,
                      label="Average",
                      help='Average value as threshold')

        form.addParam('stdFilter',
                      FloatParam,
                      allowsNull=True,
                      condition='filterOption==%d' % self.FILTER_STD,
                      label="std",
                      help='std value as threshold')

        form.addParam('thresholdDirection',
                      BooleanParam,
                      condition='filterOption',
                      label="keep greater than the threshold",
                      help='Set true if you want to keep values greater than the threshold. And set false'
                           'if the values lesser than the threshold will be discarded')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.tomos = self.inputCoordinates.get().getPrecedents()

        for tomo in self.inputCoordinates.get().getPrecedents():
            tomoId = tomo.getObjId()
            self._insertFunctionStep(self.generateSideInfo, tomoId)
            self._insertFunctionStep(self.calculatingStatisticsStep, tomoId)
        self._insertFunctionStep(self.createOutputStep)
        self._insertFunctionStep(self.closeOutputSets)

    # --------------------------- STEPS functions ----------------------------
    def generateSideInfo(self, tomoId):
        """ Generates side information and input files to feed the Xmipp filter coordinates algorithm """

        extraPrefix = self._getExtraPath(str(tomoId))
        path.makePath(extraPrefix)

        _ = utils.writeOutputCoordinates3dXmdFile(self.inputCoordinates.get(),
                                                  os.path.join(self._getExtraPath(str(tomoId)),
                                                               METADATA_INPUT_COORDINATES + XMD_EXT),
                                                  tomoId)

    def calculatingStatisticsStep(self, tomId):
        """ Given a tomogram and a set of coordinates, a ball around is considered and
         the statistic of the tomogram inside the ball with center at the coordiante
         are calculated and generated in a metadata """

        fnOut = METADATA_COORDINATES_STATS + str(tomId) + XMD_EXT
        fnInCoord = self._getExtraPath(os.path.join(str(tomId), METADATA_INPUT_COORDINATES + XMD_EXT))

        params = ' --inTomo %s' % self.retrieveMap(tomId).getFileName()
        params += ' --coordinates %s' % fnInCoord
        params += ' --radius %f' % self.radius.get()
        params += ' -o %s' % self._getExtraPath(os.path.join(str(tomId), fnOut))

        self.runJob('xmipp_tomo_filter_coordinates', params)

    def createOutputStep(self):
        """ Generates a Xmipp metadata file (.xmd) containing the statistical information associated to
        the each coordinate with information extracted from the resolution map. """

        self.getOutputSetOfCoordinates3D()

        outputMdFile = self._getExtraPath(OUTPUT_XMD_COORS)

        for tomo in self.tomos:
            tomoId = tomo.getObjId()
            inputMdFileName = METADATA_COORDINATES_STATS + str(tomoId) + XMD_EXT
            inputMdFile = self._getExtraPath(os.path.join(str(tomoId), inputMdFileName))

            mdCoor = md.MetaData()

            tom_fn = self.retrieveMap(tomoId).getFileName()
            x_pos, y_pos, z_pos, avg, std = utils.readXmdStatisticsFile(inputMdFile)

            for i in range(len(x_pos)):
                avg_i = avg[i]
                std_i = std[i]

                if (self.filterOption == self.NO_FILTER) or \
                        (
                                not self.thresholdDirection.get() and self.filterOption == self.FILTER_AVERAGE and self.averageFilter.get() > avg_i) or \
                        (
                                self.thresholdDirection.get() and self.filterOption == self.FILTER_AVERAGE and self.averageFilter.get() < avg_i) or \
                        (
                                not self.thresholdDirection.get() and self.filterOption == self.FILTER_STD and self.stdFilter.get() > std_i) or \
                        (
                                self.thresholdDirection.get() and self.filterOption == self.FILTER_STD and self.stdFilter.get() < std_i):
                    # Fill metadata
                    mdRow = md.Row()
                    mdRow.setValue(emlib.MDL_IMAGE, tom_fn)
                    mdRow.setValue(emlib.MDL_AVG, avg_i)
                    mdRow.setValue(emlib.MDL_STDDEV, std_i)
                    mdRow.setValue(emlib.MDL_XCOOR, x_pos[i])
                    mdRow.setValue(emlib.MDL_YCOOR, y_pos[i])
                    mdRow.setValue(emlib.MDL_ZCOOR, z_pos[i])
                    mdRow.writeToMd(mdCoor, mdCoor.addObject())

                    # Create output object
                    newCoord3D = Coordinate3D()
                    newCoord3D.setVolume(tomo)
                    newCoord3D.setX(x_pos[i], constants.BOTTOM_LEFT_CORNER)
                    newCoord3D.setY(y_pos[i], constants.BOTTOM_LEFT_CORNER)
                    newCoord3D.setZ(z_pos[i], constants.BOTTOM_LEFT_CORNER)

                    newCoord3D._resAvg = Float(avg_i)
                    newCoord3D._resStd = Float(std_i)

                    newCoord3D.setVolId(tomoId)
                    self.outputSetOfCoordinates3D.append(newCoord3D)
                    self.outputSetOfCoordinates3D.update(newCoord3D)

            mdCoor.write(outputMdFile)

            self.outputSetOfCoordinates3D.write()
            self._store()

    def closeOutputSets(self):
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.outputSetOfCoordinates3D.setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions --------------------------------------------
    def retrieveMap(self, tomoId):
        """ This method return a the given mask/resolution map from the input set given the correspondent tomoId. """

        found = False

        for tomo in self.inputSetOfTomograms.get():
            if tomo.getObjId() == tomoId:
                found = True
                return tomo

        if not found:
            raise Exception("Not map found in input set with tomoId with value" + tomoId)

    def getOutputSetOfCoordinates3D(self):
        if hasattr(self, "outputSetOfCoordinates3D"):
            self.outputSetOfCoordinates3D.enableAppend()

        else:
            outputSetOfCoordinates3D = self._createSetOfCoordinates3D(volSet=self.inputSetOfTomograms.get(),
                                                                      suffix='_filtered')

            outputSetOfCoordinates3D.setSamplingRate(self.inputSetOfTomograms.get().getSamplingRate())
            outputSetOfCoordinates3D.setPrecedents(self.inputSetOfTomograms.get())

            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(outputSetOfCoordinates3D=outputSetOfCoordinates3D)
            self._defineSourceRelation(self.inputSetOfTomograms.get(), outputSetOfCoordinates3D)

        return self.outputSetOfCoordinates3D

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        messages = []

        return messages

    def _summary(self):
        summary = []

        return summary
