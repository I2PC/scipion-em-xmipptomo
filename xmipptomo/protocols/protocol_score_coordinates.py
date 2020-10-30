# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
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

import numpy as np
from scipy.spatial import cKDTree

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.object import Float

from pwem.emlib.image import ImageHandler
import pwem.emlib.metadata as md

from xmipp3.convert import (openMd, readPosCoordinates, rowToCoordinate,
                            rowFromMd)

from tomo.protocols import ProtTomoPicking
from tomo.utils import delaunayTriangulation, extractVesicles


class XmippProtScoreCoordinates(ProtTomoPicking):
    '''Scoring and (optional) filtering of coordinates based on different scoring
    functions (normals angle, carbon distance, neighbour distance)'''

    _label = 'score/filter coordinates'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input 3D coordinates", important=True,
                      help='Select the set of 3D coordinates to compare')
        form.addParam('filter', params.EnumParam, choices=['score', 'filter'],
                       default=0, display=params.EnumParam.DISPLAY_HLIST,
                       label='Operation mode',
                       help='Deterimen wheter to retrieve all the coordinates scored or to '
                            'filter out unwanted coordinates based on a threshold')
        form.addParam('outliers', params.BooleanParam, default=True,
                      label="Score outluiers?")
        form.addParam('outliersThreshold', params.FloatParam, default=0.8,
                      label="Outliers distance threshold", condition='outliers == True and filter == 1',
                      help='Score value between 0 and 1')
        form.addParam('angle', params.BooleanParam, default=True,
                      label="Score normals?")
        form.addParam('angleThreshold', params.FloatParam, default=0.8,
                      label="Angle threshold", condition='angle == True and filter == 1',
                      help='Score value between 0 and 1')
        form.addParam('carbon', params.BooleanParam, default=True,
                      label="Score carbon closeness?")
        form.addParam('carbonThreshold', params.FloatParam, default=0.8,
                      label="Carbon distance threshold", condition='outliers == True and filter == 1',
                      help='Score value between 0 and 1')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        pwutils.makePath(self._getExtraPath('inputCoords'))
        pwutils.makePath(self._getExtraPath('outputCoords'))
        coordinates = self.inputCoordinates.get()
        self.tomos = coordinates.getPrecedents()
        self.tomoNames = [pwutils.removeBaseExt(tomo.getFileName()) for tomo in self.tomos]
        self._insertFunctionStep('computeParams', coordinates)
        if self.outliers.get():
            self._insertFunctionStep('detectOutliers')
        if self.carbon.get():
            self._insertFunctionStep('detectCarbonCloseness', coordinates)
        self._insertFunctionStep('createOutputStep', coordinates)

    def computeParams(self, coordinates):
        self.tomo_vesicles = extractVesicles(coordinates)

    def detectOutliers(self):
        self.scoreOutliers = {}
        # threshold = np.exp(-self.outliersThreshold.get())
        for tomoName in self.tomoNames:
            self.scoreOutliers[tomoName] = []
            for idv, vesicle in enumerate(self.tomo_vesicles[tomoName]['vesicles']):
                tree = cKDTree(vesicle)
                for idp, point in enumerate(vesicle):
                    distance, _ = tree.query(point, k=2)
                    self.scoreOutliers[tomoName].append((self.tomo_vesicles[tomoName]['ids'][idv][idp],
                                                         np.exp(-distance[1])))

    def detectCarbonCloseness(self, coordinates):
        self.scoreCarbon = {}
        for tomo, tomoName in zip(self.tomos.iterItems(), self.tomoNames):
            self.scoreCarbon[tomoName] = []
            projFile, dfactor = self.projectTomo(tomo)
            coordList = self.generateCoordList(tomo, coordinates, dfactor)
            self.writeTomoCoordinates(tomo, coordList, self._getTomoPos(projFile))
            args = '-i %s -c %s -o %s -b %d' \
                   % (projFile, self._getExtraPath('inputCoords'),
                      self._getExtraPath('outputCoords'), coordinates.getBoxSize() // dfactor)
            self.runJob('xmipp_deep_micrograph_cleaner', args)
            baseName = pwutils.removeBaseExt(projFile)
            outFile = self._getExtraPath('outputCoords', baseName + ".pos")
            posMd = readPosCoordinates(outFile)
            posMd.addLabel(md.MDL_ITEM_ID)
            for objId in posMd:
                if posMd.getValue(md.MDL_ENABLED, objId) == 0:
                    posMd.setValue(md.MDL_ENABLED, 1, objId)
                coord = rowToCoordinate(rowFromMd(posMd, objId))
                self.scoreCarbon[tomoName].append((objId,
                                                   coord._xmipp_goodRegionScore.get()))

    def createOutputStep(self, coordinates):
        outSet = self._createSetOfCoordinates3D(coordinates)
        outSet.setPrecedents(coordinates.getPrecedents())
        outSet.setBoxSize(coordinates.getBoxSize())
        outSet.setSamplingRate(coordinates.getSamplingRate())

        for tomo, tomoName in zip(self.tomos, self.tomoNames):
            scoreOutliers = sorted(self.scoreOutliers[tomoName], key=lambda tup: tup[0])
            scoreCarbon = sorted(self.scoreCarbon[tomoName], key=lambda tup: tup[0])
            for tupleCarbon, tupleOutlier in zip(scoreCarbon, scoreOutliers):
                newCoord = coordinates[tupleCarbon[0]].clone()
                newCoord.setBoxSize(outSet.getBoxSize())
                newCoord.setVolume(tomo)
                if self.filter.get():
                    if self.outliers.get() and self.outliersThreshold.get() <= tupleOutlier[1]:
                        newCoord.outlierScore = Float(tupleOutlier[1])
                    else:
                        continue
                    if self.carbon.get() and self.carbonThreshold.get() <= tupleCarbon[1]:
                        newCoord.carbonScore = Float(tupleCarbon[1])
                    else:
                        continue
                    outSet.append(newCoord)
                else:
                    if self.outliers.get():
                        newCoord.outlierScore = Float(tupleOutlier[1])
                    if self.carbon.get():
                        newCoord.carbonScore = Float(tupleCarbon[1])
                    outSet.append(newCoord)
        self._defineOutputs(outputCoordinates=outSet)
        self._defineSourceRelation(coordinates, outSet)


    # --------------------------- UTILS functions ----------------------
    def projectTomo(self, tomo):
        dfactor = None
        outFile = pwutils.removeBaseExt(tomo.getFileName()) + '_projected.mrc'
        ih = ImageHandler()
        outProjection = ih.createImage()
        tomoData = np.squeeze(ih.read(tomo.getFileName()).getData())
        projection = np.sum(tomoData, axis=0)
        outProjection.setData(projection)
        ih.write(outProjection, self._getExtraPath(outFile))
        bgDim = max(projection.shape)
        if bgDim > 200:  # Mejor 500
            dfactor = bgDim / 200
            args = '-i %s --step %f --method fourier' % \
                   (self._getExtraPath(outFile), dfactor)
            self.runJob("xmipp_transform_downsample", args)
        return self._getExtraPath(outFile), dfactor

    def _getTomoPos(self, fileName):
        """ Return the corresponding .pos file for a given tomogram. """
        baseName = pwutils.removeBaseExt(fileName)
        return self._getExtraPath('inputCoords', baseName + ".pos")

    def writeTomoCoordinates(self, tomo, coordList, outputFn, isManual=True,
                             getPosFunc=None):
        """ Write the pos file as expected by Xmipp with the coordinates
        of a given tomogram.
        Params:
            tomo: input tomogram.
            coordList: list of (x, y) pairs of the mic coordinates.
            outputFn: output filename for the pos file .
            isManual: if the coordinates are 'Manual' or 'Supervised'
            getPosFunc: a function to get the positions from the coordinate,
                it can be useful for scaling the coordinates if needed.
        """
        if getPosFunc is None:
            getPosFunc = lambda coord: coord.getPosition()

        state = 'Manual' if isManual else 'Supervised'
        f = openMd(outputFn, state)

        for coord in coordList:
            x, y, z = getPosFunc(coord)
            f.write(" %06d   1   %d  %d  %d   %06d \n"
                    % (coord.getObjId(), x, y, 1, tomo.getObjId()))

        f.close()

    def generateCoordList(self, tomo, coordinates, dfactor):
        if dfactor is not None:
            coordList = []
            for coord in coordinates.iterCoordinates(volume=tomo):
                cloned_coord = coord.clone()
                x, y, z = cloned_coord.getPosition()
                cloned_coord.setPosition(x / dfactor,
                                         y / dfactor,
                                         z / dfactor)
                coordList.append(cloned_coord)
            return coordList
        else:
            return [coord.clone() for coord in coordinates.iterCoordinates(volume=tomo)]

    # --------------------------- INFO functions ----------------------

