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
from xmipp_base import createMetaDataFromPattern

from tomo.protocols import ProtTomoPicking
from tomo3D.utils import delaunayTriangulation
from tomo.utils import extractVesicles, initDictVesicles

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
        form.addParam('outliersThreshold', params.FloatParam, default=1,
                      label="Outliers distance threshold", condition='outliers == True and filter == 1',
                      help='Z-Score value from 0 to infinite')
        form.addParam('carbon', params.BooleanParam, default=True,
                      label="Score carbon closeness?")
        form.addParam('carbonThreshold', params.FloatParam, default=0.8,
                      label="Carbon distance threshold", condition='carbon == True and filter == 1',
                      help='Score value between 0 and 1')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        pwutils.makePath(self._getExtraPath('inputCoords'))
        pwutils.makePath(self._getExtraPath('outputCoords'))
        coordinates = self.inputCoordinates.get()
        self.tomos = coordinates.getPrecedents()
        self._insertFunctionStep('computeParams', coordinates)
        if self.outliers.get():
            self._insertFunctionStep('detectOutliers')
        if self.carbon.get():
            self._insertFunctionStep('detectCarbonCloseness', coordinates)
        self._insertFunctionStep('createOutputStep', coordinates)

    def computeParams(self, coordinates):
        self.tomo_vesicles, self.tomoNames = initDictVesicles(coordinates)
        for tomoName in self.tomoNames:
            self.tomo_vesicles = extractVesicles(coordinates, self.tomo_vesicles, tomoName)

    def detectOutliers(self):
        self.scoreOutliers = {}
        # threshold = np.exp(-self.outliersThreshold.get())
        self.scoreOutliers = []
        distribution = []
        for tomoName in self.tomoNames:
            for idv, vesicle in enumerate(self.tomo_vesicles[tomoName]['vesicles']):
                tree = cKDTree(vesicle)
                for idp, point in enumerate(vesicle):
                    distance, _ = tree.query(point, k=10)
                    distribution.append(np.mean(distance[1:]))
                    # median_1 = np.abs(distance[1:] - np.mean(distance[1:]))
                    # median_2 = np.mean(median_1)
                    # distribution.append(median_2)
                    self.scoreOutliers.append([self.tomo_vesicles[tomoName]['ids'][idv][idp], 0])

        distribution = np.asarray(distribution)
        z_scores = np.abs((distribution - np.mean(distribution)) / np.std(distribution))
        for idn in range(len(self.scoreOutliers)):
            self.scoreOutliers[idn][1] = z_scores[idn]

    def detectCarbonCloseness(self, coordinates):
        # self.scoreCarbon = {}
        self.scoreCarbon = []
        for tomoName in self.tomoNames:
            idt = self.tomo_vesicles[tomoName]["volId"]
            tomo = self.tomos[idt].clone()
            projFile = self.projectTomo(tomo)
            inputTomoPathMetadataFname = self._getTmpPath("inputTomo.xmd")
            tomo_md = createMetaDataFromPattern(projFile)
            tomo_md.write(inputTomoPathMetadataFname)
            coordList = self.generateCoordList(tomo, coordinates)
            self.writeTomoCoordinates(tomo, coordList, self._getTomoPos(projFile))
            args = '-i %s -c %s -o %s -b %d' \
                   % (inputTomoPathMetadataFname, self._getExtraPath('inputCoords'),
                      self._getExtraPath('outputCoords'), coordinates.getBoxSize())
            self.runJob('xmipp_deep_micrograph_cleaner', args)
            baseName = pwutils.removeBaseExt(projFile)
            outFile = self._getExtraPath('outputCoords', baseName + ".pos")
            posMd = readPosCoordinates(outFile)
            posMd.addLabel(md.MDL_ITEM_ID)
            for objId in posMd:
                if posMd.getValue(md.MDL_ENABLED, objId) == 0:
                    posMd.setValue(md.MDL_ENABLED, 1, objId)
                coord = rowToCoordinate(rowFromMd(posMd, objId))
                self.scoreCarbon.append([posMd.getValue(md.MDL_ITEM_ID, objId),
                                         coord._xmipp_goodRegionScore.get()])
            pwutils.cleanPath(projFile)

    def createOutputStep(self, coordinates):
        outSet = self._createSetOfCoordinates3D(coordinates)
        outSet.setPrecedents(coordinates.getPrecedents())
        outSet.setBoxSize(coordinates.getBoxSize())
        outSet.setSamplingRate(coordinates.getSamplingRate())
        scoreCarbon = dict(self.scoreCarbon) if self.carbon.get() else None
        scoreOutliers = dict(self.scoreOutliers) if self.outliers.get() else None

        for coord in coordinates.iterCoordinates():
            newCoord = coord.clone()
            newCoord.setBoxSize(outSet.getBoxSize())
            newCoord.setVolume(coord.getVolume())
            scoreCoordOutlier = scoreOutliers[coord.getObjId()] if scoreOutliers else None
            scoreCoordCarbon = scoreCarbon[coord.getObjId()] if scoreCarbon else None
            if self.filter.get():
                if self.outliers.get() and not self.carbon.get():
                    if self.outliersThreshold.get() >= scoreCoordOutlier:
                        newCoord.outlierScore = Float(scoreCoordOutlier)
                        outSet.append(newCoord)
                elif not self.outliers.get() and self.carbon.get():
                    if self.carbonThreshold.get() <= scoreCoordCarbon:
                        newCoord.carbonScore = Float(scoreCoordCarbon)
                        outSet.append(newCoord)
                elif self.outliers.get() and self.carbon.get():
                    if self.outliersThreshold.get() >= scoreCoordOutlier and \
                       self.carbonThreshold.get() <= scoreCoordCarbon:
                        newCoord.outlierScore = Float(scoreCoordOutlier)
                        newCoord.carbonScore = Float(scoreCoordCarbon)
                        outSet.append(newCoord)
                else:
                    print("All scoring modes are disabled. Exting")
                    break
            else:
                if self.outliers.get():
                    newCoord.outlierScore = Float(scoreCoordOutlier)
                if self.carbon.get():
                    newCoord.carbonScore = Float(scoreCoordCarbon)
                outSet.append(newCoord)
        self._defineOutputs(outputCoordinates=outSet)
        self._defineSourceRelation(coordinates, outSet)


    # --------------------------- UTILS functions ----------------------
    def projectTomo(self, tomo):
        outFile = pwutils.removeBaseExt(tomo.getFileName()) + '_projected.mrc'
        ih = ImageHandler()
        outProjection = ih.createImage()
        tomoData = np.squeeze(ih.read(tomo.getFileName()).getData())
        projection = np.sum(tomoData, axis=0)
        outProjection.setData(projection)
        ih.write(outProjection, self._getExtraPath(outFile))
        return self._getExtraPath(outFile)

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

    def generateCoordList(self, tomo, coordinates):
        return [coord.clone() for coord in coordinates.iterCoordinates(volume=tomo)]

    # --------------------------- INFO functions ----------------------

