# **************************************************************************
# *
# * Authors:    Jose Luis Vilas (jlvilas@cnb.csic.es)
# * Authors:    Federico P. de Isidro-Gomez (fp.deisidro@cnb.csic.es)
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
from pyworkflow.protocol.params import FloatParam, PointerParam
from pyworkflow.protocol import params


from tomo.protocols import ProtTomoBase, ProtImport, ProtTomoPicking

import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol

from xmipptomo import Plugin
from xmipptomo import utils

METADATA_COORDINATES_STATS = 'coordinateStats_'
METADATA_INPUT_COORDINATES = "inputCoordinates"
XMD_EXT = '.xmd'


class XmippProtFilterCoordinatesByMap(ProtTomoPicking):
    '''Filter coordinate by map both given a mask or a resolucion map from a tomogram'''

    _label = 'Filter coordinates by map'
    _devStatus = BETA

    # def __init__(self, **args):
    #     EMProtocol.__init__(self, **args)
    #     ProtTomoBase.__init__(self)

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputCoordinates',
                      params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input 3D coordinates",
                      important=True,
                      help='Select the set of 3D coordinates to be filtered')

        form.addParam('inputSetOfTomograms',
                      params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Input Tomogram",
                      important=True,
                      help='Select the Set Of Tomograms to be used. The coordinates'
                           'make references to their corresponding tomograms, then, the'
                           'statistics of the the enviroment of each coordinates will'
                           'be calculated. Thus it is possible to associate a mean, and'
                           'a standard deviation to each coordinate.')

        form.addParam('radius',
                      params.FloatParam,
                      default=50,
                      label="Radius",
                      help='Radius of the ball with center at the coordinate')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        print("sssssssssssssssssssssssssssssssssssssssss")

        self.tomos = self.inputCoordinates.get().getPrecedents()

        print("sssssssssssssssssssssssssssssssssssssssss")

        for tomo in self.inputCoordinates.get().getPrecedents():#self.tomos:
            tomoId = tomo.getObjId()
            self._insertFunctionStep(self.generateSideInfo, tomoId)
            self._insertFunctionStep(self.calculatingStatisticsStep, tomoId)

    # --------------------------- STEPS functions ----------------------------
    def generateSideInfo(self, tomoId):
        """ Generates side information and input files to feed the Xmipp filter coordinates algorithm """

        extraPrefix = self._getExtraPath(str(tomoId))
        path.makePath(extraPrefix)


        utils.writeOutputCoordinates3dXmdFile(self.inputCoordinates.get(),
                                              os.path.join(self._getExtraPath(str(tomoId)),
                                                           METADATA_INPUT_COORDINATES + XMD_EXT),
                                              tomoId)

    def calculatingStatisticsStep(self, tomId):
        """ Given a tomogram and a set of coordinates, a ball around is considered and
         the statistic of the tomogram inside the ball with center at the coordiante
         are calculated and generated in a metadata """

        fnOut = METADATA_COORDINATES_STATS+ str(tomId) + XMD_EXT
        fnInCoord = self._getExtraPath(os.path.join(str(tomId), METADATA_INPUT_COORDINATES + XMD_EXT))

        params = ' --inTomo %s' % self.retrieveMap(tomId).getFileName()
        params += ' --coordinates %s' % fnInCoord
        params += ' --radius %f' % self.radius.get()
        params += ' -o %s' % self._getExtraPath(os.join.path(str(tomId), fnOut))

        self.runJob('xmipp_tomo_filter_coordinates', params)

    # --------------------------- UTILS functions --------------------------------------------
    def retrieveMap(self, tomoId):
        """ This method return a the given mask/resolution map from the input set given the correspondent tomoId. """

        found = False

        for tomo in self.inputSetOfTomograms.get():
            if tomo.getVolId() == tomoId:
                found = True
                return tomo

        if not found:
            raise Exception("Not map found in input set with tomoId with value" + tomoId)

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        messages = []

        return messages

    def _summary(self):
        summary = []

        return summary
