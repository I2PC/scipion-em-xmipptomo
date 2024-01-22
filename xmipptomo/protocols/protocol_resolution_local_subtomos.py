# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import numpy as np
import os

from pyworkflow import VERSION_2_0
from pyworkflow.utils import getExt
from pyworkflow.object import Set, Float
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam,
                                        LEVEL_ADVANCED)

from pwem.protocols import ProtAnalysis3D
from pwem.emlib.image import ImageHandler
import pwem.emlib.metadata as md
from pwem.protocols import EMProtocol
import tomo.constants as const
from tomo.objects import Tomogram

from pwem.emlib import lib
from tomo.protocols import ProtTomoBase

MONOTOMO_METHOD_URL = 'http://github.com/I2PC/scipion/wiki/XmippProtMonoTomo'
TOMOGRAM_RESOLUTION_FILE = 'localResolutionTomogram_'
FULL_TOMOGRAM_FILE = 'fullTomogram_'
TOMOGRAMFOLDER = 'tomo_'
HISTOGRAM_RESOLUTION_FILE = 'histogram_resolution_'
BINARY_MASK = 'binarymask'
MRCEXT = '.mrc'
XMDEXT = '.xmd'

class XmippProtMonoTomoSubtomos(EMProtocol, ProtTomoBase):
    """
    Given a tomogram the protocol assigns local resolutions to each voxel of the tomogram.
    To do that, the protocol makes use of two half tomograms, called odd and even.
    These tomograms are reconstructed with the same alignment parameter but using the
    half of the data. For instance, the odd/even-images of the tilt series, or much
    better usign the odd/even frames of the movies (recommended). The result is a
    tomogram with the values of local resolution.
    """
    _label = 'resolution MonoTomo subtomos'
    _lastUpdateVersion = VERSION_2_0

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('coords',
                      PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label='Coordinates',
                      help='3D coordinates of the subtomograms for which the local resolution will be estimated.'
                           'The coordinate denotes the center of the subtomogram.')

        form.addParam('tomogram', PointerParam, pointerClass='SetOfTomograms',
                      label="Tomogram", important=True,
                      help='Select the tomogram where the subtomograms are located. If the tomogram'
                           'has odd/even halves associated they will be used. Otherwise'
                           'the full tomogram will be used.')

        form.addParam('boxsize', IntParam, pointerClass='SetOfTomograms',
                      label="Boxsize (px)", important=True,
                      help='The subtomograms are extracted as a cube. The box size defines the edge of the cube. Local'
                           'resolution will be estimated in these cubes')

        group = form.addGroup('Extra parameters')
        line = group.addLine('Resolution Range (Å)',
                             help="Resolution range (and step in expert mode) "
                                  "to analyze the local resolution.")

        group.addParam('significance', FloatParam, default=0.95,
                       expertLevel=LEVEL_ADVANCED,
                       label="Significance",
                       help='Resolution is computed using hypothesis tests, '
                            'this value determines the significance of that test')

        line.addParam('minRes', FloatParam, default=0, label='High')
        line.addParam('maxRes', FloatParam, allowsNull=True, label='Low')
        line.addParam('stepSize', FloatParam, allowsNull=True, default=2.0,
                      expertLevel=LEVEL_ADVANCED, label='Step')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        for t in self.tomogram.get():
            tomId = t.getObjId()
            self._insertFunctionStep(self.resolutionMonoTomoStep, tomId)

        self._insertFunctionStep(self.createOutputStep)

    def writeMdCoordinates(self, tomo, tomoPath):
        """
            Returns the filename of a metadata with the coordinates.
        """
        mdCoor = lib.MetaData()

        tsid = tomo.getTsId()

        for item in self.coords.get().iterCoordinates(volume=tomo):
            coord = item

            if coord.getTomoId() == tsid:
                nRow = md.Row()
                coord.setVolume(tomo)
                nRow.setValue(lib.MDL_XCOOR, int(coord.getX(const.BOTTOM_LEFT_CORNER)))
                nRow.setValue(lib.MDL_YCOOR, int(coord.getY(const.BOTTOM_LEFT_CORNER)))
                nRow.setValue(lib.MDL_ZCOOR, int(coord.getZ(const.BOTTOM_LEFT_CORNER)))
                nRow.setValue(lib.MDL_PARTICLE_ID, int(coord.getObjId()))

                #alignmentToRow(transform, nRow, ALIGN_PROJ)
                nRow.addToMd(mdCoor)

                #newCoord = item.clone()
                #newCoord.setVolume(coord.getVolume())

        fnCoor = os.path.join(tomoPath, "%s.xmd" % tsid)
        mdCoor.write(fnCoor)

        return fnCoor

    def resolutionMonoTomoStep(self, tomId):
        '''
        This function estimates the local resolution from the oddTomo and the evenTomo.
        The output is generated in pseudo streaming. It is not a full streaming due to
        the input is not updated during the execution.
        '''
        mytomo = self.tomogram.get()[tomId]
        if self.tomogram.get().hasOddEven():
            fnTom, fnHalf = mytomo.getHalfMaps().split(',')
        else:
            fnTom = mytomo.getFileName()

        tsId = mytomo.getTsId()

        #Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        if not os.path.exists(tomoPath):
            os.mkdir(tomoPath)

        sampling = mytomo.getSamplingRate()

        #Defining outfiles
        outputlocalResTomoFn = tsId+'_localResTomo.mrc'

        # Number of frequencies
        if self.stepSize.hasValue():
            freq_step = self.stepSize.get()
        else:
            freq_step = 2.0

        fnCoor = self.writeMdCoordinates(mytomo, tomoPath)

        params = ' --tomo %s' % fnTom
        if self.tomogram.get().hasOddEven():
            params += ' --half %s' % fnHalf
        params += ' --coordinates %s ' %fnCoor

        params += ' --sampling %f' % sampling
        params += ' --boxsize %i' % self.boxsize.get()
        params += ' --lowRes %f' % self.maxRes.get()
        params += ' --highRes %f' % self.minRes.get()
        params += ' --resStep %f' % freq_step
        params += ' -o %s' % self._getExtraPath(tsId)
        params += ' --threads %f' % 4

        self.runJob('xmipp_tomo_resolution_subtomo', params)

        outputLocalResolutionSetOfTomograms = self.getOutputLocalResolutionSetOfTomograms()

        newTomogram = Tomogram()
        newTomogram.setLocation(os.path.join(self._getExtraPath(tsId), outputlocalResTomoFn))
        newTomogram.copyInfo(mytomo)
        newTomogram.copyAttributes(mytomo, '_origin')

        newTomogram.setSamplingRate(mytomo.getSamplingRate())
        outputLocalResolutionSetOfTomograms.append(newTomogram)
        outputLocalResolutionSetOfTomograms.update(newTomogram)
        outputLocalResolutionSetOfTomograms.write()
        self._store()


    def createOutputPath(self, filename, tomId, ext):
        '''
        This function takes a filename as basis, and add the id as suffix and completes
        the path with an extension. Exmaple: filename = 'tomogram_' id = 5, ext = '.mrc'
        the output will be tomogram_5.mrc
        '''
        tomoPath = self._getExtraPath(TOMOGRAMFOLDER + str(tomId))
        fnPath = os.path.join(tomoPath, filename+str(tomId)+ext)
        return fnPath

    def createHistrogram(self, tomId):
        '''
        The histogram of local resolution values of the output tomogram is computed
        '''

        ts = self.tomogram.get()[tomId]
        tsId = ts.getTsId()
        fnLocRes = self.createOutputPath(TOMOGRAM_RESOLUTION_FILE, tsId, MRCEXT)
        fnHist = self.createOutputPath(HISTOGRAM_RESOLUTION_FILE, tsId, XMDEXT)
        m, M = self.getMinMax(fnLocRes)

        freq_step = self.stepSize.get() if self.stepSize.hasValue() else 10

        range_res = round((M - m) / freq_step)

        params = ' -i %s' % fnLocRes

        if self.useMask.get():
            params += ' --mask binary_file %s' % self._getFileName(BINARY_MASK)

        params += ' --steps %f' % range_res
        params += ' --range %f %f' % (m, M-freq_step)
        params += ' -o %s' % fnHist

        self.runJob('xmipp_image_histogram', params)

    def getMinMax(self, imageFile):
        '''
        This function computes the maximum and minimum values of a tomogram given by an
        imageFile
        '''
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        min_res = round(np.amin(imgData) * 100) / 100
        max_res = round(np.amax(imgData) * 100) / 100
        return min_res, max_res

    def getOutputLocalResolutionSetOfTomograms(self):
        '''
        This function degines the output of the protocol
        '''
        if hasattr(self, "outputLocalResolutionSetOfTomograms"):
            self.outputLocalResolutionSetOfTomograms.enableAppend()
        else:
            outputLocalResolutionSetOfTomograms = self._createSetOfTomograms(suffix='LocalResolution')
            outputLocalResolutionSetOfTomograms.copyInfo(self.tomogram.get())
            samplingRate = self.tomogram.get().getSamplingRate()
            outputLocalResolutionSetOfTomograms.setSamplingRate(samplingRate)
            outputLocalResolutionSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputLocalResolutionSetOfTomograms=outputLocalResolutionSetOfTomograms)
            self._defineSourceRelation(self.tomogram, outputLocalResolutionSetOfTomograms)
        return self.outputLocalResolutionSetOfTomograms

    def createOutputStep(self):
        '''
        This function closes the generated output of the protocol
        '''
        self.getOutputLocalResolutionSetOfTomograms().setStreamState(Set.STREAM_CLOSED)
        self._store()

    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'LocalResolution_Tomogram'):
            messages.append(
                'Information about the method/article in ' + MONOTOMO_METHOD_URL)
        return messages

    def _summary(self):
        summary = []
        if hasattr(self, 'min_res_init') and hasattr(self, 'max_res_init'):
            summary.append("Highest resolution %.2f Å,   "
                           "Lowest resolution %.2f Å. \n" % (self.min_res_init,
                                                             self.max_res_init))
        return summary

    def _citations(self):
        return ['Vilas2020']
