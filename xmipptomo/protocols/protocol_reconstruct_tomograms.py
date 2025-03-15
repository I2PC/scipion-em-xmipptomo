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
from pyworkflow.object import Set, Float
from pyworkflow.protocol.params import (PointerParam, FloatParam, EnumParam, LEVEL_ADVANCED)
import pyworkflow.utils.path as path

from pwem.protocols import ProtAnalysis3D
from pwem.emlib.image import ImageHandler
import pwem.emlib.metadata as md
from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram

import pwem.emlib as emlib

MRCEXT = '.mrc'


class XmippProtReconstructTomograms(EMProtocol, ProtTomoBase):
    """
    Given a set of Tilt series with the corresponding alignment parameters. This protocol
    will reconstruct the tomograms associated to the tilt series."""
    _label = 'reconstruct tomogram'
    _lastUpdateVersion = VERSION_2_0

    RECONSTRUCT_ART = 0
    RECONSTRUCT_SIRT = 1

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputTiltSeries', PointerParam, pointerClass='SetOfTiltSeries',
                      label="Tilt Series", important=True,
                      help='Select the Set of Tilt Series that will be used to reconstruct the tomograms.')

        form.addParam('reconstructionMethod',
                      EnumParam,
                      choices=['ART', 'SIRT'],
                      default=self.RECONSTRUCT_ART,
                      label="Reconstruction Algorithm",
                      isplay=EnumParam.DISPLAY_COMBO,
                      help='Select an option to reconstruct tomograms: \n '
                           '_ART_: Arithmetic reconstruction technique. \n'
                           '_SIRT_: (only with MPI) Simultaneous Iterative Reconstruction Technique. \n ')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        for ts in self.inputTiltSeries.get():
            tomId = ts.getObjId()
            self._insertFunctionStep(self.reconstructTomogramStep, tomId)
        self._insertFunctionStep('createOutputStep')

    def createXmdFile(self, ts, fnMd):

        mdTs = md.MetaData()
        idx = 1
        for ti in ts:
            tiltAngle = ti.getTiltAngle()
            fn = ti.getFileName()
            strimg = str(idx) + '@' + fn
            idx = idx + 1

            mdRow = md.Row()
            mdRow.setValue(emlib.MDL_IMAGE, strimg)
            mdRow.setValue(emlib.MDL_ANGLE_PSI, 0.0)
            mdRow.setValue(emlib.MDL_ANGLE_TILT, tiltAngle)
            mdRow.setValue(emlib.MDL_ANGLE_PSI, 0.0)
            mdRow.setValue(emlib.MDL_SHIFT_X, 0.0)
            mdRow.setValue(emlib.MDL_SHIFT_Y, 0.0)
            mdRow.writeToMd(mdTs, mdTs.addObject())

        mdTs.write(fnMd)

        return fnMd


    def reconstructTomogramStep(self, tomId):
        '''
        This function computes the reconstructed tomogram
        '''
        ts = self.inputTiltSeries.get()[tomId]
        tsId = ts.getTsId()

        #Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        path.makePath(tomoPath)

        fnMd = os.path.join(tomoPath, 'tiltseriesAngles.xmd')
        fnMdTs = self.createXmdFile(ts, fnMd)

        #Defining outfiles
        fullTomogramName = self.createOutputPath(tsId, MRCEXT)

        params = ' -i %s' % fnMdTs
        params += ' --parallel_mode %s ' % self.defineReconstructionMethod()
        if self.reconstructionMethod.get() == self.RECONSTRUCT_ART:
            params += ' --thr %i ' % self.numberOfThreads.get()
        params += ' -o %s' % fullTomogramName

        self.runJob('xmipp_reconstruct_art', params)

        newTomogram = Tomogram()
        tomo = self.inputTiltSeries.get()[tomId]
        newTomogram.copyInfo(tomo)
        newTomogram.copyAttributes(tomo, '_origin')

        newTomogram.setLocation(fullTomogramName)

        newTomogram.setSamplingRate(tomo.getSamplingRate())

        outputSetOfTomograms = self.getOutputSetOfTomograms()
        outputSetOfTomograms.append(newTomogram)
        outputSetOfTomograms.update(newTomogram)
        outputSetOfTomograms.write()
        self._store()

    def defineReconstructionMethod(self):
        if self.reconstructionMethod.get() == self.RECONSTRUCT_ART:
            return 'ART'
        elif self.reconstructionMethod.get() == self.RECONSTRUCT_SIRT:
            return 'SIRT'


    def createOutputPath(self, tomId, ext):
        '''
        This function takes a filename as basis, and add the id as suffix and completes
        the path with an extension. Exmaple: filename = 'tomogram_' id = 5, ext = '.mrc'
        the output will be tomogram_5.mrc
        '''
        tomoPath = self._getExtraPath(str(tomId))
        fnPath = os.path.join(tomoPath, str(tomId)+ext)
        return fnPath


    def getOutputSetOfTomograms(self):
        '''
        This function degines the output of the protocol
        '''
        if hasattr(self, "outputSetOfTomograms"):
            self.outputSetOfTomograms.enableAppend()
        else:
            outputSetOfTomograms = self._createSetOfTomograms()
            outputSetOfTomograms.copyInfo(self.inputTiltSeries.get())
            samplingRate = self.inputTiltSeries.get().getSamplingRate()
            outputSetOfTomograms.setSamplingRate(samplingRate)
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTomograms=outputSetOfTomograms)
            self._defineSourceRelation(self.inputTiltSeries, outputSetOfTomograms)
        return self.outputSetOfTomograms

    def createOutputStep(self):
        '''
        This function closes the generated output of the protocol
        '''
        self.getOutputSetOfTomograms().setStreamState(Set.STREAM_CLOSED)
        self._store()

    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'outputSetOfTomograms'):
            messages.append('')
        return messages

    def _summary(self):
        summary = []

        return summary

    def _citations(self):
        return []
