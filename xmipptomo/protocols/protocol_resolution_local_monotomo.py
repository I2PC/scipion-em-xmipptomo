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
from pyworkflow.protocol.params import (PointerParam, BooleanParam, FloatParam,
                                        LEVEL_ADVANCED)

from pwem.protocols import ProtAnalysis3D
from pwem.emlib.image import ImageHandler
import pwem.emlib.metadata as md
from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram

MONOTOMO_METHOD_URL = 'http://github.com/I2PC/scipion/wiki/XmippProtMonoTomo'
TOMOGRAM_RESOLUTION_FILE = 'localResolutionTomogram_'
FULL_TOMOGRAM_FILE = 'fullTomogram_'
HISTOGRAM_RESOLUTION_FILE = 'histogram_resolution_'
BINARY_MASK = 'binarymask'
MRCEXT = '.mrc'
XMDEXT = '.xmd'


class XmippProtMonoTomo(EMProtocol, ProtTomoBase):
    """
    Given a tomogram the protocol assigns local resolutions to each voxel of the tomogram.
    To do that, thje protocol makes use of two half tomograms, called odd and even.
    These tomograms are reconstructed with the same alignment parameter but using the
    half of the data. For instance, the odd/even-images of the tilt series, or much
    better usign the odd/even frames of the movies (recommended). The result is a
    tomogram with the values of local resolution.
    """
    _label = 'local Resolution MonoTomo'
    _lastUpdateVersion = VERSION_2_0

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.min_res_init = Float()
        self.max_res_init = Float()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('oddTomograms', PointerParam, pointerClass='SetOfTomograms',
                      label="Odd tomogram", important=True,
                      help='Select the odd tomogram for determining the '
                           'local resolution tomogram.')

        form.addParam('evenTomograms', PointerParam, pointerClass='SetOfTomograms',
                      label="Even Tomogram", important=True,
                      help='Select the even tomogram for determining the  '
                           'local resolution tomogram.')

        form.addParam('useMask', BooleanParam, default=False,
                      label="Use mask?",
                      help='The mask determines which points are specimen'
                           ' and which are not.')

        form.addParam('Mask', PointerParam, pointerClass='VolumeMask',
                      condition='useMask', allowsNull=True,
                      label="Binary Mask",
                      help='The mask determines which points are specimen'
                           ' and which are not')

        group = form.addGroup('Extra parameters')
        line = group.addLine('Resolution Range (Å)',
                             help="Resolution range (and step in expert mode) "
                                  "to analyze the local resolution.")

        group.addParam('significance', FloatParam, default=0.95,
                       expertLevel=LEVEL_ADVANCED,
                       label="Significance",
                       help='Relution is computed using hypothesis tests, '
                            'this value determines the significance of that test')

        line.addParam('minRes', FloatParam, default=0, label='High')
        line.addParam('maxRes', FloatParam, allowsNull=True, label='Low')
        line.addParam('stepSize', FloatParam, allowsNull=True, default=0.5,
                      expertLevel=LEVEL_ADVANCED, label='Step')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self.min_res_init = Float(self.minRes.get())
        self.max_res_init = Float(self.maxRes.get())
        for tom_odd, tom_even in zip(self.oddTomograms.get(), self.evenTomograms.get()):
            if tom_odd.getObjId() == tom_even.getObjId():
                tomId = tom_odd.getObjId()
                #self._insertFunctionStep('convertInputStep')
                self._insertFunctionStep(self.resolutionMonoTomoStep, self.oddTomograms.get(), self.evenTomograms.get(), tomId)
                self._insertFunctionStep(self.createHistrogram, tomId)
        self._insertFunctionStep('createOutputStep')


    def convertInputStep(self):
        """
        This function takes the input tomograms and if their extension is not mrc,
         then the function convert them into mrc format.
        """
        extVol1 = getExt(self.vol1Fn)
        extVol2 = getExt(self.vol2Fn)
        if (extVol1 == '.mrc') or (extVol1 == '.map'):
            self.vol1Fn = self.vol1Fn + ':mrc'
        if (extVol2 == '.mrc') or (extVol2 == '.map'):
            self.vol2Fn = self.vol2Fn + ':mrc'


    def resolutionMonoTomoStep(self, oddTomos, evenTomos, tomId):
        '''
        This function estimates the local resolution from the oddTomo and the evenTomo.
        The output is generated in pseudo streaming. It is not a full streaming due to
        the input is not updated during the execution.
        '''
        self.vol1Fn = oddTomos[tomId].getFileName()
        self.vol2Fn = evenTomos[tomId].getFileName()

        ts = self.oddTomograms.get()[tomId]
        tsId = ts.getTsId()

        #Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        os.mkdir(tomoPath)

        #Defining outfiles
        outputlocalResTomoFn = self.createOutputPath(TOMOGRAM_RESOLUTION_FILE, tsId, MRCEXT)
        fullTomogramName = self.createOutputPath(FULL_TOMOGRAM_FILE, tsId, MRCEXT)

        # Number of frequencies
        if self.stepSize.hasValue():
            freq_step = self.stepSize.get()
        else:
            freq_step = 0.5

        params = ' --vol %s' % self.vol1Fn
        params += ' --vol2 %s' % self.vol2Fn

        params += ' --meanVol %s' % fullTomogramName
        if self.useMask.get():
            params += ' --mask %s' % self._getFileName(BINARY_MASK)
        params += ' --sampling_rate %f' % self.oddTomograms.get().getSamplingRate()
        params += ' --minRes %f' % self.minRes.get()
        params += ' --maxRes %f' % self.maxRes.get()
        params += ' --step %f' % freq_step
        params += ' -o %s' % outputlocalResTomoFn
        params += ' --significance %f' % self.significance.get()

        self.runJob('xmipp_resolution_monotomo', params)

        outputLocalResolutionSetOfTomograms = self.getOutputLocalResolutionSetOfTomograms()

        newTomogram = Tomogram()
        tomo = self.oddTomograms.get()[tomId]
        newTomogram.copyInfo(tomo)
        newTomogram.copyAttributes(tomo, '_origin')

        newTomogram.setLocation(outputlocalResTomoFn)

        newTomogram.setSamplingRate(tomo.getSamplingRate())
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
        tomoPath = self._getExtraPath(str(tomId))
        fnPath = os.path.join(tomoPath, filename+str(tomId)+ext)
        return fnPath

    def createHistrogram(self, tomId):
        '''
        The histogram of local resolution values of the output tomogram is computed
        '''
        print(tomId)
        ts = self.oddTomograms.get()[tomId]
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
            outputLocalResolutionSetOfTomograms.copyInfo(self.oddTomograms.get())
            samplingRate = self.oddTomograms.get().getSamplingRate()
            outputLocalResolutionSetOfTomograms.setSamplingRate(samplingRate)
            outputLocalResolutionSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputLocalResolutionSetOfTomograms=outputLocalResolutionSetOfTomograms)
            self._defineSourceRelation(self.oddTomograms, outputLocalResolutionSetOfTomograms)
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
