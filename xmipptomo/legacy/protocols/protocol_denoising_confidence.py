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


import os

from pyworkflow import VERSION_2_0
from pyworkflow.utils import getExt
from pyworkflow.object import Set, Float
import pyworkflow.utils.path as path
from pyworkflow.protocol.params import (PointerParam, BooleanParam, FloatParam,
                                        LEVEL_ADVANCED)

from pwem.protocols import ProtAnalysis3D
from pwem.emlib.image import ImageHandler
import pwem.emlib.metadata as md
import pwem.emlib as emlib
from pwem.protocols import EMProtocol
from ..utils import writeMdTiltSeries, xmdToTiltSeries
from pyworkflow.protocol.params import IntParam

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTomograms, Tomogram, SetOfTiltSeries, TiltSeries, TiltImage

TOMOCONF_METHOD_URL = 'http://github.com/I2PC/scipion/wiki/XmippProtMonoTomo'
TOMOGRAM_CONF_FILE = 'confTomogram_'
FULL_TOMOGRAM_FILE = 'fullTomogram_'
FN_CONFIDENCEMAP = '_confidence.mrc'
TOMOGRAMFOLDER = 'tomo_'
BINARY_MASK = 'binarymask'
MRCEXT = '.mrc'
XMDEXT = '.xmd'




class XmippProtConfTomo(EMProtocol, ProtTomoBase):
    """
    This protocol can be used with tilt series or with tomograms
    Tilt series: The algorithm estimates the local probabilities of each signal
    Tomograms:
    """
    _label = 'local confidence maps'
    _lastUpdateVersion = VERSION_2_0

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        #self.outputConfTs = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('hasAssociatedOddEven', BooleanParam,
                      label="Has odd-even associated?", important=True,
                      default=True,
                      help='(True) The odd-even SetOfTiltSeries or the SetOfTomograms are internally associated'
                           'to the input SetOfTiltSeries or SetOfTomograms. (False) The user have to provide '
                           'in the input the odd and even sets.')

        form.addParam('inputSet', PointerParam, pointerClass='SetOfTomograms, SetOfTiltSeries',
                      label="Input TiltSeres/Tomograms", important=True,
                      help='Select the odd tomogram for estimating the '
                           'confidence tomogram.')

        form.addParam('oddInput', PointerParam, pointerClass='SetOfTomograms, SetOfTiltSeries',
                      label="Odd tilt series/tomogram", important=True,
                      condition='not hasAssociatedOddEven',
                      help='Select the odd tilt series/tomogram for estimating the '
                           'confidence tomogram.')

        form.addParam('evenInput', PointerParam, pointerClass='SetOfTomograms, SetOfTiltSeries',
                      label="Even tilt series/tomogram", important=True,
                      condition='not hasAssociatedOddEven',
                      help='Select the even tilt series/tomogram for estimating the  '
                           'confidence map.')

        form.addParam('applyMedian', BooleanParam, default=True,
                      label="median filter",
                      help='The mask determines which points are specimen'
                           ' and which are not.')

        form.addParam('locality', IntParam,
                      label="Locality", default=40,
                      help='Edge of the square local windows where local distribution of noise will be measured.')

        form.addParam('sigmaGauss', FloatParam, default=2,
                      label="sigma",
                      help='The mask determines which points are specimen'
                           ' and which are not.')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):

        if self.hasAssociatedOddEven:
            for ts in self.inputSet.get():
                self._insertFunctionStep(self.confidenceMapStep, ts)
        self._insertFunctionStep(self.closeStreamStep)

    def executeConfidenceMap(self, fnTs, odir, evenfn=None, ):
        """
        This step calls the xmipp binaries
        """
        params = ' --tiltseries %s' % fnTs
        if self.hasAssociatedOddEven:
            params += ' --sampling_rate %f' % self.inputSet.get().getSamplingRate()
        else:
            params += ' --sampling_rate %f' % self.oddinput.get().getSamplingRate()
        params += ' --odir %s' % odir
        if self.applyMedian.get():
            params += ' --medianFilter '

        params += ' --sigmaGauss %f' % self.sigmaGauss.get()

        params += ' --threads %i' % self.numberOfThreads.get()

        self.runJob('xmipp_tomo_confidence_map', params)

    def confidenceMapStep(self, inTs):
        '''
        This function estimates the local confidence maps from the oddTS and the evenTS.
        The output is generated in pseudo streaming. It is not a full streaming due to
        the input is not updated during the execution.
        '''
        ts = inTs
        tsId = ts.getTsId()

        # Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        os.mkdir(tomoPath)


        # Defining the output folder
        tomoPath = self._getExtraPath(tsId)
        fnTs = writeMdTiltSeries(ts, tomoPath, '.xmd')
        odir = self.createOutputPath(TOMOGRAM_CONF_FILE, tsId, MRCEXT)

        # Defining outfiles
        self.executeConfidenceMap(fnTs, odir, None)

        outputSetOfTiltSeriesDenoised = self.getOutputSetOfTiltSeries(inTs, suffixCreate='denoised')
        outputSetOfTiltSeriesConfidence = self.getOutputSetOfTiltSeriesConf(inTs, suffixCreate='confidence')
        fnOut_denoised = os.path.join(odir, 'ts_denoised.xmd')
        fnOut_confidence = os.path.join(odir, 'ts_confidence.xmd')
        newTs_denoised = xmdToTiltSeries(outputSetOfTiltSeriesDenoised, inTs, fnOut_denoised, odir=odir, tsid=tsId, suffix='denoised')
        newTs_conf = xmdToTiltSeries(outputSetOfTiltSeriesConfidence, inTs, fnOut_confidence, odir=odir, tsid=tsId, suffix='conf')
        outputSetOfTiltSeriesDenoised.update(newTs_denoised)
        outputSetOfTiltSeriesDenoised.write()
        outputSetOfTiltSeriesConfidence.update(newTs_conf)
        outputSetOfTiltSeriesConfidence.write()
        self._store()

    def getOutputSetOfTiltSeries(self, inputObj, suffixCreate):
        '''
        This function defines the output of the protocol
        '''
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()
        else:
            outputTs = self._createSetOfTiltSeries(suffix=suffixCreate)
            outputTs.copyInfo(inputObj)
            samplingRate = inputObj.getSamplingRate()
            outputTs.setSamplingRate(samplingRate)
            outputTs.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTiltSeries=outputTs)
            self._defineSourceRelation(inputObj, outputTs)
        return self.outputSetOfTiltSeries

    def getOutputSetOfTiltSeriesConf(self, inputObj, suffixCreate):
        '''
        This function defines the output of the protocol
        '''
        if hasattr(self, "outputSetOfTiltSeriesConfidence"):
            self.outputSetOfTiltSeriesConfidence.enableAppend()
        else:
            outputTs = self._createSetOfTiltSeries(suffix=suffixCreate)
            outputTs.copyInfo(inputObj)
            samplingRate = inputObj.getSamplingRate()
            outputTs.setSamplingRate(samplingRate)
            outputTs.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTiltSeriesConfidence=outputTs)
            self._defineSourceRelation(inputObj, outputTs)
        return self.outputSetOfTiltSeriesConfidence

    def closeStreamStep(self):
        self.getOutputSetOfTiltSeries(self.inputSet.get(), suffixCreate='denoised').setStreamState(Set.STREAM_CLOSED)
        self.getOutputSetOfTiltSeriesConf(self.inputSet.get(), suffixCreate='confidence').setStreamState(Set.STREAM_CLOSED)
        self._store()

    def addMrcExt(self, fn):
        extFn = getExt(fn)
        if (extFn == '.mrc') or (extFn == '.map'):
            fn = fn + ':mrc'

        return fn

    def createOutputPath(self, filename, tomId, ext):
        '''
        This function takes a filename as basis, and add the id as suffix and completes
        the path with an extension. Exmaple: filename = 'tomogram_' id = 5, ext = '.mrc'
        the output will be tomogram_5.mrc
        '''
        tomoPath = self._getExtraPath(str(tomId))
        return tomoPath

    def createOutputStep(self, oddTomos, evenTomos):
        '''
        This function closes the generated output of the protocol
        '''

        if issubclass(type(oddTomos), SetOfTiltSeries) and issubclass(type(evenTomos), SetOfTiltSeries):
            self.getOutputConfidenceSetOfTomograms().setStreamState(Set.STREAM_CLOSED)
        else:
            self.getOutputConfidenceSetOfTiltSeries().setStreamState(Set.STREAM_CLOSED)
        self._store()

    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'confTomo'):
            messages.append(
                'Information about the method/article in ' + TOMOCONF_METHOD_URL)
        return messages

    def _summary(self):
        summary = []

        return summary
