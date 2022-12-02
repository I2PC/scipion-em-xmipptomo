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

from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram

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
    blabla
    """
    _label = 'local confidence tomogram'
    _lastUpdateVersion = VERSION_2_0

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.min_res_init = Float()
        self.max_res_init = Float()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('oddTomograms', PointerParam, pointerClass='SetOfTomograms, SetOfTiltSeries',
                      label="Odd tomogram", important=True,
                      help='Select the odd tomogram for estimating the '
                           'confidence tomogram.')

        form.addParam('evenTomograms', PointerParam, pointerClass='SetOfTomograms, SetOfTiltSeries',
                      label="Even Tomogram", important=True,
                      help='Select the even tomogram for estimating the  '
                           'confidence tomogram.')

        form.addParam('useMask', BooleanParam, default=False,
                      label="Use mask?",
                      help='The mask determines which points are specimen'
                           ' and which are not.')

        form.addParam('Mask', PointerParam, pointerClass='VolumeMask',
                      condition='useMask', allowsNull=True,
                      label="Binary Mask",
                      help='The mask determines which points are specimen'
                           ' and which are not')
                           
        form.addParam('applyMedian', BooleanParam, default=True,
                      label="median filter",
                      help='The mask determines which points are specimen'
                           ' and which are not.')

        form.addParam('applySmoothingBeforeConfidence', BooleanParam, default=True,
                      label="Gaussian filter before confidence",
                      help='The mask determines which points are specimen'
                           ' and which are not.')

        form.addParam('applySmoothingAfterConfidence', BooleanParam, default=True,
                      label="Gaussian filter after confidence",
                      help='The mask determines which points are specimen'
                           ' and which are not.')

        form.addParam('sigmaGauss', FloatParam, default=2,
                      label="sigma", condition='applySmoothingAfterConfidence or applySmoothingBeforeConfidence',
                      help='The mask determines which points are specimen'
                           ' and which are not.')

        form.addParam('estimateLocalResolution', BooleanParam, default=True,
                      label="Estimate local resolution?",
                      help='The mask determines which points are specimen'
                           ' and which are not.')
                           
        line = form.addLine('Resolution Range (Å)',
                             help="Resolution range (and step in expert mode) "
                                  "to analyze the local resolution.")

        line.addParam('lowRes', FloatParam, default=150, condition='estimateLocalResolution',
                      label="low",
                      help='low resolution')
                           
        line.addParam('resolutionStep', FloatParam, default=2, condition='estimateLocalResolution',
                      label="Step",
                      help='Resolution Steo.')

        form.addParam('significance', FloatParam, default=0.05,
                      label="significance",
                      help='significance.')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        for tom_odd, tom_even in zip(self.oddTomograms.get(), self.evenTomograms.get()):
            if tom_odd.getObjId() == tom_even.getObjId():
                tomId = tom_odd.getObjId()
                self._insertFunctionStep(self.confidenceMapStep, self.oddTomograms.get(), self.evenTomograms.get(), tomId)
        self._insertFunctionStep('createOutputStep')

    def executeConfidenceMap(self, oddfn, evenfn, odir):
        params = ' --odd %s' % oddfn
        params += ' --even %s' % evenfn

        if self.useMask.get():
            params += ' --mask %s' % self._getFileName(BINARY_MASK)
        params += ' --sampling_rate %f' % self.oddTomograms.get().getSamplingRate()
        params += ' --odir %s' % odir
        params += ' --fdr %f' % self.significance.get()
        params += ' --significance %f' % (1-self.significance.get())
        if self.applyMedian.get():
            params += ' --medianFilter '
        if self.applySmoothingBeforeConfidence.get():
            params += ' --applySmoothingBeforeConfidence'
        if self.applySmoothingAfterConfidence.get():
            params += ' --applySmoothingAfterConfidence'

        if self.applySmoothingAfterConfidence.get() or self.applySmoothingBeforeConfidence.get():
            params += ' --sigmaGauss %f' % self.sigmaGauss.get()


        if self.estimateLocalResolution.get():
            params += ' --localResolution'
            params += ' --lowRes %f' % self.lowRes.get()
            params += ' --highRes %f' % 1
            params += ' --step %f' % self.resolutionStep.get()

        params += ' --threads %i' % self.numberOfThreads.get()

        self.runJob('xmipp_tomo_confidence_map', params)


    def confidenceMapStep(self, oddTomos, evenTomos, tomId):
        '''
        This function estimates the local resolution from the oddTomo and the evenTomo.
        The output is generated in pseudo streaming. It is not a full streaming due to
        the input is not updated during the execution.
        '''
        self.vol1Fn = oddTomos[tomId].getFileName()
        self.vol2Fn = evenTomos[tomId].getFileName()

        extVol1 = getExt(self.vol1Fn)
        extVol2 = getExt(self.vol2Fn)
        if (extVol1 == '.mrc') or (extVol1 == '.map'):
            self.vol1Fn = self.vol1Fn + ':mrc'
        if (extVol2 == '.mrc') or (extVol2 == '.map'):
            self.vol2Fn = self.vol2Fn + ':mrc'

        ts = self.oddTomograms.get()[tomId]
        tsId = ts.getTsId()

        #Defining the output folder
        tomoPath = self._getExtraPath(TOMOGRAMFOLDER + tsId)
        os.mkdir(tomoPath)

        #Defining outfiles
        odir = self.createOutputPath(TOMOGRAM_CONF_FILE, tsId, MRCEXT)

        self.executeConfidenceMap(self.vol1Fn, self.vol2Fn, odir)

        outputConfidenceSetOfTomograms = self.getOutputConfidenceSetOfTomograms()

        newTomogram = Tomogram()
        tomo = self.oddTomograms.get()[tomId]
        newTomogram.copyInfo(tomo)
        newTomogram.copyAttributes(tomo, '_origin')

        newTomogram.setLocation(odir)

        newTomogram.setSamplingRate(tomo.getSamplingRate())
        outputConfidenceSetOfTomograms.append(newTomogram)
        outputConfidenceSetOfTomograms.update(newTomogram)
        outputConfidenceSetOfTomograms.write()
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

    def getOutputConfidenceSetOfTomograms(self):
        '''
        This function degines the output of the protocol
        '''
        if hasattr(self, "outputConfidenceSetOfTomograms"):
            self.outputConfidenceSetOfTomograms.enableAppend()
        else:
            outputConfidenceSetOfTomograms = self._createSetOfTomograms(suffix='confTomo')
            outputConfidenceSetOfTomograms.copyInfo(self.oddTomograms.get())
            samplingRate = self.oddTomograms.get().getSamplingRate()
            outputConfidenceSetOfTomograms.setSamplingRate(samplingRate)
            outputConfidenceSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputConfidenceSetOfTomograms=outputConfidenceSetOfTomograms)
            self._defineSourceRelation(self.oddTomograms, outputConfidenceSetOfTomograms)
        return self.outputConfidenceSetOfTomograms

    def createOutputStep(self):
        '''
        This function closes the generated output of the protocol
        '''
        self.getOutputConfidenceSetOfTomograms().setStreamState(Set.STREAM_CLOSED)
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
        if hasattr(self, 'min_res_init') and hasattr(self, 'max_res_init'):
            summary.append("Highest resolution %.2f Å,   "
                           "Lowest resolution %.2f Å. \n" % (self.min_res_init,
                                                             self.max_res_init))
        return summary

    def _citations(self):
        return ['Vilas2020']
