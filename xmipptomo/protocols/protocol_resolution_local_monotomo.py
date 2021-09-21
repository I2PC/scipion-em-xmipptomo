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
from pyworkflow.object import Float
from pyworkflow.utils import getExt
from pyworkflow.object import Set
import copy
from pyworkflow.protocol.params import (PointerParam, BooleanParam, FloatParam,
                                        LEVEL_ADVANCED)

from pwem.protocols import ProtAnalysis3D
from pwem.emlib.image import ImageHandler
from pwem.objects import Volume
import pwem.emlib.metadata as md
from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTomograms, Tomogram

MONOTOMO_METHOD_URL = 'http://github.com/I2PC/scipion/wiki/XmippProtMonoTomo'
OUTPUT_RESOLUTION_FILE = 'resolutionMap'
FN_MEAN_VOL = 'meanvol'
METADATA_MASK_FILE = 'metadataresolutions'
FN_METADATA_HISTOGRAM = 'mdhist'
BINARY_MASK = 'binarymask'
FN_GAUSSIAN_MAP = 'gaussianfilter'


class XmippProtMonoTomo(EMProtocol, ProtTomoBase):
    """
    Given a tomogram the protocol assigns local resolutions to each voxel of the tomogram.
    """
    _label = 'local Resolution MonoTomo'
    _lastUpdateVersion = VERSION_2_0

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('oddTomograms', PointerParam, pointerClass='SetOfTomograms',
                      label="Odd tomogram", important=True,
                      help='Select a volume for determining its '
                           'local resolution.')

        form.addParam('evenTomograms', PointerParam, pointerClass='SetOfTomograms',
                      label="Even Tomogram", important=True,
                      help='Select a second volume for determining a '
                           'local resolution.')

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
                                  "to evaluate the local resolution.")

        group.addParam('significance', FloatParam, default=0.95,
                       expertLevel=LEVEL_ADVANCED,
                       label="Significance",
                       help='Relution is computed using hipothesis tests, '
                            'this value determines the significance of that test')

        line.addParam('minRes', FloatParam, default=0, label='High')
        line.addParam('maxRes', FloatParam, allowsNull=True, label='Low')
        line.addParam('stepSize', FloatParam, allowsNull=True, default=0.25,
                      expertLevel=LEVEL_ADVANCED, label='Step')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called """
        myDict = {
            FN_MEAN_VOL: self._getExtraPath('mean_Tomogram.mrc'),
            OUTPUT_RESOLUTION_FILE: self._getExtraPath('TomogramLocalResolution.mrc'),
            METADATA_MASK_FILE: self._getExtraPath('mask_data.xmd'),
            FN_METADATA_HISTOGRAM: self._getExtraPath('hist.xmd'),
            BINARY_MASK: self._getExtraPath('binarized_mask.mrc'),
            FN_GAUSSIAN_MAP: self._getExtraPath('gaussianfilted.mrc'),
        }
        self._updateFilenamesDict(myDict)

    def _insertAllSteps(self):
        #self._createFilenameTemplates()

        print(self.oddTomograms.get())
        print(self.evenTomograms.get())
        counter = 0
        for tom_odd, tom_even in zip(self.oddTomograms.get(), self.evenTomograms.get()):
            counter = counter + 1
            print(counter)

        print(list(zip(self.oddTomograms.get(), self.evenTomograms.get())))

        for tom_odd, tom_even in zip(self.oddTomograms.get(), self.evenTomograms.get()):
            print('ini odd = ', tom_odd.getObjId())
            print('ini even = ', tom_even.getObjId())

            #oid = tom_odd.getObjId().copy()
            #eid = tom_even.getObjId().copy()

            if tom_odd.getObjId() == tom_even.getObjId():
                tomId = tom_odd.getObjId()
                print('id = ', tomId)
                #self.vol1Fn = self.oddTomograms.get()[tomId].getFileName()
                #self.vol2Fn = self.evenTomograms.get()[tomId].getFileName()
                #self._insertFunctionStep('convertInputStep')

                self._insertFunctionStep('resolutionMonoTomoStep', self.oddTomograms.get(), self.evenTomograms.get(), tomId)
                #self._insertFunctionStep("createHistrogram", tomId)
        #self._insertFunctionStep('createOutputStep', tomId)


    def convertInputStep(self):
        """ Read the input volume.
        """
        extVol1 = getExt(self.vol1Fn)
        extVol2 = getExt(self.vol2Fn)
        if (extVol1 == '.mrc') or (extVol1 == '.map'):
            self.vol1Fn = self.vol1Fn + ':mrc'
        if (extVol2 == '.mrc') or (extVol2 == '.map'):
            self.vol2Fn = self.vol2Fn + ':mrc'

        '''
        if self.useMask.get():
            self.maskFn = self.Mask.get().getFileName()
            extMask = getExt(self.maskFn)

            if (extMask == '.mrc') or (extMask == '.map'):
                self.maskFn = self.maskFn + ':mrc'
            #TODO: check if hasValue()
            if self.Mask.hasValue():
                params = ' -i %s' % self.maskFn
                params += ' -o %s' % self._getFileName(BINARY_MASK)
                params += ' --select below %f' % 0.5  # Mask threshold = 0.5 self.maskthreshold.get()
                params += ' --substitute binarize'

                self.runJob('xmipp_transform_threshold', params)
        '''

    def resolutionMonoTomoStep(self, oddTomos, evenTomos, tomId):
        self.vol1Fn = oddTomos[tomId].getFileName()
        self.vol2Fn = evenTomos[tomId].getFileName()
        tomoPath = self._getExtraPath('tomo_' + str(tomId))
        os.mkdir(tomoPath)
        outputlocalResTomoFn = os.path.join(tomoPath, 'localResolutionTomogram_%i.mrc' % tomId)

        # Number of frequencies
        if self.stepSize.hasValue():
            freq_step = self.stepSize.get()
        else:
            freq_step = 0.5

        params = ' --vol %s' % self.vol1Fn
        params += ' --vol2 %s' % self.vol2Fn
        params += ' --meanVol %s' % (os.path.join(tomoPath, 'fullTomogram_%i.mrc' % tomId))
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

        location = os.path.join(tomoPath, 'localResolutionTomogram_%i.mrc' % tomId)

        newTomogram.setLocation(location)

        newTomogram.setSamplingRate(tomo.getSamplingRate())
        outputLocalResolutionSetOfTomograms.append(newTomogram)
        outputLocalResolutionSetOfTomograms.update(newTomogram)
        outputLocalResolutionSetOfTomograms.write()
        self._store()

    def createHistrogram(self, tomId):

        tomoPath = self._getExtraPath('tomo_' + str(tomId))
        fnLocRes = os.path.join(tomoPath, 'localResolutionTomogram_%i.mrc' % tomId)
        m, M = self.getMinMax(fnLocRes)

        freq_step = self.stepSize.get() if self.stepSize.hasValue() else 10

        #M = float(self.max_res_init)
        #m = float(self.min_res_init)
        range_res = round((M - m) / freq_step)

        params = ' -i %s' % fnLocRes

        if self.useMask.get():
            params += ' --mask binary_file %s' % self._getFileName(BINARY_MASK)

        params += ' --steps %f' % range_res
        params += ' --range %f %f' % (m, M-freq_step)
                                      #(float(self.max_res_init.get()) - float(freq_step)))
        params += ' -o %s' % (os.path.join(tomoPath, 'hist_tomogram_%i.xmd' % tomId))

        self.runJob('xmipp_image_histogram', params)


    def readMetaDataOutput(self):
        mData = md.MetaData(self._getFileName(METADATA_MASK_FILE))
        NvoxelsOriginalMask = float(mData.getValue(md.MDL_COUNT, mData.firstObject()))
        NvoxelsOutputMask = float(mData.getValue(md.MDL_COUNT2, mData.firstObject()))
        nvox = int(round(
            ((NvoxelsOriginalMask - NvoxelsOutputMask) / NvoxelsOriginalMask) * 100))
        return nvox

    def getMinMax(self, imageFile):
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        min_res = round(np.amin(imgData) * 100) / 100
        max_res = round(np.amax(imgData) * 100) / 100
        return min_res, max_res

    def getOutputLocalResolutionSetOfTomograms(self):
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
        self.getOutputNormalizedSetOfTomograms().setStreamState(Set.STREAM_CLOSED)
        self._store()

        '''
        volume = Tomogram()
        volume.setFileName(self._getFileName(OUTPUT_RESOLUTION_FILE))

        volume.setSamplingRate(self.oddTomograms.get().getSamplingRate())
        self._defineOutputs(resolution_Volume=volume)
        self._defineSourceRelation(self.oddTomograms, volume)

        # Setting the min max for the summary
        imageFile = self._getFileName(OUTPUT_RESOLUTION_FILE)
        min_, max_ = self.getMinMax(imageFile)
        self.min_res_init.set(round(min_ * 100) / 100)
        self.max_res_init.set(round(max_ * 100) / 100)
        self._store(self.min_res_init)
        self._store(self.max_res_init)
        '''

    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'LocalResolution_Tomogram'):
            messages.append(
                'Information about the method/article in ' + MONOTOMO_METHOD_URL)
        return messages

    def _summary(self):
        summary = []
        summary.append("Highest resolution %.2f Å,   "
                       "Lowest resolution %.2f Å. \n" % (self.min_res_init,
                                                         self.max_res_init))
        return summary

    def _citations(self):
        return ['Vilas2020']
