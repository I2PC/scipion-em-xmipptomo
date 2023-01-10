# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csi.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
# *********************************************************************

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol.params import EnumParam, FloatParam, IntParam, BooleanParam, PointerParam
import pyworkflow.protocol.constants as const
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfTomograms, Tomogram
import os
import pyworkflow.utils.path as path
from pyworkflow.object import Set
import tomo.objects as tomoObj
from pwem.emlib.image import ImageHandler
from xmipptomo.protocols.protocol_crop_resize_base import XmippProtResizeBase



class XmippProtResizeTomograms(XmippProtResizeBase):
    """
    Protocol to to resize tomograms using xmipp_image_resize.
    The protocol allows to change the size of a tomogram/s by means
    of different methods
    """

    TOMOGRAMFOLDER = 'tomo_'
    SUFIXRESIZE = '_resized.mrc'

    _label = 'resize tomograms'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSet',
                      PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms',
                      help='Select a set of tomograms to be resized.')
        self._defineParamsReSize(form)

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        for tomo in self.inputSet.get():
            self._insertFunctionStep(self.resizeTomograms, self.inputSet.get(), tomo.getObjId())
        self._insertFunctionStep('closeStreamStep')

    # --------------------------- STEP functions --------------------------------
    def resizeTomograms(self, tomoSet, tomId):
        '''
        This function resize the tomograms by means of xmipp_image_resize.
        The output is create in pseudo-streaming. pseudo because the input is
        not open, but the output is updated during the execution of the protocol
        '''
        ts = self.inputSet.get()[tomId]
        tsId = ts.getTsId()

        #Defining the output folder
        os.mkdir(self._getExtraPath(str(tsId)))

        inputTomo = tomoSet[tomId].getFileName()
        outputTomo = self.outputTomoFileName(tomoSet, tsId)

        # Launching the xmipp command
        params =  ' -i %s ' % inputTomo
        params += ' -o %s ' % outputTomo
        samplingRate = self.inputSet.get().getSamplingRate()
        params += self.resizeCommonArgsResize(self, samplingRate)

        self.runJob("xmipp_image_resize", params)

        # Creating the output in pseudostreaming
        outputresizedSetOfTomograms = self.getOutputSetOfTomograms()
        newTomogram = Tomogram()
        tomo = self.inputSet.get()[tomId]
        newTomogram.copyInfo(tomo)
        newTomogram.setLocation(outputTomo)

        newTomogram.setSamplingRate(self.samplingRate)
        outputresizedSetOfTomograms.append(newTomogram)
        outputresizedSetOfTomograms.update(newTomogram)
        outputresizedSetOfTomograms.write()
        self._store()

    def outputTomoFileName(self, tomoSet, tomId):
        tomoPath = self._getExtraPath(str(tomId))
        auxfn = os.path.basename(tomoSet.getFileName())
        outputTomo = os.path.splitext(auxfn)[0] + str(tomId) + self.SUFIXRESIZE
        return os.path.join(tomoPath, outputTomo)

    def closeStreamStep(self):
        self.getOutputSetOfTomograms().setStreamState(Set.STREAM_CLOSED)
        self._store()

    def getOutputSetOfTomograms(self):
        '''
        This function defines the output of the protocol
        '''
        if hasattr(self, "outputSetOfTomograms"):
            self.outputSetOfTomograms.enableAppend()
        else:
            outputSetOfTomograms = self._createSetOfTomograms(suffix='resized')
            outputSetOfTomograms.copyInfo(self.inputSet.get())
            outputSetOfTomograms.setSamplingRate(self.samplingRate)
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTomograms=outputSetOfTomograms)
            self._defineSourceRelation(self.inputSet, outputSetOfTomograms)
        return self.outputSetOfTomograms

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfTomograms'):
            summary.append("Input Tomograms: %d.\n"
                           "Target sampling rate: %.2f A/px.\n"
                           % (self.inputSet.get().getSize(),
                              self.outputSetOfTomograms.getSamplingRate()))
        else:
            summary.append("Output not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfTomograms'):
            methods.append("%d tomograms have been resized using the xmipp_image_resize program to a %.2f A/px "
                           "target sampling rate.\n"
                           % (self.outputSetOfTomograms.getSize(),
                              self.outputSetOfTomograms.getSamplingRate()))
        else:
            methods.append("Output tomograms not ready yet.")
        return methods
