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
from xmipptomo.protocols.protocol_crop_resize_base import XmippProtCropBase3D

SUFIXCROPPED = '_cropped.mrc'

class XmippProtCropTomograms(XmippProtCropBase3D):
    """
    Protocol to crop tomograms using xmipp_transform_window.
    The protocol allows to change the size of a tomogram/s, by removing the
    borders defined by the users
    """
    _label = 'crop tomograms'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSet',
                      PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Input set of tomograms',
                      help='Select a set of tomograms to be cropped.')

        self._defineParamsReSize2D(form)
        self._defineParamsReSize3D(form)

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        for tomo in self.inputSet.get():
            self._insertFunctionStep(self.cropTomogramsStep, self.inputSet.get(), tomo.getObjId())
        self._insertFunctionStep('closeStreamStep')

    # --------------------------- STEP functions --------------------------------
    def cropTomogramsStep(self, tomoSet, tomId):
        '''
        This function resize the tomograms by means of xmipp_transform_window.
        The output is create in pseudo-streaming. pseudo because the input is
        not open, but the output is updated during the execution of the protocol
        '''
        ts = self.inputSet.get()[tomId]
        tsId = ts.getTsId()

        outputPath = self._getExtraPath(str(tsId))
        # Defining the output folder
        os.mkdir(outputPath)

        inputTomo = tomoSet[tomId].getFileName()
        outputTomo = self.outputTomoFileName(outputPath, tsId, SUFIXCROPPED)

        # Launching the xmipp command
        params = ' -i %s ' % inputTomo
        params += ' -o %s ' % outputTomo
        params += ' --corners %i %i %i %i %i %i ' %(self.xcrop0.get(), self.ycrop0.get(),
                                                 self.zcrop0.get(), self.xcropF.get(),
                                                 self.ycropF.get(), self.zcropF.get())
        params += ' --physical'
        self.runJob("xmipp_transform_window", params)

        # Creating the output in pseudostreaming
        outputresizedSetOfTomograms = self.getOutputSetOfTomograms()
        newTomogram = Tomogram()
        tomo = self.inputSet.get()[tomId]
        newTomogram.copyInfo(tomo)
        newTomogram.copyAttributes(tomo, '_origin')

        newTomogram.setLocation(outputTomo)

        newTomogram.setSamplingRate(self.inputSet.get().getSamplingRate())
        outputresizedSetOfTomograms.append(newTomogram)
        outputresizedSetOfTomograms.update(newTomogram)
        outputresizedSetOfTomograms.write()
        self._store()

    def outputTomoFileName(self, folder, tomId, ext):
        fnPath = os.path.join(folder, str(tomId) + ext)
        return fnPath

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
            outputSetOfTomograms = self._createSetOfTomograms(suffix='cropped')
            outputSetOfTomograms.copyInfo(self.inputSet.get())
            outputSetOfTomograms.setSamplingRate(self.inputSet.get().getSamplingRate())
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTomograms=outputSetOfTomograms)
            self._defineSourceRelation(self.inputSet, outputSetOfTomograms)
        return self.outputSetOfTomograms

    # --------------------------- INFO functions ----------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSetOfTomograms'):
            summary.append("%d tomograms have been resized using the xmipp_transform_window program\n"
                           % self.outputSetOfTomograms.getSize())
        else:
            summary.append("Output not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputSetOfTomograms'):
            methods.append("%d tomograms have been resized using the xmipp_transform_window program\n"
                           % self.outputSetOfTomograms.getSize())
        else:
            methods.append("Output tomograms not ready yet.")
        return methods
