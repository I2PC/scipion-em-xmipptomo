# **************************************************************************
# *
# * Authors:     JL Vilas (jlvilas@cnb.csic.es)
# *
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
"""
Protocols for tomogram/subtomogram filter operations.
"""

from pyworkflow.object import Float
from pyworkflow.protocol.params import (FloatParam, EnumParam, DigFreqParam,
                                        BooleanParam, PointerParam)

from pwem.objects import ImageDim
from pwem.constants import FILTER_LOW_PASS, FILTER_HIGH_PASS, FILTER_BAND_PASS
from pwem.protocols import ProtFilterParticles, ProtFilterVolumes

from xmipp3.constants import (FILTER_SPACE_FOURIER, FILTER_SPACE_REAL,
                              FILTER_SPACE_WAVELET)
import pyworkflow.utils as pwutils
import pyworkflow.utils.path as path
import os

from .preprocess.protocol_filter_base import (XmippProtFilterBase, FM_LOW_PASS, FM_HIGH_PASS, FM_BAND_PASS)


class XmippProtFilter3DObjects(XmippProtFilterBase):
    """
    This is the class for filtering 2D-Objects (Tomograms, or Volumes)
    """
    _label = 'filter tomos/subtomos'

    def __init__(self, **args):
        XmippProtFilterBase.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('input3D', PointerParam,
                      pointerClass='SetOfVolumes',
                      label='Tomos/subtomos', important=True)

        self._defineProcessParams(form)

    def _insertAllSteps(self):
        for tom in self.input3D.get():
            #self._insertFunctionStep('convertInputStep')
            self._insertFunctionStep('filterStep',tom.getObjId())
            #self._insertFunctionStep('createOutputStep')

    def _convertInputStep(self):
        pass

    def filterStep(self, tomObjId):
        ts = self.input3D.get()[tomObjId]
        tsId = ts.getTsId()


        fnIn = ts.getFileName()
        fnOut = pwutils.removeExt(os.path.basename(fnIn)) + '_filtered.mrc'
        path.makePath(self._getExtraPath(tsId))

        if self.filterSpace == FILTER_SPACE_FOURIER:
            lowFreq = self.lowFreqA.get()
            highFreq = self.highFreqA.get()
            freqDecay = self.freqDecayA.get()
            print("lowFreq, highFreq, freqDecay: ", lowFreq, highFreq, freqDecay)

            mode = self.filterModeFourier.get()

            if mode == FM_LOW_PASS:
                filterStr = " low_pass %f %f " % (highFreq, freqDecay)
            elif mode == FM_HIGH_PASS:
                filterStr = " high_pass %f %f " % (lowFreq, freqDecay)
            elif mode == FM_BAND_PASS:
                filterStr = " band_pass %f %f %f " % (lowFreq, highFreq, freqDecay)
            else:
                raise Exception("Unknown fourier filter mode: %d" % mode)

            args = " --fourier " + filterStr
        else:
            raise Exception("Unknown filter space: %d" % protocol.filterSpace.get())

        fnOut = self._getExtraPath(fnOut)
        args += " -o %s " % fnOut

        args += " -i %s" % fnIn
        self.runJob("xmipp_transform_filter", args)


    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        if not self.input3D.get():
            validateMsgs.append("An input volume is required.")

        return validateMsgs

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = ["The filter uses Xmipp see [Sorzano2007a]."]
        return methods

