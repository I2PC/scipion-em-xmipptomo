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


from pyworkflow import VERSION_2_0
from pyworkflow.object import Float
from pyworkflow.utils import getExt
from pyworkflow.protocol.params import (PointerParam, BooleanParam, FloatParam,
                                        LEVEL_ADVANCED, EnumParam, DigFreqParam)
from pwem.objects import Volume
from pwem.protocols import ProtAnalysis3D

from xmipp3.constants import (FILTER_SPACE_FOURIER, FILTER_SPACE_REAL,
                              FILTER_SPACE_WAVELET)
from pwem.constants import FILTER_LOW_PASS, FILTER_HIGH_PASS, FILTER_BAND_PASS

# Fourier filters
FM_LOW_PASS = FILTER_LOW_PASS  # 0
FM_HIGH_PASS = FILTER_HIGH_PASS  # 1
FM_BAND_PASS = FILTER_BAND_PASS  # 2

# Real Space Filters
FM_MEDIAN = 0
#Wavelets decomposition base
FM_DAUB4   = 0
FM_DAUB12  = 1
FM_DAUB20  = 2

class XmippProtFilterBase(ProtAnalysis3D):
    """
    This is the common class for filtering 2D-images (micrographs or particles)
    and 3D-volumes (Tomograms, or Volumes)
    """
    _lastUpdateVersion = VERSION_2_0

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)

    @classmethod
    def getModesCondition(cls, filterMode, *filterModes):
        return ' or '.join('%s==%d' % (filterMode,fm) for fm in filterModes)

    @classmethod
    def _defineProcessParams(cls, form):
        form.addParam('filterSpace', EnumParam, choices=['fourier', 'real', 'wavelets'],
                      default=FILTER_SPACE_FOURIER,
                      label="Filter space")
        form.addParam('filterModeFourier', EnumParam, choices=['low pass', 'high pass', 'band pass', 'ctf'],
                      default=FM_BAND_PASS,
                      condition='filterSpace == %d' % FILTER_SPACE_FOURIER,
                      label="Filter mode",
                      help='Depending on the filter mode some frequency (freq.) components\n'
                           'are kept and some are removed.\n '
                           '_low pass_: components below *High freq.* are preserved.\n '
                           '_high pass_: components above *Low freq.* are preserved.\n '
                           '_band pass_: components between *Low freq.* and *High freq.* '
                           'are preserved. \n'
                           'ctf: apply first CTF in CTFset to all the particles. This is normally for simulated data.\n'
                           '   : This is not a CTF correction.'
                      )

        form.addParam('filterModeReal', EnumParam, choices=['median'],
                      default=FM_MEDIAN,
                      condition='filterSpace == %d' % FILTER_SPACE_REAL,
                      label="Filter mode",
                      help='median: replace each pixel with the median of neighboring pixels.\n'
                      )
        form.addParam('filterModeWavelets', EnumParam, choices=['daub4', 'daub12', 'daub20'],
                      default=FM_DAUB4,
                      condition='filterSpace == %d' % FILTER_SPACE_WAVELET,
                      label="Filter mode",
                      help='DAUB4: filter using the DAUB4 wavelet transform.\n '
                      )

        # String that identifies filter in Fourier Space
        fourierCondition = 'filterSpace == %d and (%s)' % (FILTER_SPACE_FOURIER,
                                                           cls.getModesCondition('filterModeFourier',
                                                                                 FM_LOW_PASS,
                                                                                 FM_HIGH_PASS,
                                                                                 FM_BAND_PASS))
        # String that identifies filter in Real Space
        realCondition = 'filterSpace == %d and (%s)' % (FILTER_SPACE_REAL,
                                                        cls.getModesCondition('filterModeReal',
                                                                              FM_MEDIAN))
        # String that identifies filter in Real Space
        waveletCondition = 'filterSpace == %d and (%s)' % (FILTER_SPACE_WAVELET,
                                                           cls.getModesCondition('filterModeWavelets',
                                                                                 FM_DAUB4,
                                                                                 FM_DAUB12,
                                                                                 FM_DAUB20))
        # fourier

        form.addParam('freqInAngstrom', BooleanParam, default=True,
                      condition='filterSpace == %d' % FILTER_SPACE_FOURIER,
                      label='Provide resolution in Angstroms?',
                      help='If *Yes*, the resolution values for the filter\n'
                           'should be provided in Angstroms. If *No*, the\n'
                           'values should be in normalized frequencies (between 0 and 0.5).')
        # Resolution in Angstroms (inverse of frequencies)
        line = form.addLine('Resolution (A)',
                            condition=fourierCondition + ' and freqInAngstrom',
                            help='Range of resolutions to use in the filter')
        line.addParam('lowFreqA', FloatParam, default=60,
                      condition='(' + cls.getModesCondition('filterModeFourier',
                                                            FM_BAND_PASS,
                                                            FM_HIGH_PASS) + ') and freqInAngstrom',
                      label='Lowest')
        line.addParam('highFreqA', FloatParam, default=10,
                      condition='(' + cls.getModesCondition('filterModeFourier',
                                                            FM_BAND_PASS,
                                                            FM_LOW_PASS) + ') and freqInAngstrom',
                      label='Highest')

        form.addParam('freqDecayA', FloatParam, default=100,
                      condition=fourierCondition + ' and freqInAngstrom',
                      label='Decay length',
                      help=('Amplitude decay in a [[https://en.wikipedia.org/'
                            'wiki/Raised-cosine_filter][raised cosine]]'))

        # Normalized frequencies ("digital frequencies")
        line = form.addLine('Frequency (normalized)',
                            condition=fourierCondition + ' and (not freqInAngstrom)',
                            help='Range of frequencies to use in the filter')
        line.addParam('lowFreqDig', DigFreqParam, default=0.02,
                      condition='(' + cls.getModesCondition('filterModeFourier',
                                                            FM_BAND_PASS,
                                                            FM_HIGH_PASS) + ') and (not freqInAngstrom)',
                      label='Lowest')
        line.addParam('highFreqDig', DigFreqParam, default=0.35,
                      condition='(' + cls.getModesCondition('filterModeFourier',
                                                            FM_BAND_PASS,
                                                            FM_LOW_PASS) + ') and (not freqInAngstrom)',
                      label='Highest')

        form.addParam('freqDecayDig', FloatParam, default=0.02,
                      condition=fourierCondition + ' and (not freqInAngstrom)',
                      label='Frequency decay',
                      help=('Amplitude decay in a [[https://en.wikipedia.org/'
                            'wiki/Raised-cosine_filter][raised cosine]]'))

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called """
        myDict = {OUTPUT_3DFSC: self._getExtraPath("3dFSC.mrc"),
                  OUTPUT_DIRECTIONAL_FILTER: self._getExtraPath("filteredMap.mrc"),
                  }
        self._updateFilenamesDict(myDict)

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        # Convert input into xmipp Metadata format
        self._insertFunctionStep('convertInputStep')


    def mrc_convert(self, fileName, outputFileName):
        """Check if the extension is .mrc, if not then uses xmipp to convert it
        """
        ext = getExt(fileName)
        if (ext != '.mrc') and (ext != '.map'):
            params = ' -i "%s"' % fileName
            params += ' -o "%s"' % outputFileName
            self.runJob('xmipp_image_convert', params)
            out = outputFileName + ':mrc'
        else:
            out = fileName + ':mrc'
        return out




    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'resolution_Volume'):
            messages.append(
                'Information about the method/article in ')
        return messages

    def _summary(self):
        summary = []
        summary.append(" ")
        return summary

    def _citations(self):
        return ['']
