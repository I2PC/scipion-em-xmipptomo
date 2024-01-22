# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
# *
# * your institution
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

from pyworkflow.object import Set

from pwem.protocols import EMProtocol

from tomo.protocols.protocol_base import ProtTomoBase
from tomo.objects import (SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries,
                          CTFTomoSeries, CTFTomo, SetOfTiltSeriesCoordinates,
                          TiltSeries, TiltImage)

OUTPUT_TOMOGRAMS_NAME = "Tomograms"

class xTomoIO(EMProtocol, ProtTomoBase):
    def __init__(self, **args):

        # Possible outputs (synchronize these names with the constants)
        self.TiltSeriesCoordinates = None
        self.FiducialModelNoGaps = None
        self.FiducialModelGaps = None
        self.TiltSeries = None
        self.InterpolatedTiltSeries = None
        self.CTFTomoSeries = None
        self.FailedTiltSeries = None
        self.Tomograms = None
        self.Coordinates3D = None

    def getOutputSetOfTomograms(self, inputSet, samplingRate):

        if self.Tomograms:
            getattr(self, OUTPUT_TOMOGRAMS_NAME).enableAppend()
        else:
            outputSetOfTomograms = self._createSetOfTomograms()

            if isinstance(inputSet, SetOfTomograms):
                outputSetOfTomograms.copyInfo(inputSet)

            elif isinstance(inputSet, SetOfTiltSeries):
                outputSetOfTomograms.setAcquisition(inputSet.getAcquisition())
                outputSetOfTomograms.setSamplingRate(inputSet.getSamplingRate())

            outputSetOfTomograms.setSamplingRate(samplingRate)
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{OUTPUT_TOMOGRAMS_NAME: outputSetOfTomograms})
            self._defineSourceRelation(inputSet, outputSetOfTomograms)

        return self.Tomograms
