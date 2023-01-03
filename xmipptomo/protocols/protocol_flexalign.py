# **************************************************************************
# *
# * Authors:     David Strelak (davidstrelak@gmail.com)
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

from pyworkflow import BETA
from tomo.protocols import ProtTsCorrectMotion
from xmipp3.protocols import XmippProtFlexAlign
from pathlib import Path


class XmippProtTsFlexAlign(ProtTsCorrectMotion, XmippProtMovieCorr):
    """
    Simple protocol to average TiltSeries movies as basic
    motion correction. It is used mainly for testing purposes.
    """
    _label = 'tiltseries FlexAlign'
    _devStatus = BETA
    evenOddCapable = True
    
    def _defineParams(self, form):
        ProtTsCorrectMotion._defineParams(self, form)
        form.addSection(label="FlexAlign")
        XmippProtMovieCorr._defineAlignmentParams(self, form)

    def _getOutputMicName(self, movie):
        return self._getTiltImageMRoot(movie) + ".mrc"

    def _processTiltImageM(self, workingFolder, tiltImageM, *args):
        # XXX hack - create needed attributes and directories, so that we don't have to change the implementation of the FlexAlign
        dir = self._getOutputMovieFolder(tiltImageM)
        Path(dir).mkdir(parents=True, exist_ok=True)
        self.inputMovies = self.inputTiltSeriesM

        # Get actual user decision about keeping aligned movies
        doSaveMovie = self.doSaveMovie.get()

        # If need to do even/odd
        if self._doSplitEvenOdd():
            self.doSaveMovie.set(True)

        self._processMovie(tiltImageM)

        # Restore user decision about saving movies
        self.doSaveMovie.set(doSaveMovie)

        self.inputMovies = None
