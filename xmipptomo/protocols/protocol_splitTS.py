# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
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
# **************************************************************************

import os
import pyworkflow.utils.path as path
import pyworkflow.utils as pwutils
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase


class XmippProtSplitTiltSeries(EMProtocol, ProtTomoBase):
    """
    Wrapper protocol to Xmipp split Odd Even on tilt-series
    """
    _label = 'split tilt-series'

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series',
                      help='Select a set of tilt-series to be split into two sets (odd and even).'
                           'It means, the set of tilt-series is split in two subsets.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep('splitTiltSeries', ts.getObjId())
            self._insertFunctionStep('convertXmdToStackStep', ts.getObjId())
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def splitTiltSeries(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        tsFileName = ts.getFirstItem().getFileName()
        path.makePath(self._getExtraPath(tsId))

        tsFileNameOdd = pwutils.removeExt(os.path.basename(tsFileName)) + "_odd.xmd"
        tsFileNameEven = pwutils.removeExt(os.path.basename(tsFileName)) + "_even.xmd"

        paramsOddEven = {
            'inputImg': tsFileName,
            'outputOdd': self._getExtraPath(os.path.join(tsId, tsFileNameOdd)),
            'outputEven': self._getExtraPath(os.path.join(tsId, tsFileNameEven)),
            'type': "frames",
        }

        argsOddEven = "--img %(inputImg)s " \
                      "-o %(outputOdd)s " \
                      "-e %(outputEven)s " \
                      "--type %(type)s "

        self.runJob('xmipp_image_odd_even', argsOddEven % paramsOddEven)

    def convertXmdToStackStep(self, tsObjId):
        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        tsFileName = ts.getFirstItem().getFileName()

        tsFileNameOdd = pwutils.removeExt(os.path.basename(tsFileName)) + "_odd.xmd"
        tsFileNameEven = pwutils.removeExt(os.path.basename(tsFileName)) + "_even.xmd"

        tsFileNameOddMrc = pwutils.removeExt(os.path.basename(tsFileNameOdd)) + ".mrc"
        tsFileNameEvenMrc = pwutils.removeExt(os.path.basename(tsFileNameEven)) + ".mrc"

        paramsConvertOdd = {
            'inputXmdOdd': self._getExtraPath(os.path.join(tsId, tsFileNameOdd)),
            'outputMrcOdd': self._getExtraPath(os.path.join(tsId, tsFileNameOddMrc)),
        }

        argsConvertOdd = "-i %(inputXmdOdd)s " \
                         "-o %(outputMrcOdd)s "

        self.runJob('xmipp_image_convert', argsConvertOdd % paramsConvertOdd)

        paramsConvertEven = {
            'inputXmdEven': self._getExtraPath(os.path.join(tsId, tsFileNameEven)),
            'outputMrcEven': self._getExtraPath(os.path.join(tsId, tsFileNameEvenMrc)),
        }

        argsConvertEven = "-i %(inputXmdEven)s " \
                          "-o %(outputMrcEven)s "

        self.runJob('xmipp_image_convert', argsConvertEven % paramsConvertEven)

    def createOutputStep(self):
        oddSet = self._createSetOfMovies(suffix='odd')
        evenSet = self._createSetOfMovies(suffix='even')

        for movie in self.inputMovies.get():
            fnMovie = movie.getFileName()

            fnMovieOddMrc = self._getExtraPath(pwutils.removeExt(os.path.basename(fnMovie)) + "_odd.mrc")
            fnMovieEvenMrc = self._getExtraPath(pwutils.removeExt(os.path.basename(fnMovie)) + "_even.mrc")

            imgOutOdd = Movie()
            imgOutEven = Movie()

            imgOutOdd.setFileName(fnMovieOddMrc)
            imgOutEven.setFileName(fnMovieEvenMrc)

            imgOutOdd.setSamplingRate(movie.getSamplingRate())
            imgOutEven.setSamplingRate(movie.getSamplingRate())

            oddSet.append(imgOutOdd)
            evenSet.append(imgOutEven)

        oddSet.copyInfo(self.inputMovies.get())
        evenSet.copyInfo(self.inputMovies.get())

        oddSet.setSamplingRate(self.inputMovies.get().getSamplingRate())
        evenSet.setSamplingRate(self.inputMovies.get().getSamplingRate())

        self._defineOutputs(oddMovie=oddSet)
        self._defineOutputs(evenMovie=evenSet)

        self._defineSourceRelation(self.inputMovies, oddSet)
        self._defineSourceRelation(self.inputMovies, evenSet)

    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'resolution_Volume'):
            messages.append(
                'Information about the method/article in ')
        return messages

    def _summary(self):
        summary = []

        return summary

    def _citations(self):
        return ['']
