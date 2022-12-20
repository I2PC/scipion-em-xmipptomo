# **************************************************************************
# *
# * Authors:    Carlos Oscar Sorzano (coss@cnb.csic.es)
# *             Daniel Marchán Torres (da.marchan@cnb.csic.es)  -- streaming version
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************
"""
Consensus alignment protocol
"""

import os
from datetime import datetime
from pyworkflow.gui.plotter import Plotter
import numpy as np
from math import ceil
import re
try:
    from itertools import izip
except ImportError:
    izip = zip
from pwem.objects import SetOfMovies, SetOfMicrographs, MovieAlignment, Image

from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.protocol.constants import (STATUS_NEW)
from xmipp3.convert import getScipionObj
from pwem.constants import ALIGN_NONE

# Tomo imports
from tomo.protocols import ProtTsCorrectMotion
from pyworkflow import BETA
from tomo.objects import SetOfTiltSeries, TiltSeriesDict

ACCEPTED = 'Accepted'
DISCARDED = 'Discarded'


class XmippProtTsConsensusAlignment(ProtTsCorrectMotion):
    """
    Protocol to estimate the agreement between different movie alignment
    algorithms in the Global Shifts).
    """

    _label = 'movie ts alignment consensus'
    outputName = 'tsConsensusAlignments'
    _devStatus = BETA


    def _defineParams(self, form):
        form.addSection(label='Input Consensus')
        form.addParam('inputTiltSeriesM1', params.PointerParam,
                      pointerClass='SetOfTiltSeriesM',
                      label="Reference Aligned TS", important=True,
                      help='Select the aligned movies to evaluate (this first set will give the global shifts)')

        form.addParam('inputTiltSeriesM2', params.PointerParam,
                      pointerClass='SetOfTiltSeriesM',
                      label="Secondary Aligned TS",
                      help='Shift to be compared with reference alignment')

        # form.addParam('minConsCorrelation', params.FloatParam, default=-1,
        #               label='Minimum consensus shifts correlation',
        #               help="Minimum value for the consensus correlations between shifts trajectories."
        #                    "\nIf there are noticeable discrepancies "
        #                    "between the two estimations below this correlation, "
        #                    "it will be discarded.")
        #
        # form.addParam('trajectoryPlot', params.BooleanParam, default=False,
        #               label='Global Alignment Trajectory Plot',
        #               help="This will generate a plot for each movie where the reference and the secondary trajectory"
        #                    "will be plot in the same graph with its correlation value.")

        form.addParallelSection(threads=4, mpi=1)

# --------------------------- INSERT steps functions -------------------------
    def _insertAllSteps(self):
        self._initialize()

        inputTsM1, inputTsM2 = self._getInputTs()

        self._ciStepId = self._insertFunctionStep('convertInputStep',
                                                  inputTsM1.getObjId())

        self._ciStepId = self._insertFunctionStep('convertInputStep',
                                                  inputTsM2.getObjId())

        self._insertFunctionStep('createOutputStep', wait=True,
                                 prerequisites=[self._ciStepId])

        self._coStep = self._steps[-1]  # get last step

        self._tsDict1 = TiltSeriesDict(inputTsM1, self._getOutputSet(),
                                       newItemsCallback=self._insertNewSteps,
                                       doneItemsCallback=self._updateOutput)

        self._tsDict2 = TiltSeriesDict(inputTsM2, self._getOutputSet(),
                                       newItemsCallback=self._insertNewSteps,
                                       doneItemsCallback=self._updateOutput)

        self._tsDict1.update()
        self._tsDict2.update()


    # def _insertAllSteps(self):
    #     self.initializeParams()
    #     movieSteps = self._insertNewMovieSteps(self.allMovies1.keys(),
    #                                            self.allMovies2.keys(),
    #                                            self.insertedDict)
    #
    #
    #     self._insertFunctionStep('createOutputStep',
    #                              prerequisites=movieSteps, wait=True)

    def _getInputTs(self):
        """ Return the tiltSeries input object. """
        return self._getInput1TsPointer().get(), self._getInput2TsPointer().get()

    def _getInput1TsPointer(self):
        return self.inputTiltSeriesM1

    def _getInput2TsPointer(self):
        return self.inputTiltSeriesM2

    # def processTiltImageStep(self, tsId, tiltImageId, *args):
    #     """ To be implemented in subclasses. """
    #     pass
    #
    # def processTiltSeriesStep(self, tsId):
    #     """ To be implemented in subclasses. """
    #     pass

    def _processTiltImageM(self, workingFolder, tiltImageM, *args):
        """ This function should be implemented in subclasses to really provide
        the processing step for this TiltSeries Movie.
        Output corrected image (and DW one) should be copied to expected name.
        """
        pass


    def createOutputStep(self):
        pass

    def alignmentCorrelationMovieStep(self, movieId):
        movie1 = self.allMovies1.get(movieId)
        movie2 = self.allMovies2.get(movieId)
        fn1 = movie1.getFileName()
        fn2 = movie1.getFileName()
        movieID1 = movie1.getObjId()
        movieID2 = movie2.getObjId()
        doneFn = self._getMovieDone(movieId)

        if self.isContinued() and self._isMovieDone(movieId):
            self.info("Skipping movie with ID: %s, seems to be done" % movieId)
            return

        # Clean old finished files
        pwutils.cleanPath(doneFn)

        if (movie1 is None) or (movie2 is None):
            print('AlignmentCorrelationMovieStep movie1 or movie2 are None')
            return

        print(fn1)
        print(fn2)
        print(movieID1)
        print(movieID2)

        alignment1 = movie1.getAlignment()
        alignment2 = movie2.getAlignment()
        shiftX_1, shiftY_1 = alignment1.getShifts()
        shiftX_2, shiftY_2 = alignment2.getShifts()

        # Transformation of the shifts to calculate the shifts trajectory correlation
        S1 = np.ones([3, len(shiftX_1)])
        S2 = np.ones([3, len(shiftX_2)])

        S1[0, :] = shiftX_1
        S1[1, :] = shiftY_1
        S2[0, :] = shiftX_2
        S2[2, :] = shiftY_2

        A = np.dot(np.dot(S1, S2.T), np.linalg.inv(np.dot(S2, S2.T)))
        S2_p = np.dot(A, S2)

        S1_cart = np.array([S1[0, :]/S1[2, :], S1[1, :]/S1[2, :]])
        S2_p_cart = np.array([S2_p[0, :] / S2_p[2, :], S2_p[1, :] / S2_p[2, :]])
        rmse_cart = np.sqrt((np.square(S1_cart - S2_p_cart)).mean())
        maxe_cart = np.max(S1_cart - S2_p_cart)
        corrX_cart = np.corrcoef(S1_cart[0, :], S2_p_cart[0, :])[0, 1]
        corrY_cart = np.corrcoef(S1_cart[1, :], S2_p_cart[1, :])[0, 1]
        corr_cart = np.min([corrY_cart, corrX_cart])

        print('Root Mean Squared Error %f' %rmse_cart)
        print('Max Error %f' %maxe_cart)
        print('Correlation X %f' %corrX_cart)
        print('Correlation Y %f' %corrY_cart)
        print('General Corr (min X&Y) %f' %corr_cart)

        if corr_cart >= self.minConsCorrelation.get():
            print('Alignment shift trajectory correlated')
            fn = self._getMovieSelecFileAccepted()
            with open(fn, 'a') as f:
                f.write('%d T\n' % movieID1)

        elif corr_cart < self.minConsCorrelation.get():
            print('Discrepancy in the alignment with correlation %f' %corr_cart)
            fn = self._getMovieSelecFileDiscarded()
            with open(fn, 'a') as f:
                f.write('%d F\n' % movieID1)

        stats_loc = {'shift_corr': corr_cart, 'shift_corr_X': corrX_cart, 'shift_corr_Y': corrY_cart,
                     'max_error': maxe_cart, 'rmse_error': rmse_cart, 'S1_cart': S1_cart, 'S2_p_cart': S2_p_cart}

        self.stats[movieID1] = stats_loc
        self._store()
        # Mark this ctf as finished
        open(doneFn, 'w').close()

    # HASTA AQUI ES ALGO DE TOMO el resto SPA

    def _checkNewOutput(self):
        """ Check for already selected movies and update the output set. """
        # Load previously done items (from text file)
        doneListDiscarded = self._readCertainDoneList(DISCARDED)
        doneListAccepted = self._readCertainDoneList(ACCEPTED)
        # Check for newly done items
        movieListIdAccepted = self._readtMovieId(True)
        movieListIdDiscarded = self._readtMovieId(False)

        newDoneAccepted = [movieId for movieId in movieListIdAccepted
                           if movieId not in doneListAccepted]
        newDoneDiscarded = [movieId for movieId in movieListIdDiscarded
                            if movieId not in doneListDiscarded]

        firstTimeAccepted = len(doneListAccepted) == 0
        firstTimeDiscarded = len(doneListDiscarded) == 0

        allDone = len(doneListAccepted) + len(doneListDiscarded) +\
                  len(newDoneAccepted) + len(newDoneDiscarded)

        # We have finished when there is not more input movies (stream closed)
        # and the number of processed movies is equal to the number of inputs
        maxMovieSize = min(len(self.allMovies1), len(self.allMovies2))
        self.finished = (self.isStreamClosed and allDone == maxMovieSize)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        def readOrCreateOutputs(doneList, newDone, label=''):
            if len(doneList) > 0 or len(newDone) > 0:
                movSet = self._loadOutputSet(SetOfMovies, 'movies'+label+'.sqlite')
                micSet = self._loadOutputSet(SetOfMicrographs, 'micrographs'+label+'.sqlite')
                label = ACCEPTED if label == '' else DISCARDED
                self.fillOutput(movSet, micSet, newDone, label)
                movSet.setSamplingRate(self.samplingRate)
                micSet.setSamplingRate(self.samplingRate)
                micSet.setAcquisition(self.acquisition.clone())
                movSet.setAcquisition(self.acquisition.clone())

                return movSet, micSet
            return None, None

        movieSet, micSet = readOrCreateOutputs(doneListAccepted, newDoneAccepted)
        movieSetDiscarded, micSetDiscarded = readOrCreateOutputs(doneListDiscarded, newDoneDiscarded, DISCARDED)

        if not self.finished and not newDoneDiscarded and not newDoneAccepted:
        # If we are not finished and no new output have been produced
        # it does not make sense to proceed and updated the outputs
        # so we exit from the function here
            return

        def updateRelationsAndClose(movieSet, micSet, first, label=''):
            if os.path.exists(self._getPath('movies'+label+'.sqlite')):
                micsAttrName = 'outputMicrographs'+label
                self._updateOutputSet(micsAttrName, micSet, streamMode)
                self._updateOutputSet('outputMovies'+label, movieSet, streamMode)

                if first:
                    # We consider that Movies are 'transformed' into the Micrographs
                    # This will allow to extend the micrograph associated to a set of
                    # movies to another set of micrographs generated from a
                    # different movie alignment
                    self._defineTransformRelation(self.inputMovies1, micSet)

                micSet.close()
                movieSet.close()

        updateRelationsAndClose(movieSet, micSet, firstTimeAccepted)
        updateRelationsAndClose(movieSetDiscarded, micSetDiscarded, firstTimeDiscarded, DISCARDED)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

    def fillOutput(self, movieSet, micSet, newDone, label):
        if newDone:
            inputMovieSet = self._loadInputMovieSet(self.movieFn1)
            inputMicSet = self._loadInputMicrographSet(self.micsFn)

            for movieId in newDone:
                movie = inputMovieSet[movieId].clone()
                mic = inputMicSet[movieId].clone()

                movie.setEnabled(self._getEnable(movieId))
                mic.setEnabled(self._getEnable(movieId))
                alignment1 = movie.getAlignment()
                shiftX_1, shiftY_1 = alignment1.getShifts()
                setAttribute(mic, '_alignment_corr', self.stats[movieId]['shift_corr'])
                setAttribute(mic, '_alignment_rmse_error', self.stats[movieId]['rmse_error'])
                setAttribute(mic, '_alignment_max_error', self.stats[movieId]['max_error'])
                alignment = MovieAlignment(xshifts=shiftX_1, yshifts=shiftY_1)
                movie.setAlignment(alignment)

                self._writeCertainDoneList(movieId, label)

                if self.trajectoryPlot.get():
                    self._createAndSaveTrajectoriesPlot(movieId)
                    mic.plotCart = Image()
                    mic.plotCart.setFileName(self._getTrajectoriesPlot(movieId))

                movieSet.append(movie)
                micSet.append(mic)

            inputMovieSet.close()
            inputMicSet.close()

    def _loadOutputSet(self, SetClass, baseName):
        """
        Load the output set if it exists or create a new one.
        """
        setFile = self._getPath(baseName)
        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            if (outputSet.__len__() == 0):
                pwutils.path.cleanPath(setFile)
        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        return outputSet

    def _loadInputMovieSet(self, moviesFn):
        self.debug("Loading input db: %s" % moviesFn)
        movieSet = SetOfMovies(filename=moviesFn)
        movieSet.loadAllProperties()
        return movieSet

    def _loadInputMicrographSet(self, micsFn):
        self.debug("Loading input db: %s" % micsFn)
        micSet = SetOfMicrographs(filename=micsFn)
        micSet.loadAllProperties()
        return micSet

    def _summary(self):
        # return message
        pass

    def _validate(self):
        """ The function of this hook is to add some validation before the
        protocol is launched to be executed. It should return a list of
        errors. If the list is empty the protocol can be executed.
        """
        errors = []
        if (self.inputMovies1.get().hasAlignment() == ALIGN_NONE) or \
           (self.inputMovies2.get().hasAlignment() == ALIGN_NONE):
            errors.append("The inputs ( _Input Movies 1_ or _Input Movies 2_ must be aligned before")

        return errors


    # ------------------------------------ Utils functions ------------------------------------
    def _isMovieDone(self, id):
        """ A movie is done if the marker file exists. """
        return os.path.exists(self._getMovieDone(id))

    def _getMovieDone(self, id):
        """ Return the file that is used as a flag of termination. """
        return self._getExtraPath('DONE', 'movie_%06d.TXT' % id)

    def _readDoneList(self):
        """ Read from a file the id's of the items that have been done. """
        doneFile = self._getAllDone()
        doneList = []
        # Check what items have been previously done
        if os.path.exists(doneFile):
            with open(doneFile) as f:
                doneList += [int(line.strip()) for line in f]
        return doneList

    def _getAllDone(self):
        return self._getExtraPath('DONE_all.TXT')

    def _writeDoneList(self, partList):
        """ Write to a text file the items that have been done. """
        with open(self._getAllDone(), 'a') as f:
            for part in partList:
                f.write('%d\n' % part.getObjId())

    def _getMicsPath(self):
        prot1 = self.inputMovies1.getObjValue()  # pointer to previous protocol

        if hasattr(prot1, 'outputMicrographs'):
            path1 = prot1.outputMicrographs.getFileName()
            if os.path.getsize(path1) > 0:
                return path1
        elif hasattr(prot1, 'outputMicrographsDoseWeighted'):
            path2 = prot1.outputMicrographsDoseWeighted.getFileName()
            if os.path.getsize(path2) > 0:
                return path2
        else:
            return None

    def _readCertainDoneList(self, label):
        """ Read from a text file the id's of the items
        that have been done. """
        doneFile = self._getCertainDone(label)
        doneList = []
        # Check what items have been previously done
        if os.path.exists(doneFile):
            with open(doneFile) as f:
                doneList += [int(line.strip()) for line in f]
        return doneList

    def _writeCertainDoneList(self, movieId, label):
        """ Write to a text file the items that have been done. """
        doneFile = self._getCertainDone(label)
        with open(doneFile, 'a') as f:
            f.write('%d\n' % movieId)

    def _createAndSaveTrajectoriesPlot(self, movieId):
        """ Write to a text file the items that have been done. """
        stats = self.stats[movieId]
        fn = self._getExtraPath('global_trajectories_%d' %movieId+'_plot_cart.png')
        shift_X1 = stats['S1_cart'][0, :]
        shift_Y1 = stats['S1_cart'][1, :]
        shift_X2 = stats['S2_p_cart'][0, :]
        shift_Y2 = stats['S2_p_cart'][1, :]
        # ---------------- PLOT -----------------------
        figureSize = (8, 6)
        plotter = Plotter(*figureSize)
        figure = plotter.getFigure()
        ax = figure.add_subplot(111)
        ax.grid()
        ax.axis('equal')
        ax.set_title('Global Alignment Trajectories')
        ax.set_xlabel('X-axis (CorrX:%.3f)' %stats['shift_corr_X'])
        ax.set_ylabel('Y-axis (CorrX:%.3f)' %stats['shift_corr_Y'])
        # Max range of the plot of the two coordinates
        plotRange = max(max(shift_X1) - min(shift_X1),
                        max(shift_Y1) - min(shift_Y1))
        i = 1
        skipLabels = ceil(len(shift_X1) / 10.0)
        for x, y in izip(shift_X1, shift_Y1):
            if i % skipLabels == 0:
                ax.text(x - 0.02 * plotRange, y + 0.02 * plotRange, str(i))
            i += 1

        ax.plot(shift_X1, shift_Y1, color='b', label='reference shifts')
        ax.plot(shift_X2, shift_Y2, color='r', label='target shifts')
        ax.plot(shift_X1, shift_Y1, 'yo')

        # setting the plot windows to properly see the data
        ax.axis([min(shift_X1) - 0.1 * plotRange, max(shift_X1) + 0.1 * plotRange,
                 min(shift_Y1) - 0.1 * plotRange, max(shift_Y1) + 0.1 * plotRange])

        ax.legend()
        plotter.tightLayout()
        plotter.savefig(fn)
        plotter.close()

    def _getTrajectoriesPlot(self, movieId):
        """ Write to a text file the items that have been done. """
        return self._getExtraPath('global_trajectories_%d' %movieId+'_plot_cart.png')

    def _getCertainDone(self, label):
        return self._getExtraPath('DONE_'+label+'.TXT')

    def _getMovieSelecFileAccepted(self):
        return self._getExtraPath('selection-movie-accepted.txt')

    def _getMovieSelecFileDiscarded(self):
        return self._getExtraPath('selection-movie-discarded.txt')

    def _readtMovieId(self, accepted):
        if accepted:
            fn = self._getMovieSelecFileAccepted()
        else:
            fn = self._getMovieSelecFileDiscarded()
        moviesList = []
        # Check what items have been previously done
        if os.path.exists(fn):
            with open(fn) as f:
                moviesList += [int(line.strip().split()[0]) for line in f]
        return moviesList

    def _getEnable(self, movieId):
        fn = self._getMovieSelecFileAccepted()
        # Check what items have been previously done
        if os.path.exists(fn):
            with open(fn) as f:
                for line in f:
                    if movieId == int(line.strip().split()[0]):
                        if line.strip().split()[1] == 'T':
                            return True
                        else:
                            return False

def setAttribute(obj, label, value):
    if value is None:
        return
    setattr(obj, label, getScipionObj(value))

