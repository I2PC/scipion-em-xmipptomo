# **************************************************************************
# *
# * Authors:    Mikel Iceta Tena (miceta@cnb.csic.es)
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
# * Initial release: june 2023
# **************************************************************************

"""
Deep Consensus picking protocol suited for Cryo-ET
"""
import os, time

import xmipp3
from xmipp3 import XmippProtocol
from xmipp3.protocols.protocol_pick_noise import pickNoise_prepareInput, IN_COORDS_POS_DIR_BASENAME
from xmipp3.convert import readSetOfParticles, writeSetOfParticles
from xmipp3.convert import readSetOfCoordinates, writeSetOfCoordinates
from xmipp3.convert import setXmippAttributes

from pwem import emlib
from pwem.protocols import ProtParticlePicking, ProtUserSubSet


from pyworkflow import VERSION_3_0
import pyworkflow.utils as pwutils
from pyworkflow.protocol import params, STATUS_NEW
from pyworkflow import PROD, BETA

from threading import Semaphore

# Define some descriptors
AND = 'by_all'
OR = 'by_at_least_one'
UNION_INTERSECTIONS = 'by_at_least_two'


class XmippProtConsensusPicking3D(ProtParticlePicking, XmippProtocol):
    """ 
        This protocol receives a set of coordinates representing 3D particles
        that come from different picking sources. It then trains a model and
        uses it to predict which ones are correct.

        The sources can be either different apps/algorithms or the same one
        with different parameters. The 3D coordinates are first examined and
        processed to ensure that differently-picked-same-particles are treated
        as such (given a tolerance).

        The model is trained with subsets formed by the picked particles (pos)
        and random noise in the volumes (neg).

        There is no streaming support yet.
    """
    # Identification parameters
    _label = 'deep consensus picking 3D'
    _devStatus = BETA
    _lastUpdateVersion = VERSION_3_0
    _conda_env = 'xmipp_DLTK_v0.3'
    _stepsCheckSecs = 5 # Scipion steps check interval (in seconds)

    # Protocol-specific options/switches/etc

    # Form options: NETWORK MODEL
    MODEL_TRAIN_NEW = 0
    MODEL_TRAIN_PRETRAIN = 1
    MODEL_TRAIN_PREVRUN = 2
    FORM_MODEL_TRAIN_TYPELIST_LABELS = ["From scratch", "Existing model", "Previous run"]
    FORM_MODEL_TRAIN_TYPELIST = [MODEL_TRAIN_NEW, MODEL_TRAIN_PRETRAIN, MODEL_TRAIN_PREVRUN]

    # Form options: additional data
    ADD_DATA_TRAIN_NEW = 0
    ADD_DATA_TRAIN_PRECOMP = 1
    ADD_DATA_TRAIN_CUST = 2
    ADD_DATA_TRAIN_TYPELIST_LABELS = ["None","Precompiled","Custom"]
    ADD_DATA_TRAIN_TYPELIST = []


    # Form options: work with ST or Coords
    FORM_DATA_TRAIN_CUSTOM_OPT = ["Subtomograms","Coordinates"]
    FORM_DATA_TRAIN_CUSTOM_OPT_SUBTOM = 0
    FORM_DATA_TRAIN_CUSTOM_OPT_COORDS = 1
    FORM_DATA_TRAIN_CUSTOM_TYPELIST = [FORM_DATA_TRAIN_CUSTOM_OPT_SUBTOM, FORM_DATA_TRAIN_CUSTOM_OPT_COORDS]

    

    # Parameters for the streaming part
    cnn_semaphore = Semaphore(1)

    def __init__(self, **args):
        ProtParticlePicking.__init__(self,**args)

    #--------------- DEFINE param functions ---------------

    def _defineParams(self, form : params.Form):
        ## Multiprocessing params
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                        label="Use GPU for the model (default: Yes)",
                        help="If yes, the protocol will try to use a GPU for "
                             "model training and execution. Note that this "
                             "will greatly decrease execution time."
                        )
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                        label="GPU ID",
                        help="Your system may have several GPUs installed, "
                             " choose the one you'd like to use(default: 0)."
                        )
        form.addParallelSection(threads=1, mpi=1)

        form.addSection(label='Input')

        group_model = form.addGroup('Neural Network model')

        ## Neural Network parameters
        group_model.addParam('modelInitialization', params.EnumParam,
            choices = self.FORM_MODEL_TRAIN_TYPELIST_LABELS,
            default = self.MODEL_TRAIN_NEW,
            label = 'Select a model',
            help = 'When set to *%s*, the network will start with a fresh and randomly '
            'initialized model. The option *%s* will let you choose a previously trained '
            'model. Lastly, *%s* will utilize the same model that was used in the '
            'previous run of this protocol within this project.'
            % tuple(self.FORM_MODEL_TRAIN_TYPELIST_LABELS))
        ## Model choices
        # For previous runs
        group_model.addParam('continueRun', params.PointerParam,
            pointerClass = self.getClassName(),
            condition = 'modelInitialization == %s'%self.MODEL_TRAIN_PREVRUN, allowsNull=True,
            label = 'Select previous run',
            help = 'Choose from a previous run to continue from.'
        )
        # For NOT NEW models
        group_model.addParam('skipTraining', params.BooleanParam,
            default = False,
            condition = 'modelInitialization != %s'%self.MODEL_TRAIN_NEW,
            label = 'Skip training step',
            help = ' When set to *Yes*, the volumes will be directly fed to the model, '
            ' If set to *No*, you must provide a training set of volumes.'
        )

        group_input = form.addGroup('Input')
        ## Input
        group_input.addParam('inputCoordinates', params.MultiPointerParam,
                        pointerClass = 'SetOfCoordinates', allowsNull=False,
                        label = 'Input coordinates',
                        help = 'Select the set of 3D coordinates that represent the subtomograms to be used as input data.'  
        )

        form.addSection(label='Preprocess')
        form.addParam('consensusRadius', params.FloatParam, default=0.1,
                        label="Same-element relative radius",
                        validators=[params.Positive],
                        help='Two sets of coordinates are determined to be of '
                        'the same particle if they are within this radius. The'
                        ' radius is given in [fraction of particle size] units.'
                        
        
        )
        group_input.addParam('classThreshold', params.FloatParam, default=0.5,
                        label = 'Tolerance threshold',
                        help='Choose a threshold in the range [0,1] to adjust '
                        'the threshold used internally to determine if the '
                        'input is classified as a _particle_ or as bad _noise_'
                        '. When set to -1 all particles are considered _good_.'
        )

        form.addSection(label='Training')
        form.addParam('nEpochs', params.IntParam, default=12,
                        label = 'Cycles (total epochs)',
                        help = 'Number of process cycles that will be done '
                        'with the data in order to train the Neural Network.',
        )
        form.addParam('learningrate', params.FloatParam, default = 0.01,
                        label = 'Learning rate',
                        help = 'Hyperparameter that controls the difference '
                        'between a calculated weight and its next value. '
                        'Higher values will result in faster learning, but '
                        'at the expense of sub-optimal weight values. '
                        'Very low value cause the NN to get stuck in local '
                        'minimas.'
        )
        form.addParam('dynLearningrate', params.BooleanParam, default = True,
                        label = 'Dynamic learning rate',
                        help = 'The learning rate can be updated on runtime '
                        'depending on the evolution of the execution. '                    
        )

        form.addParam('convergstop', params.BooleanParam, default = True,
                        label = 'Stop on convergence',
                        help = 'When set to *Yes*, the protocol will stop '
                        'the training when no improvement is detected in '
                        '2 consecutive epochs. Make sure that the learning '
                        'rate is not too low or this option might get your '
                        'model stopped in a sub-optimal, local minima. This'
                        ' is not recommended for small datasets.'
        )


        #form.addSection(label='Previously labeled data')
        form.addSection(label='Streaming')
        form.addParam('doPreliminarPredict', params.BooleanParam, default=False,
                        label = 'Predict before fully trained',
                        help = 'This protocol might do predictions using '
                        'the model before it is fully trained. These '
                        'results are stored in a different output set.'
        )

        form.addParam('extractingBatch', params.IntParam, default='5',
                        label = 'Extraction batch size',
                        help = 'Amount of subtomograms in an extraction batch.'
        )

        form.addParam('trainingBatch', params.IntParam, default='5',
                        label = 'Training batch size',
                        help = 'Amount of subtomograms in a training batch. '
                        'If the provided subtomograms are not enough for the '
                        'NN, this has to be increased.'
        )

    def _validate(self):
        error_message = []
        return error_message

    #--------------- INSERT steps functions ----------------
    

    def _insertAllSteps(self):
        pass

    #--------------- STEPS functions -----------------------

    

    #--------------- INFO functions -------------------------



    def _citations(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []

    #--------------- UTILS functions -------------------------
    def cnnFree(self):
        return self.cnn_semaphore._value == 1

    def loadTrainedParams(self):
        pass
    def saveTrainedParams(self, params):
        pass
        
class XmippProtDeepConsSubSet3D():
    """
        This protocol allows a user to create subsets from the GUI for the
        Deep Consensus 3D protocol to use.
    """
    def __init__(self, **args):
        raise NotImplementedError
