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

# Tomo-specific
from tomo.protocols import ProtTomoPicking
from tomo.objects import SetOfCoordinates3D, Coordinate3D, SetOfTomograms, Tomogram

# Needed for the GUI definition and pyworkflow objects
from pyworkflow.protocol import params
from pyworkflow import BETA
from pyworkflow.object import Integer, Float
import pyworkflow.utils as pwutils


import pandas as pd
import numpy as np

from pwem import emlib
from pwem.emlib.image import ImageHandler
import tomo.constants as tconst


# TODO: probably import the program class from Xmipp3.protocols
# Import the workers from the xmipp package

# Define some descriptors
AND = 'by_all'
OR = 'by_at_least_one'
UNION_INTERSECTIONS = 'by_at_least_two'

# How many pickers need to se something to choose it as positive picking
REQUIRED_PICKERS = 2

class XmippProtPickingConsensusTomo(ProtTomoPicking):
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
    #_conda_env = 'xmipp_DLTK_v0.3'
    _conda_env = 'tfm_mikel'
    _stepsCheckSecs = 5 # Scipion steps check interval (in seconds)
    _possibleOutputs = {'output3DCoordinates' : SetOfCoordinates3D}

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
        group_input.addParam('inputSets', params.MultiPointerParam,
                        pointerClass = SetOfCoordinates3D, allowsNull=False,
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

    #--------------- INSERT steps functions ----------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep(self.preProcessStep)
        self._insertFunctionStep(self.coordConsensusStep)
        self._insertFunctionStep(self.processTrainStep)
        self._insertFunctionStep(self.processScoreStep)
        self._insertFunctionStep(self.postProcessStep)
        self._insertFunctionStep(self.createOutputStep)

    #--------------- STEPS functions -----------------------

    # BLOCK 1 - Protocol - Load information
    def preProcessStep(self):
        """
        Block 1 - Aconditioning of data structures

        Generate the picker relation table, launch the coord
        consensus, generate the representative table, fill
        the representation table. Filter representatives and
        set labels for good, dubious and bad subtomos
        """

        # The form has a parameter called inputSets that carries
        # the 3D coordinates from the input pickers
        self.inputSetsOf3DCoordinates = [item.get() for item in self.inputSets]

        # Calculate the total amount of ROIs to save resources
        # The SetOfCoordinates3D is a EMSet -> Set
        totalROIs = sum(map(len, self.inputSetsOf3DCoordinates))

        # Ahora tengo: sets de coordenadas que vienen de distintos pickers
        # Toca: juntar todo

        # Combined table
        colnames = ['pick_id','x', 'y', 'z', 'tomo_id']
        self.untreated = pd.DataFrame(index=range(totalROIs),columns=colnames)

        # Pickers table
        colnames_md = ['boxsize', 'samplingrate']
        self.nr_pickers = len(self.inputSetsOf3DCoordinates)
        self.pickerMD = pd.DataFrame(index=range(self.nr_pickers), columns=colnames_md)

        # For each of the sets selected as input in the GUI...
        pickerCoordinates : SetOfCoordinates3D
        # Index for total ROIs
        indizea = 0
        for pick_id, pickerCoordinates in enumerate(self.inputSetsOf3DCoordinates):
            # Assign incrementing ID
            id = pick_id
            # Picker parameters
            bsize = int(pickerCoordinates.getBoxSize())
            srate = pickerCoordinates.getSamplingRate()
            # Assign the corresponding line
            self.pickerMD.loc[pick_id, 'boxsize'] = bsize
            self.pickerMD.loc[pick_id, 'samplingrate'] = srate
            # Get the coordinates
            coords = pickerCoordinates.iterCoordinates()

            # For each individual coordinate in this particular set...
            coordinate : Coordinate3D
            for coordinate in coords:
                asoc_vol : Tomogram = coordinate.getVolume()
                tomo_id = asoc_vol.getFileName()
                c_x = coordinate.getX(tconst.SCIPION)
                c_y = coordinate.getY(tconst.SCIPION)
                c_z = coordinate.getZ(tconst.SCIPION)
                self.untreated.loc[indizea, 'pick_id'] = id
                self.untreated.loc[indizea, 'x'] = c_x
                self.untreated.loc[indizea, 'y'] = c_y
                self.untreated.loc[indizea, 'z'] = c_z
                self.untreated.loc[indizea, 'tomo_id'] = tomo_id
                indizea += 1

        # Join the DFs to get all of the required information in one single DF
        self.untreated = self.untreated.join(self.pickerMD, on='pick_id')
        
        # Content of the self.untreated DF
        # pick_id','x', 'y', 'z', 'tomo_id', 'boxsize', 'samplingrate'

        # Get different tomogram names
        self.allTomoIds = self.untreated['tomo_id'].unique()
        self.coordinatesByTomogram = []
        self.coordinatesByTomogramFileNames = []

        # Generate a separate folder for each tomogram's coordinates
        self.FolderPickedPerTomo = self._getExtraPath() + "/pickedpertomo/"
        pwutils.makePath(self.FolderPickedPerTomo)

        # Generate per tomogram dataframes and write to XMD
        for name in self.allTomoIds:
            print("Unique tomogram found: " + name)
            singleTomoDf : pd.DataFrame = self.untreated[self.untreated['tomo_id'] == name]
            self.coordinatesByTomogram.append(singleTomoDf)
            savedfile = self.writeCoords(self.FolderPickedPerTomo, singleTomoDf, name)
            self.coordinatesByTomogramFileNames.append(savedfile)

        # Print sizes before doing the consensus
        print(str(self.nr_pickers) + " pickers with a total of "+ str(totalROIs)+ " coordinates from "
              + str(len(self.allTomoIds)) + " individual tomograms.")
        print("\nPICKER SUMMARY")
        print(self.pickerMD)
        print("")
    
        # Do the box size consensus
        self.boxSizeConsensusStep()
        self.samplingRateConsensusStep()

        # Ahora tengo: un fichero por cada tomogram_id con todos sus pickings
        # End block
        # END STEP

    # BLOCK 1 - Protocol - write coords from DF (raw, not consensuated)
    def writeCoords(self, path, df, tomoname) -> str:
        """
        Block 1 AUX - Write coordinates into Xmipp Metadata format
        path: folder to save the data
        df: dataframe containing the picking data
        tomoname: tomogram of which this data is from
        """

        # (String) Path of the tomogram volume, get only the filename (last element of split-array)
        # Also remove .xmd
        filename = str(tomoname).split("/")[-1].split(".")[0]
        print("Saving coords for... " + filename )
        
        # Create a Xmipp MD Object
        outMD = emlib.MetaData()
        outMD.setColumnValues(emlib.MDL_REF, df['pick_id'].tolist())
        outMD.setColumnValues(emlib.MDL_X, df['x'].tolist())
        outMD.setColumnValues(emlib.MDL_Y, df['y'].tolist())
        outMD.setColumnValues(emlib.MDL_Z, df['z'].tolist())
        outMD.setColumnValues(emlib.MDL_PICKING_PARTICLE_SIZE, df['boxsize'].tolist())
        outMD.setColumnValues(emlib.MDL_SAMPLINGRATE, df['samplingrate'].tolist())
        outMD.setColumnValues(emlib.MDL_TOMOGRAM_VOLUME, df['tomo_id'].tolist())
        
        composedFileName = path + filename + "_allpickedcoords.xmd"
        outMD.write(composedFileName)
        return composedFileName
    
    # BLOCK 1 - Protocol - select box size
    def boxSizeConsensusStep(self, method="biggest"):
        """
        Block 1 AUX - Perform consensus in the box size
        
        Consensuates the box size from the different inputs according to
        a selected method.

        Methods:
          - biggest: max value amongst the pickers
          - smallest: min value amongst the pickers
          - mean: average value amongst the pickers
          - first: first in the list
        """

        # Fetch the different box sizes from pickers
        assert self.pickerMD is not None
        values = self.pickerMD['boxsize']

        result = self.valueConsensus(values, method)

        print("Determined box size: " + str(result))
        
        self.consBoxSize = Integer(result)

   # BLOCK 1 - Protocol - choose sampling rate consensus
    def samplingRateConsensusStep(self, method="smallest"):
        """
        Block 1 AUX - Perform consensus in the sampling rate (angstroms/px)
        
        Consensuates the sampling rate from the different inputs according to
        a selected method.

        Methods:
          - biggest: max value amongst the pickers
          - smallest: min value amongst the pickers
          - mean: average value amongst the pickers
          - first: first in the list
        """

        # Fetch the different box sizes from pickers
        assert self.pickerMD is not None
        values = self.pickerMD['samplingrate']
        result = self.valueConsensus(values, method)
        print("Determined sampling rate (A/px): " + str(result))
        self.consSampRate = Float()
    
    # BLOCK 1 - Protocol - Auxiliar value consensus function
    def valueConsensus(self, inputSet: pd.Series, method="biggest"):
        values = inputSet

        if method == "biggest":
            result = values.max()
        elif method == "smallest":
            result = values.min()
        elif method == "mean":
            result = values.mean()
        elif method == "first":
            result = values[0]

        return result

    # BLOCK 2 - Program - consensuate coordinates
    def coordConsensusStep(self):
        """
        Block 2 - Perform consensus in the coordinates

        This step launches a call to the associated Xmipp program, triggering
        a 3D coordinates consensus.
        """

        # Generate a separate folder for the coord consensus output
        self.FolderCoordConsensus = self._getExtraPath() + "/coordConsensus/"
        pwutils.makePath(self.FolderCoordConsensus)

        program = "xmipp_coordinates_consensus_tomo"
        tomoPickingMdFname : str
        ih = ImageHandler()
        for tomoPickingMdFname in self.coordinatesByTomogramFileNames:
            dims = ih.getDimensions(tomoPickingMdFname)
            args = ''
            args += '--input ' + tomoPickingMdFname
            args += ' --outputAll ' + self.FolderCoordConsensus + tomoPickingMdFname.split("/")[-1].split(".")[0] + "_consensus_all.xmd"
            args += ' --outputPos ' + self.FolderCoordConsensus + tomoPickingMdFname.split("/")[-1].split(".")[0] + "_consensus_pos.xmd"
            args += ' --outputDoubt ' + self.FolderCoordConsensus + tomoPickingMdFname.split("/")[-1].split(".")[0] + "_consensus_doubt.xmd"
            args += ' --boxsize ' + str(self.consBoxSize)
            args += ' --radius ' + str(float(self.consensusRadius.get()))
            args += ' --number ' + str(REQUIRED_PICKERS)
            args += ' --size ' + ' '.join(dims)
            print('\nHanding over to Xmipp program for coordinate consensus')
            self.runJob(program, args)
                     
    # BLOCK 2 - Program - Launch NN train (if needed)
    def processTrainStep(self):
        """
        Block 2 - Neural Network TRAIN operations

        Extract the good subtomos,
        generate noise for data augmentation, do the data splits
        and launch the neural network training. Save the NN
        model into a H5 file both intermediate and final.
        """
        # This protocol executes on the Xmipp program side

        # TODO: Integrate with picking noise program
        pass

    # BLOCK 2 - Program - Load model (if needed) and score 
    def processScoreStep(self):
        """
        Block 2 - Neural Network SCORE operations

        Load the selected NN model, extract the dubious
        subtomograms and score them against the model. Then
        generate the score tables.
        """
        # This protocol executes on the Xmipp program side
        pass

    # BLOCK 3 - filter tables according to thresholds
    def postProcessStep(self):
        """
        Block 3 - Handling of posterous things

        Apply the desired filters and thresholds to the 
        results.
        """
        pass
    
    # BLOCK 3 - Prepare output for Scipion GUI
    def createOutputStep(self):
        """
        Block 3 - Prepare output

        Leave everything in a format suitable for later processing
        in the CryoET workflows.
        """
        #
        tomograms : SetOfTomograms = None
        outputPath = self._getExtraPath("CBOX_3D")
        suffix = self._getOutputSuffix(SetOfCoordinates3D)

        coordinates :SetOfCoordinates3D = self._createSetOfCoordinates3D()
        coordinates.setName("")
        coordinates.setSamplingRate()


        # TODO: Xmipp metadata para salida (no guai para next protocols)
        # A ve tambien te digo mejor parametro extendido
        name = self.OUTPUT_PREFIX + suffix
        self._defineOutputs(**{name: coordinates})
        self._defineSourceRelation(tomograms, coordinates)

    #--------------- INFO functions -------------------------

    #--------------- UTILS functions -------------------------

    def _validate(self):
        errors = []
        errors += self._validateParallelProcessing()
        errors += self._validateCondaEnvironment()
        errors += self._validateXmippBinaries()
        return errors

    def _validateParallelProcessing(self):
        nGpus = len(self.getGpuList())
        nMPI = self.numberOfMpi.get()
        nThreads = self.numberOfThreads.get()
        errors = []
        if nGpus < 1:
            errors.append("A GPU is needed for this protocol to run.")
        if nThreads != 1:
            errors.append("Multithreading not yet supported. Set Threads parameter to 1.")
        if nMPI != 1:
            errors.append("Multiprocessing not yet supported. Set MPI parameter to 1")
        return errors

    def _validateCondaEnvironment(self):
        errors = []
        return errors

    def _validateXmippBinaries(self):
        errors = []
        return errors
        
class XmippProtDeepConsSubSet3D():
    """
        This protocol allows a user to create subsets from the GUI for the
        Deep Consensus 3D protocol to use.
    """
    def __init__(self, **args):
        raise NotImplementedError
