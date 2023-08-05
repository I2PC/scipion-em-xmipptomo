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
# * Initial release: september 2023
# **************************************************************************

"""
Deep Consensus picking protocol suited for Cryo-ET
"""
import os, time

# Tomo-specific
from tomo.protocols import ProtTomoPicking
from tomo.objects import SetOfCoordinates3D, Coordinate3D, SetOfTomograms, Tomogram

# Needed for the GUI definition and pyworkflow objects
from pyworkflow.protocol import params, LEVEL_ADVANCED
from pyworkflow import BETA
from pyworkflow.object import Integer, Float
import pyworkflow.utils as pwutils


import pandas as pd
import numpy as np

from pwem import emlib
from pwem.emlib.image import ImageHandler
import tomo.constants as tconst

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
    _label      = 'deep consensus picking 3D'
    _devStatus  = BETA
    #_conda_env = 'xmipp_DLTK_v0.3'
    _conda_env  = 'scipion3'
    _stepsCheckSecs  = 5 # Scipion steps check interval (in seconds)
    _possibleOutputs = {'output3DCoordinates' : SetOfCoordinates3D}

    # Protocol-specific options/switches/etc

    # Form options: NETWORK MODEL
    MODEL_TRAIN_NEW         = 0
    MODEL_TRAIN_PRETRAIN    = 1
    MODEL_TRAIN_PREVRUN     = 2
    FORM_MODEL_TRAIN_TYPELIST_LABELS = ["From scratch", "Existing model", "Previous run"]
    FORM_MODEL_TRAIN_TYPELIST        = [MODEL_TRAIN_NEW, MODEL_TRAIN_PRETRAIN]

    # Form options: additional data
    ADD_DATA_TRAIN_NEW      = 0
    ADD_DATA_TRAIN_PRECOMP  = 1
    ADD_DATA_TRAIN_CUST     = 2
    ADD_DATA_TRAIN_TYPELIST_LABELS  = ["None", "Precompiled", "Custom"]
    ADD_DATA_TRAIN_TYPELIST         = [ADD_DATA_TRAIN_NEW, ADD_DATA_TRAIN_PRECOMP, ADD_DATA_TRAIN_CUST]

    # Form options: Coordinates consensus representant election
    COORD_CONS_FIRST    = 0
    COORD_CONS_CENTROID = 1
    FORM_COORD_CONS_TYPELIST_LABELS = ["First found", "Calculate centroid"]
    FORM_COORD_CONS_TYPELIST        = [COORD_CONS_FIRST, COORD_CONS_CENTROID]

    # Form options: Values consensus
    VALUE_CONS_FIRST    = 0
    VALUE_CONS_BIG      = 1
    VALUE_CONS_SMALL    = 2
    VALUE_CONS_MEAN     = 3
    FORM_VALUE_CONS_TYPELIST_LABELS = ["First found", "Biggest", "Smallest", "Mean value"]
    FORM_VALUE_CONS_TYPELIST        = [VALUE_CONS_FIRST, VALUE_CONS_BIG, VALUE_CONS_SMALL, VALUE_CONS_MEAN]

    #--------------- DEFINE param functions ---------------

    def _defineParams(self, form : params.Form):
        ## Multiprocessing params
        # form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
        #                 label="Use GPU for the model (default: Yes)",
        #                 help="If yes, the protocol will try to use a GPU for "
        #                      "model training and execution. Note that this "
        #                      "will greatly decrease execution time."
        #                 )
        # form.addHidden(params.GPU_LIST, params.StringParam, default='0',
        #                 label="GPU ID",
        #                 help="Your system may have several GPUs installed, "
        #                      " choose the one you'd like to use(default: 0)."
        #                 )
        # form.addParallelSection(threads=1, mpi=1)

        form.addSection(label='Main')

        group_model = form.addGroup('Neural Network model')

        group_model.addParam('votingMode', params.BooleanParam,
            default = False,
            expertLevel=LEVEL_ADVANCED,
            label = 'Voting instead of DNN',
            help = 'Activating this will ignore all NN parameters and only perform '
            'a consensus based on regular voting after doing the coordinates fusion.'
        )

        group_model.addParam('votingThreshold', params.FloatParam,
            default = 0.5,
            label = "Required consensus threshold",
            condition = 'votingMode = True',
            help = 'Sets the required consensus threshold (0,1] ratio needed for '
            'a result to be considered good in simple voting mode. Bear in mind '
            'the amount of input pickers when choosing this value. For instance, '
            'a 0.7 value will not be achievable with two input pickers as the volumes '
            'will appear either in 0.5 or in 1.0 of the total pickers.'
        )

        ## Neural Network parameters
        group_model.addParam('modelInitialization', params.EnumParam,
            choices = self.FORM_MODEL_TRAIN_TYPELIST_LABELS,
            default = self.MODEL_TRAIN_NEW,
            label = 'Select a model',
            # help = 'When set to *%s*, the network will start with a fresh and randomly '
            # 'initialized model. The option *%s* will let you choose a previously trained '
            # 'model. Lastly, *%s* will utilize the same model that was used in the '
            # 'previous run of this protocol.'
            help = 'When set to *%s*, the network will start with a fresh and randomly '
            'initialized model. The option *%s* will let you choose a previously trained '
            'model.'
            % tuple(self.FORM_MODEL_TRAIN_TYPELIST_LABELS))
        ## Model choices
        # For previous runs
        # group_model.addParam('continueRun', params.PointerParam,
        #     pointerClass = self.getClassName(),
        #     condition = 'modelInitialization == %s'%self.MODEL_TRAIN_PREVRUN, allowsNull=True,
        #     label = 'Select previous run',
        #     help = 'Choose from a previous run to continue from.'
        # )
        # For NOT NEW models
        group_model.addParam('skipTraining', params.BooleanParam,
            default = False,
            condition = 'modelInitialization != %s'%self.MODEL_TRAIN_NEW,
            label = 'Skip training step',
            help = ' When set to *Yes*, the volumes will be directly fed to the model, '
            ' If set to *No*, you must provide a training set of volumes.'
        )
        form.addParam('trainingBatch', params.IntParam, default='5',
                        label = 'Training batch size',
                        help = 'Amount of subtomograms in a training batch. '
                        'If the provided subtomograms are not enough for the '
                        'NN, this has to be increased. If the machine hangs due'
                        ' to memory issues, this has to be reduced.'
        )

        group_input = form.addGroup('Input')
        ## Input
        group_input.addParam('inputSets', params.MultiPointerParam,
                        pointerClass = SetOfCoordinates3D, allowsNull=False,
                        label = 'Input coordinates',
                        help = 'Select the set of 3D coordinates that represent the subtomograms to be used as input data.'  
        )
        group_input.addParam('classThreshold', params.FloatParam, default=0.5,
                        label = 'Tolerance threshold',
                        help='Choose a threshold in the range (0,1] to adjust '
                        'the threshold used internally to determine if the '
                        'input is classified as a _particle_ or as bad _noise_'
                        '. When set to -1 all particles are considered _good_.'
        )

        form.addSection(label='Preprocess')
        form.addParam('coordConsensusRadius', params.FloatParam, default=0.1,
                        label="Same-element relative radius",
                        validators=[params.Positive],
                        help='Two sets of coordinates are determined to be of '
                        'the same particle if they are within this radius. The'
                        ' radius is given in [fraction of particle size] units.'
        )

        form.addParam('coordConsensusType', params.EnumParam, 
                      choices = self.FORM_COORD_CONS_TYPELIST_LABELS,
                      default = self.COORD_CONS_FIRST,
                      label = 'Representant choosing method',
                      help = 'When assimilating all the pickings related to the'
                      ' same ROI... *%s* will choose the number of the first '
                      'element in the list, while *%s* will calculate a mean and '
                      'later force the resize of all subtomograms to match it.' 
                      % tuple(self.FORM_COORD_CONS_TYPELIST_LABELS)
        )

        form.addParam('valueConsensusType', params.EnumParam,
                      
                      choices = self.FORM_VALUE_CONS_TYPELIST_LABELS,
                      default = self.VALUE_CONS_SMALL,
                      label = 'Boxsize choosing method',
                      help = '*%s* *%s* *%s* *%s*',
        )

        form.addParam('fracNoise', params.FloatParam, default=0.9,
                      label="Amount of noise picked for negative input",
                      help='Controls how much noise is picked and given '
                      'to the NN as negative input during training. It is'
                      ' expressed in [0..1] - fraction of the total amount'
                      ' of coordinates found on input')

        

        form.addSection(label='Training')
        form.addParam('nEpochs', params.IntParam, default=6,
                        label = 'Cycles (total epochs)',
                        help = 'Number of process cycles that will be done '
                        'with the data in order to train the Neural Network.',
        )

        form.addParam('learningRate', params.FloatParam, default = 0.01,
                        label = 'Learning rate',
                        help = 'Hyperparameter that controls the difference '
                        'between a calculated weight and its next value. '
                        'Higher values will result in faster learning, but '
                        'at the expense of sub-optimal weight values. '
                        'Very low value cause the NN to get stuck in local '
                        'minimas.'
        )

        form.addParam('dynLearningRate', params.BooleanParam, default = True,
                        label = 'Dynamic learning rate',
                        help = 'The learning rate can be updated on runtime '
                        'depending on the evolution of the execution. '                    
        )

        form.addParam('convergStop', params.BooleanParam, default = True,
                        label = 'Stop on convergence',
                        help = 'When set to *Yes*, the protocol will stop '
                        'the training when no improvement is detected in '
                        '2 consecutive epochs. Make sure that the learning '
                        'rate is not too low or this option might get your '
                        'model stopped in a sub-optimal, local minima. This'
                        ' is not recommended for small datasets.'
        )

        form.addParam('forceDataAugment', params.BooleanParam, default = False,
                      label = "Force data augmentation",
                      help = 'By default, the protocol will not try to '
                      'perform data augmentation on datasets if at least two '
                      'of the input pickers contain 900 structures detected.'
                      )


        # form.addSection(label='Streaming')
        # form.addParam('doPreliminarPredict', params.BooleanParam, default=False,
        #                 label = 'Predict before fully trained',
        #                 help = 'This protocol might do predictions using '
        #                 'the model before it is fully trained. These '
        #                 'results are stored in a different output set.'
        # )

        # form.addParam('extractingBatch', params.IntParam, default='5',
        #                 label = 'Extraction batch size',
        #                 help = 'Amount of subtomograms in an extraction batch.'
        # )

        

    #--------------- INSERT steps functions ----------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep(self.preProcessStep)
        self._insertFunctionStep(self.coordConsensusStep)
        self._insertFunctionStep(self.prepareNNStep)
        if not bool(self.skipTraining):
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

        # GET THE INFORMATION FROM THE FORM AND STORE IN SELF -----------------
        # The form has a parameter called inputSets that carries
        # the 3D coordinates from the input pickers
        self.inputSetsOf3DCoordinates = [item.get() for item in self.inputSets]

        # Calculate the total amount of ROIs to save resources
        # The SetOfCoordinates3D is a EMSet -> Set so len() can be applied
        self.totalROIs : int = sum(map(len, self.inputSetsOf3DCoordinates))

        # Get the number of input pickers given in the form
        self.nr_pickers : int = len(self.inputSetsOf3DCoordinates)

        # Get fraction for noise picking
        self.noiseFrac = float(self.fracNoise)

        # Get the method of doing the value consensus
        self.valueConsType : int = self.valueConsensusType

        # Get the method of coordinate consensus
        self.coordConsType : int = self.coordConsensusType

        # GENERATE THE NEEDED TABLES TO START ---------------------------------
        # Combined table of untreated data
        colnames = ['pick_id','x', 'y', 'z', 'tomo_id', 'boxsize', 'samplingrate']
        self.untreated = pd.DataFrame(index=range(self.totalROIs),columns=colnames)

        # Pickers data table
        colnames_md = ['boxsize', 'samplingrate']
        self.pickerMD = pd.DataFrame(index=range(self.nr_pickers), columns=colnames_md)

        # START ENTERING THE DATA INTO THE RAW DATA TABLE ---------------------
        # Index for total ROIs inside the next loop
        globalIndex = 0
        # For each of the sets selected as input in the GUI...
        pickerCoordinates : SetOfCoordinates3D
        for pick_id, pickerCoordinates in enumerate(self.inputSetsOf3DCoordinates):
            # Picker parameters
            bsize = int(pickerCoordinates.getBoxSize())
            srate = pickerCoordinates.getSamplingRate()

            # Assign the corresponding line
            self.pickerMD.loc[pick_id, 'boxsize'] = bsize
            self.pickerMD.loc[pick_id, 'samplingrate'] = srate

            # For each individual coordinate in this particular set...
            coordinate : Coordinate3D
            for coordinate in pickerCoordinates.iterCoordinates():
                asoc_vol : Tomogram = coordinate.getVolume()
                tomo_id = asoc_vol.getFileName()
                c_x = coordinate.getX(tconst.BOTTOM_LEFT_CORNER)
                c_y = coordinate.getY(tconst.BOTTOM_LEFT_CORNER)
                c_z = coordinate.getZ(tconst.BOTTOM_LEFT_CORNER)
                self.untreated.loc[globalIndex, 'pick_id'] = pick_id
                self.untreated.loc[globalIndex, 'x'] = c_x
                self.untreated.loc[globalIndex, 'y'] = c_y
                self.untreated.loc[globalIndex, 'z'] = c_z
                self.untreated.loc[globalIndex, 'tomo_id'] = tomo_id
                self.untreated.loc[globalIndex, 'boxsize'] = bsize
                self.untreated.loc[globalIndex, 'samplingrate'] = srate
                globalIndex += 1
        
        # Content of the self.untreated DF now
        # pick_id','x', 'y', 'z', 'tomo_id', 'boxsize', 'samplingrate'

        # Get different tomogram names
        self.uniqueTomoIDs = self.untreated['tomo_id'].unique()
        # self.coordinatesByTomogramFileNames = []

        # Generate a separate folder for each tomogram's coordinates
        pwutils.makePath(self._getPickedPerTomoPath())

        # Generate per tomogram dataframes and write to XMD
        for name in self.uniqueTomoIDs:
            print("Unique tomogram found: " + name)
            singleTomoDf : pd.DataFrame = self.untreated[self.untreated['tomo_id'] == name]
            savedfile = self._getAllCoordsFilename(self._stripTomoFilename(name))
            self.writeCoords(savedfile, singleTomoDf)
            # self.coordinatesByTomogramFileNames.append(savedfile)

        # Print sizes before doing the consensus
        print(str(self.nr_pickers) + " pickers with a total of "+ str(self.totalROIs)+ " coordinates from "
              + str(len(self.uniqueTomoIDs)) + " unique tomograms.")
        print("\nPICKER SUMMARY")
        print(self.pickerMD)
        print("")
    
        # Do the box size consensus
        self.BSSRConsensusStep(self.valueConsType)

        # Ahora tengo: un fichero por cada tomogram_id con todos sus pickings, asi como consenso en tamannos
        # End block
        # END STEP

    # BLOCK 1 - Protocol - write coords from DF (raw, not consensuated)
    def writeCoords(self, outpath: str, df: pd.DataFrame):
        """
        Block 1 AUX - Write coordinates into Xmipp Metadata format
        path: folder to save the data
        df: dataframe containing the picking data
        tomoname: tomogram of which this data is from
        """

        # (String) Path of the tomogram volume, get only the filename (last element of split-array)
        # Also remove .xmd
        print("Saving coords for... " + outpath )
        
        # Create a Xmipp MD Object
        outMD = emlib.MetaData()
        outMD.setColumnValues(emlib.MDL_REF, df['pick_id'].tolist())
        outMD.setColumnValues(emlib.MDL_XCOOR, list(map(int, df['x'])))
        outMD.setColumnValues(emlib.MDL_YCOOR, list(map(int, df['y'])))
        outMD.setColumnValues(emlib.MDL_ZCOOR, list(map(int, df['z'])))
        outMD.setColumnValues(emlib.MDL_PICKING_PARTICLE_SIZE, df['boxsize'].tolist())
        outMD.setColumnValues(emlib.MDL_SAMPLINGRATE, df['samplingrate'].tolist())
        outMD.setColumnValues(emlib.MDL_TOMOGRAM_VOLUME, df['tomo_id'].tolist())
        
        outMD.write(outpath)
    
    # BLOCK 1 - Protocol - select box size
    def BSSRConsensusStep(self, method):
        """
        Block 1 AUX - Perform consensus in the box size and samplingrate
        
        Consensuates the BSSR from the different inputs according to
        a selected method.

        Methods:
          - biggest: max value amongst the pickers
          - smallest: min value amongst the pickers
          - mean: average value amongst the pickers
          - first: first in the list
        """

        # Fetch the different box sizes from pickers
        assert self.pickerMD is not None

        if method ==  self.VALUE_CONS_BIG:
            result = self.pickerMD.iloc[self.pickerMD['boxsize'].argmax()]
        elif method == self.VALUE_CONS_SMALL:
            result = self.pickerMD.iloc[self.pickerMD['boxsize'].argmin()]
        elif method == self.VALUE_CONS_MEAN:
            raise NotImplemented
        elif method == self.VALUE_CONS_FIRST:
            result = self.pickerMD[0]

        print("Determined box size: " + str(result['boxsize']))
        print("Determined sampling rate (A/px): " + str(result['samplingrate']))
        
        self.consBoxSize = Integer(result['boxsize'])
        self.consSampRate = Float(result['samplingrate'])  

    # BLOCK 2 - Program - consensuate coordinates
    def coordConsensusStep(self):
        """
        Block 2 - Perform consensus in the coordinates

        This step launches a call to the associated Xmipp program, triggering
        a 3D coordinates consensus.
        """

        # Generate a separate folder for the coord consensus output
        pwutils.makePath(self._getCoordConsensusPath())

        program = "xmipp_coordinates_consensus_tomo"
        for tomo_id in self.uniqueTomoIDs:
            tomoname = self._stripTomoFilename(tomo_id)
            args = ''
            args += '--input ' + self._getAllCoordsFilename(tomoname)
            args += ' --outputAll ' + self._getConsCoordsFilename(tomoname)
            args += ' --outputPos ' + self._getPosCoordsFilename(tomoname)
            args += ' --outputDoubt ' + self._getDoubtCoordsFilename(tomoname)
            args += ' --boxsize ' + str(self.consBoxSize)
            args += ' --radius ' + str(float(self.coordConsensusRadius.get()))
            args += ' --number ' + str(self.nr_pickers)
            print('\nHanding over to Xmipp program for coordinate consensus')
            self.runJob(program, args)
    
    # BLOCK 2 - Program - Launch Noise Picking algorithm for data
    def noisePick(self, tomoPath, coordsPath, outPath):
        """
        Block 2 - Noise picker for training stages

        Generates noise for input tomogram, thus adding an
        almost certain negative input for the NN training stage,
        when training is selected from the GUI or the model is 
        built from scratch.
        """

        program = "xmipp_pick_noise_tomo"
        ih = ImageHandler()
        
        # Prepare and launch script
        dims = ih.getDimensions(tomoPath)
        args = ''
        args += ' --input ' + coordsPath # Ruta al _allpickedcoordinates
        args += ' --radius 0.8'
        args += ' --boxsize ' + str(self.consBoxSize)
        args += ' --samplingrate ' + str(self.consSampRate)
        args += ' --size ' + ' '.join(map(str, dims[0:3]))
        args += ' --limit ' + str(self.noiseFrac)
        args += ' --threads ' + str(self.numberOfThreads)
        args += ' --output ' + outPath
        print('\nHanding over to Xmipp program for noise picking')
        self.runJob(program, args)                   
  
    def prepareNNStep(self):
        """
        Block 2 - Neural Network TRAIN operations

        - Launch extraction of all needed subtomos into /dataset/pos | neg | doubt
        - Launch split generation
        - Launch NN tren (chuchú)
        - Save the NN
        model into a H5 file both intermediate and final.
        """
        # This protocol executes on the Xmipp program side

        # Check if train is needed
        trainSkip = bool(self.skipTraining)

        # Generate directory structure
        folders = [self._getDoubtSubtomogramPath()]

        # Only generate if training is needed
        if not trainSkip:
            folders.append(self._getPosSubtomogramPath())
            folders.append(self._getNegSubtomogramPath())

        for folder in folders:
            pwutils.makePath(folder)

        # Tengo: carpetas
        # Necesito: extraer
        for tomoPath in self.uniqueTomoIDs:
            tomoName = self._stripTomoFilename(tomoPath)
            if not trainSkip:
                self.tomogramExtract(tomoPath, self._getPosCoordsFilename(tomoName), self._getPosSubtomogramPath())
                self.noisePick(tomoPath, self._getAllCoordsFilename(tomoName), self._getNegCoordsFilename(tomoName))
                self.tomogramExtract(tomoPath, self._getNegCoordsFilename(tomoName), self._getNegSubtomogramPath())
            self.tomogramExtract(tomoPath, self._getDoubtCoordsFilename(tomoName), self._getDoubtSubtomogramPath())
        # Tengo: todo extraido
        # Fin de paso - next: train if needed   
            
    # BLOCK 2 - Program - Launch NN train (if needed)
    def processTrainStep(self):
        """
        Block 2 - Call for NN train

        Launches the script responsible for training the NN.
        """
        folder = self._getNnPath()
        pwutils.makePath(folder)

        program = "xmipp_deep_picking_consensus_tomo"
        args = ''
        args += ' -t ' + self.numberOfThreads
        args += ' -g ' + ','.join(map(str, self.getGpuList()))
        args += ' --mode training'
        args += ' --netpath ' + self._getNnPath()
        args += ' --consboxsize ' + str(self.consBoxSize)
        args += ' --conssamprate ' + str(self.consSampRate)
        args += ' --truevolpath ' + self._getPosSubtomogramPath()
        args += ' --falsevolpath ' + self._getNegSubtomogramPath()
        args += ' -e 5'
        args += ' -l 0.001'
        args += ' -r 0.0001'
        args += ' --ensemble 1'
        self.runJob(program, args)

        print('\nHanding over to Xmipp program for Train and Score')
             
          
    # BLOCK 2 - Program - Launch tomogram extraction step
    def tomogramExtract(self, tomoPath, coordsPath, outPath):
        """
        Block 2 - NN preparation

        This program extracts subtomograms using Xmipp.
        tomoPath - .mrc file path
        coordsPath - .xmd coordinates
        outPath - folder to leave .mrc files in
        """

        program = "xmipp_tomo_extract_subtomograms"
        args = ''
        args += '--tomogram ' + tomoPath
        args += ' --coordinates ' + coordsPath
        args += ' --boxsize ' + str(self.consBoxSize)
        args += ' -o ' + outPath
        args += ' --threads ' + str(self.numberOfThreads)

        print('\nHanding over to Xmipp program for tomogram extraction')
        self.runJob(program, args)


    # BLOCK 2 - Program - Load model (if needed) and score 
    def processScoreStep(self):
        """
        Block 2 - Neural Network SCORE operations

        Load the selected NN model, extract the dubious
        subtomograms and score them against the model. Then
        generate the score tables.
        """
        # This protocol executes on the Xmipp program side
        print('\nHanding over to Xmipp program for Train and Score')

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
        errors += self._validateNumberOfPickers()
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

    def _validateNumberOfPickers(self):
        nPickers = len(self.inputSets)
        errors = []
        if nPickers < 2:
            errors.append("More than one picker must be provided to do any kind of consensus.")
        return errors
    
        #--------------- FILENAMES functions -------------------

    def _getNnPath(self, *args):
        return self._getExtraPath('nn')

    def _stripTomoFilename(self, tomopath: str):
        return tomopath.split("/")[-1].split(".")[0]
    
    def _getCoordConsensusPath(self, *args):
        return self._getExtraPath('coordconsensus', *args)

    def _getPosCoordsFilename(self, tomo_name: str):
        return self._getCoordConsensusPath(tomo_name+"_cons_pos.xmd")
    
    def _getNegCoordsFilename(self, tomo_name: str):
        return self._getCoordConsensusPath(tomo_name+"_cons_neg.xmd")
    
    def _getDoubtCoordsFilename(self, tomo_name: str):
        return self._getCoordConsensusPath(tomo_name+"_cons_doubt.xmd")
    
    def _getConsCoordsFilename(self, tomo_name:str):
        return self._getCoordConsensusPath(tomo_name+"_cons_all.xmd")
    
    def _getAllCoordsFilename(self, tomo_name: str):
        return self._getPickedPerTomoPath(tomo_name+"_allpickedcoords.xmd")

    def _getDatasetPath(self, *args):
        return self._getExtraPath('dataset', *args)
    
    def _getPosSubtomogramPath(self, *args):
        return self._getDatasetPath('pos', *args)
    
    def _getNegSubtomogramPath(self, *args):
        return self._getDatasetPath('neg', *args)
    
    def _getDoubtSubtomogramPath(self, *args):
        return self._getDatasetPath('doubt', *args)

    def _getPickedPerTomoPath(self, *args):
        return self._getExtraPath('pickedpertomo', *args)
    
    
        
class XmippProtDeepConsSubSet3D():
    """
        This protocol allows a user to create subsets from the GUI for the
        Deep Consensus 3D protocol to use.
    """
    def __init__(self, **args):
        raise NotImplementedError
