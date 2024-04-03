# **************************************************************************
# *
# * Authors:    Mikel Iceta Tena (miceta@cnb.csic.es)
# *             Jose Luis Vilas (jlvilas@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# * v 0.2 - January 2024
# *     * Complete refactor, now Pandas is not needed
# *     * Support for irregular sets (different nº of tomos, different Srate, different BSize)
# **************************************************************************

"""
Deep Consensus picking protocol suited for Cryo-ET
"""
import os,time

# Tomo-specific
import tomo.constants as tconst
from tomo.protocols import ProtTomoPicking
from xmipp3 import XmippProtocol
from tomo.objects import SetOfCoordinates3D, Coordinate3D, SetOfTomograms, Tomogram

# Needed for the GUI definition and pyworkflow objects
from pyworkflow.protocol import params, LEVEL_ADVANCED, STEPS_PARALLEL
from pyworkflow import BETA
import pyworkflow.utils as pwutils

# Data and objects management
import numpy as np
from pyworkflow.object import Integer, Float, Boolean, List, Dict, String
from pwem.protocols import EMProtocol
from pwem import emlib
from pwem.emlib.image import ImageHandler

class XmippProtPickingConsensusTomo(ProtTomoPicking, EMProtocol, XmippProtocol):
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
    _conda_env = 'xmipp_DLTK_v1.0'
    _stepsCheckSecs  = 60 # Scipion steps check interval (in seconds)
    _possibleOutputs = {'output3DCoordinates' : SetOfCoordinates3D}

    # Protocol-specific options/switches/etc

    # Form options: NETWORK MODEL
    MODEL_TRAIN_NEW         = 0
    MODEL_TRAIN_PRETRAIN    = 1
    FORM_MODEL_TRAIN_TYPELIST_LABELS = ["From scratch", "Existing model"]
    FORM_MODEL_TRAIN_TYPELIST        = [MODEL_TRAIN_NEW, MODEL_TRAIN_PRETRAIN]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------- DEFINE param functions ---------------
    def _defineParams(self, form : params.Form):
        # Multiprocessing params
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                        label="Use GPU for the model (default: Yes)",
                        help="If yes, the protocol will try to use a GPU for "
                             "model training and execution. Note that this "
                             "will greatly decrease execution time."
                        )
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                        label="GPU ID (default: 0)",
                        help="Your system may have several GPUs installed, "
                             " choose the one you'd like to use."
                        )
        form.addParallelSection(threads=8, mpi=1)
        # Main parameters
        form.addSection(label='Main')
        form.addParam('inputSets', params.MultiPointerParam,
                        pointerClass = SetOfCoordinates3D, allowsNull=False,
                        label = 'Input coordinates',
                        help = 'Select the set of 3D coordinates that represent the subtomograms to be used as input data.'  
        )
        form.addParam('neededNumberOfPickers', params.IntParam, default=2,
                        label="Input positive consensus",
                        help='Amount of input pickers choosing a coordinate needed '
                        'to deem a coordinate as a positive input for the NN.'
        )
        form.addParam('coordConsensusRadius', params.FloatParam, default=0.5,
                        label="Same-element relative radius",
                        expertLevel=LEVEL_ADVANCED,
                        validators=[params.Positive],
                        help='Two sets of coordinates are determined to be of '
                        'the same particle if they are within this radius. The'
                        ' radius is given in [fraction of particle size] units.'
        )
        # Additional data
        form.addSection(label='Additional data')
        form.addParam('doPositiveInput', params.BooleanParam, default=False,
                             label = 'Manually insert positive inputs',
                             help = 'This option enables the input of positive-labelled data '
                             'into the NN training. For example, previously checked or hand '
                             'picked coordinates.'
        )    
        form.addParam('positiveInputSets', params.MultiPointerParam,
                        condition = 'doPositiveInput == True', 
                        pointerClass = SetOfCoordinates3D, allowsNull=True,
                        label = 'Positive references',
                        help = 'Select pickings that are presumed to be true e.g. hand-picked coordinates.'  
        )
        form.addParam('doNegativeInput', params.BooleanParam, default=False,
                             label = 'Manually insert negative inputs',
                             help = 'This option enables the input of negative-labelled data '
                             'into the NN training. For example, previously picked gold or noise.'
        )    
        form.addParam('negativeInputSets', params.MultiPointerParam,
                        condition = 'doNegativeInput == True', 
                        pointerClass = SetOfCoordinates3D, allowsNull=True,
                        label = 'Negative references',
                        help = 'Select pickings that are presumed to be negative e.g. hand-picked noise.'  
        )
        # NN Training modes
        form.addSection(label='Training')
        form.addParam('modelInitialization', params.EnumParam,
            choices = self.FORM_MODEL_TRAIN_TYPELIST_LABELS,
            default = self.MODEL_TRAIN_NEW,
            label = 'Select a model',
            help = 'When set to *%s*, the network will start with a fresh and randomly '
            'initialized model. The option *%s* will let you choose a previously trained '
            'model.'
            % tuple(self.FORM_MODEL_TRAIN_TYPELIST_LABELS)
        )
        form.addParam('skipTraining', params.BooleanParam,
            default = False,
            condition = 'modelInitialization == %s'%self.MODEL_TRAIN_PRETRAIN,
            label = 'Skip training step',
            help = ' When set to *Yes*, the volumes will be directly fed to the model, '
            ' If set to *No*, you must provide a training set of volumes.'
        )
        form.addParam('modelFile', params.FileParam,
            condition = 'modelInitialization == %s'%self.MODEL_TRAIN_PRETRAIN,
            label = 'Saved model file',
            help = 'Select the H5 format filename of the trained neural network. Note that '
            'the architecture must have been compiled using this same protocol.'
        )
        form.addParam('trainingBatch', params.IntParam, default='32',
                label = 'Training batch size',
                help = 'Amount of subtomograms in a training batch. '
                'If the provided subtomograms are not enough for the '
                'NN, this has to be increased. If the machine hangs due'
                ' to memory issues, this has to be reduced.'
        )
        form.addParam('numberEpochs', params.IntParam, default=6,
                        label = 'Cycles (total epochs)',
                        help = 'Number of process cycles that will be done '
                        'with the data in order to train the Neural Network.'
        )
        form.addParam('regulStrength', params.FloatParam, default = 0.00001,
                      label = 'L2 regularisation strength',
                      help = 'Hyperparameter that controls the extra term '
                      'added to the cost function to evade overfitting in '
                      'the trained model'
        )
        form.addParam('learningRatef', params.FloatParam, default = 0.0005,
                        label = 'Learning rate',
                        help = 'Hyperparameter that controls the difference '
                        'between a calculated weight and its next value. '
                        'Higher values will result in faster learning, but '
                        'at the expense of sub-optimal weight values. '
                        'Very low value cause the NN to get stuck in local '
                        'minimas.'
        )
        form.addParam('choiceDynLearningRate', params.BooleanParam, default = True,
                        expertLevel=LEVEL_ADVANCED,
                        label = 'Dynamic learning rate',
                        help = 'The learning rate can be updated on runtime '
                        'depending on the evolution of the execution. '                    
        )
        form.addParam('convergStop', params.BooleanParam, default = True,
                        expertLevel=LEVEL_ADVANCED,
                        label = 'Stop on convergence',
                        help = 'When set to *Yes*, the protocol will stop '
                        'the training when no improvement is detected in '
                        '2 consecutive epochs. Make sure that the learning '
                        'rate is not too low or this option might get your '
                        'model stopped in a sub-optimal, local minima. This'
                        ' is not recommended for small datasets.'
        )
        form.addParam('forceDataAugment', params.BooleanParam, default = True,
                      expertLevel=LEVEL_ADVANCED,
                      label = "Force data augmentation",
                      help = 'By default, the protocol will not try to '
                      'perform data augmentation on datasets if at least two '
                      'of the input pickers contain 900 structures detected.'
        )
        form.addParam('votingMode', params.BooleanParam,
            default = False,
            expertLevel=LEVEL_ADVANCED,
            label = 'Voting instead of DNN',
            help = 'Activating this will ignore all NN parameters and only perform '
            'a consensus based on regular voting after doing the coordinates fusion. '
            'Every NN parameter will be ignored if this option is active.'
        )
        form.addParam('votingThreshold', params.FloatParam,
            default = 0.3,
            label = "Required consensus threshold",
            condition = 'votingMode == True',
            help = 'Sets the required consensus threshold (0,1] ratio needed for '
            'a result to be considered good in simple voting mode. Bear in mind '
            'the amount of input pickers when choosing this value. For instance, '
            'a 0.7 value will not be achievable with two input pickers as the volumes '
            'will appear either in 0.5 or in 1.0 of the total pickers.'
        )
        
    #--------------- INSERT steps functions ----------------
    def _insertAllSteps(self):
        # readInputs and doBSSR consensus need to be outside of the steps system
        # Otherwise there is no global knowledge here of tsIds, etc.
        self.readInputs()  # Serial and out of the steps 
        self.doBSSRConsensus() # Serial and out of the steps
        # From here on it can be steps system
        this = self._insertFunctionStep(self.calculateScalingFactorsStep, prerequisites=[])
        this = self._insertFunctionStep(self.generateScaledCoordinatesStep, prerequisites=[this])
        this = self._insertFunctionStep(self.writeScaledCoordinatesStep, prerequisites=[this]) # Serial because of nature
        deps = []
        threadsPerTomogram = int(self.nThreads) // len(self.allTsIds)
        threadsPerTomogram = threadsPerTomogram if threadsPerTomogram > 0 else 1
        for tsId in self.allTsIds:
            thisStep = self._insertFunctionStep(self.coordConsensusStep, tsId, prerequisites=this)
            deps.append(self._insertFunctionStep(self.noisePickStep, tsId, threadsPerTomogram, prerequisites=thisStep))
        this = self._insertFunctionStep(self.tomogramScalingStep, prerequisites=deps)
        this = self._insertFunctionStep(self.extractionStep, prerequisites=[this])
        # this = self._insertFunctionStep(self.processTrainStep)
        # this = self._insertFunctionStep(self.processScoreStep)
        # this = self._insertFunctionStep(self.postProcessStep)
        # this = self._insertFunctionStep(self.createOutputStep prerequisites=this)

    #--------------- STEPS functions -----------------------

    def readInputs(self) -> None:
        """
        Reads the input form, establishes some hardcoded parameters.
        Uses Scipion datatypes for DB persistence.
        """
        # Compute
        self.nThreads : String = str(self.numberOfThreads.get())

        # Sets of coordinates
        self.inputSetsOf3DCoordinates : List = List([item.get() for item in self.inputSets])
        self.inputSetsOf3DCoordinatesPositive : List = List([item.get() for item in self.positiveInputSets])
        self.inputSetsOf3DCoordinatesNegative : List = List([item.get() for item in self.negativeInputSets])

        # Pickers
        self.nr_pickers : Integer = len(self.inputSetsOf3DCoordinates)
        self.nr_pickersPos : Integer = len(self.inputSetsOf3DCoordinatesPositive)
        self.nr_pickersNeg : Integer = len(self.inputSetsOf3DCoordinatesNegative)
        self.nr_pickersTotal : Integer = self.nr_pickersNeg + self.nr_pickersPos + self.nr_pickers

        # ROIs
        self.nr_ROIs : Integer = sum(map(len, self.inputSetsOf3DCoordinates))
        self.nr_ROIsPos : Integer = sum(map(len, self.inputSetsOf3DCoordinatesPositive))
        self.nr_ROIsNeg : Integer = sum(map(len, self.inputSetsOf3DCoordinatesNegative))
        self.nr_ROIsTotal : Integer = self.nr_ROIs + self.nr_ROIsPos + self.nr_ROIsNeg

        # Are there positive and negative examples?        
        self.havePositive : Boolean = len(self.inputSetsOf3DCoordinatesPositive) > 0        
        self.haveNegative : Boolean = len(self.inputSetsOf3DCoordinatesNegative) > 0

        # Hardcoded parameters
        self.noiseFrac = Float(0.9)
        self.valueConsType = Integer(2)
        self.coordConsType = Integer(1)
        self.valFrac = Float(0.2)
        self.globalParticleId = Integer(0) # Initially zero and incremented during extraction, a bit sketchy but it works

        # Form parameters
        self.coordConsRadius = Float(self.coordConsensusRadius.get()) # Assimilation radius
        self.trainSkip = Boolean(self.skipTraining.get()) # Skip training step?
        self.trainType = Integer(self.modelInitialization.get()) # Model initialisation
        self.nEpochs = Integer(self.numberEpochs.get()) # NN Epochs
        self.regStrength = Float(self.regulStrength.get()) #L1L2 reg str
        self.learningRate = Float(self.learningRatef.get())
        self.dynLearningRate = Boolean(self.choiceDynLearningRate.get())
        self.convergeStop = Boolean(self.convergStop.get()) 
        self.augment = Boolean(self.forceDataAugment.get())
        self.batchSize = Integer(self.trainingBatch.get())
        self.nrPickersNeeded = Integer(self.neededNumberOfPickers.get())

        # Populated in this function
        self.allTsIds = List() # Dedup list of all TSID found in the inputs
        self.allTsIds_raw = List() # Duplicated list of all TSID found in all pickers
        self.allTsIds_common = List() # Contains the intersection of all TSID in the input pickers
        self.allTsIds_filedicts = Dict() # Contains dict(key:tsID , value: dict(key:samplingrate, value: filename))
        
        # Populated in later functions
        self.consBoxSize = Integer()
        self.consSampRate = Float()
        self.inputScalingFactors = List()
        self.scaledTomoDims = List([-1,-1,-1])

        """
        allTsIds_filedicts dictionary of dictionaries
        You can see if a ts tsId is available as a sr samplingrate with -> self.allTsIds_filedicts[tsId][sr]
        { "IS_202009" : {"8.0": "file1.mrc", "13.0": "file2.mrc"},
        "IS_202012" : {"8.0": "filex.mrc", "10.0": "filey.mrc"}
        }
        """
        inputSet : SetOfCoordinates3D
        for inputSet in self.inputSetsOf3DCoordinates:
            # What Ts are found here?
            pickerTsAndFiles = self.getUniqueTsWithFilesInSet(inputSet)
            self.allTsIds_raw.append(list(pickerTsAndFiles.keys()))

            # Incorporate file as a representative of this SampRate
            for key in pickerTsAndFiles.keys():
                # Generate the list if first time with this tsId
                if key not in self.allTsIds_filedicts.keys(): 
                    self.allTsIds_filedicts[key] = Dict() 
                
                sr = inputSet.getSamplingRate()
                self.allTsIds_filedicts[key][sr] = pickerTsAndFiles[key]

        # Flatten list with method for strings list and deduplicate with list(set()) method
        self.allTsIds_raw = sum(self.allTsIds_raw, [])
        self.allTsIds = list(set(self.allTsIds_raw))
        self.allTsIds_common = list(set(tsId for tsId in self.allTsIds if self.allTsIds_raw.count(tsId) == len(self.inputSetsOf3DCoordinates)))

        # Debug print the TSID+available SRs
        print("ALL TomoID availability")
        for key in self.allTsIds_filedicts.keys():
            myList = list(self.allTsIds_filedicts[key].keys())
            print(f"TS ID {key} has these candidates: {myList}.")

        # Print the sets to console output
        print(f"Total unique tomograms: {len(self.allTsIds)} in {len(self.inputSetsOf3DCoordinates)} pickers")
        print(f"Tomograms found in all pickers: {len(self.allTsIds_common)}")

        # Print pickers information
        self.printPickersInfo()

        # Generate internal protocol folders
        folders = [ self._getScaledPath() ]
        folders.append(self._getCoordConsPath())
        folders.append(self._getNoisePath())
        folders.append(self._getScaledTomoPath())
        # folders.append(self._getCoordConsensusPath())
        for f in folders:
            pwutils.makePath(f)
        
    def doBSSRConsensus(self) -> None:
        
        # With all information loaded, now we search for the biggest SR
        allSRs = List()
        allSRs = [inputSet.getSamplingRate() for inputSet in self.inputSetsOf3DCoordinates]
        self.consSampRate = max(allSRs)        
        allBSs = [inputSet.getBoxSize() for inputSet in self.inputSetsOf3DCoordinates]
        self.consBoxSize = min(allBSs)

        self.printBSSRConsensusInfo()

    def calculateScalingFactorsStep(self) -> None:
        """
        Generates a new entry for each Picker MD dictionary: picker_sf (scale factor)
        """
        # Generate all -1 for later error detection capabilities
        self.inputScalingFactors = [-1]*len(self.inputSetsOf3DCoordinates)

        # Calculate the needed scaling factor for each picker according to consensuated SR
        # Only calculates the number, not scaling is performed in this function
        picker : SetOfCoordinates3D
        for index,picker in enumerate(self.inputSetsOf3DCoordinates):
            if picker.getSamplingRate() == self.consSampRate:
                self.inputScalingFactors[index] = 1.00
            else:
                self.inputScalingFactors[index] = picker.getSamplingRate() / self.consSampRate

        assert -1 not in self.inputScalingFactors

        self.printScalingFactorsInfo()

    def generateScaledCoordinatesStep(self) -> None:
        """
        Generates a new entry for a all Picker MD dictionary: picker_coords_scaled
        """
        inputSet : SetOfCoordinates3D
        for index, inputSet in enumerate(self.inputSetsOf3DCoordinates):
            coord : Coordinate3D
            if self.inputScalingFactors[index] != 1.0:
                for coord in inputSet.iterCoordinates():
                    coord.scale(self.inputScalingFactors[index])

        dims = [-1,-1,-1]
        for tsId in self.allTsIds:
            if self.consSampRate in self.allTsIds_filedicts[tsId].keys():
                ih = ImageHandler()
                dims = ih.getDimensions(self.allTsIds_filedicts[tsId][self.consSampRate])
                self.scaledTomoDims[0] = dims[0]
                self.scaledTomoDims[1] = dims[1]
                self.scaledTomoDims[2] = dims[2]

        assert (-1 not in dims) and (-1 not in self.scaledTomoDims)

    def writeScaledCoordinatesStep(self):
        """
        Writes an XMD file with the scaled coordinates from each picker. One file per TSID
        
        Maintains a reference to the original picker but uses tsId for reference and not tomofile
        The resultant XMD can not be previewed in tomogram mode because of that

        But it's okay we only want to save scaled coords + picker of origin
        """
        # part_id : Integer = 0
        myMdDict = dict()
        for tsId in self.allTsIds:
            myMdDict[tsId] = emlib.MetaData()

        row = emlib.metadata.Row()
        inputSet : SetOfCoordinates3D
        for pickerIndex, inputSet in enumerate(self.inputSetsOf3DCoordinates): # For each picker...
            coord : Coordinate3D
            for coord in inputSet.iterCoordinates():
                row.setValue(emlib.MDL_REF, pickerIndex) # Internally used for tracking which picker this comes from
                row.setValue(emlib.MDL_XCOOR, int(coord.getX(tconst.BOTTOM_LEFT_CORNER)))
                row.setValue(emlib.MDL_YCOOR, int(coord.getY(tconst.BOTTOM_LEFT_CORNER)))
                row.setValue(emlib.MDL_ZCOOR, int(coord.getZ(tconst.BOTTOM_LEFT_CORNER)))
                row.addToMd(myMdDict[coord.getTomoId()])

        for tsId in self.allTsIds:
            assert myMdDict[tsId] is not None
            myMdDict[tsId].write(self._getScaledFile(tsId))

    def coordConsensusStep(self, tsId) -> None:
        """
        Launches the consensus of the coordinates per tsId
        Saves into pos, doubt, neg and all coords XMD files. 
        One call per tsId
        """

        program = "xmipp_tomo_coordinates_consensus"

        args = ' --boxsize ' + str(self.consBoxSize)
        args += ' --samplingrate ' + str(self.consSampRate)
        args += ' --radius ' + str(self.coordConsRadius)
        args += ' --number ' + str(self.nrPickersNeeded)
        args += ' --constype ' + str(self.coordConsType)
        args += ' --infile ' + self._getScaledFile(tsId)
        args += ' --outall ' + self._getAllCoordsFilename(tsId)
        args += ' --outpos ' + self._getPosCoordsFilename(tsId)
        args += ' --outdoubt ' + self._getDoubtCoordsFilename(tsId)
        args += ' --outneg ' + self._getNegCoordsFilename(tsId)
        args += ' --tsid ' + str(tsId)   
        
        self.runJob(program, args)

    def noisePickStep(self, tsId: str, threads: int) -> None:
        """
        Launches the noise picking step.
        One call per tsId
        """
        program = "xmipp_tomo_pick_noise"
        
        args = ''
        args += ' --infile ' + self._getPosCoordsFilename(tsId)
        args += ' --output ' + self._getNoiseCoordsFilename(tsId)
        args += ' --boxsize ' + str(self.consBoxSize)
        args += ' --tomosize ' + ' '.join(map(str,self.scaledTomoDims[0:3]))
        args += ' --threads ' + str(threads)
        args += ' --tsid ' + str(tsId)  
        
        self.runJob(program, args)

    def tomogramScalingStep(self) -> None:
        """
        This program analyzes the current tsId panorama and decides if any of them needs to be
        generated parting from a higher-resolution one. If so, it performs this scaling.
        If not, it does nothing to that particular tsId.
        """
        count = 0
        program = "xmipp_image_resize"
        params = ''
        start = time.time()
        # Perform scaling if needed for every tsId that is present in the protocol
        for tsId in self.allTsIds:
            if self.consSampRate not in self.allTsIds_filedicts[tsId].keys():
                count += 1
                closestSR = max(self.allTsIds_filedicts[tsId].keys())
                inputTomo = self.allTsIds_filedicts[tsId][closestSR]
                outputTomo = self._getScaledTomoFilename(tsId)
                params += ' -i ' + inputTomo
                params += ' -o ' + outputTomo
                params += ' --factor ' + str(closestSR / self.consSampRate)
                self.runJob(program, params)
                # Add the newly created tomogram to the list
                self.allTsIds_filedicts[tsId][self.consSampRate] = outputTomo
        # Calculate elapsed wall time and print
        end = time.time()
        print(f"tomogramScalingStep scaled {count} tomograms in {end - start} seconds.")

    def extractionStep(self) -> None:
        """
        Performs the extraction of the subtomograms from the previously generated coords files.
        It is done in a per-tsId basis.
        """        
        for tsId in self.allTsIds:
            tomoPath = self.allTsIds_filedicts[tsId][self.consSampRate]
            posPath = self._getExtractedSubtomosPathPos()
            negPath = self._getExtractedSubtomosPathNeg()
            doubtPath = self._getExtractedSubtomosPathDoubt()
            coordsFile_pos = self._getPosCoordsFilename(tsId)
            coordsFile_neg = self._getNegCoordsFilename(tsId)
            coordsFile_doubt = self._getDoubtCoordsFilename(tsId)
            
            files_coords = [coordsFile_pos, coordsFile_doubt, coordsFile_neg]
            paths_extract = [posPath, doubtPath, negPath]

            program = "xmipp_tomo_extract_subtomograms"
            for coordsFile, outPath in zip(files_coords,paths_extract):
                args = ' '
                args += '--tomogram ' + tomoPath
                args += '--coordinates ' + coordsFile
                args += '--boxsize ' + str(self.consBoxSize)
                args += '-o ' + outPath
                args += '--threads ' + str(self.nThreads)
                args += '--normalize '
                args += '--factor 1.0 '
                self.runJob(program, args)
            
    def processTrainStep(self):
        """
        Block 2 - Call for NN train

        Launches the script responsible for training the NN.
        """
        program = "xmipp_deep_picking_consensus_tomo"
        # program = "xmipp_deep_picking_consensus_tomo"
        args = ''
        args += ' -t ' + str(self.nThreads)
        args += ' -g ' + ','.join(map(str, self.getGpuList()))
        args += ' --mode training'
        args += ' --batchsize ' + str(self.batchSize)
        args += ' --netpath ' + self._getNnPath()
        args += ' --consboxsize ' + str(self.consBoxSize)
        args += ' --conssamprate ' + str(self.consSampRate)
        args += ' --ttype ' + str(self.trainType)
        args += ' --valfrac ' + str(self.valFrac)
        args += ' --truevolpath ' + self._getPosSubtomogramPath()
        args += ' --falsevolpath ' + self._getNegSubtomogramPath()
        args += ' -e ' + str(self.nEpochs)
        args += ' -l ' + str(self.learningRate)
        args += ' -r ' + str(self.regStrength)
        args += ' --ensemble 1'
        if not self.convergeStop:
            args += ' -s'
        print('\nHanding over to Xmipp program for Train')
        self.runJob(program, args)

    def processScoreStep(self):
        """
        Block 2 - Neural Network SCORE operations

        Load the selected NN model, extract the dubious
        subtomograms and score them against the model. Then
        generate the score tables.
        """
        program = "conda run -n xmipp_DLTK_v1.0 xmipp_deep_picking_consensus_tomo"
        # program = "xmipp_deep_picking_consensus_tomo"
        args = ''
        args += ' -t ' + str(self.nThreads)
        args += ' -g ' + ','.join(map(str, self.getGpuList()))
        args += ' --mode scoring'
        args += ' --batchsize ' + str(self.batchSize)
        args += ' --netpath ' + self._getNnPath()
        args += ' --netname ' + "dpc_nn.h5"
        args += ' --consboxsize ' + str(self.consBoxSize)
        args += ' --conssamprate ' + str(self.consSampRate)
        # TODO: change this to something more sensible based on folder structure por favor
        args += ' --inputvolpath ' + self._getCombinedDatasetFile()
        args += ' --outputpath ' + self._getOutputFile()

        print('\nHanding over to Xmipp program for Score')
        self.runJob(program, args)

    def postProcessStep(self):
        """
        Block 3 - Handling of posterous things

        Apply the desired filters and thresholds to the 
        results.
        """
        # TODO: fully redo

        # Open the scored output XMD file
        md = emlib.MetaData(self._getOutputFile())
        md_filtered = emlib.MetaData()
        
        # Generate the dataframe
        cols_out = ['subtomoname','identifier','boxsize','samplingrate','tomoname','x','y','z','score']
        size_paco = md.size()

        self.scoredDf = pd.DataFrame(index=range(size_paco),columns=cols_out)

        # Load the data into memory in the protocol realm
        myIndex = 0
        for row in md:
            subtomoName = str(md.getValue(emlib.MDL_IMAGE, row))
            partId = int(md.getValue(emlib.MDL_PARTICLE_ID, row))
            bsize = int(md.getValue(emlib.MDL_PICKING_PARTICLE_SIZE, row))
            srate = float(md.getValue(emlib.MDL_SAMPLINGRATE, row))
            tomoName = md.getValue(emlib.MDL_TOMOGRAM_VOLUME, row)
            c_x = int(md.getValue(emlib.MDL_XCOOR, row))
            c_y = int(md.getValue(emlib.MDL_YCOOR, row))
            c_z = int(md.getValue(emlib.MDL_ZCOOR, row))
            score = float(md.getValue(emlib.MDL_ZSCORE, row))

            self.scoredDf.loc[myIndex, 'subtomoname'] =  subtomoName
            self.scoredDf.loc[myIndex, 'identifier'] =  partId
            self.scoredDf.loc[myIndex, 'boxsize'] =  bsize
            self.scoredDf.loc[myIndex, 'samplingrate'] =  srate
            self.scoredDf.loc[myIndex, 'tomoname'] =  tomoName
            self.scoredDf.loc[myIndex, 'x'] =  c_x
            self.scoredDf.loc[myIndex, 'y'] =  c_y
            self.scoredDf.loc[myIndex, 'z'] =  c_z
            self.scoredDf.loc[myIndex, 'score'] =  score

            # Add also to the filtered XMD if passes filter
            if score >= self.outScoreThreshold:
                where = md_filtered.addObject()
                md_filtered.setRow(md.getRow(row), where)
            
            myIndex += 1
        md_filtered.write(self._getOutputFilteredFile())
    
        # Apply the filter to the column of the zscore
        self.scoredDf_filtered : pd.DataFrame = self.scoredDf[self.scoredDf['score'] >= self.outScoreThreshold]

        # Deduplicate the ones that are literally the same
        intermediate : pd.DataFrame = self.scoredDf_filtered.drop_duplicates(subset=['x','y','z','tomoname']).reset_index()
        print("Equal dedup: from %d to %d" %(len(self.scoredDf_filtered), len(intermediate)))
        # No descomentes esto hasta haber arreglado dedup
        # self.scoredDf_filtered_dedup : pd.DataFrame = self.dedup(intermediate)
        self.scoredDf_filtered_dedup = intermediate
    
    def createOutputStep(self):
        """
        Block 3 - Prepare output

        Leave everything in a format suitable for later processing
        in the CryoET workflows.
        """
        # Get SetOfTomograms
        # tomos = [Tomogram(id) for id in self.uniqueTomoIDs]
        # tomograms : SetOfTomograms = SetOfTomograms(tomos)
        # outputPath = self._getExtraPath("CBOX_3D")
        suffix = self._getOutputSuffix(SetOfCoordinates3D)

        tomoDict = {}
        setOfCoords : SetOfCoordinates3D
        outTomos : SetOfTomograms
        for setOfCoords in self.inputSetsOf3DCoordinates:
            precedent : Tomogram
            precedents = setOfCoords.getPrecedents()
            listaguay = [elem.clone() for elem in precedents]
            # TODO: No hacer esto con outTomos en el futuro, usarlos del tomodict
            outTomos = precedents
            for precedent in listaguay:
                tsid = precedent.getTsId()
                tomoDict[tsid] = precedent
            

        print("Mi diccionario es asin de grande:")
        print(tomoDict.keys())

        coordinates : SetOfCoordinates3D = self._createSetOfCoordinates3D(outTomos, suffix)
        coordinates.setName("output coordinates - scored")
        coordinates.setSamplingRate(self.consSampRate)
        coordinates.setBoxSize(self.consBoxSize)
        for ind in self.scoredDf_filtered_dedup.index:
            c : Coordinate3D = Coordinate3D()
            tsid = self._getTsIdFromName(self.scoredDf_filtered_dedup['tomoname'][ind])
            t : Tomogram = tomoDict[tsid]
            c.setVolume(t)
            c.setX(self.scoredDf_filtered_dedup['x'][ind], tconst.BOTTOM_LEFT_CORNER)
            c.setY(self.scoredDf_filtered_dedup['y'][ind], tconst.BOTTOM_LEFT_CORNER)
            c.setZ(self.scoredDf_filtered_dedup['z'][ind], tconst.BOTTOM_LEFT_CORNER)
            # TODO: set score with new function
            # c.setScore()
            coordinates.append(c)
        
        name = self.OUTPUT_PREFIX + suffix
        self._defineOutputs(**{name: coordinates})

        inset : SetOfCoordinates3D
        for inset in self.inputSets:
            self._defineSourceRelation(inset, coordinates)
        
    #--------------- INFO functions -------------------------
    def _summary(self):
        summary = []
        summary.append(f"- Executed in {self.nThreads} threads and {len(self.getGpuList())}.")
        summary.append(f"- Found {self.nr_pickers} input pickers in total.")
        summary.append(f"- Found {len(self.allTsIds)} tsId from which {len(self.allTsIds_common)} are present in all pickers.")
        summary.append(f"- Used {self.consSampRate}A/px as the common sampling rate and {self.consBoxSize}px as box size.")
        # summary.append(f"- Total input coordinates: {self.nr_ROIsTotal}. Total output coordinates: {self.nr_ROIsTotal_out}.")
        # summary.append(f"- ")
        # summary.append("- Coordinates Consensus removed KKQLO duplicates")
        return summary

    #--------------- UTILS functions -------------------------

    def getUniqueTsWithFilesInSet(self, inSet : SetOfCoordinates3D) -> dict:
        '''
        Get the unique tomograms with filenames in a single dictionary
        '''
        # Get the list of tomograms in this set, ALL including unused ones
        precedents : SetOfTomograms = inSet.getPrecedents()
        """
        Dictionary contains
          - *tsid_name* : String - tomogram filename
        """
        res = dict()
        tomo : Tomogram
        for tomo in precedents:
            res[tomo.getTsId()] = tomo.getFileName()
        return res

    def printPickersInfo(self) -> None:
        '''
        Prints the data of the input pickers, one by one.
        Does not print the coordinates, only the amount of them.
        '''
        print("BEGIN Input pickers information")
        picker: SetOfCoordinates3D
        for index, picker in enumerate(self.inputSetsOf3DCoordinates):
            print("------------------------")
            print(f"ID: {index}")
            print(f"SamplinRate: {picker.getSamplingRate():.2f}A/px")
            print(f"BoxSize: {picker.getBoxSize()}px")
            print(f"Pickings: {picker.getSize()}")
        print("END Input pickers information")
    
    def printBSSRConsensusInfo(self) -> None:
        '''
        Prints some lines with the actual BSSR consensus.
        '''
        assert self.consBoxSize is not None
        assert self.consSampRate is not None

        print("BEGIN BSSR information")
        print("------------------------")
        if len(self.allTsIds_common) < len(self.allTsIds):
            print("=======================================================================================")
            print(f"= HEY! IRREGULAR N OF TOMOs IN PICKERS. TOMOGRAMS WILL GET DOWNSAMPLED TO {self.consSampRate:.2f}A/px   =")
            print("=======================================================================================")
        print(f"Consensus Box Size: {self.consBoxSize}px")
        print(f"Consensus Sampling Rate: {self.consSampRate:.2f}A/px")
        print("------------------------")
        print("END BSSR information")
    
    def printScalingFactorsInfo(self) -> None:
        """
        Prints scale factor for all pickers
        """
        print("BEGIN ScalingFactor information")
        print("------------------------")
        for i in range(len(self.inputSetsOf3DCoordinates)):
            sf = self.inputScalingFactors[i]
            print(f"Picker with ID {i} has scaling of {sf:.2f}")
        print("------------------------")
        print("END ScalingFactor information")
        
    #--------------- VALIDATE functions ----------------------
    def _validate(self) -> list:
        errors = []
        errors += self._validateParallelProcessing()
        errors += self._validateNrOfPickersNeeded()
        errors += self._validateExistingModel()
        return errors

    def _validateParallelProcessing(self) -> list:
        nGpus = len(self.getGpuList())
        nMPI = self.numberOfMpi.get()
        errors = []
        if nGpus < 1:
            errors.append("A GPU is needed for this protocol to run.")
        if nMPI != 1:
            errors.append("Multiprocessing not yet supported. Set MPI parameter to 1 and use Threads instead.")
        return errors
    
    def _validateNrOfPickersNeeded(self) -> list:
        howManyInputs = len(self.inputSets.get())
        howManySelected = int(self.neededNumberOfPickers.get())
        errors = []
        if howManySelected > howManyInputs:
            errors.append("Required pickers number must be lesser or equal to the number of input pickers.")
        return errors
    
    def _validateExistingModel(self) -> list:
        execMode = self.modelInitialization.get()
        errors = []
        if execMode == self.MODEL_TRAIN_PRETRAIN:
            inFile = self.modelFile.get()
            if not os.path.isfile(inFile):
                errors.append("Can not access provided trained model file. Check route and permissions.")
        return errors
        
    #--------------- FILENAMES functions -------------------

    # Extracted subtomogram file management
    def _getExtractedSubtomosPath(self, *args) -> str:
        return self._getExtraPath('extracted', *args)
    def _getExtractedSubtomosPathPos(self, *args) -> str:
        return self._getExtractedSubtomosPath('pos', *args)
    def _getExtractedSubtomosPathNeg(self, *args) -> str:
        return self._getExtractedSubtomosPath('neg', *args)
    def _getExtractedSubtomosPathDoubt(self, *args) -> str:
        return self._getExtractedSubtomosPath('doubt', *args)

    # Scaled tomograms file management
    def _getScaledTomoPath(self, *args) -> str:
        return self._getExtraPath('scaledTomos', *args)
    def _getScaledTomoFilename(self, tsId: str) -> str:
        return self._getScaledTomoPath(tsId)

    # Scaled coordinates file management
    def _getScaledPath(self, *args) -> str:
        return self._getExtraPath('scaled', *args)
    def _getScaledFile(self, tsId : str) -> str:
        return self._getScaledPath(tsId + '_scaled_coords_all.xmd')

    # Consensus coordinates file management
    def _getCoordConsPath(self, *args) -> str:
        return self._getExtraPath('coordcons', *args)
    def _getAllCoordsFilename(self, tsId: str) -> str:
        return self._getCoordConsPath(tsId+"_cons_all.xmd")
    def _getPosCoordsFilename(self, tsId: str) -> str:
        return self._getCoordConsPath(tsId+"_cons_pos.xmd")
    def _getDoubtCoordsFilename(self, tsId: str) -> str:
        return self._getCoordConsPath(tsId+"_cons_doubt.xmd")
    def _getNegCoordsFilename(self, tsId: str) -> str:
        return self._getCoordConsPath(tsId+"_cons_neg.xmd")
    
    # Noise picking file management
    def _getNoisePath(self, *args) -> str:
        return self._getExtraPath('noisepick', *args)
    def _getNoiseCoordsFilename(self, tsId: str) -> str:
        return self._getNoisePath(tsId+"_noise.xmd")
    
    # MIERDA OLD

    # def _getOutputPath(self, *args):
    #     return self._getExtraPath('out', *args)
    
    # def _getOutputFile(self, *args):
    #     return self._getOutputPath('nn_output_scored.xmd', *args)
    
    # def _getOutputFilteredFile(self, *args):
    #     return self._getOutputPath('nn_output_scored_filtered.xmd', *args)

    # def _getNnPath(self, *args):
    #     return self._getExtraPath('nn', *args)

    # def _getNnResultFilename(self, *args):
    #     return self._getNnPath("nn_res")

    # def _stripTomoFilename(self, tomopath: str):
    #     return tomopath.split("/")[-1].split(".")[0]
    
    
    
    # def _getPosCoordsScaledFilename(self, tomo_name: str):
    #     return self._getCoordConsensusPath(tomo_name+"_cons_pos_scl.xmd")
    
    # def _getNegCoordsScaledFilename(self, tomo_name: str):
    #     return self._getCoordConsensusPath(tomo_name+"_cons_neg_scl.xmd")
    
    # def _getDoubtCoordsScaledFilename(self, tomo_name: str):
    #     return self._getCoordConsensusPath(tomo_name+"_cons_doubt_scl.xmd")
    
    # def _getConsCoordsFilename(self, tomo_name:str):
    #     return self._getCoordConsensusPath(tomo_name+"_cons_all.xmd")
    
    # def _getAllCoordsFilename(self, tomo_name: str):
    #     return self._getPickedPerTomoPath(tomo_name+"_allpickedcoords.xmd")
    
    # def _getAllTruthCoordsFilename(self, tomo_name: str):
    #     return self._getPickedPerTomoPath(tomo_name+"_truth.xmd")
    
    # def _getAllLieCoordsFilename(self, tomo_name: str):
    #     return self._getPickedPerTomoPath(tomo_name+"_lie.xmd")

    # def _getDatasetPath(self, *args):
    #     return self._getExtraPath('dataset', *args)
    
    # def _getPosSubtomogramPath(self, *args):
    #     return self._getDatasetPath('pos', *args)
    
    # def _getNegSubtomogramPath(self, *args):
    #     return self._getDatasetPath('neg', *args)
    
    # def _getDoubtSubtomogramPath(self, *args):
    #     return self._getDatasetPath('doubt', *args)

    # def _getCombinedDatasetFile(self, *args):
    #     return self._getDoubtSubtomogramPath('combined.xmd', *args)

    # def _getFilenameFromTsId(self, *args):
    #     pass

    # def _getPickedPerTomoPath(self, *args):
    #     return self._getExtraPath('pickedpertomo', *args)
