# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Martín Salinas Antón (ssalinasmartin@gmail.com)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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

# General imports
import os, math
from typing import Tuple, Union, TypedDict
from emtable import Table

# Scipion em imports
from pwem.protocols import EMProtocol
from pwem.objects import SetOfParticles, CTFModel, Float, Volume
from pwem.emlib import MetaData, metadata, MDL_CTF_PHASE_SHIFT, MDL_CTF_DEFOCUSU
from pwem.emlib import MDL_CTF_DEFOCUSV, MDL_CTF_DEFOCUS_ANGLE, MDL_ANGLE_TILT
from pyworkflow import BETA
from pyworkflow.protocol import params
from pyworkflow.utils import Message

# External plugin imports
from tomo.protocols import ProtTomoBase
from tomo.objects import SubTomogram, Coordinate3D, TiltSeries, TomoAcquisition, SetOfSubTomograms
from xmipp3.convert import readSetOfParticles

# Plugin imports
from ..utils import calculateRotationAngleAndShiftsFromTM

# Protocol output variable name
OUTPUTATTRIBUTE = 'outputSetOfParticles'

class XmippProtProjectSubtomograms(EMProtocol, ProtTomoBase):
    """Extracts proyections from subtomograms"""

    # Protocol constants
    _label = 'project subtomograms'
    _devStatus = BETA
    _possibleOutputs = {OUTPUTATTRIBUTE: SetOfParticles}

    # Form constants
    METHOD_FOURIER = 0
    METHOD_REAL_SPACE = 1
    METHOD_SHEARS = 2
    TYPE_N_SAMPLES = 0
    TYPE_STEP = 1
    TYPE_TILT_SERIES = 2
    INTERPOLATION_BSPLINE = 0
    INTERPOLATION_NEAREST = 1
    INTERPOLATION_LINEAR = 2

    # Emtable star file column names (temporary)
    COLUMN_ANGLE_PSI = 'anglePsi'
    COLUMN_ANGLE_ROT = 'angleRot'
    COLUMN_ANGLE_TILT = 'angleTilt'

    # --------------------------- Class constructor --------------------------------------------
    def __init__(self, **args):
        # Calling parent class constructor
        super().__init__(**args)

        # Defining execution mode. Steps will take place in parallel now
        # Full tutorial on how to parallelize protocols can be read here:
        # https://scipion-em.github.io/docs/release-3.0.0/docs/developer/parallelization.html
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        # Defining parallel arguments
        form.addParallelSection(threads=4)

        # Generating form
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSubtomograms', params.PointerParam, pointerClass="SetOfSubTomograms,SetOfVolumes", important=True,
                      label='Set of subtomograms', help="Set of subtomograms whose projections will be generated.")
        form.addParam('hasCtfCorrected', params.BooleanParam, default=True, label='Is CTF corrected?: ',
                        help='Set this option to True if the input set of subtomograms has no CTF or the CTF has been corrected.')
        form.addParam('inputCTF', params.PointerParam, pointerClass="SetOfCTFTomoSeries", condition='not hasCtfCorrected',
                      label='CTF', help="Select the CTF corresponding to this set of subtomograms.")
        form.addParam('defocusDir', params.BooleanParam, default=True, label='Defocus increase with z positive?: ', condition='not hasCtfCorrected',
                        help='This flag must be put if the defocus increases or decreases along the z-axis. This is required to set the local CTF.')
        form.addParam('cleanTmps', params.BooleanParam, default=True, label='Clean temporary files: ', expertLevel=params.LEVEL_ADVANCED,
                        help='Clean temporary files after finishing the execution.\nThis is useful to reduce unnecessary disk usage.')
        form.addParam('transformMethod', params.EnumParam, display=params.EnumParam.DISPLAY_COMBO, default=self.METHOD_FOURIER,
                        choices=['Fourier', 'Real space', 'Shears'], label="Transform method: ", expertLevel=params.LEVEL_ADVANCED,
                        help='Select the algorithm that will be used to obtain the projections.')
        
        # Parameter group for fourier transform method
        fourierGroup = form.addGroup('Fourier parameters', condition=f"transformMethod=={self.METHOD_FOURIER}", expertLevel=params.LEVEL_ADVANCED)
        fourierGroup.addParam('pad', params.IntParam, default=2, label="Pad: ", help="Controls the padding factor.")
        fourierGroup.addParam('maxfreq', params.FloatParam, default=0.25, label="Maximum frequency: ",
                                help="Maximum frequency for the pixels.\nBy default, pixels with frequency more than 0.25 are not considered.")
        fourierGroup.addParam('interp', params.EnumParam, default=self.INTERPOLATION_BSPLINE, label="Interpolation method: ", choices=["BSpline", "Nearest", "Linear"],
                                help="Method for interpolation.\nOptions:\n\nBSpline: Cubic BSpline\nNearest: Nearest Neighborhood\nLinear: Linear BSpline")
        
        # Tilt related parameter group
        tiltGroup = form.addGroup('Tilt parameters')
        tiltGroup.addParam('tiltTypeGeneration', params.EnumParam, display=params.EnumParam.DISPLAY_COMBO, default=self.TYPE_N_SAMPLES,
                        choices=['NSamples', 'Step', 'Tilt Series'], label="Type of sample generation: ",
                        help='Select the method for generating samples:\n\n'
                            '*NSamples*: N samples are generated homogeneously across the whole tilt range.\n'
                            '*Step*: For the whole tilt range, a sample is generated every N gedrees.\n'
                            '*Tilt Series*: Given a set of Tilt Series, the angles at which each Tilt Series was taken are used.')
        tiltGroup.addParam('tiltRangeNSamples', params.IntParam, condition=f'tiltTypeGeneration=={self.TYPE_N_SAMPLES}', label='Number of samples:',
                        help='Number of samples to be produced.\nIt has to be 1 or greater.')
        tiltGroup.addParam('tiltRangeStep', params.IntParam, condition=f'tiltTypeGeneration=={self.TYPE_STEP}', label='Step:',
                        help='Number of degrees each sample will be separated from the next.\nIt has to be greater than 0.')
        tiltGroup.addParam('tiltRangeTS', params.PointerParam, pointerClass="SetOfTiltSeries",  condition=f'tiltTypeGeneration=={self.TYPE_TILT_SERIES}',
                           label='Set of Tilt Series:', help='Set of Tilt Series where the angles of each Tilt Series will be obtained for the projection.')
        tiltLine = tiltGroup.addLine("Tilt range (degrees)", condition=f'tiltTypeGeneration!={self.TYPE_TILT_SERIES}',
                        help='The initial and final values of the range of angles the projection will be produced on.\nDefaults to -60º for initial and 60º for final.')
        tiltLine.addParam('tiltRangeStart', params.IntParam, default=-60, label='Start: ')
        tiltLine.addParam('tiltRangeEnd', params.IntParam, default=60, label='End: ')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        # Defining list of function ids to be waited by the createOutput function
        paramDeps = []

        if self.tiltTypeGeneration.get() == self.TYPE_TILT_SERIES:
            # Obtaining list of angle tilts and rots per Tilt Series
            angles = self.getAngleDictionary()

            # Generating a metadata angle file and a param file for each Tilt Series
            for tsId in angles:
                # Getting metadata angle file path
                angleFile = self.getAngleFileAbsolutePath(tsId)

                # Generating metadata angle file with Tilt Series's data
                angTable = Table(columns=[self.COLUMN_ANGLE_PSI, self.COLUMN_ANGLE_ROT, self.COLUMN_ANGLE_TILT])
                for values in angles[tsId]:
                    angTable.addRow(
                        0.0,                            # anglePsi
                        values[self.COLUMN_ANGLE_ROT],  # angleRot
                        values[self.COLUMN_ANGLE_TILT]  # angleTilt
                    )
                angTable.write(angleFile, tableName='projectionAngles')
                
                # Generating param file
                paramDeps.append(self._insertFunctionStep(self.generateParamFile, tsId))

        else:
            # If type of generation is not from Tilt Series, generate single param file
            paramDeps.append(self._insertFunctionStep(self.generateParamFile))
        
        # Generating projections for each subtomogram
        generationDeps = []
        for subtomogram in self.inputSubtomograms.get():
            tsId = '' if self.tiltTypeGeneration.get() != self.TYPE_TILT_SERIES or len([ts for ts in self.tiltRangeTS.get()]) == 1 else subtomogram.getCoordinate3D().getTomoId()
            generationDeps.append(self._insertFunctionStep(self.generateSubtomogramProjections, subtomogram.getFileName(), tsId, prerequisites=paramDeps))
        
        # Conditionally removing temporary files
        if self.cleanTmps.get():
            self._insertFunctionStep(self.removeTempFiles, prerequisites=generationDeps)
        
        # Create output
        self._insertFunctionStep(self.createOutputStep, prerequisites=generationDeps)

    # --------------------------- STEPS functions --------------------------------------------
    def generateParamFile(self, tsId: str=''):
        """
        This function writes the config file for Xmipp Phantom.
        """
        # Generating file content
        content = '# XMIPP_STAR_1 *\n'
        content += '# Projection Parameters\n'
        content += 'data_block1\n'
        content += '# X and Y projection dimensions [Xdim Ydim]\n'
        content += '_dimensions2D   \'{}\'\n'.format(self.getSubtomogramDimensions())
        if self.tiltTypeGeneration.get() == self.TYPE_TILT_SERIES:
            content += "# Projection angles file. Used when projection angles are read instead of calculated\n"
            content += f"_projAngleFile {self.getAngleFileAbsolutePath(tsId)}\n"
        else:
            content += '# Rotation range and number of samples [Start Finish NSamples]\n'
            content += '_projRotRange    \'0 0 1\'\n'
            content += '# Rotation angle added noise  [noise (bias)]\n'
            content += '_projRotNoise   \'0\'\n'
            content += '# Tilt range and number of samples for Tilt [Start Finish NSamples]\n'
            content += '_projTiltRange    \'{} {} {}\'\n'.format(self.tiltRangeStart.get(), self.tiltRangeEnd.get(), self.getStepValue())
            content += '# Tilt angle added noise\n'
            content += '_projTiltNoise   \'0\'\n'
            content += '# Psi range and number of samples\n'
            content += '_projPsiRange    \'0 0 0\'\n'
            content += '# Psi added noise\n_projPsiNoise   \'0\'\n'
            content += '# Noise\n'

        # Writing content to file and closing
        with open(self.getXmippParamPath(tsId=tsId), "w") as paramFile:
            paramFile.write(content)

    def generateSubtomogramProjections(self, subtomogram: Union[SubTomogram, Volume, str], tsId: str=''):
        """
        This function generates the projection for a given input subtomogram.
        """
        params = '-i {}'.format(self.getSubtomogramAbsolutePath(subtomogram))   # Path to subtomogram
        params += ' -o {}'.format(self.getProjectionAbsolutePath(subtomogram))  # Path to output projection
        params += ' --method {}'.format(self.getMethodValue())                  # Projection algorithm
        params += ' --params {}'.format(self.getXmippParamPath(tsId=tsId))      # Path to Xmipp param file
        self.runJob("xmipp_phantom_project", params)
    
    def renameMetadataFile(self, metadataFile: str):
        """
        This function renames the given metadata file with a generic name.
        """
        os.rename(metadataFile, self.getReferenceMetadataFile())

    def removeTempFiles(self):
        """
        This function removes the temporary files of this protocol.
        """
        # Creating list of things to remove
        removeList = []

        if self.tiltTypeGeneration.get() == self.TYPE_TILT_SERIES:
            # For each Tilt Series, add angle file and param file
            for tsId in self.getAngleDictionary():
                removeList.append(self.getAngleFileAbsolutePath(tsId))
                removeList.append(self.getXmippParamPath(tsId=tsId))
        else:
            # Only one Xmipp param file
            removeList.append(self.getXmippParamPath())

        # Removing items
        self.runJob('rm -f', ' '.join(removeList))

    def createOutputStep(self):
        """
        This function generates the outputs of the protocol.
        """
        # Extracting input
        inputSubtomograms = self.inputSubtomograms.get()

        # Creating empty set of particles and setting sampling rate, alignment, and dimensions
        outputSetOfParticles = self._createSetOfParticles()
        outputSetOfParticles.setSamplingRate(inputSubtomograms.getSamplingRate())
        outputSetOfParticles.setAlignmentProj()
        dimensions = self.getSubtomogramDimensions().split(' ')
        outputSetOfParticles.setDim((int(dimensions[0]), int(dimensions[1]), 1))

        # Setting acquisition info
        acquisition = TomoAcquisition()
        acquisition.copyInfo(inputSubtomograms.getAcquisition())
        outputSetOfParticles.setAcquisition(acquisition)

        # Getting input element list
        inputList = [subtomogram.getFileName() for subtomogram in inputSubtomograms]

        # Adding projections of each subtomogram as a particle each
        for subtomogram in inputList:
            # Setting CTF for every output particle
            mdCtf = MetaData(self.getProjectionMetadataAbsolutePath(subtomogram))
            for row in metadata.iterRows(mdCtf):
                # If CTF does not exist or has been corrected, set some values to 0
                row.setValue(MDL_CTF_PHASE_SHIFT, 0.0)
                if self.hasCtfCorrected:
                    row.setValue(MDL_CTF_DEFOCUSU, 0.0)
                    row.setValue(MDL_CTF_DEFOCUSV, 0.0)
                    row.setValue(MDL_CTF_DEFOCUS_ANGLE, 0.0)
                else:
                    # Calculate subtomogram coordinates on the Tilt Series
                    defU, defV = self.getCorrectedDefocus(row.getValue(MDL_ANGLE_TILT), subtomogram.getCoordinate3D())
                    row.setValue(MDL_CTF_DEFOCUSU, defU)
                    row.setValue(MDL_CTF_DEFOCUSV, defV)
                    # TODO: This only works if the TS is aligned
                    row.setValue(MDL_CTF_DEFOCUS_ANGLE, 0.0)
                
                # Add metadata row to file
                row.addToMd(mdCtf)
            
            # Write metadata file with modified info
            mdCtf.write(self.getProjectionMetadataAbsolutePath(subtomogram))

            # Add particles to set from metadata file
            readSetOfParticles(self.getProjectionMetadataAbsolutePath(subtomogram), outputSetOfParticles)
        
        # Defining the ouput with summary and source relation
        outputSetOfParticles.setObjComment(self.getSummary(outputSetOfParticles))
        self._defineOutputs(outputSetOfParticles=outputSetOfParticles)
        self._defineSourceRelation(self.inputSubtomograms, outputSetOfParticles)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        """
        This method validates the received params and checks that they all fullfill the requirements needed to run the protocol.
        """
        errors = []

        # Checking if number of samples is greater or equal than 1
        if self.tiltTypeGeneration.get() == self.TYPE_N_SAMPLES and (self.tiltRangeNSamples.get() != None) and self.tiltRangeNSamples.get() < 1:
            errors.append('The number of samples cannot be less than 1.')
        
        # Checking if the step is greater than 0
        if self.tiltTypeGeneration.get() == self.TYPE_STEP and (self.tiltRangeStep.get() != None) and self.tiltRangeStep.get() <= 0:
            errors.append('The step must be greater than 0.')
        
        # Checking if input projectable set can be related to Tilt Series when TS is selected as angle extraction method
        angleExtractionError = self._validateGenerationType()
        if angleExtractionError:
            errors.append(angleExtractionError)
        
        # Checking if MPI is selected (only threads are allowed)
        if self.numberOfMpi > 1:
            errors.append('MPI cannot be selected, because Scipion is going to drop support for it. Select threads instead.')

        return errors
    
    # --------------------------- Auxiliary INFO functions --------------------------------------------
    def _validateGenerationType(self):
        # Initializing base error string
        baseErrorStr = 'In order to use Tilt Series as a method of angle extraction, at least one of the following conditions must be met:\n'
        baseErrorStr += '\t- There must be only one Tilt Series within the input set of Tilt Series.\n'
        baseErrorStr += '\t- The input subtomograms must be subtomograms (not volumes) and they must contain coordinates that lead to one of the input Tilt Series.\n'
        baseErrorStr += 'In this case, the following validation error occured:\n'

        # Initializing error string
        error = ''

        # Checking if input projectable set can be related to Tilt Series when TS is selected as angle extraction method
        if self.tiltTypeGeneration.get() == self.TYPE_TILT_SERIES and len([ts for ts in self.tiltRangeTS.get()]) > 1:
            # Checking if input subtomograms are of class SetOfSubTomograms (SetOfVolumes can't have a Coordinate3D to link them to a TiltSeries)
            if not isinstance(self.inputSubtomograms.get(), SetOfSubTomograms):
                error = 'Input subtomograms are not of a real set of subtomograms.\n'
            # Checking if input subtomograms do have a Coordinate3D and such Coordinate3D has a tomoId that matches a tsId
            else:
                # Getting list of Tilt Series ids
                tsIdList = [ts.getTsId() for ts in self.tiltRangeTS.get()]

                for subtomogram in self.inputSubtomograms.get():
                    if not subtomogram.getCoordinate3D():
                        error = 'At least one input subtomogram does not have a 3D coordinate assigned, so it cannot be traced back to a Tilt Series.\n'
                        break
                    elif subtomogram.getCoordinate3D().getTomoId() not in tsIdList:
                        error = 'At least one input subtomogram\'s 3D coordinate does not match any Tilt Series id.\n'
                        break
            
        return baseErrorStr + error if error else ''

    # --------------------------- UTILS functions --------------------------------------------
    def scapePath(self, path: str) -> str:
        """
        This function returns the given path with all the spaces in folder names scaped to avoid errors.
        """
        # os.path.abspath adds '\\' when finding a foldername with '\ ', so '\\\' needs to be replaced with ''
        # Then, '\' is inserted before every space again, to include now possible folders with spaces in the absolute path
        return path.replace('\\\ ', ' ').replace(' ', '\ ')
    
    def getSubtomogramRelativePath(self, subtomogram: Union[SubTomogram, Volume, str]) -> str:
        """
        This method returns a the subtomogram path relative to current directory.
        Path is scaped to support spaces.
        Example:
            if a file is in /home/username/documents/test/import_file.mrc
            and current directory is /home/username/documents
            this will return '/test/import_file.mrc'
        """
        return self.scapePath(subtomogram if isinstance(subtomogram, str) else subtomogram.getFileName())

    def getSubtomogramAbsolutePath(self, subtomogram: Union[SubTomogram, Volume, str]) -> str:
        """
        This method returns a the absolute path for the subtomogram.
        Path is scaped to support spaces.
        Example: '/home/username/documents/test/import_file.mrc'
        """
        return self.scapePath(os.path.abspath(self.getSubtomogramRelativePath(subtomogram)))

    def getSubtomogramName(self, filename: str) -> str:
        """
        This method returns the full name of the given subtomogram files.
        Example: import_file.mrc
        """
        return os.path.basename(filename)
    
    def getCleanSubtomogramName(self, filename: str) -> str:
        """
        This method returns the full name of the given subtomogram file without the 'import_' prefix.
        Example:
            if filename is 'import_file.mrc'
            this will return 'file.mrc'
        """
        return self.getSubtomogramName(filename).replace('import_', '')

    def getProjectionName(self, subtomogram: Union[SubTomogram, Volume, str]) -> str:
        """
        This function returns the name of the projection file for a given input subtomogram.
        """
        name = os.path.splitext(self.getCleanSubtomogramName(self.getSubtomogramAbsolutePath(subtomogram)))[0]
        return '{}_image.mrcs'.format(name)
    
    def getProjectionMetadataName(self, subtomogram: Union[SubTomogram, Volume, str]) -> str:
        """
        This function returns the filename of the metadata file for a given input subtomogram.
        """
        return os.path.splitext(self.getProjectionName(subtomogram))[0] + '.xmd'

    def getProjectionAbsolutePath(self, subtomogram: Union[SubTomogram, Volume, str]) -> str:
        """
        This function returns the full path of a given subtomogram.
        """
        return self.scapePath(os.path.abspath(self._getExtraPath(self.getProjectionName(subtomogram))))
    
    def getProjectionMetadataAbsolutePath(self, subtomogram: Union[SubTomogram, Volume, str]) -> str:
        """
        This function returns the full path of a given subtomogram's metadata file.
        """
        return os.path.splitext(self.getProjectionAbsolutePath(subtomogram))[0] + '.xmd'
    
    def getAngleFileAbsolutePath(self, tsId: str) -> str:
        """
        This function returns the full path of a given tilt series's metadata angle file.
        """
        return os.path.abspath(self._getExtraPath(f'{tsId}.xmd'))

    def getStepValue(self) -> float:
        """
        This function translates the provided sample generation input to number of samples for Xmipp phantom project.
        """
        if self.tiltTypeGeneration.get() == self.TYPE_N_SAMPLES:
            return self.tiltRangeNSamples.get()
        else:
            # Converting step to number of samples
            return ((self.tiltRangeStart.get() % 360) - (self.tiltRangeEnd.get() % 360)) / self.tiltRangeStep.get()
    
    def getMethodValue(self) -> str:
        """
        This function returns the string value associated to the form value provided by the user regarding transform method.
        """
        if self.transformMethod.get() == self.METHOD_FOURIER:
            # Obtaining string and params for fourier method
            methodString = f'fourier {self.pad.get()} {self.maxfreq.get()} '

            # Adding interpolation method
            if self.interp.get() == self.INTERPOLATION_BSPLINE:
                methodString += 'bspline'
            elif self.interp.get() == self.INTERPOLATION_NEAREST:
                methodString += 'nearest'
            else:
                methodString += 'linear'
            
            # Returning complete fourier params
            return methodString
        elif self.transformMethod.get() == self.METHOD_REAL_SPACE:
            return 'real_space'
        else:
            return 'shears'

    def getSubtomogramDimensions(self) -> str:
        """
        This function retuns the first two dimensions of the subtomograms.
        """
        try:
            dimensions = self.inputSubtomograms.get().getFirstItem().getDimensions()
            return '{} {}'.format(dimensions[0], dimensions[1])
        except TypeError:
            errorMessage = " ".join(["No subtomograms were received. Check the output of the previous protocol.",
                                     "If you are using the integrated test, run the extract subtomos's test first."])
            raise TypeError(errorMessage)

    def getXmippParamPath(self, tsId: str='') -> str:
        """
        This function returns the path for the config file for Xmipp Phantom.
        """
        filename = 'projections' + (f'_{tsId}' if tsId else '') + '.param'
        return self.scapePath(os.path.abspath(os.path.join(self._getExtraPath(''), filename)))
    
    def getReferenceMetadataFile(self) -> str:
        """
        This function returns the path for the reference metadata file used when generating projections
        based on a Tilt Series.
        """
        return os.path.join(os.path.dirname(self.getXmippParamPath()), 'reference.xmd')
    
    def getSummary(self, setOfParticles: SetOfParticles) -> str:
        """
        Returns the summary of a given set of particles.
        The summary consists of a text including the number of particles in the set.
        """
        return "Number of projections generated: {}".format(setOfParticles.getSize())

    def getCorrectedDefocus(self, tiltAngle: Float, coordinates: Coordinate3D) -> Tuple[Float, Float]:
        """
        This function returns the corrected defocusU and defocusV given a tilt angle and
        the coordinates of a subtomogram.
        """
        # Converting tilt angle to radians (received in degrees)
        radiansTiltAngle = math.radians(tiltAngle)
        
        # Calculating defocus direction
        defocusDir = -1 if self.defocusDir.get() else 1

        for ctf in self.inputCTF.get():
            # From the input set of CTFs, get the Tilt Series of each CTF
            ts = ctf.getTiltSeries()

            # If the Tilt Series id matches the id of the tomogram where the
            # subtomogram containing the given coordinates comes from, correct
            # the defocus of the closest CTF to tilt angle from current Tilt Series
            if (ts.getTsId() == coordinates.getTomoId()):
                # Obtaining closest CTF to tilt angle from current Tilt Series
                closestCTF = self.getClosestCTF(ts, tiltAngle)

                # Obtain CTF's defocus
                defocusU, defocusV = closestCTF.getDefocusU(), closestCTF.getDefocusV()
                
                # Obtain and return corrected defocus
                generalDefocus = (coordinates.getX() * math.cos(radiansTiltAngle) + coordinates.getZ() * math.sin(radiansTiltAngle)) * ts.getSamplingRate() * math.sin(radiansTiltAngle)
                correctedDefU = defocusU + defocusDir * generalDefocus
                correctedDefV = defocusV + defocusDir * generalDefocus
                return correctedDefU, correctedDefV

    def getClosestCTF(self, tiltSeries: TiltSeries, tiltAngle: Float) -> CTFModel:
        """
        This function returns the Tilt Image's CTF inside a given
        Tilt Series which is closest to the given tilt angle.
        """
        # Initial distance
        distance = 360

        # Find , 
        for tiltImage in tiltSeries:
            # Obtaining difference in tilt between given angle and tilt image angle
            tiltDistance = abs(tiltAngle - tiltImage.getTiltAngle())

            # If tilt difference is less than current distance, update distance and
            # keep current image's CTF as current closest CTF
            if tiltDistance < distance:
                distance = tiltDistance
                outputCTF = tiltImage.getCTF()
        
        # Returning closest CTF
        return outputCTF
    
    def getAngleDictionary(self) -> TypedDict:
        """
        This function returs a dictionary containing all the angles of each Tilt Series of the input set
        """
        # Initializing empty dictionary
        angleDict = {}
        
        if self.tiltTypeGeneration.get() == self.TYPE_TILT_SERIES:
            for ts in self.tiltRangeTS.get():
                tsId = ts.getTsId()
                # Adding angle list to dictionary
                angleDict[tsId] = []
                for ti in ts:
                    # Add each angle tilt and rotation to the list
                    angleDict[tsId].append({self.COLUMN_ANGLE_ROT: calculateRotationAngleAndShiftsFromTM(ti)[0], self.COLUMN_ANGLE_TILT: ti.getTiltAngle()})

        # Returning result dictionary
        return angleDict
