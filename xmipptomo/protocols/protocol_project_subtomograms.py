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
import os

from pwem.protocols import EMProtocol
from pwem.objects import SetOfParticles

from pyworkflow import BETA
from pyworkflow.protocol import params

from tomo.protocols import ProtTomoBase
from xmipp3.convert import readSetOfParticles

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
    TYPE_N_SAMPLES = 0
    TYPE_STEP = 1

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('inputSubtomograms', params.PointerParam, pointerClass="SetOfSubTomograms,SetOfVolumes",
                      label='Set of subtomograms', help="Set of subtomograms to be carried into SPA")
        form.addParam('cleanTmps', params.BooleanParam, default='True', label='Clean temporary files: ', expertLevel=params.LEVEL_ADVANCED,
                        help='Clean temporary files after finishing the execution.\nThis is useful to reduce unnecessary disk usage.')
        form.addParam('transformMethod', params.EnumParam, display=params.EnumParam.DISPLAY_COMBO, default=self.METHOD_FOURIER,
                        choices=['Fourier', 'Real space'], label="Transform method: ", expertLevel=params.LEVEL_ADVANCED,
                        help='Select the algorithm that will be used to obtain the projections.')
        tiltGroup = form.addGroup('Tilt range')
        tiltGroup.addParam('tiltRangeStart', params.IntParam, default=-60, label='Tilt range start:',
                        help='The initial value of the range of angles the projection will be produced on.\nDefaults to -60º.')
        tiltGroup.addParam('tiltRangeEnd', params.IntParam, default=60, label='Tilt range end:',
                        help='The final value of the range of angles the projection will be produced on.\nDefaults to 60º.')
        tiltGroup.addParam('tiltTypeGeneration', params.EnumParam, display=params.EnumParam.DISPLAY_COMBO, default=self.TYPE_N_SAMPLES,
                        choices=['NSamples', 'Step'], label="Type of sample generation: ",
                        help='Select either the number of samples to be taken or the separation in degrees between each sample.')
        tiltGroup.addParam('tiltRangeNSamples', params.IntParam, condition='tiltTypeGeneration==%d' % self.TYPE_N_SAMPLES, label='Number of samples:',
                        help='Number of samples to be produced.\nIt has to be 1 or greater.')
        tiltGroup.addParam('tiltRangeStep', params.IntParam, condition='tiltTypeGeneration==%d' % self.TYPE_STEP, label='Step:',
                        help='Number of degrees each sample will be separated from the next.\nIt has to be greater than 0.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        # Writing config file
        self._insertFunctionStep('generateConfigFile')

        # Generating projections for each subtomogram
        self._insertFunctionStep('generateSubtomogramProjections')
        
        # Conditionally removing temporary files
        if self.cleanTmps.get():
            self._insertFunctionStep('removeTempFiles')
        
        # Create output
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def generateConfigFile(self):
        """
        This function writes the config file for Xmipp Phantom.
        """
        confFile = open(self.getXmippParamPath(), "w")

        # Generating file content
        content = '# XMIPP_STAR_1 *\n'
        content += '# Projection Parameters\n'
        content += 'data_block1\n'
        content += '# X and Y projection dimensions [Xdim Ydim]\n'
        content += '_dimensions2D   \'{}\'\n'.format(self.getSubtomogramDimensions())
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
        confFile.write(content)
        confFile.close()

    def generateSubtomogramProjections(self):
        """
        This function generates the projection for all input subtomograms.
        """
        for subtomogram in self.inputSubtomograms.get():
            params = '-i {}'.format(self.getSubtomogramAbsolutePath(subtomogram))   # Path to subtomogram
            params += ' -o {}'.format(self.getProjectionAbsolutePath(subtomogram))  # Path to output projection
            params += ' --method {}'.format(self.getMethodValue())                  # Projection algorithm
            params += ' --params {}'.format(self.getXmippParamPath())               # Path to Xmipp phantom param file
            self.runJob("xmipp_phantom_project", params, cwd='/home')

    def removeTempFiles(self):
        """
        This function removes the temporary files of this protocol.
        """
        # Removing Xmipp Phantom config file
        self.runJob('rm', self.getXmippParamPath())

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

        # Adding projections of each subtomogram as a particle each
        for subtomogram in inputSubtomograms.iterItems():
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
        if self.tiltTypeGeneration.get() == self.TYPE_N_SAMPLES and (not self.tiltRangeNSamples.get() == None) and self.tiltRangeNSamples.get() < 1:
            errors.append('The number of samples cannot be less than 1.')
        
        # Checking if the step is greater than 0
        if self.tiltTypeGeneration.get() == self.TYPE_STEP and (not self.tiltRangeStep.get() == None) and self.tiltRangeStep.get() <= 0:
            errors.append('The step must be greater than 0.')

        return errors

    def _summary(self):
        """
        This method usually returns a summary of the text provided by '_methods'.
        """
        return []

    def _methods(self):
        """
        This method returns a text intended to be copied and pasted in the paper section 'materials & methods'.
        """
        return []

    # --------------------------- UTILS functions --------------------------------------------
    def scapePath(self, path):
        """
        This function returns the given path with all the spaces in folder names scaped to avoid errors.
        """
        # os.path.baspath adds '\\' when finding a foldername with '\ ', so '\\\' needs to be replaced with ''
        # Then, '\' is inserted before every space again, to include now possible folders with spaces in the absolute path
        return path.replace('\\\ ', ' ').replace(' ', '\ ')
    
    def getSubtomogramRelativePath(self, subtomogram):
        """
        This method returns a the subtomogram path relative to current directory.
        Path is scaped to support spaces.
        Example:
            if a file is in /home/username/documents/test/import_file.mrc
            and current directory is /home/username/documents
            this will return '/test/import_file.mrc'
        """
        return self.scapePath(subtomogram.getFileName())

    def getSubtomogramAbsolutePath(self, subtomogram):
        """
        This method returns a the absolute path for the subtomogram.
        Path is scaped to support spaces.
        Example: '/home/username/documents/test/import_file.mrc'
        """
        return self.scapePath(os.path.abspath(self.getSubtomogramRelativePath(subtomogram)))

    def getSubtomogramName(self, filename):
        """
        This method returns the full name of the given subtomogram files.
        Example: import_file.mrc
        """
        return os.path.basename(filename)
    
    def getCleanSubtomogramName(self, filename):
        """
        This method returns the full name of the given subtomogram file without the 'import_' prefix.
        Example:
            if filename is 'import_file.mrc'
            this will return 'file.mrc'
        """
        return self.getSubtomogramName(filename).replace('import_', '')

    def getProjectionName(self, subtomogram):
        """
        This function returns the name of the projection file for a given input subtomogram.
        """
        name = os.path.splitext(self.getCleanSubtomogramName(self.getSubtomogramAbsolutePath(subtomogram)))[0]
        return '{}_image.mrc'.format(name)
    
    def getProjectionMetadataName(self, subtomogram):
        """
        This function returns the filename of the metadata file for a given input subtomogram.
        """
        return self.getProjectionName(subtomogram).replace('.mrc', '.xmd')

    def getProjectionAbsolutePath(self, subtomogram):
        """
        This function returns the full path of a given subtomogram.
        """
        return self.scapePath(os.path.abspath(self._getExtraPath(self.getProjectionName(subtomogram))))
    
    def getProjectionMetadataAbsolutePath(self, subtomogram):
        """
        This function returns the full path of a given subtomogram's metadata file.
        """
        return self.getProjectionAbsolutePath(subtomogram).replace('.mrc', '.xmd')

    def getStepValue(self):
        """
        This function translates the provided sample generation input to number of samples for Xmipp phantom.
        """
        if self.tiltTypeGeneration.get() == self.TYPE_N_SAMPLES:
            return self.tiltRangeNSamples.get()
        else:
            # Converting step to number of samples
            return ((self.tiltRangeStart.get() % 360) - (self.tiltRangeEnd.get() % 360)) / self.tiltRangeStep.get()
    
    def getMethodValue(self):
        """
        This function returns the string value associated to the form value provided by the user regarding transform method.
        """
        if self.transformMethod.get() == self.METHOD_FOURIER:
            return 'fourier'
        else:
            return 'real_space'

    def getSubtomogramDimensions(self):
        """
        This function retuns the first two dimensions of the subtomograms.
        """
        dimensions = self.inputSubtomograms.get().getFirstItem().getDimensions()
        return '{} {}'.format(dimensions[0], dimensions[1])

    def getXmippParamPath(self):
        """
        This function returns the path for the config file for Xmipp Phantom.
        """
        return self.scapePath(os.path.abspath(os.path.join(self._getExtraPath(''), 'xmippPhantom.param')))
    
    def getSummary(self, setOfParticles):
        """
        Returns the summary of a given set of particles.
        The summary consists of a text including the number of particles in the set.
        """
        return "Number of projections generated: {}".format(setOfParticles.getSize())
