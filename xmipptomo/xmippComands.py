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
from pyworkflow.object import Integer, Float, String, Pointer, Boolean, CsvList
from tomo.objects import TiltImage, SetOfTiltSeries, TiltSeries, TomoAcquisition, Transform, SubTomogram


def xmipp_phantom_project(inputdata, outputdata, method, projparams):

    params = '-i %s ' % inputdata  # Path to subtomogram
    params += ' -o %s ' % outputdata  # Path to output projection
    params += ' --method %s ' % method  # Projection algorithm
    params += ' --params %s ' % projparams  # Path to Xmipp phantom param file

    return params

def generateParamsFileProjection(paramsFn, dimX, dimY, tiltRangeStart, tiltRangeEnd, tiltSamples):
    """
    This function writes the config file for Xmipp Phantom.
    """
    confFile = open(paramsFn, "w")

    # Generating file content
    content = '# XMIPP_STAR_1 *\n'
    content += '# Projection Parameters\n'
    content += 'data_block1\n'
    content += '# X and Y projection dimensions [Xdim Ydim]\n'
    content += '_dimensions2D %i %i \n' % (dimX, dimY)
    content += '# Rotation range and number of samples [Start Finish NSamples]\n'
    content += '_projRotRange    \'0 0 1\'\n'
    content += '# Rotation angle added noise  [noise (bias)]\n'
    content += '_projRotNoise   \'0\'\n'
    content += '# Tilt range and number of samples for Tilt [Start Finish NSamples]\n'
    content += '_projTiltRange    \'%f %f %f\'\n' % (tiltRangeStart, tiltRangeEnd, tiltSamples)
    content += '# Tilt angle added noise\n'
    content += '_projTiltNoise   \'0\'\n'
    content += '# Psi range and number of samples\n'
    content += '_projPsiRange    \'0 0 0\'\n'
    content += '# Psi added noise\n_projPsiNoise   \'0\'\n'
    content += '# Noise\n'

    # Writing content to file and closing
    confFile.write(content)
    confFile.close()

