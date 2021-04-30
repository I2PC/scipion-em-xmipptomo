# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
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


import numpy as np
import mrcfile

from pwem.emlib.image import ImageHandler

from pyworkflow.protocol import params
import pyworkflow.utils as pwutils
from pwem.objects import Transform

from tomo.objects import Coordinate3D
from tomo.protocols import ProtTomoPicking


class XmippProtAlignTomos(ProtTomoPicking):
    '''Protocol to align a series of Tomograms to a common reference defined by a
    Tomogram. The alignment is perform only on a central portion of the Tomograms to avoid
    working with such big volumes.'''

    _label = 'align tomos'

    def __init__(self, **args):
        ProtTomoPicking.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label='Coordinates to be aligned', important=True)
        form.addParam('referenceTomograms', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label='Reference Tomogram', important=True)
        form.addParam('square', params.BooleanParam,
                      default=True,
                      label="Square Subregion?",
                      help="If True, we will try to define a square region based on the "
                           "dimension provided to align the Tomograms. Otherwise, all the dimensions "
                           "have to be provided to define the geometry of the subregion.")
        form.addParam('x', params.IntParam, default=512,
                      label="X Dimension of Subregion")
        form.addParam('y', params.IntParam, default=512,
                      condition="square==False",
                      label="Y Dimension of Subregion")
        form.addParam('z', params.IntParam, default=512,
                      condition="square==False",
                      label="Z Dimension of Subregion")

    # --------------------------- INSERT steps functions ---------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('extractSubregions')
        self._insertFunctionStep('alignToReference')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ---------------------------
    def extractSubregions(self):
        # Define ImageHandler instance to read Tomograms and get input parameters
        ih = ImageHandler()
        dim_x = self.x.get()
        dim_y = dim_x if self.square.get() else self.y.get()
        dim_z = dim_x if self.square.get() else self.z.get()

        # Read input Tomograms, extract subregions and write them to tmp folder
        inputTomograms = self.inputCoordinates.get().getPrecedents()
        for tomo_input_file, tomo_reference_file in zip(inputTomograms.getFiles(),
                                                        self.referenceTomograms.get().getFiles()):
            # Extract subregion from input
            subregion_file = pwutils.removeBaseExt(tomo_input_file) + "_input_subregion.mrc"
            tomo_data = np.squeeze(ih.read(tomo_input_file).getData())
            tomo_subregion_data = self.crop_center(tomo_data, dim_x, dim_y, dim_z)
            img_subregion = ih.createImage()
            img_subregion.setData(tomo_subregion_data.astype(np.float64))
            ih.write(img_subregion, self._getTmpPath(subregion_file))

            # Extract subregion from reference
            subregion_file = pwutils.removeBaseExt(tomo_reference_file) + "_reference_subregion.mrc"
            tomo_data = np.squeeze(ih.read(tomo_reference_file).getData())
            tomo_subregion_data = self.crop_center(tomo_data, dim_x, dim_y, dim_z)
            img_subregion = ih.createImage()
            img_subregion.setData(tomo_subregion_data.astype(np.float64))
            ih.write(img_subregion, self._getTmpPath(subregion_file))

    def alignToReference(self):

        # Align Input Tomos to Reference Tomos and save Transformation Matrix
        inputTomograms = self.inputCoordinates.get().getPrecedents()
        for tomo_input, tomo_reference in zip(inputTomograms.getFiles(),
                                              self.referenceTomograms.get().getFiles()):
            reference_file = self._getTmpPath(pwutils.removeBaseExt(tomo_reference) +
                                              "_reference_subregion.mrc")
            input_file = self._getTmpPath(pwutils.removeBaseExt(tomo_input) +
                                          "_input_subregion.mrc")
            matrix_file = self._getExtraPath(pwutils.removeBaseExt(tomo_input) +
                                             "_input_align_matrix.txt")
            params = f"--i1 {reference_file} --i2 {input_file} --dontScale --local --copyGeo {matrix_file}"
            self.runJob("xmipp_volume_align", params)

    def createOutputStep(self):
        output_set = self._createSetOfCoordinates3D(self.referenceTomograms.get())
        inputCoordinates = self.inputCoordinates.get()
        inputTomograms = inputCoordinates.getPrecedents()
        referenceTomograms = self.referenceTomograms.get()
        output_set.setPrecedents(referenceTomograms)
        output_set.setBoxSize(inputCoordinates.getBoxSize())
        for tomo_input in inputTomograms.iterItems():
            tomo_reference = referenceTomograms[tomo_input.getObjId()]
            tomo_input_dim = np.asarray(tomo_input.getDimensions()) / 2
            tomo_reference_dim = np.asarray(tomo_reference.getDimensions()) / 2
            tomo_input_tr, tomo_reference_tr = self.createTransformation(t=-tomo_input_dim), \
                                               self.createTransformation(t=tomo_reference_dim)
            matrix_file = self._getExtraPath(pwutils.removeBaseExt(tomo_input.getFileName()) + "_input_align_matrix.txt")
            matrix = np.loadtxt(matrix_file).reshape(4, 4)
            for coordinate in inputCoordinates.iterCoordinates(volume=tomo_input):
                moved_Coord = Coordinate3D()
                moved_Coord.setVolume(tomo_reference)
                original_pos, original_transformation = np.asarray(coordinate.getPosition() + (1,)), \
                                                        coordinate.getMatrix()
                tr_coord = tomo_reference_tr @ matrix @ tomo_input_tr
                moved_pos, moved_transformation = np.dot(tr_coord, original_pos), np.dot(tr_coord, original_transformation)
                moved_Coord.setPosition(moved_pos[0], moved_pos[1], moved_pos[2])
                moved_Coord.setMatrix(moved_transformation)
                output_set.append(moved_Coord)
        self._defineOutputs(alignedSetOfCoordinates3D=output_set)
        self._defineSourceRelation(self.inputCoordinates, output_set)
        self._defineSourceRelation(self.referenceTomograms, output_set)

    # --------------------------- UTILS functions ---------------------------
    def crop_center(self, img, dim_x, dim_y, dim_z):
        img = np.swapaxes(img, 0, 2)
        x, y, z = img.shape
        start_x = x // 2 - (dim_x // 2)
        start_y = y // 2 - (dim_y // 2)
        start_z = z // 2 - (dim_z // 2)

        # Check if dimension lie outside of the Tomogram and correct them
        if start_x < 0:
            start_x = 0
            dim_x = x
            print("Subregion X dimension exceed Tomogram shape. We will take all X dimension "
                  "of Tomogram to define the subregion...")
        if start_y < 0:
            start_y = 0
            dim_y = y
            print("Subregion Y dimension exceed Tomogram shape. We will take all Y dimension "
                  "of Tomogram to define the subregion...")
        if start_z < 0:
            start_z = 0
            dim_z = z
            print("Subregion Z dimension exceed Tomogram shape. We will take all Z dimension "
                  "of Tomogram to define the subregion...")

        img = img[start_x:start_x + dim_x, start_y:start_y + dim_y, start_z:start_z + dim_z]
        img = np.swapaxes(img, 0, 2)
        return img

    def createTransformation(self, rot=np.eye(3), t=np.array([0, 0, 0])):
        transformation = np.eye(4)
        transformation[:-1, :-1] = rot
        transformation[:-1, -1] = t
        return transformation

    # --------------------------- INFO functions ---------------------------
    def _summary(self):
        summary = []

        if self.getOutputsSize() >= 1:
            summary.append('A total of *%d Coordinates* have been aligned '
                           'to *Tomograms %s*' % (self.inputCoordinates.get().getSize(),
                                                  self.referenceTomograms.get().getSize()))
        else:
            summary.append('Output not ready yet')

        return summary

    def _validate(self):
        validateMsgs = []

        if self.inputCoordinates.get().getPrecedents().getSize() != self.referenceTomograms.get().getSize():
            validateMsgs.append("The number of Tomograms associated to the input coordinates and the number of "
                                "Reference Tomograms are different. This protocol assumes a direct correspondence "
                                "in the Tomograms in the two sets mentioned before to align the coordinates. Please, "
                                "provide two sets of equal size and with the aforementioned correspondence.")

        return validateMsgs
