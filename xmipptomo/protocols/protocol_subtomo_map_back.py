# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
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
import enum
import os

from pwem.emlib import lib
import pwem.emlib.metadata as md
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pwem.protocols import EMProtocol

from pyworkflow import BETA
from pyworkflow.protocol import STEPS_PARALLEL
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, FloatParam
from pyworkflow.utils import createLink

from tomo.objects import Tomogram, SetOfCoordinates3D, SetOfSubTomograms, SetOfClassesSubTomograms, SubTomogram, \
    MATRIX_CONVERSION, SetOfTomograms
from tomo.protocols import ProtTomoBase
from xmipp3.convert import alignmentToRow
import tomo.constants as const

REFERENCE = 'Reference'


# Painting types
class PAINTING_TYPES:
    COPY=0
    AVERAGE=1
    HIGHLIGHT=2
    BINARIZE=3


class MapBackOutputs(enum.Enum):
    tomograms = SetOfTomograms

class XmippProtSubtomoMapBack(EMProtocol, ProtTomoBase):
    """ This protocol takes a tomogram, a reference subtomogram and a metadata with geometrical parameters
   (x,y,z) and places the reference subtomogram on the tomogram at the designated locations (map back).
   It has different representation options."""

    _label = 'map back subtomos'
    _devStatus = BETA
    _possibleOutputs = MapBackOutputs

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
        self.tomos = None

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('selection', EnumParam,
                      choices=['Class', 'Subtomograms'], default=0, display=EnumParam.DISPLAY_HLIST,
                      label='Select input type',
                      help='Class: input is a class generated previously in Scipion.\nSubtomograms: input is a set of'
                           ' subtomograms previously aligned.')
        form.addParam('inputClasses', PointerParam, pointerClass=SetOfClassesSubTomograms, label='Class',
                      condition='selection==0', allowsNull=True,
                      help="Subtomogram class from which the coordinates of the subtomograms and the reference will be "
                           "used. It should be a SetOfClassesSubTomograms with just 1 item.")
        form.addParam('inputSubtomos', PointerParam, pointerClass=[SetOfSubTomograms, SetOfCoordinates3D], label='Subtomograms/coordinates',
                      condition='selection==1', allowsNull=True,
                      help="Subtomograms to be mapped back, they should have alignment and coordinates.")
        form.addParam('inputRef', PointerParam, pointerClass="Volume, SubTomogram, AverageSubTomogram",
                      label=REFERENCE, condition='selection==1', allowsNull=True,
                      help="Subtomogram reference, average, representative or initial model of the subtomograms.")
        form.addParam('inputTomograms', PointerParam, pointerClass="SetOfTomograms",
                      label='Original tomograms', help="Original tomograms from which the subtomograms were extracted", allowsNull=True)
        form.addParam('invertContrast', BooleanParam, default=False, label='Invert reference contrast',
                      help="Invert the contrast if the reference is black over a white background.  Xmipp, Spider, "
                           "Relion and Eman require white particles over a black background. ")
        form.addParam('paintingType', EnumParam,
                      choices=['Copy', 'Average', 'Highlight', 'Binarize'],
                      default=PAINTING_TYPES.HIGHLIGHT, important=True,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Painting mode',
                      help='The program has several painting options:\n*Copy*: Copying the reference onto the tomogram.'
                           '\n*Average*: Setting the region occupied by the reference in the tomogram to the average '
                           'value of that region.\n*Highlight*: Add the reference multiplied by a constant to the '
                           'location specified.\n*Binarize*: Copy a binarized version of the reference onto the '
                           'tomogram.')
        form.addParam('removeBackground', BooleanParam, default=False, label='Remove background',
                      help="Set tomogram to 0", condition="paintingType == %s or paintingType == %s" % (PAINTING_TYPES.COPY, PAINTING_TYPES.BINARIZE))
        form.addParam('threshold', FloatParam, default=0.5, label='Threshold',
                      help="threshold applied to tomogram", condition="paintingType == %s or paintingType == %s" % (PAINTING_TYPES.AVERAGE, PAINTING_TYPES.BINARIZE))
        form.addParam('constant', FloatParam, default=2, label='Multiplier',
                      help="constant to multiply the reference", condition="paintingType == %s or paintingType == %s" % (PAINTING_TYPES.COPY, PAINTING_TYPES.HIGHLIGHT))

        form.addParallelSection(threads=4, mpi=1)
    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        prepRefId = self._insertFunctionStep(self.prepareReference, self.invertContrast.get())


        mapBacksStepIds=[]
        # For each tomogram
        for key, value in self._getTomogramsInvolved().items():
            mapBackStepId = self._insertFunctionStep(self.runMapBack, value.getTsId(),
                                     self.paintingType.get(),
                                     self.removeBackground.get(),
                                     self.threshold.get(),
                                     prerequisites=prepRefId)

            mapBacksStepIds.append(mapBackStepId)

        self._insertFunctionStep(self.createOutput, prerequisites=mapBacksStepIds)

    # --------------------------- STEPS functions -------------------------------


    def prepareReference(self, invertContrast):
        """ Prepares the reference for the mapback"""
        fnRef = self.getFinalRefName()
        sourceRef = self.getSourceReferenceFn()

        # Do we need this conversion!!
        # img = ImageHandler()
        # img.convert(sourceRef, fnRef)

        refSampling = self.inputRef.get().getSamplingRate()
        tomoSampling = self.getInputSetOfTomograms().getSamplingRate()
        # for xmipp 0.5 means halving, 2 means doubling
        factor = refSampling/tomoSampling
        if invertContrast:
            self.runJob("xmipp_image_operate", " -i %s  --mult -1 -o %s" % (sourceRef, fnRef))
            sourceRef = fnRef

        if factor != 1:
            self.runJob("xmipp_image_resize", " -i %s  --factor %0.2f -o %s" % (sourceRef, factor, fnRef))


        if not os.path.exists(fnRef):
            createLink(sourceRef, fnRef)

    def getSourceReferenceFn(self):
        """ Returns the source reference file name: representative from the first class or the reference param"""
        if self._useClasses():
            firstClass = self.inputClasses.get().getFirstItem()
            return firstClass.getRepresentative()
        else:
            return self.inputRef.get().getFileName()

    def getFinalRefName(self):
        """ returns the final path of the reference"""
        return self._getExtraPath('reference.mrc')

    def getFinalTomoName(self, tomo):
        """ Returns the final tomogram name having a tomogram. Uses the Tilt Series Id"""

        return self._getExtraPath('tomogram_%s.mrc' % tomo.getTsId())

    def removeTomogramBackground(self, tomo):
        """ Converts the tomogram into a black box by multiplying al voxels by 0"""
        img = ImageHandler()
        fnTomo = self.getFinalTomoName(tomo)
        img.convert(tomo, fnTomo)

        if self.paintingType.get() == PAINTING_TYPES.COPY or \
                self.paintingType.get() == PAINTING_TYPES.BINARIZE:
            if self.removeBackground.get() == True:
                self.runJob("xmipp_image_operate", " -i %s  --mult 0" % fnTomo)

    def _getTomogramsInvolved(self):

        """ Returns only the tomograms involved in the input set"""
        if self.tomos is None:

            input = self._getInput()

            if isinstance(input, SetOfSubTomograms):
                self.tomos = input.getTomograms()
            else:
                self.tomos = input.getPrecedentsInvolved()

        return self.tomos
    def useOtherSetOfTomograms(self):
        """ Returns true if tomograms to be used are those in inputTomograms"""
        return self.inputTomograms.get() is not None

    def getInputSetOfTomograms(self):
        """ Returns the set of tomograms from the input"""

        if self.useOtherSetOfTomograms():
            return self.inputTomograms.get()
        else:
            input = self._getInput()

            if isinstance(input, SetOfSubTomograms):
                input = input.getCoordinates3D()

            return input.getPrecedents()

    def _getInput(self):
        """ Returns the iterator on the input set that could be on classes, subtomograms or coordinates"""
        if self._useClasses():
            return self.inputClasses.get().getFirstItem()
        else:
            return self.inputSubtomos.get()

    def _useClasses(self):
        """ Returns true if inputClasses attribute has to be used as input"""
        return self.selection == 0

    def _isInputA3DClass(self):
        """ Returns true if input is a 3D class"""
        return isinstance(self.inputClasses.get(), SetOfClassesSubTomograms)

    def _isInputASetOfSubtomograms(self):
        """ Returns true if the input is a set of subtomograms"""
        return isinstance(self.inputSubtomos.get(), SetOfSubTomograms)

    def _isInputASetOfCoordinates(self):
        """ Returns true if the input is a set of coordinates"""
        return isinstance(self.inputSubtomos.get(), SetOfCoordinates3D)

    def getTomogram(self, tsId):
        """ Returns a tomogram object based on its tsId"""
        return self._getTomogramsInvolved()[tsId]

    def runMapBack(self, tsId, paintingType, removeBackground, threshold):

        tomo = self.getTomogram(tsId)
        self.debug("TomoInvolved: %s" % tomo)

        # Removing background
        self.removeTomogramBackground(tomo)

        input = self._getInput()

        inputSR = input.getSamplingRate()
        tomoSR = tomo.getSamplingRate()
        scaleFactor = inputSR/tomoSR
        mdGeometry = lib.MetaData()

        ref = self.getFinalRefName()

        if self.paintingType.get() == PAINTING_TYPES.COPY and self.constant.get() != 1:
            initialref = self.inputRef.get().getFileName()
            if initialref.endswith('.mrc'):
                initialref += ':mrc'
            ref = self._getExtraPath("ref_mult.mrc")
            self.runJob("xmipp_image_operate", " -i %s  --mult %d -o %s" %
                        (initialref, self.constant.get(), ref))

        # Using subtomograms
        usingSubtomograms = isinstance(input, SetOfSubTomograms)

        where = "_coordinate._tomoId='%s'" if usingSubtomograms else "_tomoId='%s'"
        where = where % tsId

        for item in input.iterItems(where=where):
            self.debug("Mapping back %s" % item)
            if usingSubtomograms:
                coord = item.getCoordinate3D()
                # A coordinte does not have an objId in this case, we set it
                coord.setObjId(item.getObjId())
                transform = item.getTransform(convention=MATRIX_CONVERSION.XMIPP)
            else: # Coordinate
                coord = item
                transform = Transform(matrix=item.getMatrix(convention=MATRIX_CONVERSION.XMIPP))

            if coord.getVolId() == tomo.getObjId() \
                    or coord.getTomoId() == tomo.getTsId():
                nRow = md.Row()
                nRow.setValue(lib.MDL_ITEM_ID, int(coord.getObjId()))
                coord.setVolume(tomo)
                nRow.setValue(lib.MDL_XCOOR, int(coord.getX(const.BOTTOM_LEFT_CORNER)*scaleFactor))
                nRow.setValue(lib.MDL_YCOOR, int(coord.getY(const.BOTTOM_LEFT_CORNER)*scaleFactor))
                nRow.setValue(lib.MDL_ZCOOR, int(coord.getZ(const.BOTTOM_LEFT_CORNER)*scaleFactor))
                # Compute inverse matrix
                #A = subtomo.getTransform().getMatrix()
                #subtomo.getTransform().setMatrix(np.linalg.inv(A))
                # Convert transform matrix to Euler Angles (rot, tilt, psi)
                from pwem import ALIGN_PROJ
                alignmentToRow(transform, nRow, ALIGN_PROJ)
                nRow.addToMd(mdGeometry)
        fnGeometry = self._getExtraPath("geometry%s.xmd" % tsId)
        mdGeometry.write(fnGeometry)

        if scaleFactor != 1:
            args = "-i %s -o %s --scale %d" % (ref, ref, scaleFactor)
            self.runJob('xmipp_transform_geometry', args)

        if self.paintingType.get() == PAINTING_TYPES.COPY:
            painting = 'copy'
        elif self.paintingType.get() == PAINTING_TYPES.AVERAGE:
            painting = 'avg %d' % self.threshold.get()
        elif self.paintingType.get() == PAINTING_TYPES.HIGHLIGHT:
            painting = 'highlight %d' % self.constant.get()
        elif self.paintingType.get() == PAINTING_TYPES.BINARIZE:
            painting = 'copy_binary %f' % self.threshold.get()

        tomogram = self.getFinalTomoName(tomo)
        args = " -i %s -o %s --geom %s --ref %s --method %s" % (tomogram, tomogram,
                                                                self._getExtraPath("geometry%s.xmd" % tsId),
                                                                ref, painting)
        self.runJob("xmipp_tomo_map_back", args)

    def createOutput(self):

        inputTomos = self._getTomogramsInvolved()
        inputTomosSet = self.getInputSetOfTomograms()

        outputTomos = self._createSetOfTomograms()
        outputTomos.copyInfo(inputTomosSet)
        for inputTomo in inputTomos.values():
            tomo = Tomogram()
            tomo.setLocation(self.getFinalTomoName(inputTomo))
            outputTomos.append(tomo)
        self._defineOutputs(**{MapBackOutputs.tomograms.name: outputTomos})
        if self._useClasses():
            self._defineSourceRelation(self.inputClasses, outputTomos)
        else:
            self._defineSourceRelation(inputTomosSet, outputTomos)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        if self._useClasses():
            subtomo = self.inputClasses.get().getFirstItem().getFirstItem()
            if not subtomo.hasCoordinate3D():
                validateMsgs.append('Please provide a class which contains subtomograms with 3D coordinates.')

            if not subtomo.hasTransform():
                validateMsgs.append('Please provide a class which contains subtomograms with alignment.')

        else:
            if not self.inputRef.get():
                validateMsgs.append("When using coordinates or subtomograms the %s is mandatory." % REFERENCE)

            if self._isInputASetOfSubtomograms():
                subtomo = self.inputSubtomos.get().getFirstItem()
                if not subtomo.hasCoordinate3D():
                    validateMsgs.append('Please provide a set of subtomograms which contains subtomograms with 3D '
                                            'coordinates.')
        return validateMsgs

    def _summary(self):
        summary = []
        if self._useClasses():
            summary.append("Using a class with its reference.")
        elif self._isInputASetOfSubtomograms():
            summary.append("Using subtomograms.")
        else:
            summary.append("Using 3D coordinates.")

        return summary

    def _methods(self):
        methods = []
        if self._useClasses():
            setSize = len(self.inputClasses.get().getFirstItem())
        else:
            setSize = len(self.inputSubtomos.get())

        if hasattr(self, MapBackOutputs.tomograms.name):
            methods.append("A reference was mapped back %d times into %s." %
                       (setSize, self.getObjectTag(getattr(self, MapBackOutputs.tomograms.name))))
        return methods
