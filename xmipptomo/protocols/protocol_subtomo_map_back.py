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

from os.path import basename

from pwem import ALIGN_3D
from pwem.emlib import lib
import pwem.emlib.metadata as md
from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, FloatParam

from tomo.objects import Tomogram
from tomo.protocols import ProtTomoBase
from xmipp3.convert import alignmentToRow


class XmippProtSubtomoMapBack(EMProtocol, ProtTomoBase):
    """ This protocol takes a tomogram, a reference subtomogram and a metadata with geometrical parameters
   (x,y,z) and places the reference subtomogram on the tomogram at the designated locations (map back).
   It has different representation options."""

    _label = 'map back subtomos'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('selection', EnumParam,
                      choices=['Class', 'Subtomograms'], default=0, display=EnumParam.DISPLAY_HLIST,
                      label='Select input type',
                      help='Class: input is a class generated previously in Scipion.\nSubtomograms: input is a set of'
                           ' subtomograms previously aligned.')
        form.addParam('inputClasses', PointerParam, pointerClass="SetOfClassesSubTomograms", label='Class',
                      condition='selection==0', allowsNull=True,
                      help="Subtomogram class from which the coordinates of the subtomograms and the reference will be "
                           "used. It should be a SetOfClassesSubTomograms with just 1 item.")
        form.addParam('inputSubtomos', PointerParam, pointerClass="SetOfSubTomograms", label='Subtomograms',
                      condition='selection==1', allowsNull=True,
                      help="Subtomograms to be mapped back, they should have alignment and coordinates.")
        form.addParam('inputRef', PointerParam, pointerClass="Volume, SubTomogram, AverageSubTomogram",
                      label='Reference', condition='selection==1', allowsNull=True,
                      help="Subtomogram reference, average, representative or initial model of the subtomograms.")
        form.addParam('inputTomograms', PointerParam, pointerClass="SetOfTomograms",
                      label='Original tomograms', help="Original tomograms from which the subtomograms were extracted")
        form.addParam('invertContrast', BooleanParam, default=False, label='Invert reference contrast',
                      help="Invert the contrast if the reference is black over a white background.  Xmipp, Spider, "
                           "Relion and Eman require white particles over a black background. ")
        form.addParam('paintingType', EnumParam,
                      choices=['Copy', 'Average', 'Highlight', 'Binarize'],
                      default=0, important=True,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Painting mode',
                      help='The program has several painting options:\n*Copy*: Copying the reference onto the tomogram.'
                           '\n*Average*: Setting the region occupied by the reference in the tomogram to the average '
                           'value of that region.\n*Highlight*: Add the reference multiplied by a constant to the '
                           'location specified.\n*Binarize*: Copy a binarized version of the reference onto the '
                           'tomogram.')
        form.addParam('removeBackground', BooleanParam, default=False, label='Remove background',
                      help="Set tomogram to 0", condition="paintingType == 0 or paintingType == 3")
        form.addParam('threshold', FloatParam, default=0.5, label='Threshold',
                      help="threshold applied to tomogram", condition="paintingType == 1 or paintingType == 3")
        form.addParam('constant', FloatParam, default=2, label='Multiplier',
                      help="constant to multiply the reference",
                      condition="paintingType == 2")

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInput')
        if self.selection == 0:
            for subtomoClass in self.inputClasses.get():
                self._insertFunctionStep('runMapBack', subtomoClass.getObjId())
        else:
            self._insertFunctionStep('runMapBack', 0)
        self._insertFunctionStep('createOutput')

    # --------------------------- STEPS functions -------------------------------
    def convertInput(self):
        for tomo in self.inputTomograms.get().iterItems():
            img = ImageHandler()
            fnTomo = self._getExtraPath('tomogram_%d.mrc' % tomo.getObjId())
            img.convert(tomo, fnTomo)
            if self.selection == 0:
                for classSubt in self.inputClasses.get().iterItems():
                    cId = classSubt.getFirstItem().getClassId()
                    fnRef = self._getExtraPath('reference%d.mrc' % cId)
                    img.convert(classSubt.getRepresentative(), fnRef)
                    if self.invertContrast.get() == True:
                        self.runJob("xmipp_image_operate", " -i %s  --mult -1" % fnRef)
            else:
                if self.invertContrast.get() == True:
                    self.runJob("xmipp_image_operate", " -i %s  --mult -1" % self.inputRef.get().getFileName())

            if self.paintingType.get() == 0 or self.paintingType.get() == 3:
                if self.removeBackground.get() == True:
                    self.runJob("xmipp_image_operate", " -i %s  --mult 0" % fnTomo)

    def runMapBack(self, classId):
        for tomo in self.inputTomograms.get().iterItems():
            if self.selection == 0:
                TsSubtomo = self.inputClasses.get().getSamplingRate()
            else:
                TsSubtomo = self.inputRef.get().getSamplingRate()
            TsTomo = tomo.getSamplingRate()
            scaleFactor = TsSubtomo/TsTomo
            mdGeometry = lib.MetaData()

            if self.selection == 0:
                inputSet = self.inputClasses.get().getFirstItem()
                ref = self._getExtraPath("reference%d.mrc" % classId)
            else:
                inputSet = self.inputSubtomos.get()
                ref = self.inputRef.get().getFileName()

            for subtomo in inputSet.iterItems():
                if subtomo.getCoordinate3D().getVolId() == tomo.getObjId() \
                        or basename(subtomo.getVolName()) == tomo.getBaseName().partition('_')[2]:
                    nRow = md.Row()
                    nRow.setValue(lib.MDL_ITEM_ID, int(subtomo.getObjId()))
                    nRow.setValue(lib.MDL_XCOOR, int(subtomo.getCoordinate3D().getX()*scaleFactor))
                    nRow.setValue(lib.MDL_YCOOR, int(subtomo.getCoordinate3D().getY()*scaleFactor))
                    nRow.setValue(lib.MDL_ZCOOR, int(subtomo.getCoordinate3D().getZ()*scaleFactor))
                    # Convert transform matrix to Euler Angles (rot, tilt, psi)
                    alignmentToRow(subtomo.getTransform(), nRow, ALIGN_3D)
                    nRow.addToMd(mdGeometry)
            fnGeometry = self._getExtraPath("geometry%d.xmd" % classId)
            mdGeometry.write(fnGeometry)

            if TsSubtomo != TsTomo:
                factor = TsSubtomo/TsTomo
                args = "-i %s -o %s --scale %d" % (ref, ref, factor)
                self.runJob('xmipp_transform_geometry', args)

            if self.paintingType.get() == 0:
                painting = 'copy'
            elif self.paintingType.get() == 1:
                painting = 'avg %d' % self.threshold.get()
            elif self.paintingType.get() == 2:
                painting = 'highlight %d' % self.constant.get()
            elif self.paintingType.get() == 3:
                painting = 'copy_binary %f' % self.threshold.get()

            tomogram = self._getExtraPath("tomogram_%d.mrc" % tomo.getObjId())
            args = " -i %s -o %s --geom %s --ref %s --method %s" % (tomogram, tomogram,
                                                                    self._getExtraPath("geometry%d.xmd" % classId),
                                                                    ref, painting)
            self.runJob("xmipp_tomo_map_back", args)

    def createOutput(self):
        inputTomos = self.inputTomograms.get()
        outputTomos = self._createSetOfTomograms()
        outputTomos.copyInfo(inputTomos)
        for j, _ in enumerate(inputTomos):
            tomo = Tomogram()
            tomo.setLocation(self._getExtraPath("tomogram_%d.mrc" % int(j+1)))
            outputTomos.append(tomo)
        self._defineOutputs(outputTomograms=outputTomos)
        self._defineSourceRelation(self.inputTomograms, outputTomos)
        if self.selection == 0:
            self._defineSourceRelation(self.inputClasses, outputTomos)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        if self.selection == 0:
            for subtomo in self.inputClasses.get().getFirstItem().iterItems():
                if not subtomo.hasCoordinate3D():
                    validateMsgs.append('Please provide a class which contains subtomograms with 3D coordinates.')
                    break
                if not subtomo.hasTransform():
                    validateMsgs.append('Please provide a class which contains subtomograms with alignment.')
                    break
        else:
            for subtomo in self.inputSubtomos.get().iterItems():
                if not subtomo.hasCoordinate3D():
                    validateMsgs.append('Please provide a set of subtomograms which contains subtomograms with 3D '
                                        'coordinates.')
                    break
                if not subtomo.hasTransform():
                    validateMsgs.append('Please provide a set of subtomograms which contains subtomograms with '
                                        'alignment.')
                    break
        return validateMsgs

    def _summary(self):
        summary = []
        if self.selection == 0:
            refSize = len(self.inputClasses.get())
            setSize = len(self.inputClasses.get().getFirstItem())
        else:
            refSize = 1
            setSize = len(self.inputSubtomos.get())
        summary.append("%d subtomogram reference mapped back %d times to original tomograms" % (refSize, setSize))
        return summary

    def _methods(self):
        methods = []
        if self.selection == 0:
            refSize = len(self.inputClasses.get())
            setSize = len(self.inputClasses.get().getFirstItem())
        else:
            refSize = 1
            setSize = len(self.inputSubtomos.get())
        methods.append("References from %d subtomogram classes mapped back %d times to original tomograms %s" %
                       (refSize, setSize, self.getObjectTag('inputTomograms')))
        return methods
