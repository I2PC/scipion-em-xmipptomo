# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from pwem import emlib
import pwem.emlib.metadata as md
from pwem.emlib.image import ImageHandler
from pwem.constants import (NO_INDEX, ALIGN_NONE, ALIGN_PROJ, ALIGN_2D, ALIGN_3D)
import numpy as np
from collections import OrderedDict
from xmipp3.base import getLabelPythonType, iterMdRows

def prefixAttribute(attribute):
    return '_xmipp_%s' % attribute

if not getattr(emlib, "GHOST_ACTIVATED", False):
    """ Some of MDL may not exist when Ghost is activated
    """
    # This dictionary will be used to map
    # between CTFModel properties and Xmipp labels
    ACQUISITION_DICT = OrderedDict([
           ("_amplitudeContrast", emlib.MDL_CTF_Q0),
           ("_sphericalAberration", emlib.MDL_CTF_CS),
           ("_voltage", emlib.MDL_CTF_VOLTAGE)
           ])

    COOR_DICT = OrderedDict([
        ("_x", emlib.MDL_XCOOR),
        ("_y", emlib.MDL_YCOOR),
        ("_z", emlib.MDL_ZCOOR)
    ])

    ALIGNMENT_DICT = OrderedDict([
           (prefixAttribute("shiftX"), emlib.MDL_SHIFT_X),
           (prefixAttribute("shiftY"), emlib.MDL_SHIFT_Y),
           (prefixAttribute("shiftZ"), emlib.MDL_SHIFT_Z),
           (prefixAttribute("flip"), emlib.MDL_FLIP),
           (prefixAttribute("anglePsi"), emlib.MDL_ANGLE_PSI),
           (prefixAttribute("angleRot"), emlib.MDL_ANGLE_ROT),
           (prefixAttribute("angleTilt"), emlib.MDL_ANGLE_TILT),
           ])

def writeSetOfSubtomograms(imgSet, filename, blockName='Particles', **kwargs):
    writeSetOfImages(imgSet, filename, particleToRow, blockName, **kwargs)

def writeSetOfImages(imgSet, filename, imgToFunc,
                     blockName='Images', **kwargs):
    """ This function will write a SetOfImages as a Xmipp metadata.
    Params:
        imgSet: the set of images to be written (particles,
        micrographs or volumes)
        filename: the filename where to write the metadata.
        rowFunc: this function can be used to setup the row before
            adding to metadata.
    """
    mdSet = emlib.MetaData()

    setOfImagesToMd(imgSet, mdSet, imgToFunc, **kwargs)
    mdSet.write('%s@%s' % (blockName, filename))

def particleToRow(part, partRow, **kwargs):
    """ Set labels values from Particle to md row. """
    imageToRow(part, partRow, emlib.MDL_IMAGE, **kwargs)
    coord = part.getCoordinate3D()
    if coord is not None:
        coordinateToRow(coord, partRow, copyId=False)

def setOfImagesToMd(imgSet, mdIn, imgToFunc, **kwargs):
    """ This function will fill Xmipp metadata from a SetOfMicrographs
    Params:
        imgSet: the set of images to be converted to metadata
        mdIn: metadata to be filled
        rowFunc: this function can be used to setup the row before
            adding to metadata.
    """

    if 'alignType' not in kwargs:
        kwargs['alignType'] = imgSet.getAlignment()

    if 'where' in kwargs:
        where = kwargs['where']
        for img in imgSet.iterItems(where=where):
            objId = mdIn.addObject()
            imgRow = md.Row()
            imgToFunc(img, imgRow, **kwargs)
            imgRow.writeToMd(mdIn, objId)
    else:
        for img in imgSet:
            objId = mdIn.addObject()
            imgRow = md.Row()
            imgToFunc(img, imgRow, **kwargs)
            imgRow.writeToMd(mdIn, objId)

def imageToRow(img, imgRow, imgLabel, **kwargs):
    # Provide a hook to be used if something is needed to be
    # done for special cases before converting image to row
    preprocessImageRow = kwargs.get('preprocessImageRow', None)
    if preprocessImageRow:
        preprocessImageRow(img, imgRow)

    setRowId(imgRow, img)  # Set the id in the metadata as MDL_ITEM_ID
    index, filename = img.getLocation()
    fn = locationToXmipp(index, filename)
    imgRow.setValue(imgLabel, fn)

    # alignment is mandatory at this point, it shoud be check
    # and detected defaults if not passed at readSetOf.. level
    alignType = kwargs.get('alignType')

    if alignType != ALIGN_NONE:
        alignmentToRow(img.getTransform(), imgRow, alignType)

    if kwargs.get('writeAcquisition', True) and img.hasAcquisition():
        acquisitionToRow(img.getAcquisition(), imgRow)

    # Write all extra labels to the row
    objectToRow(img, imgRow, {})

    # Provide a hook to be used if something is needed to be
    # done for special cases before converting image to row
    postprocessImageRow = kwargs.get('postprocessImageRow', None)
    if postprocessImageRow:
        postprocessImageRow(img, imgRow)

def coordinateToRow(coord, coordRow, copyId=True):
    """ Set labels values from Coordinate coord to md row. """
    if copyId:
        setRowId(coordRow, coord)
    objectToRow(coord, coordRow, COOR_DICT)
    if coord.getTomoId():
        coordRow.setValue(emlib.MDL_MICROGRAPH, str(coord.getTomoId()))

def setRowId(mdRow, obj, label=emlib.MDL_ITEM_ID):
    mdRow.setValue(label, int(obj.getObjId()))


def objectToRow(obj, row, attrDict, extraLabels=[]):
    """ This function will convert an EMObject into a Row.
    Params:
        obj: the EMObject instance (input)
        row: the Row instance (output)  -see emlib.metadata.utils.Row()-
        attrDict: dictionary with the map between obj attributes(keys) and
            row MDLabels in Xmipp (values).
        extraLabels: a list with extra labels that could be included
            as _xmipp_labelName
    """
    if obj.isEnabled():
        enabled = 1
    else:
        enabled = -1
    row.setValue(emlib.MDL_ENABLED, enabled)

    for attr, label in attrDict.items():
        if hasattr(obj, attr):
            valueType = getLabelPythonType(label)
            value = getattr(obj, attr).get()
            try:
                row.setValue(label, valueType(value))
            except Exception as e:
                print(e)
                print("Problems found converting metadata: ")
                print("Label id = %s" % label)
                print("Attribute = %s" % attr)
                print("Value = %s" % value)
                print("Value type = %s" % valueType)
                raise e
            row.setValue(label, valueType(getattr(obj, attr).get()))

    attrLabels = attrDict.values()

    for label in extraLabels:
        attrName = prefixAttribute(emlib.label2Str(label))
        if label not in attrLabels and hasattr(obj, attrName):
            value = obj.getAttributeValue(attrName)
            row.setValue(label, value)

def locationToXmipp(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Xmipp.
    """
    return ImageHandler.locationToXmipp((index, filename))

def alignmentToRow(alignment, alignmentRow, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
                          -> for xmipp implies alignment
    """
    if alignment is None:
        return

    is2D = alignType == ALIGN_2D
    inverseTransform = alignType == ALIGN_PROJ
    # only flip is meaninfull if 2D case
    # in that case the 2x2 determinant is negative
    flip = False
    matrix = alignment.getMatrix()
    if alignType == ALIGN_2D:
        # get 2x2 matrix and check if negative
        flip = bool(np.linalg.det(matrix[0:2, 0:2]) < 0)
        if flip:
            matrix[0, :2] *= -1.  # invert only the first two columns keep x
            matrix[2, 2] = 1.  # set 3D rot
        else:
            pass

    elif alignType == ALIGN_3D:
        flip = bool(np.linalg.det(matrix[0:3, 0:3]) < 0)
        if flip:
            matrix[0, :4] *= -1.  # now, invert first line including x
            matrix[3, 3] = 1.  # set 3D rot
        else:
            pass

    else:
        flip = bool(np.linalg.det(matrix[0:3, 0:3]) < 0)
        if flip:
            raise Exception("the det of the transformation matrix is "
                            "negative. This is not a valid transformation "
                            "matrix for Scipion.")
    shifts, angles = geometryFromMatrix(matrix, inverseTransform)
    alignmentRow.setValue(emlib.MDL_SHIFT_X, shifts[0])
    alignmentRow.setValue(emlib.MDL_SHIFT_Y, shifts[1])

    if is2D:
        angle = angles[0] + angles[2]
        alignmentRow.setValue(emlib.MDL_ANGLE_PSI,  angle)
    else:
        # if alignType == ALIGN_3D:
        alignmentRow.setValue(emlib.MDL_SHIFT_Z, shifts[2])
        alignmentRow.setValue(emlib.MDL_ANGLE_ROT,  angles[0])
        alignmentRow.setValue(emlib.MDL_ANGLE_TILT, angles[1])
        alignmentRow.setValue(emlib.MDL_ANGLE_PSI,  angles[2])
    alignmentRow.setValue(emlib.MDL_FLIP, flip)

def acquisitionToRow(acquisition, ctfRow):
    """ Set labels values from acquisition to md row. """
    objectToRow(acquisition, ctfRow, ACQUISITION_DICT)

def geometryFromMatrix(matrix, inverseTransform):
    from pwem.convert.transformations import (translation_from_matrix,
                                              euler_from_matrix)
    if inverseTransform:
        matrix = np.linalg.inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    angles = -np.rad2deg(euler_from_matrix(matrix, axes='szyz'))
    return shifts, angles