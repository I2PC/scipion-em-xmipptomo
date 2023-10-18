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
from pwem.constants import (NO_INDEX, ALIGN_NONE, ALIGN_PROJ, ALIGN_2D,
                            ALIGN_3D)
from pwem.objects import Acquisition, Transform
from pyworkflow.object import ObjectWrap, String
import numpy as np
from collections import OrderedDict
from xmipp3.base import getLabelPythonType, iterMdRows
from tomo.objects import SubTomogram, Coordinate3D, CTFTomo

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

    CTF_DICT = OrderedDict([
           ("_defocusU", emlib.MDL_CTF_DEFOCUSU),
           ("_defocusV", emlib.MDL_CTF_DEFOCUSV),
           ("_defocusAngle", emlib.MDL_CTF_DEFOCUS_ANGLE),
           ("_resolution", emlib.MDL_CTF_CRIT_MAXFREQ),
           ("_fitQuality", emlib.MDL_CTF_CRIT_FITTINGSCORE)
           ])

    # TODO: remove next dictionary when all
    # cTFmodel has resolution and fitQuality
    CTF_DICT_NORESOLUTION = OrderedDict([
            ("_defocusU", emlib.MDL_CTF_DEFOCUSU),
            ("_defocusV", emlib.MDL_CTF_DEFOCUSV),
            ("_defocusAngle", emlib.MDL_CTF_DEFOCUS_ANGLE)
            ])

    CTF_PSD_DICT = OrderedDict([
           ("_psdFile", emlib.MDL_PSD),
           (prefixAttribute("enhanced_psd"), emlib.MDL_PSD_ENHANCED),
           (prefixAttribute("ctfmodel_quadrant"), emlib.MDL_IMAGE1),
           (prefixAttribute("ctfmodel_halfplane"), emlib.MDL_IMAGE1)
           ])

    CTF_EXTRA_LABELS = [
        emlib.MDL_CTF_CA,
        emlib.MDL_CTF_ENERGY_LOSS,
        emlib.MDL_CTF_LENS_STABILITY,
        emlib.MDL_CTF_CONVERGENCE_CONE,
        emlib.MDL_CTF_LONGITUDINAL_DISPLACEMENT,
        emlib.MDL_CTF_TRANSVERSAL_DISPLACEMENT,
        emlib.MDL_CTF_K,
        emlib.MDL_CTF_BG_GAUSSIAN_K,
        emlib.MDL_CTF_BG_GAUSSIAN_SIGMAU,
        emlib.MDL_CTF_BG_GAUSSIAN_SIGMAV,
        emlib.MDL_CTF_BG_GAUSSIAN_CU,
        emlib.MDL_CTF_BG_GAUSSIAN_CV,
        emlib.MDL_CTF_BG_SQRT_K,
        emlib.MDL_CTF_BG_SQRT_U,
        emlib.MDL_CTF_BG_SQRT_V,
        emlib.MDL_CTF_BG_SQRT_ANGLE,
        emlib.MDL_CTF_BG_BASELINE,
        emlib.MDL_CTF_BG_GAUSSIAN2_K,
        emlib.MDL_CTF_BG_GAUSSIAN2_SIGMAU,
        emlib.MDL_CTF_BG_GAUSSIAN2_SIGMAV,
        emlib.MDL_CTF_BG_GAUSSIAN2_CU,
        emlib.MDL_CTF_BG_GAUSSIAN2_CV,
        emlib.MDL_CTF_BG_GAUSSIAN2_ANGLE,
        emlib.MDL_CTF_CRIT_FITTINGCORR13,
        emlib.MDL_CTF_CRIT_ICENESS,
        emlib.MDL_CTF_VPP_RADIUS,
        emlib.MDL_CTF_DOWNSAMPLE_PERFORMED,
        emlib.MDL_CTF_CRIT_PSDVARIANCE,
        emlib.MDL_CTF_CRIT_PSDPCA1VARIANCE,
        emlib.MDL_CTF_CRIT_PSDPCARUNSTEST,
        emlib.MDL_CTF_CRIT_FIRSTZEROAVG,
        emlib.MDL_CTF_CRIT_DAMPING,
        emlib.MDL_CTF_CRIT_FIRSTZERORATIO,
        emlib.MDL_CTF_CRIT_PSDCORRELATION90,
        emlib.MDL_CTF_CRIT_PSDRADIALINTEGRAL,
        emlib.MDL_CTF_CRIT_NORMALITY,
        # In xmipp the ctf also contains acquisition information
        emlib.MDL_CTF_Q0,
        emlib.MDL_CTF_CS,
        emlib.MDL_CTF_VOLTAGE,
        emlib.MDL_CTF_SAMPLING_RATE
        ]

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
    alignmentRow.setValue(emlib.MDL_SHIFT_Z, shifts[2])

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

def readSetOfSubtomograms(filename, partSet, **kwargs):
    readSetOfImages(filename, partSet, rowToParticle, **kwargs)

def readSetOfImages(filename, imgSet, rowToFunc, **kwargs):
    """read from Xmipp image metadata.
        filename: The metadata filename where the image are.
        imgSet: the SetOfParticles that will be populated.
        rowToFunc: this function will be used to convert the row to Object
    """
    imgMd = emlib.MetaData(filename)

    # By default remove disabled items from metadata
    # be careful if you need to preserve the original number of items
    if kwargs.get('removeDisabled', True):
        imgMd.removeDisabled()

    # If the type of alignment is not sent through the kwargs
    # try to deduced from the metadata labels
    if 'alignType' not in kwargs:
        imgRow = rowFromMd(imgMd, imgMd.firstObject())
        if _containsAny(imgRow, ALIGNMENT_DICT):
            if imgRow.containsLabel(emlib.MDL_ANGLE_TILT):
                kwargs['alignType'] = ALIGN_PROJ
            else:
                kwargs['alignType'] = ALIGN_2D
        else:
            kwargs['alignType'] = ALIGN_NONE

    if imgMd.size() > 0:
        for objId in imgMd:
            imgRow = rowFromMd(imgMd, objId)
            img = rowToFunc(imgRow, **kwargs)
            imgSet.append(img)

        imgSet.setHasCTF(img.hasCTF())
        imgSet.setAlignment(kwargs['alignType'])

def rowToParticle(partRow, **kwargs):
    return _rowToParticle(partRow, SubTomogram, **kwargs)

def _rowToParticle(partRow, particleClass, **kwargs):
    """ Create a Particle from a row of a metadata. """
    # Since postprocessImage is intended to be after the object is
    # setup, we need to intercept it here and call it at the end
    postprocessImageRow = kwargs.get('postprocessImageRow', None)
    if postprocessImageRow:
        del kwargs['postprocessImageRow']

    img = rowToImage(partRow, emlib.MDL_IMAGE, particleClass, **kwargs)
    img.setCoordinate3D(rowToCoordinate(partRow))
    # copy micId if available
    # if not copy micrograph name if available
    try:
        if partRow.hasLabel(emlib.MDL_MICROGRAPH_ID):
            img.setMicId(partRow.getValue(emlib.MDL_MICROGRAPH_ID))
#        elif partRow.hasLabel(emlib.MDL_MICROGRAPH):
#            micName = partRow.getValue(emlib.MDL_MICROGRAPH)
#            img._micrograph = micName
#            print("setting micname as %s" % micName)
#            img.printAll()
#            print("getAttributes1 %s" % img._micrograph)
#            print("getAttributes2 %s" % getattr(img, "_micrograph", 'kk')
#        else:
#            print("WARNING: No micname")
    except Exception as e:
        print("Warning:", e.message)

    if postprocessImageRow:
        postprocessImageRow(img, partRow)

    return img

def rowFromMd(mdIn, objId):
    row = md.Row()
    row.readFromMd(mdIn, objId)
    return row

def _containsAny(row, labels):
    """ Check if the labels (values) in labelsDict
    are present in the row.
    """
    values = labels.values() if isinstance(labels, dict) else labels
    return any(row.containsLabel(l) for l in values)

def rowToImage(imgRow, imgLabel, imgClass, **kwargs):
    """ Create an Image from a row of a metadata. """
    img = imgClass()

    # Provide a hook to be used if something is needed to be
    # done for special cases before converting image to row
    preprocessImageRow = kwargs.get('preprocessImageRow', None)
    if preprocessImageRow:
        preprocessImageRow(img, imgRow)

    # Decompose Xmipp filename
    index, filename = xmippToLocation(imgRow.getValue(imgLabel))
    img.setLocation(index, filename)

    if imgRow.containsLabel(emlib.MDL_REF):
        img.setClassId(imgRow.getValue(emlib.MDL_REF))
    elif imgRow.containsLabel(emlib.MDL_REF3D):
        img.setClassId(imgRow.getValue(emlib.MDL_REF3D))

    if kwargs.get('readCtf', True):
        img.setCTF(rowToCtfModel(imgRow))

    # alignment is mandatory at this point, it shoud be check
    # and detected defaults if not passed at readSetOf.. level
    alignType = kwargs.get('alignType')

    if alignType != ALIGN_NONE:
        img.setTransform(rowToAlignment(imgRow, alignType))

    if kwargs.get('readAcquisition', True):
        img.setAcquisition(rowToAcquisition(imgRow))

    if kwargs.get('magnification', None):
        img.getAcquisition().setMagnification(kwargs.get("magnification"))

    setObjId(img, imgRow)
    # Read some extra labels
    rowToObject(imgRow, img, {}, [])

    # Provide a hook to be used if something is needed to be
    # done for special cases before converting image to row
    postprocessImageRow = kwargs.get('postprocessImageRow', None)
    if postprocessImageRow:
        postprocessImageRow(img, imgRow)

    return img

def rowToCoordinate(coordRow):
    """ Create a Coordinate from a row of a metadata. """
    # Check that all required labels are present in the row
    if _containsAll(coordRow, COOR_DICT):
        coord = Coordinate3D()
        rowToObject(coordRow, coord, COOR_DICT)

        # Setup the micId if is integer value
        try:
            coord.setTomoId(int(coordRow.getValue(emlib.MDL_MICROGRAPH_ID)))
        except Exception:
            pass
    else:
        coord = None

    return coord

def xmippToLocation(xmippFilename):
    """ Return a location (index, filename) given
    a Xmipp filename with the index@filename structure. """
    if '@' in xmippFilename:
        return emlib.FileName(xmippFilename).decompose()
    else:
        return NO_INDEX, str(xmippFilename)

def rowToCtfModel(ctfRow):
    """ Create a CTFModel from a row of a metadata. """
    # Check if the row has CTF values, this could be called from a xmipp
    # particles metadata
    if _containsAll(ctfRow, CTF_DICT_NORESOLUTION):

        # for compatibility reason ignore resolution and fitQuality
        # Instantiate Scipion CTF Model
        ctfModel = CTFTomo()

        # Case for metadata coming with Xmipp resolution label
        # Populate Scipion CTF from metadata row (using mapping dictionary
        # plus extra labels
        if ctfRow.hasLabel(md.MDL_CTF_PHASE_SHIFT):
            ctfModel.setPhaseShift(ctfRow.getValue(md.MDL_CTF_PHASE_SHIFT, 0))
        if ctfRow.containsLabel(emlib.label2Str(emlib.MDL_CTF_CRIT_MAXFREQ)):
            rowToObject(ctfRow, ctfModel, CTF_DICT,
                        extraLabels=CTF_EXTRA_LABELS)
        else:
            rowToObject(ctfRow, ctfModel, CTF_DICT_NORESOLUTION)

        # Standarize defocus values
        ctfModel.standardize()
        # Set psd file names
        setPsdFiles(ctfModel, ctfRow)
        # ctfModel.setPhaseShift(0.0)  # for consistency with ctfModel

    else:
        ctfModel = None

    return ctfModel

def rowToAlignment(alignmentRow, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
        """
    is2D = alignType == ALIGN_2D
    inverseTransform = alignType == ALIGN_PROJ

    if _containsAny(alignmentRow, ALIGNMENT_DICT):
        alignment = Transform()
        angles = np.zeros(3)
        shifts = np.zeros(3)
        flip = alignmentRow.getValue(emlib.MDL_FLIP)

        shifts[0] = alignmentRow.getValue(emlib.MDL_SHIFT_X, 0.)
        shifts[1] = alignmentRow.getValue(emlib.MDL_SHIFT_Y, 0.)
        if not is2D:
            angles[0] = alignmentRow.getValue(emlib.MDL_ANGLE_ROT, 0.)
            angles[1] = alignmentRow.getValue(emlib.MDL_ANGLE_TILT, 0.)
            shifts[2] = alignmentRow.getValue(emlib.MDL_SHIFT_Z, 0.)
            angles[2] = alignmentRow.getValue(emlib.MDL_ANGLE_PSI, 0.)
            if flip:
                angles[1] = angles[1]+180  # tilt + 180
                angles[2] = - angles[2]    # - psi, COSS: this is mirroring X
                shifts[0] = -shifts[0]     # -x
        else:
            psi = alignmentRow.getValue(emlib.MDL_ANGLE_PSI, 0.)
            rot = alignmentRow.getValue(emlib.MDL_ANGLE_ROT, 0.)
            if rot != 0. and psi != 0:
                print("HORROR rot and psi are different from zero in 2D case")
            angles[0] = \
                alignmentRow.getValue(emlib.MDL_ANGLE_PSI, 0.)\
                + alignmentRow.getValue(emlib.MDL_ANGLE_ROT, 0.)

        matrix = matrixFromGeometry(shifts, angles, inverseTransform)

        if flip:
            if alignType == ALIGN_2D:
                matrix[0, :2] *= -1.  # invert only the first two columns
                # keep x
                matrix[2, 2] = -1.  # set 3D rot
            elif alignType == ALIGN_3D:
                matrix[0, :3] *= -1.  # now, invert first line excluding x
                matrix[3, 3] *= -1.
            elif alignType == ALIGN_PROJ:
                pass

        alignment.setMatrix(matrix)

        # FIXME: now are also storing the alignment parameters since
        # the conversions to the Transform matrix have not been extensively
        # tested.
        # After this, we should only keep the matrix
        # for paramName, label in ALIGNMENT_DICT.iter():
        #    if alignmentRow.hasLabel(label):
        #        setattr(alignment, paramName,
        #                alignmentRow.getValueAsObject(label))
    else:
        alignment = None

    return alignment

def rowToAcquisition(acquisitionRow):
    """ Create an acquisition from a row of a metadata. """
    if _containsAll(acquisitionRow, ACQUISITION_DICT):
        acquisition = Acquisition()
        rowToObject(acquisitionRow, acquisition, ACQUISITION_DICT)
    else:
        acquisition = None

    return acquisition

def rowToObject(row, obj, attrDict, extraLabels=[]):
    """ This function will convert from a Row to an EMObject.
    Params:
        row: the Row instance (input)  -see emlib.metadata.utils.Row()-
        obj: the EMObject instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and
            row MDLabels in Xmipp (values).
        extraLabels: a list with extra labels that could be included
            as _xmipp_labelName
    """
    obj.setEnabled(row.getValue(emlib.MDL_ENABLED, 1) > 0)

    for attr, label in attrDict.items():
        value = row.getValue(label)
        if not hasattr(obj, attr):
            setattr(obj, attr, ObjectWrap(value))
        else:
            getattr(obj, attr).set(value)

    attrLabels = attrDict.values()

    for label in extraLabels:
        if label not in attrLabels and row.hasLabel(label):
            labelStr = emlib.label2Str(label)
            setattr(obj, prefixAttribute(labelStr), row.getValueAsObject(label))

def setObjId(obj, mdRow, label=emlib.MDL_ITEM_ID):
    if mdRow.containsLabel(label):
        obj.setObjId(mdRow.getValue(label))
    else:
        obj.setObjId(None)

def _containsAll(row, labels):
    """ Check if the labels (values) in labelsDict
    are present in the row.
    """
    values = labels.values() if isinstance(labels, dict) else labels
    return all(row.containsLabel(l) for l in values)

def setPsdFiles(ctfModel, ctfRow):
    """ Set the PSD files of CTF estimation related
    to this ctfModel. The values will be read from
    the ctfRow if present.
    """
    for attr, label in CTF_PSD_DICT.items():
        if ctfRow.containsLabel(label):
            setattr(ctfModel, attr, String(ctfRow.getValue(label)))

def matrixFromGeometry(shifts, angles, inverseTransform):
    """ Create the transformation matrix from a given
    2D shifts in X and Y...and the 3 euler angles.
    """
    from pwem.convert import euler_matrix
    radAngles = -np.deg2rad(angles)

    M = euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    if inverseTransform:
        M[:3, 3] = -shifts[:3]
        M = np.linalg.inv(M)
    else:
        M[:3, 3] = shifts[:3]

    return M