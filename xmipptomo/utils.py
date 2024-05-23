# **************************************************************************
# *
# * Authors:     Federico P. de Isidro Gomez (fp.deisidro@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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
"""
This module contains utils functions for xmipp tomo protocols
"""

# General imports
import math
import csv
import os
import shutil

# Scipion em imports
import emtable
from pwem import ALIGN_PROJ
from pwem.objects import Integer, CTFModel, Transform
from pwem.emlib import lib
import pwem.emlib.metadata as md
from pwem.emlib.image import ImageHandler
import pyworkflow as pw

# External plugin imports
from tomo.objects import TiltSeries, TiltImage, SetOfCTFTomoSeries
from tomo.objects import MATRIX_CONVERSION, TiltSeries, TiltImage
from tomo.constants import BOTTOM_LEFT_CORNER
from xmipp3.convert import alignmentToRow

OUTPUT_TILTSERIES_NAME = "TiltSeries"
OUTPUT_TS_INTERPOLATED_NAME = "InterpolatedTiltSeries"

def calculateRotationAngleAndShiftsFromTM(ti):
    """ This method calculates the rot and shifts of a tilt image from its associated transformation matrix."""
    transform = ti.getTransform()
    if transform is None:
        rotationAngle = 0.0
        sx = 0.0
        sy = 0.0
    else:
        tm = transform.getMatrix()
        cosRotationAngle = tm[0][0]
        sinRotationAngle = tm[1][0]
        rotationAngle = float(math.degrees(math.atan(sinRotationAngle / cosRotationAngle)))
        sx = tm[0][2]
        sy = tm[1][2]

    return rotationAngle, sx, sy

def readXmdStatisticsFile(fnmd):
    xPos = []
    yPos = []
    zPos = []
    avg = []
    std = []

    table = emtable.Table(fileName=fnmd)

    for row in table.iterRows(fileName='noname@' + fnmd):
        avg.append(row.get('avg'))
        std.append(row.get('stddev'))
        xPos.append(row.get('xcoor'))
        yPos.append(row.get('ycoor'))
        zPos.append(row.get('zcoor'))

    return xPos, yPos, zPos, avg, std


def tiltSeriesParticleToXmd(tsParticle):
    mdtsp = lib.MetaData()
    for ti in tsParticle:
        tm = ti.getTransformationMatrix()
        fn = ti.parseFileName()
        nRow = md.Row()
        nRow.setValue(lib.MDL_IMAGE, fn)
        alignmentToRow(tm, nRow, ALIGN_PROJ)
        nRow.addToMd(mdtsp)


def writeOutputCoordinates3dXmdFile(soc, filePath, tomoId=None):
    """ Generates a 3D coordinates xmd file from the set of coordinates associated to a given tomogram (identified by
     its tomo tomoId). If no tomoId is input the xmd output file will contain all the coordinates belonging to the
     set. """

    xmdHeader = "# XMIPP_STAR_1 *\n" \
                "#\n" \
                "data_noname\n" \
                "loop_\n" \
                " _xcoor\n" \
                " _ycoor\n" \
                " _zcoor\n"

    coordinatesInfo = []
    fieldNames = ['x', 'y', 'z']

    for coord in soc.iterCoordinates(tomoId):
        coordinatesInfo.append([coord.getX(BOTTOM_LEFT_CORNER),
                                coord.getY(BOTTOM_LEFT_CORNER),
                                coord.getZ(BOTTOM_LEFT_CORNER)])

    with open(filePath, 'w') as f:
        f.write(xmdHeader)
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldNames)

        for ci in coordinatesInfo:
            writer.writerow({'x': ci[0],
                             'y': ci[1],
                             'z': ci[2]})


def xmdToTiltSeries(outputSetOfTs, inTs, fnXmd, sampling=1, odir='', tsid='defaulttsId', suffix=''):
    """
    This function takes a metadata files as input and stores the Tilt Series
    """
    mdts = md.MetaData(fnXmd)
    counter = 1

    ih = ImageHandler()
    newTs = TiltSeries(tsId=tsid)
    newTs.copyInfo(inTs, copyId=True)
    outputSetOfTs.append(newTs)
    fnStack = os.path.join(odir, tsid + suffix + '.mrcs')

    for objId in mdts:
        fnImg = os.path.join(odir, mdts.getValue(lib.MDL_IMAGE, objId))
        mdts.getValue(lib.MDL_ANGLE_TILT, objId)

        originalTi = inTs[counter]
        newTi = TiltImage()
        newTi.copyInfo(originalTi, copyId=True, copyTM=True)
        newTi.setOddEven([])

        newLocation = (counter, fnStack)
        ih.convert(fnImg, newLocation)
        pw.utils.cleanPath(fnImg)
        newTi.setLocation(newLocation)

        newTs.append(newTi)
        newTs.setSamplingRate(sampling)
        counter = counter + 1
    return newTs


def writeMdTiltSeries(ts, tomoPath):
    """
        Returns a metadata with the tilt series information, TsID, filename and tilt angle.
    """

    mdts = lib.MetaData()
    tsid = ts.getTsId()

    fnts = os.path.join(tomoPath, "%s_ts.xmd" % tsid)
    for _, item in enumerate(ts):

        rot, sx, sy = calculateRotationAngleAndShiftsFromTM(item)

        tiIndex = item.getLocation()[0]
        fn = str(tiIndex) + "@" + item.getFileName()
        nRow = md.Row()
        nRow.setValue(lib.MDL_IMAGE, fn)

        '''
        if item.hasCTF():
            ctf = item.getCTF()
            defU = ctf.getCTF().getDefocusU()
            defV = ctf.getCTF().getDefocusV()
            defAng = ctf.getCTF().getDefocusAngle()
            nRow.setValue(lib.MDL_CTF_DEFOCUSU, defU)
            nRow.setValue(lib.MDL_CTF_DEFOCUSU, defV)
            nRow.setValue(lib.MDL_CTF_DEFOCUS_ANGLE, defAng)
        '''
        '''
        if ts.hasOddEven():
            fnOdd = item.getOdd()
            fnEven = item.getEven()
            nRow.setValue(lib.MDL_HALF1, fnOdd)
            nRow.setValue(lib.MDL_HALF2, fnEven)
        '''

        nRow.setValue(lib.MDL_TSID, tsid)
        tilt = item.getTiltAngle()
        nRow.setValue(lib.MDL_ANGLE_TILT, tilt)
        nRow.setValue(lib.MDL_ANGLE_ROT, rot)
        nRow.setValue(lib.MDL_SHIFT_X, sx)
        nRow.setValue(lib.MDL_SHIFT_Y, sy)
        nRow.addToMd(mdts)
    mdts.write(fnts)

    return fnts

def getCTFfromId(setOfCTFs: SetOfCTFTomoSeries, targetTsId: Integer) -> CTFModel:
    """
    This function returns the CTF from the set with the given target Tilt series id. 
    """
    # Iterate CTF set looking for the one with targetTsId
    for ctf in setOfCTFs:
        # If ctf id matches target TS id, return such CTF
        if targetTsId == ctf.getTsId():
            return ctf

def setGeometricalParametersToRow(row, fnZero, rot, tilt, psi, sx, sy, defU, defV, defAng):
    if defAng is None:
        defAng = 0.0
    row.setValue(lib.MDL_IMAGE, fnZero)
    row.setValue(lib.MDL_CTF_DEFOCUSU, defU)
    row.setValue(lib.MDL_CTF_DEFOCUSV, defV)
    row.setValue(lib.MDL_CTF_DEFOCUS_ANGLE, defAng)
    row.setValue(lib.MDL_ANGLE_TILT, tilt)
    row.setValue(lib.MDL_ANGLE_ROT, rot)
    row.setValue(lib.MDL_ANGLE_PSI, psi)
    row.setValue(lib.MDL_SHIFT_X, sx)
    row.setValue(lib.MDL_SHIFT_Y, sy)
    return row

def removeTmpElements(tmpElements):
    """ This function removes all given temporary files and directories. """
    # Removing selected elements
    for item in tmpElements:
        if os.path.exists(item):
            if os.path.isdir(item):
                shutil.rmtree(item)
            else:
                os.remove(item)


def retrieveXmipp3dCoordinatesIntoList(coordFilePath, xmdFormat=0):
    """ This method takes a xmipp metadata (xmd) 3D coordinates file path and returns a list of tuples containing
    every coordinate. This method also transform the coordinates into the Scipion convention. This method allows
    different xmd formats containing coordinates information:
        format=0: plain coordinates, xmd files only contains (x, y, z) values.
        format=1: coordinates with alignment information, xmd files contains also shifts and angle values."""

    coorList = []

    with open(coordFilePath) as f:
        inputLines = f.readlines()

    if xmdFormat == 0:
        for line in inputLines[7:]:
            vector = line.split()

            coorList.append([float(vector[0]),
                             float(vector[1]),
                             float(vector[2])])

    if xmdFormat == 1:
        for line in inputLines[15:]:
            vector = line.split()

            coorList.append([float(vector[-3]),
                             float(vector[-2]),
                             float(vector[-1])])

    return coorList


def writeMdCoordinates(setOfCoordinates, tomo, fnCoor):
    """
        Write the xmd file containing the set of coordinates corresponding to the given tomogram at the specified
        location
    """
    mdCoor = lib.MetaData()

    tsid = tomo.getTsId()

    coordDict = []
    lines = []

    fnCoorDirectory = os.path.dirname(fnCoor)
    if not os.path.exists(fnCoorDirectory):
        os.makedirs(fnCoorDirectory)

    for item in setOfCoordinates.iterCoordinates(volume=tomo):
        coord = item
        transform = Transform(matrix=item.getMatrix(convention=MATRIX_CONVERSION.XMIPP))

        if coord.getTomoId() == tsid:
            nRow = md.Row()
            nRow.setValue(lib.MDL_ITEM_ID, int(coord.getObjId()))
            coord.setVolume(tomo)

            nRow.setValue(lib.MDL_XCOOR, int(coord.getX(BOTTOM_LEFT_CORNER)))
            nRow.setValue(lib.MDL_YCOOR, int(coord.getY(BOTTOM_LEFT_CORNER)))
            nRow.setValue(lib.MDL_ZCOOR, int(coord.getZ(BOTTOM_LEFT_CORNER)))

            alignmentToRow(transform, nRow, ALIGN_PROJ)
            nRow.addToMd(mdCoor)

            newCoord = item.clone()
            newCoord.setVolume(coord.getVolume())
            coordDict.append(newCoord)
            lines.append(coordDict)

    mdCoor.write(fnCoor)

    return fnCoor
