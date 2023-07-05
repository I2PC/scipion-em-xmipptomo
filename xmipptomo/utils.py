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

import shutil
import math
import csv
import os.path

import os
import emtable

from pwem import ALIGN_PROJ
from pwem.emlib import lib
import pwem.emlib.metadata as md
from pwem.objects import Transform

import pyworkflow as pw
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import MATRIX_CONVERSION
from xmipp3.convert import alignmentToRow
from pwem.emlib import lib
import pwem.emlib.metadata as md
from pwem.emlib.image import ImageHandler
from tomo.objects import TiltSeries, TiltImage

OUTPUT_TILTSERIES_NAME = "TiltSeries"
OUTPUT_TS_INTERPOLATED_NAME = "InterpolatedTiltSeries"

def calculateRotationAngleFromTM(ti):
    """ This method calculates que tilt image rotation angle from its associated transformation matrix."""

    tm = ti.getTransform().getMatrix()
    cosRotationAngle = tm[0][0]
    sinRotationAngle = tm[1][0]
    rotationAngle = math.degrees(math.atan(sinRotationAngle / cosRotationAngle))

    return rotationAngle


def readXmdStatisticsFile(fnmd):
    x_pos = []
    y_pos = []
    z_pos = []
    avg = []
    std = []

    table = emtable.Table(fileName=fnmd)

    for row in table.iterRows(fileName='noname@' + fnmd):
        avg.append(row.get('avg'))
        std.append(row.get('stddev'))
        x_pos.append(row.get('xcoor'))
        y_pos.append(row.get('ycoor'))
        z_pos.append(row.get('zcoor'))

    return x_pos, y_pos, z_pos, avg, std


def writeOutputCoordinates3dXmdFile(soc, filePath, tomoId=None):
    """ Generates a 3D coordinates xmd file from the set of coordinates associated to a given tomogram (identified by
     its tomo tomoId). If no tomoId is input the the xmd output file will contain all the coordinates belonging to the
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
    fnStack = os.path.join(odir, tsid + suffix +'.mrcs')

    for objId in mdts:
        fnImg = os.path.join(odir, mdts.getValue(lib.MDL_IMAGE, objId))
        tilt = mdts.getValue(lib.MDL_ANGLE_TILT, objId)

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

def writeMdTiltSeries(ts, tomoPath, fnXmd=None):
    """
        Returns a metadata with the tilt series information, TsID, filename and tilt angle.
    """
    mdts = lib.MetaData()

    tsid = ts.getTsId()

    for index, item in enumerate(ts):

        #transform = item.getTransform()
        #if transform is None:
        #    tm = convertMatrix(np.eye(4))
        #else:
        #    tm = transform.getMatrix(convention=MATRIX_CONVERSION.XMIPP)
        #Maq = Transform(matrix=tm)

        tiIndex = item.getLocation()[0]
        fn = str(tiIndex) + "@" + item.getFileName()
        nRow = md.Row()
        nRow.setValue(lib.MDL_IMAGE, fn)
        if ts.hasOddEven():
            fnOdd  = item.getOdd()
            fnEven = item.getEven()
            nRow.setValue(lib.MDL_HALF1, fnOdd)
            nRow.setValue(lib.MDL_HALF2, fnEven)
        nRow.setValue(lib.MDL_TSID, tsid)
        nRow.setValue(lib.MDL_ANGLE_TILT, item.getTiltAngle())
        # nRow.setValue(lib.MDL_ANGLE_ROT, int(coord.getY(const.BOTTOM_LEFT_CORNER)))
        # alignmentToRow(Maq, nRow, ALIGN_PROJ)
        nRow.addToMd(mdts)

        fnts = os.path.join(tomoPath, "%s_ts.xmd" % tsid)

    mdts.write(fnts)

    return fnts

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

    fnCoor_directory = os.path.dirname(fnCoor)
    if not os.path.exists(fnCoor_directory):
        os.makedirs(fnCoor_directory)

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
