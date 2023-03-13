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

import math
import csv
import os.path

import emtable

from pwem import ALIGN_PROJ
from pwem.emlib import lib
import pwem.emlib.metadata as md
from pwem.objects import Transform
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import MATRIX_CONVERSION
from xmipp3.convert import alignmentToRow


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


def retrieveXmipp3dCoordinatesIntoList(coordFilePath):
    """ This method takes an xmipp metadata (xmd) 3D coordinates file path and returns a list of tuples containing
    every coordinate. This method also transform the coordinates into the Scipion convention. """

    coorList = []

    with open(coordFilePath) as f:
        inputLines = f.readlines()

    for line in inputLines[7:]:
        vector = line.split()

        coorList.append([float(vector[0]),
                         float(vector[1]),
                         float(vector[2])])

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
