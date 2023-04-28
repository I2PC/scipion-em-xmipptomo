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
import os

import emtable
from tomo.constants import BOTTOM_LEFT_CORNER
from pwem.emlib import lib
import numpy as np
import pwem.emlib.metadata as md
from tomo.objects import MATRIX_CONVERSION, convertMatrix, TiltSeries, TiltImage
from pwem.objects import Transform


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

def xmdToTiltSeries(fnXmd, sampling=1, tsid='defaulttsId'):
    """
    This function takes a metadata files as input and stores the Tilt Series in the SQLite database
    """
    mdts = md.MetaData(fnXmd)
    counter = 0
    print('entro')
    for objId in md:
        fnImg = mdts.getValue(lib.MDL_IMAGE, objId)
        tilt = mdts.getValue(lib.MDL_ANGLE_TILT, objId)
        #tsid = mdts.getValue(lib.MDL_TSID, objId)

        if counter == 0:
            newTs = TiltSeries(tsId=tsid)

        newTs.copyInfo(tsid)
        newTi = TiltImage()
        newTi.setLocation(fnImg)
        newTi.setTiltAngle(tilt)
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

    for item in ts:
        transform = item.getTransform()

        #if transform is None:
        #    tm = convertMatrix(np.eye(4))
        #else:
        #    tm = transform.getMatrix(convention=MATRIX_CONVERSION.XMIPP)
        #Maq = Transform(matrix=tm)

        tiIndex = item.getLocation()[0]
        fn = str(tiIndex) + "@" + item.getFileName()
        nRow = md.Row()
        nRow.setValue(lib.MDL_IMAGE, fn)
        nRow.setValue(lib.MDL_TSID, tsid)
        nRow.setValue(lib.MDL_ANGLE_TILT, item.getTiltAngle())
        # nRow.setValue(lib.MDL_ANGLE_ROT, int(coord.getY(const.BOTTOM_LEFT_CORNER)))
        # alignmentToRow(Maq, nRow, ALIGN_PROJ)
        nRow.addToMd(mdts)

    if ts is None:
        fnts = os.path.join(tomoPath, "%s_ts.xmd" % tsid)
    else:
        fnts = os.path.join(tomoPath, fnXmd)

    print(fnts)
    mdts.write(fnts)

    return fnts

