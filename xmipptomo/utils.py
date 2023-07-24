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
import numpy as np
from pwem.emlib import lib
import pwem.emlib.metadata as md
from pwem.emlib.image import ImageHandler
import pyworkflow as pw
from pyworkflow.object import Set
from tomo.objects import MATRIX_CONVERSION, convertMatrix, TiltSeries, TiltImage
from pwem import ALIGN_PROJ
from xmipp3.convert import alignmentToRow

OUTPUT_TILTSERIES_NAME = "TiltSeries"
OUTPUT_TS_INTERPOLATED_NAME = "InterpolatedTiltSeries"

from pwem.objects import Transform


def calculateRotationAngleAndShiftsFromTM(ti):
    """ This method calculates que tilt image rotation angle from its associated transformation matrix."""

    tm = ti.getTransform().getMatrix()
    cosRotationAngle = tm[0][0]
    sinRotationAngle = tm[1][0]
    rotationAngle = math.degrees(math.atan(sinRotationAngle / cosRotationAngle))
    Sx = tm[0][2]
    Sy = tm[1][2]

    return rotationAngle, Sx, Sy


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

        transform = item.getTransform()
        if transform is None:
            rot = 0
            Sx = 0
            Sy = 0
        else:
            rot, Sx, Sy = calculateRotationAngleAndShiftsFromTM(item)

        tiIndex = item.getLocation()[0]
        fn = str(tiIndex) + "@" + item.getFileName()
        nRow = md.Row()
        nRow.setValue(lib.MDL_IMAGE, fn)

        if item.hasCTF():
            defU = item.getCTF().getDefocusU()
            defV = item.getCTF().getDefocusV()
            defAng = item.getCTF().getDefocusAngle()
            nRow.setValue(lib.MDL_CTF_DEFOCUSU, defU)
            nRow.setValue(lib.MDL_CTF_DEFOCUSU, defV)
            nRow.setValue(lib.MDL_CTF_DEFOCUS_ANGLE, defAng)

        if ts.hasOddEven():
            fnOdd  = item.getOdd()
            fnEven = item.getEven()
            nRow.setValue(lib.MDL_HALF1, fnOdd)
            nRow.setValue(lib.MDL_HALF2, fnEven)
        nRow.setValue(lib.MDL_TSID, tsid)
        tilt = item.getTiltAngle()
        nRow.setValue(lib.MDL_ANGLE_TILT, tilt)
        nRow.setValue(lib.MDL_ANGLE_ROT, rot)
        nRow.setValue(lib.MDL_SHIFT_X, Sx)
        nRow.setValue(lib.MDL_SHIFT_Y, Sy)
        nRow.addToMd(mdts)

        fnts = os.path.join(tomoPath, "%s_ts.xmd" % tsid)

    mdts.write(fnts)

    return fnts

