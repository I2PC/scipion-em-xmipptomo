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
import os
import emtable

import pyworkflow as pw
from tomo.constants import BOTTOM_LEFT_CORNER
from pwem.emlib import lib
import pwem.emlib.metadata as md
from pwem.emlib.image import ImageHandler
from tomo.objects import TiltSeries, TiltImage


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

OUTPUT_TILTSERIES_NAME = "TiltSeries"
OUTPUT_TS_INTERPOLATED_NAME = "InterpolatedTiltSeries"

def calculateRotationAngleFromTM(ti):
    """ This method calculates que tilt image rotation angle from its associated transformation matrix."""

    tm = ti.getTransform().getMatrix()
    cosRotationAngle = tm[0][0]
    sinRotationAngle = tm[1][0]
    rotationAngle = math.degrees(math.atan(sinRotationAngle / cosRotationAngle))

    return rotationAngle


def writeXmippMetadataTiltAngleList(ts, angleFilePath):
    """ This method takes a Scipion tilt-series object and return a Xmipp metadata (xmd) tilt angle file containing
    every angle of each tilt-image. """

    header = "# XMIPP_STAR_1 * \n" \
             "#\n" \
             "data_noname\n" \
             "loop_\n" \
             "_angleTilt\n"

    angleList = []

    for ti in ts:
        angleList.append(ti.getTiltAngle())

    with open(angleFilePath, 'w') as f:
        f.write(header)
        f.writelines("%s\n" % angle for angle in angleList)


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


def readXmippMetadataEnabledTiltImages(xmdPath):
    """ This method takes an Xmipp metadata (xmd) file containing the enabled images from a tilt series and retrieves a
     matrix containing the enable info and the tilt image location. """

    enableInfoList = []

    with open(xmdPath) as f:
        enableInfoText = f.read().splitlines()

    for line in enableInfoText[6:]:
        # Split enable and location
        vectorLine = line.split()

        # Split location in index and path
        locationInfo = vectorLine[1].split("@")
        enableInfoList.append([vectorLine[0], int(locationInfo[0]), locationInfo[1]])

    return enableInfoList


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

    if len(coordinatesInfo) == 0:
        return False

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
    return True


def writeOutputTiltSeriesCoordinates3dXmdFile(soc, filePath, sr, halfX, halfY, tsId=None):
    """ Generates a 3D coordinates xmd file from the set of coordinates associated to a given tilt-series (identified by
     its tomo tsId). If no tsId is input the xmd output file will contain all the coordinates belonging to the
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

    if tsId is None:
        for coord in soc:
            coordinatesInfo.append([(coord.getX()/sr)-halfX,
                                    (coord.getY()/sr)-halfY,
                                    (coord.getZ()/sr)])
    else:
        for coord in soc:
            if coord.getTsId() == tsId:
                coordinatesInfo.append([(coord.getX()/sr)+halfX,
                                        (coord.getY()/sr)+halfY,
                                        (coord.getZ()/sr)])

    if len(coordinatesInfo) == 0:
        return False

    with open(filePath, 'w') as f:
        f.write(xmdHeader)
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldNames)

        for ci in coordinatesInfo:
            writer.writerow({'x': ci[0],
                             'y': ci[1],
                             'z': ci[2]})

    return True


def readResidualStatisticsXmdFile(xmdFilePath):
    """ This method takes the file path of a Xmipp metadata file (.xmd) and generates a dictionary with all the
    information associated to the residuals from each landmark model: convex hull area and perimeter, statistical
    tests passed and failed, and its associated coordinate. """

    def paramType(string):
        """
            0: convex hull area
            1: convex hull parameter
            2: statistical test
        """

        if string == "chArea":
            return 0
        elif string == "chPerim":
            return 1
        else:
            return 2

    statisticsInfoTable = {}

    table = emtable.Table(fileName=xmdFilePath)

    for row in table.iterRows(fileName='noname@'+xmdFilePath):
        en = row.get('enabled')
        name = str(row.get('image'))
        min = row.get('min')  # convex hull area/parameter or p-value
        _ = row.get('max')  # convex hull area/parameter or p-value pondered by FDR
        xCoor = row.get('xcoor')
        yCoor = row.get('ycoor')
        zCoor = row.get('zcoor')

        key, test = name.split('_')
        parType = paramType(test)

        if key in statisticsInfoTable.keys():
            # Convex hull area
            if parType == 0:
                statisticsInfoTable[key][0] = min

            # Convex hull perimeter
            elif parType == 1:
                statisticsInfoTable[key][1] = min

            # Passed tests
            elif parType == 2 and en == 1:
                statisticsInfoTable[key][2].append(test)

            # Failed tests
            elif parType == 2 and en == -1:
                statisticsInfoTable[key][3].append(test)

        else:
            # Convex hull area
            if parType == 0:
                statisticsInfoTable[key] = [min, 0, [], [], [xCoor, yCoor, zCoor]]

            # Convex hull perimeter
            elif parType == 1:
                statisticsInfoTable[key] = [0, min, [], [], [xCoor, yCoor, zCoor]]

            # Passed tests
            elif parType == 2 and en == 1:
                statisticsInfoTable[key] = [0, 0, [test], [], [xCoor, yCoor, zCoor]]

            # Failed tests
            elif parType == 2 and en == -1:
                statisticsInfoTable[key] = [0, 0, [], [test], [xCoor, yCoor, zCoor]]

    return statisticsInfoTable


def calculateAverageRotationAngleFromTM(ts):
    """ This method calculates que average tilt image rotation angle from its associated transformation matrix."""
    avgRotationAngle = 0

    if not ts.getFirstItem().hasTransform():
        return avgRotationAngle

    for ti in ts:
        tm = ti.getTransform().getMatrix()
        cosRotationAngle = tm[0][0]
        sinRotationAngle = tm[1][0]
        avgRotationAngle += math.degrees(math.atan(sinRotationAngle/cosRotationAngle))

    avgRotationAngle = avgRotationAngle / ts.getSize()

    return avgRotationAngle

