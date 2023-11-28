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


def readXmippMetadataEnabledTiltImages(xmdPath):
    """ This method takes a Xmipp metadata (xmd) file containing the enabled images from a tilt series and retrieves a
     matrix containing the enable label and the tilt image location. """

    enableInfoList = []

    mdEnable = md.MetaData(xmdPath)

    for objId in mdEnable:
        imageName = mdEnable.getValue(lib.MDL_IMAGE, objId)
        enable = mdEnable.getValue(lib.MDL_ENABLED, objId)
        imgNumber, imgPath = imageName.split("@")

        enableInfoList.append([enable, imgNumber, imgPath])

    return enableInfoList

    # with open(xmdPath) as f:
    #     enableInfoText = f.read().splitlines()
    #
    # for line in enableInfoText[6:]:
    #     # Split enable and location
    #     vectorLine = line.split()
    #
    #     # Split location in index and path
    #     locationInfo = vectorLine[1].split("@")
    #     enableInfoList.append([vectorLine[0], int(locationInfo[0]), locationInfo[1]])
    #
    # return enableInfoList


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


def writeMdTiltSeries(ts, tomoPath, fnXmd=None):
    """
        Returns a metadata with the tilt series information, TsID, filename and tilt angle.
    """

    mdts = lib.MetaData()
    tsid = ts.getTsId()

    for _, item in enumerate(ts):

        transform = item.getTransform()
        if transform is None:
            rot = 0
            sx = 0
            sy = 0
        else:
            rot, sx, sy = calculateRotationAngleAndShiftsFromTM(item)

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
            fnOdd = item.getOdd()
            fnEven = item.getEven()
            nRow.setValue(lib.MDL_HALF1, fnOdd)
            nRow.setValue(lib.MDL_HALF2, fnEven)
        nRow.setValue(lib.MDL_TSID, tsid)
        tilt = item.getTiltAngle()
        nRow.setValue(lib.MDL_ANGLE_TILT, tilt)
        nRow.setValue(lib.MDL_ANGLE_ROT, rot)
        nRow.setValue(lib.MDL_SHIFT_X, sx)
        nRow.setValue(lib.MDL_SHIFT_Y, sy)
        nRow.addToMd(mdts)

        fnts = os.path.join(tomoPath, "%s_ts.xmd" % tsid)

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


def parseLandmarkCoordinatesFile(lmFile):
    """ This function retrive a list of landmark coordinates form xmd file as generated by xmipp program
    xmipp_tomo_detect_landmarks"""

    lmInfo = []
    table = emtable.Table(fileName=lmFile)

    for row in table.iterRows(fileName='noname@' + lmFile):
        xCoor = row.get('xcoor')
        yCoor = row.get('ycoor')
        tiltIm = row.get('zcoor')

        lmInfo.append([xCoor, yCoor, tiltIm])

    return lmInfo


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
