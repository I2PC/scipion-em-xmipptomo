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

from tomo.constants import BOTTOM_LEFT_CORNER


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

    with open(filePath, 'w') as f:
        f.write(xmdHeader)
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldNames)

        for ci in coordinatesInfo:
            writer.writerow({'x': ci[0],
                             'y': ci[1],
                             'z': ci[2]})
