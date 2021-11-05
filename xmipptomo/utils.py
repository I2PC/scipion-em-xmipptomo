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


def calculateRotationAngleFromTM(ti):
    """ This method calculates que tilt image rotation angle from its associated transformation matrix."""

    tm = ti.getTransform().getMatrix()
    cosRotationAngle = tm[0][0]
    sinRotationAngle = tm[1][0]
    rotationAngle = math.degrees(math.atan(sinRotationAngle/cosRotationAngle))

    return rotationAngle


def retrieveXmipp3dCoordinatesIntoList(ts, angleFilePath):
    """ This method takes a Scipion tilt-series object and return a xmipp metadata (xmd) tilt angle file contaiing
    every angle of each tilt-image. """

    header = "# XMIPP_STAR_1 * " \
             "#" \
             "data_noname" \
             "loop_" \
             "_angleTilt"

    angleList = []

    for ti in ts:
        angleList.append(ti.getTiltAngle())

    with open(angleFilePath, 'w') as f:
        f.write(header)
        f.writelines("%s\n" % angle for angle in angleList)
