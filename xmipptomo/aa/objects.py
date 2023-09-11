# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
# *
# * your institution
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
from pyworkflow.object import Integer, Float, String, Pointer, Boolean, CsvList
from tomo.objects import TiltImage, SetOfTiltSeries, TiltSeries, TomoAcquisition, Transform, SubTomogram


class TiltParticle(TiltImage):
    """
    Represents a single 2D image of the set of images that defines a subtomogram
    """

    def __init__(self, **kwargs):
        # Tilt particles have the same information than the tilt images
        # This means the same transformation matrix (except the shifts that are corrected)
        TiltImage.__init__(self, **kwargs)
        # Tilt particles are assume square with edge equal to the boxsize
        self._boxSize = Integer()
        # The position of the tilt image is given by 2D coordinates
        self._coord2DX = Float()
        self._coord2DY = Float()
        self._tsId = String()
        # The Euler angles of the tilt particle are described by the transformation matrix (_transfom)
        # of the tilt image. However, in the refinement the refined angles and shift will be slightly
        # different from the measured in the tilt image. The transformation defined by these deviations is
        # store in _deltaTransform. The transformation matrix of the tilt particle will be the combination
        # of D (deltaTransform), M (tilt image transformation matrix), and R (subtomogram transformation matrix)
        # The final Euler matrix will be RMD
        self._deltaTransform = None

        # Note the CTF is defined in the TiltImage class, but is local and belongs to each tilt particle

    def getTsId(self):
        return self._tsId

    def setTsId(self, value):
        self._tsId.set(value)

    def setCoord2DX(self, xvalue):
        self._coord2DX.set(xvalue)

    def getCoord2DX(self):
        return self._coord2DX

    def setCoord2DY(self, yvalue):
        self._coord2DY.set(yvalue)

    def getCoord2DY(self):
        return self._coord2DY

    def getBoxSize(self):
        return self._boxSize

    def setBoxSize(self, boxsize):
        self._boxSize.set(boxsize)

    def getDeltaTransform(self):
        return self._deltaTransform

    def setDeltaTransform(self, matrix):
        self._deltaTransform = matrix

    def __str__(self):

        s = super().__str__()

        return s


class TiltSeriesParticle(TiltSeries):
    ITEM_TYPE = TiltParticle

    def __init__(self,  **kwargs):
        TiltSeries.__init__(self, **kwargs)
        #self._mapperPath.setStore(False)
        self._originalTs = String('ts')
        # The boxsize defines the length of the edge of the particle
        self._boxsize = Integer()

        # This is the tilt series Id to identify the tilt series where the particle was picked
        self._tsId = String()

        # Each tilt Series particle represents a subtomogram. The transformation matrix to align the subtomogram
        # is given by the next matrix
        self._transformSubtomo = None

        # The next set of coordinates define the position of the particle in the tomogram
        self._coord3DX = Integer()
        self._coord3DY = Integer()
        self._coord3DZ = Integer()



    def getOriginalTs(self):
        return self._originalTs

    def setOriginalTs(self, value):
        self._originalTs.set(value)

    '''
    def getFileName(self):
        return self._fileName

    def setFileName(self, fn):
        self._fileName = fn
    '''

    def getTsId(self):
        return self._tsId

    def setTsId(self, value):
        self._tsId.set(value)

    def getBoxSize(self):
        return self._boxsize

    def setBoxSize(self, value):
        self._boxsize.set(value)

    def getTransformSubtomo(self):
        return self._transformSubtomo

    def setTransformSubtomo(self, matrix):
        self._transformSubtomo.set(matrix)

    def setCoord3DX(self, xvalue):
        self._coord3DX.set(xvalue)

    def getCoord3DX(self):
        return self._coord3DX

    def setCoord3DY(self, yvalue):
        self._coord3DY.set(yvalue)

    def getCoord3DY(self):
        return self._coord3DY

    def setCoord3DZ(self, zvalue):
        self._coord3DY.set(zvalue)

    def getCoord3DZ(self):
        return self._coord3DZ




class SetOfTiltSeriesParticle(SetOfTiltSeries):
    ITEM_TYPE = TiltSeriesParticle

    def __init__(self, **kwargs):
        SetOfTiltSeries.__init__(self, **kwargs)

    def __str__(self):
        """ String representation of a set of coordinates. """
        return "%s (%d items, %s x %s x %s, %0.2f Å/px)" % ('SetOfTiltSeriesParticle', self.getSize(), self._anglesCount,
                                                         self.getDim()[0], self.getDim()[1], self.getSamplingRate())

'''
    def __str__(self):

        s = super().__str__()

        return s
'''