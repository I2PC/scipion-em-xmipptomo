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
import pwem.objects.data as data
from tomo.objects import TiltImage, TiltImageBase, SetOfTiltSeries, SetOfTiltSeriesBase, TomoAcquisition, Transform, SubTomogram


class TiltParticle(TiltImage, data.CTFModel):
    """
    Represents a single 2D image of the set of images that defines a subtomogram
    """
    def __init__(self, **kwargs):
        TiltImage.__init__(self, **kwargs)
        self._boxSize = Integer()
        self._coordX = Float()
        self._coordY = Float()
        self._ctfCorrected = None
        #TODO: use transform of tilt image

    def setCoordX(self, xvalue):
        self._coordX = xvalue

    def getCoordX(self):
        return self._coordX

    def setCoordY(self, yvalue):
        self._coordY = yvalue

    def getCoordY(self):
        return self._coordY

    def setTransformationMatrix(self, tm):
        self._transfMatrix = tm

    def getTransformationMatrix(self):
        return self._transfMatrix

    def setctfCorrected(self, isCorrected):
        self._ctfCorrected = isCorrected

    def getctfCorrected(self):
        return self._ctfCorrected


class TiltSeriesParticle(data.SetOfImages, SubTomogram):
    ITEM_TYPE = TiltParticle

    def __init__(self, **kwargs):
        data.SetOfImages.__init__(self, **kwargs)
        self._ctfCorrected = Boolean(False)


class SetOfTiltSeriesParticle(data.SetOfImages):
    ITEM_TYPE = TiltSeriesParticle

    def __init__(self, **kwargs):
        data.SetOfImages.__init__(self, **kwargs)
        self._tsId = String(kwargs.get('tsId', None))
        # TiltSeries will always be used inside a SetOfTiltSeries
        # so, let's do not store the mapper path by default
        self._mapperPath.setStore(False)
        self._acquisition = TomoAcquisition()
        self._anglesCount = Integer()
        self._hasAlignment = Boolean(False)
        self._ctfCorrected = Boolean(False)

    def setAnglesCount(self, value):
        self._anglesCount.set(value)

    def getAnglesCount(self):
        self._anglesCount.get()

    def setHasAlignment(self, booleanvalue):
        self._hasAlignment.set(booleanvalue)

    def hasAlignment(self):
        self._hasAlignment.get()

    def setCtfCorrected(self, booleanvalue):
        self._ctfCorrected.set(booleanvalue)

    def hasCtfCorrected(self):
        self._ctfCorrected.get()

    def copyInfo(self, other):
        """ Copy basic information (sampling rate and ctf)
        from other set of images to current one"""
        self.copyAttributes(other, '_samplingRate', '_hasCtf',
                            '_alignment', '_isPhaseFlipped', '_isAmplitudeCorrected')

    def __str__(self):
        """ String representation of a set of coordinates. """
        return "%s (%d items, %s x %s x %s, %0.2f Å/px)" % ('SetOfTiltSeriesParticle', self.getSize(), self._anglesCount,
                                                         self.getDim()[0], self.getDim()[1], self.getSamplingRate())





