# **************************************************************************
# *
# * Authors:       Jose Luis Vilas Prieto (jlvilas@cnb.csic.es) [1]
# *                Federico P. de Isidro Gomez (fp.deisidro@cnb.csi.es) [1]
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
# **************************************************************************ç

import os

import numpy as np

from pwem import emlib
from pyworkflow import BETA
from pyworkflow.object import Set
from pyworkflow.protocol.params import PointerParam
import pyworkflow.utils.path as path
from pwem.protocols import EMProtocol
from tomo.protocols import ProtTomoBase
import pwem.emlib.metadata as md
from tomo.objects import SetOfTiltSeriesCoordinates, TiltSeriesCoordinate


class XmippProtCalculate3dCoordinatesFromTS(EMProtocol, ProtTomoBase):
    """
    Scipion protocol to calculate the 3D coordinate based on a set of 2D coordinates obtained from the tilt-series
    using an SPA 2D picker.
    """

    _label = 'Tilt-series calculate coordinates 3D'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        ProtTomoBase.__init__(self)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series associated to the coordinates')

        form.addParam('inputSetOfCoordinates',
                      PointerParam,
                      pointerClass='SetOfCoordinates',
                      important=True,
                      label='Input set of 2D coordinates form tilt-series')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.assignTiltPairs)
        self._insertFunctionStep(self.calculateCoordinates3D)
        self._insertFunctionStep(self._closeOutputSet)

    # --------------------------- STEPS functions ----------------------------

    def convertInputStep(self):
        for tsId in self.getTsIdList():
            extraPrefix = self._getExtraPath(tsId)
            tmpPrefix = self._getTmpPath(tsId)

            path.makePath(tmpPrefix)
            path.makePath(extraPrefix)

        # Write coordinates info .xmd files
        for m in self.inputSetOfCoordinates.get().iterMicrographs():
            mObjId = m.getObjId()

            coordList = []

            for c in self.inputSetOfCoordinates.get():
                if c.getMicId() == mObjId:  # We are doing this because interating sets in scipion is horrible and exasperating
                    coordList.append([c.getObjId(), c.getX(), c.getY()])

            tsId = m._tsId.get()
            avgAngle = m._avgAngle.get()

            tmpPrefix = self._getTmpPath(tsId)
            outputFilePath = os.path.join(tmpPrefix, "%s_%0.f.xmd" % (tsId, avgAngle))

            print(outputFilePath)

            self.writeMicCoordinates(mObjId, coordList, outputFilePath)

    def assignTiltPairs(self):

        for m in self.inputSetOfCoordinates.get().iterMicrographs():
            tsId = m._tsId.get()

            tmpPrefix = self._getTmpPath(tsId)
            extraPrefix = self._getExtraPath(tsId)

            # Tilt pair particles
            fnuntilt = "particles@" + os.path.join(tmpPrefix, "%s_%0.f.xmd" % (tsId, 15.0))
            fntilt = "particles@" + os.path.join(tmpPrefix, "%s_%0.f.xmd" % (tsId, -15.0))
            fnmicsize = m.getFileName()
            maxShift = 50
            threshold = 0.25
            odir = extraPrefix

            params = ' --untiltcoor %s' % fnuntilt
            params += ' --tiltcoor %s' % fntilt
            params += ' --tiltmicsize %s' % fnmicsize
            params += ' --maxshift %f' % maxShift
            params += ' --particlesize %d' % self.getBoxSize()
            params += ' --threshold %f' % threshold
            params += ' --odir %s' % odir
            self.runJob('xmipp_image_assignment_tilt_pair', params)

            # Estimate the tilt axis
            fnposUntilt = "particles@" + os.path.join(extraPrefix, "%s_%0.f.pos" % (tsId, 15.0))
            fnposTilt = "particles@" + os.path.join(extraPrefix, "%s_%0.f.pos" % (tsId, -15.0))
            fnO = os.path.join(extraPrefix, 'tiltAngle.xmd')

            params = ' --untilted %s' % fnposUntilt
            params += ' --tilted %s' % fnposTilt
            params += ' -o %s' % fnO
            self.runJob('xmipp_angular_estimate_tilt_axis', params)

            break  # Esto es para comprobar unicamente una pareja

    def calculateCoordinates3D(self):
        for m in self.inputSetOfCoordinates.get().iterMicrographs():
            tsId = m._tsId.get()

            extraPrefix = self._getExtraPath(tsId)

            fnposUntilt = "particles@" + os.path.join(extraPrefix, "%s_%0.f.pos" % (tsId, 15.0))
            fnposTilt = "particles@" + os.path.join(extraPrefix, "%s_%0.f.pos" % (tsId, -15.0))

            mdCoor1 = md.MetaData()
            mdCoor2 = md.MetaData()

            mdCoor1.read(fnposUntilt)
            mdCoor2.read(fnposTilt)

            # *** calculate from sqlite
            a1 = np.radians(-15)
            a2 = np.radians(15)

            xdim, ydim, _ = self.inputSetOfTiltSeries.get().getDimensions()
            xdim2 = xdim/2
            ydim2 = ydim/2

            self.getOutputSetOfTiltSeriesCoordinates()

            print("check1")

            for objId in mdCoor1:

                print(objId)

                x1 = mdCoor1.getValue(emlib.MDL_XCOOR, objId) - xdim2
                x2 = mdCoor2.getValue(emlib.MDL_XCOOR, objId) - xdim2
                y1 = mdCoor1.getValue(emlib.MDL_YCOOR, objId) - ydim2
                y2 = mdCoor2.getValue(emlib.MDL_YCOOR, objId) - ydim2

                # x1 = x*cos(a1)-z*sin(a1)
                # x2 = x*cos(a2)-z*sin(a2)

                m = np.matrix([[np.cos(a1), -np.sin(a1)], [np.cos(a2), -np.sin(a2)]])
                m_inv = np.linalg.pinv(m)

                x = m_inv[0, 0] * x1 + m_inv[0, 1] * x2
                y = (y1 + y2) / 2
                z = m_inv[1, 0] * x1 + m_inv[1, 1] * x2

                print(x)
                print(y)
                print(z)

                newCoord3D = TiltSeriesCoordinate()
                newCoord3D.setTsId(tsId)
                newCoord3D.setPosition(x,
                                       y,
                                       z,
                                       sampling_rate=self.inputSetOfTiltSeries.get().getSamplingRate())

                self.tiltSeriesCoordinates.append(newCoord3D)

            self.tiltSeriesCoordinates.write()
            self._store()

            break

    def CloseOutputStep(self):
        self.getOutputSetOfTiltSeriesCoordinates().setStreamState(Set.STREAM_CLOSED)

        self._store()

    # --------------------------- UTILS functions ----------------------------
    def getBoxSize(self):
        return self.inputSetOfCoordinates.get().getBoxSize()

    def getTsIdList(self):
        tsIdList = []

        for m in self.inputSetOfCoordinates.get().getMicrographs():
            tsIdList.append(m._tsId.get())

        # Remove tsId duplicated in list
        tsIdList = list(dict.fromkeys(tsIdList))

        return tsIdList

    def writeMicCoordinates(self, micObjId, coordList, outputFn, isManual=True):
        """ Write the pos file as expected by Xmipp with the coordinates
        of a given micrograph.
        Params:
            mic: input micrograph.
            coordList: list of (x, y) coordinates.
            outputFn: output filename for the pos file .
            isManual: if the coordinates are 'Manual' or 'Supervised'
            getPosFunc: a function to get the positions from the coordinate,
                it can be useful for scaling the coordinates if needed.
        """
        state = 'Manual' if isManual else 'Supervised'
        f = self.openMd(outputFn, state)

        for coord in coordList:
            coordObjId = coord[0]
            x = coord[1]
            y = coord[2]
            f.write(" %06d   1   %d  %d  %d   %06d \n"
                    % (coordObjId, x, y, 1, micObjId))

        f.close()

    def openMd(self, fn, state='Manual'):
        # We are going to write metadata directy to file to do it faster
        f = open(fn, 'w')
        ismanual = state == 'Manual'
        block = 'data_particles' if ismanual else 'data_particles_auto'
        s = """# XMIPP_STAR_1 *
#
data_header
loop_
 _pickingMicrographState
%s
%s
loop_
 _itemId
 _enabled
 _xcoor
 _ycoor
 _cost
 _micrographId
    """ % (state, block)
        f.write(s)
        return f

    def getOutputSetOfTiltSeriesCoordinates(self):
        if hasattr(self, "tiltSeriesCoordinates"):
            self.tiltSeriesCoordinates.enableAppend()

        else:
            outputSetOfCoordinates3D = SetOfTiltSeriesCoordinates.create(self._getPath(), suffix='Fiducials3D')

            outputSetOfCoordinates3D.setSetOfTiltSeries(self.inputSetOfTiltSeries)
            outputSetOfCoordinates3D.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(tiltSeriesCoordinates=outputSetOfCoordinates3D)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfCoordinates3D)

        return self.tiltSeriesCoordinates
