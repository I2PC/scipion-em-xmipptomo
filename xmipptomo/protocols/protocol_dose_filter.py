# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csi.es) [1]
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

import os
import imod.utils as utils
import pwem.objects as data

from pyworkflow import utils as pwutils
from pyworkflow import BETA
from pyworkflow.object import Set
import pyworkflow.protocol.params as params
import pyworkflow.utils.path as path

import pwem.emlib.metadata as md
import pwem.emlib as emlib
from pwem.protocols import EMProtocol

from tomo.protocols.protocol_base import ProtTomoImportFiles
from tomo.protocols import ProtTomoBase
import tomo.objects as tomoObj

from xmipptomo import utils

SCIPION_IMPORT = 0
FIXED_DOSE = 1


class XmippProtDoseFilter(ProtTomoImportFiles, EMProtocol, ProtTomoBase):
    """
    Tilt-series' dose filtering based on  T. Grant, N. Grigorieff, eLife 2015
    More info:
        https://doi.org/10.7554/eLife.06980
    """

    _label = 'Dose filter'
    _devStatus = BETA

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection('Input')

        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of tilt-series to be filtered.')

        form.addParam('initialDose',
                      params.FloatParam,
                      default=0.0,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Initial dose (e/sq A)',
                      help='Dose applied before any of the images in the input file were taken; this value will be '
                           'added to all the prior dose values, however they were obtained.')

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        for ts in self.inputSetOfTiltSeries.get():
            self._insertFunctionStep(self.doseFilterStep, ts.getObjId())
            self._insertFunctionStep(self.createOutputStep, ts.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions ----------------------------
    def doseFilterStep(self, tsObjId):
        """Apply the dose fitler to every tilt series"""

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        path.makePath(tmpPrefix)
        path.makePath(extraPrefix)

        firstItem = ts.getFirstItem()

        fnMd = 'image_and_dose.xmd'
        mdDose = md.MetaData()
        idx = 1
        for ti in ts:
            doseValue = ti.getAcquisition().getAccumDose()
            fn = ti.getFileName()
            strimg = str(idx) + '@' + fn
            idx = idx + 1

            mdRow = md.Row()
            mdRow.setValue(emlib.MDL_IMAGE, strimg)
            mdRow.setValue(emlib.MDL_DOSE, doseValue)
            mdRow.writeToMd(mdDose, mdDose.addObject())

        mdDose.write(fnMd)

        params = ' -i %s '          % fnMd
        params += ' -o %s '         % (os.path.join(extraPrefix, os.path.splitext(os.path.basename(ts.getFileName()))[0]  + '.mrcs'))
        params += ' --sampling %s' % self.inputSetOfTiltSeries.get().getSamplingRate()
        params += ' --voltage %f '  % ts.getAcquisition().getVoltage()

        self.runJob('xmipp_tomo_tiltseries_dose_filter', params)


    def createOutputStep(self, tsObjId):
        """Generate output filtered tilt series"""

        ts = self.inputSetOfTiltSeries.get()[tsObjId]
        tsId = ts.getTsId()

        extraPrefix = self._getExtraPath(tsId)

        self.getOutputSetOfTiltSeries()

        newTs = tomoObj.TiltSeries(tsId=tsId)
        newTs.copyInfo(ts)

        self.outputSetOfTiltSeries.append(newTs)

        for index, tiltImage in enumerate(ts):
            newTi = tomoObj.TiltImage()
            newTi.copyInfo(tiltImage, copyId=True, copyTM=True)
            newTi.setAcquisition(tiltImage.getAcquisition())


            #To be reviewed in hunting day
            pathTi = os.path.join(extraPrefix, os.path.splitext(os.path.basename(ts.getFileName()))[0] + '.mrcs')
            newTi.setLocation(index + 1, pathTi) #(os.path.join(extraPrefix, tiltImage.parseFileName())))

            newTs.append(newTi)

        newTs.write(properties=False)

        self.outputSetOfTiltSeries.update(newTs)
        self.outputSetOfTiltSeries.write()

        self._store()

    def getOutputSetOfTiltSeries(self):
        if hasattr(self, "outputSetOfTiltSeries"):
            self.outputSetOfTiltSeries.enableAppend()
        else:
            outputSetOfTiltSeries = self._createSetOfTiltSeries()
            outputSetOfTiltSeries.copyInfo(self.inputSetOfTiltSeries.get())
            outputSetOfTiltSeries.setDim(self.inputSetOfTiltSeries.get().getDim())
            outputSetOfTiltSeries.setSamplingRate(self.inputSetOfTiltSeries.get().getSamplingRate())
            outputSetOfTiltSeries.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(outputSetOfTiltSeries=outputSetOfTiltSeries)
            self._defineSourceRelation(self.inputSetOfTiltSeries, outputSetOfTiltSeries)
        return self.outputSetOfTiltSeries

    def closeOutputSetsStep(self):
        self.outputSetOfTiltSeries.setStreamState(Set.STREAM_CLOSED)
        self.outputSetOfTiltSeries.write()
        self._store()

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        validateMsgs = []

        if self.inputDoseType.get() == SCIPION_IMPORT:
            for ts in self.inputSetOfTiltSeries.get():
                if ts.getFirstItem().getAcquisition().getDosePerFrame() == None:
                    validateMsgs.append("%s has no dose information stored in Scipion Metadata. To solve this import "
                                        "the tilt-series with the mdoc option." % ts.getTsId())

        return validateMsgs

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           % (self.inputSetOfTiltSeries.get().getSize(),
                              self.outputSetOfTiltSeries.getSize()))
        elif hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            summary.append("Input Tilt-Series: %d.\nTransformation matrices calculated: %d.\n"
                           "Interpolated Tilt-Series: %d.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputSetOfTiltSeries.getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The transformation matrix has been calculated for %d "
                           "Tilt-series using the IMOD procedure.\n"
                           % (self.outputSetOfTiltSeries.getSize()))
        elif hasattr(self, 'outputInterpolatedSetOfTiltSeries'):
            methods.append("The transformation matrix has been calculated for %d "
                           "Tilt-series using the IMOD procedure.\n"
                           "Also, interpolation has been completed for %d Tilt-series.\n"
                           % (self.outputSetOfTiltSeries.getSize(),
                              self.outputInterpolatedSetOfTiltSeries.getSize()))
        else:
            methods.append("Output classes not ready yet.")
        return methods
