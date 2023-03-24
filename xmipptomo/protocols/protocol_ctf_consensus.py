# **************************************************************************
# *
# * Authors: Daniel Marchan Torres    (da.marchan@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import logging
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow import BETA
import pyworkflow.object as pwobj

from tomo.objects import SetOfCTFTomoSeries
from pwem import emlib
import xmipp3
from xmipp3.convert import setXmippAttribute, getScipionObj

from cmath import rect, phase
from math import radians, degrees

logger = logging.getLogger(__name__)

OUTPUT_BAD_CTF_SERIE = "badCTFTomoSeries"
OUTPUT_GOOD_CTF_SERIE = "goodCTFTomoSeries"

DEFOCUS = 'defocus'
ASTIGMATISM = 'astigmatism'
RESOLUTION = 'resolution'


class XmippProtCTFTomoSeriesConsensus(EMProtocol):
    """
     Validate a set of CTF tomo series and separate into two sets (good and
     bad tomo series )
    """
    _label = 'ctf consensus tomo'
    _devStatus = BETA

    _countDict = {}
    _dfCons = {}
    _astigCons = {}
    _freqResol = {}

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        """ Define input parameters from this program into the given form. """
        form.addSection(label='Input')
        form.addParam('inputCtfTomoSeries', params.PointerParam, important=True,
                      pointerClass='SetOfCTFTomoSeries',
                      label='Input ctf tomo series')

        form.addParam('inputCtfTomoSeries2', params.PointerParam, important=True,
                      pointerClass='SetOfCTFTomoSeries',
                      label='Input secondary ctf tomo series',
                      help='CTF tomo series to be compared with reference CTF')
        form.addSection(label='Consensus')
        form.addParam('validationType', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=["Global", "Per tilt"],
                      label="Validation type",
                      default=0,
                      help="Global mode: the series with at least "
                           "one image that does not satisfy the criteria are "
                           "rejected \n"
                           "Per tilt mode: the series containing a certain "
                           "number of images that does not satisfy the "
                           "criteria are rejected")
        form.addParam('numberImages', params.IntParam,
                      label='Number of images to rejected',
                      condition='validationType==1', default=None,
                      allowsNull=True,
                      help="Number of images taking into account to rejected a "
                           "ctf series")
        self.addCriteriaParams(form)
        form.addParam('averageDefocus', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Average equivalent metadata?',
                      help='If *Yes*, making an average of those metadata present '
                           'in both CTF estimations (defocus, astigmatism angle...)\n '
                           'If *No*, the primary estimation metadata will persist.')
        form.addParam('includeSecondary', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Include all secondary metadata?',
                      help='If *Yes*, all metadata in the *Secondary CTF* will '
                           'be included in the resulting CTF.\n '
                           'If *No*, only the primary metadata (plus consensus '
                           'scores) will be in the resulting CTF.')

    def addCriteriaParams(self, form):
        form.addParam('defocusCriteria', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=["Yes", "No"],
                      label="Defocus tolerance",
                      default=1, help="Validate the defocus deviation taking "
                                      "into account a threshold(tolerance) "
                                      "respect to a defocus expected value.")
        form.addParam('defocusTolerance', params.FloatParam,
                      label='Tolerance value defocus',
                      condition='defocusCriteria==0', default=0.1,
                      help="Defocus tolerance value")

        form.addParam('astigmatismCriteria', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=["Yes", "No"],
                      label="Astigmatism",
                      default=1, help="Validate the astigmatism taking into "
                                      "account a tolerance value.")
        form.addParam('astigmatismTolerance', params.FloatParam,
                      label='Tolerance value astigmatism',
                      condition='astigmatismCriteria==0', default=0.1,
                      help="Astigmatism tolerance value")

        form.addParam('resolutionCriteria', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=["Yes", "No"],
                      label="Resolution",
                      default=1, help="Validate the resolution taking into "
                                      "account a expected resolution.")
        form.addParam('resolutionTolerance', params.FloatParam,
                      label='Minimum consensus resolution (A)',
                      allowsNull=True,
                      condition='resolutionCriteria==0', default=None,
                      help="Minimum value for the consensus resolution in "
                           "Angstroms.\nIf there are noticeable discrepancies "
                           "between the two estimations below this resolution, "
                           "it will be discarded. 'Option for calculating consensus resolution. "
                           "The algorithm assumes that two CTF are "
                           "consistent if the phase (wave aberration function) "
                           "of the two CTFs are closer than 90 degrees.\n"
                           "The reported consensusResolution is the resolution "
                           "at which the two CTF phases differ in 90 degrees.'")

# ---------------------INSERT ALL STEPS --------------------------------------
    def _insertAllSteps(self):
        self.initializeParams()
        self._insertFunctionStep(self.ctfConsensusStep)
        self._insertFunctionStep(self.createOutputStep)

# ----------------------------STEPS --------------------------------------
    def initializeParams(self):
        self.setSecondaryAttributes()

    def setSecondaryAttributes(self):
        if self.includeSecondary:
            item = self.inputCtfTomoSeries.get().getFirstItem()
            ctf1Attr = set(item.getObjDict().keys())

            item = self.inputCtfTomoSeries.get().getFirstItem()
            ctf2Attr = set(item.getObjDict().keys())
            self.secondaryAttributes = ctf2Attr - ctf1Attr
        else:
            self.secondaryAttributes = set()

    def ctfConsensusStep(self):
        """
        Validate all ctf tomo series and separate into two sets(good and bad
        following the selected criteria)
        """
        self.goodCTFTomoSeries = None
        self.badCTFTomoSeries = None
        self.validateDict = {}
        imagesToReject = 1  # Case of global validation type

        ctfSeriesSet1 = self.inputCtfTomoSeries.get()
        ctfSeriesSetDict1 = {series.getTsId(): series.getObjId() for series in ctfSeriesSet1.iterItems()}
        ctfSeriesSetIds1 = ctfSeriesSetDict1.keys()
        ctfSeriesSet2 = self.inputCtfTomoSeries2.get()
        ctfSeriesSetDict2 = {series.getTsId(): series.getObjId() for series in ctfSeriesSet2.iterItems()}
        ctfSeriesSetIds2 = ctfSeriesSetDict2.keys()
        sharedSeriesIDs = list(set(ctfSeriesSetIds1).intersection(set(ctfSeriesSetIds2)))
        print('Shared TS')
        print(sharedSeriesIDs)

        if self.validationType.get() == 1:
            imagesToReject = self.numberImages.get()

        for tsId in sharedSeriesIDs:
            newCTFTomoSeries = ctfSeriesSet1[ctfSeriesSetDict1[tsId]].clone()
            ctfEstItems = []
            failedDefocusCriteria = 0
            failedAstigmatismCriteria = 0
            failedResolutionCriteria = 0

            ctfSeriesDict1 = {ctf.getObjId(): ctf.clone() for ctf in ctfSeriesSet1[ctfSeriesSetDict1[tsId]].iterItems()}
            ctfSeriesIds1 = ctfSeriesDict1.keys()
            ctfSeriesDict2 = {ctf.getObjId(): ctf.clone() for ctf in ctfSeriesSet2[ctfSeriesSetDict2[tsId]].iterItems()}
            ctfSeriesIds2 = ctfSeriesDict2.keys()
            sharedCtfIDs = list(set(ctfSeriesIds1).intersection(set(ctfSeriesIds2)))
            print('Shared tilt images')
            print(sharedCtfIDs)

            dfDict = {}
            astigDict = {}
            resolDict = {}
            # TODO: un dict q guarde el angulo al que se saca todo esto
            md1 = emlib.MetaData()
            md2 = emlib.MetaData()
            for ctfId in sharedCtfIDs:
                ctf1 = ctfSeriesDict1.get(ctfId)
                ctf2 = ctfSeriesDict2.get(ctfId)
                # Defocus criteria
                if self.defocusCriteria.get() == 0:
                    df = self._calculateConsensusDefocus(ctf1, ctf2)
                    dfDict[ctfId] = df
                    if df < self.defocusTolerance.get():
                        failedDefocusCriteria += 1
                # Astigmatism criteria
                if self.astigmatismCriteria.get() == 0:
                    astig = self._calculateConsensusAstigmatism(ctf1, ctf2)
                    astigDict[ctfId] = astig
                    if astig < self.astigmatismTolerance.get():
                        failedAstigmatismCriteria += 1
                # Resolution criteria
                if self.resolutionCriteria.get() == 0:
                    res = self._calculateConsensusResolution(md1, md2, ctf1, ctf2)
                    setAttribute(ctf1, '_consensus_resolution', res)
                    resolDict[ctfId] = res
                    if res < self.resolutionTolerance.get():
                        failedResolutionCriteria += 1

                if self.averageDefocus.get():
                    ctf = self._fillCTFDefocus(ctf1, ctf2)
                else:
                    ctf = ctf1

                # TODO: - Secondary data
                if self.includeSecondary.get():
                    for attr in self.secondaryAttributes:
                        copyAttribute(ctf2, ctf, attr)

                ctfEstItems.append(ctf)

            restDict = {}
            if self.defocusCriteria.get() == 0:
                self._dfCons[tsId] = dfDict
                restDict.update({DEFOCUS: failedDefocusCriteria})
            if self.astigmatismCriteria.get() == 0:
                self._astigCons[tsId] = astigDict
                restDict.update({ASTIGMATISM: failedAstigmatismCriteria})
            if self.resolutionCriteria.get() == 0:
                self._freqResol[tsId] = resolDict
                restDict.update({RESOLUTION: failedResolutionCriteria})

            self._countDict[tsId] = restDict

            if self.validateConsensusCtf(tsId, imagesToReject):
                newCTFTomoSeries.setEnabled(True)
                output = self.getOutputSetOfCTFTomoSeries(OUTPUT_GOOD_CTF_SERIE)
                output.append(newCTFTomoSeries)
            else:
                newCTFTomoSeries.setEnabled(False)
                output = self.getOutputSetOfCTFTomoSeries(OUTPUT_BAD_CTF_SERIE)
                output.append(newCTFTomoSeries)

            for ctfItem in ctfEstItems:
                newCTFTomoSeries.append(ctfItem)

        self.info(self._countDict)


    def _ctfToMd(self, ctf, ctfMd):
        """ Write the proper metadata for Xmipp from a given CTF """
        ctfMd.clear()
        ctfRow = emlib.metadata.Row()
        xmipp3.convert.ctfModelToRow(ctf, ctfRow)
        ctfRow.addToMd(ctfMd)

    def _calculateConsensusDefocus(self, ctfTomo1, ctfTomo2):
        avg1 = (ctfTomo1.getDefocusU() + ctfTomo1.getDefocusV())/2
        avg2 = (ctfTomo2.getDefocusU() + ctfTomo2.getDefocusV())/2
        return _calculateDiff(avg1, avg2)


    def _calculateConsensusAstigmatism(self, ctfTomo1, ctfTomo2):
        astig1 = abs(ctfTomo1.getDefocusU() - ctfTomo1.getDefocusV())
        astig2 = abs(ctfTomo2.getDefocusU() - ctfTomo2.getDefocusV())
        return _calculateDiff(astig1, astig2)

    def _calculateConsensusResolution(self, md1, md2, ctf1, ctf2):
        # Podras definir los metadata aqui pero puedo perjudicar en tiempo
        try:
            self._ctfToMd(ctf1, md1)
            self._ctfToMd(ctf2, md2)
            res = emlib.errorMaxFreqCTFs2D(md1, md2)  # TODO relacionar el angulo a la cons_res
        except TypeError as exc:
            print("Error reading ctf for id:%s. %s" % (ctf1.getObjId(), exc))
            res = 0  # more coherent number
        return res

    def validateConsensusCtf(self, tsId, imagesToReject=1):
        """
         Validate the ctftomo serie.
         Return true if  the image surpass the threshold limits
        """
        criteria = True
        if self.defocusCriteria.get() == 0:
            criteria = True if self._countDict[tsId][DEFOCUS] < imagesToReject else False
        if self.astigmatismCriteria.get() == 0:
            criteria = True if self._countDict[tsId][ASTIGMATISM] < imagesToReject else False
        if self.resolutionCriteria.get() == 0:
            criteria = True if self._countDict[tsId][RESOLUTION] < imagesToReject else False

        return criteria

    def _getSetOfTiltSeries(self, pointer=False):
        return self.inputCtfTomoSeries.get().getSetOfTiltSeries(pointer=pointer)

    def _fillCTFDefocus(self, ctf1, ctf2):
        newDefocusU = 0.5 * (ctf1.getDefocusU() + ctf2.getDefocusU())
        newDefocusV = 0.5 * (ctf1.getDefocusV() + ctf2.getDefocusV())
        newDefocusAngle = averageAngles(ctf1.getDefocusAngle(), ctf2.getDefocusAngle())
        ctf1.setStandardDefocus(newDefocusU, newDefocusV, newDefocusAngle)
        return ctf1

    def getOutputSetOfCTFTomoSeries(self, outputSetName):
        outputSetOfCTFTomoSeries = getattr(self, outputSetName, None)

        if outputSetOfCTFTomoSeries:
            outputSetOfCTFTomoSeries.enableAppend()
        else:
            outputSetOfCTFTomoSeries = SetOfCTFTomoSeries.create(self._getPath(),
                                                                 prefix=outputSetName)
            outputSetOfCTFTomoSeries.setSetOfTiltSeries(self._getSetOfTiltSeries())
            outputSetOfCTFTomoSeries.setStreamState(pwobj.Set.STREAM_OPEN)
            self._defineOutputs(**{outputSetName: outputSetOfCTFTomoSeries})

        return outputSetOfCTFTomoSeries

    def createOutputStep(self):
        if self.goodCTFTomoSeries is not None:
            self.goodCTFTomoSeries.setStreamState(pwobj.Set.STREAM_CLOSED)
            self.goodCTFTomoSeries.write()
            self._store()
        if self.badCTFTomoSeries is not None:
            self.badCTFTomoSeries.setStreamState(pwobj.Set.STREAM_CLOSED)
            self.badCTFTomoSeries.write()
            self._store()

    def _summary(self):
        return self._countDict

# -------------------------------- UTILS -----------------------------------------
def _calculateDiff(value1, value2):
    return abs(value1 - value2)/(0.5 * (value1 + value2))

def setAttribute(obj, label, value):
    if value is None:
        return
    setattr(obj, label, getScipionObj(value))

def averageAngles(angle1, angle2):
    c1 = rect(1, radians(angle1*2))
    c2 = rect(1, radians(angle2*2))
    return degrees(phase((c1 + c2)*0.5))/2


def copyAttribute(src, dst, label, default=None):
    setAttribute(dst, label, getattr(src, label, default))
