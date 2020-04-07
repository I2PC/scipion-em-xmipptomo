# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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
import numpy as np
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import IntParam
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms


class XmippProtPhantomSubtomo(EMProtocol, ProtTomoBase):
    """ Create subtomogram phantoms """

    _label = 'phantom create subtomo'

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('nsubtomos', IntParam, label='Number of subtomograms', help="How many phantom subtomograms")
        form.addParam('dim', IntParam, label='Dimension of subtomograms', help="dimx = dimy = dimz")
        form.addParam('rotmin', IntParam, label='Min rot angle')
        form.addParam('rotmax', IntParam, label='Max rot angle')
        form.addParam('tiltmin', IntParam, label='Min tilt angle')
        form.addParam('tiltmax', IntParam, label='Max tilt angle')
        form.addParam('psimin', IntParam, label='Min psi angle')
        form.addParam('psimax', IntParam, label='Max psi angle')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createPhantomsStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def createPhantomsStep(self):
        fnDescr = self._getExtraPath("phantom.descr")
        fhDescr = open(fnDescr, 'w')
        dim = self.dim.get()
        fhDescr.write("%s %s %s 0\ncyl + 1 0 0 0 15 15 2 0 0 0\nsph + 1 0 0 5 2\ncyl + 1 0 0 -5 2 2 10 0 90 0"
                      % (dim, dim, dim))
        fhDescr.close()
        fnVol = self._getExtraPath("phantom.vol")
        self.runJob("xmipp_phantom_create", " -i %s -o %s" % (fnDescr, fnVol))
        for i in range(self.nsubtomos.get()):
            fnPhantomi = self._getExtraPath("phantom%03d.vol" % i)
            rot = np.random.randint(self.rotmin.get(), self.rotmax.get())
            tilt = np.random.randint(self.tiltmin.get(), self.tiltmax.get())
            psi = np.random.randint(self.psimin.get(), self.psimax.get())
            self.runJob("xmipp_transform_geometry",
                        " -i %s -o %s --rotate_volume euler %d %d %d "
                        % (fnVol, fnPhantomi, rot, tilt, psi))

    def createOutputStep(self):
        self.outputSubTomogramsSet = self._createSetOfSubTomograms(self._getOutputSuffix(SetOfSubTomograms))
        for item in self.getInputTomograms().iterItems():
            for ind, tomoFile in enumerate(self.tomoFiles):
                if os.path.basename(tomoFile) == os.path.basename(item.getFileName()):
                    coordSet = self.lines[ind]
                    outputSet = self.readSetOfSubTomograms(
                        self._getExtraPath(pwutils.replaceBaseExt(tomoFile, "hdf")),
                        self.outputSubTomogramsSet, coordSet)

        self._defineOutputs(outputSetOfSubtomogram=outputSet)
        self._defineSourceRelation(self.inputCoordinates, outputSet)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputSubtomograms'):
            summary.append("Output phantom not ready yet.")
        else:
            summary.append("%s phantoms created with random orientations" % self.nsubtomos.get())
        return summary

    def _methods(self):
        if not hasattr(self, 'outputSubtomograms'):
            return ["Output phantoms not ready yet."]
        else:
            return ["%s phantoms created with random orientations" % self.nsubtomos.get()]