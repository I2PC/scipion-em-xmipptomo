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
from pwem.convert.transformations import euler_matrix
from pwem.objects.data import Transform
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import IntParam, FloatParam, EnumParam, PointerParam, StringParam, BooleanParam
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms, SubTomogram, TomoAcquisition


class XmippProtPhantomSubtomo(EMProtocol, ProtTomoBase):
    """ Create subtomogram phantoms """

    _label = 'phantom create subtomo'

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('option', EnumParam, choices=['Import volume', 'Create'], default=0,
                      display=EnumParam.DISPLAY_HLIST, label=' ',
                      help="Import a volume or create 'base' phantom manually")
        form.addParam('inputVolume', PointerParam, pointerClass="Volume", label='Input volume',
                      condition="option==0", help="Volume used as 'base' phantom")
        form.addParam('create', StringParam, label='Create phantom', condition="option==1",
                      default='40 40 40 0\ncyl + 1 0 0 0 15 15 2 0 0 0\nsph + 1 0 0 5 2\ncyl + 1 0 0 -5 2 2 10 0 90 0\n'
                              'sph + 1 0 -5 5 2', help="create a phantom description.")
        form.addParam('sampling', FloatParam, label='Sampling rate', default=4)
        form.addParam('mwfilter', BooleanParam, label='Apply missing wedge?', default=False,
                      help='Apply a filter to simulate the missing wedge along Y axis.')
        form.addParam('mwangle', IntParam, label='Missing wedge angle', default=60, condition='mwfilter==True',
                      help='Missing wedge (along y) for data between +- this angle.')

        form.addSection(label='Transform')
        form.addParam('nsubtomos', IntParam, label='Number of subtomograms', default=50,
                      help="How many phantom subtomograms")
        form.addParam('rotmin', IntParam, label='Min rot angle', default=0)
        form.addParam('rotmax', IntParam, label='Max rot angle', default=60)
        form.addParam('tiltmin', IntParam, label='Min tilt angle', default=0)
        form.addParam('tiltmax', IntParam, label='Max tilt angle', default=60)
        form.addParam('psimin', IntParam, label='Min psi angle', default=0)
        form.addParam('psimax', IntParam, label='Max psi angle', default=60)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createPhantomsStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def createPhantomsStep(self):
        mwfilter = self.mwfilter.get()
        fnVol = self._getExtraPath("phantom000.vol")
        if self.option.get() == 0:
            inputVol = self.inputVolume.get()
            fnInVol = inputVol.getLocation()[1]
            dim = inputVol.getDim()
            if mwfilter:
                mwangle = self.mwangle.get()
                self.runJob("xmipp_transform_filter", " --fourier wedge -%d %d 0 0 0 -i %s -o %s"
                            % (mwangle, mwangle, fnInVol, fnVol))
            else:
                self.runJob("xmipp_image_convert", " -i %s -o %s" % (fnInVol, fnVol))
        else:
            desc = self.create.get()
            fnDescr = self._getExtraPath("phantom.descr")
            fhDescr = open(fnDescr, 'w')
            fhDescr.write(desc)
            fhDescr.close()
            dim = [desc.split()[0], desc.split()[1], desc.split()[2]]
            self.runJob("xmipp_phantom_create", " -i %s -o %s" % (fnDescr, fnVol))
            if mwfilter:
                mwangle = self.mwangle.get()
                self.runJob("xmipp_transform_filter", " --fourier wedge -%d %d 0 0 0 -i %s -o %s"
                            % (mwangle, mwangle, fnVol, fnVol))

        self.outputSet = self._createSetOfSubTomograms(self._getOutputSuffix(SetOfSubTomograms))
        self.outputSet.setDim(dim)
        self.outputSet.setSamplingRate(self.sampling.get())

        subtomobase = SubTomogram()
        acq = TomoAcquisition()
        subtomobase.setAcquisition(acq)
        subtomobase.setLocation(fnVol)
        subtomobase.setSamplingRate(self.sampling.get())
        transformBase = Transform()
        transformBase.setMatrix(np.identity(4))
        subtomobase.setTransform(transformBase)
        self.outputSet.append(subtomobase)

        for i in range(int(self.nsubtomos.get())-1):
            fnPhantomi = self._getExtraPath("phantom%03d.vol" % int(i+1))
            rot = np.random.randint(self.rotmin.get(), self.rotmax.get())
            tilt = np.random.randint(self.tiltmin.get(), self.tiltmax.get())
            psi = np.random.randint(self.psimin.get(), self.psimax.get())
            self.runJob("xmipp_transform_geometry", " -i %s -o %s --rotate_volume euler %d %d %d "
                        % (fnVol, fnPhantomi, rot, tilt, psi))

            if mwfilter:
                self.runJob("xmipp_transform_filter", " --fourier wedge -%d %d 0 0 0 -i %s -o %s"
                            % (mwangle, mwangle, fnPhantomi, fnPhantomi))

            subtomo = SubTomogram()
            subtomo.setAcquisition(acq)
            subtomo.setLocation(fnPhantomi)
            subtomo.setSamplingRate(self.sampling.get())
            A = euler_matrix(np.deg2rad(rot), np.deg2rad(tilt), np.deg2rad(psi), 'szyz')
            transform = Transform()
            transform.setMatrix(A)
            subtomo.setTransform(transform)
            self.outputSet.append(subtomo)

    def createOutputStep(self):
        self._defineOutputs(outputSubtomograms=self.outputSet)
        if self.option.get() == 0:
            self._defineSourceRelation(self.inputVolume.get(), self.outputSet)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if self.rotmin.get() >= self.rotmax.get():
            errors.append("rot max must be bigger than rot min")
        if self.tiltmin.get() >= self.tiltmax.get():
            errors.append("tilt max must be bigger than tilt min")
        if self.psimin.get() >= self.psimax.get():
            errors.append("psi max must be bigger than psi min")
        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputSubtomograms'):
            summary.append("Output phantom not ready yet.")
        else:
            summary.append("%s phantoms created with random orientations" % self.nsubtomos.get())
            if self.mwfilter.get():
                summary.append("Missing wedge applied between +-%d along Y axis" % self.mwangle.get())
        return summary

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputSubtomograms'):
            methods.append("Output phantoms not ready yet.")
            return methods
        else:
            methods.append("%s phantoms created with random orientations." % self.nsubtomos.get())
            if self.mwfilter.get():
                methods.append("Missing wedge applied between +-%d along Y axis." % self.mwangle.get())
            return methods
