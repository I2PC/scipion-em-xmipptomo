# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *              Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
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
import enum

import numpy as np
from pwem.emlib.image import ImageHandler as ih
from pwem.emlib import lib
from pwem.objects import Particle, Volume, Transform, String, SetOfVolumes, SetOfParticles
from pwem.protocols import ProtAnalysis3D
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, EnumParam, IntParam, BooleanParam
from tomo.objects import SetOfSubTomograms

class SubtomoProjectOutput(enum.Enum):
    particles = SetOfParticles
    average = Particle

class XmippProtSubtomoProject(ProtAnalysis3D):
    """
    Project a set of volumes or subtomograms to obtain their X, Y or Z projection of the desired range of slices.
    """
    _label = 'subtomo projection'
    _devStatus = BETA
    _possibleOutputs = SubtomoProjectOutput

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='General parameters')
        form.addParam('input', PointerParam, pointerClass=[SetOfSubTomograms, SetOfVolumes],
                      label='Input Volumes', help='This protocol can *not* work with .em files *if* the input is a set'
                                                  ' of tomograms or a set of volumes, ')
        form.addParam('radAvg', BooleanParam, default=False, label='Compute radial average?',
                      help='Compute the radial average with respect to Z of the input volumes and from them, '
                           'it computes their projections in the desired direction')
        form.addParam('dirParam', EnumParam, choices=['X', 'Y', 'Z'], default=2, display=EnumParam.DISPLAY_HLIST,
                      label='Projection direction', condition='radAvg == False')
        form.addParam('rangeParam', EnumParam, choices=['All', 'Range'], default=0, display=EnumParam.DISPLAY_HLIST,
                      label='Range of slices', condition='radAvg == False',
                      help='Range of slices used to compute the projection, where 0 is the central slice.')
        form.addParam('cropParam', IntParam, default=10, label='Slices', condition="rangeParam == 1",
                      help='Crop this amount of voxels in each side of the selected direction.')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):

        self._insertFunctionStep('projectZStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def projectZStep(self):
        input = self.input.get()
        x, y, z = input.getDim()
        dir = self.dirParam.get()
        if self.rangeParam.get() == 1:
            cropParam = self.cropParam.get()

        fnProj = self._getExtraPath("projections.mrcs")
        lib.createEmptyFile(fnProj, x, y, 1, input.getSize())

        for subtomo in input.iterItems():
            fn = subtomo.getLocation()
            if fn[1].endswith('.mrc'):
                fn = list(fn)
                fn[1] += ':mrc'
                fn = tuple(fn)
                subtomo.setFileName(fn[1])
            vol = Volume()
            vol.setLocation('%d@%s' % fn)
            vol = ih().read(vol.getLocation())
            img = ih().createImage()
            if self.radAvg.get():
                img = vol.radialAverageAxis()
            else:
                volData = vol.getData()
                proj = np.empty([x, y])
                if dir == 0:
                    if self.rangeParam.get() == 1:
                        volData = volData[:, :, int(x/2 - cropParam):int(x/2 + cropParam):1]
                    for zi in range(z):
                        for yi in range(y):
                            proj[zi, yi] = np.sum(volData[zi, yi, :])
                elif dir == 1:
                    if self.rangeParam.get() == 1:
                        volData = volData[:, int(x/2 - cropParam):int(x/2 + cropParam):1, :]
                    for zi in range(z):
                        for xi in range(x):
                            proj[zi, xi] = np.sum(volData[zi, :, xi])
                elif dir == 2:
                    if self.rangeParam.get() == 1:
                        volData = volData[int(x/2 - cropParam):int(x/2 + cropParam):1, :, :]
                    for xi in range(x):
                        for yi in range(y):
                            proj[xi, yi] = np.sum(volData[:, yi, xi])
                img.setData(proj)

            img.write('%d@%s' % (subtomo.getObjId(), fnProj))

    def createOutputStep(self):
        input = self.input.get()
        imgSetOut = self._createSetOfParticles()
        imgSetOut.setSamplingRate(input.getSamplingRate())
        imgSetOut.setAlignmentProj()

        # Input could be SetOfVolumes or SetOfSubtomograms
        for item in input.iterItems():
            idx = item.getObjId()
            p = Particle()
            p.setLocation(ih._convertToLocation((idx, self._getExtraPath("projections.mrcs"))))
            p._subtomogramID = String(idx)

            if item.hasTransform():
                transform = Transform()
                transform.setMatrix(item.getTransform().getMatrix())
                p.setTransform(transform)
            imgSetOut.append(p)

        imgSetOut.setObjComment(self.getSummary(imgSetOut))

        self._defineOutputs(**{SubtomoProjectOutput.particles.name:imgSetOut})
        self._defineSourceRelation(self.input, imgSetOut)

        if self.radAvg.get():
            avgFile = self._getExtraPath("average.xmp")
            imgh = ih()
            avgImage = imgh.computeAverage(imgSetOut)
            avgImage.write(avgFile)
            avg = Particle()
            avg.setLocation(1, avgFile)
            avg.copyInfo(imgSetOut)
            self._defineOutputs(**{SubtomoProjectOutput.average.name:avg})
            self._defineSourceRelation(self.input, avg)

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        vols = self.input.get()
        return ["Projection of %d volumes with dimensions %s obtained."
                % (vols.getSize(), vols.getDimensions())]

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Output views not ready yet.")

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                summary.append("*%s:* \n %s " % (key, output.getObjComment()))
        return summary

    def getSummary(self, imgSetOut):
        summary = []
        summary.append("Number of projections generated: %s" % imgSetOut.getSize())
        return "\n".join(summary)
