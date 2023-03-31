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

    _dirChoices = ['X', 'Y', 'Z']

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='General parameters')
        form.addParam('input', PointerParam, pointerClass=[SetOfSubTomograms, SetOfVolumes],
                      label='Input Volumes', help='This protocol can *not* work with .em files *if* the input is a set'
                                                  ' of tomograms or a set of volumes, ')
        form.addParam('radAvg', BooleanParam, default=False, label='Compute radial average?',
                      help='Compute the radial average with respect to Z of the input volumes and from them, '
                           'it computes their projections in the desired direction')
        form.addParam('dirParam', EnumParam, choices=self._dirChoices, default=2, display=EnumParam.DISPLAY_HLIST,
                      label='Projection direction', condition='radAvg == False')
        form.addParam('rangeParam', EnumParam, choices=['All', 'Range'], default=0, display=EnumParam.DISPLAY_HLIST,
                      label='Range of slices', condition='radAvg == False',
                      help='Range of slices used to compute the projection, where 0 is the central slice.')
        form.addParam('cropParam', IntParam, default=10, label='Slices', condition="rangeParam == 1",
                      help='Crop this amount of voxels in each side of the selected direction.')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):

        self._insertFunctionStep(self.projectStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -------------------------------
    def projectStep(self):
        input = self.input.get()
        x, y, z = input.getDim()

        print("Dimensions (x,y,z) are: %s, %s, %s." % (x,y,z))


        dir = self.dirParam.get()
        partialProjection = self.rangeParam.get() == 1

        if partialProjection:
            bottomSlice = int(x/2 - self.cropParam.get())
            topSlice = int(x/2 + self.cropParam.get())


        fnProj = self._getExtraPath("projections.mrcs")
        lib.createEmptyFile(fnProj, x, y, 1, input.getSize())

        tmpOrientedSubtomoFn = self._getExtraPath("oriented.mrc")

        for subtomo in input.iterItems():

            fn = "%s@%s" % subtomo.getLocation()
            if fn.endswith('.mrc'):
                fn += ':mrc'

            if subtomo.hasTransform():

                ih().rotateVolume(fn, tmpOrientedSubtomoFn, subtomo.getTransform())
                fn = tmpOrientedSubtomoFn

            vol = ih().read(fn)
            img = ih().createImage()

            if self.radAvg.get():
                img = vol.radialAverageAxis()
            else:
                volData = vol.getData()

                proj = np.empty([x, y])

                # X axis
                if dir == 0:
                    if partialProjection:
                        volData = volData[:, :, bottomSlice:topSlice]
                    for zi in range(z):
                        for yi in range(y):
                            proj[yi, zi] = np.sum(volData[zi, yi, :])
                # Y axis
                elif dir == 1:
                    if partialProjection:
                        volData = volData[:, bottomSlice:topSlice, :]
                    for zi in range(z):
                        for xi in range(x):
                            proj[zi, xi] = np.sum(volData[zi, :, xi])
                # Z axis
                elif dir == 2:
                    if partialProjection:
                        volData = volData[bottomSlice:topSlice, :, :]
                    for xi in range(x):
                        for yi in range(y):
                            proj[yi, xi] = np.sum(volData[:, yi, xi])

                # Make the projection to be the image data
                img.setData(proj)

            # Write the image at a specific slice
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

        methods=[]
        if self.radAvg:
            methods.append("Radial average computed for %d volumes with dimensions %s were obtained." % (vols.getSize(), vols.getDimensions()))
        else:
            methods.append("Projections on the %s axis of %d volumes with dimensions %s were obtained using %s slices."
                       % (self._dirChoices[self.dirParam.get()], vols.getSize(), vols.getDimensions(),
                          "all" if self.rangeParam.get() == 0 else self.cropParam.get()))
        return methods


    def _summary(self):
        summary = []

        if self.radAvg:
            summary.append("Radial average: %s" % self.radAvg)
        else:

            summary.append("Projections direction: %s" % self._dirChoices[self.dirParam.get()])
            summary.append("Number: %s" % ("All" if self.rangeParam.get() == 0 else self.cropParam.get()))

        return summary

