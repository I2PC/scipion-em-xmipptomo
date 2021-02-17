# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************


import numpy as np

from pwem.objects import Transform

from pyworkflow.protocol import params

from tomo.protocols import ProtTomoPicking
from tomo.objects import SubTomogram, TomoAcquisition


class XmippProtAlignTransform(ProtTomoPicking):
    '''Protocol to rotate a series of alignments to a common reference defined by a
    Subtomogram Average'''

    _label = 'align transformations'

    def __init__(self, **args):
        ProtTomoPicking.__init__(self, **args)
        # self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('firstAverage', params.PointerParam,
                      pointerClass='AverageSubTomogram',
                      label='Reference Subtomogram Average', important=True)
        form.addParam('secondAverage', params.PointerParam,
                      pointerClass='AverageSubTomogram',
                      label='Moving Subtomogram Average', important=True)
        form.addParam('secondSubtomos', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      label="Moving Subtomograms", important=True)

    # --------------------------- INSERT steps functions ---------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('findAlignment')
        self._insertFunctionStep('alignToReference')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ---------------------------
    def findAlignment(self):
        first_average = self.firstAverage.get()
        second_average = self.secondAverage.get()

        # Find transformation needed to align the Transformations to a common axis
        params = ' --i1 %s --i2 %s --local --dontScale ' \
                 '--copyGeo %s' % \
                 (first_average.getFileName(), second_average.getFileName(),
                  self._getExtraPath("AlignMatricesAxis.txt"))
        self.runJob("xmipp_volume_align", params)

    def alignToReference(self):
        self.second_subtomos = self.secondSubtomos.get()

        # Extract Transformation Matrices from input SubTomograms
        second_matrices = self.queryMatrices(self.second_subtomos)

        # Align Transformations matrices to a common axis using previously computed Matrix
        self.alignAxisMatrix = np.loadtxt(self._getExtraPath("AlignMatricesAxis.txt")).reshape(4, 4)
        self.second_matrices = np.asarray([(tr[0], self.alignAxisMatrix @ tr[1]) for tr in second_matrices])

    def createOutputStep(self):
        outSubtomos = self._createSetOfSubTomograms()
        outSubtomos.setSamplingRate(self.second_subtomos.getSamplingRate())
        outSubtomos.setCoordinates3D(self.second_subtomos.getCoordinates3D())
        # if self.second_subtomos.getAcquisition() is not None:
        #     acquisition = TomoAcquisition()
        #     acquisition.setAngleMin(self.second_subtomos.getFirstItem().getAcquisition().getAngleMin())
        #     acquisition.setAngleMax(self.second_subtomos.getFirstItem().getAcquisition().getAngleMax())
        #     acquisition.setStep(self.second_subtomos.getFirstItem().getAcquisition().getStep())
        #     outSubtomos.setAcquisition(acquisition)
        for ids, inSubtomo in enumerate(self.second_subtomos.iterItems()):
            subtomogram = SubTomogram()
            subtomogram.setObjId(self.second_matrices[ids][0])
            subtomogram.setLocation(inSubtomo.getLocation())
            subtomogram.setCoordinate3D(inSubtomo.getCoordinate3D())
            subtomogram.setTransform(Transform(self.second_matrices[ids][1]))
            subtomogram.setVolName(inSubtomo.getVolName())
            outSubtomos.append(subtomogram)
        self._defineOutputs(outputSetOfSubtomogram=outSubtomos)
        self._defineSourceRelation(self.firstAverage, outSubtomos)
        self._defineSourceRelation(self.secondAverage, outSubtomos)

    # --------------------------- UTILS functions ---------------------------
    def queryMatrices(self, subtomos):
        matrices = [(subtomo.getObjId(), subtomo.getTransform().getMatrix())
                    for subtomo in subtomos.iterItems()]
        return matrices

    # --------------------------- INFO functions ---------------------------
    def _summary(self):
        summary = []

        if self.getOutputsSize() >= 1:
            summary.append('A total of *%d Subtomogram Transformations* have been aligned '
                           'to *Subtomogram Average %s*' % (self.secondSubtomos.get().getSize(),
                                                              self.firstAverage.get().getFileName()))
        else:
            summary.append('Output not ready yet')

        return summary