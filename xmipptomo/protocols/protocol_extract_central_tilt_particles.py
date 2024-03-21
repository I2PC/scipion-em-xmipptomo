# **************************************************************************
# *
# * Authors:    Oier Lauzirika Zarrabeitia (olauzirika@cnb.csic.es)
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


from pyworkflow import BETA

from pyworkflow.protocol import params

from pwem.protocols import EMProtocol
from pwem.objects import SetOfParticles, Particle
from pwem import emlib

from tomo.protocols import ProtTomoBase
from xmipptomo.objects import SetOfTiltSeriesParticle, TiltSeriesParticle, TiltParticle

import math

class ExtractCentralTiltParticles(ProtTomoBase, EMProtocol):
    '''Protocol to extract the central image of a set of tilt series particles
    '''

    _label = 'extract central tilt particles'
    _devStatus = BETA

    def __init__(self, *args, **kwargs):
        EMProtocol.__init__(self, *args, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass=SetOfTiltSeriesParticle,
                      label='Input particles', important=True)

    # --------------------------- INSERT steps functions ---------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ---------------------------
    def createOutputStep(self):
        inputParticles: SetOfTiltSeriesParticle = self.inputParticles.get()
        outputParticles: SetOfParticles = self._createSetOfParticles()
        outputParticles.enableAppend()
        outputParticles.setAcquisition(outputParticles.getAcquisition())
        outputParticles.setSamplingRate(inputParticles.getSamplingRate())
        
        tiltSeriesParticle: TiltSeriesParticle
        for tiltSeriesParticle in inputParticles:
            particle = Particle(objId=tiltSeriesParticle.getObjId())

            smallestTiltAngle = math.inf
            tiltParticle: TiltParticle
            for tiltParticle in tiltSeriesParticle:
                absoluteTiltAngle = abs(tiltParticle.getTiltAngle())
                
                if absoluteTiltAngle < smallestTiltAngle:
                    smallestTiltAngle = absoluteTiltAngle
                    
                    particle.setLocation(tiltParticle.getLocation())
                    particle.setCTF(tiltParticle.getCTF())
                    
            outputParticles.append(particle)
                    
        self._defineOutputs(output=outputParticles)  
        self._defineTransformRelation(inputParticles, outputParticles)
        
    # --------------------------- UTILS functions ---------------------------


    # --------------------------- INFO functions ---------------------------
