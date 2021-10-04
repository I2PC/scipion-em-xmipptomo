# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# * [1] Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from pwem.wizards import ColorScaleWizardBase
from pyworkflow.wizard import Wizard
from .protocols import XmippProtConnectedComponents
from xmipptomo.viewers import XmippMonoTomoViewer
from .viewers.monotomo_tree_provider import MonotomoViewer


class XmippConnectedCompWizard(Wizard):
    _targets = ([(XmippProtConnectedComponents, ['distance'])])

    def show(self, form):
        tomoCCProt = form.protocol
        inputCoordinates = tomoCCProt.inputCoordinates.get()
        if not inputCoordinates:
            print('You must specify input coordinates')
            return
        boxSize = inputCoordinates.getBoxSize()
        if not boxSize:
            print('These coordinates do not have box size. Please, enter distance manually.')
            return
        distance = boxSize * 3
        form.setVar('distance', distance)


class ColorScaleWizard(ColorScaleWizardBase):
    _targets = ColorScaleWizardBase.defineTargets(XmippMonoTomoViewer)


class SelectTomogramWizard(Wizard):
    _targets = [(XmippMonoTomoViewer, ['inputTomogram'])]

    def _selectTomogram(self, form):
        localResolutionTomos = form.protocol.outputLocalResolutionSetOfTomograms
        view = MonotomoViewer(form.getTkRoot(), form.protocol,
                              localResolutionTomos).show()
        return view

    def show(self, form, *args):
        view = self._selectTomogram(form.protocol)
        if len(view.values):
            tomogram = view.values[0]
            form.protocol.selectedTomogram = tomogram
            tomogramId = (tomogram.getTsId() if tomogram.getTsId() is not None
                          else tomogram.getObjId())
            form.setVar('inputTomogram', str(tomogramId) +
                        "  ("+str(tomogram) + ")")
