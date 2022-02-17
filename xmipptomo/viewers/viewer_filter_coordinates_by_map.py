# -*- coding: utf-8 -*-
# **************************************************************************
# * Authors:    Jose Luis Vilas (jlvilas@cnb.csic.es)
# * Authors:    Federico P. de Isidro-Gomez (fp.deisidro@cnb.csic.es)
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

import matplotlib.pyplot as plt
from matplotlib import cm
from pyworkflow.gui.dialog import showError

import os

from pyworkflow.protocol.params import (LabelParam, EnumParam,
                                        IntParam, LEVEL_ADVANCED, StringParam)
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pwem.viewers import DataView, EmPlotter
from pwem.viewers.viewer_localres import LocalResolutionViewer
from pwem.emlib.metadata import MetaData, MDL_X, MDL_COUNT

from xmipp3.viewers.plotter import XmippPlotter

from xmipptomo.protocols.protocol_filter_coordinates_by_map import XmippProtFilterCoordinatesByMap


class XmippProtFilterCoordinatesByMapViewer(ProtocolViewer):
    """
    Visualization tools for filtering a set of coordinates by tomogram values

    """
    _label = 'viewer filter coordinates by map'
    _targets = [XmippProtFilterCoordinatesByMap]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, *args, **kwargs):
        self.selectedTomogram = None
        ProtocolViewer.__init__(self, *args, **kwargs)

    def _defineParams(self, form):

        form.addSection(label='Visualization')

        form.addParam('inputTomogram', StringParam,
                      label="Select a tomogram",
                      important=True)

        form.addParam('doShowCreateHistogram', LabelParam,
                      label="Show histogram")


    def _getVisualizeDict(self):
        return {'doShowOriginalVolumeSlices': self._showOriginalVolumeSlices,
                'doShowCreateHistogram': self._showHistogram,
                }

    def _showResolutionSlices(self, param=None):
        cm = DataView(STATISTICS_METADATA)
        return [cm]

    def _showHistogram(self, param=None):
        cm = DataView(STATISTICS_METADATA)
        return [cm]

    def _plotHistogram(self, param=None):

        fn = self.createPath(STATISTICS_METADATA)
        md = MetaData()
        md.read(fn)
        x_axis = []
        y_axis = []

        i = 0
        for idx in md:
            x_axis_ = md.getValue(MDL_X, idx)
            if i == 0:
                x0 = x_axis_
            elif i == 1:
                x1 = x_axis_
            y_axis_ = md.getValue(MDL_COUNT, idx)

            i += 1
            x_axis.append(x_axis_)
            y_axis.append(y_axis_)

        plotter = EmPlotter()
        plotter.createSubPlot("Histogram",
                              "Tomogram Value (a.u.)", "# of Counts")

    barwidth = x1 - x0

    plotter.plotDataBar(x_axis, y_axis, barwidth)

    return [plotter]
