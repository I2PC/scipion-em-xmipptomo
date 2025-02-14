# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
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
from pwem.wizards import ColorScaleWizardBase
from pyworkflow.gui.dialog import showError

import os

from pyworkflow.protocol.params import (LabelParam, EnumParam,
                                        IntParam, LEVEL_ADVANCED, StringParam)
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pwem.viewers import ChimeraView, DataView, EmPlotter
from pwem.viewers.viewer_localres import LocalResolutionViewer
from pwem.emlib.metadata import MetaData, MDL_X, MDL_COUNT
from pwem.constants import AX_Z

from xmipp3.viewers.plotter import XmippPlotter

from xmipptomo.protocols.protocol_resolution_local_monotomo import (XmippProtMonoTomo,
                                                           TOMOGRAM_RESOLUTION_FILE, FULL_TOMOGRAM_FILE,
                                                           HISTOGRAM_RESOLUTION_FILE, MRCEXT, XMDEXT)


class XmippMonoTomoViewer(LocalResolutionViewer):
    """
    Visualization tools for MonoTomo results.

    MonoTomo is a Xmipp package for estimating the local resolution of 
    tomograms, primarily by cryo-electrons tomography (cryo-EM).
    """
    _label = 'viewer MonoTomo'
    _targets = [XmippProtMonoTomo]
    _environments = [DESKTOP_TKINTER]

    @staticmethod
    def getColorMapChoices():
        return plt.colormaps()

    def __init__(self, *args, **kwargs):
        self.selectedTomogram = None
        ProtocolViewer.__init__(self, *args, **kwargs)

    def _defineParams(self, form):

        form.addSection(label='Visualization')

        form.addParam('inputTomogram', StringParam,
                      label="Select a tomogram",
                      important=True)

        form.addParam('doShowResolutionSlices', LabelParam,
                      label="Show resolution slices")

        form.addParam('doShowOriginalVolumeSlices', LabelParam,
                      label="Show original volume slices")

        form.addParam('doShowResHistogram', LabelParam,
                      label="Show resolution histogram")

        group = form.addGroup('Colored resolution Slices and Volumes')

        group.addParam('sliceAxis', EnumParam, default=AX_Z,
                       choices=['x', 'y', 'z'],
                       display=EnumParam.DISPLAY_HLIST,
                       label='Slice axis')

        group.addParam('doShowVolumeColorSlices', LabelParam,
                       label="Show colored resolution slices")

        group.addParam('doShowOneColorslice', LabelParam,
                       expertLevel=LEVEL_ADVANCED,
                       label='Show selected slice')
        group.addParam('sliceNumber', IntParam, default=-1,
                       expertLevel=LEVEL_ADVANCED,
                       label='Show slice number')

        group.addParam('doShowChimera', LabelParam,
                       label="Show Resolution Tomogram in Chimera")

        ColorScaleWizardBase.defineColorScaleParams(group, defaultHighest=self.protocol.max_res_init,
                                                    defaultLowest=self.protocol.min_res_init)

    def _getVisualizeDict(self):
        return {'doShowOriginalVolumeSlices': self._showOriginalVolumeSlices,
                'doShowResolutionSlices': self._showResolutionSlices,
                'doShowVolumeColorSlices': self._showVolumeColorSlicesResolution,
                'doShowOneColorslice': self._showOneColorslice,
                'doShowResHistogram': self._plotHistogram,
                'doShowChimera': self._showChimera,
                }

    def _showResolutionSlices(self, param=None):
        if self.validateTomogram():
            cm = DataView(self.createPath(TOMOGRAM_RESOLUTION_FILE, MRCEXT))
            return [cm]

    def _showOriginalVolumeSlices(self, param=None):
        if self.validateTomogram():
            cm = DataView(self.createPath(FULL_TOMOGRAM_FILE, MRCEXT))
            return [cm]

    def _showVolumeColorSlicesResolution(self, param=None):
        if self.validateTomogram():
            self._showVolumeColorSlices(self.createPath(TOMOGRAM_RESOLUTION_FILE, MRCEXT))

    def _showVolumeColorSlices(self, mapFile):
        imgData, min_Res, max_Res, voldim__ = self.getImgData(mapFile)

        xplotter = XmippPlotter(x=2, y=2, mainTitle="Local Resolution Slices "
                                                    "along %s-axis."
                                                    % self._getAxis())

        # The slices to be shown are close to the center. Volume size is divided in
        # 9 segments, the fouth central ones are selected i.e. 3,4,5,6
        for i in range(3, 7):
            sliceNumber = self.getSlice(i, imgData)
            plot = self._createSlicePlot(imgData, sliceNumber, xplotter)
        xplotter.getColorBar(plot)

        return [plt.show()]

    def _showOneColorslice(self, param=None):
        imageFile = self.createPath(TOMOGRAM_RESOLUTION_FILE, MRCEXT)
        imgData, _, _, volDim  = self.getImgData(imageFile)

        xplotter = XmippPlotter(x=1, y=1, mainTitle="Local Resolution Slices "
                                                    "along %s-axis."
                                                    % self._getAxis())
        sliceNumber = self.sliceNumber.get()
        if sliceNumber < 0:
            sliceNumber = int(volDim[0] / 2)
        else:
            sliceNumber -= 1
        # sliceNumber has no sense to start in zero
        plot = self._createSlicePlot(imgData, sliceNumber, xplotter)
        xplotter.getColorBar(plot)

        return [plt.show()]

    def _createSlicePlot(self, imgData, sliceNumber, xplotter):
        a = xplotter.createSubPlot("Slice %s" % (sliceNumber + 1), '', '')
        matrix = self.getSliceImage(imgData, sliceNumber, self._getAxis())
        plot = xplotter.plotMatrix(a, matrix, self.lowest.get(), self.highest.get(),
                                   cmap=self.getColorMap(),
                                   interpolation="nearest")
        return plot

    def _plotHistogram(self, param=None):
        if self.validateTomogram():
            fn = self.createPath(HISTOGRAM_RESOLUTION_FILE, XMDEXT)
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
            plotter.createSubPlot("Resolutions Histogram",
                                  "Resolution (A)", "# of Counts")

        barwidth = x1 - x0

        plotter.plotDataBar(x_axis, y_axis, barwidth)

        return [plotter]

    def _getIdFromString(self, tomoDataStr):
        '''
        This function takes the String with the information of the tomogram, like
        '25 Tomogram (1000x1000x200)' and returns the id = 25
        '''
        listInfo = str(tomoDataStr).split()
        return listInfo[0]

    def _getAxis(self):
        return self.getEnumText('sliceAxis')

    def _showChimera(self, param=None):
        if self.validateTomogram():
            fnResVol = self.createPath(TOMOGRAM_RESOLUTION_FILE, MRCEXT)

            fnOrigMap = self.createPath(FULL_TOMOGRAM_FILE, MRCEXT)
            sampRate = 1 #sampling is not set in the binary

            cmdFile = self.protocol._getExtraPath('chimera_resolution_map.py')
            self.createChimeraScript(cmdFile, fnResVol, fnOrigMap, sampRate,
                                     numColors=self.intervals.get(),
                                     lowResLimit=self.highest.get(),
                                     highResLimit=self.lowest.get())
            view = ChimeraView(cmdFile)
            return [view]

    def getColorMap(self):
        cmap = cm.get_cmap(self.colorMap.get())
        if cmap is None:
            cmap = cm.jet
        return cmap

    def createPath(self, fn, ext):
        '''
        This function returns the output path corresponding to a tomogram id, tomId, with
        a filename (fn) and an extension, ext
        '''

        tomId = self._getIdFromString(self.inputTomogram)
        return self.protocol.createOutputPath(fn, tomId, ext)

    def validateTomogram(self):
        '''
        This function returns true or false if there is a tomogram selected
        or not.
        '''
        if not self.inputTomogram.empty():
            return True
        else:
            showError("Tomogram missing", "You need to select one tomogram.",
                      self.getTkRoot())
            return False
