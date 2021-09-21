# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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

import pyworkflow.viewer as pwviewer
import pwem.viewers.views as vi

from .monotomo_tree_provider import MonoTomoViewerProvider
from ..protocols import XmippProtMonoTomo


class MonotomoDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    using pyvista
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        XmippProtMonoTomo
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return vi.ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        localResolutionTomos = self.protocol.outputLocalResolutionSetOfTomograms

        tomoList = []
        pathList = []
        from os.path import basename
        for tom in localResolutionTomos.iterItems():
            tomoList.append(tom.clone())
            #pathList.append(basename(tom.getFileName()))
        #[tom.clone() for tom in localResolutionTomos.iterItems()]
        print('1')
        print(tomoList)
        print(pathList)
        print('2')
        #print(os.path.basename(self.protocol._getExtraPath())
        vesicleProvider = MonoTomoViewerProvider(tomoList, None, None)

        return views