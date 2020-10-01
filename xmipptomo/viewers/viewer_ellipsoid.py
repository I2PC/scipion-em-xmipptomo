# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez
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

import numpy as np
from os.path import exists
from pyworkflow.viewer import (ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO)
from pyworkflow.protocol.params import LabelParam
from xmipptomo.protocols import XmippProtFitEllipsoid
import matplotlib.pyplot as plt


class XmippFitEllipsoidViewer(ProtocolViewer):
    """ Visualization of the output of fit_ellipsoid protocol
    """
    _label = 'viewer fit ellipsoid'
    _targets = [XmippProtFitEllipsoid]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('displayEllipsoids', LabelParam, default=False, label='Display adjusted ellipsoids:')

    def _getVisualizeDict(self):
        return {'displayEllipsoids': self._visualizeEllipsoids}

    def _visualizeEllipsoids(self, obj, **args):
        """display a 3D view of the adjusted ellipsoids"""
        views = []
        fn = "%s" % (self.protocol._getExtraPath("ellipsoid.txt"))
        if exists(fn):
            try:
                c = np.empty(10)
                fh = open(fn)
                for i in range(len(c)):
                    nline = fh.readline()
                    nline = nline.rstrip()
                    c[i] = nline
                    print(c[i])

                x = np.linspace(-9, 9, 400)
                y = np.linspace(-5, 5, 400)
                z = np.linspace(-5, 5, 400)

                plotter = plt.contour(x, y, z, (c[0]*x*x + c[1]*y*y + c[2]*z*z + 2*c[3]*x*y + 2*c[4]*x*z + 2*c[5]*y*z
                                                + 2*c[6]*x + 2*c[7]*y + 2*c[8]*z + c[9]))
                print("------g-------")
                views.append(plotter)
                print("------h------")
            except Exception as e:
                print(e)
        return views
