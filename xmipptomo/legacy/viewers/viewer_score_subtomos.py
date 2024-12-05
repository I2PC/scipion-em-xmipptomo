# **************************************************************************
# *
# * Authors:  Pablo Conesa (pconesa@cnb.csic.es),
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
from tomo.objects import SetOfSubTomograms
from xmipptomo.protocols.protocol_score_transform import ScoreTransformOutputs
from pyworkflow.viewer import Viewer

from xmipptomo.protocols import XmippProtScoreTransform
from pwem.viewers.plotter import EmPlotter, plt


class XmippTomoScoreSubtomoViewer(Viewer):
    """ Visualize the output of protocol reconstruct swarm """
    _label = 'Score transformation viewer'
    _targets = [XmippProtScoreTransform]

    def _visualize(self, scoreProtocol: XmippProtScoreTransform, **kwargs):
        # Keep input object in case we need to launch
        # a new protocol and set dependencies

        views =[]
        subtomos = getattr(scoreProtocol, ScoreTransformOutputs.Subtomograms.name)

        views.append(self.plotScoreHistogram(subtomos, **kwargs))

        views.append(self.plotAngularDistibution(subtomos,**kwargs))

        return views

    def plotAngularDistibution(self, subtomos:SetOfSubTomograms, **kwargs):

        plotter = EmPlotter(x=1, y=1, windowTitle='Angular distribution')

        plotter.plotAngularDistributionFromSet(subtomos, "Angular distribution")

        return plotter

    def plotScoreHistogram(self, subtomos:SetOfSubTomograms, **kwargs):

        plotter = EmPlotter(x=1, y=1, windowTitle='Angular distance histogram')


        results = subtomos.getUniqueValues(["id", XmippProtScoreTransform.SCORE_ATTR])

        scores = results[XmippProtScoreTransform.SCORE_ATTR]

        plotter.createSubPlot("Angular distance histogram", "angle", "count")

        plotter.plotHist(scores, 30)

        return plotter
