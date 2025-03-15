# **************************************************************************
# *
# * Authors:    Federico P. de Isidro-Gomez
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
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


from pyworkflow.tests import BaseTest, setupTestProject

from xmipptomo.protocols import XmippProtPhantomTomo

from xmipptomo.protocols.protocol_phantom_tomo import OutputPhantomTomos


class TestXmippTomoPhantom(BaseTest):
    """This class check if the protocol create phantom tomograms works properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def test_phantom(self):
        phantom = self.newProtocol(XmippProtPhantomTomo,
                                   ntomos=2,
                                   nparticles=10)
        self.launchProtocol(phantom)
        tomograms = getattr(phantom, OutputPhantomTomos.tomograms.name)
        self.assertSetSize(tomograms, 2, "There was a problem with subtomograms output")

        self.assertIsNotNone(tomograms.getAcquisition().getMagnification(),
                             "Acquisition magnification not populated properly")

        coordinates = getattr(phantom, OutputPhantomTomos.coordinates3D.name)
        self.assertSetSize(coordinates, 20, "There was a problem with subtomograms output")

        return phantom
