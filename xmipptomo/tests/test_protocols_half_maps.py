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


import numpy as np

from pwem.emlib.image import ImageHandler
from pyworkflow.tests import BaseTest, setupTestProject

from xmipptomo.protocols import XmippProtPhantomSubtomo, XmippProtHalfMapsSubtomo


class TestXmipptomoHalfMaps(BaseTest):
    """This class check if the protocol to create half maps from a SetOfSubtomograms
    works properly"""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _phantom(self):
        phantom = self.newProtocol(XmippProtPhantomSubtomo, option=1)
        self.launchProtocol(phantom)
        self.assertIsNotNone(phantom.outputSubtomograms,
                             "There was a problem with subtomograms output")
        return phantom

    def _half_maps(self, subtomos):
        half_maps = self.newProtocol(XmippProtHalfMapsSubtomo, inputSubtomograms=subtomos)
        self.launchProtocol(half_maps)
        self.assertIsNotNone(half_maps.halfMaps,
                             "There was a problem with half maps output")
        return half_maps

    def test_half_maps(self):
        phantom = self._phantom()
        half_maps = self._half_maps(phantom.outputSubtomograms)
        map_even = np.squeeze(ImageHandler().read(half_maps.halfMaps[1].getFileName()).getData())
        map_odd = np.squeeze(ImageHandler().read(half_maps.halfMaps[2].getFileName()).getData())
        error = np.sqrt(np.sum((map_even - map_odd) ** 2) / map_even.size)
        self.assertAlmostEqual(error, 0.0, delta=0.1, msg="Unexpected half maps")
        return half_maps
