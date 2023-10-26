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
        halfMaps = self.newProtocol(XmippProtHalfMapsSubtomo, inputSubtomograms=subtomos)
        self.launchProtocol(halfMaps)
        self.assertIsNotNone(halfMaps.halfMaps,
                             "There was a problem with half maps output")
        return halfMaps

    def test_half_maps(self):
        phantom = self._phantom()
        halfMaps = self._half_maps(phantom.outputSubtomograms)
        mapEven = np.squeeze(ImageHandler().read(halfMaps.halfMaps[1].getFileName()).getData())
        mapOdd = np.squeeze(ImageHandler().read(halfMaps.halfMaps[2].getFileName()).getData())
        error = np.sqrt(np.sum((mapEven - mapOdd) ** 2) / mapEven.size)
        self.assertAlmostEqual(error, 0.0, delta=0.1, msg="Unexpected half maps")
        return halfMaps
