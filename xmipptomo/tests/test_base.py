# **************************************************************************
# *
# * Authors:    J.L. Vilas (jlvilas@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from tomo.protocols.protocol_import_coordinates import IMPORT_FROM_EMAN
import tomo.protocols


class XmippBaseTest(BaseTest):
    def _runImportCoordinatesAndTomograms(self):
        protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
                                                   auto=IMPORT_FROM_EMAN,
                                                   filesPath=self.coords3D,
                                                   importTomograms=protImportTomogram.Tomograms,
                                                   filesPattern='', boxSize=32,
                                                   samplingRate=5)

        self.launchProtocol(protImportCoordinates3d)
        self.assertSetSize(protImportTomogram.Tomograms, 1,
                           "There was a problem with tomogram output")
        self.assertSetSize(protImportCoordinates3d.outputCoordinates, 5,
                           "There was a problem with coordinates 3d output")

        return protImportTomogram, protImportCoordinates3d
