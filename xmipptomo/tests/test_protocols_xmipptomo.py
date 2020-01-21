# **************************************************************************
# *
# * Authors:    Estrella Fernandez Gimenez [1]
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
from tomo.protocols import ProtImportSubTomograms
from tomo.tests import DataSet
from xmipptomo.protocols import XmippProtConnectedComponents, XmippProtUndoAlignSubtomo, \
    XmippProtApplyTransformSubtomo, XmippProtRoiIJ, XmippProtSubtomoMapBack, XmippProtSubtomoProject

# class TestXmippProtProjectZ(TestXmippBase):
#     """This class check if the protocol project top in Xmipp works properly."""
#
#     @classmethod
#     def setUpClass(cls):
#         setupTestProject(cls)
#         cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
#
#     def importVolumes(self):
#         args = {'filesPath': self.dsXmipp.getFile('volumes/'),
#                 'filesPattern': '*mrc',
#                 'samplingRate': 1
#                 }
#         prot = self.newProtocol(ProtImportVolumes, **args)
#         prot.setObjLabel('import volume')
#         self.launchProtocol(prot)
#
#         return prot
#
#     def _createProjZ(self, inputVolumes):
#         XmippProtSubtomoProject = importFromPlugin("xmipp3.protocols",
#                                                    "XmippProtSubtomoProject",
#                                                    errorMsg=TOMO_IMPORT_ERROR,
#                                                    doRaise=True)
#         prot = self.newProtocol(XmippProtSubtomoProject)
#         prot.input.set(inputVolumes)
#         self.launchProtocol(prot)
#
#         return prot
#
#     def test_top(self):
#         protImportVols = self.importVolumes()
#         prot = self._createProjZ(protImportVols.outputVolumes)