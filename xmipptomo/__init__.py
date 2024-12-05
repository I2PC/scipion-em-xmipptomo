# **************************************************************************
# *
# * Authors:     you (you@yourinstitution.email)
# *
# * your institution
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
import xmipp3
import subprocess, os
import pyworkflow.utils as pwutils

_logo = "xmipp_logo.png"
_references = ['delaRosaTrevin2013', 'Jimenez2022']
__version__ = "3.24.06.2" #X.YY.MM.sv
        # X.Y.M = version of the xmipp release associated.
        # sv = Set this to ".0" on each xmipp  release.
        # For not release version (hotfix) increase it --> ".1", ".2", ...

class Plugin(xmipp3.Plugin):

    @classmethod
    def getTensorFlowEnviron(cls):
        """ Create the needed environment for XmippTomo programs. """
        environ = xmipp3.Plugin.getEnviron()
        environ.update({
            "TF_FORCE_GPU_ALLOW_GROWTH": "'true'"
        }, position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def defineBinaries(cls, env):

        pass
