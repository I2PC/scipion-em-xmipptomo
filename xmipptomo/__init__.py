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
_references = ['delaRosaTrevin2013']
__version__ = "3.1.0"

class Plugin(xmipp3.Plugin):

    @classmethod
    def getTensorFlowEnviron(cls):
        """ Create the needed environment for XmippTomo programs. """
        environ = pwutils.Environ(os.environ)
        environ.update({
            "TF_FORCE_GPU_ALLOW_GROWTH": "'true'"
        }, position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def defineBinaries(cls, env):
        preMsgs = []
        cudaMsgs = []
        nvidiaDriverVer = None
        if os.environ.get('CUDA', 'True') == 'True':
            try:
                nvidiaDriverVer = subprocess.Popen(["nvidia-smi",
                                                    "--query-gpu=driver_version",
                                                    "--format=csv,noheader"],
                                                   env=xmipp3.Plugin.getEnviron(),
                                                   stdout=subprocess.PIPE
                                                   ).stdout.read().decode('utf-8').split(".")[0]
                if int(nvidiaDriverVer) < 390:
                    preMsgs.append("Incompatible driver %s" % nvidiaDriverVer)
                    cudaMsgs.append("Your NVIDIA drivers are too old (<390). "
                                    "Tensorflow was installed without GPU support. "
                                    "Just CPU computations enabled (slow computations).")
                    nvidiaDriverVer = None
            except Exception as e:
                preMsgs.append(str(e))

        if nvidiaDriverVer is not None:
            preMsgs.append("CUDA support find. Driver version: %s" % nvidiaDriverVer)
            msg = "Tensorflow installed with CUDA SUPPORT."
            cudaMsgs.append(msg)
            useGpu = True
        else:
            preMsgs.append("CUDA will NOT be USED. (not found or incompatible)")
            msg = ("Tensorflow installed without GPU. Just CPU computations "
                   "enabled (slow computations). To enable CUDA (drivers>390 needed), "
                   "set CUDA=True in 'scipion.conf' file")
            cudaMsgs.append(msg)
            useGpu = False

        # commands  = [(command, target), (cmd, tgt), ...]
        cmdsInstall = [(cmd, envName + ".yml") for cmd, envName in
                       xmipp3.CondaEnvManager.yieldInstallAllCmds(useGpu=useGpu)]

        now = xmipp3.datetime.now()
        installDLvars = {'modelsUrl': "http://scipion.cnb.csic.es/downloads/scipion/software/em",
                         'syncBin': cls.getHome('bin/xmipp_sync_data'),
                         'modelsDir': cls.getHome('models'),
                         'modelsPrefix': "models_UPDATED_on",
                         'xmippLibToken': 'xmippLibToken',
                         'libXmipp': cls.getHome('lib/libXmipp.so'),
                         'preMsgsStr': ' ; '.join(preMsgs),
                         'afterMsgs': "\n > ".join(cudaMsgs)}

        installDLvars.update({'modelsTarget': "%s_%s_%s_%s"
                                              % (installDLvars['modelsPrefix'],
                                                 now.day, now.month, now.year)})

        modelsDownloadCmd = ("rm %(modelsPrefix)s_* %(xmippLibToken)s 2>/dev/null ; "
                             "echo 'Downloading pre-trained models...' ; "
                             "%(syncBin)s update %(modelsDir)s %(modelsUrl)s DLmodels && "
                             "touch %(modelsTarget)s && echo ' > %(afterMsgs)s'"
                             % installDLvars,  # End of command
                             installDLvars['modelsTarget'])  # Target

        xmippInstallCheck = ("if ls %(libXmipp)s > /dev/null ; "
                             "then touch %(xmippLibToken)s; echo ' > %(preMsgsStr)s' ; "
                             "else echo ; echo ' > Xmipp installation not found, "
                             "please install it first (xmippSrc or xmippBin*).';echo;"
                             " fi" % installDLvars,  # End of command
                             installDLvars['xmippLibToken'])  # Target

        env.addPackage(xmipp3.constants.XMIPP_DLTK_NAME, version='0.2', urlSuffix='external',
                       commands=[xmippInstallCheck] + cmdsInstall + [modelsDownloadCmd],
                       deps=[], tar=xmipp3.constants.XMIPP_DLTK_NAME + '.tgz',
                       default=True)
