from pwem.viewers.viewers_data import DataViewer, ObjectView
from pwem.viewers import showj
from xmipptomo.protocols.protocol_phantom_subtomo import XmippProtPhantomSubtomo
class XmippPhantomSubtomoViewer(DataViewer):
    """ Wrapper to visualize Subtomo phantom outputs """
    _label = 'Subtomogram phantom create viewer'
    _targets = [XmippProtPhantomSubtomo]

    def _visualize(self, obj, **args):
        if hasattr(self.protocol, XmippProtPhantomSubtomo._possibleOutputs.outputSubtomograms.name):
            obj = self.protocol.outputSubtomograms
            fn = obj.getFileName()
            labels = 'id enabled comment _filename _transform._matrix phantom_rot phantom_tilt phantom_psi phantom_shiftX phantom_shiftY phantom_shiftZ'
            self._views.append(ObjectView(self._project, obj.strId(), fn,
                                          viewParams={showj.ORDER: labels,
                                                      showj.VISIBLE: labels,
                                                      showj.MODE: showj.MODE_MD,
                                                      showj.RENDER: '_filename'}))
            return self._views