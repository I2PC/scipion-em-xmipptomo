# **************************************************************************
# *
# * Authors:     Scipion Team
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

import pyworkflow.viewer as pwviewer
from pyworkflow.gui import *

from pyworkflow.gui.dialog import ListDialog


class MonoTomoTreeProvider(TreeProvider):
    """ Model class that will retrieve the information from
        Tomogram and  prepare the  columns/rows models required by the TreeDialog GUI.
    """
    COL_TS = 'Tomogram'
    COL_INFO = 'Info'

    ORDER_DICT = {COL_TS: 'id'}

    def __init__(self, parent, protocol, objs):
        self.title = 'Tomograms display'
        self.parent = parent
        self.protocol = protocol
        self.objs = objs
        TreeProvider.__init__(self, sortingColumnName=self.COL_TS)
        self.selectedDict = {}
        self.mapper = protocol.mapper
        self.maxNum = 200

    def getObjects(self):
        # Retrieve all objects of type className
        objects = []

        orderBy = self.ORDER_DICT.get(self.getSortingColumnName(), 'id')
        direction = 'ASC' if self.isSortingAscending() else 'DESC'

        for obj in self.objs.iterItems(orderBy=orderBy, direction=direction):
            item = obj.clone()
            item._allowsSelection = True
            item._parentObject = None
            objects.append(item)

        return objects

    def getColumns(self):
        cols = [
            (self.COL_TS, 100),
            (self.COL_INFO, 350)]
        return cols

    def getObjectInfo(self, obj):
        itemId = obj.getTsId()
        if itemId is None:
            itemId = str(obj.getObjId())

        key = obj.getObjId()
        text = itemId
        values = [str(obj)]
        tags = ''
        opened = True

        item = {
            'key': key, 'text': text,
            'values': tuple(values),
            'open': opened,
            'selected': False,
            'tags': tags
        }
        return item

    def getObjectActions(self, obj):
        pass


class MonoTomoListDialog(ListDialog):
    def __init__(self, parent, title, provider, itemDoubleClick=False, **kwargs):
        self.parent = parent
        self.provider = provider
        self._itemDoubleClick = itemDoubleClick
        ListDialog.__init__(self, parent, title, provider, message=None,
                            allowSelect=True,  cancelButton=True,
                            selectOnDoubleClick=True,
                            **kwargs)

    def body(self, bodyFrame):
        bodyFrame.config()
        gui.configureWeigths(bodyFrame)
        dialogFrame = tk.Frame(bodyFrame)
        dialogFrame.grid(row=0, column=0, sticky='news', padx=5, pady=5)
        dialogFrame.config()
        gui.configureWeigths(dialogFrame, row=1)
        self._createFilterBox(dialogFrame)
        self._col = 0
        self._createTree(dialogFrame)
        self.initial_focus = self.tree
        if self._itemDoubleClick:
            self.tree.itemDoubleClick = self.doubleClickOnItem
        self.tree.focus(self.tree.getFirst())

    def doubleClickOnItem(self, e=None):
       pass


class MonotomoViewer(pwviewer.View):
    """ This class implements a view using Tkinter ListDialog
    and the ImodTreeProvider.
    """
    def __init__(self, parent, protocol, objs,  **kwargs):
        self._tkParent = parent
        self._protocol = protocol
        self._provider = MonoTomoTreeProvider(parent, protocol, objs)
        self.title = self._provider.title

    def show(self):
        return MonoTomoListDialog(self._tkParent, self.title, self._provider)
