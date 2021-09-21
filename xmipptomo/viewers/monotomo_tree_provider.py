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
from os.path import join, isfile

from pyworkflow.utils import removeBaseExt
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider


class MonoTomoViewerProvider(TomogramsTreeProvider):

    def getObjectInfo(self, tomo):
        print(tomo.get().getFileName())
        tomoName = removeBaseExt(tomo.getFileName())
        return {'key': tomoName, 'parent': None}

    def getColumns(self):
        return [('tomo', 400)]