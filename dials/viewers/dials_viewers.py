# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

from pyworkflow.viewer import CommandView, Viewer, DESKTOP_TKINTER
from pwed.objects import SetOfDiffractionImages, SetOfSpots
from ..protocols import DialsProtFindSpots

# Base on https://github.com/scipion-em/scipion-em-bsoft/blob/devel/bsoft/viewers/bsoft_viewers.py
# ProtocolView in Gautomat

class DialsImageView(CommandView):
    def __init__(self, modelFile, reflectionFile=None, **kwargs):
        
        if reflectionFile is None:
            CommandView.__init__(self, 'dials.image_viewer "%s" &' % modelFile, **kwargs)
        else:
            CommandView.__init__(self, 'dials.image_viewer "%s" "%s" &' % (modelFile, reflectionFile), **kwargs)

class DialsReciprocalLatticeView(CommandView):
    def __init__(self, modelFile, reflectionFile=None, **kwargs):
        CommandView.__init__(self, 'dials.reciprocal_lattice_viewer "%s" "%s" &' % (modelFile, reflectionFile), **kwargs)

class DialsImageViewer(Viewer):
    _environments = [DESKTOP_TKINTER]
    _targets = [DialsProtFindSpots]
    
    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)

    def visualize(self, obj, **kwargs):
        cls = type(obj)

        modelFn = obj.getModelFile()
        DialsImageView(modelFn).show()
