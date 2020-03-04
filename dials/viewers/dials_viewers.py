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

import pyworkflow.protocol.params as params
from pyworkflow.viewer import CommandView, ProtocolViewer, DESKTOP_TKINTER
from ..protocols import DialsProtImportDiffractionImages, DialsProtFindSpots, DialsProtIndexSpots, DialsProtRefineSpots, DialsProtIntegrateSpots

# Base on https://github.com/scipion-em/scipion-em-bsoft/blob/devel/bsoft/viewers/bsoft_viewers.py
# ProtocolView in Gautomat


class DialsImageView(CommandView):
    def __init__(self, modelFile, reflectionFile=None, **kwargs):

        if reflectionFile is None:
            CommandView.__init__(
                self, 'dials.image_viewer "%s" &' % modelFile, **kwargs)
        else:
            CommandView.__init__(self, 'dials.image_viewer "%s" "%s" &' % (
                modelFile, reflectionFile), **kwargs)


class DialsReciprocalLatticeView(CommandView):
    def __init__(self, modelFile, reflectionFile=None, **kwargs):
        CommandView.__init__(self, 'dials.reciprocal_lattice_viewer "%s" "%s" &' % (
            modelFile, reflectionFile), **kwargs)


IMAGE_VIEWER = 0
RECIPROCAL_VIEWER = 1


class DialsFoundSpotsViewer(ProtocolViewer):
    ''' Vizualisation of Dials imported images and results from spotfinding and indexing '''

    _environments = [DESKTOP_TKINTER]
    _targets = [DialsProtImportDiffractionImages, DialsProtFindSpots, DialsProtIndexSpots,
                DialsProtRefineSpots, DialsProtIntegrateSpots]

    def _defineParams(self, form):
        form.addSection(label='Pick viewer')

        form.addParam('viewSelection', params.EnumParam, choices=['image viewer', 'reciprocal lattice viewer'],
                      default=IMAGE_VIEWER, label='Display data with', display=params.EnumParam.DISPLAY_HLIST,
                      help='*image viewer*: Display the images used in spotfinding. Option to show the found spots on the images\n'
                      '*reciprocal viewer*: View the found spots in reciprocal space. Does not work if the protocol is importing diffraction images.')

        form.addParam('viewSpotsOnImage', params.BooleanParam,
                      default=True, label='View the found spots on the images?', condition='viewSelection==%d' % IMAGE_VIEWER)

    def _getVisualizeDict(self):
        visualizeDict = {
            'viewSelection': self._viewerSelection
        }
        return visualizeDict

    def _viewerSelection(self, paramName=None):
        if self.viewSelection == IMAGE_VIEWER:
            return self._viewImages()
        elif self.viewSelection == RECIPROCAL_VIEWER:
            return self._viewReciprocal()

    def _viewImages(self, reflFn=None, **kwargs):
        if self.viewSpotsOnImage:
            reflFn = self._getRefls()
        DialsImageView(self._getModel(), reflectionFile=reflFn).show()

    def _viewReciprocal(self, **kwargs):
        DialsReciprocalLatticeView(
            self._getModel(), reflectionFile=self._getRefls()).show()

    def _getModel(self):
        try:
            return self.protocol.getModelFile()
        except AttributeError:
            try:
                return self.protocol.getOutputModelFile()
            except AttributeError:
                return self.protocol.getInputModelFile()

    def _getRefls(self):
        try:
            return self.protocol.getReflFile()
        except AttributeError:
            pass
        try:
            return self.protocol.getOutputReflFile()
        except AttributeError:
            pass
        try:
            return self.protocol.getInputReflFile()
        except AttributeError:
            pass
        return None
