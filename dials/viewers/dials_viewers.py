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

import tkinter as tk
import textwrap

from pyworkflow.gui.widgets import Button, HotButton
import pyworkflow.gui as pwgui
import pyworkflow.gui.text as pwtext
import pyworkflow.utils as pwutils
import pyworkflow.viewer as pwviewer
import pyworkflow.protocol.params as params

import dials.protocols as dialsProt
import dials.utils as dutils


class DialsImageView(pwviewer.CommandView):
    def __init__(self, modelFile, reflectionFile=None, **kwargs):

        if reflectionFile is None:
            pwviewer.CommandView.__init__(
                self, 'dials.image_viewer "%s" &' % modelFile, **kwargs)
        else:
            pwviewer.CommandView.__init__(self, 'dials.image_viewer experiments="%s" reflections="%s" &' % (
                modelFile, reflectionFile), **kwargs)


class DialsReciprocalLatticeView(pwviewer.CommandView):
    def __init__(self, modelFile, reflectionFile=None, **kwargs):
        pwviewer.CommandView.__init__(self, 'dials.reciprocal_lattice_viewer "%s" "%s" &' % (
            modelFile, reflectionFile), **kwargs)


IMAGE_VIEWER = 0
RECIPROCAL_VIEWER = 1


class DialsFoundSpotsViewer(pwviewer.ProtocolViewer):
    ''' Vizualisation of Dials imported images and results from spotfinding and indexing '''

    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [dialsProt.DialsProtImportDiffractionImages, dialsProt.DialsProtFindSpots, dialsProt.DialsProtIndexSpots,
                dialsProt.DialsProtRefineSpots, dialsProt.DialsProtIntegrateSpots]

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


class DialsResultsViewer(pwviewer.Viewer):
    """ Viewer for analyzing results from dials processing in general  """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [dialsProt.DialsProtImportDiffractionImages,
                dialsProt.DialsProtFindSpots,
                dialsProt.DialsProtIndexSpots,
                dialsProt.DialsProtRefineSpots,
                dialsProt.DialsProtIntegrateSpots,
                dialsProt.DialsProtSymmetry,
                dialsProt.DialsProtScaling,
                dialsProt.DialsProtExport]

    def visualize(self, obj, **kwargs):
        self.resultsWindow = self.tkWindow(ResultsWindow,
                                           title='Results',
                                           protocol=obj
                                           )
        self.resultsWindow.show()


class ResultsWindow(pwgui.Window):

    def __init__(self, **kwargs):
        pwgui.Window.__init__(self, **kwargs)

        self.protocol = kwargs.get('protocol')

        content = tk.Frame(self.root)
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        content.columnconfigure(0, weight=1)

    def _createContent(self, content):
        topFrame = tk.Frame(content)
        content.columnconfigure(0, weight=1)
        topFrame.grid(row=0, column=0, sticky='new', padx=5, pady=5)

        content.rowconfigure(1, weight=1)
        resultsFrame = tk.Frame(content)
        # Descide how much relative space the frames should fill
        resultsFrame.rowconfigure(0, weight=100)
        resultsFrame.columnconfigure(0, weight=2)
        resultsFrame.grid(row=1, column=0, sticky='news', padx=5, pady=5)

        buttonsFrame = tk.Frame(content)
        buttonsFrame.grid(row=2, column=0, sticky='new', padx=5, pady=5)

        self._fillTopFrame(topFrame)
        self._fillResultsFrame(resultsFrame)
        self._fillButtonsFrame(buttonsFrame)

    def _fillTopFrame(self, frame):
        p1 = tk.Label(frame, text='Project: ')
        p1.grid(row=0, column=0, sticky='nw', padx=5, pady=5)
        projName = self.protocol.getProject().getShortName()
        p2 = tk.Label(frame, text=projName, font=self.fontBold)
        p2.grid(row=0, column=1, sticky='nw', padx=5, pady=5)
        dataSource = self.protocol.getDatasets()
        splitSource = dataSource.split('\n')
        expl = splitSource[0]
        source = "\n".join(splitSource[1:])
        p3 = tk.Label(frame, text=expl)
        p3.grid(row=1, column=0, sticky='nw', padx=5, pady=5)
        p4 = tk.Label(frame, text=source)
        p4.grid(row=1, column=1, sticky='nw', padx=5, pady=5)

    def _fillResultsFrame(self, frame):
        try:
            results = self.protocol.getLogOutput()
        except Exception as e:
            results = e
        txt = textwrap.dedent(results)
        # TODO: add monospace font
        s = pwtext.Text(frame)
        s.addText(txt)
        s.configure(state='disabled')
        s.grid(row=0, column=0, sticky='news')

    def _fillButtonsFrame(self, frame):
        subframe = tk.Frame(frame)
        subframe.grid(row=0, column=0, sticky='nw')
        frame.columnconfigure(1, weight=1)

        imgPlainBtn = Button(subframe, "View plain images",
                             command=self._viewPlainImages)
        imgPlainBtn.grid(row=0, column=0, sticky='nw', padx=(0, 5))
        if self._getModel() is None:
            imgPlainBtn['state'] = 'disabled'

        imgOverlaidBtn = Button(subframe, "View images with spots",
                                command=self._viewOverlaidImages)
        imgOverlaidBtn.grid(row=0, column=1, sticky='nw', padx=(0, 5))
        if None in {self._getModel(), self._getRefls()}:
            imgOverlaidBtn['state'] = 'disabled'

        reciprocalBtn = Button(subframe, "View reciprocal lattice",
                               command=self._viewReciprocal)
        reciprocalBtn.grid(row=0, column=2, sticky='nw', padx=(0, 5))
        if None in {self._getModel(), self._getRefls()}:
            reciprocalBtn['state'] = 'disabled'

        htmlBtn = HotButton(subframe, 'Open HTML Report',
                            command=self._openHTML)
        htmlBtn.grid(row=0, column=3, sticky='nw', padx=(0, 5))
        if self._getHtml() is None:
            htmlBtn['state'] = 'disabled'

        closeBtn = self.createCloseButton(frame)
        closeBtn.grid(row=0, column=1, sticky='ne')

    def _getModel(self):
        try:
            return self.protocol.getModelFile()
        except AttributeError:
            pass
        try:
            return self.protocol.getOutputModelFile()
        except AttributeError:
            pass
        try:
            return self.protocol.getInputModelFile()
        except AttributeError:
            pass
        return None

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

    def _getHtml(self):
        try:
            reportPath = self.protocol.getOutputHtmlFile()
        except AttributeError:
            return None

        if pwutils.exists(reportPath):
            return reportPath
        else:
            return None

    def _viewPlainImages(self, e=None):
        DialsImageView(self._getModel(), reflectionFile=None).show()

    def _viewOverlaidImages(self, e=None):
        DialsImageView(self._getModel(),
                       reflectionFile=self._getRefls()).show()

    def _viewReciprocal(self, e=None):
        DialsReciprocalLatticeView(
            self._getModel(), reflectionFile=self._getRefls()).show()

    def _openHTML(self, e=None):
        if self._getHtml() is not None:
            dutils._showHtmlReport(self._getHtml())
        else:
            self.showInfo("Could not find an HTML file to open")
