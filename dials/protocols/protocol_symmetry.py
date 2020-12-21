# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              V. E.G: Bengtsson (viktor.bengtsson@mmk.su.se) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Department of Materials and Environmental Chemistry, Stockholm University
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

import os
import re
from glob import glob
from pathlib import Path

import pyworkflow.protocol as pwprot
import dials.utils as dutils

from pwed.objects import IndexedSpot, SetOfIndexedSpots
from pwed.protocols import EdBaseProtocol
from dials.convert import writeJson, readRefl, writeRefl, copyDialsFile
from dials.constants import *


class DialsProtSymmetry(EdBaseProtocol):
    """ Protocol for checking symmetry of integrated spots using the POINTLESS algorithm as implemented in DIALS
    """
    _label = 'symmetry'

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        # EdProtIndexSpots._defineParams(self, form)

        # The start of the actually relevant part.
        form.addSection(label='Input')

        form.addParam('inputSet', pwprot.PointerParam,
                      pointerClass='SetOfIndexedSpots',
                      label="Indexed spots to symmetry check",
                      help="")

        # Allow adding anything else with command line syntax
        group = form.addGroup('Raw command line input parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED)
        group.addParam('commandLineInput', pwprot.StringParam,
                       default='',
                       help="Anything added here will be added at the end of the command line")

        # Add a section for creating an html report
        form.addSection('HTML report')
        form.addParam('makeReport', pwprot.BooleanParam,
                      label='Do you want to create an HTML report for the output?', default=False,
                      help="",
                      )

        form.addParam('showReport', pwprot.BooleanParam,
                      label='Do you want to open the report as soon as the protocol is done?', default=False,
                      help="",
                      condition="makeReport",
                      )

        group = form.addGroup('Parameters',
                              condition="makeReport",)

        self.extDepOptions = ['remote', 'local', 'embed']
        group.addParam('externalDependencies', pwprot.EnumParam,
                       label='External dependencies: ',
                       choices=self.extDepOptions,
                       default=REMOTE,
                       help="Whether to use remote external dependencies (files relocatable but requires an internet connection), local (does not require internet connection but files may not be relocatable) or embed all external dependencies (inflates the html file size).",
                       )

        group.addParam('pixelsPerBin', pwprot.IntParam,
                       label='Pixels per bin',
                       default=40,
                       GE=1,
                       allowsNull=True,
                       )

        group.addParam('centroidDiffMax', pwprot.FloatParam,
                       label='Centroid diff max',
                       default=None,
                       allowsNull=True,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       help="Magnitude in pixels of shifts mapped to the extreme colours in the heatmap plots centroid_diff_x and centroid_diff_y",
                       )

        # Allow adding anything else with command line syntax
        group = form.addGroup('HTML report command line parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED,
                              condition="makeReport",)
        group.addParam('commandLineInputReport', pwprot.StringParam,
                       default='',
                       help="Anything added here will be added at the end of the command line")

   # -------------------------- INSERT functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(
            'convertInputStep', self.inputSet.getObjId())
        self._insertFunctionStep('symmetryStep')
        if self.makeReport:
            self._insertFunctionStep('makeHtmlReportStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, inputSpotId):
        inputSet = self.inputSet.get()
        self.info("Number of images: %s" % inputSet.getSize())
        self.info("Number of spots: %s" % inputSet.getSpots())
        # Write new model and/or reflection file if no was supplied from set
        if self._checkWriteModel():
            self.writeJson(inputSet, self.getInputModelFile())
        if self._checkWriteRefl():
            self.writeRefl(inputSet, self.getInputReflFile())

    def symmetryStep(self):
        program = 'dials.symmetry'
        arguments = self._prepCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    def makeHtmlReportStep(self):
        prog = 'dials.report'
        arguments = self._prepCommandlineReport()
        self.runJob(prog, arguments)
        if self.showReport:
            dutils._showHtmlReport(self.getOutputHtmlFile())

    def createOutputStep(self):
        # Check that the indexing created proper output
        assert(os.path.exists(self.getOutputReflFile()))
        assert(os.path.exists(self.getOutputModelFile()))

        outputSet = self._createSetOfIndexedSpots()
        outputSet.setDialsModel(self.getOutputModelFile())
        outputSet.setDialsRefl(self.getOutputReflFile())

        outputSet.write()

        self._defineOutputs(outputSymmetrizedSpots=outputSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []

        if self.getDatasets() not in (None, ''):
            summary.append(self.getDatasets())
            summary.append("\n")

        if self.getLogOutput() not in (None, ''):
            summary.append(self.getLogOutput())

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputModelFile(self):
        if self.getSetModel():
            return self.getSetModel()
        else:
            return self._getExtraPath('integrated.expt')

    def getInputReflFile(self):
        if self.getSetRefl():
            return self.getSetRefl()
        else:
            return self._getExtraPath('integrated.refl')

    def getOutputModelFile(self):
        return self._getExtraPath('symmetrized.expt')

    def getOutputReflFile(self):
        return self._getExtraPath('symmetrized.refl')

    def getOutputHtmlFile(self):
        return self._getExtraPath('dials.report.html')

    def getPhilPath(self):
        return self._getTmpPath('symmetry.phil')

    def existsPath(self, path):
        return os.path.exists(path)

    def getSetModel(self):
        inputSet = self.inputSet.get()
        if self.existsPath(inputSet.getDialsModel()):
            return inputSet.getDialsModel()
        else:
            return None

    def getSetRefl(self):
        inputSet = self.inputSet.get()
        if self.existsPath(inputSet.getDialsRefl()):
            return inputSet.getDialsRefl()
        else:
            return None

    def getLogFilePath(self, program='dials.symmetry'):
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        return logPath

    def getDatasets(self):
        return dutils.getDatasets(self.getInputModelFile())

    def getLogOutput(self):
        logOutput = dutils.readLog(
            self.getLogFilePath(),
            'Recommended',
            'Saving')
        return logOutput.strip()

    def _checkWriteModel(self):
        return self.getSetModel() != self.getInputModelFile()

    def _checkWriteRefl(self):
        return self.getSetRefl() != self.getInputReflFile()

    def _prepCommandline(self, program):
        "Create the command line input to run dials programs"

        # Input basic parameters
        logPath = self.getLogFilePath(program)
        params = "{} {} output.log={} output.experiments={} output.reflections={}".format(
            self.getInputModelFile(),
            self.getInputReflFile(),
            logPath,
            self.getOutputModelFile(),
            self.getOutputReflFile(),
        )

        # Update the command line with additional parameters

        if self.commandLineInput.get():
            params += " {}".format(self.commandLineInput.get())

        return params

    def _prepCommandlineReport(self):
        "Create the command line input to run dials programs"
        # Input basic parameters
        params = "{} {} output.html={} output.external_dependencies={}".format(
            self.getOutputModelFile(),
            self.getOutputReflFile(),
            self.getOutputHtmlFile(),
            self.extDepOptions[self.externalDependencies.get()]
        )

        if self.pixelsPerBin.get():
            params += " pixels_per_bin={}".format(self.pixelsPerBin.get())

        if self.centroidDiffMax.get():
            params += " centroid_diff_max={}".format(
                self.centroidDiffMax.get())

        if self.commandLineInputReport.get() not in (None, ''):
            params += " {}".format(self.commandLineInputReport.get())

        return params
