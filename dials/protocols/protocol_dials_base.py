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

from _typeshed import Self
import os
from glob import glob
from pathlib import Path
import json
import textwrap

import pyworkflow.protocol as pwprot
import dials.utils as dutils

from pwed.protocols import EdBaseProtocol
from pwed.convert import find_subranges
from dials.convert import writeJson, readRefl, writeRefl, writeRefinementPhil, copyDialsFile
from dials.constants import *


class DialsProtBase(EdBaseProtocol):
    """ Base protocol for DIALS
    """

    # -------------------------- UTILS functions ------------------------------

    INPUT_EXPT_FILENAME = 'input.expt'
    OUTPUT_EXPT_FILENAME = 'output.expt'
    INPUT_REFL_FILENAME = 'input.refl'
    OUTPUT_REFL_FILENAME = 'output.refl'

    def getInputModelFile(self):
        if self.getSetModel():
            return self.getSetModel()
        else:
            return self._getExtraPath(self.INPUT_EXPT_FILENAME)

    def getInputReflFile(self):
        if self.getSetRefl():
            return self.getSetRefl()
        else:
            return self._getExtraPath(self.INPUT_REFL_FILENAME)

    def getDatasets(self):
        return dutils.getDatasets(self.getInputModelFile())

    def getOutputModelFile(self):
        return self._getExtraPath(self.OUTPUT_EXPT_FILENAME)

    def getOutputReflFile(self):
        return self._getExtraPath(self.OUTPUT_REFL_FILENAME)

    def _getModelSource(self):
        try:
            return self.inputImages.get()
        except AttributeError:
            pass
        try:
            return self.inputSpots.get()
        except AttributeError:
            pass
        try:
            return self.inputSet.get()
        except AttributeError:
            pass
        return None

    def getSetModel(self):
        if self._getModelSource() is None:
            return None
        elif dutils.existsPath(self._getModelSource().getDialsModel()):
            return self._getModelSource().getDialsModel()
        else:
            return None

    def getSetRefl(self):
        if self._getModelSource() is None:
            return None
        elif dutils.existsPath(self._getModelSource().getDialsRefl()):
            return self._getModelSource().getDialsRefl()
        else:
            return None

    def getLogFilePath(self, program='dials.*'):
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        return logPath

    def getLogOutput(self):
        return ''

    def _checkWriteModel(self):
        return self.getSetModel() != self.getInputModelFile()

    def _checkWriteRefl(self):
        return self.getSetRefl() != self.getInputReflFile()

    def _extraParams(self):
        params = ""
        return params

    def _getExtraPhilsPath(self):
        return PhilBase._addPhilPath()

    def _getCLI(self):
        return CliBase._getCommandLineInput().rstrip()

    def _prepareCommandline(self, program):
        "Create the command line input to run dials programs"

        # Input basic parameters
        logPath = self.getLogFilePath(program)
        params = f"{self.getInputModelFile()} {self.getInputReflFile()} output.log={logPath} output.experiments={self.getOutputModelFile()} output.reflections={self.getOutputReflFile()}"

        # Update the command line with additional parameters

        params += self._extraParams()

        params += self._getExtraPhilsPath()

        params += self._getCLI()

        return params


class CliBase(EdBaseProtocol):
    def _defineCliParams(self, form):
        # Allow adding anything else with command line syntax
        group = form.addGroup('Raw command line input parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED)
        group.addParam('commandLineInput', pwprot.StringParam,
                       default='',
                       help="Anything added here will be added at the end of the command line")

    def _getCommandLineInput(self):
        if self.commandLineInput.get():
            return f" {self.commandLineInput.get()}"
        else:
            return ""


class RefineParamsBase(EdBaseProtocol):

    # -------------------------- DEFINE param functions -----------------------

    def _defineParametrisations(self, form):
        group = form.addGroup('Model parametrisation')

        group.addParam('beamFixInSpindlePlane', pwprot.BooleanParam,
                       label='Fix beam in spindle plane?', default=True,
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted to a"
                       "direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow refinement of"
                       "the unit cell volume.",
                       )

        group.addParam('beamFixOutSpindlePlane', pwprot.BooleanParam,
                       label='Fix beam out of spindle plane?', default=False,
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted to"
                       "a direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow refinement"
                       "of the unit cell volume.",
                       )

        group.addParam('beamFixWavelength', pwprot.BooleanParam,
                       label='Fix beam wavelength?', default=True,
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted"
                       "to a direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow"
                       "refinement of the unit cell volume.",
                       )

        group.addParam('crystalFixCell', pwprot.BooleanParam,
                       label='Crystal: Fix cell?', default=False,
                       help="Fix crystal parameters",
                       )

        group.addParam('crystalFixOrientation', pwprot.BooleanParam,
                       label='Crystal: Fix orientation?', default=False,
                       help="Fix crystal parameters",
                       )

        group.addParam('detectorFixPosition', pwprot.BooleanParam,
                       label='Fix detector position?', default=False,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.",
                       )

        group.addParam('detectorFixOrientation', pwprot.BooleanParam,
                       label='Fix detector orientation?', default=False,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.",
                       )

        group.addParam('detectorFixDistance', pwprot.BooleanParam,
                       label='Fix detector distance?', default=True,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.",
                       )

        group.addParam('goniometerFixInBeamPlane', pwprot.BooleanParam,
                       label='Fix goniometer in beam plane?', default=True,
                       help="Whether to fix goniometer parameters. By default, fix all."
                       "Alternatively the setting matrix can be constrained to allow"
                       "rotation only within the spindle-beam plane or to allow"
                       "rotation only around an axis that lies in that plane. Set to"
                       "None to refine the in two orthogonal directions.",
                       )

        group.addParam('goniometerFixOutBeamPlane', pwprot.BooleanParam,
                       label='Fix goniometer out of beam plane?', default=True,
                       help="Whether to fix goniometer parameters. By default, fix all."
                       "Alternatively the setting matrix can be constrained to allow"
                       "rotation only within the spindle-beam plane or to allow"
                       "rotation only around an axis that lies in that plane. Set to"
                       "None to refine the in two orthogonal directions.",
                       )

    def getBeamFixParams(self):
        beamfix = []
        if self.beamFixInSpindlePlane and self.beamFixOutSpindlePlane and self.beamFixWavelength:
            beamfix += "'*all in_spindle_plane out_spindle_plane wavelength'"
        else:
            beamfix += "'all "
            if self.beamFixInSpindlePlane:
                beamfix += "*"
            beamfix += "in_spindle_plane "
            if self.beamFixOutSpindlePlane:
                beamfix += "*"
            beamfix += "out_spindle_plane "
            if self.beamFixWavelength:
                beamfix += "*"
            beamfix += "wavelength'"
        beamfixparams = " refinement.parameterisation.beam.fix={}".format(
            "".join(beamfix))
        return beamfixparams

    def getCrystalFixParams(self):
        crystalfix = []
        if self.crystalFixCell and self.crystalFixOrientation:
            crystalfix += "'*all cell orientation'"
        else:
            crystalfix += "'all "
            if self.crystalFixCell:
                crystalfix += "*"
            crystalfix += "cell "
            if self.crystalFixOrientation:
                crystalfix += "*"
            crystalfix += "orientation'"
        crystalfixparams = " refinement.parameterisation.crystal.fix={}".format(
            ''.join(crystalfix))
        return crystalfixparams

    def getDetectorFixParams(self):
        detectorfix = []
        if self.detectorFixPosition and self.detectorFixOrientation and self.detectorFixDistance:
            detectorfix += "'*all position orientation distance'"
        else:
            detectorfix += "'all "
            if self.detectorFixPosition:
                detectorfix += "*"
            detectorfix += "position "
            if self.detectorFixOrientation:
                detectorfix += "*"
            detectorfix += "orientation "
            if self.detectorFixDistance:
                detectorfix += "*"
            detectorfix += "distance'"
        detectorfixparams = " refinement.parameterisation.detector.fix={}".format(
            ''.join(detectorfix))
        return detectorfixparams

    def getGonioFixParams(self):
        goniofix = []
        if self.goniometerFixInBeamPlane and self.goniometerFixOutBeamPlane:
            goniofix += "'*all in_beam_plane out_beam_plane'"
        else:
            goniofix += "'all "
            if self.goniometerFixInBeamPlane:
                goniofix += "*"
            goniofix += "in_beam_plane "
            if self.goniometerFixOutBeamPlane:
                goniofix += "*"
            goniofix += "out_beam_plane'"
        goniofixparams = " refinement.parameterisation.goniometer.fix={}".format(
            "".join(goniofix))
        return goniofixparams


class HtmlBase(EdBaseProtocol):
    def _defineHtmlParams(self, form):
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

    # -------------------------- STEPS functions -------------------------------

    def makeHtmlReportStep(self):
        prog = 'dials.report'
        arguments = self._prepCommandlineReport()
        self.runJob(prog, arguments)
        if self.showReport:
            dutils._showHtmlReport(self.getOutputHtmlFile())

        # -------------------------- UTILS functions ------------------------------

    OUTPUT_HTML_FILENAME = 'dials.report.html'

    def getOutputHtmlFile(self):
        return self._getExtraPath(self.OUTPUT_HTML_FILENAME)

    def _prepCommandlineReport(self):
        "Create the command line input to run dials programs"
        # Input basic parameters
        params = "{} {} output.html={} output.external_dependencies={}".format(
            DialsProtBase.getOutputModelFile(),
            DialsProtBase.getOutputReflFile(),
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


class PhilBase(EdBaseProtocol):
    def _definePhilParams(self, form):
        # Allow an easy way to import a phil file with parameters
        form.addParam('extraPhilPath', pwprot.PathParam,
                      expertLevel=pwprot.LEVEL_ADVANCED,
                      allowsNull=True,
                      default=None,
                      label="Add phil file",
                      help="Enter the path to a phil file that you want to add to include.")

    def _addPhilPath(self):
        if self.extraPhilPath.get():
            return f" {self.extraPhilPath.get('').strip()}"
        else:
            return ""
