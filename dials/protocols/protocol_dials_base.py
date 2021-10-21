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


import pyworkflow.protocol as pwprot
import dials.utils as dutils

from pwed.protocols import EdBaseProtocol
from dials.constants import *


class DialsProtBase(EdBaseProtocol):
    """ Base protocol for DIALS
    """
    # Define default filenames
    INPUT_EXPT_FILENAME = 'input.expt'
    OUTPUT_EXPT_FILENAME = 'output.expt'
    INPUT_REFL_FILENAME = 'input.refl'
    OUTPUT_REFL_FILENAME = 'output.refl'
    OUTPUT_HTML_FILENAME = 'dials.report.html'
    OUTPUT_JSON_FILENAME = "dials.program.json"

    # -------------------------- STEPS functions -----------------------------

    def makeHtmlReportStep(self):
        HtmlBase.makeHtmlReportStep(self)

    # -------------------------- UTILS functions -----------------------------

    def getInputModelFile(self, inputSource=None):
        if self.getSetModel(inputSource):
            return self.getSetModel(inputSource)
        else:
            return self._getExtraPath(self.INPUT_EXPT_FILENAME)

    def getInputReflFile(self, inputSource=None):
        if self.getSetRefl(inputSource):
            return self.getSetRefl(inputSource)
        else:
            return self._getExtraPath(self.INPUT_REFL_FILENAME)

    def getDatasets(self):
        return dutils.getDatasets(self.getInputModelFile())

    def getOutputModelFile(self):
        return self._getExtraPath(self.OUTPUT_EXPT_FILENAME)

    def getOutputReflFile(self):
        return self._getExtraPath(self.OUTPUT_REFL_FILENAME)

    def getOutputJsonFile(self):
        return self._getExtraPath(self.OUTPUT_JSON_FILENAME)

    def getProjectName(self, newProjectName=None):
        # Function to get the name of the overall project.
        # Useful for the metadata of exported files. Could also be used
        # in summaries.
        if newProjectName:
            return newProjectName
        else:
            return self.getProject().getShortName()

    def getCrystalName(self, newCrystalName=None):
        # Function to return crystal name. Can be connected to pwed.objects
        # properties at a later point to define a name during import and use
        # it throughout the processing. Until then it provides "XTAL" as
        # default and can be overwritten in protocols.
        if newCrystalName:
            return newCrystalName
        else:
            return "XTAL"

    def _getModelSources(self, inputSource=None):
        sources = []
        if inputSource is not None:
            sources.append(inputSource)
        try:
            sources.append(self.inputImages.get())
        except AttributeError:
            pass
        try:
            sources.append(self.inputSpots.get())
        except AttributeError:
            pass
        try:
            sources.append(self.inputSet.get())
        except AttributeError:
            pass

        return sources

    def getSetModel(self, inputSource=None):
        for source in self._getModelSources(inputSource):
            try:
                if dutils.existsPath(source.getDialsModel()):
                    return source.getDialsModel()
            except TypeError:
                pass

    def getSetRefl(self, inputSource=None):
        for source in self._getModelSources(inputSource):
            try:
                if dutils.existsPath(source.getDialsRefl()):
                    return source.getDialsRefl()
            except TypeError:
                pass

    def getLogFilePath(self, program='dials.*'):
        logPath = f"{self._getLogsPath()}/{program}.log"
        return logPath

    def getLogOutput(self):
        return ''

    def _checkWriteModel(self, inputSource=None):
        return self.getSetModel(inputSource) != self.getInputModelFile(inputSource)

    def _checkWriteRefl(self, inputSource=None):
        return self.getSetRefl(inputSource) != self.getInputReflFile(inputSource)

    def _initialParams(self, program):
        # Base method that can more easily be overridden when needed
        params = (f"{self.getInputModelFile()} {self.getInputReflFile()} "
                  f"output.log={self.getLogFilePath(program)} "
                  f"output.experiments={self.getOutputModelFile()} "
                  f"output.reflections={self.getOutputReflFile()}")

        return params

    def _extraParams(self):
        params = ""
        return params

    def _getExtraPhilsPath(self):
        return PhilBase._addPhilPath(self)

    def _getCLI(self):
        return CliBase._getCommandLineInput(self).rstrip()

    def _prepareCommandline(self, program=None):
        """Create the command line input to run dials programs"""

        # Input basic parameters
        self.info(f"Program is {program}")
        params = self._initialParams(program)

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
                       help="Anything added here will be added at the "
                       "end of the command line")

    def _getCommandLineInput(self):
        if self.commandLineInput.get():
            return f" {self.commandLineInput.get()}"
        else:
            return ""


class RefineParamsBase(EdBaseProtocol):

    # -------------------------- DEFINE param functions ----------------------

    def _defineParametrisations(self, form):
        group = form.addGroup('Model parametrisation')

        group.addHidden('beamFixAll', pwprot.BooleanParam,
                        label='Fix all beam parameters?', default=False,
                        help="Whether to fix beam parameters. By default, "
                        "in_spindle_plane is selected, and one of the two "
                        "parameters is fixed. If a goniometer is present "
                        "this leads to the beam orientation being restricted"
                        " to a direction in the initial spindle-beam plane. "
                        "Wavelength is also fixed by default, to allow "
                        "refinement of the unit cell volume.",
                        )

        group.addParam('beamFixInSpindlePlane', pwprot.BooleanParam,
                       label='Fix beam in spindle plane?', default=True,
                       condition="beamFixAll==False",
                       help="Whether to fix beam parameters. By default, "
                       "in_spindle_plane is selected, and one of the two "
                       "parameters is fixed. If a goniometer is present this "
                       "leads to the beam orientation being restricted to a "
                       "direction in the initial spindle-beam plane. "
                       "Wavelength is also fixed by default, to allow "
                       "refinement of the unit cell volume.",
                       )

        group.addParam('beamFixOutSpindlePlane', pwprot.BooleanParam,
                       label='Fix beam out of spindle plane?', default=False,
                       condition="beamFixAll==False",
                       help="Whether to fix beam parameters. By default, "
                       "in_spindle_plane is selected, and one of the two "
                       "parameters is fixed. If a goniometer is present "
                       "this leads to the beam orientation being restricted "
                       "to a direction in the initial spindle-beam plane. "
                       "Wavelength is also fixed by default, to allow "
                       "refinement of the unit cell volume.",
                       )

        group.addParam('beamFixWavelength', pwprot.BooleanParam,
                       label='Fix beam wavelength?', default=True,
                       condition="beamFixAll==False",
                       help="Whether to fix beam parameters. By default, "
                       "in_spindle_plane is selected, and one of the two "
                       "parameters is fixed. If a goniometer is present this "
                       "leads to the beam orientation being restricted to a "
                       "direction in the initial spindle-beam plane. "
                       "Wavelength is also fixed by default, to allow "
                       "refinement of the unit cell volume.",
                       )

        group.addParam('beamForceStatic', pwprot.BooleanParam,
                       label="Force static parametrisation for the beam? "
                       "(only applies to scan-varying refinement)",
                       default=True,
                       help="Force a static parametrisation for the beam "
                       "when doing scan-varying refinement",
                       )

        group.addParam('crystalFixCell', pwprot.BooleanParam,
                       label='Crystal: Fix cell?', default=False,
                       help="Fix crystal parameters",
                       )

        group.addParam('crystalFixOrientation', pwprot.BooleanParam,
                       label='Crystal: Fix orientation?', default=False,
                       help="Fix crystal parameters",
                       )

        group.addHidden('detectorFixAll', pwprot.BooleanParam,
                        label='Fix all detector parameters?', default=False,
                        help="Fix detector parameters. The translational "
                        "parameters (position) may be set"
                        "separately to the orientation.",
                        )

        group.addParam('detectorFixPosition', pwprot.BooleanParam,
                       label='Fix detector position?', default=False,
                       help="Fix detector parameters. The translational "
                       "parameters (position) may be set"
                       "separately to the orientation.",
                       condition="detectorFixAll==False",
                       )

        group.addParam('detectorFixOrientation', pwprot.BooleanParam,
                       label='Fix detector orientation?', default=False,
                       help="Fix detector parameters. The translational "
                       "parameters (position) may be set"
                       "separately to the orientation.",
                       condition="detectorFixAll==False",
                       )

        group.addParam('detectorFixDistance', pwprot.BooleanParam,
                       label='Fix detector distance?', default=True,
                       help="Fix detector parameters. The translational "
                       "parameters (position) may be set"
                       "separately to the orientation.",
                       condition="detectorFixAll==False",
                       )

        group.addParam('goniometerFixInBeamPlane', pwprot.BooleanParam,
                       label='Fix goniometer in beam plane?', default=True,
                       help="Whether to fix goniometer parameters. By default,"
                       " fix all. Alternatively the setting matrix can be "
                       "constrained to allow rotation only within the spindle-"
                       "beam plane or to allow rotation only around an axis "
                       "that lies in that plane. Set to None to refine the in "
                       "two orthogonal directions.",
                       )

        group.addParam('goniometerFixOutBeamPlane', pwprot.BooleanParam,
                       label='Fix goniometer out of beam plane?',
                       default=True,
                       help="Whether to fix goniometer parameters. By default,"
                       " fix all. Alternatively the setting matrix can be "
                       "constrained to allow rotation only within the spindle-"
                       "beam plane or to allow rotation only around an axis "
                       "that lies in that plane. Set to None to refine the in "
                       "two orthogonal directions.",
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
        beamfixparams = f" refinement.parameterisation.beam.fix={''.join(beamfix)}"
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
        crystalfixparams = f" refinement.parameterisation.crystal.fix={''.join(crystalfix)}"
        return crystalfixparams

    def getDetectorFixParams(self):
        detectorfix = []
        if self.detectorFixAll or (self.detectorFixPosition and self.detectorFixOrientation and self.detectorFixDistance):
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
        detectorfixparams = f" refinement.parameterisation.detector.fix={''.join(detectorfix)}"
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
        goniofixparams = f" refinement.parameterisation.goniometer.fix={''.join(goniofix)}"
        return goniofixparams


class HtmlBase(EdBaseProtocol):
    def _defineHtmlParams(self, form):
        # Add a section for creating an html report
        form.addSection('HTML report')
        form.addParam('makeReport', pwprot.BooleanParam,
                      label='Do you want to create an HTML report for the output?',
                      default=False,
                      help="",
                      )

        form.addParam('showReport', pwprot.BooleanParam,
                      label='Do you want to open the report as soon as the protocol is done?',
                      default=False,
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
                       help="Whether to use remote external dependencies "
                       "(files relocatable but requires an internet "
                       "connection), local (does not require internet "
                       "connection but files may not be relocatable) or "
                       "embed all external dependencies (inflates the html"
                       " file size).",
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
                       help="Magnitude in pixels of shifts mapped to the "
                       "extreme colours in the heatmap plots centroid_diff_x"
                       " and centroid_diff_y",
                       )

        # Allow adding anything else with command line syntax
        group = form.addGroup('HTML report command line parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED,
                              condition="makeReport",)
        group.addParam('commandLineInputReport', pwprot.StringParam,
                       default='',
                       help="Anything added here will be added at the end"
                       " of the command line")

    # -------------------------- STEPS functions -------------------------------

    def makeHtmlReportStep(self):
        prog = 'dials.report'
        arguments = HtmlBase._prepCommandlineReport(self)
        self.runJob(prog, arguments)
        if self.showReport:
            dutils._showHtmlReport(HtmlBase.getOutputHtmlFile(self))

    # -------------------------- UTILS functions ------------------------------
    def getOutputHtmlFile(self):
        if self.OUTPUT_HTML_FILENAME:
            return self._getExtraPath(self.OUTPUT_HTML_FILENAME)
        else:
            return self._getExtraPath('dials.report.html')

    def _prepCommandlineReport(self):
        "Create the command line input to run dials programs"
        # Input basic parameters
        params = (f"{DialsProtBase.getOutputModelFile(self)} "
                  f"{DialsProtBase.getOutputReflFile(self)} "
                  f"output.html={HtmlBase.getOutputHtmlFile(self)} "
                  f"output.external_dependencies="
                  f"{self.extDepOptions[self.externalDependencies.get()]}")

        if self.pixelsPerBin.get():
            params += f" pixels_per_bin={self.pixelsPerBin.get()}"

        if self.centroidDiffMax.get():
            params += f" centroid_diff_max={self.centroidDiffMax.get()}"

        if self.commandLineInputReport.get() not in (None, ''):
            params += f" {self.commandLineInputReport.get()}"

        return params


class PhilBase(EdBaseProtocol):
    def _definePhilParams(self, form):
        # Allow an easy way to import a phil file with parameters
        form.addParam('extraPhilPath', pwprot.PathParam,
                      expertLevel=pwprot.LEVEL_ADVANCED,
                      allowsNull=True,
                      default=None,
                      label="Add phil file",
                      help="Enter the path to a phil file that you want "
                      "to add to include.")

    def _addPhilPath(self):
        if self.extraPhilPath.get():
            return f" {self.extraPhilPath.get('').strip()}"
        else:
            return ""
