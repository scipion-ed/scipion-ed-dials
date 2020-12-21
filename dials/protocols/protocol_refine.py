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
from pwed.protocols import EdProtRefineSpots
from pwed.convert import find_subranges
from dials.convert import writeJson, readRefl, writeRefl, writeRefinementPhil, copyDialsFile
from dials.constants import *


class DialsProtRefineSpots(EdProtRefineSpots):
    """ Protocol for refining spots using Dials
    """
    _label = 'refine'

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        # EdProtRefineSpots._defineParams(self, form)

        # The start of the actually relevant part.
        form.addSection(label='Input')

        form.addParam('inputSet', pwprot.PointerParam,
                      pointerClass='SetOfIndexedSpots',
                      label="Input indexed spots",
                      help="")

        # TODO: Add auto option
        form.addParam('scanVarying', pwprot.BooleanParam,
                      label='Perform scan-varying refinement?', default=False,
                      help="Allow models that are not forced to be static to vary during the scan, Auto will run one macrocycle with static then scan varying refinement for the crystal",
                      )

        # Help messages are copied from the DIALS documentation at
        # https://dials.github.io/documentation/programs/dials_index.html
        form.addParam('Nproc', pwprot.IntParam,
                      label="How many processors do you want to use?",
                      default=1, help='The number of processes to use. Not all choices of refinement engine support nproc > 1.'
                      'Where multiprocessing is possible, it is helpful only in certain circumstances,'
                      'so this is not recommended for typical use.',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        group = form.addGroup('Parametrisation')

        group.addParam('beamFixAll', pwprot.BooleanParam,
                       label='Fix all beam parameters?', default=False,
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted to a"
                       "direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow refinement of"
                       "the unit cell volume.",
                       )

        group.addParam('beamFixInSpindlePlane', pwprot.BooleanParam,
                       label='Fix beam in spindle plane?', default=True, condition="beamFixAll==False",
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted to a"
                       "direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow refinement of"
                       "the unit cell volume.",
                       )

        group.addParam('beamFixOutSpindlePlane', pwprot.BooleanParam,
                       label='Fix beam out of spindle plane?', default=False, condition="beamFixAll==False",
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted to"
                       "a direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow refinement"
                       "of the unit cell volume.",
                       )

        group.addParam('beamFixWavelength', pwprot.BooleanParam,
                       label='Fix beam wavelength?', default=True, condition="beamFixAll==False",
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted"
                       "to a direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow"
                       "refinement of the unit cell volume.",
                       )

        group.addParam('beamForceStatic', pwprot.BooleanParam,
                       label='Force static parametrisation for the beam?', default=True, condition="scanVarying==True",
                       help="Force a static parameterisation for the beam when doing scan-varying refinement",
                       )

        group.addParam('crystalFixCell', pwprot.BooleanParam,
                       label='Crystal: Fix cell?', default=True,
                       help="Fix crystal parameters",
                       )

        group.addParam('crystalFixOrientation', pwprot.BooleanParam,
                       label='Crystal: Fix orientation?', default=True,
                       help="Fix crystal parameters",
                       )

        group.addParam('detectorFixAll', pwprot.BooleanParam,
                       label='Fix all detector parameters?', default=False,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.",
                       )

        group.addParam('detectorFixPosition', pwprot.BooleanParam,
                       label='Fix detector position?', default=False,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.", condition="detectorFixAll==False",
                       )

        group.addParam('detectorFixOrientation', pwprot.BooleanParam,
                       label='Fix detector orientation?', default=False,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.", condition="detectorFixAll==False",
                       )

        group.addParam('detectorFixdistance', pwprot.BooleanParam,
                       label='Fix detector distance?', default=False,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.", condition="detectorFixAll==False",
                       )

        group = form.addGroup('Refinery')

        group.addParam('doSetMaxIterations', pwprot.BooleanParam,
                       label='Do you want to set the maximum number of iterations?', default=False,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       help="Maximum number of iterations in refinement before termination.",
                       )

        group.addParam('refineryMaxIterations', pwprot.IntParam,
                       default=None,
                       allowsNull=True, help="Maximum number of iterations in refinement before termination."
                       "None implies the engine supplies its own default.",
                       label='Max iterations', condition="doSetMaxIterations",
                       )

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
        self._insertFunctionStep('refineStep')
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

    def refineStep(self):
        program = 'dials.refine'
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

        reflectionData = readRefl(self.getOutputReflFile())
        iSpot = IndexedSpot()
        numberOfSpots = reflectionData[2]
        reflDict = reflectionData[4]

        outputSet.setSpots(numberOfSpots)

        for i in range(0, numberOfSpots):
            iSpot.setObjId(i+1)
            iSpot.setSpotId(reflDict['id'][i])
            iSpot.setBbox(reflDict['bbox'][i])
            iSpot.setFlag(reflDict['flags'][i])
            iSpot.setIntensitySumValue(reflDict['intensity.sum.value'][i])
            iSpot.setIntensitySumVariance(
                reflDict['intensity.sum.variance'][i])
            iSpot.setNSignal(reflDict['n_signal'][i])
            iSpot.setPanel(reflDict['panel'][i])
            try:
                iSpot.setShoebox(reflDict['shoebox'][i])
            except IndexError:
                pass
            iSpot.setXyzobsPxValue(reflDict['xyzobs.px.value'][i])
            iSpot.setXyzobsPxVariance(reflDict['xyzobs.px.variance'][i])
            outputSet.append(iSpot)

        outputSet.write()

        self._defineOutputs(outputRefinedSpots=outputSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        if self.getDatasets() not in (None, ''):
            summary.append(self.getDatasets())

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputModelFile(self):
        if self.getSetModel():
            return self.getSetModel()
        else:
            return self._getExtraPath('indexed.expt')

    def getInputReflFile(self):
        if self.getSetRefl():
            return self.getSetRefl()
        else:
            return self._getExtraPath('indexed.refl')

    def getOutputModelFile(self):
        return self._getExtraPath('refined.expt')

    def getOutputReflFile(self):
        return self._getExtraPath('refined.refl')

    def getOutputHtmlFile(self):
        return self._getExtraPath('dials.report.html')

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

    def getDatasets(self):
        return dutils.getDatasets(self.getInputModelFile())

    def getLogOutput(self):
        return ''

    def _checkWriteModel(self):
        return self.getSetModel() != self.getInputModelFile()

    def _checkWriteRefl(self):
        return self.getSetRefl() != self.getInputReflFile()

    def _prepCommandline(self, program):
        "Create the command line input to run dials programs"

        # Input basic parameters
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        params = "{} {} output.log={} output.experiments={} output.reflections={}".format(
            self.getInputModelFile(),
            self.getInputReflFile(),
            logPath,
            self.getOutputModelFile(),
            self.getOutputReflFile(),
        )

        # Update the command line with additional parameters

        if self.Nproc.get() not in (None, 1):
            params += " nproc={}".format(self.Nproc.get())

        if self.scanVarying:
            params += " scan_varying=True"

        if self.beamFixAll:
            params += " beam.fix='*all in_spindle_plane out_spindle_plane wavelength'"
        else:
            beamfix = []
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
            params += " beam.fix={}".format("".join(beamfix))

        if self.beamForceStatic.get() not in (None, True):
            params += " beam.force_static={}".format(
                self.beamForceStatic.get())

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
            params += " crystal.fix={}".format(''.join(crystalfix))

        detectorfix = []
        if self.detectorFixAll:
            detectorfix += "'*all position orientation distance'"
        else:
            detectorfix += "'all "
            if self.detectorFixPosition:
                detectorfix += "*"
            detectorfix += "position "
            if self.detectorFixOrientation:
                detectorfix += "*"
            detectorfix += "orientation "
            if self.detectorFixdistance:
                detectorfix += "*"
            detectorfix += "distance'"
        params += " detector.fix={}".format(''.join(detectorfix))

        if self.refineryMaxIterations.get() is not None:
            params += " refinery.max_iterations={}".format(
                self.refineryMaxIterations.get())

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
