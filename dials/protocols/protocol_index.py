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
import json
import textwrap

import pyworkflow.protocol as pwprot
import dials.utils as dutils

from pwed.objects import DiffractionImage, SetOfDiffractionImages, DiffractionSpot, SetOfSpots, IndexedSpot, SetOfIndexedSpots
from pwed.protocols import EdProtIndexSpots
from pwed.convert import find_subranges
from dials.convert import writeJson, readRefl, writeRefl, writeRefinementPhil, copyDialsFile
from dials.constants import *


class DialsProtIndexSpots(EdProtIndexSpots):
    """ Protocol for indexing spots using Dials
    """
    _label = 'index'

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        # EdProtIndexSpots._defineParams(self, form)

        # Check which parts of indexing to perform. Reindexing is probably to small to warrant
        # its own protocol.
        form.addSection(label='Input')

        form.addParam('doIndex', pwprot.BooleanParam,
                      default=True,
                      label='Do you want to index from start?')

        form.addParam('doRefineBravaisSettings', pwprot.BooleanParam,
                      default=False,
                      label='Do you want to refine the Bravais settings?')

        form.addParam('doReindex', pwprot.BooleanParam,
                      default=False,
                      label='Do you want to reindex after refining Bravais settings?',
                      condition='doRefineBravaisSettings')

        # Allow some options that are only relevant for reindexing
        group = form.addGroup('Reindex',
                              condition='doReindex')

        group.addParam('doReindexModel', pwprot.BooleanParam,
                       default=False, label="Do you want to reindex the experimental model?")

        group.addParam('doReindexReflections', pwprot.BooleanParam,
                       default=False, label="Do you want to reindex the experimental reflections?")

        # The start typical inputs.

        form.addParam('inputImages', pwprot.PointerParam,
                      pointerClass='SetOfDiffractionImages',
                      label="Input diffraction images",
                      help="")

        form.addParam('inputSpots', pwprot.PointerParam,
                      pointerClass='SetOfSpots',
                      label='Input strong spots',
                      help="")

        # Help messages are copied from the DIALS documentation at
        # https://dials.github.io/documentation/programs/dials_index.html
        form.addParam('indexNproc', pwprot.IntParam,
                      label="How many processors do you want to use?",
                      default=1,
                      help="The number of processes to use.")

        form.addParam('enterSpaceGroup', pwprot.BooleanParam,
                      default=False, label='Use a known space group?')

        form.addParam('knownSpaceGroup', pwprot.StringParam,
                      condition='enterSpaceGroup', default='', label='Space group:')

        form.addParam('enterUnitCell', pwprot.BooleanParam,
                      default=False, label='Use a known unit cell?')

        form.addParam('knownUnitCell', pwprot.StringParam,
                      condition='enterUnitCell', default='', label='Unit cell:')

        group = form.addGroup('Other indexing parameters',
                              condition='doIndex', expertLevel=pwprot.LEVEL_ADVANCED,)

        group.addParam('indexMmSearchScope', pwprot.FloatParam, default=4.0,
                       help="Global radius of origin offset search.",
                       label='mm search scope',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('indexWideSearchBinning', pwprot.FloatParam, default=2,
                       help="Modify the coarseness of the wide grid search for the beam centre.",
                       label='Wide search binning',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('indexMinCellVolume', pwprot.FloatParam, default=25,
                       help="Minimum unit cell volume (in Angstrom^3).",
                       label='Min cell volume',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('indexMinCell', pwprot.FloatParam, default=3,
                       help="Minimum length of candidate unit cell basis vectors (in Angstrom).",
                       label='Min_cell',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('indexMaxCell', pwprot.FloatParam, default=None,
                       label='Max_cell', allowsNull=True,
                       help="Maximum length of candidate unit cell basis vectors (in Angstrom).",
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('misindexCheckGridScope', pwprot.IntParam, default=0,
                       help="Search scope for testing misindexing on h, k, l.",
                       label='Misindexing check grid scope',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('doFilter_ice', pwprot.BooleanParam, default=False,
                       label='Filter ice?', expertLevel=pwprot.LEVEL_ADVANCED,
                       help="Filter out reflections at typical ice ring resolutions before max_cell estimation.")

        group = form.addGroup('Refinement parameter configuration',
                              condition='doIndex', expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('refineNproc', pwprot.IntParam,
                       default=1, label='nproc',
                       help='The number of processes to use. Not all choices of refinement engine support nproc > 1.'
                       'Where multiprocessing is possible, it is helpful only in certain circumstances,'
                       'so this is not recommended for typical use.'
                       )

        group = form.addGroup('Model parametrisation',
                              condition='doIndex')

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

        group = form.addGroup('Refinery',
                              expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('doSetMaxIterations', pwprot.BooleanParam,
                       label='Do you want to set the maximum number of iterations?', default=False,
                       help="Maximum number of iterations in refinement before termination.",
                       )

        group.addParam('refineryMaxIterations', pwprot.IntParam,
                       default=None,
                       allowsNull=True, help="Maximum number of iterations in refinement before termination."
                       "None implies the engine supplies its own default.",
                       label='Max iterations', condition="doSetMaxIterations",
                       )

        # Allow some options if the Bravais settings are to be refined
        group = form.addGroup('Refine Bravais settings',
                              condition='doRefineBravaisSettings')

        group.addParam('refineBravNproc', pwprot.IntParam,
                       default=4, label="How many processors do you want to use?",
                       help="The number of processes to use.")

        # Allow adding anything else with command line syntax
        group = form.addGroup('Raw command line input parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED)
        group.addParam('commandLineInputIndexing', pwprot.StringParam,
                       default='', condition="doIndex==True",
                       label='Indexing command line',
                       help="Anything added here will be added at the end of the command line for indexing")
        group.addParam('commandLineInputBravais', pwprot.StringParam,
                       default='', condition="doRefineBravaisSettings",
                       label='Bravais setting command line',
                       help="Anything added here will be added at the end of the command line for Bravais settings refinement")
        group.addParam('commandLineInputReindexing', pwprot.StringParam,
                       default='', condition="doReindex",
                       label='Reindexing command line',
                       help="Anything added here will be added at the end of the command line for reindexing")

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
            'convertInputStep', self.inputImages.getObjId(), self.inputSpots.getObjId())
        if self.doIndex:
            self._insertFunctionStep('indexStep')
        if self.doRefineBravaisSettings:
            self._insertFunctionStep('refineBravaisStep')
        if self.doReindex:
            self._insertFunctionStep('reindexStep')
        self._insertFunctionStep('createOutputStep')
        if self.makeReport:
            self._insertFunctionStep('makeHtmlReportStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, inputImgId, inputSpotId):
        inputImages = self.inputImages.get()
        inputSpots = self.inputSpots.get()
        self.info("Number of images: %s" % inputImages.getSize())
        self.info("Number of spots: %s" % inputSpots.getSpots())
        # Write new model and/or reflection file if no was supplied from set
        if self._checkWriteModel():
            self.writeJson(inputImages, self.getInputModelFile())
        if self._checkWriteRefl():
            self.writeRefl(inputSpots, self.getInputReflFile())

    def indexStep(self):
        program = 'dials.index'
        arguments = self._prepIndexCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    def refineBravaisStep(self):
        program = 'dials.refine_bravais_settings'
        arguments = self._prepBravaisCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())
            return
        assert(self.getBravaisSummary() is not None)

    def reindexStep(self):
        program = 'dials.reindex'
        arguments = self._prepReindexCommandline(program)
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
        # Find the most processed model file and reflection file and copy to output
        if self.existsPath(self.getReindexedModelFile()):
            copyDialsFile(self.getReindexedModelFile(),
                          self.getOutputModelFile())
        elif self.existsPath(self.getBravaisModelFile(self.getBravaisId())):
            copyDialsFile(self.getBravaisModelFile(self.getBravaisId()),
                          self.getOutputModelFile())
        elif self.existsPath(self.getIndexedModelFile()):
            copyDialsFile(self.getIndexedModelFile(),
                          self.getOutputModelFile())

        if self.existsPath(self.getReindexedReflFile()):
            copyDialsFile(self.getReindexedReflFile(),
                          self.getOutputReflFile())
        elif self.existsPath(self.getIndexedReflFile()):
            copyDialsFile(self.getIndexedReflFile(),
                          self.getOutputReflFile())

            # Check that the indexing created proper output
        assert(self.existsPath(self.getOutputReflFile()))
        assert(self.existsPath(self.getOutputModelFile()))

        # TODO: Add Diffraction images as well
        outputSet = self._createSetOfIndexedSpots()
        outputSet.setDialsModel(self.getOutputModelFile())
        outputSet.setDialsRefl(self.getOutputReflFile())
        try:
            # FIXME: readRefl does not work when reading reindexed.refl. Complains about "Extra data"
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
        except Exception as e:
            self.info(e)

        outputSet.write()

        self._defineOutputs(outputIndexedSpots=outputSet)

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
            return self._getExtraPath('imported.expt')

    def getInputReflFile(self):
        if self.getSetRefl():
            return self.getSetRefl()
        else:
            return self._getExtraPath('strong.refl')

    def getIndexedModelFile(self):
        return self._getTmpPath('indexed.expt')

    def getIndexedReflFile(self):
        return self._getTmpPath('indexed.refl')

    def getDatasets(self):
        return dutils.getDatasets(self.getInputModelFile())

    def getLogOutput(self):
        logOutput = ''
        if self.getIndexLogOutput() not in (None, ''):
            logOutput += self.getIndexLogOutput()
        if self.getBravaisLogOutput() not in (None, ''):
            logOutput += self.getBravaisLogOutput()
        return logOutput

    def getIndexLogOutput(self):
        try:
            indexOutput = dutils.readLog(
                self.getLogFilePath('dials.index'),
                'crystal models:',
                'Saving',
                flush='###'
            )
        except FileNotFoundError:
            indexOutput = None
        if indexOutput not in (None, ''):
            indexOut = "\n{}".format(textwrap.dedent(indexOutput))
        else:
            indexOut = indexOutput
        return indexOut

    def getBravaisLogOutput(self):
        # Try-except to avoid problems when there is no log file to read
        try:
            bravaisOutput = dutils.readLog(
                self.getLogFilePath('dials.refine_bravais_settings'),
                'Chiral',
                'Saving')
        except FileNotFoundError:
            bravaisOutput = None
        if bravaisOutput not in (None, ''):
            bravaisOut = "\n{}".format(textwrap.dedent(bravaisOutput))
        else:
            bravaisOut = bravaisOutput
        return bravaisOut

    def getBravaisPath(self, fn=None):
        if fn is None:
            return self._getTmpPath()
        else:
            return self._getTmpPath(fn)

    def getBravaisId(self):
        summary = self.getBravaisSummary()
        if summary is None:
            return None
        else:
            # Use highest ranked suggestion
            keys = list(summary.keys())
            return keys[0]

    def getBravaisModelFile(self, fileId):
        fn = "bravais_setting_{}.expt".format(fileId)
        return self.getBravaisPath(fn)

    def getChangeOfBasisOp(self, fileId):
        summary = self.getBravaisSummary()
        if summary is None:
            # TODO: Add parameter to manually supply default
            change_of_basis_op = 'a,b,c'
        elif fileId is None:
            cbop = summary[self.getBravaisId()]["cb_op"]
            change_of_basis_op = cbop
        else:
            # TODO: Optionally output table from find_bravais_settings
            # and let the user choose which to use
            cbop = summary[fileId]["cb_op"]
            change_of_basis_op = cbop
        return change_of_basis_op

    def getBravaisSummary(self):
        fn = self.getBravaisPath('bravais_summary.json')
        if self.existsPath(fn):
            with open(fn) as f:
                summary = json.load(f)
            return summary
        else:
            return None

    def getReindexedModelFile(self):
        return self._getTmpPath('reindexed.expt')

    def getReindexedReflFile(self):
        return self._getTmpPath('reindexed.refl')

    def getOutputModelFile(self):
        return self._getExtraPath('indexed.expt')

    def getOutputReflFile(self):
        return self._getExtraPath('indexed.refl')

    def getOutputHtmlFile(self):
        return self._getExtraPath('dials.report.html')

    def getPhilPath(self):
        return self._getTmpPath('index.phil')

    def existsPath(self, path):
        return os.path.exists(path)

    def getSetModel(self):
        inputImages = self.inputImages.get()
        inputSpots = self.inputSpots.get()
        if self.existsPath(inputSpots.getDialsModel()):
            return inputSpots.getDialsModel()
        elif self.existsPath(inputImages.getDialsModel()):
            return inputImages.getDialsModel()
        else:
            return None

    def getSetRefl(self):
        inputImages = self.inputImages.get()
        inputSpots = self.inputSpots.get()
        if self.existsPath(inputSpots.getDialsRefl()):
            return inputSpots.getDialsRefl()
        elif self.existsPath(inputImages.getDialsRefl()):
            return inputImages.getDialsRefl()
        else:
            return None

    def getLogFilePath(self, program='dials.index'):
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        return logPath

    def _checkWriteModel(self):
        return self.getSetModel() != self.getInputModelFile()

    def _checkWriteRefl(self):
        return self.getSetRefl() != self.getInputReflFile()

    def _prepIndexCommandline(self, program):
        "Create the command line input to run dials programs"

        # Input basic parameters
        logPath = self.getLogFilePath(program)
        params = "{} {} output.log={} output.experiments={} output.reflections={}".format(
            self.getInputModelFile(),
            self.getInputReflFile(),
            logPath,
            self.getIndexedModelFile(),
            self.getIndexedReflFile(),
        )

        # Update the command line with additional parameters

        if self.indexNproc.get() not in (None, 1):
            params += " indexing.nproc={}".format(self.indexNproc.get())

        if self.enterSpaceGroup.get():
            params += " indexing.known_symmetry.space_group={}".format(
                self.knownSpaceGroup.get())

        if self.enterUnitCell.get():
            params += " indexing.known_symmetry.unit_cell={}".format(
                self.knownUnitCell.get())

        if self.indexMmSearchScope.get() not in (None, 4.0):
            params += " indexing.mm_search_scope={}".format(
                self.indexMmSearchScope.get())

        if self.indexWideSearchBinning.get() not in (None, 2):
            params += " indexing.wide_search_binning={}".format(
                self.indexWideSearchBinning.get())

        if self.indexMinCellVolume.get() not in (None, 25):
            params += " indexing.min_cell_volume={}".format(
                self.indexMinCellVolume.get())

        if self.indexMinCell.get() not in (None, 3.0):
            params += " indexing.min_cell={}".format(self.indexMinCell.get())

        if self.indexMaxCell.get() is not None:
            params += " indexing.max_cell={}".format(self.indexMaxCell.get())

        if self.misindexCheckGridScope.get() not in (None, 0):
            params += " check_misindexing.grid_search_scope={}".format(
                self.misindexCheckGridScope.get())

        if self.doFilter_ice.get():
            params += " indexing.max_cell_estimation.filter_ice={}".format(
                self.doFilter_ice.get())

        if self.refineNproc.get() not in (None, 1):
            params += " refinement.nproc={}".format(self.refineNproc.get())

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
        params += " refinement.parameterisation.beam.fix={}".format(
            "".join(beamfix))

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
            params += " refinement.parameterisation.crystal.fix={}".format(
                ''.join(crystalfix))

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
        params += " refinement.parameterisation.detector.fix={}".format(
            ''.join(detectorfix))

        if self.refineryMaxIterations.get() is not None:
            params += " refinery.max_iterations={}".format(
                self.refineryMaxIterations.get())

        if self.commandLineInputIndexing.get():
            params += " {}".format(self.commandLineInputIndexing.get())

        return params

    def _prepBravaisCommandline(self, program):
        "Create the command line input to run dials programs"
        # Input basic parameters
        logPath = self.getLogFilePath(program)
        params = "{} {} output.log={} output.directory={}".format(
            self.getIndexedModelFile(), self.getIndexedReflFile(), logPath, self.getBravaisPath())

        # Update the command line with additional parameters

        if self.refineBravNproc.get() not in (None, 4):
            params += " nproc={}".format(self.refineBravNproc.get())

        if self.commandLineInputBravais.get():
            params += " {}".format(self.commandLineInputBravais.get())

        return params

    def _prepReindexCommandline(self, program):
        "Create the command line input to run dials programs"
        # Input basic parameters
        params = "change_of_basis_op={}".format(
            self.getChangeOfBasisOp(self.getBravaisId()))

        if self.doReindexModel.get():
            params += " {} output.experiments={}".format(
                self.getIndexedModelFile(), self.getReindexedModelFile())

        if self.doReindexReflections.get():
            params += " {} output.reflections={}".format(
                self.getIndexedReflFile(), self.getReindexedReflFile())

        if self.commandLineInputReindexing.get():
            params += " {}".format(self.commandLineInputReindexing.get())

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
