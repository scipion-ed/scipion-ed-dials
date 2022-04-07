# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              V. E.G: Bengtsson (viktor.bengtsson@mmk.su.se) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Department of Materials and Environmental Chemistry,
# *     Stockholm University
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

import json
import textwrap
from dials.objects import IsNoneError

import pyworkflow.protocol as pwprot
import dials.utils as dutils

from pwed.objects import IndexedSpot
from pwed.protocols import EdProtIndexSpots
from dials.protocols import DialsProtBase, HtmlBase, RefineParamsBase
from dials.convert import readRefl, copyDialsFile
from dials.constants import *


class DialsProtIndexSpots(EdProtIndexSpots, DialsProtBase):
    """ Protocol for indexing spots using Dials
    """
    _label = 'index'

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        # EdProtIndexSpots._defineParams(self, form)

        # Check which parts of indexing to perform. Reindexing is probably too small to warrant
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
                      label="Do you want to reindex after refining "
                      "Bravais settings?",
                      condition='doRefineBravaisSettings')

        # Keep options to maintain compatibility with workflows etc
        form.addHidden('doReindexModel', pwprot.BooleanParam,
                       default=False,
                       label="Do you want to reindex the experimental model?")

        form.addHidden('doReindexReflections', pwprot.BooleanParam,
                       default=True,
                       label="Do you want to reindex the experimental "
                       "reflections?")

        # The start of typical inputs.

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
                      label="How many processes do you want to use?",
                      default=1,
                      help="The number of processes to use.")

        form.addParam('enterSpaceGroup', pwprot.BooleanParam,
                      default=False,
                      label='Use a known space group?')

        form.addParam('knownSpaceGroup', pwprot.StringParam,
                      condition='enterSpaceGroup',
                      default='',
                      label='Space group:')

        form.addParam('enterUnitCell', pwprot.BooleanParam,
                      default=False,
                      label='Use a known unit cell?')

        form.addParam('knownUnitCell', pwprot.StringParam,
                      condition='enterUnitCell',
                      default='',
                      label='Unit cell:')

        group = form.addGroup('Other indexing parameters',
                              condition='doIndex',
                              expertLevel=pwprot.LEVEL_ADVANCED,)

        group.addParam('indexMmSearchScope', pwprot.FloatParam,
                       default=4.0,
                       help="Global radius of origin offset search.",
                       label='mm search scope',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('indexWideSearchBinning', pwprot.FloatParam,
                       default=2,
                       help="Modify the coarseness of the wide grid search "
                       "for the beam centre.",
                       label='Wide search binning',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('indexMinCellVolume', pwprot.FloatParam,
                       default=25,
                       help="Minimum unit cell volume (in Angstrom^3).",
                       label='Min cell volume',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('indexMinCell', pwprot.FloatParam,
                       default=3,
                       help="Minimum length of candidate unit cell basis "
                       "vectors (in Angstrom).",
                       label='Min_cell',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('indexMaxCell', pwprot.FloatParam,
                       default=None,
                       label='Max_cell',
                       allowsNull=True,
                       help="Maximum length of candidate unit cell basis "
                       "vectors (in Angstrom).",
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('misindexCheckGridScope', pwprot.IntParam,
                       default=0,
                       help="Search scope for testing misindexing "
                       "on h, k, l.",
                       label='Misindexing check grid scope',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('doFilter_ice', pwprot.BooleanParam,
                       default=False,
                       label='Filter ice?',
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       help="Filter out reflections at typical ice ring "
                       "resolutions before max_cell estimation.")

        group = form.addGroup('Refinement parameter configuration',
                              condition='doIndex',
                              expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('refineNproc', pwprot.IntParam,
                       default=1,
                       label='nproc',
                       help="The number of processes to use. Not all choices "
                       "of refinement engine support nproc > 1. Where "
                       "multiprocessing is possible, it is helpful only in "
                       "certain circumstances, so this is not recommended for"
                       " typical use."
                       )

        RefineParamsBase._defineParametrisations(self, form)

        group = form.addGroup('Refinery',
                              expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('doSetMaxIterations', pwprot.BooleanParam,
                       label="Do you want to set the maximum number of "
                       "iterations?",
                       default=False,
                       help="Maximum number of iterations in refinement "
                       "before termination.",
                       )

        group.addParam('refineryMaxIterations', pwprot.IntParam,
                       default=None,
                       allowsNull=True,
                       help="Maximum number of iterations in refinement "
                       "before termination."
                       "None implies the engine supplies its own default.",
                       label='Max iterations',
                       condition="doSetMaxIterations",
                       )

        # Allow some options if the Bravais settings are to be refined
        group = form.addGroup('Refine Bravais settings',
                              condition='doRefineBravaisSettings')

        group.addParam('refineBravNproc', pwprot.IntParam,
                       default=4,
                       label="How many processors do you want to use?",
                       help="The number of processes to use.")

        group.addParam('copyBeamFix', pwprot.BooleanParam,
                       default=True,
                       label="Copy beam parametrisation from indexing?",
                       help="Do you want to use the indexing parametrisation "
                       "of the beam instead of the default parametrisation "
                       "for Bravais setting refinement?",)

        group.addParam('copyCrystalFix', pwprot.BooleanParam,
                       default=True,
                       label="Copy crystal parametrisation from indexing?",
                       help="Do you want to use the indexing parametrisation "
                       "of the crystal instead of the default parametrisation"
                       " for Bravais setting refinement?",)

        group.addParam('copyDetectorFix', pwprot.BooleanParam,
                       default=True,
                       label="Copy detector parametrisation from indexing?",
                       help="Do you want to use the indexing parametrisation "
                       "of the detector instead of the default parametrisation"
                       " for Bravais setting refinement?",)

        group.addParam('copyGonioFix', pwprot.BooleanParam,
                       default=True,
                       label="Copy goniometer parametrisation from indexing?",
                       help="Do you want to use the indexing parametrisation "
                       "of the goniometer instead of the default "
                       "parametrisation for Bravais setting refinement?",)

        # Allow an easy way to import a phil file with parameters
        # Not reusing from base protocols to allow different ones for
        # each function
        group = form.addGroup('Add parameters with phil files',
                              expertLevel=pwprot.LEVEL_ADVANCED,)
        group.addParam('extraPhilPathIndexing', pwprot.PathParam,
                       allowsNull=True,
                       default=None,
                       condition="doIndex==True",
                       label="Additional indexing phil file",
                       help="Enter the path to a phil file that you want to "
                       "add to include.")
        group.addParam('extraPhilPathBravais', pwprot.PathParam,
                       allowsNull=True,
                       default=None,
                       condition="doRefineBravaisSettings",
                       label="Additional bravais settings phil file",
                       help="Enter the path to a phil file that you want to "
                       "add to include.")
        group.addParam('extraPhilPathReindexing', pwprot.PathParam,
                       allowsNull=True,
                       default=None,
                       condition="doReindex",
                       label="Additional reindexing phil file",
                       help="Enter the path to a phil file that you want to "
                       "add to include.")

        # Allow adding anything else with command line syntax
        group = form.addGroup('Raw command line input parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED)
        group.addParam('commandLineInputIndexing', pwprot.StringParam,
                       default='',
                       condition="doIndex==True",
                       label='Indexing command line',
                       help="Anything added here will be added at the end of"
                       " the command line for indexing")
        group.addParam('commandLineInputBravais', pwprot.StringParam,
                       default='',
                       condition="doRefineBravaisSettings",
                       label='Bravais setting command line',
                       help="Anything added here will be added at the end of"
                       " the command line for Bravais settings refinement")
        group.addParam('commandLineInputReindexing', pwprot.StringParam,
                       default='',
                       condition="doReindex",
                       label='Reindexing command line',
                       help="Anything added here will be added at the end of "
                       "the command line for reindexing")

        # Add a section for creating an html report
        HtmlBase._defineHtmlParams(self, form)

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
        self.info(f"Number of images: {inputImages.getSize()}")
        self.info(f"Number of spots: {inputSpots.getSpots()}")
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
        if self.getBravaisSummary() is None:
            raise IsNoneError

    def reindexStep(self):
        program = 'dials.reindex'
        arguments = self._prepReindexCommandline()
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    def createOutputStep(self):
        # Find the most processed model file and reflection file and copy
        # to output
        if dutils.existsPath(self.getReindexedModelFile()):
            copyDialsFile(self.getReindexedModelFile(),
                          self.getOutputModelFile())
        elif dutils.existsPath(self.getBravaisModelFile(self.getBravaisId())):
            copyDialsFile(self.getBravaisModelFile(self.getBravaisId()),
                          self.getOutputModelFile())
        elif dutils.existsPath(self.getIndexedModelFile()):
            copyDialsFile(self.getIndexedModelFile(),
                          self.getOutputModelFile())

        if dutils.existsPath(self.getReindexedReflFile()):
            copyDialsFile(self.getReindexedReflFile(),
                          self.getOutputReflFile())
        elif dutils.existsPath(self.getIndexedReflFile()):
            copyDialsFile(self.getIndexedReflFile(),
                          self.getOutputReflFile())

        # Check that the indexing created proper output
        dutils.verifyPathExistence(self.getOutputReflFile(),
                                   self.getOutputModelFile())

        # TODO: Add Diffraction images as well
        outputSet = self._createSetOfIndexedSpots()
        outputSet.setDialsModel(self.getOutputModelFile())
        outputSet.setDialsRefl(self.getOutputReflFile())
        try:
            # FIXME: readRefl does not work when reading
            # reindexed.refl. Complains about "Extra data"
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

    # -------------------------- BASE methods to be overridden -----------------

    INPUT_EXPT_FILENAME = 'imported.expt'
    OUTPUT_EXPT_FILENAME = 'indexed.expt'
    INPUT_REFL_FILENAME = 'strong.refl'
    OUTPUT_REFL_FILENAME = 'indexed.refl'

    # -------------------------- UTILS functions ------------------------------

    def getIndexedModelFile(self):
        return self._getTmpPath('indexed.expt')

    def getIndexedReflFile(self):
        return self._getTmpPath('indexed.refl')

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
            indexOut = f"\n{textwrap.dedent(indexOutput)}"
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
            bravaisOut = f"\n{textwrap.dedent(bravaisOutput)}"
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
        fn = f"bravais_setting_{fileId}.expt"
        return self.getBravaisPath(fn)

    def getChangeOfBasisOp(self, fileId):
        # This function is likely the problem in issue #10
        summary = self.getBravaisSummary()
        if summary is None:
            # TODO: Add parameter to manually supply default
            change_of_basis_op = 'a,b,c'
        elif fileId is None:
            cbop = summary[self.getBravaisId()]["cb_op"]
            change_of_basis_op = cbop
        else:
            # TODO: #20 Allow manual selection of Bravais Setting.
            cbop = summary[fileId]["cb_op"]
            change_of_basis_op = cbop
        return change_of_basis_op

    def getBravaisSummary(self):
        fn = self.getBravaisPath('bravais_summary.json')
        if dutils.existsPath(fn):
            with open(fn) as f:
                summary = json.load(f)
            return summary
        else:
            return None

    def getReindexedModelFile(self):
        return self._getTmpPath('reindexed.expt')

    def getReindexedReflFile(self):
        return self._getTmpPath('reindexed.refl')

    def getKnownUnitCell(self):
        return self.fixString(self.knownUnitCell.get())

    # Placeholder for using phils as default

    def getPhilPath(self):
        return self._getTmpPath('index.phil')

    # Get the path of additional phil files
    def getExtraPhilsPathIndexing(self):
        return self.extraPhilPathIndexing.get('').strip()

    def getExtraPhilsPathBravais(self):
        return self.extraPhilPathBravais.get('').strip()

    def getExtraPhilsPathReindexing(self):
        return self.extraPhilPathReindexing.get('').strip()

    def _checkWriteModel(self):
        return self.getSetModel() != self.getInputModelFile()

    def _checkWriteRefl(self):
        return self.getSetRefl() != self.getInputReflFile()

    def _prepIndexCommandline(self, program):
        "Create the command line input to run dials programs"

        # Input basic parameters
        logPath = self.getLogFilePath(program)
        params = (f"{self.getInputModelFile()} {self.getInputReflFile()} "
                  f"output.log={logPath} "
                  f"output.experiments={self.getIndexedModelFile()} "
                  f"output.reflections={self.getIndexedReflFile()}")

        # Update the command line with additional parameters

        if self.indexNproc.get() not in (None, 1):
            params += f" indexing.nproc={self.indexNproc.get()}"

        if self.enterSpaceGroup.get():
            params += (f" indexing.known_symmetry.space_group="
                       f"{self.knownSpaceGroup.get()}")

        if self.enterUnitCell.get():
            params += (f" indexing.known_symmetry.unit_cell="
                       f"{self.getKnownUnitCell()}")

        if self.indexMmSearchScope.get() not in (None, 4.0):
            params += (f" indexing.mm_search_scope="
                       f"{self.indexMmSearchScope.get()}")

        if self.indexWideSearchBinning.get() not in (None, 2):
            params += (f" indexing.wide_search_binning="
                       f"{self.indexWideSearchBinning.get()}")

        if self.indexMinCellVolume.get() not in (None, 25):
            params += (f" indexing.min_cell_volume="
                       f"{self.indexMinCellVolume.get()}")

        if self.indexMinCell.get() not in (None, 3.0):
            params += f" indexing.min_cell={self.indexMinCell.get()}"

        if self.indexMaxCell.get() is not None:
            params += f" indexing.max_cell={self.indexMaxCell.get()}"

        if self.misindexCheckGridScope.get() not in (None, 0):
            params += (f" check_misindexing.grid_search_scope="
                       f"{self.misindexCheckGridScope.get()}")

        if self.doFilter_ice.get():
            params += (f" indexing.max_cell_estimation.filter_ice="
                       f"{self.doFilter_ice.get()}")

        if self.refineNproc.get() not in (None, 1):
            params += f" refinement.nproc={self.refineNproc.get()}"

        params += RefineParamsBase.getBeamFixParams(self)

        params += RefineParamsBase.getCrystalFixParams(self)

        params += RefineParamsBase.getDetectorFixParams(self)

        params += RefineParamsBase.getGonioFixParams(self)

        if self.refineryMaxIterations.get() is not None:
            params += (f" refinery.max_iterations="
                       f"{self.refineryMaxIterations.get()}")

        if self.extraPhilPathIndexing.get():
            params += f" {self.getExtraPhilsPathIndexing()}"

        if self.commandLineInputIndexing.get():
            params += f" {self.commandLineInputIndexing.get()}"

        return params

    def _prepBravaisCommandline(self, program):
        "Create the command line input to run dials programs"
        # Input basic parameters
        logPath = self.getLogFilePath(program)
        params = (f"{self.getIndexedModelFile()} {self.getIndexedReflFile()} "
                  f"output.log={logPath} "
                  f"output.directory={self.getBravaisPath()}")

        # Update the command line with additional parameters

        if self.refineBravNproc.get() not in (None, 4):
            params += f" nproc={self.refineBravNproc.get()}"

        if self.copyBeamFix:
            params += RefineParamsBase.getBeamFixParams(self)

        if self.copyCrystalFix:
            params += RefineParamsBase.getCrystalFixParams(self)

        if self.copyDetectorFix:
            params += RefineParamsBase.getDetectorFixParams(self)

        if self.copyGonioFix:
            params += RefineParamsBase.getGonioFixParams(self)

        if self.extraPhilPathBravais.get():
            params += f" {self.getExtraPhilsPathBravais()}"

        if self.commandLineInputBravais.get():
            params += f" {self.commandLineInputBravais.get()}"

        return params

    def _prepReindexCommandline(self):
        "Create the command line input to run dials programs"
        # Input basic parameters
        # FIXME: Fix issue #10
        params = (f"change_of_basis_op="
                  f"{self.getChangeOfBasisOp(self.getBravaisId())}")

        if self.doReindexModel.get():
            params += (f" {self.getIndexedModelFile()} "
                       f"output.experiments={self.getReindexedModelFile()}")

        if self.doReindexReflections.get():
            params += (f" {self.getIndexedReflFile()} "
                       f"output.reflections={self.getReindexedReflFile()}")

        if self.extraPhilPathReindexing.get():
            params += f" {self.getExtraPhilsPathReindexing()}"

        if self.commandLineInputReindexing.get():
            params += f" {self.commandLineInputReindexing.get()}"

        return params
