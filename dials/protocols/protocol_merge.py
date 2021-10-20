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

import textwrap
from pwed.protocols.protocol_base import EdProtMerge

import pyworkflow.protocol as pwprot
import pyworkflow.utils as pwutils
import dials.utils as dutils
import dials.convert as dconv

from pwed.protocols import EdProtMerge
from pwed.objects import ExportFile
from pwed.utils import CutRes
from dials.protocols import DialsProtBase, PhilBase, CliBase
from dials.constants import *


class DialsProtMerge(EdProtMerge, DialsProtBase, CutRes):
    """ Protocol for merging spots using Dials
    """
    _label = 'merge'

    # -------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        # EdProtIndexSpots._defineParams(self, form)

        # The start of the actually relevant part.
        form.addSection(label='Input')

        form.addParam('inputSets', pwprot.MultiPointerParam,
                      pointerClass='SetOfIndexedSpots',
                      label='Sets to merge',
                      minNumObjects=1,
                      maxNumObjects=0,
                      )

        self._defineResolutionParams(form)

        group = form.addGroup("MTZ parameters")
        group.addParam("mtzName", pwprot.StringParam,
                       label="Name of output file",
                       default="merged.mtz",
                       )
        group.addParam("xtalName", pwprot.StringParam,
                       label="Crystal name",
                       default="XTAL",
                       )
        group.addParam("datasetNames", pwprot.StringParam,
                       label="Dataset name(s)",
                       allowsNull=True,
                       default=None,
                       help="Dataset name to be used in MTZ file output "
                       "(multiple names allowed for MAD datasets)",
                       )

        form.addParam("assessSpaceGroup", pwprot.BooleanParam,
                      label="Asses space group?",
                      default=True,
                      help="Option to assess space group by testing presence "
                      "of axial reflections",
                      )

        form.addParam("anomalous", pwprot.BooleanParam,
                      label="Output anomalous intensities?",
                      default=True,
                      help="Output anomalous as well as mean intensities.",
                      )

        form.addParam("truncate", pwprot.BooleanParam,
                      label="Truncate merged data?",
                      default=True,)

        form.addParam("wavelengthTolerance", pwprot.FloatParam,
                      label="Wavelength tolerace",
                      default=1e-4,
                      allowsNull=True,
                      help="Absolute tolerance for determining wavelength "
                            "grouping for merging."
                      )

        form.addParam("combinePartials", pwprot.BooleanParam,
                      label="Combine partials?",
                      default=True,
                      help="Combine partials that have the same partial id "
                            "into one reflection, with an updated partiality "
                            "given by the sum of the individual partialities."
                      )

        form.addParam("partialityThreshold", pwprot.FloatParam,
                      label="Partiality threshold",
                      default=0.4,
                      allowsNull=True,
                      help="All reflections with partiality values above the"
                            " partiality threshold will be retained. This is "
                            "done after any combination of partials if "
                            "applicable.",
                      )

        form.addParam("bestUnitCell", pwprot.StringParam,
                      label="Best unit cell",
                      allowsNull=True,
                      default=None,
                      help="Best unit cell value, to use when performing "
                            "resolution cutting, and as the overall unit cell"
                            " in the merged mtz. If undefined, the median cell"
                            " will be used.",
                      )

        form.addParam("nResidues", pwprot.IntParam,
                      allowsNull=True,
                      default=200,
                      help="Number of residues to use in Wilson scaling",
                      )

        group = form.addGroup("Merging controls")
        group.addParam("useInternalVariance", pwprot.BooleanParam,
                       label="Use internal variance?",
                       default=False,
                       )

        group.addParam("nBins", pwprot.IntParam,
                       label="Number of bins (minimum 5)",
                       default=20,
                       allowsNull=True,
                       )

        group.addParam("mergingAnomalous", pwprot.BooleanParam,
                       label="Make reported merging stats anomalous?",
                       default=False,
                       help="Option to control whether reported merging stats"
                             " are anomalous.",
                       )

        # Allow an easy way to import a phil file with parameters
        PhilBase._definePhilParams(self, form)

        # Allow adding anything else with command line syntax
        CliBase._defineCliParams(self, form)

    # -------------------------- INSERT functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(
            'convertInputStep',
            [inputSet.get().getObjId() for inputSet in self.inputSets])
        self._insertFunctionStep('mergeStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, inputSpotId):
        for iS in self.inputSets:
            inputSet = iS.get()
            self.info(f"ObjId: {inputSet.getObjId()}")
            self.info(f"Number of images: {inputSet.getSize()}")
            self.info(f"Number of spots: {inputSet.getSpots()}")
            # Write new model and/or reflection file if no was supplied from set
            if self._checkWriteModel(inputSet):
                dconv.writeJson(inputSet, self.getInputModelFile(inputSet))
            if self._checkWriteRefl(inputSet):
                dconv.writeRefl(inputSet, self.getInputReflFile(inputSet))

    def mergeStep(self):
        program = 'dials.merge'
        arguments = self._prepareCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    def createOutputStep(self):
        outputSet = self._createSetOfExportFiles()
        eFile = ExportFile()
        eFile.setFilePath(self.getExport())
        eFile.setFileType(self.getFileType())
        outputSet.append(eFile)
        outputSet.write()
        self._defineOutputs(exportedFileSet=outputSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        errors.append(self.resolutionValidation())
        return errors

    def _summary(self):
        summary = []

        if self.getDatasets() not in (None, ''):
            summary.append(self.getDatasets())

        nSets = len(self.inputSets)
        if nSets > 1:
            summary.append(f'\nMerged {nSets} different datasets together')
        elif nSets == 1:
            summary.append('\nMerged a single dataset')

        if self.getDMin() is not None:
            summary.append(f'High resolution cutoff at {self.getDMin()} Å')

        if self.getDMax() is not None:
            summary.append(f'Low resolution cutoff at {self.getDMax()} Å')

        if self._getCLI() != '':
            summary.append(f"Additional command line input:\n"
                           f"{self._getCLI().strip()}")

        return summary
    # -------------------------- BASE methods to be overridden -----------------

    INPUT_EXPT_FILENAME = 'scaled.expt'
    INPUT_REFL_FILENAME = 'scaled.refl'
    OUTPUT_HTML_FILENAME = 'dials.merge.html'

    def getOutputHtmlFile(self):
        return self._getExtraPath(self.OUTPUT_HTML_FILENAME)

    def getLogOutput(self):
        logOutput = ''

        if self.getMergingStatisticsLogOutput() not in (None, ''):
            logOutput += self.getMergingStatisticsLogOutput()
        return logOutput

    def _initialParams(self, program):
        # Base method that can more easily be overridden when needed
        params = (f"{self.getAllInputFiles()} "
                  f"output.log={self.getLogFilePath(program)} "
                  f"output.html={self.getOutputHtmlFile()} "
                  f"output.mtz={self.getMtzName()} "
                  f"output.crystal_names={self.crystalName()} "
                  f"output.project_name={self.getProjectName()}"
                  f"assess_space_group={self.assesSpaceGroup.get()} "
                  f"anomalous={self.anomalous.get()} "
                  f"truncate={self.truncate.get()} "
                  f"wavelength_tolerance={self.wavelengthTolerance.get()} "
                  f"combine_partials={self.combinePartials.get()} "
                  f"partiality_threshold={self.partialityThreshold.get()} "
                  )

        return params

    def _extraParams(self):
        params = f"{self.getDatasetNameLine()}"
        if self.getDMin():
            params += f" d_min={self.getDMin()}"

        if self.getDMax():
            params += f" d_max={self.getDMax()}"

        if self.bestUnitCell.get():
            params += f" best_unit_cell={self.bestUnitCell.get()}"

        if self.nResidues.get():
            params += f" n_residues={self.nResidues.get()}"

        params += (f" merging.use_internal_variance="
                   f"{self.useInternalVariance.get()}")

        if self.nBins.get() and self.nBins.get() >= 5:
            params += f" merging.n_bins={self.nBins.get()}"

        params += f" merging.anomalous={self.mergingAnomalous.get()}"

        return params

    # -------------------------- UTILS functions ------------------------------
    def crystalName(self):
        return self.getCrystalName(self.xtalName.get())

    def getMtzName(self):
        return self.mtzName.get()

    def getDatasetNameLine(self):
        if self.datasetNames.get():
            return f" output.dataset_names={self.datasetNames.get()}"
        else:
            return ""

    def getAllInputFiles(self):
        files = ""
        for iS in self.inputSets:
            files += f"{self.getInputModelFile(iS.get())} {self.getInputReflFile(iS.get())} "
        return files.strip()

    def getMergingStatisticsLogOutput(self):
        mergingStats = dutils.readLog(
            self.getLogFilePath(program='dials.merge'),
            'Merging statistics',
            'Writing html report')
        if mergingStats not in (None, ''):
            mergeStats = f"\n{textwrap.dedent(mergingStats)}"
        else:
            mergeStats = mergingStats
        return mergeStats
