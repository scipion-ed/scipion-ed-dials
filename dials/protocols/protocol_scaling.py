# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              V. E.G: Bengtsson (viktor.bengtsson@mmk.su.se) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Department of Materials and Environmental Chemistry, Stockholm
# *     University
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
# *  e-mail address "scipion@cnb.csic.es"
# *
# **************************************************************************

import textwrap

import pyworkflow.protocol as pwprot
import pyworkflow.utils as pwutils
import dials.utils as dutils
import dials.convert as dconv

from pwed.protocols import EdProtScaling
from pwed.utils import CutRes
from dials.protocols import DialsProtBase, PhilBase, CliBase
from dials.constants import *


class DialsProtScaling(EdProtScaling, DialsProtBase, CutRes):
    """ Protocol for scaling spots using Dials
    """
    _label = "scale"

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        # EdProtIndexSpots._defineParams(self, form)

        # The start of the actually relevant part.
        form.addSection(label="Input")

        inputsetsLabel = "Spots to scale"
        form.addParam("inputSets", pwprot.MultiPointerParam,
                      pointerClass="SetOfIndexedSpots",
                      label=inputsetsLabel,
                      minNumObjects=1,
                      maxNumObjects=0,
                      help="")

        form.addParam("showReport", pwprot.BooleanParam,
                      label="Do you want to view the HTML report after the"
                      " processing?",
                      default=False,
                      )

        group = form.addGroup("Export")
        group.addParam("exportMergedMtz", pwprot.BooleanParam,
                       default=False,
                       label="Export a merged mtz file?",)

        group.addParam("mergedMtzName", pwprot.StringParam,
                       label="Name of merged mtz:",
                       default="merged.mtz",
                       condition="exportMergedMtz")

        group.addParam("exportUnmergedMtz", pwprot.BooleanParam,
                       default=False,
                       label="Export an unmerged mtz file?",)

        group.addParam("unmergedMtzName", pwprot.StringParam,
                       label="Name of unmerged mtz:",
                       default="unmerged.mtz",
                       condition="exportUnmergedMtz")

        group.addParam("specifyExportPath", pwprot.BooleanParam,
                       label="Do you want to specify a directory for the exported file(s)?",
                       default=False,
                       condition="exportMergedMtz or exportUnmergedMtz")

        group.addParam("exportPath", pwprot.PathParam,
                       label="Target directory for export",
                       default="",
                       condition="specifyExportPath")

        group.addParam("crystalName", pwprot.StringParam,
                       label="Crystal name for metadata",
                       default="XTAL",
                       condition="exportMergedMtz or exportUnmergedMtz")

        group = form.addGroup("Cut data")

        self._defineResolutionParams(form)

        group.addParam("partialityCutoff", pwprot.FloatParam,
                       label="Partiality cutoff",
                       default=0.4,
                       allowsNull=True,
                       help="Value below which reflections are removed "
                       "from the dataset due to low partiality.",
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       )

        group.addParam("minIsigi", pwprot.FloatParam,
                       label="min I/sigma(I)",
                       default=-5,
                       allowsNull=True,
                       help="Value below which reflections are removed "
                       "from the dataset due to low I/sigI in either "
                       "profile or summation intensity estimates",
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       )

        group = form.addGroup("Scaling options")
        group.addParam("checkConsistentIndexing", pwprot.BooleanParam,
                       default=False,
                       label="Check indexing consistency between datasets?",
                       help="If True, run dials.cosym on all data in the "
                       "data preparation step, to ensure consistent indexing.",
                       )

        group.addParam("outlierRejection", pwprot.EnumParam,
                       label="Outlier rejection",
                       choices=["standard", "simple"],
                       default=STANDARD,
                       help="Choice of outlier rejection routine. Standard "
                       "may take a significant amount of time to run for "
                       "large datasets or high multiplicities, whereas simple"
                       " should be quick for these datasets.",
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       )

        group.addParam("outlierZmax", pwprot.FloatParam,
                       label="Outlier z-max",
                       default=6.0,
                       allowsNull=True,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       )

        group = form.addGroup("Filtering")

        group.addParam("filteringMethod", pwprot.EnumParam,
                       label="Filtering method",
                       choices=["None", "delta CC(1/2)"],
                       default=NONE,
                       display=pwprot.EnumParam.DISPLAY_HLIST,
                       help="Optionally use a filter to remove some "
                       "datasets or groups of images based on the contribution"
                       " to CC(1/2). CC(1/2) will be calculated for all "
                       "datasets or groups of images combined. Each dataset "
                       "or group will then be removed and a new CC(1/2) "
                       "calculated and compared to the original, giving a "
                       "ranking of all individual contributions to the overall"
                       " CC(1/2). The mean and the standard deviation of these"
                       " individual contributions are calculated. All datasets "
                       "or image groups that reduce CC(1/2) with more than a "
                       "number of standard deviations (default 4.0) from the "
                       "mean will be removed. The cycle is then repeated by "
                       "comparing to the overall CC(1/2) of the remaining "
                       "datasets until no datasets are removed or further "
                       "removals would violate the limits set by the user "
                       "(max cycles, max percent removed, and minimum "
                       "completeness).",
                       )

        group.addParam("ccHalfMaxCycles", pwprot.IntParam,
                       label="Max cycles",
                       default=6,
                       GE=1,
                       allowsNull=True,
                       condition=f"filteringMethod=={DELTA_CC_HALF}",
                       )

        group.addParam("ccHalfMaxPercentRemoved", pwprot.FloatParam,
                       label="Max percent removed",
                       default=10,
                       allowsNull=True,
                       condition=f"filteringMethod=={DELTA_CC_HALF}",
                       )

        group.addParam("ccHalfMinCompleteness", pwprot.FloatParam,
                       label="Minimum completeness in percent",
                       default=None,
                       GE=1,
                       LE=100,
                       allowsNull=True,
                       condition=f"filteringMethod=={DELTA_CC_HALF}",
                       )

        group.addParam("ccHalfMode", pwprot.EnumParam,
                       label="Mode",
                       choices=["dataset", "image group"],
                       default=DATASET,
                       display=pwprot.EnumParam.DISPLAY_HLIST,
                       help="Perform analysis on whole datasets or "
                       "batch groups",
                       condition="filteringMethod==1",
                       )

        group.addParam("ccHalfGroupSize", pwprot.IntParam,
                       label="Group size",
                       default=10,
                       GE=1,
                       allowsNull=True,
                       condition=f"ccHalfMode=={IMAGE_GROUP}",
                       help="The number of images to group together "
                       "when calculating delta CC(1/2) in image_group "
                       "mode",
                       )

        group.addParam("ccHalfStdcutoff", pwprot.FloatParam,
                       label="Std cutoff",
                       default=4.0,
                       allowsNull=True,
                       help="Datasets with a ΔCC½ below (mean(ΔCC½) - "
                       "(std cutoff)*standard_deviation(ΔCC½)) are removed",
                       condition="filteringMethod==1",
                       )

        # Select which images to exclude
        ImageExclusions._defineExcludeParams(self, form, inputsetsLabel)

        # Allow an easy way to import a phil file with parameters
        PhilBase._definePhilParams(self, form)

        # Allow adding anything else with command line syntax
        CliBase._defineCliParams(self, form)

   # -------------------------- INSERT functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(
            "convertInputStep",
            [inputSet.get().getObjId() for inputSet in self.inputSets])
        self._insertFunctionStep("scaleStep")
        if self.showReport:
            self._insertFunctionStep("showHtmlReportStep")
        self._insertFunctionStep("createOutputStep")

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

    def scaleStep(self):
        program = "dials.scale"
        arguments = self._prepareCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    def showHtmlReportStep(self):
        try:
            dutils._showHtmlReport(self.getOutputHtmlFile())
        except:
            self.info(self.getError())

    def createOutputStep(self):
        # Check that the indexing created proper output
        dutils.verifyPathExistence(self.getOutputReflFile(),
                                   self.getOutputModelFile())

        outputSet = self._createSetOfIndexedSpots()
        outputSet.setDialsModel(self.getOutputModelFile())
        outputSet.setDialsRefl(self.getOutputReflFile())
        try:
            outputSet.setDialsHtml(self.getOutputHtmlFile())
        except FileNotFoundError:
            self.info(self.getError())

        outputSet.write()

        self._defineOutputs(outputScaledSpots=outputSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if self.swappedResolution():
            errors.append(self.resSwapMsg())
        return errors

    def _summary(self):
        summary = []

        if self.getDatasets() not in (None, ""):
            summary.append(self.getDatasets())

        nSets = len(self.inputSets)
        if nSets > 1:
            summary.append(f"\nScaled {nSets} different datasets together")
        elif nSets == 1:
            summary.append("\nScaled a single dataset")

        if self.getDMin() is not None:
            summary.append(f"High resolution cutoff at {self.getDMin()} Å")

        if self.getDMax() is not None:
            summary.append(f"Low resolution cutoff at {self.getDMax()} Å")

        if self.checkConsistentIndexing:
            summary.append(
                "Reindexed all datasets with dials.cosym before scaling")

        if self.filteringMethod.get() is DELTA_CC_HALF:
            if self.ccHalfMode.get() is DATASET:
                mode = "datasets"
            elif self.ccHalfMode.get() is IMAGE_GROUP:
                mode = "image groups"

            summary.append(f"Filtered {mode} based on ΔCC½ with std "
                           f"cutoff {self.ccHalfStdcutoff.get()}")
        if self.excludeImages:
            for iG in ImageExclusions.getImageExclusions():
                summary.append(f"Excluded images {iG.get()}")
        if self._getCLI() != "":
            summary.append(f"Additional command line input:\n"
                           f"{self._getCLI().strip()}")

        if self.getSpaceGroupLogOutput() not in (None, ""):
            summary.append(self.getSpaceGroupLogOutput())

        if self.mtzExport("merged_mtz"):
            summary.append(
                f"Exported merged mtz file {self.getMergedMtzPath()}")

        if self.mtzExport("unmerged_mtz"):
            summary.append(
                f"Exported unmerged mtz file {self.getUnmergedMtzPath()}")

        return summary
    # -------------------------- BASE methods to be overridden -----------------

    INPUT_EXPT_FILENAME = "symmetrized.expt"
    OUTPUT_EXPT_FILENAME = "scaled.expt"
    INPUT_REFL_FILENAME = "symmetrized.refl"
    OUTPUT_REFL_FILENAME = "scaled.refl"
    OUTPUT_HTML_FILENAME = "dials.scale.html"

    def getOutputHtmlFile(self):
        return self._getExtraPath(self.OUTPUT_HTML_FILENAME)

    def getDatasets(self):
        return dutils.getDatasets(self.getOutputModelFile())

    def getCrystalName(self):
        return self.crystalName.get()

    def getLogOutput(self):
        logOutput = ""
        if self.getSpaceGroupLogOutput() not in (None, ""):
            logOutput += self.getSpaceGroupLogOutput()

        if self.getMergingStatisticsLogOutput() not in (None, ""):
            logOutput += self.getMergingStatisticsLogOutput()
        return logOutput

    def _initialParams(self, program):
        # Base method that can more easily be overridden when needed
        params = (f"{self.getAllInputFiles()} "
                  f"output.log={self.getLogFilePath(program)} "
                  f"output.experiments={self.getOutputModelFile()} "
                  f"output.reflections={self.getOutputReflFile()} "
                  f"output.html={self.getOutputHtmlFile()} "
                  f"{self.getMtzLine()}"
                  f"filtering.output.scale_and_filter_results="
                  f"{self.getOutputScaleJson()}")

        return params

    def _extraParams(self):
        params = ""
        if self.getDMin():
            params += f" cut_data.d_min={self.getDMin()}"

        if self.getDMax():
            params += f" cut_data.d_max={self.getDMax()}"

        if self.partialityCutoff.get():
            params += (f" cut_data.partiality_cutoff="
                       f"{self.partialityCutoff.get()}")

        if self.minIsigi.get():
            params += f" cut_data.min_isigi={self.minIsigi.get()}"
        # Scaling options

        if self.checkConsistentIndexing.get():
            params += (f" scaling_options.check_consistent_indexing="
                       f"{self.checkConsistentIndexing.get()}")

        if self.outlierRejection.get() is STANDARD:
            params += " outlier_rejection=standard"
        elif self.outlierRejection.get() is SIMPLE:
            params += " outlier_rejection=simple"

        if self.outlierZmax.get():
            params += f" outlier_zmax={self.outlierZmax.get()}"

        # Filtering

        if self.filteringMethod.get() is DELTA_CC_HALF:
            params += " filtering.method=deltacchalf"
        elif self.filteringMethod.get() is NONE:
            params += " filtering.method=None"

        if self.ccHalfMaxCycles.get():
            params += (f" filtering.deltacchalf.max_cycles="
                       f"{self.ccHalfMaxCycles.get()}")

        if self.ccHalfMaxPercentRemoved.get():
            params += (f" filtering.deltacchalf.max_percent_removed="
                       f"{self.ccHalfMaxPercentRemoved.get()}")

        if self.ccHalfMinCompleteness.get():
            params += (f" filtering.deltacchalf.min_completeness="
                       f"{self.ccHalfMinCompleteness.get()}")

        if self.ccHalfMode.get() is DATASET:
            params += " filtering.deltacchalf.mode=dataset"
        elif self.ccHalfMode.get() is IMAGE_GROUP:
            params += " filtering.deltacchalf.mode=image_group"

        if self.ccHalfGroupSize.get():
            params += (f" filtering.deltacchalf.group_size="
                       f"{self.ccHalfGroupSize.get()}")

        if self.ccHalfStdcutoff.get():
            params += (f" filtering.deltacchalf.stdcutoff="
                       f"{self.ccHalfStdcutoff.get()}")

        if self.excludeImages:
            for iG in ImageExclusions.getImageExclusions():
                params += f" exclude_images={iG.get()}"
        return params

    # -------------------------- UTILS functions ------------------------------

    def getOutputScaleJson(self):
        return self._getExtraPath("scale_and_filter_results.json")

    def getPhilPath(self):
        return self._getTmpPath("scale.phil")

    def getAllInputFiles(self):
        files = ""
        for iS in self.inputSets:
            files += (f"{self.getInputModelFile(iS.get())} "
                      f"{self.getInputReflFile(iS.get())} ")
        return files.strip()

    def getSpaceGroupLogOutput(self):
        spaceGroup = dutils.readLog(
            self.getLogFilePath(program="dials.scale"),
            "Space group being used",
            "Scaling models have been initialised")
        return spaceGroup

    def getMergingStatisticsLogOutput(self):
        mergingStats = dutils.readLog(
            self.getLogFilePath(program="dials.scale"),
            "Merging statistics",
            "Writing html report")
        if mergingStats not in (None, ""):
            mergeStats = f"\n{textwrap.dedent(mergingStats)}"
        else:
            mergeStats = mergingStats
        return mergeStats

    def getMtzDir(self):
        # Sanitize directory
        if self.specifyExportPath.get():
            exportPath = self.exportPath.get()
            if exportPath in (None, ""):
                return self._getExtraPath()
            elif dutils.existsPath(exportPath):
                return exportPath
        return self._getExtraPath()

    def getMtzPath(self, fn=None):
        if fn:
            return dutils.joinPath(self.getMtzDir(), fn)
        return self.getMtzDir()

    def getMergedMtzPath(self):
        # Useful for getMergedMtzLine()
        # Separate function to allow use in testing and createOutputStep
        return self.getMtzPath(self.mergedMtzName.get())

    def getUnmergedMtzPath(self):
        # Useful for getUnmergedMtzLine()
        # Separate function to allow use in testing and createOutputStep
        return self.getMtzPath(self.unmergedMtzName.get())

    def getMergedMtzLine(self):
        # Should have been part of _extraParams(), but want all output options
        # to be part of _initialParams()
        if self.mtzExport("merged_mtz"):
            return f"output.merged_mtz={self.getMergedMtzPath()} "
        else:
            return ""

    def getUnmergedMtzLine(self):
        # Should have been part of _extraParams(), but want all output options
        # to be part of _initialParams()
        if self.mtzExport("unmerged_mtz"):
            return f"output.unmerged_mtz={self.getUnmergedMtzPath()} "
        else:
            return ""

    def getMtzLine(self):
        exportMtz = (self.mtzExport("merged_mtz")
                     or self.mtzExport("unmerged_mtz"))
        mtzLine = (f"{self.getMergedMtzLine()}"
                   f"{self.getUnmergedMtzLine()}")
        if exportMtz:
            mtzLine += f"output.crystal_name={self.getCrystalName()} "
            mtzLine += f"output.project_name={self.getProjectName()} "

        return mtzLine

    def mtzExport(self, key=None):
        exportStatus = {"merged_mtz": self.exportMergedMtz.get(),
                        "unmerged_mtz": self.exportUnmergedMtz.get()}
        if key in exportStatus:
            return exportStatus[key]
        else:
            return exportStatus

    def getDMax(self):
        return self.dMax.get()

    def getDMin(self):
        return self.dMin.get()

    def swappedResolution(self):
        # d_min (high resolution) should always be smaller than d_max
        # (low resolution).
        if self.getDMin() is not None and self.getDMax() is not None:
            # Check for the case where both d_min and d_max are set
            # and have wrong relative values
            return self.getDMin() > self.getDMax()
        else:
            # If at least one value is None, then no swap is possible
            return False


class ImageExclusions(DialsProtBase):
    def _defineExcludeParams(self, form, inputsetsLabel="the input set"):
        group = form.addGroup("Selections")

        group.addParam("excludeImages", pwprot.BooleanParam,
                       label="Do you want to exclude images from a dataset?",
                       default=False,
                       help="",
                       )

        group.addParam("numberOfExclusions", pwprot.IntParam,
                       label=("How many groups of images do you want to "
                              "exclude?"),
                       default=0,
                       condition="excludeImages",
                       help="If you want to use more than 20 groups, you will "
                       "need to add it as command line parameters under "
                       "advanced options",
                       )

        group.addParam("imageGroup1", pwprot.StringParam,
                       label="Image group 1",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions "
                       "in range(1,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup2", pwprot.StringParam,
                       label="Image group 2",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(2,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup3", pwprot.StringParam,
                       label="Image group 3",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(3,21)",
                       help=f"Input in the format exp:start:end\nExclude "
                       f"a range of images (start,stop) from the dataset "
                       f"with experiment identifier exp  (inclusive of "
                       f"frames start, stop). For the first dataset listed "
                       f"in {inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset listed, "
                       f"the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup4", pwprot.StringParam,
                       label="Image group 4",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(4,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically "
                       f"0. For the next it is 1, and so on.\nTo exclude images"
                       f" 22, 23 and 24 from the second dataset listed, the "
                       f"syntax is 1:22:24.",
                       )

        group.addParam("imageGroup5", pwprot.StringParam,
                       label="Image group 5",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range"
                       "(5,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup6", pwprot.StringParam,
                       label="Image group 6",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range"
                       "(6,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup7", pwprot.StringParam,
                       label="Image group 7",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range"
                       "(7, 21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup8", pwprot.StringParam,
                       label="Image group 8",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(8, 21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup9", pwprot.StringParam,
                       label="Image group 9",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range"
                       "(9,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup10", pwprot.StringParam,
                       label="Image group 10",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(10,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup11", pwprot.StringParam,
                       label="Image group 11",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(11,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup12", pwprot.StringParam,
                       label="Image group 12",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(12,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup13", pwprot.StringParam,
                       label="Image group 13",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(13,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup14", pwprot.StringParam,
                       label="Image group 14",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(14,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup15", pwprot.StringParam,
                       label="Image group 15",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(15,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup16", pwprot.StringParam,
                       label="Image group 16",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(16,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup17", pwprot.StringParam,
                       label="Image group 17",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(17,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup18", pwprot.StringParam,
                       label="Image group 18",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(18,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup19", pwprot.StringParam,
                       label="Image group 19",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(19,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

        group.addParam("imageGroup20", pwprot.StringParam,
                       label="Image group 20",
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in "
                       "range(20,21)",
                       help=f"Input in the format exp:start:end\nExclude a "
                       f"range of images (start,stop) from the dataset with "
                       f"experiment identifier exp  (inclusive of frames "
                       f"start, stop). For the first dataset listed in "
                       f"{inputsetsLabel}, the identifier exp is typically"
                       f" 0. For the next it is 1, and so on.\nTo exclude "
                       f"images 22, 23 and 24 from the second dataset "
                       f"listed, the syntax is 1:22:24.",
                       )

    def getExclusions(self):
        return self.numberOfExclusions.get()

    def getImageExclusions(self):
        imageGroups = []
        if self.imageGroup1.get() is not None and self.getExclusions() >= 1:
            imageGroups.append(self.imageGroup1)
        if self.imageGroup2.get() is not None and self.getExclusions() >= 2:
            imageGroups.append(self.imageGroup2)
        if self.imageGroup3.get() is not None and self.getExclusions() >= 3:
            imageGroups.append(self.imageGroup3)
        if self.imageGroup4.get() is not None and self.getExclusions() >= 4:
            imageGroups.append(self.imageGroup4)
        if self.imageGroup5.get() is not None and self.getExclusions() >= 5:
            imageGroups.append(self.imageGroup5)
        if self.imageGroup6.get() is not None and self.getExclusions() >= 6:
            imageGroups.append(self.imageGroup6)
        if self.imageGroup7.get() is not None and self.getExclusions() >= 7:
            imageGroups.append(self.imageGroup7)
        if self.imageGroup8.get() is not None and self.getExclusions() >= 8:
            imageGroups.append(self.imageGroup8)
        if self.imageGroup9.get() is not None and self.getExclusions() >= 9:
            imageGroups.append(self.imageGroup9)
        if self.imageGroup10.get() is not None and self.getExclusions() >= 10:
            imageGroups.append(self.imageGroup10)
        if self.imageGroup11.get() is not None and self.getExclusions() >= 11:
            imageGroups.append(self.imageGroup11)
        if self.imageGroup12.get() is not None and self.getExclusions() >= 12:
            imageGroups.append(self.imageGroup12)
        if self.imageGroup13.get() is not None and self.getExclusions() >= 13:
            imageGroups.append(self.imageGroup13)
        if self.imageGroup14.get() is not None and self.getExclusions() >= 14:
            imageGroups.append(self.imageGroup14)
        if self.imageGroup15.get() is not None and self.getExclusions() >= 15:
            imageGroups.append(self.imageGroup15)
        if self.imageGroup16.get() is not None and self.getExclusions() >= 16:
            imageGroups.append(self.imageGroup16)
        if self.imageGroup17.get() is not None and self.getExclusions() >= 17:
            imageGroups.append(self.imageGroup17)
        if self.imageGroup18.get() is not None and self.getExclusions() >= 18:
            imageGroups.append(self.imageGroup18)
        if self.imageGroup19.get() is not None and self.getExclusions() >= 19:
            imageGroups.append(self.imageGroup19)
        if self.imageGroup20.get() is not None and self.getExclusions() >= 20:
            imageGroups.append(self.imageGroup20)

        return imageGroups
