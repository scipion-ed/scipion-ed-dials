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
import textwrap

import pyworkflow.protocol as pwprot
import pyworkflow.utils as pwutils
import dials.utils as dutils

from pwed.objects import IndexedSpot, SetOfIndexedSpots
from pwed.protocols import EdBaseProtocol
from dials.convert import writeJson, readRefl, writeRefl, copyDialsFile
from dials.constants import *


class DialsProtScaling(EdBaseProtocol):
    """ Protocol for integrating spots using Dials
    """
    _label = 'scale'

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        # EdProtIndexSpots._defineParams(self, form)

        # The start of the actually relevant part.
        form.addSection(label='Input')

        inputsetsLabel = 'Spots to scale'
        form.addParam('inputSets', pwprot.MultiPointerParam,
                      pointerClass='SetOfIndexedSpots',
                      label=inputsetsLabel,
                      minNumObjects=1,
                      maxNumObjects=0,
                      help="")

        form.addParam('showReport', pwprot.BooleanParam,
                      label='Do you want to view the HTML report after the processing?',
                      default=False,
                      )

        group = form.addGroup('Cut data')

        group.addParam('dMin', pwprot.FloatParam,
                       default=None,
                       allowsNull=True,
                       label="High resolution limit",
                       help="The maximum resolution limit")

        group.addParam('dMax', pwprot.FloatParam,
                       default=None,
                       allowsNull=True,
                       label="Low resolution limit",
                       help="The minimum resolution limit")

        group.addParam('partialityCutoff', pwprot.FloatParam,
                       label='Partiality cutoff',
                       default=0.4,
                       allowsNull=True,
                       help="Value below which reflections are removed from the dataset due to low partiality.",
                       expertLevel=pwprot.LEVEL_ADVANCED,)

        group.addParam('minIsigi', pwprot.FloatParam,
                       label='min I/sigma(I)',
                       default=-5,
                       allowsNull=True,
                       help="Value below which reflections are removed from the dataset due to low I/sigI in either profile or summation intensity estimates",
                       expertLevel=pwprot.LEVEL_ADVANCED,)

        group = form.addGroup('Scaling options')
        group.addParam('checkConsistentIndexing', pwprot.BooleanParam,
                       default=False,
                       label='Check indexing consistency between datasets?',
                       help="If True, run dials.cosym on all data in the data preparation step, to ensure consistent indexing.",)

        group.addParam('outlierRejection', pwprot.EnumParam,
                       label='Outlier rejection',
                       choices=['standard', 'simple'],
                       default=STANDARD,
                       help="Choice of outlier rejection routine. Standard may take a significant amount of time to run for large datasets or high multiplicities, whereas simple should be quick for these datasets.",
                       expertLevel=pwprot.LEVEL_ADVANCED,)

        group.addParam('outlierZmax', pwprot.FloatParam,
                       label="Outlier z-max",
                       default=6.0,
                       allowsNull=True,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       )

        group = form.addGroup('Filtering')

        group.addParam("filteringMethod", pwprot.EnumParam,
                       label="Filtring method",
                       choices=['None', 'delta CC(1/2)'],
                       default=NONE,
                       help="Choice of whether to do any filtering cycles, default None.",
                       )

        group.addParam("ccHalfMaxCycles", pwprot.IntParam,
                       label="Max cycles",
                       default=6,
                       GE=1,
                       allowsNull=True,
                       condition="filteringMethod=={}".format(DELTA_CC_HALF),
                       )

        group.addParam("ccHalfMaxPercentRemoved", pwprot.FloatParam,
                       label="Max percent removed",
                       default=10,
                       allowsNull=True,
                       condition="filteringMethod=={}".format(DELTA_CC_HALF),
                       )

        group.addParam("ccHalfMinCompleteness", pwprot.FloatParam,
                       label="Minimum completeness in percent",
                       default=None,
                       GE=1,
                       LE=100,
                       allowsNull=True,
                       condition="filteringMethod=={}".format(DELTA_CC_HALF),
                       )

        group.addParam("ccHalfMode", pwprot.EnumParam,
                       label="Mode",
                       choices=["dataset", "image group"],
                       default=DATASET,
                       help="Perform analysis on whole datasets or batch groups",
                       condition="filteringMethod==1",
                       )

        group.addParam("ccHalfGroupSize", pwprot.IntParam,
                       label="Group size",
                       default=10,
                       GE=1,
                       allowsNull=True,
                       condition="ccHalfMode=={}".format(IMAGE_GROUP),
                       help="The number of images to group together when calculating delta CC(1/2) in image_group mode",
                       )

        group.addParam("ccHalfStdcutoff", pwprot.FloatParam,
                       label="Std cutoff",
                       default=4.0,
                       allowsNull=True,
                       help="Datasets with a ΔCC½ below (mean - stdcutoff*std) are removed",
                       condition="filteringMethod==1",
                       )

        group = form.addGroup('Selections')

        group.addParam('excludeImages', pwprot.BooleanParam,
                       label='Do you want to exclude images from a dataset?', default=False,
                       help="",
                       )

        group.addParam('numberOfExclusions', pwprot.IntParam,
                       label='How many groups of images do you want to exclude?',
                       default=0,
                       condition="excludeImages",
                       help="If you want to use more than 20 groups, you will need to add it as command line parameters under advanced options",
                       )

        group.addParam('imageGroup1', pwprot.StringParam,
                       label='Image group 1',
                       default=None,
                       allowsNull=True,
                       condition="excludeImages and numberOfExclusions in range(1,21)",
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup2', pwprot.StringParam,
                       label='Image group 2',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(2,21)',
                       help="Input in the format exp: start: end\nExclude a range of images(start, stop) from the dataset with experiment identifier exp(inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1: 22: 24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup3', pwprot.StringParam,
                       label='Image group 3',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(3,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup4', pwprot.StringParam,
                       label='Image group 4',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(4,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup5', pwprot.StringParam,
                       label='Image group 5',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(5,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup6', pwprot.StringParam,
                       label='Image group 6',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(6,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup7', pwprot.StringParam,
                       label='Image group 7',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(7,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup8', pwprot.StringParam,
                       label='Image group 8',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(8,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup9', pwprot.StringParam,
                       label='Image group 9',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(9,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup10', pwprot.StringParam,
                       label='Image group 10',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(10,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup11', pwprot.StringParam,
                       label='Image group 11',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(11,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup12', pwprot.StringParam,
                       label='Image group 12',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(12,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup13', pwprot.StringParam,
                       label='Image group 13',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(13,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup14', pwprot.StringParam,
                       label='Image group 14',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(14,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup15', pwprot.StringParam,
                       label='Image group 15',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(15,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup16', pwprot.StringParam,
                       label='Image group 16',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(16,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup17', pwprot.StringParam,
                       label='Image group 17',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(17,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup18', pwprot.StringParam,
                       label='Image group 18',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(18,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup19', pwprot.StringParam,
                       label='Image group 19',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(19,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        group.addParam('imageGroup20', pwprot.StringParam,
                       label='Image group 20',
                       default=None,
                       allowsNull=True,
                       condition='excludeImages and numberOfExclusions in range(20,21)',
                       help="Input in the format exp:start:end\nExclude a range of images (start,stop) from the dataset with experiment identifier exp  (inclusive of frames start, stop). For the first dataset listed in {}, the identifier exp is typically 0. For the next it is 1, and so on.\nTo exclude images 22, 23 and 24 from the second dataset listed, the syntax is 1:22:24.".format(
                           inputsetsLabel),
                       )

        # Allow adding anything else with command line syntax
        group = form.addGroup('Raw command line input parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED)
        group.addParam('commandLineInput', pwprot.StringParam,
                       default='',
                       help="Anything added here will be added at the end of the command line")

   # -------------------------- INSERT functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(
            'convertInputStep', [inputSet.get().getObjId() for inputSet in self.inputSets])
        self._insertFunctionStep('scaleStep')
        if self.showReport:
            self._insertFunctionStep('showHtmlReportStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, inputSpotId):
        for iS in self.inputSets:
            inputSet = iS.get()
            self.info("ObjId: %s" %
                      inputSet.getObjId())
            self.info("Number of images: %s" % inputSet.getSize())
            self.info("Number of spots: %s" % inputSet.getSpots())
            # Write new model and/or reflection file if no was supplied from set
            if self._checkWriteModel(inputSet):
                self.writeJson(inputSet, self.getInputModelFile(inputSet))
            if self._checkWriteRefl(inputSet):
                self.writeRefl(inputSet, self.getInputReflFile(inputSet))

    def scaleStep(self):
        program = 'dials.scale'
        arguments = self._prepCommandline(program)
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
        assert(pwutils.exists(self.getOutputReflFile()))
        assert(pwutils.exists(self.getOutputModelFile()))

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
        return errors

    def _summary(self):
        summary = []

        if self.getDatasets() not in (None, ''):
            summary.append(self.getDatasets())

        nSets = len(self.inputSets)
        if nSets > 1:
            summary.append(
                '\nScaled {} different datasets together'.format(nSets))
        elif nSets == 1:
            summary.append('\nScaled a single dataset')

        if self.dMin.get() is not None:
            summary.append(
                'High resolution cutoff at {} Å'.format(self.dMin.get()))
        if self.dMax.get() is not None:
            summary.append(
                'Low resolution cutoff at {} Å'.format(self.dMax.get()))
        if self.checkConsistentIndexing:
            summary.append(
                'Reindexed all datasets with dials.cosym before scaling')
        if self.filteringMethod.get() is DELTA_CC_HALF:
            if self.ccHalfMode.get() is DATASET:
                mode = 'datasets'
            elif self.ccHalfMode.get() is IMAGE_GROUP:
                mode = 'image groups'
            summary.append('Filtered {} based on ΔCC½ with std cutoff {}'.format(
                mode, self.ccHalfStdcutoff.get()))
        if self.excludeImages:
            for iG in self.getImageExclusions():
                summary.append('Excluded images {}'.format(iG.get()))
        if self.commandLineInput.get() != '':
            summary.append('Additional command line input:\n{}'.format(
                self.commandLineInput.get()))

        if self.getSpaceGroupLogOutput() not in (None, ''):
            summary.append(self.getSpaceGroupLogOutput())

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputModelFile(self, inputSet):
        if self.getSetModel(inputSet):
            return self.getSetModel(inputSet)
        else:
            return self._getExtraPath('symmetrized.expt')

    def getInputReflFile(self, inputSet):
        if self.getSetRefl(inputSet):
            return self.getSetRefl(inputSet)
        else:
            return self._getExtraPath('symmetrized.refl')

    def getOutputModelFile(self):
        return self._getExtraPath('scaled.expt')

    def getOutputReflFile(self):
        return self._getExtraPath('scaled.refl')

    def getOutputHtmlFile(self):
        return self._getExtraPath('dials.scale.html')

    def getOutputScaleJson(self):
        return self._getExtraPath('scale_and_filter_results.json')

    def getLogFilePath(self, program='dials.scale'):
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        return logPath

    def getPhilPath(self):
        return self._getTmpPath('scale.phil')

    def getSetModel(self, inputSet):
        if pwutils.exists(inputSet.getDialsModel()):
            return inputSet.getDialsModel()
        else:
            return None

    def getSetRefl(self, inputSet):
        if pwutils.exists(inputSet.getDialsRefl()):
            return inputSet.getDialsRefl()
        else:
            return None

    def getAllInputFiles(self):
        files = ""
        for iS in self.inputSets:
            files += "{} {} ".format(self.getInputModelFile(iS.get()),
                                     self.getInputReflFile(iS.get()))
        return files.strip()

    def getDatasets(self):
        return dutils.getDatasets(self.getOutputModelFile())

    def getSpaceGroupLogOutput(self):
        spaceGroup = dutils.readLog(
            self.getLogFilePath(),
            'Space group being used',
            'Scaling models have been initialised')
        return spaceGroup

    def getMergingStatisticsLogOutput(self):
        mergingStats = dutils.readLog(
            self.getLogFilePath(),
            'Merging statistics',
            'Writing html report')
        if mergingStats not in (None, ''):
            mergeStats = "\n{}".format(textwrap.dedent(mergingStats))
        else:
            mergeStats = mergingStats
        return mergeStats

    def getLogOutput(self):
        logOutput = ''
        if self.getSpaceGroupLogOutput() not in (None, ''):
            logOutput += self.getSpaceGroupLogOutput()

        if self.getMergingStatisticsLogOutput() not in (None, ''):
            logOutput += self.getMergingStatisticsLogOutput()
        return logOutput

    def getImageExclusions(self):
        imageGroups = []
        if self.imageGroup1.get() is not None and self.numberOfExclusions.get() >= 1:
            imageGroups.append(self.imageGroup1)
        if self.imageGroup2.get() is not None and self.numberOfExclusions.get() >= 2:
            imageGroups.append(self.imageGroup2)
        if self.imageGroup3.get() is not None and self.numberOfExclusions.get() >= 3:
            imageGroups.append(self.imageGroup3)
        if self.imageGroup4.get() is not None and self.numberOfExclusions.get() >= 4:
            imageGroups.append(self.imageGroup4)
        if self.imageGroup5.get() is not None and self.numberOfExclusions.get() >= 5:
            imageGroups.append(self.imageGroup5)
        if self.imageGroup6.get() is not None and self.numberOfExclusions.get() >= 6:
            imageGroups.append(self.imageGroup6)
        if self.imageGroup7.get() is not None and self.numberOfExclusions.get() >= 7:
            imageGroups.append(self.imageGroup7)
        if self.imageGroup8.get() is not None and self.numberOfExclusions.get() >= 8:
            imageGroups.append(self.imageGroup8)
        if self.imageGroup9.get() is not None and self.numberOfExclusions.get() >= 9:
            imageGroups.append(self.imageGroup9)
        if self.imageGroup10.get() is not None and self.numberOfExclusions.get() >= 10:
            imageGroups.append(self.imageGroup10)
        if self.imageGroup11.get() is not None and self.numberOfExclusions.get() >= 11:
            imageGroups.append(self.imageGroup11)
        if self.imageGroup12.get() is not None and self.numberOfExclusions.get() >= 12:
            imageGroups.append(self.imageGroup12)
        if self.imageGroup13.get() is not None and self.numberOfExclusions.get() >= 13:
            imageGroups.append(self.imageGroup13)
        if self.imageGroup14.get() is not None and self.numberOfExclusions.get() >= 14:
            imageGroups.append(self.imageGroup14)
        if self.imageGroup15.get() is not None and self.numberOfExclusions.get() >= 15:
            imageGroups.append(self.imageGroup15)
        if self.imageGroup16.get() is not None and self.numberOfExclusions.get() >= 16:
            imageGroups.append(self.imageGroup16)
        if self.imageGroup17.get() is not None and self.numberOfExclusions.get() >= 17:
            imageGroups.append(self.imageGroup17)
        if self.imageGroup18.get() is not None and self.numberOfExclusions.get() >= 18:
            imageGroups.append(self.imageGroup18)
        if self.imageGroup19.get() is not None and self.numberOfExclusions.get() >= 19:
            imageGroups.append(self.imageGroup19)
        if self.imageGroup20.get() is not None and self.numberOfExclusions.get() >= 20:
            imageGroups.append(self.imageGroup20)

        return imageGroups

    def _checkWriteModel(self, inputSet):
        return self.getSetModel(inputSet) != self.getInputModelFile(inputSet)

    def _checkWriteRefl(self, inputSet):
        return self.getSetRefl(inputSet) != self.getInputReflFile(inputSet)

    def _prepCommandline(self, program):
        "Create the command line input to run dials programs"

        # Input basic parameters
        logPath = self.getLogFilePath(program)
        params = "{} output.log={} output.experiments={} output.reflections={} output.html={} filtering.output.scale_and_filter_results={}".format(
            self.getAllInputFiles(),
            logPath,
            self.getOutputModelFile(),
            self.getOutputReflFile(),
            self.getOutputHtmlFile(),
            self.getOutputScaleJson()
        )

        # Update the command line with additional parameters

        # Cut data

        if self.dMin.get():
            params += " cut_data.d_min={}".format(self.dMin.get())

        if self.dMax.get():
            params += " cut_data.d_max={}".format(self.dMax.get())

        if self.partialityCutoff.get():
            params += " cut_data.partiality_cutoff={}".format(
                self.partialityCutoff.get())

        if self.minIsigi.get():
            params += " cut_data.min_isigi={}".format(self.minIsigi.get())

        # Scaling options

        if self.checkConsistentIndexing.get():
            params += " scaling_options.check_consistent_indexing={}".format(
                self.checkConsistentIndexing.get())

        if self.outlierRejection.get() is STANDARD:
            params += " outlier_rejection=standard"
        elif self.outlierRejection.get() is SIMPLE:
            params += " outlier_rejection=simple"

        if self.outlierZmax.get():
            params += " outlier_zmax={}".format(self.outlierZmax.get())

        # Filtering

        if self.filteringMethod.get() is DELTA_CC_HALF:
            params += " filtering.method=deltacchalf"
        elif self.filteringMethod.get() is NONE:
            params += " filtering.method=None"

        if self.ccHalfMaxCycles.get():
            params += " filtering.deltacchalf.max_cycles={}".format(
                self.ccHalfMaxCycles.get())

        if self.ccHalfMaxPercentRemoved.get():
            params += " filtering.deltacchalf.max_percent_removed={}".format(
                self.ccHalfMaxPercentRemoved.get())

        if self.ccHalfMinCompleteness.get():
            params += " filtering.deltacchalf.min_completeness={}".format(
                self.ccHalfMinCompleteness.get())

        if self.ccHalfMode.get() is DATASET:
            params += " filtering.deltacchalf.mode=dataset"
        elif self.ccHalfMode.get() is IMAGE_GROUP:
            params += " filtering.deltacchalf.mode=image_group"

        if self.ccHalfGroupSize.get():
            params += " filtering.deltacchalf.group_size={}".format(
                self.ccHalfGroupSize.get())

        if self.ccHalfStdcutoff.get():
            params += " filtering.deltacchalf.stdcutoff={}".format(
                self.ccHalfStdcutoff.get())

        if self.excludeImages:
            for iG in self.getImageExclusions():
                params += " exclude_images={}".format(iG.get())

        # Unrestrained command line input

        if self.commandLineInput.get():
            params += " {}".format(self.commandLineInput.get())

        return params
