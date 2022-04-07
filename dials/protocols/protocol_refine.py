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

import os

import pyworkflow.protocol as pwprot
import dials.utils as dutils

from pwed.objects import IndexedSpot
from pwed.protocols import EdProtRefineSpots
from dials.protocols import (
    DialsProtBase, PhilBase, CliBase, HtmlBase, RefineParamsBase)
from dials.convert import readRefl, writeRestraintsPhil
from dials.constants import *


class DialsProtRefineSpots(EdProtRefineSpots, DialsProtBase):
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

        # Keep old option and syntax for compatibility, but do not show it.
        form.addHidden('scanVarying', pwprot.BooleanParam,
                       label='Perform scan-varying refinement?', default=None,
                       allowsNull=True)

        form.addParam('useScanVaryingFromWorkflow', pwprot.BooleanParam,
                      label='Use scan varying option from imported workflow?',
                      default=True,
                      condition='scanVarying==True',
                      help="It appears that a workflow from an older version"
                      " used scan-varying refinement. Do you want to keep "
                      "using that option?")

        form.addParam('scanVaryingNew', pwprot.EnumParam,
                      label='Scan varying or static?',
                      choices=['Auto', 'Static',
                               'Scan-varying', 'Dials default'],
                      default=UNSET,
                      condition="scanVarying!=True or "
                      "useScanVaryingFromWorkflow==False",
                      help="Allow models that are not forced to be static to "
                      "vary during the scan, Auto will run one macrocycle with"
                      " static then scan varying refinement for the crystal. "
                      "The option \"Dials default\" will use the default as "
                      "indicated in https://dials.github.io/documentation/"
                      "programs/dials_refine.html.",
                      )

        form.addParam('useRestraint', pwprot.BooleanParam,
                      label="Do you want to restrain the unit cell?",
                      default=False,)

        form.addParam('targetUnitCell', pwprot.StringParam,
                      label="Unit cell values",
                      allowsNull=True, default=None,
                      condition='useRestraint==True',
                      help="Target unit cell parameters for the restraint for"
                      " this parameterisation")

        form.addParam('targetSigmas', pwprot.StringParam,
                      label="Sigmas for determining restraint weight",
                      allowsNull=True, default=None,
                      condition='useRestraint==True',
                      help="The unit cell target values are associated with "
                      "sigmas which are used to determine the weight of each "
                      "restraint. A sigma of zero will remove the restraint at"
                      " that position. If symmetry constrains two cell "
                      "dimensions to be equal then only the smaller of the two"
                      " sigmas will be kept")

        # Help messages are copied from the DIALS documentation at
        # https://dials.github.io/documentation/programs/dials_index.html
        form.addParam('Nproc', pwprot.IntParam,
                      label="How many processors do you want to use?",
                      default=1,
                      help="The number of processes to use. Not all choices of"
                      " refinement engine support nproc > 1. Where "
                      "multiprocessing is possible, it is helpful only in "
                      "certain circumstances, so this is not recommended for "
                      "typical use.",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        RefineParamsBase._defineParametrisations(self, form)

        group = form.addGroup('Refinery')

        group.addParam('doSetMaxIterations', pwprot.BooleanParam,
                       label="Do you want to set the maximum number of "
                       "iterations?",
                       default=False,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       help="Maximum number of iterations in refinement before"
                       " termination.",
                       )

        group.addParam('refineryMaxIterations', pwprot.IntParam,
                       default=None,
                       allowsNull=True,
                       help="Maximum number of iterations in refinement before"
                       " termination."
                       "None implies the engine supplies its own default.",
                       label='Max iterations', condition="doSetMaxIterations",
                       )

        # Allow an easy way to import a phil file with parameters
        PhilBase._definePhilParams(self, form)

        # Allow adding anything else with command line syntax
        CliBase._defineCliParams(self, form)

        # Add a section for creating an html report
        HtmlBase._defineHtmlParams(self, form)

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
        self.info(f"Number of images: {inputSet.getSize()}")
        self.info(f"Number of spots: {inputSet.getSpots()}")
        # Write new model and/or reflection file if no was supplied from set
        if self._checkWriteModel():
            self.writeJson(inputSet, self.getInputModelFile())
        if self._checkWriteRefl():
            self.writeRefl(inputSet, self.getInputReflFile())

    def refineStep(self):
        program = 'dials.refine'
        arguments = self._prepareCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    def createOutputStep(self):
        # Check that the indexing created proper output
        dutils.verifyPathExistence(self.getOutputModelFile(),
                                   self.getOutputReflFile())

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

    # -------------------------- BASE methods to be overridden -----------------

    INPUT_EXPT_FILENAME = 'indexed.expt'
    OUTPUT_EXPT_FILENAME = 'refined.expt'
    INPUT_REFL_FILENAME = 'indexed.refl'
    OUTPUT_REFL_FILENAME = 'refined.refl'

    def _extraParams(self):
        params = ""
        if self.Nproc.get() not in (None, 1):
            params += f" nproc={self.Nproc.get()}"

        params += self.getScanVaryingCommand()

        params += RefineParamsBase.getBeamFixParams(self)

        if self.beamForceStatic.get() not in (None, True):
            if self.getScanVaryingStatus():
                params += f" beam.force_static={self.beamForceStatic.get()}"

        params += RefineParamsBase.getCrystalFixParams(self)

        params += RefineParamsBase.getDetectorFixParams(self)

        params += RefineParamsBase.getGonioFixParams(self)

        if self.refineryMaxIterations.get() is not None:
            params += (f" refinery.max_iterations="
                       f"{self.refineryMaxIterations.get()}")

        if self.useRestraint:
            # Make the phil when it is known that it will be used
            self.makeRestraintsPhil()
            params += f" {self.getRestraintsPhil()}"

        return params

    # -------------------------- UTILS functions ------------------------------

    def getScanVaryingCommand(self):
        if self.scanVarying.get() is True and self.useScanVaryingFromWorkflow.get() is True:
            return " scan_varying=True"
        elif self.scanVaryingNew.get() == SCAN_VARYING:
            return " scan_varying=True"
        elif self.scanVaryingNew.get() == STATIC:
            return " scan_varying=False"
        elif self.scanVaryingNew.get() == AUTO:
            return " scan_varying=Auto"
        elif self.scanVaryingNew.get() == UNSET:
            return ""
        return ""

    def getScanVaryingStatus(self):
        try:
            return self.getScanVaryingCommand().split("=")[-1].lower() == "true"
        except:
            return None

    def getRestraintsPhil(self):
        return self._getExtraPath("restraints.phil")

    def makeRestraintsPhil(self):
        writeRestraintsPhil(
            fn=self.getRestraintsPhil(),
            values=self.targetUnitCell.get(),
            sigmas=self.targetSigmas.get()
        )

    def _checkWriteModel(self):
        return self.getSetModel() != self.getInputModelFile()

    def _checkWriteRefl(self):
        return self.getSetRefl() != self.getInputReflFile()

    def _prepCommandline(self, program):
        "Create the command line input to run dials programs"

        # Input basic parameters
        logPath = f"{self._getLogsPath()}/{program}.log"
        params = (f"{self.getInputModelFile()} {self.getInputReflFile()} "
                  f"output.log={logPath} "
                  f"output.experiments={self.getOutputModelFile()}"
                  f" output.reflections={self.getOutputReflFile()}")

        # Update the command line with additional parameters

        if self.Nproc.get() not in (None, 1):
            params += f" nproc={self.Nproc.get()}"

        params += self.getScanVaryingCommand()

        if self.beamFixAll:
            params += (" beam.fix='*all in_spindle_plane"
                       " out_spindle_plane wavelength'")
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
            params += f" beam.fix={''.join(beamfix)}"

        if self.beamForceStatic.get() not in (None, True):
            if self.getScanVaryingStatus():
                params += f" beam.force_static={self.beamForceStatic.get()}"

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
            params += f" crystal.fix={''.join(crystalfix)}"

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
        params += f" detector.fix={''.join(detectorfix)}"

        if self.refineryMaxIterations.get() is not None:
            params += (f" refinery.max_iterations="
                       f"{self.refineryMaxIterations.get()}")

        if self.useRestraint:
            # Make the phil when it is known that it will be used
            self.makeRestraintsPhil()
            params += f" {self.getRestraintsPhil()}"

        if self.extraPhilPath.get():
            params += f" {self.getExtraPhilsPath()}"

        if self.commandLineInput.get():
            params += f" {self.commandLineInput.get()}"

        return params
