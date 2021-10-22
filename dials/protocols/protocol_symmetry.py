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

import pyworkflow.protocol as pwprot
import dials.utils as dutils

from dials.protocols import DialsProtBase, PhilBase, CliBase, HtmlBase
from dials.constants import *
import dials.convert as dconv


class DialsProtSymmetry(DialsProtBase, HtmlBase):
    """ Protocol for checking symmetry of integrated spots using the
    POINTLESS algorithm as implemented in DIALS
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

        # Allow an easy way to import a phil file with parameters
        PhilBase._definePhilParams(self, form)

        # Allow adding anything else with command line syntax
        CliBase._defineCliParams(self, form)

        # Add a section for creating an html report
        self._defineHtmlParams(form)

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
        self.info(f"Number of images: {inputSet.getSize()}")
        self.info(f"Number of spots: {inputSet.getSpots()}")
        # Write new model and/or reflection file if no was supplied from set
        if self._checkWriteModel():
            dconv.writeJson(inputSet, self.getInputModelFile())
        if self._checkWriteRefl():
            dconv.writeRefl(inputSet, self.getInputReflFile())

    def symmetryStep(self):
        program = 'dials.symmetry'
        arguments = self._prepareCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    def createOutputStep(self):
        # Check that the indexing created proper output
        dutils.verifyPathExistence(self.getOutputReflFile(),
                                   self.getOutputModelFile())

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
    # -------------------------- BASE methods to be overridden -----------------

    INPUT_EXPT_FILENAME = 'integrated.expt'
    OUTPUT_EXPT_FILENAME = 'symmetrized.expt'
    INPUT_REFL_FILENAME = 'integrated.refl'
    OUTPUT_REFL_FILENAME = 'symmetrized.refl'
    OUTPUT_HTML_FILENAME = "dials.symmetry.html"
    OUTPUT_JSON_FILENAME = "dials.symmetry.json"

    def getLogOutput(self):
        logOutput = dutils.readLog(
            self.getLogFilePath("dials.symmetry"),
            'Recommended',
            'Saving')
        return logOutput.strip()

    def _extraParams(self):
        params = ""
        params += f" output.html={self.getOutputHtmlFile()}"
        params += f" output.json={self.getOutputJsonFile()}"
        return params

    # -------------------------- UTILS functions ------------------------------

    # Prepare to use phils as default files
    def getPhilPath(self):
        return self._getTmpPath('symmetry.phil')
