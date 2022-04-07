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

from pwed.objects import IndexedSpot
from pwed.protocols import EdProtIntegrateSpots
from dials.protocols import DialsProtBase, PhilBase, CliBase, HtmlBase
from pwed.convert import find_subranges
from pwed.utils import CutRes
from dials.convert import readRefl
from dials.constants import *


class DialsProtIntegrateSpots(EdProtIntegrateSpots, DialsProtBase, CutRes):
    """ Protocol for integrating spots using Dials
    """
    _label = 'integrate'

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        # EdProtIndexSpots._defineParams(self, form)

        # The start of the actually relevant part.
        form.addSection(label='Input')

        form.addParam('inputSet', pwprot.PointerParam,
                      pointerClass='SetOfIndexedSpots',
                      label="Input spots",
                      help="")

        # Help messages are copied from the DIALS documentation at
        # https://dials.github.io/documentation/programs/dials_integrate.html
        form.addParam('nproc', pwprot.IntParam,
                      label="How many processors do you want to use?",
                      default=1,
                      help="The number of processes to use.")

        form.addParam('doFilter_ice', pwprot.BooleanParam, default=False,
                      label='Filter ice?',
                      expertLevel=pwprot.LEVEL_ADVANCED,
                      help="Filter out reflections at typical ice ring "
                      "resolutions before max_cell estimation.")

        form.addParam('useScanRanges', pwprot.BooleanParam,
                      label='Cut out some images with scan_ranges?',
                      default=False,
                      help="Explicitly specify the images to be processed. "
                      "Only applicable when experiment list contains a single"
                      " imageset.",
                      expertLevel=pwprot.LEVEL_ADVANCED,
                      )

        # Define d_min and d_max
        self._defineResolutionParams(form)

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
        self._insertFunctionStep('integrateStep')
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

    def integrateStep(self):
        program = 'dials.integrate'
        arguments = self._prepareCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())
    # TODO: Create a temporary "SetOfIndexedSpotsFile" that
    # only saves the file location

    def createOutputStep(self):
        # Check that the indexing created proper output
        dutils.verifyPathExistence(self.getOutputReflFile(),
                                   self.getOutputModelFile())

        outputSet = self._createSetOfIndexedSpots()
        outputSet.setDialsModel(self.getOutputModelFile())
        outputSet.setDialsRefl(self.getOutputReflFile())

        try:
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
                iSpot.setPanel(reflDict['panel'][i])
                iSpot.setXyzobsPxValue(reflDict['xyzobs.px.value'][i])
                iSpot.setXyzobsPxVariance(reflDict['xyzobs.px.variance'][i])
                outputSet.append(iSpot)
        except Exception as e:
            self.info(
                f"createOutputStep created an exception with the message {e}")

        outputSet.write()

        self._defineOutputs(outputIntegratedSpots=outputSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if self.swappedResolution():
            errors.append(self.resSwapMsg())
        return errors

    def _summary(self):
        summary = []
        if self.getDatasets() not in (None, ''):
            summary.append(self.getDatasets())
            summary.append("\n")

        return summary

    # -------------------------- BASE methods to be overridden -----------------

    INPUT_EXPT_FILENAME = 'sv_refined.expt'
    OUTPUT_EXPT_FILENAME = 'integrated_model.expt'
    INPUT_REFL_FILENAME = 'sv_refined.refl'
    OUTPUT_REFL_FILENAME = 'integrated_reflections.refl'

    def getLogOutput(self):
        logOutput = dutils.readLog(
            self.getLogFilePath(program='dials.integrate'),
            'Summary vs resolution',
            'Timing')
        return logOutput

    def _initialParams(self, program):
        # Add output.phil parameter
        params = (f"{self.getInputModelFile()} {self.getInputReflFile()} "
                  f"output.log={self.getLogFilePath(program)} "
                  f"output.experiments={self.getOutputModelFile()} "
                  f"output.reflections={self.getOutputReflFile()} "
                  f"output.phil={self.getOutputPhilFile()}")

        return params

    def _extraParams(self):
        params = ""
        if self.useScanRanges.get() is True:
            params += f" {self._createScanRanges()}"

        if self.nproc.get() not in (None, 1):
            params += f" nproc={self.nproc.get()}"

        if self.doFilter_ice.get():
            params += f" filter.ice_rings={self.doFilter_ice.get()}"

        if self.getDMin():
            params += f" prediction.d_min={self.getDMin()}"

        if self.getDMax():
            params += f" prediction.d_max={self.getDMax()}"
        return params

    # -------------------------- UTILS functions ------------------------------

    # Placeholder for defaulting to creating phil files
    def getPhilPath(self):
        return self._getTmpPath('integrate.phil')

    def getOutputPhilFile(self):
        return self._getExtraPath('dials.integrate.phil')

    def _createScanRanges(self):
        # Go through the
        images = [image.getObjId() for image in self.inputImages.get()
                  if image.getIgnore() is not True]
        scanranges = find_subranges(images)
        scanrange = ' '.join(f'spotfinder.scan_range={i},{j}'
                             for i, j in scanranges)
        return scanrange
