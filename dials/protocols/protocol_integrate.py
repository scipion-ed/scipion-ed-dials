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

from pwed.objects import IndexedSpot, SetOfIndexedSpots, IndexedSpot, SetOfIndexedSpots
from pwed.protocols import EdProtIntegrateSpots
from pwed.convert import find_subranges
from dials.convert import writeJson, readRefl, writeRefl, writeRefinementPhil, copyDialsFile


class DialsProtIntegrateSpots(EdProtIntegrateSpots):
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
                      label='Filter ice?', expertLevel=pwprot.LEVEL_ADVANCED,
                      help="Filter out reflections at typical ice ring resolutions before max_cell estimation.")

        form.addParam('useScanRanges', pwprot.BooleanParam,
                      label='Cut out some images with scan_ranges?', default=False,
                      help="Explicitly specify the images to be processed. Only applicable when experiment list contains a single imageset.", expertLevel=pwprot.LEVEL_ADVANCED,
                      )

        form.addParam('dMin', pwprot.FloatParam,
                      default=None,
                      allowsNull=True,
                      label="High resolution limit",
                      help="The maximum resolution limit")

        form.addParam('dMax', pwprot.FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Low resolution limit",
                      help="The minimum resolution limit")

        # Allow adding anything else with command line syntax
        group = form.addGroup('Raw command line input parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED)
        group.addParam('commandLineInput', pwprot.StringParam,
                       default='',
                       help="Anything added here will be added at the end of the command line")

   # -------------------------- INSERT functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(
            'convertInputStep', self.inputSet.getObjId())
        self._insertFunctionStep('integrateStep')
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

    def integrateStep(self):
        program = 'dials.integrate'
        arguments = self._prepCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())
    # TODO: Create a temporary "SetOfIndexedSpotsFile" that only saves the file location

    def createOutputStep(self):
        # Check that the indexing created proper output
        assert(os.path.exists(self.getOutputReflFile()))
        assert(os.path.exists(self.getOutputModelFile()))

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
                "createOutputStep created an exception with the message {}".format(e))

        outputSet.write()

        self._defineOutputs(outputIntegratedSpots=outputSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        datasets = dutils.getModelDataPath(self.getInputModelFile())
        if len(datasets) >= 1:
            newlineDatasets = '\n'.join(
                '{}'.format(item) for item in datasets)
            summary.append(
                "Source of data:\n{}".format(newlineDatasets))

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputModelFile(self):
        if self.getSetModel():
            return self.getSetModel()
        else:
            return self._getExtraPath('sv_refined.expt')

    def getInputReflFile(self):
        if self.getSetRefl():
            return self.getSetRefl()
        else:
            return self._getExtraPath('sv_refined.refl')

    def getOutputModelFile(self):
        return self._getExtraPath('integrated_model.expt')

    def getOutputReflFile(self):
        return self._getExtraPath('integrated_reflections.refl')

    def getPhilPath(self):
        return self._getTmpPath('integrate.phil')

    def existsPath(self, path):
        return os.path.exists(path)

    def getSetModel(self):
        inputSet = self.inputSet.get()
        inputSet = self.inputSet.get()
        if self.existsPath(inputSet.getDialsModel()):
            return inputSet.getDialsModel()
        elif self.existsPath(inputSet.getDialsModel()):
            return inputSet.getDialsModel()
        else:
            return None

    def getSetRefl(self):
        inputSet = self.inputSet.get()
        inputSet = self.inputSet.get()
        if self.existsPath(inputSet.getDialsRefl()):
            return inputSet.getDialsRefl()
        elif self.existsPath(inputSet.getDialsRefl()):
            return inputSet.getDialsRefl()
        else:
            return None

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
        if self.useScanRanges.get() is True:
            params += " {}".format(self._createScanRanges())

        if self.nproc.get() not in (None, 1):
            params += " nproc={}".format(self.nproc.get())

        if self.doFilter_ice.get():
            params += " filter.ice_rings={}".format(
                self.doFilter_ice.get())

        if self.dMin.get():
            params += " prediction.d_min={}".format(self.dMin.get())

        if self.dMax.get():
            params += " prediction.d_max={}".format(self.dMax.get())

        if self.commandLineInput.get():
            params += " {}".format(self.commandLineInput.get())

        return params

    def _createScanRanges(self):
        # Go through the
        images = [image.getObjId() for image in self.inputImages.get()
                  if image.getIgnore() is not True]
        scanranges = find_subranges(images)
        scanrange = ' '.join('spotfinder.scan_range={},{}'.format(i, j)
                             for i, j in scanranges)
        return scanrange
