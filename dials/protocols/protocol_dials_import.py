# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              V. E.G. Bengtsson (viktor.bengtsson@mmk.su.se) [2]
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
import pathlib
from glob import glob

import pyworkflow as pw
import pyworkflow.protocol as pwprot
import dials.utils as dutils

from pwed.objects import DiffractionImage, SetOfDiffractionImages
from pwed.protocols import ProtImportDiffractionImages


class DialsProtImportDiffractionImages(ProtImportDiffractionImages):
    """ Base class for other Import protocols.
    All imports protocols will have:
    1) Several options to import from (_getImportOptions function)
    2) First option will always be "from files". (for this option
      files with a given pattern will be retrieved  and the ### will
      be used to mark an ID part from the filename.
      - For each file a function to process it will be called
        (_importFile(fileName, fileId))
    """

    _label = 'import diffraction images'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtImportDiffractionImages._defineParams(self, form)

        # Allow adding anything else with command line syntax
        group = form.addGroup('Raw command line input parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED)
        group.addParam('commandLineInput', pwprot.StringParam,
                       default='',
                       help="Anything added here will be added at the end of the command line")

    # -------------------------- INSERT functions ------------------------------
    def _insertAllSteps(self):
        self.loadPatterns()
        super()._insertFunctionStep(
            'convertInputStep', self._pattern)
        self._insertFunctionStep('importStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def importStep(self):
        # Run dials import on the images
        self.info("Rotation axis is {}".format(self.getRotationAxis()))
        program = 'dials.import'
        arguments = self._prepareCommandLineArguments(program)
        self.runJob(program, arguments)
        assert(os.path.exists(self.getOutputModelFile()))

    def createOutputStep(self):
        super().createOutputStep(dialsModel=self.getOutputModelFile())

    # -------------------------- INFO functions -------------------------------
    # def _validate(self)

    # -------------------------- BASE methods to be overridden -----------------
    # def here

    # -------------------------- UTILS functions ------------------------------

    def getCmdparamStart(self):
        self.loadPatterns()
        if self.useTemplate:
            return "template={}".format(self._templatePattern)
        else:
            fileString = " ".join([i[0] for i in self.getMatchingFiles()])
            return fileString

    def getOutputModelFile(self):
        return self._getExtraPath('imported.expt')

    def getDatasets(self):
        return dutils.getDatasets(self.getOutputModelFile())

    def getLogOutput(self):
        return ''

    def _prepareCommandLineArguments(self, program):
        # Make a string to append to
        cmdparams = "{}".format(self.getCmdparamStart())

        # Add standard output paths
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        cmdparams += " output.log={} output.experiments={}".format(
            logPath, self.getOutputModelFile())

        if self.getRotationAxis():
            cmdparams += " goniometer.axes={}".format(
                ",".join(map(str, self.getRotationAxis())))

        if self.overwriteDetectorDistance.get() is not None:
            cmdparams += " distance={}".format(
                self.overwriteDetectorDistance.get())

        if self.commandLineInput.get():
            cmdparams += " {}".format(self.commandLineInput.get())

        return cmdparams
