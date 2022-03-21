# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              V. E.G. Bengtsson (viktor.bengtsson@mmk.su.se) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Department of Materials and Environmental Chemistry,
# * Stockholm University
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

from dials.protocols import DialsProtBase, CliBase, PhilBase
import dials.utils as dutils

from pwed.protocols import ProtImportDiffractionImages


class DialsProtImportDiffractionImages(
        ProtImportDiffractionImages,
        DialsProtBase):
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

        PhilBase._definePhilParams(self, form)

        CliBase._defineCliParams(self, form)

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
        self.info(f"Rotation axis is {self.getRotationAxis()}")
        program = 'dials.import'
        arguments = self._prepareCommandline(program)
        self.runJob(program, arguments)
        dutils.verifyPathExistence(self.getOutputModelFile())

    def createOutputStep(self):
        super().createOutputStep(dialsModel=self.getOutputModelFile())

    # -------------------------- INFO functions -------------------------------
    # def _validate(self)

    # -------------------------- BASE methods to be overridden -----------------
    OUTPUT_EXPT_FILENAME = 'imported.expt'

    def getDatasets(self):
        return dutils.getDatasets(self.getOutputModelFile())

    def _initialParams(self, program):
        params = (f"{self.getCmdparamStart()} "
                  f"output.log={self.getLogFilePath(program)} "
                  f"output.experiments={self.getOutputModelFile()}"
                  )
        return params

    def _extraParams(self):
        params = self._getDialsOverwrites()
        return params

    # -------------------------- UTILS functions ------------------------------

    def getCmdparamStart(self):
        self.loadPatterns()
        if self.useTemplate:
            return f"template={self._templatePattern}"
        else:
            fileString = " ".join([i[0] for i in self.getMatchingFiles()])
            return fileString

    def _getDialsOverwrites(self):
        params = ""
        if self.getRotationAxis():
            params += (
                f" goniometer.axes={self.list2str(self.getRotationAxis())}")

        if self.getNewDetectorDistance() is not None:
            params += f" distance={self.getNewDetectorDistance()}"
        return params
