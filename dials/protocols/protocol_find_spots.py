# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.protocol as pwprot

from pwed.objects import DiffractionImage, SetOfDiffractionImages
from pwed.protocols import EdProtFindSpots


class DialsProtFindSpots(EdProtFindSpots):
    """ Base class to all EM protocols.
    It will contains some common functionalities.
    """
    _base = True

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputImages', pwprot.PointerParam,
                      pointerCllass='SetOfDiffractionImages',
                      label="Input diffraction images",
                      help="")

        # form.addParam('filesPattern', pwprot.StringParam,
        #               label='Pattern',
        #               help="Pattern of the tilt series\n\n"
        #                    "The pattern can contain standard wildcards such as\n"
        #                    "*, ?, etc.\n\n"
        #                    "It should also contains the following special tags:"
        #                    "   {TS}: tilt series identifier "
        #                    "         (can be any UNIQUE part of the path).\n"
        #                    "   {TI}: image identifier "
        #                    "         (an integer value, unique within a tilt-series).\n"
        #                    "Examples:\n"
        #                    "")
        # form.addParam('importAction', pwprot.EnumParam,
        #               default=self.IMPORT_LINK_REL,
        #               choices=['Copy files',
        #                        'Absolute symlink',
        #                        'Relative symlink'],
        #               display=pwprot.EnumParam.DISPLAY_HLIST,
        #               expertLevel=pwprot.LEVEL_ADVANCED,
        #               label="Import action on files",
        #               help="By default ...")

    # -------------------------- INSERT functions ------------------------------
    def _insertAllSteps(self):
        self.loadPatterns()
        self._insertFunctionStep('convertInputStep', inputId)
        self._insertFunctionStep('findSpotsStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, inputId):
        pass

    def findSpotsStep(self):
        pass

    def createOutputStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        return errors

    # -------------------------- UTILS functions ------------------------------


