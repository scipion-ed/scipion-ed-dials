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
from dials.convert import writeJson


class DialsProtFindSpots(EdProtFindSpots):
    """ Base class to all EM protocols.
    It will contains some common functionalities.
    """
    _label = 'find spots'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        #EdProtFindSpots._defineParams(self, form)
        form.addSection(label='Input')

        form.addParam('inputImages', pwprot.PointerParam,
                      pointerClass='SetOfDiffractionImages',
                      label="Input diffraction images",
                      help="")

        # Help messages are copied from the DIALS documentation at
        # https://dials.github.io/documentation/programs/dials_find_spots.html
        form.addParam('dMin', pwprot.FloatParam,
                      default=None,
                      allowsNull=True,
                      label="High resolution limit",
                      help="The high resolution limit in Angstrom for a pixel "
                      "to be accepted by the filtering algorithm.")

        form.addParam('dMax', pwprot.FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Low resolution limit",
                      help="The low resolution limit in Angstrom for a pixel to"
                      " be accepted by the filtering algorithm.")

        form.addSection(label='Filtering')

        form.addParam('gain', pwprot.FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Gain")

        form.addParam('kernelSize', pwprot.IntParam,
                      default=3,
                      label="Kernel size",
                      help="The size of the local area around the spot in which to"
                      "calculate the mean and variance. The kernel is given as a box"
                      "of size (2 * nx + 1, 2 * ny + 1) centred at the pixel.")

        form.addParam('sigmaBackground', pwprot.FloatParam,
                      default=6,
                      label='sigma background',
                      help="The number of standard deviations of the index of dispersion"
                      "(variance / mean) in the local area below which the pixel"
                      "will be classified as background.")

        form.addParam('sigmaStrong', pwprot.FloatParam,
                      default=3,
                      label="sigma strong",
                      help="The number of standard deviations above the mean in the local"
                      "area above which the pixel will be classified as strong.")

        form.addParam('iceRings', pwprot.BooleanParam,
                      default=False, label='Filter out ice rings? ')

        form.addParam('untrustedAreas', pwprot.BooleanParam,
                      default=False, label='Are there untrusted areas? ',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('untrustedCircle', pwprot.StringParam,
                      condition='untrustedAreas', default='', help="An untrusted circle (xc, yc, r)")

        form.addParam('untrustedRectangle_1', pwprot.StringParam,
                      condition='untrustedAreas', default='', help="An untrusted rectangle (x0, x1, y0, y1)")

        form.addParam('untrustedRectangle_2', pwprot.StringParam,
                      condition='untrustedAreas', default='', help="A second untrusted rectangle (x0, x1, y0, y1)")
        form.addParam('minSpotSize', pwprot.IntParam,
                      default=None,
                      allowsNull=True,
                      label="Minimum spot size (pixels)", help="The minimum "
                      "number of contiguous pixels for a spot to be accepted by the filtering algorithm.",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('maxSpotSize', pwprot.IntParam,
                      default=1000,
                      label="Maximum spot size (pixels)", help="The maximum "
                      "number of contiguous pixels for a spot to be accepted by the filtering algorithm.",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('maxStrongPixelFraction', pwprot.FloatParam,
                      default=0.25,
                      label='Max fraction strong pixels',
                      help="If the fraction of pixels in an image marked as strong is"
                      "greater than this value, throw an exception",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('thresholdIntensity', pwprot.FloatParam,
                      default=0,
                      label='Minimum pixel intensity',
                      help='All pixels with a lower value will be considered part of the background',
                      expertLevel=pwprot.LEVEL_ADVANCED)

    # -------------------------- INSERT functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(
            'convertInputStep', self.inputImages.getObjId())
        self._insertFunctionStep('findSpotsStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, inputId):
        inputImages = self.inputImages.get()
        self.info("Number of images: %s" % inputImages.getSize())
        self.model = writeJson(inputImages)

    def findSpotsStep(self):
        pass

    def createOutputStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        return errors

    # -------------------------- UTILS functions ------------------------------
    def _prepCommandline(self):
        "Create the command line input to run dials programs"
        # Input experiment model file
        params = f" {self.model}"

        # Update the command line with additional parameters
        if self.dMin:
            params += f" filter.d_min={self.dMin}"

        if self.dMax:
            params += f" filter.d_max={self.dMax}"

        if self.iceRings:
            params += f" filter.ice_rings.filter={self.iceRings}"

        if self.minSpotSize:
            params += f" filter.min_spot_size={self.minSpotSize}"

        if self.maxSpotSize:
            params += f" filter.min_spot_size={self.maxSpotSize}"

        if self.maxStrongPixelFraction:
            params += f" filter.max_strong_pixel_fraction={self.maxStrongPixelFraction}"

        if self.untrustedAreas:
            params += f" filter.untrusted.circle={self.untrustedCircle}"
            params += f" filter.untrusted.rectangle={self.untrustedRectangle_1}"
            params += f" filter.untrusted.rectangle={self.untrustedRectangle_2}"

        if self.thresholdIntensity:
            params += f" threshold.global_threshold={self.thresholdIntensity}"

        if self.gain:
            params += f" threshold.gain={self.gain}"

        if self.sigmaBackground:
            params += f" threshold.sigma_background={self.sigmaBackground}"

        if self.sigmaStrong:
            params += f" threshold.sigma_strong={self.sigmaStrong}"

        if self.kernelSize:
            params += f" threshold.kernel_size={self.kernelSize},{self.kernelSize}"

        return params
