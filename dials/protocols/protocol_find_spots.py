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
from pwed.convert import find_subranges
from dials.convert import writeJson, readRefl


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
        jsonName = self._getExtraPath('imported.expt')
        self.model = writeJson(inputImages, fn=jsonName)

    def findSpotsStep(self):
        self.program = 'dials.find_spots'
        arguments = self._prepCommandline()
        self.runJob(self.program, arguments)
        self.reflectionData = readRefl(self.reflectionsPath)

    def createOutputStep(self):
        pass

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        return errors

    # -------------------------- UTILS functions ------------------------------
    def _prepCommandline(self):
        "Create the command line input to run dials programs"
        # Input basic parameters
        logPath = "{}/{}.log".format(self._getLogsPath(), self.program)
        self.reflectionsPath = "{}/strong.refl".format(self._getExtraPath())
        params = "{} output.log={} output.reflections={} {}".format(
            self.model, logPath, self.reflectionsPath, self._createScanRanges())

        # Update the command line with additional parameters
        if self.dMin.get():
            params += " spotfinder.filter.d_min={}".format(self.dMin.get())

        if self.dMax.get():
            params += " spotfinder.filter.d_max={}".format(self.dMax.get())

        if self.iceRings.get():
            params += " spotfinder.filter.ice_rings.filter={}".format(
                self.iceRings.get())

        if self.minSpotSize.get():
            params += " spotfinder.filter.min_spot_size={}".format(
                self.minSpotSize.get())

        if self.maxSpotSize.get():
            params += " spotfinder.filter.max_spot_size={}".format(
                self.maxSpotSize.get())

        if self.maxStrongPixelFraction.get():
            params += " spotfinder.filter.max_strong_pixel_fraction={}".format(
                self.maxStrongPixelFraction.get())

        if self.untrustedAreas.get():
            params += " spotfinder.filter.untrusted.circle={}".format(
                self.untrustedCircle.get())
            params += " spotfinder.filter.untrusted.rectangle={}".format(
                self.untrustedRectangle_1.get())
            params += " spotfinder.filter.untrusted.rectangle={}".format(
                self.untrustedRectangle_2.get())

        if self.thresholdIntensity.get():
            params += " spotfinder.threshold.dispersion.global_threshold={}".format(
                self.thresholdIntensity.get())

        if self.gain.get():
            params += " spotfinder.threshold.dispersion.gain={}".format(
                self.gain.get())

        if self.sigmaBackground.get():
            params += " spotfinder.threshold.dispersion.sigma_background={}".format(
                self.sigmaBackground.get())

        if self.sigmaStrong.get():
            params += " spotfinder.threshold.dispersion.sigma_strong={}".format(
                self.sigmaStrong.get())

        if self.kernelSize.get():
            params += " spotfinder.threshold.dispersion.kernel_size={},{}".format(
                self.kernelSize.get(), self.kernelSize.get())

        return params

    def _createScanRanges(self):
        # Go through the
        images = [image.getObjId() for image in self.inputImages.get()
                  if image.getIgnore() is not True]
        scanranges = find_subranges(images)
        scanrange = ' '.join('spotfinder.scan_range={},{}'.format(i, j)
                             for i, j in scanranges)
        return scanrange
