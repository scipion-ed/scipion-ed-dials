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
from dials.protocols import DialsProtBase, CliBase, PhilBase, HtmlBase
import dials.utils as dutils

from pwed.objects import DiffractionSpot
from pwed.protocols import EdProtFindSpots
from pwed.utils import CutRes
from pwed.convert import find_subranges
from dials.convert import writeJson, readRefl, copyDialsFile
from dials.constants import *


class DialsProtFindSpots(EdProtFindSpots, DialsProtBase, CutRes):
    """ Protocol for performing spot finding using dials.find_spots
    """

    _label = 'find spots'

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        # EdProtFindSpots._defineParams(self, form)
        form.addSection(label='Input')

        form.addParam('inputImages', pwprot.PointerParam,
                      pointerClass='SetOfDiffractionImages',
                      label="Input diffraction images",
                      help="")

        # Help messages are copied from the DIALS documentation at
        # https://dials.github.io/documentation/programs/dials_find_spots.html
        dminHelp = ("The high resolution limit in Angstrom for a pixel "
                    "to be accepted by the filtering algorithm.")
        dmaxHelp = ("The low resolution limit in Angstrom for a pixel to"
                    " be accepted by the filtering algorithm.")
        self._defineResolutionParams(form, dminHelp, dmaxHelp)

        form.addParam('minImage', pwprot.IntParam,
                      label='First image to use',
                      default=None,
                      allowsNull=True,
                      help="Do not use images with lower index",
                      )

        form.addParam('maxImage', pwprot.IntParam,
                      label='Last image to use',
                      default=None,
                      allowsNull=True,
                      help="Cut images after this index. Useful for "
                      "removing frames with too much beam damage.",
                      )

        group = form.addGroup('Filtering')

        group.addParam('gain', pwprot.FloatParam,
                       default=None,
                       allowsNull=True,
                       label="Gain")

        group.addParam('kernelSize', pwprot.IntParam,
                       default=3,
                       label="Kernel size",
                       help="The size of the local area around the spot in "
                       "which to calculate the mean and variance. The kernel "
                       "is given as a box of size (2 * nx + 1, 2 * ny + 1) "
                       "centred at the pixel."
                       )

        group.addParam('sigmaBackground', pwprot.FloatParam,
                       default=6,
                       label='sigma background',
                       help="The number of standard deviations of the index "
                       "of dispersion (variance / mean) in the local area "
                       "below which the pixelwill be classified as background."
                       )

        group.addParam('sigmaStrong', pwprot.FloatParam,
                       default=3,
                       label="sigma strong",
                       help="The number of standard deviations above the mean"
                       " in the localarea above which the pixel will be "
                       "classified as strong.")

        group.addParam('iceRings', pwprot.BooleanParam,
                       default=False, label='Filter out ice rings? ')

        group.addParam('untrustedAreas', pwprot.BooleanParam,
                       default=False, label='Are there untrusted areas? ',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('untrustedCircle', pwprot.StringParam,
                       condition='untrustedAreas',
                       label='Untrusted circle',
                       default='',
                       help="An untrusted circle (xc, yc, r)",
                       )

        group.addParam('untrustedRectangle_1', pwprot.StringParam,
                       condition='untrustedAreas',
                       label='Untrusted rectangle',
                       default='',
                       help="An untrusted rectangle (x0, x1, y0, y1)",
                       )

        group.addParam('untrustedRectangle_2', pwprot.StringParam,
                       condition='untrustedAreas',
                       label='Untrusted rectangle',
                       default='',
                       help="An untrusted rectangle (x0, x1, y0, y1)",
                       )

        group.addParam('minSpotSize', pwprot.IntParam,
                       default=None,
                       allowsNull=True,
                       label="Minimum spot size (pixels)",
                       help="The minimum number of contiguous pixels for a "
                       "spot to be accepted by the filtering algorithm.",
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('maxSpotSize', pwprot.IntParam,
                       default=1000,
                       label="Maximum spot size (pixels)",
                       help="The maximum "
                       "number of contiguous pixels for a spot to be accepted"
                       " by the filtering algorithm.",
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('maxStrongPixelFraction', pwprot.FloatParam,
                       default=0.25,
                       label='Max fraction strong pixels',
                       help="If the fraction of pixels in an image marked as"
                       " strong is greater than this value, throw an "
                       "exception",
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('thresholdIntensity', pwprot.FloatParam,
                       default=0,
                       label='Minimum pixel intensity',
                       help="All pixels with a lower value will be considered"
                       " part of the background",
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('maxSeparation', pwprot.FloatParam,
                       label='Max separation',
                       default=2,
                       help="The maximum peak-to-centroid separation (in "
                       "pixels) for a spot to be accepted by the filtering "
                       "algorithm.",
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       )

        group.addParam('thresholdAlgorithm', pwprot.EnumParam,
                       label='threshold algorithm',
                       choices=['dispersion', 'dispersion extended'],
                       default=DISPERSION_EXTENDED,
                       help="",
                       )

        PhilBase._definePhilParams(self, form)

        CliBase._defineCliParams(self, form)

        HtmlBase._defineHtmlParams(self, form)

    # -------------------------- INSERT functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(
            'convertInputStep', self.inputImages.getObjId())
        self._insertFunctionStep('findSpotsStep')
        if self.makeReport:
            self._insertFunctionStep('makeHtmlReportStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, inputId):
        inputImages = self.inputImages.get()
        self.info(f"Number of images: {inputImages.getSize()}")
        fileName = self.getInputModelFile()
        try:
            if os.path.exists(inputImages.getDialsModel()):
                copyDialsFile(inputImages.getDialsModel(),
                              fn=fileName)
            else:
                self.info(f"Writing new input model file {fileName}")
                writeJson(inputImages, fn=fileName)
        except TypeError as e:
            self.info(e)
            self.info(f"Writing new input model file {fileName}")
            writeJson(inputImages, fn=fileName)

    def findSpotsStep(self):
        program = 'dials.find_spots'
        arguments = self._prepareCommandline(program)
        self.runJob(program, arguments)

    def createOutputStep(self):
        reflectionData = readRefl(self.getOutputReflFile())
        outputSet = self._createSetOfSpots()
        dSpot = DiffractionSpot()
        numberOfSpots = reflectionData[2]
        reflDict = reflectionData[4]

        outputSet.setSpots(numberOfSpots)
        outputSet.setDialsModel(self.getInputModelFile())
        outputSet.setDialsRefl(self.getOutputReflFile())

        for i in range(0, numberOfSpots):
            dSpot.setObjId(i+1)
            dSpot.setSpotId(reflDict['id'][i])
            dSpot.setBbox(reflDict['bbox'][i])
            dSpot.setFlag(reflDict['flags'][i])
            dSpot.setIntensitySumValue(reflDict['intensity.sum.value'][i])
            dSpot.setIntensitySumVariance(
                reflDict['intensity.sum.variance'][i])
            dSpot.setNSignal(reflDict['n_signal'][i])
            dSpot.setPanel(reflDict['panel'][i])
            try:
                dSpot.setShoebox(reflDict['shoebox'][i])
            except IndexError:
                pass
            dSpot.setXyzobsPxValue(reflDict['xyzobs.px.value'][i])
            dSpot.setXyzobsPxVariance(reflDict['xyzobs.px.variance'][i])
            outputSet.append(dSpot)

        outputSet.write()

        self._defineOutputs(outputDiffractionSpots=outputSet)

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

        if self.getLogOutput() not in (None, ''):
            summary.append(self.getLogOutput())

        return summary
    # -------------------------- BASE methods to be overridden -----------------

    INPUT_EXPT_FILENAME = 'imported.expt'
    OUTPUT_REFL_FILENAME = 'strong.refl'

    def getLogOutput(self):
        logOutput = dutils.readLog(
            self.getLogFilePath('dials.find_spots'),
            'Histogram',
            '---')
        return logOutput

    def _initialParams(self, program):
        # Input basic parameters
        params = (f"{self.getInputModelFile()} "
                  f"output.log={self.getLogFilePath(program)} "
                  f"output.reflections={self.getOutputReflFile()}")

        return params

    def _extraParams(self):
        params = f" {self.getScanRanges()}"
        # Update the command line with additional parameters
        if self.getDMin():
            params += f" spotfinder.filter.d_min={self.getDMin()}"

        if self.getDMax():
            params += f" spotfinder.filter.d_max={self.getDMax()}"

        if self.iceRings.get():
            params += (f" spotfinder.filter.ice_rings.filter="
                       f"{self.iceRings.get()}")

        if self.minSpotSize.get():
            params += (f" spotfinder.filter.min_spot_size="
                       f"{self.minSpotSize.get()}")

        if self.maxSpotSize.get():
            params += (f" spotfinder.filter.max_spot_size="
                       f"{self.maxSpotSize.get()}")

        if self.maxStrongPixelFraction.get():
            params += (f" spotfinder.filter.max_strong_pixel_fraction="
                       f"{self.maxStrongPixelFraction.get()}")

        if self.maxSeparation.get():
            params += (f" spotfinder.filter.max_separation="
                       f"{self.maxSeparation.get()}")

        if self.untrustedAreas.get():
            if self.untrustedCircle.get() != '':
                circle = self.fixString(self.untrustedCircle.get())
                params += (f" spotfinder.filter.untrusted.circle="
                           f"{circle}")
            if self.untrustedRectangle_1.get() != '':
                rectangle1 = self.fixString(self.untrustedRectangle_1.get())
                params += (f" spotfinder.filter.untrusted.rectangle="
                           f"{rectangle1}")
            if self.untrustedRectangle_2.get() != '':
                rectangle2 = self.fixString(self.untrustedRectangle_2.get())
                params += (f" spotfinder.filter.untrusted.rectangle="
                           f"{rectangle2}")

        if self.thresholdAlgorithm.get() is DISPERSION:
            params += f" spotfinder.threshold.algorithm=dispersion"
        elif self.thresholdAlgorithm.get() is DISPERSION_EXTENDED:
            params += f" spotfinder.threshold.algorithm=dispersion_extended"

        if self.thresholdIntensity.get():
            params += (f" spotfinder.threshold.dispersion.global_threshold="
                       f"{self.thresholdIntensity.get()}")

        if self.gain.get():
            params += (f" spotfinder.threshold.dispersion.gain="
                       f"{self.gain.get()}")

        if self.sigmaBackground.get():
            params += (f" spotfinder.threshold.dispersion.sigma_background="
                       f"{self.sigmaBackground.get()}")

        if self.sigmaStrong.get():
            params += (f" spotfinder.threshold.dispersion.sigma_strong="
                       f"{self.sigmaStrong.get()}")

        if self.kernelSize.get():
            params += (f" spotfinder.threshold.dispersion.kernel_size="
                       f"{self.kernelSize.get()},{self.kernelSize.get()}")
        return params
    # -------------------------- UTILS functions ------------------------------

    def getScanRanges(self):
        # Get a list of object IDs for good images
        objIds = [image.getObjId() for image in self.inputImages.get()
                  if image.getIgnore() is not True]
        # Get where to start the scan range
        if self.minImage.get() is None:
            first = min(objIds)
        else:
            first = self.minImage.get()
        # Get where to stop the scan range
        if self.maxImage.get() is None:
            last = max(objIds)
        else:
            last = self.maxImage.get()
        # Make a list of all images to include in the processing
        images = []
        for i in objIds:
            if first <= i <= last:
                images.append(i)
        # Exclude skipped images from the scan range
        scanranges = find_subranges(images)
        scanrange = ' '.join(f'spotfinder.scan_range={i},{j}'
                             for i, j in scanranges)
        return scanrange
