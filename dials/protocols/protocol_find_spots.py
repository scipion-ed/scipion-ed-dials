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

import pyworkflow.protocol as pwprot
import dials.utils as dutils

from pwed.objects import DiffractionImage, SetOfDiffractionImages, DiffractionSpot, SetOfSpots
from pwed.protocols import EdProtFindSpots
from pwed.convert import find_subranges
from dials.convert import writeJson, readRefl, copyDialsFile
from dials.constants import *


class DialsProtFindSpots(EdProtFindSpots):
    """ Base class to all EM protocols.
    It will contains some common functionalities.
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
                      help="Cut images after this index. Useful for removing frames with too much beam damage.",
                      )

        group = form.addGroup('Filtering')

        group.addParam('gain', pwprot.FloatParam,
                       default=None,
                       allowsNull=True,
                       label="Gain")

        group.addParam('kernelSize', pwprot.IntParam,
                       default=3,
                       label="Kernel size",
                       help="The size of the local area around the spot in which to"
                       "calculate the mean and variance. The kernel is given as a box"
                       "of size (2 * nx + 1, 2 * ny + 1) centred at the pixel.")

        group.addParam('sigmaBackground', pwprot.FloatParam,
                       default=6,
                       label='sigma background',
                       help="The number of standard deviations of the index of dispersion"
                       "(variance / mean) in the local area below which the pixel"
                       "will be classified as background.")

        group.addParam('sigmaStrong', pwprot.FloatParam,
                       default=3,
                       label="sigma strong",
                       help="The number of standard deviations above the mean in the local"
                       "area above which the pixel will be classified as strong.")

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
                       label="Minimum spot size (pixels)", help="The minimum "
                       "number of contiguous pixels for a spot to be accepted by the filtering algorithm.",
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('maxSpotSize', pwprot.IntParam,
                       default=1000,
                       label="Maximum spot size (pixels)", help="The maximum "
                       "number of contiguous pixels for a spot to be accepted by the filtering algorithm.",
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('maxStrongPixelFraction', pwprot.FloatParam,
                       default=0.25,
                       label='Max fraction strong pixels',
                       help="If the fraction of pixels in an image marked as strong is"
                       "greater than this value, throw an exception",
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('thresholdIntensity', pwprot.FloatParam,
                       default=0,
                       label='Minimum pixel intensity',
                       help='All pixels with a lower value will be considered part of the background',
                       expertLevel=pwprot.LEVEL_ADVANCED)

        group.addParam('maxSeparation', pwprot.FloatParam,
                       label='Max separation',
                       default=2,
                       help="The maximum peak-to-centroid separation (in pixels) for a spot to be accepted by the filtering algorithm.",
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       )

        group.addParam('thresholdAlgorithm', pwprot.EnumParam,
                       label='threshold algorithm',
                       choices=['dispersion', 'dispersion extended'], default=DISPERSION_EXTENDED,
                       help="",
                       )

        # Allow adding anything else with command line syntax
        group = form.addGroup('Raw command line input parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED)
        group.addParam('commandLineInput', pwprot.StringParam,
                       default='',
                       help="Anything added here will be added at the end of the command line")

        # Add a section for creating an html report
        form.addSection('HTML report')
        form.addParam('makeReport', pwprot.BooleanParam,
                      label='Do you want to create an HTML report for the output?', default=False,
                      help="",
                      )

        form.addParam('showReport', pwprot.BooleanParam,
                      label='Do you want to open the report as soon as the protocol is done?', default=False,
                      help="",
                      condition="makeReport",
                      )

        group = form.addGroup('Parameters',
                              condition="makeReport",)

        self.extDepOptions = ['remote', 'local', 'embed']
        group.addParam('externalDependencies', pwprot.EnumParam,
                       label='External dependencies: ',
                       choices=self.extDepOptions,
                       default=REMOTE,
                       help="Whether to use remote external dependencies (files relocatable but requires an internet connection), local (does not require internet connection but files may not be relocatable) or embed all external dependencies (inflates the html file size).",
                       )

        group.addParam('pixelsPerBin', pwprot.IntParam,
                       label='Pixels per bin',
                       default=40,
                       GE=1,
                       allowsNull=True,
                       )

        group.addParam('centroidDiffMax', pwprot.FloatParam,
                       label='Centroid diff max',
                       default=None,
                       allowsNull=True,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       help="Magnitude in pixels of shifts mapped to the extreme colours in the heatmap plots centroid_diff_x and centroid_diff_y",
                       )

        # Allow adding anything else with command line syntax
        group = form.addGroup('HTML report command line parameters',
                              expertLevel=pwprot.LEVEL_ADVANCED,
                              condition="makeReport",)
        group.addParam('commandLineInputReport', pwprot.StringParam,
                       default='',
                       help="Anything added here will be added at the end of the command line")

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
        self.info("Number of images: %s" % inputImages.getSize())
        fileName = self.getInputModelFile()
        try:
            if os.path.exists(inputImages.getDialsModel()):
                copyDialsFile(inputImages.getDialsModel(),
                              fn=fileName)
            else:
                self.info("Writing new input model file {}".format(fileName))
                writeJson(inputImages, fn=fileName)
        except TypeError as e:
            self.info(e)
            self.info("Writing new input model file {}".format(fileName))
            writeJson(inputImages, fn=fileName)

    def findSpotsStep(self):
        self.program = 'dials.find_spots'
        arguments = self._prepCommandline()
        self.runJob(self.program, arguments)

    def makeHtmlReportStep(self):
        prog = 'dials.report'
        arguments = self._prepCommandlineReport()
        self.runJob(prog, arguments)
        if self.showReport:
            dutils._showHtmlReport(self.getOutputHtmlFile())

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
        return errors

    def _summary(self):
        summary = []
        if self.getDatasets() not in (None, ''):
            summary.append(self.getDatasets())
            summary.append("\n")

        if self.getLogOutput() not in (None, ''):
            pass
        #    summary.append(self.getLogOutput())

        return summary

    # -------------------------- UTILS functions ------------------------------
    def getInputModelFile(self):
        return self._getExtraPath('imported.expt')

    def getOutputReflFile(self):
        return self._getExtraPath('strong.refl')

    def getOutputHtmlFile(self):
        return self._getExtraPath('dials.report.html')

    def getDatasets(self):
        return dutils.getDatasets(self.getInputModelFile())

    def getLogOutput(self):
        logOutput = dutils.readLog(
            self.getLogFilePath(),
            'Histogram',
            '---')
        return logOutput

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
        scanrange = ' '.join('spotfinder.scan_range={},{}'.format(i, j)
                             for i, j in scanranges)
        return scanrange

    def getLogFilePath(self, program='dials.find_spots'):
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        return logPath

    def _prepCommandline(self):
        "Create the command line input to run dials programs"
        # Input basic parameters
        logPath = "{}/{}.log".format(self._getLogsPath(), "dials.find_spots")
        params = "{} output.log={} output.reflections={} {}".format(
            self.getInputModelFile(), logPath, self.getOutputReflFile(), self.getScanRanges())

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

        if self.maxSeparation.get():
            params += " spotfinder.filter.max_separation={}".format(
                self.maxSeparation.get())

        if self.untrustedAreas.get():
            if self.untrustedCircle.get() != '':
                params += " spotfinder.filter.untrusted.circle={}".format(
                    self.untrustedCircle.get())
            if self.untrustedRectangle_1.get() != '':
                params += " spotfinder.filter.untrusted.rectangle={}".format(
                    self.untrustedRectangle_1.get())
            if self.untrustedRectangle_2.get() != '':
                params += " spotfinder.filter.untrusted.rectangle={}".format(
                    self.untrustedRectangle_2.get())

        if self.thresholdAlgorithm.get() is DISPERSION:
            params += " spotfinder.threshold.algorithm=dispersion"
        elif self.thresholdAlgorithm.get() is DISPERSION_EXTENDED:
            params += " spotfinder.threshold.algorithm=dispersion_extended"

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

        if self.commandLineInput.get():
            params += " {}".format(self.commandLineInput.get())

        return params

    def _createScanRanges(self):
        # FIXME: Remove if getScanRanges() works
        images = [image.getObjId() for image in self.inputImages.get()
                  if image.getIgnore() is not True]
        scanranges = find_subranges(images)
        scanrange = ' '.join('spotfinder.scan_range={},{}'.format(i, j)
                             for i, j in scanranges)
        return scanrange

    def _prepCommandlineReport(self):
        "Create the command line input to run dials programs"
        # Input basic parameters
        params = "{} output.html={} output.external_dependencies={}".format(
            self.getOutputReflFile(),
            self.getOutputHtmlFile(),
            self.extDepOptions[self.externalDependencies.get()]
        )

        if self.pixelsPerBin.get():
            params += " pixels_per_bin={}".format(self.pixelsPerBin.get())

        if self.centroidDiffMax.get():
            params += " centroid_diff_max={}".format(
                self.centroidDiffMax.get())

        if self.commandLineInputReport.get() not in (None, ''):
            params += " {}".format(self.commandLineInputReport.get())

        return params
