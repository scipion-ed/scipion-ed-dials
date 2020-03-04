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

from pwed.objects import DiffractionImage, SetOfDiffractionImages, DiffractionSpot, SetOfSpots, IndexedSpot, SetOfIndexedSpots
from pwed.protocols import EdProtIndexSpots
from pwed.convert import find_subranges
from dials.convert import writeJson, readRefl, writeRefl, writeRefinementPhil, copyInput


class DialsProtIndexSpots(EdProtIndexSpots):
    """ Protocol for indexing spots using Dials
    """
    _label = 'index'

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        # EdProtIndexSpots._defineParams(self, form)

        # Check which parts of indexing to perform. Reindexing is probably to small to warrant
        # its own protocol.
        form.addSection(label='Processing steps')

        form.addParam('doIndex', pwprot.BooleanParam,
                      default=True, label='Do you want to index from start?')
        form.addParam('doRefineBravaisSettings', pwprot.BooleanParam,
                      default=True, label='Do you want to refine the Bravais settings?')
        form.addParam('doReindex', pwprot.BooleanParam,
                      default=True, label='Do you want to reindex after refining Bravais settings?',
                      condition='doRefineBravaisSettings')

        # The start of the actually relevant part.
        form.addSection(label='Indexing basics', condition='doIndex')

        form.addParam('inputImages', pwprot.PointerParam,
                      pointerClass='SetOfDiffractionImages',
                      label="Input diffraction images",
                      help="")

        form.addParam('inputSpots', pwprot.PointerParam,
                      pointerClass='SetOfSpots',
                      label='Input strong spots',
                      help="")

        form.addParam('inputStrongPath', pwprot.StringParam,
                      label='Give the absolute path to the directory of strong.refl and import.expt',
                      default='',
                      help="",
                      )

        # Help messages are copied from the DIALS documentation at
        # https://dials.github.io/documentation/programs/dials_index.html
        form.addParam('indexNproc', pwprot.IntParam,
                      label="How many processors do you want to use?",
                      default=1,
                      help="The number of processes to use.")

        form.addParam('enterSpaceGroup', pwprot.BooleanParam,
                      default=False, label='Use a known space group?')

        form.addParam('knownSpaceGroup', pwprot.StringParam,
                      condition='enterSpaceGroup', default='', label='Space group:')

        form.addParam('enterUnitCell', pwprot.BooleanParam,
                      default=False, label='Use a known unit cell?')

        form.addParam('knownUnitCell', pwprot.StringParam,
                      condition='enterUnitCell', default='', label='Unit cell:')

        form.addParam('indexMmSearchScope', pwprot.FloatParam, default=4.0,
                      help="Global radius of origin offset search.",
                      label='mm search scope',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('indexWideSearchBinning', pwprot.FloatParam, default=2,
                      help="Modify the coarseness of the wide grid search for the beam centre.",
                      label='Wide search binning',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('indexMinCellVolume', pwprot.FloatParam, default=25,
                      help="Minimum unit cell volume (in Angstrom^3).",
                      label='Min cell volume',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('indexMinCell', pwprot.FloatParam, default=3,
                      help="Minimum length of candidate unit cell basis vectors (in Angstrom).",
                      label='Min_cell',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('indexMaxCell', pwprot.FloatParam, default=None,
                      label='Max_cell', allowsNull=True,
                      help="Maximum length of candidate unit cell basis vectors (in Angstrom).",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('misindexCheckGridScope', pwprot.IntParam, default=0,
                      help="Search scope for testing misindexing on h, k, l.",
                      label='Misindexing check grid scope',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('doFilter_ice', pwprot.BooleanParam, default=False,
                      label='Filter ice?', expertLevel=pwprot.LEVEL_ADVANCED,
                      help="Filter out reflections at typical ice ring resolutions before max_cell estimation.")

        form.addSection('Refinement parameter configuration')

        form.addParam('refineNproc', pwprot.IntParam,
                      default=1, label='nproc',
                      help='The number of processes to use. Not all choices of refinement engine support nproc > 1.'
                      'Where multiprocessing is possible, it is helpful only in certain circumstances,'
                      'so this is not recommended for typical use.',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        group = form.addGroup('Parametrisation')

        group.addParam('beamFixInSpindlePlane', pwprot.BooleanParam,
                       label='Fix beam in spindle plane?', default=True,
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted to a"
                       "direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow refinement of"
                       "the unit cell volume.",
                       )

        group.addParam('beamFixOutSpindlePlane', pwprot.BooleanParam,
                       label='Fix beam out of spindle plane?', default=False,
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted to"
                       "a direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow refinement"
                       "of the unit cell volume.",
                       )

        group.addParam('beamFixWavelength', pwprot.BooleanParam,
                       label='Fix beam wavelength?', default=True,
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted"
                       "to a direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow"
                       "refinement of the unit cell volume.",
                       )

        group.addParam('crystalFixCell', pwprot.BooleanParam,
                       label='Crystal: Fix cell?', default=True,
                       help="Fix crystal parameters",
                       )

        group.addParam('crystalFixOrientation', pwprot.BooleanParam,
                       label='Crystal: Fix orientation?', default=True,
                       help="Fix crystal parameters",
                       )

        group.addParam('detectorFixPosition', pwprot.BooleanParam,
                       label='Fix detector position?', default=True,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.",
                       )

        group.addParam('detectorFixOrientation', pwprot.BooleanParam,
                       label='Fix detector orientation?', default=True,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.",
                       )

        group.addParam('detectorFixdistance', pwprot.BooleanParam,
                       label='Fix detector distance?', default=True,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.",
                       )

        group = form.addGroup('Refinery')

        group.addParam('doSetMaxIterations', pwprot.BooleanParam,
                       label='Do you want to set the maximum number of iterations?', default=False,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       help="Maximum number of iterations in refinement before termination.",
                       )

        group.addParam('refineryMaxIterations', pwprot.IntParam,
                       default=None,
                       allowsNull=True, help="Maximum number of iterations in refinement before termination."
                       "None implies the engine supplies its own default.",
                       label='Max iterations', condition="doSetMaxIterations",
                       )

        # Allow some options if the Bravais settings are to be refined
        form.addSection(label='Refine Bravais settings',
                        condition='doRefineBravaisSettings')

        form.addParam('refineBravNproc', pwprot.IntParam,
                      default=4, label="How many processors do you want to use?",
                      help="The number of processes to use.")

        # Allow some options that are only relevant for reindexing
        form.addSection(label='Reindex',
                        condition='doReindex')

        form.addParam('doReindexModel', pwprot.BooleanParam,
                      default=False, label="Do you want to reindex the experimental model?")

        form.addParam('doReindexReflections', pwprot.BooleanParam,
                      default=False, label="Do you want to reindex the experimental reflections?")

        # Allow adding anything else with command line syntax
        form.addSection('Plain command line input')
        form.addParam('commandLineInputIndexing', pwprot.StringParam,
                      expertLevel=pwprot.LEVEL_ADVANCED,
                      default='', condition="doIndex==True",
                      help="Anything added here will be added at the end of the command line")
        form.addParam('commandLineInputBravais', pwprot.StringParam,
                      expertLevel=pwprot.LEVEL_ADVANCED,
                      default='', condition="doReindexReflections",
                      help="Anything added here will be added at the end of the command line")
        form.addParam('commandLineInputReindexing', pwprot.StringParam,
                      expertLevel=pwprot.LEVEL_ADVANCED,
                      default='', condition="doReindexModel==True",
                      help="Anything added here will be added at the end of the command line")

   # -------------------------- INSERT functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(
            'convertInputStep', self.inputImages.getObjId(), self.inputSpots.getObjId())
        if self.doIndex:
            self._insertFunctionStep('indexStep')
        if self.doRefineBravaisSettings:
            self._insertFunctionStep('refineBravaisStep')
        if self.doReindex:
            self._insertFunctionStep('reindexStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, inputImgId, inputSpotId):
        inputImages = self.inputImages.get()
        inputSpots = self.inputSpots.get()
        self.info("Number of images: %s" % inputImages.getSize())
        self.info("Number of spots: %s" % inputSpots.getSpots())
        # Write new model and/or reflection file if no was supplied from set
        if self._checkWriteModel():
            self.writeJson(inputImages, self.getInputModelFile())
        if self._checkWriteRefl():
            self.writeRefl(inputSpots, self.getInputReflFile())

    def indexStep(self):
        program = 'dials.index'
        arguments = self._prepIndexCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    def refineBravaisStep(self):
        program = 'dials.refine_bravais_settings'
        arguments = self._prepBravaisCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    def reindexStep(self):
        program = 'dials.reindex'
        arguments = self._prepReindexCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    # TODO: Create a temporary "SetOfIndexedSpotsFile" that only saves the file location
    def createOutputStep(self):
        # Check that the indexing created proper output
        assert(os.path.exists(self.getOutputReflFile()))
        assert(os.path.exists(self.getOutputModelFile()))

        reflectionData = readRefl(self.getOutputReflFile())
        outputSet = self._createSetOfIndexedSpots()
        iSpot = IndexedSpot()
        numberOfSpots = reflectionData[2]
        reflDict = reflectionData[4]

        outputSet.setSpots(numberOfSpots)
        outputSet.setDialsModel(self.getOutputModelFile())
        outputSet.setDialsRefl(self.getOutputReflFile())

        for i in range(0, numberOfSpots):
            iSpot.setObjId(i+1)
            iSpot.setSpotId(reflDict['id'][i])
            iSpot.setBbox(reflDict['bbox'][i])
            iSpot.setFlag(reflDict['flags'][i])
            iSpot.setIntensitySumValue(reflDict['intensity.sum.value'][i])
            iSpot.setIntensitySumVariance(
                reflDict['intensity.sum.variance'][i])
            iSpot.setNSignal(reflDict['n_signal'][i])
            iSpot.setPanel(reflDict['panel'][i])
            try:
                iSpot.setShoebox(reflDict['shoebox'][i])
            except IndexError:
                pass
            iSpot.setXyzobsPxValue(reflDict['xyzobs.px.value'][i])
            iSpot.setXyzobsPxVariance(reflDict['xyzobs.px.variance'][i])
            outputSet.append(iSpot)

        outputSet.write()

        self._defineOutputs(outputIndexedSpots=outputSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        return errors

    # -------------------------- UTILS functions ------------------------------
    def getInputModelFile(self):
        if self.getSetModel():
            return self.getSetModel()
        else:
            return self._getExtraPath('imported.expt')

    def getInputReflFile(self):
        if self.getSetRefl():
            return self.getSetRefl()
        else:
            return self._getExtraPath('strong.refl')

    def getOutputModelFile(self):
        return self._getExtraPath('indexed_model.expt')

    def getOutputReflFile(self):
        return self._getExtraPath('indexed_reflections.refl')

    def getPhilPath(self):
        return self._getTmpPath('index.phil')

    def existsPath(self, path):
        return os.path.exists(path)

    def getSetModel(self):
        inputImages = self.inputImages.get()
        inputSpots = self.inputSpots.get()
        if self.existsPath(inputSpots.getDialsModel()):
            return inputSpots.getDialsModel()
        elif self.existsPath(inputImages.getDialsModel()):
            return inputImages.getDialsModel()
        else:
            return None

    def getSetRefl(self):
        inputImages = self.inputImages.get()
        inputSpots = self.inputSpots.get()
        if self.existsPath(inputSpots.getDialsRefl()):
            return inputSpots.getDialsRefl()
        elif self.existsPath(inputImages.getDialsRefl()):
            return inputImages.getDialsRefl()
        else:
            return None

    def getChangeOfBasisOp(self):
        # TODO: extract change of basis op from result in refineBravaisStep
        change_of_basis_op = 'a,b,c'
        return change_of_basis_op

    def _checkWriteModel(self):
        return self.getSetModel() != self.getInputModelFile()

    def _checkWriteRefl(self):
        return self.getSetRefl() != self.getInputReflFile()

    def _prepIndexCommandline(self, program):
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

        if self.indexNproc.get() not in (None, 1):
            params += " indexing.nproc={}".format(self.indexNproc.get())

        if self.enterSpaceGroup.get():
            params += " indexing.known_symmetry.space_group={}".format(
                self.knownSpaceGroup.get())

        if self.enterUnitCell.get():
            params += " indexing.known_symmetry.unit_cell={}".format(
                self.knownUnitCell.get())

        if self.indexMmSearchScope.get() not in (None, 4.0):
            params += " index.mm_search_scope={}".format(
                self.indexMmSearchScope.get())

        if self.indexWideSearchBinning.get() not in (None, 2):
            params += " index.wide_search_binning={}".format(
                self.indexWideSearchBinning.get())

        if self.indexMinCellVolume.get() not in (None, 25):
            params += " index.min_cell_volume={}".format(
                self.indexMinCellVolume.get())

        if self.indexMinCell.get() not in (None, 3.0):
            params += " index.min_cell={}".format(self.indexMinCell.get())

        if self.indexMaxCell.get() is not None:
            params += " index.max_cell={}".format(self.indexMaxCell.get())

        if self.misindexCheckGridScope.get() not in (None, 0):
            params += " check_misindexing.grid_search_scope={}".format(
                self.misindexCheckGridScope.get())

        if self.doFilter_ice.get():
            params += " indexing.max_cell_estimation.filter_ice={}".format(
                self.doFilter_ice.get())

        if self.refineNproc.get() not in (None, 1):
            params += " refinement.nproc={}".format(self.refineNproc.get())

        """ # FIXME: make all beamfix one block
        if True in (self.beamFixInSpindlePlane.get(), self.beamFixOutSpindlePlane.get(), self.beamFixWavelength.get()):
            beamfix = []
            if self.beamFixInSpindlePlane.get() is True:
                beamfix += ["in_spindle_plane"]
            if self.beamFixOutSpindlePlane.get() is True:
                beamfix += ["out_spindle_plane"]
            if self.beamFixWavelength.get() is True:
                beamfix += ["wavelength"]
            params += " refinement.parameterisation.beam.fix={}".format(
                " ".join(beamfix)) """

        # FIXME: Combine in one line
        if self.crystalFixCell.get() not in (None, True):
            params += " refinement.parameterisation.crystal.fix.cell={}".format(
                self.crystalFixCell.get())

        if self.crystalFixOrientation.get() not in (None, True):
            params += " refinement.parameterisation.crystal.fix.orientation={}".format(
                self.crystalFixOrientation.get())

        # FIXME: Convert to one line
        if self.detectorFixPosition.get() not in (None, False):
            params += " detector.fix=position"

        if self.detectorFixOrientation.get() not in (None, False):
            params += " refinement.parameterisation.detector.fix=orientation"

        if self.detectorFixdistance.get() not in (None, False):
            params += " refinement.parameterisation.detector.fix=distance"

        if self.refineryMaxIterations.get() is not None:
            params += " refinery.max_iterations={}".format(
                self.refineryMaxIterations.get())

        if self.commandLineInputIndexing.get():
            params += " {}".format(self.commandLineInputIndexing.get())

        return params

    def _prepBravaisCommandline(self, program):
        "Create the command line input to run dials programs"
        # Input basic parameters
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        params = "{} {} output.log={}".format(
            self.getOutputModelFile(), self.getOutputReflFile(), logPath)

        # Update the command line with additional parameters

        if self.refineBravNproc.get() not in (None, 4):
            params += " nproc={}".format(self.refineBravNproc.get())

        if self.commandLineInputBravais.get():
            params += " {}".format(self.commandLineInputBravais.get())

        return params

    def _prepReindexCommandline(self, program):
        "Create the command line input to run dials programs"
        # Input basic parameters
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        params = "change_of_basis_op={}".format(
            self.getChangeOfBasisOp(), logPath)

        if self.doReindexModel.get():
            params += " {}".format(self.getOutputModelFile())

        if self.doReindexReflections.get():
            params += " {}".format(self.getOutputReflFile())

        if self.commandLineInputReindexing.get():
            params += " {}".format(self.commandLineInputReindexing.get())

        return params

    def _createScanRanges(self):
        # Go through the
        images = [image.getObjId() for image in self.inputImages.get()
                  if image.getIgnore() is not True]
        scanranges = find_subranges(images)
        scanrange = ' '.join('spotfinder.scan_range={},{}'.format(i, j)
                             for i, j in scanranges)
        return scanrange
