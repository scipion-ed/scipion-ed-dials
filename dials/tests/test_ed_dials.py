# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              V E.G. Bengtsson       (viktor.bengtsson@mmk.su.se)   [2]
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
from dials.protocols.protocol_merge import DialsProtMerge

import pyworkflow as pw
import pyworkflow.tests as pwtests

import pwed
from pwed.objects import *
from pwed.protocols import ProtImportDiffractionImages

from dials.protocols import *
from dials.convert import writeRestraintsPhil
from dials.constants import *


pw.Config.setDomain(pwed)
if not pw.Config.debugOn():
    pw.Config.toggleDebug()

# Create toggles for skipping some tests
SKIP_PIPELINES = False
SKIP_LYSO = False
SKIP_GARNET = False
if SKIP_LYSO and SKIP_GARNET:
    SKIP_PIPELINES = True
SKIP_UTILS = False
KEEP_PROTOCOL_TEST_OUTPUT = True
KEEP_UTILS_TEST_OUTPUT = False


class TestEdDialsProtocols(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        if SKIP_PIPELINES:
            cls.skipTest(cls, "Skipping pipelines")
        pwtests.setupTestProject(cls, writeLocalConfig=True)
        pwtests.setupTestOutput(cls)
        cls.dataPath = os.environ.get("SCIPION_ED_TESTDATA")
        cls.PROJECT_NAME = cls.__name__

        if not os.path.exists(cls.dataPath):
            raise Exception(
                f"Can not run DIALS tests, missing file:\n  {cls.dataPath}")

    @classmethod
    def tearDownClass(cls):
        if not KEEP_PROTOCOL_TEST_OUTPUT:
            # Clean up all output files from the test
            pw.utils.cleanPath(cls.getOutputPath())

    # Functions for running protocols

    def _runImportImages(self,
                         filesPath,
                         filesPattern='SMV/data/00{TI}.img',
                         **kwargs):
        protImport = self.newProtocol(
            ProtImportDiffractionImages,
            filesPath=filesPath,
            filesPattern=filesPattern,
            **kwargs)
        self.launchProtocol(protImport)
        return protImport

    def _runDialsImportImages(self,
                              filesPath,
                              filesPattern='SMV/data/00{TI}.img',
                              **kwargs):
        protImport = self.newProtocol(
            DialsProtImportDiffractionImages,
            filesPath=filesPath,
            filesPattern=filesPattern,
            **kwargs)
        self.launchProtocol(protImport)
        return protImport

    def _runFindSpots(self, inputImages, **kwargs):
        protFindSpots = self.newProtocol(DialsProtFindSpots,
                                         inputImages=inputImages,
                                         **kwargs)
        self.launchProtocol(protFindSpots)
        return protFindSpots

    def _runIndex(self, inputImages, inputSpots, **kwargs):
        protIndex = self.newProtocol(DialsProtIndexSpots,
                                     inputImages=inputImages,
                                     inputSpots=inputSpots,
                                     **kwargs)
        self.launchProtocol(protIndex)
        return protIndex

    def _runRefine(self, inputSet, **kwargs):
        protRefine = self.newProtocol(DialsProtRefineSpots,
                                      inputSet=inputSet,
                                      **kwargs)
        self.launchProtocol(protRefine)
        return protRefine

    def _runIntegrate(self, inputSet, **kwargs):
        protIntegrate = self.newProtocol(DialsProtIntegrateSpots,
                                         inputSet=inputSet,
                                         **kwargs)
        self.launchProtocol(protIntegrate)
        return protIntegrate

    def _runSymmetry(self, inputSet, **kwargs):
        protSymmetry = self.newProtocol(DialsProtSymmetry,
                                        inputSet=inputSet,
                                        **kwargs)
        self.launchProtocol(protSymmetry)
        return protSymmetry

    def _runScaling(self, inputSets, **kwargs):
        protScaling = self.newProtocol(DialsProtScaling,
                                       inputSets=inputSets,
                                       **kwargs)
        self.launchProtocol(protScaling)
        return protScaling

    def _runMerging(self, inputSet, **kwargs):
        protMerging = self.newProtocol(DialsProtMerge,
                                       inputSet=inputSet,
                                       **kwargs)
        self.launchProtocol(protMerging)
        return protMerging

    def _runExport(self, inputSet, **kwargs):
        protExport = self.newProtocol(DialsProtExport,
                                      inputSet=inputSet,
                                      **kwargs)
        self.launchProtocol(protExport)
        return protExport

    # Helper functions
    def makePathList(self, scanRange, digits=3,
                     location=None,
                     pattern='SMV/data/00{TI}.img',
                     wildcard='{TI}'):
        first, last = scanRange.split(',')
        numbers = range(int(first), int(last)+1)
        replacements = [str(num).zfill(digits) for num in numbers]
        images = [pattern.replace(wildcard, replacement)
                  for replacement in replacements]
        imageList = [os.path.join(location, i) for i in images]
        return imageList

    def getLysoTestExperiments(self):
        experiments = []
        lyso_experiment_14 = {
            "location": "lysozyme/experiment_14",
            "files_pattern": "SMV/data/00{TI}.img",
            "d_min": 2.5,
            "found_spots": 1322,
            "distance": 1480.56,
            "rotation_axis": "-0.6204,-0.7843,0.0",
            "lyso": True,
            "sample": "lyso",
            "scan_range": "1,49",
            "space_group": "P 4 2 2",
            "unit_cell": "77.76,77.76,40.5,90,90,90",
            "unit_cell_sigmas": "0.05,0.05,0.05,0.05,0.05,0.05",
            "cb_op": "b,c,a",
        }
        lyso_experiment_24 = {
            "location": "lysozyme/experiment_24",
            "files_pattern": "SMV/data/00{TI}.img",
            "d_min": 3.5,
            "found_spots": 847,
            "distance": None,
            "rotation_axis": "-0.6204,-0.7843,0.0",
            "lyso": True,
            "sample": "lysozyme",
            "scan_range": "1,46",
            "space_group": "P 4 2 2",
            "unit_cell": "78.24,78.24,39.75,90,90,90",
            "unit_cell_sigmas": "0.05,0.05,0.05,0.05,0.05,0.05",
            "cb_op": "b,c,a",
        }
        experiments.append(lyso_experiment_14)
        experiments.append(lyso_experiment_24)
        return experiments

    def getGarnetExperiment(self):
        experiment = {
            "location": "garnet/experiment_1",
            "files_pattern": "Garnet_Cont_{TI}.img",
            "sample": "garnet",
            "replace_rotation_axis": True,
            "rotation_axis": "-0.755,-0.656,0.0",
            "useTemplate": True,
            "tsReplacement": "0###",
            "min_spot_size": 10,
            "max_spot_size": 1500000,
            "max_separation": 7.5,
            "d_min": 0.5,
            "d_max": 3.0,
            "sigma_background": 1.0,
            "found_spots": 6223,
            "nproc": 8,
            "bravais_setting_reference": "bravais_setting_22.expt",
            "cb_op": "b+c,a+c,a+b",
        }
        return experiment

    def checkLogDataset(self, protocol, dataset, logOutput=''):
        datasetString = (f"Source of data:\n"
                         f"{os.path.join(dataset, 'SMV/data')}")
        self.assertEqual(protocol.getDatasets(), datasetString)
        outputCompare = protocol.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(), logOutput.strip())

    def _getCommand(self, protocol, program=None):
        try:
            return protocol._prepareCommandline(program)
        except AttributeError:
            pass
        try:
            return protocol._prepCommandline(program)
        except AttributeError:
            pass
        try:
            return protocol._prepareCommandLineArguments(program)
        except AttributeError:
            pass

    def assertCommand(self, protocol, commandString, program=None):
        CL = self._getCommand(protocol, program)
        self.assertEqual(CL, commandString)

    def assertFileExists(self, file):
        self.assertTrue(os.path.exists(file))

    def comparePhils(self, goodPhil='restraints.phil',
                     testPhil=None):
        self.assertIsNotNone(goodPhil)
        self.assertIsNotNone(testPhil)
        with open(goodPhil, 'r') as f1:
            with open(testPhil, 'r') as f2:
                self.assertEqual(f1.read(), f2.read())

    # Pipelines
    def test_lyso_pipeline(self):
        if SKIP_LYSO:
            self.skipTest("Skipping lyso pipeline test")

        scaledSets = []
        scaleProt = []
        multiDataset = []

        # Define parameter reference strings
        beamParams = ("refinement.parameterisation.beam.fix='all "
                      "*in_spindle_plane out_spindle_plane *wavelength'")
        beamParamsAll = ("refinement.parameterisation.beam.fix='*all "
                         "in_spindle_plane out_spindle_plane wavelength'")
        beamParamsSv = ("refinement.parameterisation.beam.fix='all "
                        "in_spindle_plane out_spindle_plane *wavelength'")
        crystalParams = ("refinement.parameterisation.crystal.fix='all "
                         "cell orientation'")
        detectorParamsAll = ("refinement.parameterisation.detector.fix='*all "
                             "position orientation distance'")
        detectorParamsNoFix = ("refinement.parameterisation.detector.fix='all "
                               "position orientation distance'")
        gonioParams = ("refinement.parameterisation.goniometer.fix='*all "
                       "in_beam_plane out_beam_plane'")

        # Start the run
        for experiment in self.getLysoTestExperiments():
            exptId = experiment['location']
            with self.subTest(msg=f"Pipeline using {exptId}", experiment=exptId):

                dataset = os.path.join(self.dataPath, exptId)
                multiDataset.append(dataset)

                # Run import
                protImport = self._runDialsImportImages(
                    filesPath=dataset,
                    filesPattern=experiment['files_pattern'],
                    objLabel=f"dials - import diffraction images\n{exptId}",
                    replaceRotationAxis=True,
                    rotationAxis=experiment['rotation_axis'],
                    overwriteDetectorDistance=experiment['distance'],
                )
                if experiment['distance'] != None:
                    distance = f" distance={experiment['distance']}"
                else:
                    distance = ''

            pathlist = self.makePathList(experiment['scan_range'],
                                         location=dataset,
                                         pattern=experiment['files_pattern'])
            importCL = (
                f"{' '.join(pathlist)} "
                f"output.log={protImport._getLogsPath()}/dials.import.log "
                f"output.experiments={protImport._getExtraPath()}/imported.expt "
                f"goniometer.axes={experiment['rotation_axis']}"
                f"{distance.rstrip()}")

            self.assertCommand(protImport, importCL,
                               program='dials.import')
            outputModel = protImport.getOutputModelFile()
            self.assertFileExists(outputModel)
            importedset = getattr(
                protImport, 'outputDiffractionImages', None)
            self.assertIsNotNone(importedset)
            setModel = importedset.getDialsModel()
            self.assertFileExists(setModel)
            self.checkLogDataset(
                protImport, dataset, '')

            # Run find spots
            protFindSpots = self._runFindSpots(
                protImport.outputDiffractionImages,
                dMin=experiment['d_min'],
                untrustedAreas=True,
                untrustedRectangle_1='0,516,255,261',
                untrustedRectangle_2='255,261,0,516',
                thresholdAlgorithm=DISPERSION_EXTENDED)
            inputModel = protFindSpots.getInputModelFile()
            self.assertFileExists(inputModel)
            findSpotCL = (
                f"{protImport._getExtraPath()}/imported.expt "
                f"output.log={protFindSpots._getLogsPath()}/dials.find_spots.log "
                f"output.reflections={protFindSpots._getExtraPath()}/strong.refl "
                f"spotfinder.scan_range={experiment['scan_range']} "
                f"spotfinder.filter.d_min={experiment['d_min']} "
                f"spotfinder.filter.max_spot_size=1000 "
                f"spotfinder.filter.max_strong_pixel_fraction=0.25 "
                f"spotfinder.filter.max_separation=2.0 "
                f"spotfinder.filter.untrusted.rectangle=0,516,255,261 "
                f"spotfinder.filter.untrusted.rectangle=255,261,0,516 "
                f"spotfinder.threshold.algorithm=dispersion_extended "
                f"spotfinder.threshold.dispersion.sigma_background=6.0 "
                f"spotfinder.threshold.dispersion.sigma_strong=3.0 "
                f"spotfinder.threshold.dispersion.kernel_size=3,3")
            self.assertCommand(protFindSpots, findSpotCL,
                               program='dials.find_spots')
            foundspotset = getattr(
                protFindSpots, 'outputDiffractionSpots', None)
            self.assertIsNotNone(foundspotset)
            self.assertEqual(foundspotset.getSpots(),
                             experiment['found_spots'])
            self.assertFileExists(foundspotset.getDialsModel())
            self.assertFileExists(foundspotset.getDialsRefl())
            self.checkLogDataset(protFindSpots, dataset,
                                 'Histogram of per-image spot count for imageset 0:')

            # Run indexing
            protIndex = self._runIndex(
                objLabel="dials - index and reindex",
                inputImages=protImport.outputDiffractionImages,
                inputSpots=protFindSpots.outputDiffractionSpots,
                doRefineBravaisSettings=True,
                doReindex=True,
                detectorFixPosition=True,
                detectorFixOrientation=True,
                detectorFixDistance=True,
                beamFixInSpindlePlane=True,
                beamFixOutSpindlePlane=True,
                beamFixWavelength=True,
                copyBeamFix=False,
                copyCrystalFix=False,
                copyDetectorFix=False,
                copyGonioFix=False,
            )
            indexTmp = protIndex._getTmpPath()
            indexLogs = protIndex._getLogsPath()
            indexExtra = protIndex._getExtraPath()

            indexCL = (
                f"{protImport._getExtraPath()}/imported.expt "
                f"{protFindSpots._getExtraPath()}/strong.refl "
                f"output.log={indexLogs}/dials.index.log "
                f"output.experiments={indexTmp}/indexed.expt "
                f"output.reflections={indexTmp}/indexed.refl "
                f"{beamParamsAll} {crystalParams} "
                f"{detectorParamsAll} {gonioParams}")
            refBravCL = (
                f"{indexTmp}/indexed.expt "
                f"{indexTmp}/indexed.refl "
                f"output.log={indexLogs}/dials.refine_bravais_settings.log "
                f"output.directory={indexTmp}")
            reindexCL = (
                f"change_of_basis_op={experiment['cb_op']} {indexTmp}/indexed.refl "
                f"output.reflections={indexTmp}/reindexed.refl")
            self.assertEqual(protIndex._prepIndexCommandline(
                'dials.index'), indexCL)
            self.assertEqual(protIndex._prepBravaisCommandline(
                'dials.refine_bravais_settings'), refBravCL)
            self.assertEqual(protIndex._prepReindexCommandline(), reindexCL)
            indexedset = getattr(protIndex, 'outputIndexedSpots', None)
            self.assertIsNotNone(protIndex.outputIndexedSpots)
            self.assertFileExists(indexedset.getDialsModel())
            self.assertFileExists(indexedset.getDialsRefl())
            self.checkLogDataset(protIndex, dataset)

            with self.subTest(msg='Testing with restraints in phil file'):
                spaceGroup = experiment['space_group'].replace(' ', '')
                protIndexPhil = self._runIndex(
                    objLabel="dials - index known space group",
                    inputImages=protImport.outputDiffractionImages,
                    inputSpots=protFindSpots.outputDiffractionSpots,
                    detectorFixPosition=True,
                    detectorFixOrientation=True,
                    detectorFixDistance=True,
                    beamFixInSpindlePlane=True,
                    beamFixOutSpindlePlane=True,
                    beamFixWavelength=True,
                    enterSpaceGroup=True,
                    knownSpaceGroup=spaceGroup,
                )
                indexTmpPhil = protIndexPhil._getTmpPath()
                indexLogsPhil = protIndexPhil._getLogsPath()
                indexExtraPhil = protIndexPhil._getExtraPath()
                indexCLPhil = (
                    f"{protImport._getExtraPath()}/imported.expt "
                    f"{protFindSpots._getExtraPath()}/strong.refl "
                    f"output.log={indexLogsPhil}/dials.index.log "
                    f"output.experiments={indexTmpPhil}/indexed.expt "
                    f"output.reflections={indexTmpPhil}/indexed.refl "
                    f"indexing.known_symmetry.space_group={spaceGroup} "
                    f"{beamParamsAll} {crystalParams} "
                    f"{detectorParamsAll} {gonioParams}")
                self.assertEqual(protIndexPhil._prepIndexCommandline(
                    'dials.index'), indexCLPhil)
                indexedSetPhil = getattr(
                    protIndexPhil, 'outputIndexedSpots', None)
                self.assertIsNotNone(protIndexPhil.outputIndexedSpots)
                self.assertFileExists(indexedSetPhil.getDialsModel())
                self.assertFileExists(indexedSetPhil.getDialsRefl())
                self.checkLogDataset(protIndexPhil, dataset)

                # Refinement with target restraints
                protRefinePhil = self._runRefine(
                    objLabel="dials - static refinement with restraints",
                    inputSet=protIndexPhil.outputIndexedSpots,
                    scanVaryingNew=STATIC,
                    detectorFixDistance=False,
                    useRestraint=True,
                    targetUnitCell=experiment['unit_cell'],
                    targetSigmas=experiment['unit_cell_sigmas'],
                )
                refineCLstaticPhil = (
                    f"{indexExtraPhil}/indexed.expt "
                    f"{indexExtraPhil}/indexed.refl "
                    f"output.log={protRefinePhil._getLogsPath()}/dials.refine.log "
                    f"output.experiments={protRefinePhil._getExtraPath()}/refined.expt "
                    f"output.reflections={protRefinePhil._getExtraPath()}/refined.refl "
                    f"scan_varying=False {beamParams} {crystalParams} "
                    f"{detectorParamsNoFix} {gonioParams} "
                    f"{protRefinePhil._getExtraPath()}/restraints.phil")
                self.assertCommand(protRefinePhil, refineCLstaticPhil,
                                   program='dials.refine')
                refinedsetPhil = getattr(
                    protRefinePhil, 'outputRefinedSpots', None)
                self.assertIsNotNone(protRefinePhil.outputRefinedSpots)
                self.assertFileExists(refinedsetPhil.getDialsModel())
                self.assertFileExists(refinedsetPhil.getDialsRefl())
                # Check that the correct phil file is actually created
                self.assertFileExists(protRefinePhil.getRestraintsPhil())
                referenceFile = os.path.abspath("restraints.phil")
                referencePhil = writeRestraintsPhil(
                    fn=referenceFile,
                    values=experiment["unit_cell"],
                    sigmas=experiment["unit_cell_sigmas"]
                )
                self.assertFileExists(referenceFile)
                self.comparePhils(goodPhil=referencePhil,
                                  testPhil=protRefinePhil.getRestraintsPhil())
                self.checkLogDataset(protRefinePhil, dataset)

            # Run refinement
            protRefine = self._runRefine(
                objLabel="dials - static refinement",
                inputSet=protIndex.outputIndexedSpots,
                scanVaryingNew=UNSET,
                detectorFixAll=True,
            )
            refineCLstatic = (
                f"{indexExtra}/indexed.expt "
                f"{indexExtra}/indexed.refl "
                f"output.log={protRefine._getLogsPath()}/dials.refine.log "
                f"output.experiments={protRefine._getExtraPath()}/refined.expt "
                f"output.reflections={protRefine._getExtraPath()}/refined.refl "
                f"{beamParams} {crystalParams} "
                f"{detectorParamsAll} {gonioParams}")
            self.assertCommand(protRefine, refineCLstatic,
                               program='dials.refine')
            refinedset = getattr(protRefine, 'outputRefinedSpots', None)
            self.assertIsNotNone(protRefine.outputRefinedSpots)
            self.assertFileExists(refinedset.getDialsModel())
            self.assertFileExists(refinedset.getDialsRefl())
            self.checkLogDataset(protRefine, dataset)

            # Run scan-varying refinement
            protSvRefine = self._runRefine(
                objLabel="dials - scan-varying refinement",
                inputSet=protRefine.outputRefinedSpots,
                scanVaryingNew=SCAN_VARYING,
                beamFixAll=False,
                beamFixInSpindlePlane=False,
                beamFixOutSpindlePlane=False,
                beamFixWavelength=True,
                beamForceStatic=False,
                detectorFixAll=True,
            )
            refineCLsv = (
                f"{protRefine._getExtraPath()}/refined.expt "
                f"{protRefine._getExtraPath()}/refined.refl "
                f"output.log={protSvRefine._getLogsPath()}/dials.refine.log "
                f"output.experiments={protSvRefine._getExtraPath()}/refined.expt "
                f"output.reflections={protSvRefine._getExtraPath()}/refined.refl "
                f"scan_varying=True {beamParamsSv} "
                f"beam.force_static=False {crystalParams} "
                f"{detectorParamsAll} {gonioParams}")
            self.assertCommand(protSvRefine, refineCLsv,
                               program='dials.refine')
            svRefinedset = getattr(
                protSvRefine, 'outputRefinedSpots', None)
            self.assertIsNotNone(protSvRefine.outputRefinedSpots)
            self.assertFileExists(svRefinedset.getDialsModel())
            self.assertFileExists(svRefinedset.getDialsRefl())
            self.checkLogDataset(protSvRefine, dataset)

            # Run scan-varying refinement based on old workflow
            protSvRefineOld = self._runRefine(
                objLabel="dials - scan-varying refinement",
                inputSet=protRefine.outputRefinedSpots,
                scanVarying=True,
                beamFixAll=False,
                beamFixInSpindlePlane=False,
                beamFixOutSpindlePlane=False,
                beamFixWavelength=True,
                beamForceStatic=False,
                detectorFixAll=True,
            )
            refineCLsvOld = (
                f"{protRefine._getExtraPath()}/refined.expt "
                f"{protRefine._getExtraPath()}/refined.refl "
                f"output.log={protSvRefineOld._getLogsPath()}/dials.refine.log "
                f"output.experiments={protSvRefineOld._getExtraPath()}/refined.expt "
                f"output.reflections={protSvRefineOld._getExtraPath()}/refined.refl "
                f"scan_varying=True {beamParamsSv} "
                f"beam.force_static=False {crystalParams} "
                f"{detectorParamsAll} {gonioParams}")
            self.assertCommand(protSvRefineOld, refineCLsvOld,
                               program='dials.refine')
            svRefinedsetOld = getattr(
                protSvRefine, 'outputRefinedSpots', None)
            self.assertIsNotNone(protSvRefineOld.outputRefinedSpots)
            self.assertFileExists(svRefinedsetOld.getDialsModel())
            self.assertFileExists(svRefinedsetOld.getDialsRefl())
            self.checkLogDataset(protSvRefineOld, dataset)

            # Run integration
            protIntegrate = self._runIntegrate(
                # objLabel=f"Dials integration: {exptId}",
                inputSet=protSvRefine.outputRefinedSpots,
                nproc=8,
                commandLineInput=f"prediction.d_min={experiment['d_min']}",
            )
            integrateCL = (
                f"{protSvRefine._getExtraPath()}/refined.expt "
                f"{protSvRefine._getExtraPath()}/refined.refl "
                f"output.log={protIntegrate._getLogsPath()}/dials.integrate.log "
                f"output.experiments={protIntegrate._getExtraPath()}/"
                f"integrated_model.expt "
                f"output.reflections={protIntegrate._getExtraPath()}/"
                f"integrated_reflections.refl "
                f"output.phil={protIntegrate._getExtraPath()}/dials.integrate.phil "
                f"nproc=8 prediction.d_min={experiment['d_min']}")
            self.assertCommand(protIntegrate, integrateCL,
                               program='dials.integrate')
            integratedset = getattr(
                protIntegrate, 'outputIntegratedSpots', None)
            self.assertIsNotNone(protIntegrate.outputIntegratedSpots)
            self.assertFileExists(integratedset.getDialsModel())
            self.assertFileExists(integratedset.getDialsRefl())
            self.checkLogDataset(protIntegrate, dataset,
                                 'Summary vs resolution')

            # Check symmetry and scale
            protSymmetry = self._runSymmetry(
                objLabel="dials - symmetry check",
                inputSet=protIntegrate.outputIntegratedSpots,
            )
            symmCL = (
                f"{protIntegrate._getExtraPath()}/integrated_model.expt "
                f"{protIntegrate._getExtraPath()}/integrated_reflections.refl "
                f"output.log={protSymmetry._getLogsPath()}/dials.symmetry.log "
                f"output.experiments={protSymmetry._getExtraPath()}/symmetrized.expt "
                f"output.reflections={protSymmetry._getExtraPath()}/symmetrized.refl "
                f"output.html={protSymmetry._getExtraPath()}/dials.symmetry.html "
                f"output.json={protSymmetry._getExtraPath()}/dials.symmetry.json"
            )

            self.assertCommand(protSymmetry, symmCL, 'dials.symmetry')
            symmetrizedset = getattr(
                protSymmetry, 'outputSymmetrizedSpots', None)
            self.assertIsNotNone(protSymmetry.outputSymmetrizedSpots)
            self.assertFileExists(symmetrizedset.getDialsModel())
            self.assertFileExists(symmetrizedset.getDialsRefl())
            self.checkLogDataset(
                protSymmetry, dataset,
                f"Recommended space group: {experiment['space_group']}")

            protScaling = self._runScaling(
                objLabel="dials - scaling",
                inputSets=[protIntegrate.outputIntegratedSpots],
            )
            scaleCL = (
                f"{protIntegrate._getExtraPath()}/integrated_model.expt "
                f"{protIntegrate._getExtraPath()}/integrated_reflections.refl "
                f"output.log={protScaling._getLogsPath()}/dials.scale.log "
                f"output.experiments={protScaling._getExtraPath()}/scaled.expt "
                f"output.reflections={protScaling._getExtraPath()}/scaled.refl "
                f"output.html={protScaling._getExtraPath()}/dials.scale.html "
                f"filtering.output.scale_and_filter_results="
                f"{protScaling._getExtraPath()}/scale_and_filter_results.json "
                f"cut_data.partiality_cutoff=0.4 cut_data.min_isigi=-5.0 "
                f"outlier_rejection=standard outlier_zmax=6.0 filtering.method=None "
                f"filtering.deltacchalf.max_cycles=6 "
                f"filtering.deltacchalf.max_percent_removed=10.0 "
                f"filtering.deltacchalf.mode=dataset "
                f"filtering.deltacchalf.group_size=10 "
                f"filtering.deltacchalf.stdcutoff=4.0")
            self.assertCommand(protScaling, scaleCL, 'dials.scale')
            scaledset = getattr(
                protScaling, 'outputScaledSpots', None)
            self.assertIsNotNone(protScaling.outputScaledSpots)
            self.assertFileExists(scaledset.getDialsModel())
            self.assertFileExists(scaledset.getDialsRefl())
            self.checkLogDataset(
                protScaling, dataset, 'Space group being used during scaling is P 4')

            scaledSets.append(protScaling.outputScaledSpots)
            scaleProt.append(protScaling)

        with self.subTest(msg='Scaling multiple datasets'):
            while len(scaledSets) < 2:
                self.skipTest('Not enough datasets to test scaling multiple')
            mergedFn = "merged_test.mtz"
            unmergedFn = "unmerged.mtz"
            crystalName = "Lysozyme"
            protMultiScaling = self._runScaling(
                objLabel="dials - scaling multiple",
                inputSets=scaledSets,
                checkConsistentIndexing=True,
                exportMergedMtz=True,
                mergedMtzName=mergedFn,
                exportUnmergedMtz=True,
                specifyExportPath=True,
                exportPath=self.getOutputPath(),
                crystalName=crystalName,
            )
            scaledExtra0 = scaleProt[0]._getExtraPath()
            scaledExtra1 = scaleProt[1]._getExtraPath()
            mergedMtzFile = f"{self.getOutputPath()}/{mergedFn}"
            unmergedMtzFile = f"{self.getOutputPath()}/{unmergedFn}"
            multiScaleCL = (
                f"{scaledExtra0}/scaled.expt "
                f"{scaledExtra0}/scaled.refl "
                f"{scaledExtra1}/scaled.expt "
                f"{scaledExtra1}/scaled.refl "
                f"output.log={protMultiScaling._getLogsPath()}/dials.scale.log "
                f"output.experiments={protMultiScaling._getExtraPath()}/scaled.expt "
                f"output.reflections={protMultiScaling._getExtraPath()}/scaled.refl "
                f"output.html={protMultiScaling._getExtraPath()}/dials.scale.html "
                f"output.merged_mtz={mergedMtzFile} "
                f"output.unmerged_mtz={unmergedMtzFile} "
                f"output.crystal_name={crystalName} "
                f"output.project_name={self.PROJECT_NAME} "
                f"filtering.output.scale_and_filter_results="
                f"{protMultiScaling._getExtraPath()}/scale_and_filter_results.json "
                f"cut_data.partiality_cutoff=0.4 cut_data.min_isigi=-5.0 "
                f"scaling_options.check_consistent_indexing=True "
                f"outlier_rejection=standard outlier_zmax=6.0 filtering.method=None "
                f"filtering.deltacchalf.max_cycles=6 "
                f"filtering.deltacchalf.max_percent_removed=10.0 "
                f"filtering.deltacchalf.mode=dataset "
                f"filtering.deltacchalf.group_size=10 "
                f"filtering.deltacchalf.stdcutoff=4.0")
            self.assertCommand(protMultiScaling, multiScaleCL, 'dials.scale')
            multiscaledset = getattr(
                protMultiScaling, 'outputScaledSpots', None)
            self.assertIsNotNone(protMultiScaling.outputScaledSpots)
            self.assertFileExists(mergedMtzFile)
            self.assertFileExists(unmergedMtzFile)
            self.assertFileExists(multiscaledset.getDialsModel())
            self.assertFileExists(multiscaledset.getDialsRefl())
            compareDatasets = 'Source of data:'
            for ds in multiDataset:
                compareDatasets += f"\n{os.path.join(self.dataPath, ds, 'SMV/data')}"
            self.assertEqual(protMultiScaling.getDatasets(), compareDatasets)

            # Test merging protocol
            mergeBins = 5
            protMerge = self._runMerging(
                inputSet=protMultiScaling.outputScaledSpots,
                xtalName=crystalName,
                nBins=mergeBins,
            )

            mergeCL = (
                f"{protMultiScaling._getExtraPath()}/scaled.expt "
                f"{protMultiScaling._getExtraPath()}/scaled.refl "
                f"output.log={protMerge._getLogsPath()}/dials.merge.log "
                f"output.html={protMerge._getExtraPath()}/dials.merge.html "
                f"output.mtz={protMerge._getExtraPath()}/merged.mtz "
                f"output.crystal_names={crystalName} "
                f"output.project_name={self.PROJECT_NAME}"
                f"assess_space_group=True "
                f"anomalous=True "
                f"truncate=True "
                f"wavelength_tolerance={1e-4} "
                f"combine_partials=True "
                f"partiality_threshold=0.4 "
                f"n_residues=200 "
                f"merging.use_internal_variance=False "
                f"merging.n_bins={mergeBins} "
                f"merging.anomalous=False"
            )

            self.assertCommand(protMerge, mergeCL, "dials.merge")
            mergedset = getattr(protMerge, "exportedFileSet", None)
            self.assertIsNotNone(protMerge.exportedFileSet)
            for exportFile in mergedset:
                self.assertFileExists(exportFile.getFilePath())
            self.assertEqual(protMerge.getDatasets(), compareDatasets)

            # Test export of output scaled together
            protExportMtz = self._runExport(
                inputSet=protMultiScaling.outputScaledSpots,
                exportFormat=MTZ,
            )

            exportMtzCL = (
                f"{protMultiScaling._getExtraPath()}/scaled.expt "
                f"{protMultiScaling._getExtraPath()}/scaled.refl "
                f"format=mtz mtz.hklout={protExportMtz._getExtraPath()}/"
                f"integrated_{protExportMtz.getObjId()}.mtz "
                f"output.log={protExportMtz._getLogsPath()}/dials.export.log "
                f"mtz.combine_partials=True mtz.partiality_threshold=0.99 "
                f"mtz.min_isigi=-5.0 mtz.crystal_name=XTAL "
                f"mtz.project_name={self.PROJECT_NAME}")
            self.assertCommand(protExportMtz, exportMtzCL, "dials.export")
            exportedmtzset = getattr(protExportMtz, 'exportedFileSet', None)
            self.assertIsNotNone(protExportMtz.exportedFileSet)
            for exportFile in exportedmtzset:
                self.assertFileExists(exportFile.getFilePath())
            self.assertEqual(protExportMtz.getDatasets(), compareDatasets)

        with self.subTest(msg='Excluding images'):
            exclusions = ["0:1:2", "0:10:12"]
            protScalingExclude = self._runScaling(
                objLabel="dials - scaling with exclusions",
                inputSets=[protIntegrate.outputIntegratedSpots],
                excludeImages=True,
                numberOfExclusions=2,
                imageGroup1=exclusions[0],
                imageGroup2=exclusions[1],
            )
            scaleExclusionCL = (
                f"{protIntegrate._getExtraPath()}/integrated_model.expt "
                f"{protIntegrate._getExtraPath()}/integrated_reflections.refl "
                f"output.log={protScalingExclude._getLogsPath()}/dials.scale.log "
                f"output.experiments={protScalingExclude._getExtraPath()}/scaled.expt "
                f"output.reflections={protScalingExclude._getExtraPath()}/scaled.refl "
                f"output.html={protScalingExclude._getExtraPath()}/dials.scale.html "
                f"filtering.output.scale_and_filter_results="
                f"{protScalingExclude._getExtraPath()}/scale_and_filter_results.json "
                f"cut_data.partiality_cutoff=0.4 cut_data.min_isigi=-5.0 "
                f"outlier_rejection=standard outlier_zmax=6.0 filtering.method=None "
                f"filtering.deltacchalf.max_cycles=6 "
                f"filtering.deltacchalf.max_percent_removed=10.0 "
                f"filtering.deltacchalf.mode=dataset "
                f"filtering.deltacchalf.group_size=10 "
                f"filtering.deltacchalf.stdcutoff=4.0 "
                f"exclude_images={exclusions[0]} "
                f"exclude_images={exclusions[1]}"
            )
            self.assertCommand(protScalingExclude,
                               scaleExclusionCL, 'dials.scale')
            scaledset = getattr(
                protScalingExclude, 'outputScaledSpots', None)
            self.assertIsNotNone(protScalingExclude.outputScaledSpots)
            self.assertFileExists(scaledset.getDialsModel())
            self.assertFileExists(scaledset.getDialsRefl())
            self.checkLogDataset(
                protScalingExclude, dataset, 'Space group being used during scaling is P 4')

    def test_garnet_pipeline(self):
        if SKIP_GARNET:
            self.skipTest("Skipping garnet pipeline test")

        # Define all experiment variables in one place
        experiment = self.getGarnetExperiment()
        exptId = experiment['location']
        dataset = os.path.join(self.dataPath, exptId)

        # Define parameter strings
        beamParams = ("refinement.parameterisation.beam.fix="
                      "'all *in_spindle_plane out_spindle_plane *wavelength'")
        crystalParams = ("refinement.parameterisation.crystal.fix="
                         "'all cell orientation'")
        detectorParamsDistance = ("refinement.parameterisation.detector.fix="
                                  "'all position orientation *distance'")
        gonioParams = ("refinement.parameterisation.goniometer.fix="
                       "'*all in_beam_plane out_beam_plane'")

        # Run import
        protImport = self._runDialsImportImages(
            filesPath=dataset,
            filesPattern=experiment['files_pattern'],
            objLabel=f"dials - import diffraction images\n{exptId}",
            replaceRotationAxis=experiment['replace_rotation_axis'],
            rotationAxis=experiment['rotation_axis'],
            useTemplate=experiment['useTemplate'],
            tsReplacement=experiment['tsReplacement'],
        )
        importCL = (
            f"template={dataset}/Garnet_Cont_{experiment['tsReplacement']}.img "
            f"output.log={protImport._getLogsPath()}/dials.import.log "
            f"output.experiments={protImport._getExtraPath()}/imported.expt "
            f"goniometer.axes={experiment['rotation_axis']}")

        self.assertCommand(protImport, importCL,
                           program='dials.import')
        outputModel = protImport.getOutputModelFile()
        self.assertFileExists(outputModel)
        importedset = getattr(
            protImport, 'outputDiffractionImages', None)
        self.assertIsNotNone(importedset)
        setModel = importedset.getDialsModel()
        self.assertFileExists(setModel)
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protImport.getDatasets(), datasetString)
        outputCompare = protImport.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(), ''.strip())

        # Run find spots
        protFindSpots = self._runFindSpots(
            protImport.outputDiffractionImages,
            minSpotSize=experiment['min_spot_size'],
            maxSpotSize=experiment['max_spot_size'],
            maxSeparation=experiment['max_separation'],
            dMin=experiment['d_min'],
            dMax=experiment['d_max'],
            sigmaBackground=experiment['sigma_background'],
            thresholdAlgorithm=DISPERSION,
        )
        inputModel = protFindSpots.getInputModelFile()
        self.assertFileExists(inputModel)
        findSpotCL = (f"{protImport._getExtraPath()}/imported.expt "
                      f"output.log={protFindSpots._getLogsPath()}/dials.find_spots.log "
                      f"output.reflections={protFindSpots._getExtraPath()}/strong.refl "
                      f"spotfinder.scan_range=1,566 "
                      f"spotfinder.filter.d_min={experiment['d_min']} "
                      f"spotfinder.filter.d_max={experiment['d_max']} "
                      f"spotfinder.filter.min_spot_size={experiment['min_spot_size']} "
                      f"spotfinder.filter.max_spot_size={experiment['max_spot_size']} "
                      f"spotfinder.filter.max_strong_pixel_fraction=0.25 "
                      f"spotfinder.filter.max_separation={experiment['max_separation']} "
                      f"spotfinder.threshold.algorithm=dispersion "
                      f"spotfinder.threshold.dispersion.sigma_background="
                      f"{experiment['sigma_background']} "
                      f"spotfinder.threshold.dispersion.sigma_strong=3.0 "
                      f"spotfinder.threshold.dispersion.kernel_size=3,3")
        self.assertCommand(protFindSpots, findSpotCL,
                           program="dials.find_spots")
        foundspotset = getattr(
            protFindSpots, 'outputDiffractionSpots', None)
        self.assertIsNotNone(foundspotset)
        self.assertEqual(foundspotset.getSpots(),
                         experiment['found_spots'])
        self.assertFileExists(foundspotset.getDialsModel())
        self.assertFileExists(foundspotset.getDialsRefl())
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protFindSpots.getDatasets(), datasetString)
        outputCompare = protFindSpots.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(),
                         'Histogram of per-image spot count for imageset 0:'.strip())
        with self.subTest(
                msg="Testing validation of resolution limits in Find Spots protocol"):
            with self.assertRaises(Exception):
                self._runFindSpots(
                    protImport.outputDiffractionImages,
                    dMin=experiment['d_max'],
                    dMax=experiment['d_min'],
                )

        # Run indexing
        protIndex = self._runIndex(
            objLabel="dials - index and reindex",
            inputImages=protImport.outputDiffractionImages,
            inputSpots=protFindSpots.outputDiffractionSpots,
            doRefineBravaisSettings=True,
            doReindex=True,
            copyBeamFix=False,
            copyCrystalFix=False,
            copyDetectorFix=False,
            copyGonioFix=False,
        )
        indexTmp = protIndex._getTmpPath()
        indexLogs = protIndex._getLogsPath()
        indexExtra = protIndex._getExtraPath()

        indexCL = (
            f"{protImport._getExtraPath()}/imported.expt "
            f"{protFindSpots._getExtraPath()}/strong.refl "
            f"output.log={indexLogs}/dials.index.log "
            f"output.experiments={indexTmp}/indexed.expt "
            f"output.reflections={indexTmp}/indexed.refl "
            f"{beamParams} {crystalParams} "
            f"{detectorParamsDistance} {gonioParams}")
        refBravCL = (
            f"{indexTmp}/indexed.expt {indexTmp}/indexed.refl "
            f"output.log={indexLogs}/dials.refine_bravais_settings.log "
            f"output.directory={indexTmp}")
        reindexCL = (
            f"change_of_basis_op={experiment['cb_op']} {indexTmp}/indexed.refl "
            f"output.reflections={indexTmp}/reindexed.refl")
        self.assertEqual(protIndex._prepIndexCommandline(
            'dials.index'), indexCL)
        self.assertEqual(protIndex._prepBravaisCommandline(
            'dials.refine_bravais_settings'), refBravCL)
        self.assertEqual(protIndex._prepReindexCommandline(), reindexCL)
        indexedset = getattr(protIndex, 'outputIndexedSpots', None)
        self.assertIsNotNone(protIndex.outputIndexedSpots)
        self.assertFileExists(indexedset.getDialsModel())
        self.assertFileExists(indexedset.getDialsRefl())
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protIndex.getDatasets(), datasetString)
        outputCompare = protIndex.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(), ''.strip())
        with self.subTest(
                msg="Testing refine_bravais_settings with copied parameters"):
            protIndexCopy = self._runIndex(
                objLabel="dials - index and refine bravais setting with copied parameters",
                inputImages=protImport.outputDiffractionImages,
                inputSpots=protFindSpots.outputDiffractionSpots,
                doRefineBravaisSettings=True,
                doReindex=False,
            )
            indexTmpCopy = protIndexCopy._getTmpPath()
            indexLogsCopy = protIndexCopy._getLogsPath()
            refBravCLCopy = (
                f"{indexTmpCopy}/indexed.expt {indexTmpCopy}/indexed.refl "
                f"output.log={indexLogsCopy}/dials.refine_bravais_settings.log "
                f"output.directory={indexTmpCopy} "
                f"{beamParams} {crystalParams} {detectorParamsDistance} {gonioParams}")
            self.assertEqual(protIndexCopy._prepBravaisCommandline(
                'dials.refine_bravais_settings'), refBravCLCopy)

        # Run refinement
        protRefine = self._runRefine(
            inputSet=protIndex.outputIndexedSpots,
            scanVaryingNew=UNSET,
        )
        refineCLstatic = (
            f"{indexExtra}/indexed.expt {indexExtra}/indexed.refl "
            f"output.log={protRefine._getLogsPath()}/dials.refine.log "
            f"output.experiments={protRefine._getExtraPath()}/refined.expt "
            f"output.reflections={protRefine._getExtraPath()}/refined.refl "
            f"{beamParams} {crystalParams} {detectorParamsDistance} {gonioParams}")
        self.assertCommand(protRefine, refineCLstatic,
                           program='dials.refine')
        refinedset = getattr(protRefine, 'outputRefinedSpots', None)
        self.assertIsNotNone(protRefine.outputRefinedSpots)
        self.assertFileExists(refinedset.getDialsModel())
        self.assertFileExists(refinedset.getDialsRefl())
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protRefine.getDatasets(), datasetString)
        outputCompare = protRefine.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(), ''.strip())

        # Run integration
        protIntegrate = self._runIntegrate(
            inputSet=protRefine.outputRefinedSpots,
            nproc=experiment['nproc'],
            dMin=experiment['d_min'],
            makeReport=True,
        )
        integrateCL = (
            f"{protRefine._getExtraPath()}/refined.expt "
            f"{protRefine._getExtraPath()}/refined.refl "
            f"output.log={protIntegrate._getLogsPath()}/dials.integrate.log "
            f"output.experiments={protIntegrate._getExtraPath()}/"
            f"integrated_model.expt "
            f"output.reflections={protIntegrate._getExtraPath()}/"
            f"integrated_reflections.refl "
            f"output.phil={protIntegrate._getExtraPath()}/dials.integrate.phil "
            f"nproc=8 prediction.d_min={experiment['d_min']}")
        self.assertCommand(protIntegrate, integrateCL,
                           program='dials.integrate')
        integratedset = getattr(
            protIntegrate, 'outputIntegratedSpots', None)
        self.assertIsNotNone(protIntegrate.outputIntegratedSpots)
        self.assertFileExists(integratedset.getDialsModel())
        self.assertFileExists(integratedset.getDialsRefl())
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protIntegrate.getDatasets(), datasetString)
        outputCompare = protIntegrate.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(),
                         'Summary vs resolution'.strip())
        with self.subTest(
                msg="Testing validation of resolution limits in Integrate protocol"):
            with self.assertRaises(Exception):
                self._runIntegrate(
                    inputSet=protRefine.outputRefinedSpots,
                    dMin=experiment['d_max'],
                    dMax=experiment['d_min'],
                )

        # Check symmetry and scale
        protSymmetry = self._runSymmetry(
            objLabel="dials - symmetry check",
            inputSet=protIntegrate.outputIntegratedSpots,
        )
        symmCL = (
            f"{protIntegrate._getExtraPath()}/integrated_model.expt "
            f"{protIntegrate._getExtraPath()}/integrated_reflections.refl "
            f"output.log={protSymmetry._getLogsPath()}/dials.symmetry.log "
            f"output.experiments={protSymmetry._getExtraPath()}/symmetrized.expt "
            f"output.reflections={protSymmetry._getExtraPath()}/symmetrized.refl "
            f"output.html={protSymmetry._getExtraPath()}/dials.symmetry.html "
            f"output.json={protSymmetry._getExtraPath()}/dials.symmetry.json"
        )
        self.assertCommand(protSymmetry, symmCL, 'dials.symmetry')
        symmetrizedset = getattr(
            protSymmetry, 'outputSymmetrizedSpots', None)
        self.assertIsNotNone(protSymmetry.outputSymmetrizedSpots)
        self.assertFileExists(symmetrizedset.getDialsModel())
        self.assertFileExists(symmetrizedset.getDialsRefl())
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protSymmetry.getDatasets(), datasetString)
        outputCompare = protSymmetry.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(),
                         'Recommended space group: I 4 3 2'.strip())

        protScaling = self._runScaling(
            objLabel="dials - scaling",
            inputSets=[protIntegrate.outputIntegratedSpots],
        )
        scaleCL = (
            f"{protIntegrate._getExtraPath()}/integrated_model.expt "
            f"{protIntegrate._getExtraPath()}/integrated_reflections.refl "
            f"output.log={protScaling._getLogsPath()}/dials.scale.log "
            f"output.experiments={protScaling._getExtraPath()}/scaled.expt "
            f"output.reflections={protScaling._getExtraPath()}/scaled.refl "
            f"output.html={protScaling._getExtraPath()}/dials.scale.html "
            f"filtering.output.scale_and_filter_results={protScaling._getExtraPath()}/"
            f"scale_and_filter_results.json cut_data.partiality_cutoff=0.4 "
            f"cut_data.min_isigi=-5.0 outlier_rejection=standard outlier_zmax=6.0 "
            f"filtering.method=None filtering.deltacchalf.max_cycles=6 "
            f"filtering.deltacchalf.max_percent_removed=10.0 "
            f"filtering.deltacchalf.mode=dataset filtering.deltacchalf.group_size=10 "
            f"filtering.deltacchalf.stdcutoff=4.0")
        self.assertCommand(protScaling, scaleCL, 'dials.scale')
        scaledset = getattr(
            protScaling, 'outputScaledSpots', None)
        self.assertIsNotNone(protScaling.outputScaledSpots)
        self.assertFileExists(scaledset.getDialsModel())
        self.assertFileExists(scaledset.getDialsRefl())
        self.assertEqual(protScaling.getDatasets(), datasetString)
        outputCompare = protScaling.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(),
                         'Space group being used during scaling is I 2 3'.strip())
        with self.subTest(msg="Testing validation of resolution limits in Scale protocol"):
            with self.assertRaises(Exception):
                self._runScaling(
                    inputSets=[protIntegrate.outputIntegratedSpots],
                    dMin=experiment['d_max'],
                    dMax=experiment['d_min'],
                )

        protExportMtz = self._runExport(
            inputSet=protScaling.outputScaledSpots,
            exportFormat=MTZ,
            mtzCrystalName="Garnet"
        )

        exportMtzCL = (
            f"{protScaling._getExtraPath()}/scaled.expt "
            f"{protScaling._getExtraPath()}/scaled.refl format=mtz "
            f"mtz.hklout={protExportMtz._getExtraPath()}/"
            f"integrated_{protExportMtz.getObjId()}.mtz "
            f"output.log={protExportMtz._getLogsPath()}/dials.export.log "
            f"mtz.combine_partials=True mtz.partiality_threshold=0.99 "
            f"mtz.min_isigi=-5.0 mtz.crystal_name=Garnet "
            f"mtz.project_name={self.PROJECT_NAME}")
        self.assertCommand(protExportMtz, exportMtzCL, "dials.export")
        exportedmtzset = getattr(protExportMtz, 'exportedFileSet', None)
        self.assertIsNotNone(protExportMtz.exportedFileSet)
        for exportFile in exportedmtzset:
            self.assertFileExists(exportFile.getFilePath())
        self.assertEqual(protExportMtz.getDatasets(), datasetString)


class TestEdDialsUtils(pwtests.BaseTest):
    @ classmethod
    def setUpClass(cls):
        if SKIP_UTILS:
            cls.skipTest(cls, "Skipping utils")
        pwtests.setupTestOutput(cls)
        cls.dataPath = os.environ.get("SCIPION_ED_TESTDATA")

        if not os.path.exists(cls.dataPath):
            raise Exception(
                f"Can not run utils tests, missing file:\n  {cls.dataPath}")

    @ classmethod
    def tearDownClass(cls):
        if not KEEP_UTILS_TEST_OUTPUT:
            # Clean up all output files from the test
            pw.utils.cleanPath(cls.getOutputPath())

    def comparePhils(self, goodPhil='restraints.phil', testPhil=None):
        self.assertIsNotNone(testPhil)
        with open(os.path.join(self.dataPath, "utils", goodPhil), 'r') as f1:
            with open(testPhil, 'r') as f2:
                self.assertEqual(f1.read(), f2.read())

    def test_write_restraints(self):
        # Input and output are as expected
        setFn = self.getOutputPath("restraints.phil")
        pw.utils.cleanPath(setFn)

        # Add some values as a start
        values = "10,20,30,90,90,90"
        sigmas = "0.05,0.05,0.05,0.05,0.05,0.05"

        outFn = writeRestraintsPhil(fn=setFn, values=values, sigmas=sigmas)
        self.comparePhils(goodPhil='restraints.phil', testPhil=outFn)

    def test_write_restraints_bad_input(self):
        # Check if the functions to fix bad input and typos work
        fnShort = self.getOutputPath("restraints_bad_input_short.phil")
        fnLong = self.getOutputPath("restraints_bad_input_long.phil")
        pw.utils.cleanPath(fnShort)
        pw.utils.cleanPath(fnLong)

        # Test different separators and do not add all values
        values = "10, 20 30,90 ,90"
        sigmas = "0.05,0.05"
        # Test cutting away extra values added to the end
        values2 = "10,20,30,90,90,90,90"
        sigmas2 = "0.05,0.05,0.05,0.05,0.05,0.05,0.05"

        writeRestraintsPhil(fn=fnShort, values=values, sigmas=sigmas)
        writeRestraintsPhil(fn=fnLong, values=values2, sigmas=sigmas2)

        self.comparePhils(goodPhil='restraints.phil', testPhil=fnShort)
        self.comparePhils(goodPhil='restraints.phil', testPhil=fnLong)

    def test_write_restraints_no_sigmas(self):
        setFn = self.getOutputPath("restraints_no_sigmas.phil")
        pw.utils.cleanPath(setFn)

        values = "10,20,30,90,90,90"
        writeRestraintsPhil(fn=setFn, values=values)
        self.comparePhils(goodPhil='restraints_no_sigmas.phil',
                          testPhil=setFn)
