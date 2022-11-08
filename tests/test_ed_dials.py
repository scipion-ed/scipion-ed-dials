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

import pwed
import pyworkflow as pw
import pyworkflow.tests as pwtests

import dials.constants as dconst
import tests.constants_cases as cc
import tests.testing_utils as tutils
from dials.convert import writeRestraintsPhil
from dials.objects import MissingPathException

pw.Config.setDomain(pwed)
if not pw.Config.debugOn():
    pw.Config.toggleDebug()


class TestEdDialsProtocols(tutils.ProtocolRunner, tutils.HelperCollection):
    @classmethod
    def setUpClass(cls):
        if cc.SKIP_PIPELINES:
            cls.skipTest(cls, "Skipping pipelines")
        pwtests.setupTestProject(cls, writeLocalConfig=True)
        pwtests.setupTestOutput(cls)
        cls.dataPath = os.environ.get("SCIPION_ED_TESTDATA")
        cls.PROJECT_NAME = cls.__name__

        if not os.path.exists(cls.dataPath):
            raise MissingPathException(
                f"Can not run DIALS tests, missing file:\n  {cls.dataPath}"
            )

    @classmethod
    def tearDownClass(cls):
        if not cc.KEEP_PROTOCOL_TEST_OUTPUT:
            # Clean up all output files from the test
            pw.utils.cleanPath(cls.getOutputPath())

    # Pipelines
    def test_lyso_pipeline(self):
        if cc.SKIP_LYSO:
            self.skipTest("Skipping lyso pipeline test")

        scaledSets = []
        scaleProt = []
        multiDataset = []

        # Start the run
        for experiment in [cc.lyso_experiment_14, cc.lyso_experiment_24]:
            exptId = experiment["location"]
            with self.subTest(
                msg=f"Pipeline using {exptId}", experiment=exptId
            ):

                dataset = os.path.join(self.dataPath, exptId)
                multiDataset.append(dataset)

                # Run import
                protImport = self._runDialsImportImages(
                    filesPath=dataset,
                    filesPattern=experiment["files_pattern"],
                    objLabel=f"dials - import diffraction images\n{exptId}",
                    replaceRotationAxis=True,
                    rotationAxis=experiment["rotation_axis"],
                    overwriteDetectorDistance=experiment["distance"],
                )
            self.assertCommand(
                protImport,
                cc.lysoImportCommandLine(
                    location=dataset,
                    extraPath=protImport._getExtraPath(),
                    logPath=protImport._getLogsPath(),
                    experiment=experiment,
                ),
                program="dials.import",
            )
            outputModel = protImport.getOutputModelFile()
            self.assertFileExists(outputModel)
            importedset = getattr(protImport, "outputDiffractionImages", None)
            self.assertIsNotNone(importedset)
            setModel = importedset.getDialsModel()
            self.assertFileExists(setModel)
            self.checkLogDataset(protImport, dataset, "")

            # Run find spots
            protFindSpots = self._runFindSpots(
                protImport.outputDiffractionImages,
                dMin=experiment["d_min"],
                untrustedAreas=True,
                untrustedRectangle_1="0,516,255,261",
                untrustedRectangle_2="255,261,0,516",
                thresholdAlgorithm=dconst.DISPERSION_EXTENDED,
            )
            inputModel = protFindSpots.getInputModelFile()
            self.assertFileExists(inputModel)
            self.assertCommand(
                protFindSpots,
                cc.lysoStandardFindSpotCommandLine(
                    inputExtraPath=protImport._getExtraPath(),
                    outputExtraPath=protFindSpots._getExtraPath(),
                    logPath=protFindSpots._getLogsPath(),
                    experiment=experiment,
                ),
                program="dials.find_spots",
            )
            foundspotset = getattr(
                protFindSpots, "outputDiffractionSpots", None
            )
            self.assertIsNotNone(foundspotset)
            self.assertEqual(
                foundspotset.getSpots(), experiment["found_spots"]
            )
            self.assertFileExists(foundspotset.getDialsModel())
            self.assertFileExists(foundspotset.getDialsRefl())
            self.checkLogDataset(
                protFindSpots,
                dataset,
                "Histogram of per-image spot count for imageset 0:",
            )

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
            self.assertEqual(
                protIndex._prepIndexCommandline("dials.index"),
                cc.lysoIndexCommandLine(
                    modelPath=protImport._getExtraPath(),
                    spotsPath=protFindSpots._getExtraPath(),
                    logPath=indexLogs,
                    indexTmp=indexTmp,
                ),
            )
            self.assertEqual(
                protIndex._prepBravaisCommandline(
                    "dials.refine_bravais_settings"
                ),
                cc.refineBravaisLatticeCommandLine(
                    indexTmp=indexTmp, indexLogs=indexLogs
                ),
            )
            self.assertEqual(
                protIndex._prepReindexCommandline(),
                cc.reindexCommandLine(
                    experiment=experiment, indexTmp=indexTmp
                ),
            )
            indexedset = getattr(protIndex, "outputIndexedSpots", None)
            self.assertIsNotNone(protIndex.outputIndexedSpots)
            self.assertFileExists(indexedset.getDialsModel())
            self.assertFileExists(indexedset.getDialsRefl())
            self.checkLogDataset(protIndex, dataset)

            with self.subTest(msg="Testing with restraints in phil file"):
                spaceGroup = experiment["space_group"].replace(" ", "")
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
                self.assertEqual(
                    protIndexPhil._prepIndexCommandline("dials.index"),
                    cc.lysoIndexCommandLine(
                        modelPath=protImport._getExtraPath(),
                        spotsPath=protFindSpots._getExtraPath(),
                        logPath=protIndexPhil._getLogsPath(),
                        indexTmp=protIndexPhil._getTmpPath(),
                        spaceGroup=spaceGroup,
                    ),
                )
                indexedSetPhil = getattr(
                    protIndexPhil, "outputIndexedSpots", None
                )
                self.assertIsNotNone(protIndexPhil.outputIndexedSpots)
                self.assertFileExists(indexedSetPhil.getDialsModel())
                self.assertFileExists(indexedSetPhil.getDialsRefl())
                self.checkLogDataset(protIndexPhil, dataset)

                # Refinement with target restraints
                protRefinePhil = self._runRefine(
                    objLabel="dials - static refinement with restraints",
                    inputSet=protIndexPhil.outputIndexedSpots,
                    scanVaryingNew=dconst.STATIC,
                    detectorFixDistance=False,
                    useRestraint=True,
                    targetUnitCell=experiment["unit_cell"],
                    targetSigmas=experiment["unit_cell_sigmas"],
                )
                self.assertCommand(
                    protRefinePhil,
                    cc.lysoRefineCommandLine(
                        inputExtraPath=protIndexPhil._getExtraPath(),
                        outputExtraPath=protRefinePhil._getExtraPath(),
                        logPath=protRefinePhil._getLogsPath(),
                        static=True,
                    ),
                    program="dials.refine",
                )
                refinedsetPhil = getattr(
                    protRefinePhil, "outputRefinedSpots", None
                )
                self.assertIsNotNone(protRefinePhil.outputRefinedSpots)
                self.assertFileExists(refinedsetPhil.getDialsModel())
                self.assertFileExists(refinedsetPhil.getDialsRefl())
                # Check that the correct phil file is actually created
                self.assertFileExists(protRefinePhil.getRestraintsPhil())
                referenceFile = os.path.abspath("restraints.phil")
                referencePhil = writeRestraintsPhil(
                    fn=referenceFile,
                    values=experiment["unit_cell"],
                    sigmas=experiment["unit_cell_sigmas"],
                )
                self.assertFileExists(referenceFile)
                self.comparePhils(
                    goodPhil=referencePhil,
                    testPhil=protRefinePhil.getRestraintsPhil(),
                )
                self.checkLogDataset(protRefinePhil, dataset)

            # Run refinement
            protRefine = self._runRefine(
                objLabel="dials - static refinement",
                inputSet=protIndex.outputIndexedSpots,
                scanVaryingNew=dconst.UNSET,
                detectorFixAll=True,
            )
            self.assertCommand(
                protRefine,
                cc.lysoRefineCommandLine(
                    inputExtraPath=indexExtra,
                    outputExtraPath=protRefine._getExtraPath(),
                    logPath=protRefine._getLogsPath(),
                ),
                program="dials.refine",
            )
            refinedset = getattr(protRefine, "outputRefinedSpots", None)
            self.assertIsNotNone(protRefine.outputRefinedSpots)
            self.assertFileExists(refinedset.getDialsModel())
            self.assertFileExists(refinedset.getDialsRefl())
            self.checkLogDataset(protRefine, dataset)

            # Run scan-varying refinement
            protSvRefine = self._runRefine(
                objLabel="dials - scan-varying refinement",
                inputSet=protRefine.outputRefinedSpots,
                scanVaryingNew=dconst.SCAN_VARYING,
                beamFixAll=False,
                beamFixInSpindlePlane=False,
                beamFixOutSpindlePlane=False,
                beamFixWavelength=True,
                beamForceStatic=False,
                detectorFixAll=True,
            )
            self.assertCommand(
                protSvRefine,
                cc.lysoRefineCommandLine(
                    inputExtraPath=protRefine._getExtraPath(),
                    outputExtraPath=protSvRefine._getExtraPath(),
                    logPath=protSvRefine._getLogsPath(),
                    scan_varying=True,
                ),
                program="dials.refine",
            )
            svRefinedset = getattr(protSvRefine, "outputRefinedSpots", None)
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
            self.assertCommand(
                protSvRefineOld,
                cc.lysoRefineCommandLine(
                    inputExtraPath=protRefine._getExtraPath(),
                    outputExtraPath=protSvRefineOld._getExtraPath(),
                    logPath=protSvRefineOld._getLogsPath(),
                    scan_varying=True,
                ),
                program="dials.refine",
            )
            svRefinedsetOld = getattr(
                protSvRefineOld, "outputRefinedSpots", None
            )
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
            self.assertCommand(
                protIntegrate,
                cc.integrateCommandLine(
                    inputExtraPath=protSvRefine._getExtraPath(),
                    outputExtraPath=protIntegrate._getExtraPath(),
                    logPath=protIntegrate._getLogsPath(),
                    experiment=experiment,
                ),
                program="dials.integrate",
            )
            integratedset = getattr(
                protIntegrate, "outputIntegratedSpots", None
            )
            self.assertIsNotNone(protIntegrate.outputIntegratedSpots)
            self.assertFileExists(integratedset.getDialsModel())
            self.assertFileExists(integratedset.getDialsRefl())
            self.checkLogDataset(
                protIntegrate, dataset, "Summary vs resolution"
            )

            # Check symmetry and scale
            protSymmetry = self._runSymmetry(
                objLabel="dials - symmetry check",
                inputSet=protIntegrate.outputIntegratedSpots,
            )

            self.assertCommand(
                protSymmetry,
                cc.symmetryCommandLine(
                    inputExtraPath=protIntegrate._getExtraPath(),
                    outputExtraPath=protSymmetry._getExtraPath(),
                    logPath=protSymmetry._getLogsPath(),
                ),
                "dials.symmetry",
            )
            symmetrizedset = getattr(
                protSymmetry, "outputSymmetrizedSpots", None
            )
            self.assertIsNotNone(protSymmetry.outputSymmetrizedSpots)
            self.assertFileExists(symmetrizedset.getDialsModel())
            self.assertFileExists(symmetrizedset.getDialsRefl())
            self.checkLogDataset(
                protSymmetry,
                dataset,
                f"Recommended space group: {experiment['space_group']}",
            )

            protScaling = self._runScaling(
                objLabel="dials - scaling",
                inputSets=[protIntegrate.outputIntegratedSpots],
            )
            self.assertCommand(
                protScaling,
                cc.scalingCommandLine(
                    inputExtraPath=protIntegrate._getExtraPath(),
                    outputExtraPath=protScaling._getExtraPath(),
                    logPath=protScaling._getLogsPath(),
                ),
                "dials.scale",
            )
            scaledset = getattr(protScaling, "outputScaledSpots", None)
            self.assertIsNotNone(protScaling.outputScaledSpots)
            self.assertFileExists(scaledset.getDialsModel())
            self.assertFileExists(scaledset.getDialsRefl())
            self.checkLogDataset(
                protScaling,
                dataset,
                "Space group being used during scaling is P 4",
            )

            scaledSets.append(protScaling.outputScaledSpots)
            scaleProt.append(protScaling)

        with self.subTest(msg="Scaling multiple datasets"):
            while len(scaledSets) < 2:
                self.skipTest("Not enough datasets to test scaling multiple")
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
            scaledExtra = [prot._getExtraPath() for prot in scaleProt]
            mergedMtzFile = f"{self.getOutputPath()}/{mergedFn}"
            unmergedMtzFile = f"{self.getOutputPath()}/{unmergedFn}"
            self.assertCommand(
                protMultiScaling,
                cc.lysoMultiScalingCommandLine(
                    inputExtraList=scaledExtra,
                    outputExtraPath=protMultiScaling._getExtraPath(),
                    logPath=protMultiScaling._getLogsPath(),
                    mergedFile=mergedMtzFile,
                    unmergedFile=unmergedMtzFile,
                    crystalName=crystalName,
                    projectName=self.PROJECT_NAME,
                ),
                "dials.scale",
            )
            multiscaledset = getattr(
                protMultiScaling, "outputScaledSpots", None
            )
            self.assertIsNotNone(protMultiScaling.outputScaledSpots)
            self.assertFileExists(mergedMtzFile)
            self.assertFileExists(unmergedMtzFile)
            self.assertFileExists(multiscaledset.getDialsModel())
            self.assertFileExists(multiscaledset.getDialsRefl())
            compareDatasets = "Source of data:"
            for ds in multiDataset:
                compareDatasets += (
                    f"\n{os.path.join(self.dataPath, ds, 'SMV/data')}"
                )
            self.assertEqual(protMultiScaling.getDatasets(), compareDatasets)

            # Test merging protocol
            mergeBins = 5
            protMerge = self._runMerging(
                inputSet=protMultiScaling.outputScaledSpots,
                xtalName=crystalName,
                nBins=mergeBins,
            )

            self.assertCommand(
                protMerge,
                cc.lysoMergeCommandLine(
                    inputExtraPath=protMultiScaling._getExtraPath(),
                    outputExtraPath=protMerge._getExtraPath(),
                    logPath=protMerge._getLogsPath(),
                    crystalName=crystalName,
                    projectName=self.PROJECT_NAME,
                    mergeBins=mergeBins,
                ),
                "dials.merge",
            )
            mergedset = getattr(protMerge, "exportedFileSet", None)
            self.assertIsNotNone(protMerge.exportedFileSet)
            for exportFile in mergedset:
                self.assertFileExists(exportFile.getFilePath())
            self.assertEqual(protMerge.getDatasets(), compareDatasets)

            # Test export of output scaled together
            protExportMtz = self._runExport(
                inputSet=protMultiScaling.outputScaledSpots,
                exportFormat=dconst.MTZ,
            )
            self.assertCommand(
                protExportMtz,
                cc.exportMtzCommandLine(
                    inputExtraPath=protMultiScaling._getExtraPath(),
                    outputExtraPath=protExportMtz._getExtraPath(),
                    logPath=protExportMtz._getLogsPath(),
                    objectID=protExportMtz.getObjId(),
                    projectName=self.PROJECT_NAME,
                ),
                "dials.export",
            )
            exportedmtzset = getattr(protExportMtz, "exportedFileSet", None)
            self.assertIsNotNone(protExportMtz.exportedFileSet)
            for exportFile in exportedmtzset:
                self.assertFileExists(exportFile.getFilePath())
            self.assertEqual(protExportMtz.getDatasets(), compareDatasets)

        with self.subTest(msg="Excluding images"):
            exclusions = ["0:1:2", "0:10:12"]
            protScalingExclude = self._runScaling(
                objLabel="dials - scaling with exclusions",
                inputSets=[protIntegrate.outputIntegratedSpots],
                excludeImages=True,
                numberOfExclusions=2,
                imageGroup1=exclusions[0],
                imageGroup2=exclusions[1],
            )
            self.assertCommand(
                protScalingExclude,
                cc.lysoScaleExclusionCommandLine(
                    inputExtraPath=protIntegrate._getExtraPath(),
                    outputExtraPath=protScalingExclude._getExtraPath(),
                    logPath=protScalingExclude._getLogsPath(),
                    exclusionList=exclusions,
                ),
                "dials.scale",
            )
            scaledset = getattr(protScalingExclude, "outputScaledSpots", None)
            self.assertIsNotNone(protScalingExclude.outputScaledSpots)
            self.assertFileExists(scaledset.getDialsModel())
            self.assertFileExists(scaledset.getDialsRefl())
            self.checkLogDataset(
                protScalingExclude,
                dataset,
                "Space group being used during scaling is P 4",
            )

    def test_garnet_pipeline(self):
        if cc.SKIP_GARNET:
            self.skipTest("Skipping garnet pipeline test")

        # Define all experiment variables in one place
        experiment = cc.garnet_experiment
        exptId = experiment["location"]
        dataset = os.path.join(self.dataPath, exptId)

        # Run import
        protImport = self._runDialsImportImages(
            filesPath=dataset,
            filesPattern=experiment["files_pattern"],
            objLabel=f"dials - import diffraction images\n{exptId}",
            replaceRotationAxis=experiment["replace_rotation_axis"],
            rotationAxis=experiment["rotation_axis"],
            useTemplate=experiment["useTemplate"],
            tsReplacement=experiment["tsReplacement"],
        )

        self.assertCommand(
            protImport,
            cc.garnetImportCommandLine(
                dataset=dataset,
                extraPath=protImport._getExtraPath(),
                logPath=protImport._getLogsPath(),
                experiment=experiment,
            ),
            program="dials.import",
        )
        outputModel = protImport.getOutputModelFile()
        self.assertFileExists(outputModel)
        importedset = getattr(protImport, "outputDiffractionImages", None)
        self.assertIsNotNone(importedset)
        setModel = importedset.getDialsModel()
        self.assertFileExists(setModel)
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protImport.getDatasets(), datasetString)
        outputCompare = protImport.getLogOutput().split("\n")[0]
        self.assertEqual(outputCompare.strip(), "".strip())

        # Run find spots
        protFindSpots = self._runFindSpots(
            protImport.outputDiffractionImages,
            minSpotSize=experiment["min_spot_size"],
            maxSpotSize=experiment["max_spot_size"],
            maxSeparation=experiment["max_separation"],
            dMin=experiment["d_min"],
            dMax=experiment["d_max"],
            sigmaBackground=experiment["sigma_background"],
            thresholdAlgorithm=dconst.DISPERSION,
        )
        inputModel = protFindSpots.getInputModelFile()
        self.assertFileExists(inputModel)
        self.assertCommand(
            protFindSpots,
            cc.garnetStandardFindSpotCommandLine(
                inputExtraPath=protImport._getExtraPath(),
                outputExtraPath=protFindSpots._getExtraPath(),
                logPath=protFindSpots._getLogsPath(),
                experiment=experiment,
            ),
            program="dials.find_spots",
        )
        foundspotset = getattr(protFindSpots, "outputDiffractionSpots", None)
        self.assertIsNotNone(foundspotset)
        self.assertEqual(foundspotset.getSpots(), experiment["found_spots"])
        self.assertFileExists(foundspotset.getDialsModel())
        self.assertFileExists(foundspotset.getDialsRefl())
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protFindSpots.getDatasets(), datasetString)
        outputCompare = protFindSpots.getLogOutput().split("\n")[0]
        self.assertEqual(
            outputCompare.strip(),
            "Histogram of per-image spot count for imageset 0:".strip(),
        )
        with self.subTest(
            msg="Testing validation of resolution limits in Find Spots protocol"
        ):
            with self.assertRaises(Exception):
                self._runFindSpots(
                    protImport.outputDiffractionImages,
                    dMin=experiment["d_max"],
                    dMax=experiment["d_min"],
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

        self.assertEqual(
            protIndex._prepIndexCommandline("dials.index"),
            cc.garnetIndexCommandLine(
                modelPath=protImport._getExtraPath(),
                spotsPath=protFindSpots._getExtraPath(),
                logPath=indexLogs,
                indexTmp=indexTmp,
            ),
        )
        self.assertEqual(
            protIndex._prepBravaisCommandline("dials.refine_bravais_settings"),
            cc.refineBravaisLatticeCommandLine(
                indexTmp=indexTmp, indexLogs=indexLogs
            ),
        )
        self.assertEqual(
            protIndex._prepReindexCommandline(),
            cc.reindexCommandLine(experiment=experiment, indexTmp=indexTmp),
        )
        indexedset = getattr(protIndex, "outputIndexedSpots", None)
        self.assertIsNotNone(protIndex.outputIndexedSpots)
        self.assertFileExists(indexedset.getDialsModel())
        self.assertFileExists(indexedset.getDialsRefl())
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protIndex.getDatasets(), datasetString)
        outputCompare = protIndex.getLogOutput().split("\n")[0]
        self.assertEqual(outputCompare.strip(), "".strip())
        with self.subTest(
            msg="Testing refine_bravais_settings with copied parameters"
        ):
            protIndexCopy = self._runIndex(
                objLabel="dials - index and refine bravais setting with copied parameters",
                inputImages=protImport.outputDiffractionImages,
                inputSpots=protFindSpots.outputDiffractionSpots,
                doRefineBravaisSettings=True,
                doReindex=False,
            )
            self.assertEqual(
                protIndexCopy._prepBravaisCommandline(
                    "dials.refine_bravais_settings"
                ),
                cc.refineBravaisLatticeCommandLine(
                    indexTmp=protIndexCopy._getTmpPath(),
                    indexLogs=protIndexCopy._getLogsPath(),
                    copied=True,
                ),
            )

        # Run refinement
        protRefine = self._runRefine(
            inputSet=protIndex.outputIndexedSpots,
            scanVaryingNew=dconst.UNSET,
        )
        self.assertCommand(
            protRefine,
            cc.garnetRefineCommandLine(
                inputExtraPath=protIndex._getExtraPath(),
                outputExtraPath=protRefine._getExtraPath(),
                logPath=protRefine._getLogsPath(),
            ),
            program="dials.refine",
        )
        refinedset = getattr(protRefine, "outputRefinedSpots", None)
        self.assertIsNotNone(protRefine.outputRefinedSpots)
        self.assertFileExists(refinedset.getDialsModel())
        self.assertFileExists(refinedset.getDialsRefl())
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protRefine.getDatasets(), datasetString)
        outputCompare = protRefine.getLogOutput().split("\n")[0]
        self.assertEqual(outputCompare.strip(), "".strip())

        # Run integration
        protIntegrate = self._runIntegrate(
            inputSet=protRefine.outputRefinedSpots,
            nproc=experiment["nproc"],
            dMin=experiment["d_min"],
            makeReport=True,
        )
        self.assertCommand(
            protIntegrate,
            cc.integrateCommandLine(
                inputExtraPath=protRefine._getExtraPath(),
                outputExtraPath=protIntegrate._getExtraPath(),
                logPath=protIntegrate._getLogsPath(),
                experiment=experiment,
            ),
            program="dials.integrate",
        )
        integratedset = getattr(protIntegrate, "outputIntegratedSpots", None)
        self.assertIsNotNone(protIntegrate.outputIntegratedSpots)
        self.assertFileExists(integratedset.getDialsModel())
        self.assertFileExists(integratedset.getDialsRefl())
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protIntegrate.getDatasets(), datasetString)
        outputCompare = protIntegrate.getLogOutput().split("\n")[0]
        self.assertEqual(
            outputCompare.strip(), "Summary vs resolution".strip()
        )
        with self.subTest(
            msg="Testing validation of resolution limits in Integrate protocol"
        ):
            with self.assertRaises(Exception):
                self._runIntegrate(
                    inputSet=protRefine.outputRefinedSpots,
                    dMin=experiment["d_max"],
                    dMax=experiment["d_min"],
                )

        # Check symmetry and scale
        protSymmetry = self._runSymmetry(
            objLabel="dials - symmetry check",
            inputSet=protIntegrate.outputIntegratedSpots,
        )
        self.assertCommand(
            protSymmetry,
            cc.symmetryCommandLine(
                inputExtraPath=protIntegrate._getExtraPath(),
                outputExtraPath=protSymmetry._getExtraPath(),
                logPath=protSymmetry._getLogsPath(),
            ),
            program="dials.symmetry",
        )
        symmetrizedset = getattr(protSymmetry, "outputSymmetrizedSpots", None)
        self.assertIsNotNone(protSymmetry.outputSymmetrizedSpots)
        self.assertFileExists(symmetrizedset.getDialsModel())
        self.assertFileExists(symmetrizedset.getDialsRefl())
        datasetString = f"Source of data:\n{dataset}"
        self.assertEqual(protSymmetry.getDatasets(), datasetString)
        outputCompare = protSymmetry.getLogOutput().split("\n")[0]
        self.assertEqual(
            outputCompare.strip(), "Recommended space group: I 4 3 2".strip()
        )

        protScaling = self._runScaling(
            objLabel="dials - scaling",
            inputSets=[protIntegrate.outputIntegratedSpots],
        )
        self.assertCommand(
            protScaling,
            cc.scalingCommandLine(
                inputExtraPath=protIntegrate._getExtraPath(),
                outputExtraPath=protScaling._getExtraPath(),
                logPath=protScaling._getLogsPath(),
            ),
            program="dials.scale",
        )
        scaledset = getattr(protScaling, "outputScaledSpots", None)
        self.assertIsNotNone(protScaling.outputScaledSpots)
        self.assertFileExists(scaledset.getDialsModel())
        self.assertFileExists(scaledset.getDialsRefl())
        self.assertEqual(protScaling.getDatasets(), datasetString)
        outputCompare = protScaling.getLogOutput().split("\n")[0]
        self.assertEqual(
            outputCompare.strip(),
            "Space group being used during scaling is I 2 3".strip(),
        )
        with self.subTest(
            msg="Testing validation of resolution limits in Scale protocol"
        ):
            with self.assertRaises(Exception):
                self._runScaling(
                    inputSets=[protIntegrate.outputIntegratedSpots],
                    dMin=experiment["d_max"],
                    dMax=experiment["d_min"],
                )
        crystalName = "Garnet"
        protExportMtz = self._runExport(
            inputSet=protScaling.outputScaledSpots,
            exportFormat=dconst.MTZ,
            mtzCrystalName=crystalName,
        )
        self.assertCommand(
            protExportMtz,
            cc.exportMtzCommandLine(
                inputExtraPath=protScaling._getExtraPath(),
                outputExtraPath=protExportMtz._getExtraPath(),
                logPath=protExportMtz._getLogsPath(),
                objectID=protExportMtz.getObjId(),
                projectName=self.PROJECT_NAME,
                crystalName=crystalName,
            ),
            program="dials.export",
        )
        exportedmtzset = getattr(protExportMtz, "exportedFileSet", None)
        self.assertIsNotNone(protExportMtz.exportedFileSet)
        for exportFile in exportedmtzset:
            self.assertFileExists(exportFile.getFilePath())
        self.assertEqual(protExportMtz.getDatasets(), datasetString)


class TestEdDialsUtils(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        if cc.SKIP_UTILS:
            cls.skipTest(cls, "Skipping utils")
        pwtests.setupTestOutput(cls)
        cls.dataPath = os.environ.get("SCIPION_ED_TESTDATA")

        if not os.path.exists(cls.dataPath):
            raise MissingPathException(
                f"Can not run utils tests, missing file:\n  {cls.dataPath}"
            )

    @classmethod
    def tearDownClass(cls):
        if not cc.KEEP_UTILS_TEST_OUTPUT:
            # Clean up all output files from the test
            pw.utils.cleanPath(cls.getOutputPath())

    def comparePhils(self, goodPhil="restraints.phil", testPhil=None):
        self.assertIsNotNone(testPhil)
        with open(os.path.join(self.dataPath, "utils", goodPhil)) as f1:
            with open(testPhil) as f2:
                self.assertEqual(f1.read(), f2.read())

    def test_write_restraints(self):
        # Input and output are as expected
        setFn = self.getOutputPath("restraints.phil")
        pw.utils.cleanPath(setFn)

        # Add some values as a start
        values = "10,20,30,90,90,90"
        sigmas = "0.05,0.05,0.05,0.05,0.05,0.05"

        outFn = writeRestraintsPhil(fn=setFn, values=values, sigmas=sigmas)
        self.comparePhils(goodPhil="restraints.phil", testPhil=outFn)

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

        self.comparePhils(goodPhil="restraints.phil", testPhil=fnShort)
        self.comparePhils(goodPhil="restraints.phil", testPhil=fnLong)

    def test_write_restraints_no_sigmas(self):
        setFn = self.getOutputPath("restraints_no_sigmas.phil")
        pw.utils.cleanPath(setFn)

        values = "10,20,30,90,90,90"
        writeRestraintsPhil(fn=setFn, values=values)
        self.comparePhils(goodPhil="restraints_no_sigmas.phil", testPhil=setFn)
