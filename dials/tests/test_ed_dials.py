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
import io
import json

import msgpack

import pyworkflow as pw
import pyworkflow.tests as pwtests

import pwed
from pwed.objects import *
from pwed.protocols import ProtImportDiffractionImages

from dials.protocols import *
from dials.convert import writeJson, readRefl, writeRefl
from dials.constants import *
import dials.utils as dutils


pw.Config.setDomain(pwed)


class TestEdDialsProtocols(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls, writeLocalConfig=True)
        cls.dataPath = os.environ.get("SCIPION_ED_TESTDATA")

        if not os.path.exists(cls.dataPath):
            raise Exception("Can not run DIALS tests, missing file:\n  %s"
                            % cls.dataPath)

    # Functions for running protocols

    def _runImportImages(self, filesPath, filesPattern='SMV/data/00{TI}.img', **kwargs):
        protImport = self.newProtocol(
            ProtImportDiffractionImages,
            filesPath=filesPath,
            filesPattern=filesPattern,
            **kwargs)
        self.launchProtocol(protImport)
        return protImport

    def _runDialsImportImages(self, filesPath, filesPattern='SMV/data/00{TI}.img', **kwargs):
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

    def _runExport(self, inputSet, **kwargs):
        protExport = self.newProtocol(DialsProtExport,
                                      inputSet=inputSet,
                                      **kwargs)
        self.launchProtocol(protExport)
        return protExport

    # Helper functions
    def makePathList(self, scanRange, digits=3, location=None, pattern='SMV/data/00{TI}.img', wildcard='{TI}'):
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
            'location': 'lysozyme/experiment_14',
            'files_pattern': 'SMV/data/00{TI}.img',
            'd_min': 2.5,
            'found_spots': 1322,
            'distance': 1480.56,
            'rotation_axis': '-0.6204,-0.7843,0.0',
            'lyso': True,
            'sample': 'lyso',
            'scan_range': '1,49',
            'space_group': 'P 2',
        }
        lyso_experiment_24 = {
            'location': 'lysozyme/experiment_24',
            'files_pattern': 'SMV/data/00{TI}.img',
            'd_min': 3.5,
            'found_spots': 847,
            'distance': None,
            'rotation_axis': '-0.6204,-0.7843,0.0',
            'lyso': True,
            'sample': 'lysozyme',
            'scan_range': '1,46',
            'space_group': 'P 4 2 2',
        }
        experiments.append(lyso_experiment_14)
        experiments.append(lyso_experiment_24)
        return experiments

    def getGarnetExperiment(self):
        experiment = {
            'location': 'garnet/experiment_1',
            'files_pattern': 'Garnet_Cont_{TI}.img',
            'sample': 'garnet',
            'replace_rotation_axis': True,
            'rotation_axis': '-0.755,-0.656,0.0',
            'useTemplate': True,
            'tsReplacement': '0###',
            'min_spot_size': 10,
            'max_spot_size': 1500000,
            'max_separation': 7.5,
            'd_min': 0.5,
            'd_max': 3.0,
            'sigma_background': 1.0,
            'found_spots': 6223,
            'nproc': 8,
            'bravais_setting_reference': 'bravais_setting_22.expt'
        }
        return experiment

    def checkLogDataset(self, protocol, dataset, logOutput=''):
        datasetString = "Source of data:\n{}".format(
            os.path.join(dataset, 'SMV/data'))
        self.assertEqual(protocol.getDatasets(), datasetString)
        outputCompare = protocol.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(), logOutput.strip())

    def assertCommand(self, protocol, commandString, program=None):
        if program != None:
            try:
                CL = protocol._prepCommandline(program)
            except AttributeError:
                CL = protocol._prepareCommandLineArguments(program)
        else:
            try:
                CL = protocol._prepCommandline()
            except AttributeError:
                CL = protocol._prepareCommandLineArguments()
        self.assertEqual(CL, commandString)

    # Pipelines
    def test_lyso_pipeline(self):

        scaledSets = []
        scaleProt = []
        multiDataset = []
        for experiment in self.getLysoTestExperiments():
            exptId = experiment['location']
            with self.subTest(msg='Pipeline using {}'.format(exptId), experiment=exptId):

                dataset = os.path.join(self.dataPath, exptId)
                multiDataset.append(dataset)

                # Run import
                protImport = self._runDialsImportImages(
                    filesPath=dataset,
                    filesPattern=experiment['files_pattern'],
                    objLabel="dials - import diffraction images\n"
                    "{}".format(exptId),
                    replaceRotationAxis=True,
                    rotationAxis=experiment['rotation_axis'],
                    overwriteDetectorDistance=experiment['distance'],
                )
                if experiment['distance'] != None:
                    distance = " distance={}".format(experiment['distance'])
                else:
                    distance = ''

            pathlist = self.makePathList(
                experiment['scan_range'], location=dataset, pattern=experiment['files_pattern'])
            importCL = '{} output.log={}/dials.import.log output.experiments={}/imported.expt goniometer.axes={}{}'.format(
                " ".join(pathlist),
                protImport._getLogsPath(),
                protImport._getExtraPath(),
                experiment['rotation_axis'],
                distance
            )

            self.assertCommand(protImport, importCL,
                               program='dials.import')
            outputModel = protImport.getOutputModelFile()
            self.assertTrue(os.path.exists(outputModel))
            importedset = getattr(
                protImport, 'outputDiffractionImages', None)
            self.assertIsNotNone(importedset)
            setModel = importedset.getDialsModel()
            self.assertTrue(os.path.exists(setModel))
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
            self.assertTrue(os.path.exists(inputModel))
            findSpotCL = '{}/imported.expt output.log={}/dials.find_spots.log output.reflections={}/strong.refl spotfinder.scan_range={} spotfinder.filter.d_min={} spotfinder.filter.max_spot_size=1000 spotfinder.filter.max_strong_pixel_fraction=0.25 spotfinder.filter.max_separation=2.0 spotfinder.filter.untrusted.rectangle=0,516,255,261 spotfinder.filter.untrusted.rectangle=255,261,0,516 spotfinder.threshold.algorithm=dispersion_extended spotfinder.threshold.dispersion.sigma_background=6.0 spotfinder.threshold.dispersion.sigma_strong=3.0 spotfinder.threshold.dispersion.kernel_size=3,3'.format(
                protFindSpots._getExtraPath(),
                protFindSpots._getLogsPath(),
                protFindSpots._getExtraPath(),
                experiment['scan_range'],
                experiment['d_min'],)
            self.assertCommand(protFindSpots, findSpotCL)
            foundspotset = getattr(
                protFindSpots, 'outputDiffractionSpots', None)
            self.assertIsNotNone(foundspotset)
            self.assertEqual(foundspotset.getSpots(),
                             experiment['found_spots'])
            self.assertTrue(os.path.exists(foundspotset.getDialsModel()))
            self.assertTrue(os.path.exists(foundspotset.getDialsRefl()))
            self.checkLogDataset(protFindSpots, dataset,
                                 'Histogram of per-image spot count for imageset 0:')

            # Run indexing
            protIndex = self._runIndex(
                objLabel="dials - index and reindex",
                inputImages=protImport.outputDiffractionImages,
                inputSpots=protFindSpots.outputDiffractionSpots,
                doRefineBravaisSettings=True,
                doReindex=True,
                doReindexReflections=True,
                detectorFixPosition=True,
                detectorFixOrientation=True,
                detectorFixDistance=True,
                beamFixInSpindlePlane=True,
                beamFixOutSpindlePlane=True,
                beamFixWavelength=True,
                # FIXME: Command line input should not be needed
                commandLineInputBravais='beam.fix=all detector.fix=all',
            )
            indexTmp = protIndex._getTmpPath()
            indexLogs = protIndex._getLogsPath()
            indexExtra = protIndex._getExtraPath()
            indexCL = "{}/imported.expt {}/strong.refl output.log={}/dials.index.log output.experiments={}/indexed.expt output.reflections={}/indexed.refl refinement.parameterisation.beam.fix='*all in_spindle_plane out_spindle_plane wavelength' refinement.parameterisation.crystal.fix='all cell orientation' refinement.parameterisation.detector.fix='*all position orientation distance'".format(
                protFindSpots._getExtraPath(),
                protFindSpots._getExtraPath(),
                indexLogs,
                indexTmp,
                indexTmp,
            )
            refBravCL = "{}/indexed.expt {}/indexed.refl output.log={}/dials.refine_bravais_settings.log output.directory={} beam.fix=all detector.fix=all".format(
                indexTmp,
                indexTmp,
                indexLogs,
                indexTmp,
            )
            reindexCL = 'change_of_basis_op=a,b,c {}/indexed.refl output.reflections={}/reindexed.refl'.format(
                indexTmp,
                indexTmp,
            )
            self.assertEqual(protIndex._prepIndexCommandline(
                'dials.index'), indexCL)
            self.assertEqual(protIndex._prepBravaisCommandline(
                'dials.refine_bravais_settings'), refBravCL)
            self.assertEqual(protIndex._prepReindexCommandline(
                'dials.reindex'), reindexCL)
            indexedset = getattr(protIndex, 'outputIndexedSpots', None)
            self.assertIsNotNone(protIndex.outputIndexedSpots)
            self.assertTrue(os.path.exists(indexedset.getDialsModel()))
            self.assertTrue(os.path.exists(indexedset.getDialsRefl()))
            self.checkLogDataset(protIndex, dataset)

            # Run refinement
            protRefine = self._runRefine(
                objLabel="dials - static refinement",
                inputSet=protIndex.outputIndexedSpots,
                scanVarying=False,
                detectorFixAll=True,
            )
            refineCLstatic = "{}/indexed.expt {}/indexed.refl output.log={}/dials.refine.log output.experiments={}/refined.expt output.reflections={}/refined.refl beam.fix='all *in_spindle_plane out_spindle_plane *wavelength' detector.fix='*all position orientation distance'".format(
                indexExtra,
                indexExtra,
                protRefine._getLogsPath(),
                protRefine._getExtraPath(),
                protRefine._getExtraPath(),
            )
            self.assertCommand(protRefine, refineCLstatic,
                               program='dials.refine')
            refinedset = getattr(protRefine, 'outputRefinedSpots', None)
            self.assertIsNotNone(protRefine.outputRefinedSpots)
            self.assertTrue(os.path.exists(refinedset.getDialsModel()))
            self.assertTrue(os.path.exists(refinedset.getDialsRefl()))
            self.checkLogDataset(protRefine, dataset)

            # Run scan-varying refinement
            protSvRefine = self._runRefine(
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
            refineCLsv = "{}/refined.expt {}/refined.refl output.log={}/dials.refine.log output.experiments={}/refined.expt output.reflections={}/refined.refl scan_varying=True beam.fix='all in_spindle_plane out_spindle_plane *wavelength' beam.force_static=False detector.fix='*all position orientation distance'".format(
                protRefine._getExtraPath(),
                protRefine._getExtraPath(),
                protSvRefine._getLogsPath(),
                protSvRefine._getExtraPath(),
                protSvRefine._getExtraPath(),
            )
            self.assertCommand(protSvRefine, refineCLsv,
                               program='dials.refine')
            svRefinedset = getattr(
                protSvRefine, 'outputRefinedSpots', None)
            self.assertIsNotNone(protSvRefine.outputRefinedSpots)
            self.assertTrue(os.path.exists(svRefinedset.getDialsModel()))
            self.assertTrue(os.path.exists(svRefinedset.getDialsRefl()))
            self.checkLogDataset(protSvRefine, dataset)

            # Run integration
            protIntegrate = self._runIntegrate(
                # objLabel="Dials integration: {}".format(exptId),
                inputSet=protSvRefine.outputRefinedSpots,
                nproc=8,
                commandLineInput='prediction.d_min={}'.format(
                    experiment['d_min']),
            )
            integrateCL = "{}/refined.expt {}/refined.refl output.log={}/dials.integrate.log output.experiments={}/integrated_model.expt output.reflections={}/integrated_reflections.refl nproc=8 prediction.d_min={}".format(
                protSvRefine._getExtraPath(),
                protSvRefine._getExtraPath(),
                protIntegrate._getLogsPath(),
                protIntegrate._getExtraPath(),
                protIntegrate._getExtraPath(),
                experiment['d_min']
            )
            self.assertCommand(protIntegrate, integrateCL,
                               program='dials.integrate')
            integratedset = getattr(
                protIntegrate, 'outputIntegratedSpots', None)
            self.assertIsNotNone(protIntegrate.outputIntegratedSpots)
            self.assertTrue(os.path.exists(integratedset.getDialsModel()))
            self.assertTrue(os.path.exists(integratedset.getDialsRefl()))
            self.checkLogDataset(protIntegrate, dataset,
                                 'Summary vs resolution')

            # Check symmetry and scale
            protSymmetry = self._runSymmetry(
                objLabel="dials - symmetry check",
                inputSet=protIntegrate.outputIntegratedSpots,
            )
            symmCL = "{}/integrated_model.expt {}/integrated_reflections.refl output.log={}/dials.symmetry.log output.experiments={}/symmetrized.expt output.reflections={}/symmetrized.refl".format(
                protIntegrate._getExtraPath(),
                protIntegrate._getExtraPath(),
                protSymmetry._getLogsPath(),
                protSymmetry._getExtraPath(),
                protSymmetry._getExtraPath(),
            )
            self.assertCommand(protSymmetry, symmCL, 'dials.symmetry')
            symmetrizedset = getattr(
                protSymmetry, 'outputSymmetrizedSpots', None)
            self.assertIsNotNone(protSymmetry.outputSymmetrizedSpots)
            self.assertTrue(os.path.exists(symmetrizedset.getDialsModel()))
            self.assertTrue(os.path.exists(symmetrizedset.getDialsRefl()))
            self.checkLogDataset(protSymmetry, dataset,
                                 'Recommended space group: {}'.format(experiment['space_group']))

            protScaling = self._runScaling(
                objLabel="dials - scaling",
                inputSets=[protIntegrate.outputIntegratedSpots],
            )
            scaleCL = "{}/integrated_model.expt {}/integrated_reflections.refl output.log={}/dials.scale.log output.experiments={}/scaled.expt output.reflections={}/scaled.refl output.html={}/dials.scale.html filtering.output.scale_and_filter_results={}/scale_and_filter_results.json cut_data.partiality_cutoff=0.4 cut_data.min_isigi=-5.0 outlier_rejection=standard outlier_zmax=6.0 filtering.method=None filtering.deltacchalf.max_cycles=6 filtering.deltacchalf.max_percent_removed=10.0 filtering.deltacchalf.mode=dataset filtering.deltacchalf.group_size=10 filtering.deltacchalf.stdcutoff=4.0".format(
                protIntegrate._getExtraPath(),
                protIntegrate._getExtraPath(),
                protScaling._getLogsPath(),
                protScaling._getExtraPath(),
                protScaling._getExtraPath(),
                protScaling._getExtraPath(),
                protScaling._getExtraPath(),
            )
            self.assertCommand(protScaling, scaleCL, 'dials.scale')
            scaledset = getattr(
                protScaling, 'outputScaledSpots', None)
            self.assertIsNotNone(protScaling.outputScaledSpots)
            self.assertTrue(os.path.exists(scaledset.getDialsModel()))
            self.assertTrue(os.path.exists(scaledset.getDialsRefl()))
            self.checkLogDataset(
                protScaling, dataset, 'Space group being used during scaling is P 4')

            scaledSets.append(protScaling.outputScaledSpots)
            scaleProt.append(protScaling)

        with self.subTest(msg='Scaling multiple datasets'):
            while len(scaledSets) < 2:
                self.skipTest('Not enough datasets to test scaling multiple')
            protMultiScaling = self._runScaling(
                objLabel="dials - scaling multiple",
                inputSets=scaledSets,
                checkConsistentIndexing=True,
            )
            scaledExtra0 = scaleProt[0]._getExtraPath()
            scaledExtra1 = scaleProt[1]._getExtraPath()
            multiScaleCL = "{}/scaled.expt {}/scaled.refl {}/scaled.expt {}/scaled.refl output.log={}/dials.scale.log output.experiments={}/scaled.expt output.reflections={}/scaled.refl output.html={}/dials.scale.html filtering.output.scale_and_filter_results={}/scale_and_filter_results.json cut_data.partiality_cutoff=0.4 cut_data.min_isigi=-5.0 scaling_options.check_consistent_indexing=True outlier_rejection=standard outlier_zmax=6.0 filtering.method=None filtering.deltacchalf.max_cycles=6 filtering.deltacchalf.max_percent_removed=10.0 filtering.deltacchalf.mode=dataset filtering.deltacchalf.group_size=10 filtering.deltacchalf.stdcutoff=4.0".format(
                scaledExtra0,
                scaledExtra0,
                scaledExtra1,
                scaledExtra1,
                protMultiScaling._getLogsPath(),
                protMultiScaling._getExtraPath(),
                protMultiScaling._getExtraPath(),
                protMultiScaling._getExtraPath(),
                protMultiScaling._getExtraPath(),
            )
            self.assertCommand(protMultiScaling, multiScaleCL, 'dials.scale')
            multiscaledset = getattr(
                protMultiScaling, 'outputScaledSpots', None)
            self.assertIsNotNone(protMultiScaling.outputScaledSpots)
            self.assertTrue(os.path.exists(multiscaledset.getDialsModel()))
            self.assertTrue(os.path.exists(multiscaledset.getDialsRefl()))
            compareDatasets = 'Source of data:'
            for ds in multiDataset:
                compareDatasets += "\n{}".format(
                    os.path.join(self.dataPath, ds, 'SMV/data'))
            self.assertEqual(protMultiScaling.getDatasets(), compareDatasets)

            protExportMtz = self._runExport(
                inputSet=protMultiScaling.outputScaledSpots,
                exportFormat=MTZ,
            )

            exportMtzCL = "{}/scaled.expt {}/scaled.refl format=mtz mtz.hklout={}/integrated_{}.mtz output.log={}/dials.export.log mtz.combine_partials=True mtz.partiality_threshold=0.99 mtz.min_isigi=-5.0 mtz.crystal_name=XTAL".format(
                protMultiScaling._getExtraPath(),
                protMultiScaling._getExtraPath(),
                protExportMtz._getExtraPath(),
                protExportMtz.getObjId(),
                protExportMtz._getLogsPath()
            )
            self.assertCommand(protExportMtz, exportMtzCL, "dials.export")
            exportedmtzset = getattr(protExportMtz, 'exportedFileSet', None)
            self.assertIsNotNone(protExportMtz.exportedFileSet)
            for exportFile in exportedmtzset:
                self.assertTrue(os.path.exists(exportFile.getFilePath()))
            self.assertEqual(protExportMtz.getDatasets(), compareDatasets)

    def test_garnet_pipeline(self):
        # Define all experiment variables in one place
        experiment = self.getGarnetExperiment()
        exptId = experiment['location']
        dataset = os.path.join(self.dataPath, exptId)

        # Run import
        protImport = self._runDialsImportImages(
            filesPath=dataset,
            filesPattern=experiment['files_pattern'],
            objLabel="dials - import diffraction images\n"
            "{}".format(exptId),
            replaceRotationAxis=experiment['replace_rotation_axis'],
            rotationAxis=experiment['rotation_axis'],
            useTemplate=experiment['useTemplate'],
            tsReplacement=experiment['tsReplacement'],
        )
        importCL = 'template={}/Garnet_Cont_{}.img output.log={}/dials.import.log output.experiments={}/imported.expt goniometer.axes={}'.format(
            dataset,
            experiment['tsReplacement'],
            protImport._getLogsPath(),
            protImport._getExtraPath(),
            experiment['rotation_axis'],
        )

        self.assertCommand(protImport, importCL,
                           program='dials.import')
        outputModel = protImport.getOutputModelFile()
        self.assertTrue(os.path.exists(outputModel))
        importedset = getattr(
            protImport, 'outputDiffractionImages', None)
        self.assertIsNotNone(importedset)
        setModel = importedset.getDialsModel()
        self.assertTrue(os.path.exists(setModel))
        datasetString = "Source of data:\n{}".format(dataset)
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
        self.assertTrue(os.path.exists(inputModel))
        findSpotCL = '{}/imported.expt output.log={}/dials.find_spots.log output.reflections={}/strong.refl spotfinder.scan_range=1,566 spotfinder.filter.d_min={} spotfinder.filter.d_max={} spotfinder.filter.min_spot_size={} spotfinder.filter.max_spot_size={} spotfinder.filter.max_strong_pixel_fraction=0.25 spotfinder.filter.max_separation={} spotfinder.threshold.algorithm=dispersion spotfinder.threshold.dispersion.sigma_background={} spotfinder.threshold.dispersion.sigma_strong=3.0 spotfinder.threshold.dispersion.kernel_size=3,3'.format(
            protFindSpots._getExtraPath(),
            protFindSpots._getLogsPath(),
            protFindSpots._getExtraPath(),
            experiment['d_min'],
            experiment['d_max'],
            experiment['min_spot_size'],
            experiment['max_spot_size'],
            experiment['max_separation'],
            experiment['sigma_background'],
        )
        self.assertCommand(protFindSpots, findSpotCL)
        foundspotset = getattr(
            protFindSpots, 'outputDiffractionSpots', None)
        self.assertIsNotNone(foundspotset)
        self.assertEqual(foundspotset.getSpots(),
                         experiment['found_spots'])
        self.assertTrue(os.path.exists(foundspotset.getDialsModel()))
        self.assertTrue(os.path.exists(foundspotset.getDialsRefl()))
        datasetString = "Source of data:\n{}".format(dataset)
        self.assertEqual(protFindSpots.getDatasets(), datasetString)
        outputCompare = protFindSpots.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(),
                         'Histogram of per-image spot count for imageset 0:'.strip())

        # Run indexing
        protIndex = self._runIndex(
            objLabel="dials - index and reindex",
            inputImages=protImport.outputDiffractionImages,
            inputSpots=protFindSpots.outputDiffractionSpots,
            doRefineBravaisSettings=True,
            doReindex=True,
            doReindexReflections=True,
        )
        indexTmp = protIndex._getTmpPath()
        indexLogs = protIndex._getLogsPath()
        indexExtra = protIndex._getExtraPath()
        indexCL = "{}/imported.expt {}/strong.refl output.log={}/dials.index.log output.experiments={}/indexed.expt output.reflections={}/indexed.refl refinement.parameterisation.beam.fix='all *in_spindle_plane out_spindle_plane *wavelength' refinement.parameterisation.crystal.fix='all cell orientation' refinement.parameterisation.detector.fix='all position orientation *distance'".format(
            protFindSpots._getExtraPath(),
            protFindSpots._getExtraPath(),
            indexLogs,
            indexTmp,
            indexTmp,
        )
        refBravCL = "{}/indexed.expt {}/indexed.refl output.log={}/dials.refine_bravais_settings.log output.directory={}".format(
            indexTmp,
            indexTmp,
            indexLogs,
            indexTmp,
        )
        reindexCL = 'change_of_basis_op=a,b,c {}/indexed.refl output.reflections={}/reindexed.refl'.format(
            indexTmp,
            indexTmp,
        )
        self.assertEqual(protIndex._prepIndexCommandline(
            'dials.index'), indexCL)
        self.assertEqual(protIndex._prepBravaisCommandline(
            'dials.refine_bravais_settings'), refBravCL)
        self.assertEqual(protIndex._prepReindexCommandline(
            'dials.reindex'), reindexCL)
        indexedset = getattr(protIndex, 'outputIndexedSpots', None)
        self.assertIsNotNone(protIndex.outputIndexedSpots)
        self.assertTrue(os.path.exists(indexedset.getDialsModel()))
        self.assertTrue(os.path.exists(indexedset.getDialsRefl()))
        datasetString = "Source of data:\n{}".format(dataset)
        self.assertEqual(protIndex.getDatasets(), datasetString)
        outputCompare = protIndex.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(), ''.strip())

        # Run refinement
        protRefine = self._runRefine(
            inputSet=protIndex.outputIndexedSpots,
            scanVarying=False,
        )
        refineCLstatic = "{}/indexed.expt {}/indexed.refl output.log={}/dials.refine.log output.experiments={}/refined.expt output.reflections={}/refined.refl beam.fix='all *in_spindle_plane out_spindle_plane *wavelength' detector.fix='all position orientation distance'".format(
            indexExtra,
            indexExtra,
            protRefine._getLogsPath(),
            protRefine._getExtraPath(),
            protRefine._getExtraPath(),
        )
        self.assertCommand(protRefine, refineCLstatic,
                           program='dials.refine')
        refinedset = getattr(protRefine, 'outputRefinedSpots', None)
        self.assertIsNotNone(protRefine.outputRefinedSpots)
        self.assertTrue(os.path.exists(refinedset.getDialsModel()))
        self.assertTrue(os.path.exists(refinedset.getDialsRefl()))
        datasetString = "Source of data:\n{}".format(dataset)
        self.assertEqual(protRefine.getDatasets(), datasetString)
        outputCompare = protRefine.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(), ''.strip())

        # Run integration
        protIntegrate = self._runIntegrate(
            inputSet=protRefine.outputRefinedSpots,
            nproc=experiment['nproc'],
            dMin=experiment['d_min'],
        )
        integrateCL = "{}/refined.expt {}/refined.refl output.log={}/dials.integrate.log output.experiments={}/integrated_model.expt output.reflections={}/integrated_reflections.refl nproc=8 prediction.d_min={}".format(
            protRefine._getExtraPath(),
            protRefine._getExtraPath(),
            protIntegrate._getLogsPath(),
            protIntegrate._getExtraPath(),
            protIntegrate._getExtraPath(),
            experiment['d_min']
        )
        self.assertCommand(protIntegrate, integrateCL,
                           program='dials.integrate')
        integratedset = getattr(
            protIntegrate, 'outputIntegratedSpots', None)
        self.assertIsNotNone(protIntegrate.outputIntegratedSpots)
        self.assertTrue(os.path.exists(integratedset.getDialsModel()))
        self.assertTrue(os.path.exists(integratedset.getDialsRefl()))
        datasetString = "Source of data:\n{}".format(dataset)
        self.assertEqual(protIntegrate.getDatasets(), datasetString)
        outputCompare = protIntegrate.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(),
                         'Summary vs resolution'.strip())

        # Check symmetry and scale
        protSymmetry = self._runSymmetry(
            objLabel="dials - symmetry check",
            inputSet=protIntegrate.outputIntegratedSpots,
        )
        symmCL = "{}/integrated_model.expt {}/integrated_reflections.refl output.log={}/dials.symmetry.log output.experiments={}/symmetrized.expt output.reflections={}/symmetrized.refl".format(
            protIntegrate._getExtraPath(),
            protIntegrate._getExtraPath(),
            protSymmetry._getLogsPath(),
            protSymmetry._getExtraPath(),
            protSymmetry._getExtraPath(),
        )
        self.assertCommand(protSymmetry, symmCL, 'dials.symmetry')
        symmetrizedset = getattr(
            protSymmetry, 'outputSymmetrizedSpots', None)
        self.assertIsNotNone(protSymmetry.outputSymmetrizedSpots)
        self.assertTrue(os.path.exists(symmetrizedset.getDialsModel()))
        self.assertTrue(os.path.exists(symmetrizedset.getDialsRefl()))
        datasetString = "Source of data:\n{}".format(dataset)
        self.assertEqual(protSymmetry.getDatasets(), datasetString)
        outputCompare = protSymmetry.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(),
                         'Recommended space group: I 4 3 2'.strip())

        protScaling = self._runScaling(
            objLabel="dials - scaling",
            inputSets=[protIntegrate.outputIntegratedSpots],
        )
        scaleCL = "{}/integrated_model.expt {}/integrated_reflections.refl output.log={}/dials.scale.log output.experiments={}/scaled.expt output.reflections={}/scaled.refl output.html={}/dials.scale.html filtering.output.scale_and_filter_results={}/scale_and_filter_results.json cut_data.partiality_cutoff=0.4 cut_data.min_isigi=-5.0 outlier_rejection=standard outlier_zmax=6.0 filtering.method=None filtering.deltacchalf.max_cycles=6 filtering.deltacchalf.max_percent_removed=10.0 filtering.deltacchalf.mode=dataset filtering.deltacchalf.group_size=10 filtering.deltacchalf.stdcutoff=4.0".format(
            protIntegrate._getExtraPath(),
            protIntegrate._getExtraPath(),
            protScaling._getLogsPath(),
            protScaling._getExtraPath(),
            protScaling._getExtraPath(),
            protScaling._getExtraPath(),
            protScaling._getExtraPath(),
        )
        self.assertCommand(protScaling, scaleCL, 'dials.scale')
        scaledset = getattr(
            protScaling, 'outputScaledSpots', None)
        self.assertIsNotNone(protScaling.outputScaledSpots)
        self.assertTrue(os.path.exists(scaledset.getDialsModel()))
        self.assertTrue(os.path.exists(scaledset.getDialsRefl()))
        self.assertEqual(protScaling.getDatasets(), datasetString)
        outputCompare = protScaling.getLogOutput().split('\n')[0]
        self.assertEqual(outputCompare.strip(),
                         'Space group being used during scaling is I 2 3'.strip())

        protExportMtz = self._runExport(
            inputSet=protScaling.outputScaledSpots,
            exportFormat=MTZ,
        )

        exportMtzCL = "{}/scaled.expt {}/scaled.refl format=mtz mtz.hklout={}/integrated_{}.mtz output.log={}/dials.export.log mtz.combine_partials=True mtz.partiality_threshold=0.99 mtz.min_isigi=-5.0 mtz.crystal_name=XTAL".format(
            protScaling._getExtraPath(),
            protScaling._getExtraPath(),
            protExportMtz._getExtraPath(),
            protExportMtz.getObjId(),
            protExportMtz._getLogsPath()
        )
        self.assertCommand(protExportMtz, exportMtzCL, "dials.export")
        exportedmtzset = getattr(protExportMtz, 'exportedFileSet', None)
        self.assertIsNotNone(protExportMtz.exportedFileSet)
        for exportFile in exportedmtzset:
            self.assertTrue(os.path.exists(exportFile.getFilePath()))
        self.assertEqual(protExportMtz.getDatasets(), datasetString)
