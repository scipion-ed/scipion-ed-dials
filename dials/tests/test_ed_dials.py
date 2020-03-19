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
import io
import json

import msgpack

import pyworkflow as pw
import pyworkflow.tests as pwtests

import pwed
from pwed.objects import DiffractionImage, SetOfDiffractionImages, DiffractionSpot, SetOfSpots, IndexedSpot, SetOfIndexedSpots, RefinedSpot, SetOfRefinedSpots, IntegratedSpot, SetOfIntegratedSpots
from pwed.protocols import ProtImportDiffractionImages

from dials.protocols import DialsProtImportDiffractionImages, DialsProtFindSpots, DialsProtIndexSpots, DialsProtRefineSpots, DialsProtIntegrateSpots, DialsProtExport
from dials.convert import writeJson, readRefl, writeRefl


pw.Config.setDomain(pwed)


class TestEdDialsProtocols(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls, writeLocalConfig=True)
        cls.dataPath = os.path.join(pwed.Config.SCIPION_ED_TESTDATA,
                                    '190503')

        if not os.path.exists(cls.dataPath):
            raise Exception("Can not run DIALS tests, missing file:\n  %s"
                            % cls.dataPath)

    # Functions for running protocols
    def _runImportImages(self, filesPattern, **kwargs):
        protImport = self.newProtocol(
            ProtImportDiffractionImages,
            filesPath=os.path.join(self.dataPath),
            filesPattern=filesPattern,
            **kwargs)
        self.launchProtocol(protImport)
        return protImport

    def _runDialsImportImages(self, filesPattern, **kwargs):
        protImport = self.newProtocol(
            DialsProtImportDiffractionImages,
            filesPath=os.path.join(self.dataPath),
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

    def _runExport(self, inputSet, **kwargs):
        protExport = self.newProtocol(DialsProtExport,
                                      inputSet=inputSet,
                                      **kwargs)
        self.launchProtocol(protExport)
        return protExport

    # Helper functions
    def assertSameModel(self, model1, model2):
        with io.open(model1) as m1:
            m1m = [s[s.find('extra/'):] for s in list(m1)]
        with io.open(model2) as m2:
            m2m = [s[s.find('extra/'):] for s in list(m2)]
        return self.assertEqual(m1m, m2m)

    def assertSameRefl(self, refl1, refl2):
        with io.open(refl1, 'rb') as r1:
            with io.open(refl2, 'rb') as r2:
                return self.assertEqual(r1, r2)

    def updateImagePath(self, modelPath):
        with open(modelPath) as mf:
            m = json.load(mf)
            p = m['imageset'][0]['template']
            m['imageset'][0]['template'] = p.replace(pwed.Config.SCIPION_ED_TESTDATA,
                                                     "$SCIPION_ED_TESTDATA")
        return m

    def getReferenceFile(self, reffile=None, lyso=False):
        if reffile != None:
            if lyso:
                reference = os.path.join(
                    self.dataPath, 'lyso', 'reference', reffile)
            else:
                reference = os.path.join(self.dataPath, 'manual-test', reffile)
            try:
                self.updateImagePath(reference)
            except UnicodeDecodeError:
                pass
            return reference

    # Tests of protocols
    def test_find_spots(self):
        self.skipTest("Skip this test until it can be made independent")
        protImport = self._runDialsImportImages(
            '{TS}/SMV/data/{TI}.img', skipImages=10, rotationAxis='-0.6204,-0.7843,0')
        protFindSpots = self._runFindSpots(protImport.outputDiffractionImages)
        self.assertIsNotNone(protFindSpots.getOutputReflFile())
        outputset = getattr(protFindSpots, 'outputDiffractionSpots', None)
        self.assertIsNotNone(outputset)
        self.assertEqual(outputset.getSpots(), 626)
        self.assertIsNotNone(outputset.getDialsModel())
        self.assertIsNotNone(outputset.getDialsRefl())
        with self.subTest(msg="Testing reflections in SetOfDiffractionSpots"):
            self.skipTest("Need to fix errors from paths")
            self.assertSameRefl(outputset.getDialsRefl(),
                                self.getReferenceFile('strong.refl'))
        # TODO: Add confirmation step that SetOfSpots format and values are correct

    def test_standard_pipeline(self):
        # Run import
        protImport = self._runDialsImportImages(
            '{TS}/SMV/data/{TI}.img', skipImages=10, rotationAxis='-0.6204,-0.7843,0')
        self.assertIsNotNone(protImport.getOutputModelFile())
        importedset = getattr(protImport, 'outputDiffractionImages', None)
        self.assertIsNotNone(importedset)
        self.assertIsNotNone(importedset.getDialsModel())
        with self.subTest(msg="Testing model in SetOfDiffractionImages"):
            self.assertSameModel(importedset.getDialsModel(),
                                 self.getReferenceFile('imported.expt'))
        # Run find spots
        protFindSpots = self._runFindSpots(protImport.outputDiffractionImages)
        foundspotset = getattr(protFindSpots, 'outputDiffractionSpots', None)
        self.assertIsNotNone(foundspotset)
        self.assertEqual(foundspotset.getSpots(), 626)
        self.assertIsNotNone(foundspotset.getDialsModel())
        with self.subTest(msg="Testing model in SetOfDiffractionSpots"):
            self.assertSameModel(foundspotset.getDialsModel(),
                                 self.getReferenceFile('imported.expt'))
        self.assertIsNotNone(foundspotset.getDialsRefl())
        with self.subTest(msg="Testing reflections in SetOfDiffractionSpots"):
            self.skipTest("Need to fix errors from paths")
            self.assertSameRefl(foundspotset.getDialsRefl(),
                                self.getReferenceFile('strong.refl'))
        # Run indexing
        with self.subTest(msg="Do not refine Bravais settings"):
            protIndex = self._runIndex(
                inputImages=protImport.outputDiffractionImages,
                inputSpots=protFindSpots.outputDiffractionSpots,
                doRefineBravaisSettings=False,
                doReindex=False,
                detectorFixPosition=True,
                detectorFixOrientation=False,
                detectorFixdistance=False,
            )
            indexedset = getattr(protIndex, 'outputIndexedSpots', None)
            self.assertIsNotNone(protIndex.outputIndexedSpots)
            self.assertIsNotNone(indexedset.getDialsModel())
            self.assertIsNotNone(indexedset.getDialsRefl())
            with self.subTest(msg="Testing model in SetOfIndexedSpots"):
                self.assertSameModel(indexedset.getDialsModel(),
                                     self.getReferenceFile('indexed.expt'))
            with self.subTest(msg="Testing reflections in SetOfIndexedSpots"):
                self.skipTest("Need to fix errors from paths")
                self.assertSameRefl(indexedset.getDialsRefl(),
                                    self.getReferenceFile('indexed.refl'))

        with self.subTest(msg="Refine Bravais settings and reindex"):
            protIndex = self._runIndex(
                objLabel="index and reindex",
                inputImages=protImport.outputDiffractionImages,
                inputSpots=protFindSpots.outputDiffractionSpots,
                doRefineBravaisSettings=True,
                doReindex=True,
                doReindexReflections=True,
                detectorFixPosition=True,
                detectorFixOrientation=False,
                detectorFixdistance=False,
            )
            indexedset = getattr(protIndex, 'outputIndexedSpots', None)
            self.assertIsNotNone(protIndex.outputIndexedSpots)
            self.assertIsNotNone(indexedset.getDialsModel())
            self.assertIsNotNone(indexedset.getDialsRefl())
            with self.subTest(msg="Testing model after reindexing"):
                self.assertSameModel(indexedset.getDialsModel(),
                                     self.getReferenceFile('bravais_setting_12.expt'))
            with self.subTest(msg="Testing reflections after reindexing"):
                self.skipTest("Need to fix errors from paths")
                self.assertSameRefl(indexedset.getDialsRefl(),
                                    self.getReferenceFile('reindexed.refl'))

        # Run refinement
        with self.subTest(msg="Static refinement"):
            # self.skipTest("Not implemented")
            protRefine = self._runRefine(
                objLabel="Static refinement",
                inputSet=protIndex.outputIndexedSpots,
                scanVarying=False,
            )
            refinedset = getattr(protRefine, 'outputRefinedSpots', None)
            self.assertIsNotNone(protRefine.outputRefinedSpots)
            self.assertIsNotNone(refinedset.getDialsModel())
            self.assertIsNotNone(refinedset.getDialsRefl())
            with self.subTest(msg="Testing model after refinement"):
                self.assertSameModel(refinedset.getDialsModel(),
                                     self.getReferenceFile('refined.expt'))
            with self.subTest(msg="Testing reflections after refinement"):
                self.skipTest("Need to fix errors from paths")
                self.assertSameRefl(refinedset.getDialsRefl(),
                                    self.getReferenceFile('refined.refl'))

        # Run scan-varying refinement
        with self.subTest(msg="Scan-varying refinement"):
            # self.skipTest("Not implemented")
            protSvRefine = self._runRefine(
                objLabel="Scan-varying refinement",
                inputSet=protRefine.outputRefinedSpots,
                scanVarying=True,
                beamFixAll=False,
                beamFixInSpindlePlane=False,
                beamFixOutSpindlePlane=False,
                beamFixWavelength=True,
                beamForceStatic=False,
                detectorFixAll=True,
            )
            svrefinedset = getattr(protSvRefine, 'outputRefinedSpots', None)
            self.assertIsNotNone(protSvRefine.outputRefinedSpots)
            self.assertIsNotNone(svrefinedset.getDialsModel())
            self.assertIsNotNone(svrefinedset.getDialsRefl())
            with self.subTest(msg="Testing model after scan-varying refinement"):
                self.assertSameModel(svrefinedset.getDialsModel(),
                                     self.getReferenceFile('sv_refined.expt'))
            with self.subTest(msg="Testing reflections after scan-varying refinement"):
                self.skipTest("Need to fix errors from paths")
                self.assertSameRefl(svrefinedset.getDialsRefl(),
                                    self.getReferenceFile('sv_refined.refl'))

        # Run integration
        with self.subTest(msg="Integration"):
            # self.skipTest("Not implemented yet")
            protIntegrate = self._runIntegrate(
                inputSet=protSvRefine.outputRefinedSpots,
                nproc=8,
            )
            integratedset = getattr(
                protIntegrate, 'outputIntegratedSpots', None)
            self.assertIsNotNone(protIntegrate.outputIntegratedSpots)
            self.assertIsNotNone(integratedset.getDialsModel())
            self.assertIsNotNone(integratedset.getDialsRefl())
            with self.subTest(msg="Testing model after integration"):
                self.assertSameModel(integratedset.getDialsModel(),
                                     self.getReferenceFile('integrated.expt'))
            with self.subTest(msg="Testing reflections after integration"):
                self.skipTest("Need to fix errors from paths")
                self.assertSameRefl(integratedset.getDialsRefl(),
                                    self.getReferenceFile('integrated.refl'))

        MTZ = 0
        SADABS = 1
        NXS = 2
        MMCIF = 3
        MOSFLM = 4
        XDS = 5
        XDS_ASCII = 6
        JSON = 7

        with self.subTest(msg="Export mtz"):
            protMtzExport = self._runExport(
                inputSet=protIntegrate.outputIntegratedSpots,
                exportFormat=MTZ
            )
        with self.subTest(msg="Export xds_ascii"):
            protXdsExport = self._runExport(
                inputSet=protIntegrate.outputIntegratedSpots,
                exportFormat=XDS_ASCII
            )

    def test_lyso_pipeline(self):
        # Run import
        protImport = self._runDialsImportImages(
            'lyso/{TS}/SMV/data/{TI}.img',
            objLabel="Lyso: dials import",
            rotationAxis='-0.6204,-0.7843,0',
            commandLineInput='geometry.detector.distance=2562')
        self.assertIsNotNone(protImport.getOutputModelFile())
        importedset = getattr(protImport, 'outputDiffractionImages', None)
        self.assertIsNotNone(importedset)
        self.assertIsNotNone(importedset.getDialsModel())
        with self.subTest(msg="Testing model in SetOfDiffractionImages"):
            self.assertSameModel(importedset.getDialsModel(),
                                 self.getReferenceFile('imported.expt', lyso=True))
        # Run find spots
        protFindSpots = self._runFindSpots(protImport.outputDiffractionImages,
                                           objLabel="Lyso: find spots",
                                           dMin=3.5,
                                           untrustedAreas=True,
                                           untrustedRectangle_1='0,516,255,261',
                                           untrustedRectangle_2='255,261,0,516',
                                           )
        foundspotset = getattr(protFindSpots, 'outputDiffractionSpots', None)
        self.assertIsNotNone(foundspotset)
        self.assertEqual(foundspotset.getSpots(), 785)
        self.assertIsNotNone(foundspotset.getDialsModel())
        with self.subTest(msg="Testing model in SetOfDiffractionSpots"):
            self.assertSameModel(foundspotset.getDialsModel(),
                                 self.getReferenceFile('imported.expt', lyso=True))
        self.assertIsNotNone(foundspotset.getDialsRefl())
        with self.subTest(msg="Testing reflections in SetOfDiffractionSpots"):
            self.skipTest("Need to fix errors from paths")
            self.assertSameRefl(foundspotset.getDialsRefl(),
                                self.getReferenceFile('strong.refl', lyso=True))
        # Run indexing
        with self.subTest(msg="Refine Bravais settings and reindex"):
            protIndex = self._runIndex(
                objLabel="Lyso: index and reindex",
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
            )
            indexedset = getattr(protIndex, 'outputIndexedSpots', None)
            self.assertIsNotNone(protIndex.outputIndexedSpots)
            self.assertIsNotNone(indexedset.getDialsModel())
            self.assertIsNotNone(indexedset.getDialsRefl())
            with self.subTest(msg="Testing model after reindexing"):
                self.assertSameModel(indexedset.getDialsModel(),
                                     self.getReferenceFile('bravais_setting_9.expt', lyso=True))
            with self.subTest(msg="Testing reflections after reindexing"):
                self.skipTest("Need to fix errors from paths")
                self.assertSameRefl(indexedset.getDialsRefl(),
                                    self.getReferenceFile('reindexed.refl', lyso=True))

        # Run refinement
        with self.subTest(msg="Static refinement"):
            # self.skipTest("Not implemented")
            protRefine = self._runRefine(
                objLabel="Lyso: static refinement",
                inputSet=protIndex.outputIndexedSpots,
                scanVarying=False,
                detectorFixAll=True,
            )
            refinedset = getattr(protRefine, 'outputRefinedSpots', None)
            self.assertIsNotNone(protRefine.outputRefinedSpots)
            self.assertIsNotNone(refinedset.getDialsModel())
            self.assertIsNotNone(refinedset.getDialsRefl())
            with self.subTest(msg="Testing model after refinement"):
                self.assertSameModel(refinedset.getDialsModel(),
                                     self.getReferenceFile('refined.expt', lyso=True))
            with self.subTest(msg="Testing reflections after refinement"):
                self.skipTest("Need to fix errors from paths")
                self.assertSameRefl(refinedset.getDialsRefl(),
                                    self.getReferenceFile('refined.refl', lyso=True))

        # Run scan-varying refinement
        with self.subTest(msg="Scan-varying refinement"):
            # self.skipTest("Not implemented")
            protSvRefine = self._runRefine(
                objLabel="Lyso: scan-varying refinement",
                inputSet=protRefine.outputRefinedSpots,
                scanVarying=True,
                beamFixAll=False,
                beamFixInSpindlePlane=False,
                beamFixOutSpindlePlane=False,
                beamFixWavelength=True,
                beamForceStatic=False,
                detectorFixAll=True,
            )
            svrefinedset = getattr(protSvRefine, 'outputRefinedSpots', None)
            self.assertIsNotNone(protSvRefine.outputRefinedSpots)
            self.assertIsNotNone(svrefinedset.getDialsModel())
            self.assertIsNotNone(svrefinedset.getDialsRefl())
            with self.subTest(msg="Testing model after scan-varying refinement"):
                self.assertSameModel(svrefinedset.getDialsModel(),
                                     self.getReferenceFile('sv_refined.expt', lyso=True))
            with self.subTest(msg="Testing reflections after scan-varying refinement"):
                self.skipTest("Need to fix errors from paths")
                self.assertSameRefl(svrefinedset.getDialsRefl(),
                                    self.getReferenceFile('sv_refined.refl', lyso=True))

        # Run integration
        with self.subTest(msg="Integration"):
            # self.skipTest("Not implemented yet")
            protIntegrate = self._runIntegrate(
                objLabel="Lyso: integration",
                inputSet=protSvRefine.outputRefinedSpots,
                nproc=8,
                commandLineInput='prediction.d_min=3.5',
            )
            integratedset = getattr(
                protIntegrate, 'outputIntegratedSpots', None)
            self.assertIsNotNone(protIntegrate.outputIntegratedSpots)
            self.assertIsNotNone(integratedset.getDialsModel())
            self.assertIsNotNone(integratedset.getDialsRefl())
            with self.subTest(msg="Testing model after integration"):
                self.assertSameModel(integratedset.getDialsModel(),
                                     self.getReferenceFile('integrated.expt', lyso=True))
            with self.subTest(msg="Testing reflections after integration"):
                self.skipTest("Need to fix errors from paths")
                self.assertSameRefl(integratedset.getDialsRefl(),
                                    self.getReferenceFile('integrated.refl', lyso=True))

        MTZ = 0
        SADABS = 1
        NXS = 2
        MMCIF = 3
        MOSFLM = 4
        XDS = 5
        XDS_ASCII = 6
        JSON = 7

        with self.subTest(msg="Export mtz"):
            protMtzExport = self._runExport(
                objLabel="Lyso: export mtz",
                inputSet=protIntegrate.outputIntegratedSpots,
                exportFormat=MTZ
            )
        with self.subTest(msg="Export xds_ascii"):
            protXdsExport = self._runExport(
                objLabel="Lyso: export XDS_ASCII",
                inputSet=protIntegrate.outputIntegratedSpots,
                exportFormat=XDS_ASCII
            )

    # Test of I/O utilities
    def test_writeJson(self):

        template = os.path.join(self.dataPath, 'IO-test', 'imported.expt')
        self.assertIsNotNone(template)
        protImport = self._runImportImages('{TS}/SMV/data/{TI}.img',
                                           rotationAxis='-0.6204,-0.7843,0')
        inputImages = protImport.outputDiffractionImages
        modelPath = self.getOutputPath('model.expt')
        writeJson(inputImages, fn=modelPath)
        with open(template) as tf:
            t = json.load(tf)
        m = self.updateImagePath(modelPath)
        """ with open(modelPath) as mf:
            m = json.load(mf)

        # Replace the path prefix to match with template value
        p = m['imageset'][0]['template']
        m['imageset'][0]['template'] = p.replace(pwed.Config.SCIPION_ED_TESTDATA,
                                                 "$SCIPION_ED_TESTDATA") """

        self.assertEqual(type(t), type(m))
        if type(t) == dict:
            self.assertEqual(len(t), len(m))
            for key in t:
                self.assertIn(key, m)
                self.assertEqual(t[key], m[key])
        elif type(t) == list:
            self.assertEqual(len(t), len(m))
            for itemA, itemB in zip(t, m):
                self.assertEqual(itemA, itemB)
        else:
            self.assertEqual(t, m)

    def test_readRefl(self):
        spotfile = os.path.join(self.dataPath, 'IO-test', 'strong.refl')
        result = readRefl(spotfile)
        self.assertIsNotNone(result)
