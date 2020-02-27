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
import json

import pyworkflow as pw
import pyworkflow.tests as pwtests

import pwed
from pwed.objects import DiffractionImage, SetOfDiffractionImages, DiffractionSpot, SetOfSpots, IndexedSpot, SetOfIndexedSpots
from pwed.protocols import ProtImportDiffractionImages

from dials.protocols import DialsProtFindSpots, DialsProtIndexSpots
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

    # Tests of protocols
    # TODO: Make tests independent of each other
    def test_find_spots(self):
        protImport = self._runImportImages(
            '{TS}/SMV/data/{TI}.img', skipImages=10)
        protFindSpots = self._runFindSpots(protImport.outputDiffractionImages)
        self.assertIsNotNone(protFindSpots.getOutputReflFile())
        outputset = getattr(protFindSpots, 'outputDiffractionSpots', None)
        self.assertIsNotNone(outputset)
        self.assertEqual(outputset.getSpots(), 626)
        self.exampleSetOfSpots = outputset
        # TODO: Add confirmation step that SetOfSpots format and values are correct

    def test_index(self):
        protImport = self._runImportImages(
            '{TS}/SMV/data/{TI}.img', skipImages=10)
        protFindSpots = self._runFindSpots(protImport.outputDiffractionImages)
        protIndex = self._runIndex(
            inputImages=protImport.outputDiffractionImages,
            inputSpots=protFindSpots.outputDiffractionSpots,
            inputStrongPath=os.path.join(
                self.dataPath, 'manual-test'),
            doRefineBravaisSettings=False,
            doReindex=False,
            detectorFixPosition=True,
            detectorFixOrientation=False,
            detectorFixdistance=False,
        )
        self.assertIsNotNone(protIndex.outputIndexedSpots)

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
        with open(modelPath) as mf:
            m = json.load(mf)

        # Replace the path prefix to match with template value
        p = m['imageset'][0]['template']
        m['imageset'][0]['template'] = p.replace(pwed.Config.SCIPION_ED_TESTDATA,
                                                 "$SCIPION_ED_TESTDATA")

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

    """ def test_writeRefl(self):
        spotfile = os.path.join(self.dataPath, 'IO-test', 'strong.refl')
        reflPath = self.getOutputPath('strong.refl')
        readVar = readRefl(spotfile)
        writeRefl(readVar,
                  fn=reflPath)
        with open(reflPath, 'rb') as rp:
            r = rp.read()
        with open(spotfile, 'rb') as sf:
            s = sf.read()
        self.assertEqual(len(r), len(s)) """
