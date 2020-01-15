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

import pyworkflow as pw
import pyworkflow.tests as pwtests

import pwed
from pwed.objects import DiffractionImage, SetOfDiffractionImages
from pwed.protocols import ProtImportDiffractionImages

from dials.protocols import DialsProtFindSpots
from dials.convert import writeJson


pw.Config.setDomain(pwed)


class TestEdDialsProtocols(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls, writeLocalConfig=True)
        cls.dataPath = os.environ.get('SCIPION_TEST_ED',
                                      '/data/work_software/scipion-ed/')

        if not os.path.exists(cls.dataPath):
            raise Exception("Can not run ED tests, "
                            "SCIPION_TEST_ED variable not defined. ")

    def _runImportImages(self, filesPattern, **kwargs):
        protImport = self.newProtocol(
            ProtImportDiffractionImages,
            filesPath=os.path.join(self.dataPath),
            filesPattern=filesPattern,
            **kwargs)
        self.launchProtocol(protImport)
        return protImport
    
    def _runFindSpots(self, inputImages,**kwargs):
        protFindSpots = self.newProtocol(DialsProtFindSpots,
            inputImages=inputImages,
            **kwargs)
        self.launchProtocol(protFindSpots)
        return protFindSpots

    def test_find_spots(self):

        protImport = self._runImportImages('{TS}/SMV/data/{TI}.img',skipImages=10)
        protFindSpots = self._runFindSpots(protImport.outputDiffractionImages)

        outputset=getattr(protFindSpots,'SetOfSpots',None)
        outputstats=getattr(protFindSpots,'Statistics',None)
        #self.assertIsNotNone(outputset)
        #self.assertIsNotNone(outputstats)
        # TODO: Add confirmation step that SetOfSpots format and values are correct (after defining the set)
        # TODO: Add confirmation step that Statistics format and values are correct (after defining them)

    def test_writeJson(self):
        import json
        template = os.path.join(self.dataPath,'IO-test','imported.expt')
        self.assertIsNotNone(template)
        protImport = self._runImportImages('{TS}/SMV/data/{TI}.img')
        inputImages = protImport.outputDiffractionImages
        modelPath = os.path.join(self.dataPath,'IO-test','testoutput','model.expt')
        if os.path.exists(modelPath):
            os.remove(modelPath)
        model = writeJson(inputImages,fn=modelPath)
        self.assertIsNotNone(model)
        with open(template) as tf:
            t = json.load(tf)
        with open(model) as mf:
            m = json.load(mf)
        self.assertEqual(type(t), type(m))
        if type(t) == dict:
            self.assertEqual(len(t),len(m))
            for key in t:
                self.assertIn(key,m)
                self.assertEqual(t[key],m[key])
        elif type(t) == list:
            self.assertEqual(len(t),len(m))
            for itemA, itemB in zip(t, m):
                self.assertEqual(itemA,itemB)
        else:
            self.assertEqual(t,m)
