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

    def _runImportImages(self, filesPattern):
        protImport = self.newProtocol(
            ProtImportDiffractionImages,
            filesPath=os.path.join(self.dataPath),
            filesPattern=filesPattern)
        self.launchProtocol(protImport)
        return protImport

    def test_find_spots(self):

        protImport = self._runImportImages('{TS}/SMV/data/{TI}.img')

        findSpotsProt = self.newProtocol(DialsProtFindSpots)

        findSpotsProt.inputImages.set(protImport.outputDiffractionImages)
        self.launchProtocol(findSpotsProt)
