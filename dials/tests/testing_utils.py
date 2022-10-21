# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se)
# *              V E.G. Bengtsson       (viktor.e.g.bengtsson@gmail.com)
# *
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
from typing import Union

import pyworkflow.tests as pwtests
from pwed.protocols import ProtImportDiffractionImages

from dials.protocols import (
    DialsProtBase,
    DialsProtExport,
    DialsProtFindSpots,
    DialsProtImportDiffractionImages,
    DialsProtIndexSpots,
    DialsProtIntegrateSpots,
    DialsProtMerge,
    DialsProtRefineSpots,
    DialsProtScaling,
    DialsProtSymmetry,
    HtmlBase,
)


def makePathList(
    scanRange,
    digits=3,
    location=None,
    pattern="SMV/data/00{TI}.img",
    wildcard="{TI}",
) -> list:
    first, last = scanRange.split(",")
    numbers = range(int(first), int(last) + 1)
    replacements = [str(num).zfill(digits) for num in numbers]
    images = [
        pattern.replace(wildcard, replacement) for replacement in replacements
    ]
    imageList = [os.path.join(location, i) for i in images]
    return imageList


class ProtocolRunner(pwtests.BaseTest):
    """
    A base class to define the functions used
    for running protocols during testing.
    """

    def _runImportImages(
        self, filesPath, filesPattern="SMV/data/00{TI}.img", **kwargs
    ) -> ProtImportDiffractionImages:
        protImport = self.newProtocol(
            ProtImportDiffractionImages,
            filesPath=filesPath,
            filesPattern=filesPattern,
            **kwargs,
        )
        self.launchProtocol(protImport)
        return protImport

    def _runDialsImportImages(
        self, filesPath, filesPattern="SMV/data/00{TI}.img", **kwargs
    ) -> DialsProtImportDiffractionImages:
        protImport = self.newProtocol(
            DialsProtImportDiffractionImages,
            filesPath=filesPath,
            filesPattern=filesPattern,
            **kwargs,
        )
        self.launchProtocol(protImport)
        return protImport

    def _runFindSpots(self, inputImages, **kwargs) -> DialsProtFindSpots:
        protFindSpots = self.newProtocol(
            DialsProtFindSpots, inputImages=inputImages, **kwargs
        )
        self.launchProtocol(protFindSpots)
        return protFindSpots

    def _runIndex(
        self, inputImages, inputSpots, **kwargs
    ) -> DialsProtIndexSpots:
        protIndex = self.newProtocol(
            DialsProtIndexSpots,
            inputImages=inputImages,
            inputSpots=inputSpots,
            **kwargs,
        )
        self.launchProtocol(protIndex)
        return protIndex

    def _runRefine(self, inputSet, **kwargs) -> DialsProtRefineSpots:
        protRefine = self.newProtocol(
            DialsProtRefineSpots, inputSet=inputSet, **kwargs
        )
        self.launchProtocol(protRefine)
        return protRefine

    def _runIntegrate(self, inputSet, **kwargs) -> DialsProtIntegrateSpots:
        protIntegrate = self.newProtocol(
            DialsProtIntegrateSpots, inputSet=inputSet, **kwargs
        )
        self.launchProtocol(protIntegrate)
        return protIntegrate

    def _runSymmetry(self, inputSet, **kwargs) -> DialsProtSymmetry:
        protSymmetry = self.newProtocol(
            DialsProtSymmetry, inputSet=inputSet, **kwargs
        )
        self.launchProtocol(protSymmetry)
        return protSymmetry

    def _runScaling(self, inputSets, **kwargs) -> DialsProtScaling:
        protScaling = self.newProtocol(
            DialsProtScaling, inputSets=inputSets, **kwargs
        )
        self.launchProtocol(protScaling)
        return protScaling

    def _runMerging(self, inputSet, **kwargs) -> DialsProtMerge:
        protMerging = self.newProtocol(
            DialsProtMerge, inputSet=inputSet, **kwargs
        )
        self.launchProtocol(protMerging)
        return protMerging

    def _runExport(self, inputSet, **kwargs) -> DialsProtExport:
        protExport = self.newProtocol(
            DialsProtExport, inputSet=inputSet, **kwargs
        )
        self.launchProtocol(protExport)
        return protExport


class HelperCollection(pwtests.BaseTest):
    def checkLogDataset(
        self, protocol: DialsProtBase, dataset, logOutput=""
    ) -> None:
        datasetString = (
            f"Source of data:\n" f"{os.path.join(dataset, 'SMV/data')}"
        )
        self.assertEqual(protocol.getDatasets(), datasetString)
        outputCompare = protocol.getLogOutput().split("\n")[0]
        self.assertEqual(outputCompare.strip(), logOutput.strip())

    def _getCommand(
        self,
        protocol: Union[DialsProtBase, DialsProtRefineSpots, HtmlBase],
        program=None,
    ) -> str:
        try:
            return protocol._prepareCommandline(program)
        except AttributeError:
            pass
        try:
            return protocol._prepCommandline(program)
        except AttributeError:
            pass

    def assertCommand(
        self,
        protocol: DialsProtBase,
        commandString: str,
        program: Union[str, None] = None,
    ):
        CL = self._getCommand(protocol, program)
        self.assertEqual(CL, commandString)

    def assertFileExists(self, file):
        self.assertTrue(os.path.exists(file))

    def comparePhils(self, goodPhil="restraints.phil", testPhil=None):
        self.assertIsNotNone(goodPhil)
        self.assertIsNotNone(testPhil)
        with open(goodPhil) as f1:
            with open(testPhil) as f2:
                self.assertEqual(f1.read(), f2.read())
