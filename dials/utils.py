# **************************************************************************
# *
# * Authors:     V. E.G: Bengtsson (viktor.bengtsson@mmk.su.se) [1]
# *
# * [1] Department of Materials and Environmental Chemistry, Stockholm University
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
import json
import os.path as p
import subprocess
from packaging import version
import pyworkflow.gui.text as text
import pyworkflow.utils as pwutils

from dials.objects import *


def _showHtmlReport(reportPath):
    if pwutils.exists(reportPath):
        text._open_cmd(reportPath)


def getModelDataPath(exptFile):
    dataPaths = []
    with open(exptFile) as json_file:
        data = json.load(json_file)
        for imageset in data['imageset']:
            template = imageset['template']
            dataPaths.append(p.dirname(template))
    return dataPaths


def getDatasets(modelFile):
    datasets = getModelDataPath(modelFile)
    if len(datasets) >= 1:
        newlineDatasets = '\n'.join(f'{item}' for item in datasets)
        return f"Source of data:\n{newlineDatasets}"
    else:
        return ''


def existsPath(path):
    return p.exists(path)


def joinPath(*args):
    return p.join(*args)


def readLog(logfile, start, stop, flush=None):
    # based on https://stackoverflow.com/a/18865133
    contentList = []
    content = ''
    with open(logfile) as infile:
        append = False
        for line in infile:
            if start in line.strip():
                append = True
                contentList.append(line)
                continue
            elif stop in line.strip():
                append = False
                continue
            elif flush is not None and flush in line.strip():
                append = False
                contentList.clear()
                continue
            elif append:
                contentList.append(line)

    newlineList = ''.join(f'{item}' for item in contentList)
    content = f"{newlineList}"

    return content


def verifyPathExistence(*requiredPaths):
    for reqPath in requiredPaths:
        if not existsPath(reqPath):
            raise MissingPathException
        else:
            continue


def getDialsVersion():
    # Run subprocess.run and save output object
    output = subprocess.run("dials.version", capture_output=True)
    outputString = output.stdout.decode("utf-8").strip()
    # Get the first line of the input with the full version name
    dialsVersion = outputString.split("\n")[0].strip("DIALS")
    # Get only the semantic versioning part of the version number
    versionNumber = dialsVersion.split("-")[0]
    return versionNumber


def isMinDialsVersion(minversion):
    return version.parse(getDialsVersion()) >= version.parse(minversion)
