# **************************************************************************
# *
# * Authors:     Viktor E. G. Bengtsson (viktor.bengtsson@mmk.su.se)   [1]
# *
# * [1] MMK, Stockholm University
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

import dials.utils as dutils


def _getModel(protocol):
    try:
        modelFile = protocol.getModelFile()
        if dutils.existsPath(modelFile):
            return modelFile
    except AttributeError:
        pass
    try:
        modelFile = protocol.getOutputModelFile()
        if dutils.existsPath(modelFile):
            return modelFile
    except AttributeError:
        pass
    try:
        modelFile = protocol.getInputModelFile()
        if dutils.existsPath(modelFile):
            return modelFile
    except AttributeError:
        pass
    return None


def _getRefls(protocol):
    try:
        reflFile = protocol.getReflFile()
        if dutils.existsPath(reflFile):
            return reflFile
    except AttributeError:
        pass
    try:
        reflFile = protocol.getOutputReflFile()
        if dutils.existsPath(reflFile):
            return reflFile
    except AttributeError:
        pass
    try:
        reflFile = protocol.getInputReflFile()
        if dutils.existsPath(reflFile):
            return reflFile
    except AttributeError:
        pass
    return None


def _getHtml(protocol):
    try:
        htmlFile = protocol.getOutputHtmlFile()
        if dutils.existsPath(htmlFile):
            return htmlFile
    except AttributeError:
        pass
    return None
