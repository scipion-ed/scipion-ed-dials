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

import pyworkflow.gui.text as text
import pyworkflow.utils as pwutils


def _showHtmlReport(self, reportPath):
    if pwutils.exists(reportPath):
        text._open_cmd(reportPath)
    else:
        self.showInfo('Could not find an html report')
