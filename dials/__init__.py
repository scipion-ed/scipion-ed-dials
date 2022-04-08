# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Viktor E. G. Bengtsson (viktor.bengtsson@mmk.su.se)   [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] MMK, Stockholm University
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
"""
This modules contains classes related with ED
"""

import pyworkflow.plugin as pwplugin
import pwed
from .constants import *

# Epoch indicates compatible main Scipion version
# major.minor.micro versioning starting with 1.0.0 in the new epoch
__version__ = '3!1.0.1'
_logo = "logo.jpeg"


class Plugin(pwplugin.Plugin):
    pass


pwed.Domain.registerPlugin(__name__)
