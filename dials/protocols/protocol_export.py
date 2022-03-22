# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              V. E.G: Bengtsson (viktor.bengtsson@mmk.su.se) [2]
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

import pyworkflow.protocol as pwprot
import dials.convert as dconv
import dials.utils as dutils

from pwed.objects import ExportFile
from pwed.protocols import EdProtExport
from dials.protocols import DialsProtBase, PhilBase, CliBase
from dials.constants import exportFormatChoices
from dials.constants import *


class DialsProtExport(EdProtExport, DialsProtBase):
    """ Protocol for exporting results using Dials
    """

    _label = 'export'

    # -------------------------- DEFINE param functions -----------------------

    def _defineParams(self, form):
        # EdProtIndexSpots._defineParams(self, form)

        # The start of the actually relevant part.
        form.addSection(label='Input')

        form.addParam('inputSet', pwprot.PointerParam,
                      pointerClass='SetOfIndexedSpots',
                      label="Integrated spots to export",
                      help="")

        # Help messages are copied from the DIALS documentation at
        # https://dials.github.io/documentation/programs/dials_export.html

        form.addParam('exportFormat', pwprot.EnumParam,
                      label='Which format do you want to export?',
                      choices=exportFormatChoices,
                      default=MTZ,
                      display=pwprot.EnumParam.DISPLAY_HLIST,
                      help="The output file format. Please note that "
                      "XDS_ASCII is incompatible with scaled data."
                      )

        group = form.addGroup('mtz', condition=f"exportFormat=={MTZ}")

        group.addParam('mtzCombinePartials', pwprot.BooleanParam,
                       label='Combine partial reflections?', default=True,
                       help="Combine partials that have the same partial "
                       "id into one reflection, with an updated partiality "
                       "given by the sum of the individual partialities.",
                       )

        group.addParam('mtzPartialityThreshold', pwprot.FloatParam,
                       label='Partiality threshold',
                       default=0.99, allowsNull=True,
                       condition="mtzCombinePartials==True",
                       help="All reflections with partiality values above "
                       "the partiality threshold will be retained. This is "
                       "done after any combination of partials if applicable.",
                       )

        group.addParam('mtzMinIsigi', pwprot.FloatParam,
                       label='min I/Sig(I)',
                       default=-5, allowsNull=True,
                       help="Exclude reflections with unfeasible "
                       "values of I/Sig(I)",
                       )

        group.addParam('mtzForceStaticModel', pwprot.BooleanParam,
                       label='Force static model?', default=False,
                       help="Force program to use static model even if "
                       "scan varying is present",
                       )

        group.addParam('mtzFilter_ice', pwprot.BooleanParam, default=False,
                       label='Filter ice?',
                       help="Filter out reflections at typical ice ring "
                       "resolutions before max_cell estimation.")

        group.addParam('mtzDMin', pwprot.FloatParam,
                       label='d_min',
                       default=None, allowsNull=True,
                       help="Filter out reflections with d-spacing "
                       "below d_min",
                       )

        group.addParam('mtzHklout', pwprot.StringParam,
                       label='Output file',
                       default="",
                       help="The output MTZ filename, defaults to "
                       "integrated_<jobID>.mtz",
                       )

        group.addParam('mtzCrystalName', pwprot.StringParam,
                       label='Crystal name',
                       default='XTAL',
                       help="The name of the crystal "
                       "for the mtz file metadata",
                       )

        group = form.addGroup(
            'sadabs', condition=f"exportFormat=={SADABS}")

        group.addParam('sadabsHklout', pwprot.StringParam,
                       label='Output filename',
                       default='integrated.sad',
                       help="The output raw sadabs file",
                       )

        group.addParam('sadabsRun', pwprot.IntParam,
                       label='Batch or run number',
                       default=1, allowsNull=True,
                       help="Batch number / run number for output file",
                       )

        group.addParam('sadabsPredict', pwprot.BooleanParam,
                       label='Predict from static model', default=False,
                       help="Compute centroids with static model, "
                       "not observations",
                       )

        group = form.addGroup(
            'Nexus', condition=f"exportFormat=={NXS}")
        group.addParam('nxsHklout', pwprot.StringParam,
                       label='Output filename',
                       default='integrated.nxs',
                       help="The output Nexus file",
                       )

        group.addParam('nxsInstrumentName', pwprot.StringParam,
                       label='Instrument name',
                       default='Unknown',
                       help="Name of the instrument/beamline",
                       )

        group.addParam('nxsInstrumentShortName', pwprot.StringParam,
                       label='Short instrument name',
                       default='Unknown',
                       help="Short name for instrument/beamline, "
                       "perhaps the acronym",
                       )

        group.addParam('nxsSourceName', pwprot.StringParam,
                       label='Source name',
                       default='Unknown',
                       help="Name of the source/facility",
                       )

        group.addParam('nxsSourceShortName', pwprot.StringParam,
                       label='Short source name',
                       default='Unknown',
                       help="Short name for source, perhaps the acronym",
                       )

        group = form.addGroup(
            'mmcif', condition=f"exportFormat=={MMCIF}")

        group.addParam('mmcifHklout', pwprot.StringParam,
                       label='Output name',
                       default='',
                       help="The output CIF file, defaults to "
                       "<jobID>_integrated.cif.",
                       )

        group = form.addGroup(
            'XDS_ASCII', condition=f"exportFormat=={XDS_ASCII}")

        group.addParam('xdsAsciiHklout', pwprot.StringParam,
                       label='Output name',
                       default='DIALS.HKL',
                       help="The output raw hkl file",
                       )

        group = form.addGroup(
            'json', condition=f"exportFormat=={JSON}")
        group.addParam('jsonFilename', pwprot.StringParam,
                       label='Filename',
                       default='rlp.json',
                       help="",
                       )

        group.addParam('jsonCompact', pwprot.BooleanParam,
                       label='Compact?', default=True,
                       )
        group.addParam('jsonNDigits', pwprot.IntParam,
                       label='Number of decimal places?',
                       default=6, allowsNull=True,
                       help="Number of decimal places to be used for "
                       "representing the reciprocal lattice points.",
                       )

        if dutils.isMinDialsVersion("3.9.0"):
            group = form.addGroup(
                "shelx", condition=f"exportFormat=={SHELX}")
            group.addParam("shelxName", pwprot.StringParam,
                           label="Shelx base name",
                           default="dials",
                           )
            group.addParam("shelxScale", pwprot.BooleanParam,
                           label="Scale reflections?",
                           help=("Scale reflections to maximise output "
                                 "precision in SHELX 8.2f format"),
                           default=True,
                           )
            group.addParam("shelxScaleMin", pwprot.FloatParam,
                           condition="shelxScale",
                           label="Minimum intensity value after scaling",
                           default=-9999.0
                           )
            group.addParam("shelxScaleMax", pwprot.FloatParam,
                           condition="shelxScale",
                           label="Maximum intensity value after scaling",
                           default=9999.0
                           )
            group = form.addGroup(
                "pets", condition=f"exportFormat=={PETS}")
            group.addParam("petsPrefix", pwprot.StringParam,
                           label="Filename prefix",
                           default="dials_dyn",
                           help=("The prefix for output files, where the "
                                 "default will produce dials_dyn.cif_pets")
                           )
            group.addParam("petsId", pwprot.IntParam,
                           label=("ID of experiment to export from"
                                  " multi-experiment list:"),
                           default=None,
                           allowsNull=True,
                           expertLevel=pwprot.LEVEL_ADVANCED
                           )
            group.addParam("petsPartialityCutoff", pwprot.FloatParam,
                           label="Partiality cutoff",
                           default=0.99,
                           allowsNull=True,
                           help=("Cutoff for determining which reflections are"
                                 " deemed to be fully recorded")
                           )
            group.addParam("petsFlagFiltering", pwprot.BooleanParam,
                           label="Flag filtering",
                           default=False,
                           expertLevel=pwprot.LEVEL_ADVANCED,
                           help=("If true, keep only the reflections where the"
                                 " relevant `integrated` flag is set (either "
                                 "`integrated_sum` or `integrated_prf`). This"
                                 " seems to be quite restrictive compared to"
                                 " PETS, so is not set by default.")
                           )
            group.addParam("petsVFExcitationErrorCutoff", pwprot.FloatParam,
                           label="Excitation error cutoff",
                           defualt=0.04,
                           allowsNull=True,
                           help=("Excitation error cutoff determining which "
                                 "reflections are included in virtual frames")
                           )
            group.addParam("petsVFNMerged", pwprot.IntParam,
                           label=("Number of frames to merge in a virtual"
                                  " frame"),
                           default=1,
                           allowsNull=True,
                           )
            group.addParam("petsVFStep", pwprot.IntParam,
                           label=("Steps between frames"),
                           default=1,
                           allowsNull=True,
                           )

        # Allow an easy way to import a phil file with parameters
        PhilBase._definePhilParams(self, form)

        # Allow adding anything else with command line syntax
        CliBase._defineCliParams(self, form)

    # -------------------------- INSERT functions ----------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(
            'convertInputStep', self.inputSet.getObjId())
        self._insertFunctionStep('exportStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -----------------------------
    def convertInputStep(self, inputSpotId):
        inputSet = self.inputSet.get()
        self.info(f"Number of spots: {inputSet.getSize()}")
        # Write new model and/or reflection file if no was supplied from set
        if self._checkWriteModel():
            dconv.writeJson(inputSet, self.getInputModelFile())
        if self._checkWriteRefl():
            dconv.writeRefl(inputSet, self.getInputReflFile())

    def exportStep(self):
        program = 'dials.export'
        arguments = self._prepareCommandline(program)
        try:
            self.runJob(program, arguments)
        except:
            self.info(self.getError())

    def createOutputStep(self):
        outputSet = self._createSetOfExportFiles()
        eFile = ExportFile()
        eFile.setFilePath(self.getExport())
        eFile.setFileType(self.getFileType())
        outputSet.append(eFile)
        outputSet.write()
        self._defineOutputs(exportedFileSet=outputSet)

    # -------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        if self.getDatasets() not in (None, ''):
            summary.append(self.getDatasets())

        return summary
    # -------------------------- BASE methods to be overridden ---------------

    INPUT_EXPT_FILENAME = 'integrated_model.expt'
    INPUT_REFL_FILENAME = 'integrated_reflections.refl'

    def _initialParams(self, program):
        # Base method that can more easily be overridden when needed
        params = (f"{self.getInputModelFile()} {self.getInputReflFile()} "
                  f"{self.getOutput()} "
                  f"output.log={self.getLogFilePath(program)}")

        return params

    def _extraParams(self):
        params = ""
        if self.getFormat() is MTZ:
            params += self.mtzExtraParams()

        elif self.getFormat() is SADABS:
            params += self.sadabsExtraParams()

        elif self.getFormat() is NXS:
            params += self.nxsExtraParams()

        elif self.getFormat() is JSON:
            params += self.jsonExtraParams()

        elif self.getFormat() is SHELX:
            params += self.shelxExtraParams()

        elif self.getFormat() is PETS:
            params += self.petsExtraParams()

        return params

    # -------------------------- UTILS functions ------------------------------

    def getFormat(self):
        return self.exportFormat.get()

    def getFileType(self):
        if self.getFormat() is MTZ:
            filetype = "mtz"

        if self.getFormat() is SADABS:
            filetype = "sad"

        if self.getFormat() is NXS:
            filetype = "nxs"

        if self.getFormat() is MMCIF:
            filetype = "cif"

        if self.getFormat() is XDS_ASCII:
            filetype = "XDS_ASCII"

        if self.getFormat() is JSON:
            filetype = "json"

        if self.getFormat() is SHELX:
            filetype = "shelx"

        if self.getFormat() is PETS:
            filetype = "cif_pets"

        return filetype

    def getExport(self):
        name = None
        if self.getFormat() is MTZ:
            if self.mtzHklout.get() == "":
                name = f"integrated_{self.getObjId()}.mtz"
            else:
                name = self.mtzHklout.get()

        if self.getFormat() is SADABS:
            name = self.sadabsHklout.get()

        if self.getFormat() is NXS:
            name = self.nxsHklout.get()

        if self.getFormat() is MMCIF:
            if self.mmcifHklout.get() == "":
                name = f"integrated_{self.getObjId()}.cif"
            else:
                name = self.mmcifHklout.get()

        if self.getFormat() is XDS_ASCII:
            name = self.xdsAsciiHklout.get()

        if self.getFormat() is JSON:
            name = self.jsonFilename.get()

        if self.getFormat() is SHELX:
            name = self.shelxName.get()

        if self.getFormat() is PETS:
            name = self.petsPrefix.get()

        return self.outDir(name)

    def outDir(self, fn=None):
        if fn is None:
            return self._getExtraPath()
        else:
            return self._getExtraPath(fn)

    def getOutput(self):
        mtzStr = f"mtz.hklout={self.getExport()}"

        sadabsStr = f"sadabs.hklout={self.getExport()}"

        nxsStr = f"nxs.hklout={self.getExport()}"

        mmcifStr = f"mmcif.hklout={self.getExport()}"

        xdsAsciiStr = f"xds_ascii.hklout={self.getExport()}"

        jsonStr = f"json.filename={self.getExport()}"

        shelxStr = (f"shelx.hklout={self.getExport()}.hkl "
                    f"shelx.ins={self.getExport()}.ins")

        petsStr = f"pets.filename_prefix={self.getExport()}"

        idx = self.getFormat()

        formats = exportFormatChoices

        nameStr = [mtzStr, sadabsStr, nxsStr, mmcifStr,
                   xdsAsciiStr, jsonStr, shelxStr, petsStr]

        outputString = f"format={formats[idx]} {nameStr[idx]}"
        return outputString

    def mtzExtraParams(self):
        params = ""
        if self.mtzCombinePartials:
            params += " mtz.combine_partials=True"

            params += (f" mtz.partiality_threshold="
                       f"{self.mtzPartialityThreshold.get()}")

            params += f" mtz.min_isigi={self.mtzMinIsigi.get()}"

        if self.mtzForceStaticModel:
            params += " mtz.force_static_model=True"

        if self.mtzFilter_ice:
            params += " mtz.filter_ice_rings=True"

        if self.mtzDMin.get() is not None:
            params += f" mtz.d_min={self.mtzDMin.get()}"

        params += (f" mtz.crystal_name="
                   f"{self.getCrystalName(self.mtzCrystalName.get())}")
        params += f" mtz.project_name={self.getProjectName()}"
        return params

    def sadabsExtraParams(self):
        params = ""
        if self.sadabsRun.get() != 1:
            params += f" sadabs.run={self.sadabsRun.get()}"

        if self.sadabsPredict:
            params += " sadabs.predict=True"
        return params

    def nxsExtraParams(self):
        params = ""
        params += (f" nxs.instrument_name="
                   f"{self.nxsInstrumentName.get()}")

        params += (f" nxs.instrument_short_name="
                   f"{self.nxsInstrumentShortName.get()}")

        params += f" nxs.source_name={self.nxsSourceName.get()}"

        params += (f" nxs.source_short_name="
                   f"{self.nxsSourceShortName.get()}")
        return params

    def jsonExtraParams(self):
        params = ""
        if self.jsonCompact is False:
            params += " json.compact=False"

        if self.jsonNDigits.get() != 6:
            params += f" json.n_digits={self.jsonNDigits.get()}"
        return params

    def shelxExtraParams(self):
        params = ""
        if self.shelxScale.get():
            params += (f" scale=True scale_range="
                       f"{self.shelxScaleMin.get()},"
                       f"{self.shelxScaleMax.get()}")
        else:
            params += f" scale=False"

        return params

    def petsExtraParams(self):
        params = ""
        petsID = self.petsId.get()
        eeCutoff = self.petsVFExcitationErrorCutoff.get()
        nMerged = self.petsVFNMerged.get()
        step = self.petsVFStep.get()

        if petsID not in ("", "None", None):
            params += f" id={petsID}"
        params += f" partiality_cutoff={self.petsPartialityCutoff.get()}"
        if self.petsFlagFiltering.get():
            params += f" flag_filtering=True"
        else:
            params += f" flag_filtering=False"
        if eeCutoff is not None:
            params += f" excitation_error_cutoff={eeCutoff}"
        if nMerged is not None:
            params += f" n_merged={nMerged}"
        if step is not None:
            params += f" step={step}"

        return params
