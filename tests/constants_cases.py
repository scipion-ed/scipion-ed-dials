# **************************************************************************
# *
# * Author(s):     Viktor E. G. Bengtsson (viktor.e.g.bengtsson@gmail.com)
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
from __future__ import annotations

import tests.testing_utils as tutils

# Create toggles for skipping some tests
SKIP_PIPELINES = False
SKIP_LYSO = False
SKIP_GARNET = False
if SKIP_LYSO and SKIP_GARNET:
    SKIP_PIPELINES = True
SKIP_UTILS = False
KEEP_PROTOCOL_TEST_OUTPUT = True
KEEP_UTILS_TEST_OUTPUT = False

# -------------------------- Test case dictionaries --------------------------

lyso_experiment_14 = {
    "location": "lysozyme/experiment_14",
    "files_pattern": "SMV/data/00{TI}.img",
    "d_min": 2.5,
    "found_spots": 1322,
    "distance": 1480.56,
    "rotation_axis": "-0.6204,-0.7843,0.0",
    "lyso": True,
    "sample": "lyso",
    "scan_range": "1,49",
    "space_group": "P 4 2 2",
    "unit_cell": "77.76,77.76,40.5,90,90,90",
    "unit_cell_sigmas": "0.05,0.05,0.05,0.05,0.05,0.05",
    "cb_op": "b,c,a",
}

lyso_experiment_24 = {
    "location": "lysozyme/experiment_24",
    "files_pattern": "SMV/data/00{TI}.img",
    "d_min": 3.5,
    "found_spots": 847,
    "distance": None,
    "rotation_axis": "-0.6204,-0.7843,0.0",
    "lyso": True,
    "sample": "lysozyme",
    "scan_range": "1,46",
    "space_group": "P 4 2 2",
    "unit_cell": "78.24,78.24,39.75,90,90,90",
    "unit_cell_sigmas": "0.05,0.05,0.05,0.05,0.05,0.05",
    "cb_op": "b,c,a",
}

garnet_experiment = {
    "location": "garnet/experiment_1",
    "files_pattern": "Garnet_Cont_{TI}.img",
    "sample": "garnet",
    "replace_rotation_axis": True,
    "rotation_axis": "-0.755,-0.656,0.0",
    "useTemplate": True,
    "tsReplacement": "0###",
    "min_spot_size": 10,
    "max_spot_size": 1500000,
    "max_separation": 7.5,
    "d_min": 0.5,
    "d_max": 3.0,
    "sigma_background": 1.0,
    "found_spots": 6223,
    "nproc": 8,
    "bravais_setting_reference": "bravais_setting_22.expt",
    "cb_op": "b+c,a+c,a+b",
}

# --------------------------- Pre-defined strings ---------------------------

beamParams = (
    "refinement.parameterisation.beam.fix='all "
    "*in_spindle_plane out_spindle_plane *wavelength'"
)
beamParamsAll = (
    "refinement.parameterisation.beam.fix='*all "
    "in_spindle_plane out_spindle_plane wavelength'"
)
beamParamsSv = (
    "refinement.parameterisation.beam.fix='all "
    "in_spindle_plane out_spindle_plane *wavelength'"
)
crystalParams = (
    "refinement.parameterisation.crystal.fix='all cell orientation'"
)
detectorParamsAll = (
    "refinement.parameterisation.detector.fix='*all "
    "position orientation distance'"
)
detectorParamsDistance = (
    "refinement.parameterisation.detector.fix='all "
    "position orientation *distance'"
)
detectorParamsNoFix = (
    "refinement.parameterisation.detector.fix='all "
    "position orientation distance'"
)
gonioParams = (
    "refinement.parameterisation.goniometer.fix='*all "
    "in_beam_plane out_beam_plane'"
)

# -------------------------- UTILS functions ------------------------------


def lysoImportCommandLine(
    location: str,
    extraPath: str,
    logPath: str,
    experiment: dict,
) -> str:
    pathlist = tutils.makePathList(
        experiment["scan_range"],
        location=location,
        pattern=experiment["files_pattern"],
    )
    if experiment["distance"] is not None:
        distance = f" distance={experiment['distance']}"
    else:
        distance = ""
    command_line = (
        f"{' '.join(pathlist)} "
        f"output.log={logPath}/dials.import.log "
        f"output.experiments={extraPath}/imported.expt "
        f"goniometer.axes={experiment['rotation_axis']}"
        f"{distance.rstrip()}"
    )
    return command_line


def garnetImportCommandLine(
    dataset: str,
    extraPath: str,
    logPath: str,
    experiment: dict,
) -> str:
    command_line = (
        f"template={dataset}/Garnet_Cont_{experiment['tsReplacement']}.img "
        f"output.log={logPath}/dials.import.log "
        f"output.experiments={extraPath}/imported.expt "
        f"goniometer.axes={experiment['rotation_axis']}"
    )
    return command_line


def lysoStandardFindSpotCommandLine(
    inputExtraPath: str, outputExtraPath: str, logPath: str, experiment: dict
) -> str:
    command_line = (
        f"{inputExtraPath}/imported.expt "
        f"output.log={logPath}/dials.find_spots.log "
        f"output.reflections={outputExtraPath}/strong.refl "
        f"spotfinder.scan_range={experiment['scan_range']} "
        f"spotfinder.filter.d_min={experiment['d_min']} "
        f"spotfinder.filter.max_spot_size=1000 "
        f"spotfinder.filter.max_strong_pixel_fraction=0.25 "
        f"spotfinder.filter.max_separation=2.0 "
        f"spotfinder.filter.untrusted.rectangle=0,516,255,261 "
        f"spotfinder.filter.untrusted.rectangle=255,261,0,516 "
        f"spotfinder.threshold.algorithm=dispersion_extended "
        f"spotfinder.threshold.dispersion.sigma_background=6.0 "
        f"spotfinder.threshold.dispersion.sigma_strong=3.0 "
        f"spotfinder.threshold.dispersion.kernel_size=3,3"
    )
    return command_line


def garnetStandardFindSpotCommandLine(
    inputExtraPath: str, outputExtraPath: str, logPath: str, experiment: dict
) -> str:
    command_line = (
        f"{inputExtraPath}/imported.expt "
        f"output.log={logPath}/dials.find_spots.log "
        f"output.reflections={outputExtraPath}/strong.refl "
        f"spotfinder.scan_range=1,566 "
        f"spotfinder.filter.d_min={experiment['d_min']} "
        f"spotfinder.filter.d_max={experiment['d_max']} "
        f"spotfinder.filter.min_spot_size={experiment['min_spot_size']} "
        f"spotfinder.filter.max_spot_size={experiment['max_spot_size']} "
        f"spotfinder.filter.max_strong_pixel_fraction=0.25 "
        f"spotfinder.filter.max_separation={experiment['max_separation']} "
        f"spotfinder.threshold.algorithm=dispersion "
        f"spotfinder.threshold.dispersion.sigma_background="
        f"{experiment['sigma_background']} "
        f"spotfinder.threshold.dispersion.sigma_strong=3.0 "
        f"spotfinder.threshold.dispersion.kernel_size=3,3"
    )
    return command_line


def lysoIndexCommandLine(
    modelPath: str,
    spotsPath: str,
    logPath: str,
    indexTmp: str,
    spaceGroup: str | None = None,
) -> str:
    if spaceGroup:
        spacegroupstring = f"indexing.known_symmetry.space_group={spaceGroup} "
    else:
        spacegroupstring = ""
    indexCL = (
        f"{modelPath}/imported.expt "
        f"{spotsPath}/strong.refl "
        f"output.log={logPath}/dials.index.log "
        f"output.experiments={indexTmp}/indexed.expt "
        f"output.reflections={indexTmp}/indexed.refl "
        f"{spacegroupstring}"
        f"{beamParamsAll} {crystalParams} "
        f"{detectorParamsAll} {gonioParams}"
    )
    return indexCL


def garnetIndexCommandLine(
    modelPath: str,
    spotsPath: str,
    logPath: str,
    indexTmp: str,
) -> str:
    indexCL = (
        f"{modelPath}/imported.expt "
        f"{spotsPath}/strong.refl "
        f"output.log={logPath}/dials.index.log "
        f"output.experiments={indexTmp}/indexed.expt "
        f"output.reflections={indexTmp}/indexed.refl "
        f"{beamParams} {crystalParams} "
        f"{detectorParamsDistance} {gonioParams}"
    )
    return indexCL


def refineBravaisLatticeCommandLine(
    indexTmp: str, indexLogs: str, copied: bool = False
) -> str:
    if copied:
        copied_parameters = (
            f" {beamParams} {crystalParams}"
            f" {detectorParamsDistance} {gonioParams}"
        )
    else:
        copied_parameters = ""
    refBravCL = (
        f"{indexTmp}/indexed.expt "
        f"{indexTmp}/indexed.refl "
        f"output.log={indexLogs}/dials.refine_bravais_settings.log "
        f"output.directory={indexTmp}"
        f"{copied_parameters}"
    )
    return refBravCL


def reindexCommandLine(experiment: dict, indexTmp: str) -> str:
    reindexCL = (
        f"change_of_basis_op={experiment['cb_op']} {indexTmp}/indexed.refl "
        f"output.reflections={indexTmp}/reindexed.refl"
    )
    return reindexCL


def lysoRefineCommandLine(
    inputExtraPath: str,
    outputExtraPath: str,
    logPath: str,
    scan_varying: bool = False,
    static: bool = False,
) -> str:

    restraints = ""
    detector_parameters = detectorParamsAll
    file_base = "indexed"
    beam_params = beamParams

    if scan_varying:
        beam_params = (
            f"scan_varying=True {beamParamsSv} beam.force_static=False"
        )
        file_base = "refined"
    elif static:
        beam_params = f"scan_varying=False {beamParams}"
        restraints = f" {outputExtraPath}/restraints.phil"
        detector_parameters = detectorParamsNoFix

    command_line = (
        f"{inputExtraPath}/{file_base}.expt "
        f"{inputExtraPath}/{file_base}.refl "
        f"output.log={logPath}/dials.refine.log "
        f"output.experiments={outputExtraPath}/refined.expt "
        f"output.reflections={outputExtraPath}/refined.refl "
        f"{beam_params} {crystalParams} "
        f"{detector_parameters} {gonioParams}"
        f"{restraints}"
    )
    return command_line


def garnetRefineCommandLine(
    inputExtraPath: str,
    outputExtraPath: str,
    logPath: str,
) -> str:

    command_line = (
        f"{inputExtraPath}/indexed.expt "
        f"{inputExtraPath}/indexed.refl "
        f"output.log={logPath}/dials.refine.log "
        f"output.experiments={outputExtraPath}/refined.expt "
        f"output.reflections={outputExtraPath}/refined.refl "
        f"{beamParams} {crystalParams} "
        f"{detectorParamsDistance} {gonioParams}"
    )
    return command_line


def integrateCommandLine(
    inputExtraPath: str, outputExtraPath: str, logPath: str, experiment: dict
) -> str:
    command_line = (
        f"{inputExtraPath}/refined.expt "
        f"{inputExtraPath}/refined.refl "
        f"output.log={logPath}/dials.integrate.log "
        f"output.experiments={outputExtraPath}/"
        f"integrated_model.expt "
        f"output.reflections={outputExtraPath}/"
        f"integrated_reflections.refl "
        f"output.phil={outputExtraPath}/dials.integrate.phil "
        f"nproc=8 prediction.d_min={experiment['d_min']}"
    )
    return command_line


def symmetryCommandLine(
    inputExtraPath: str, outputExtraPath: str, logPath: str
) -> str:
    command_line = (
        f"{inputExtraPath}/integrated_model.expt "
        f"{inputExtraPath}/integrated_reflections.refl "
        f"output.log={logPath}/dials.symmetry.log "
        f"output.experiments={outputExtraPath}/symmetrized.expt "
        f"output.reflections={outputExtraPath}/symmetrized.refl "
        f"output.html={outputExtraPath}/dials.symmetry.html "
        f"output.json={outputExtraPath}/dials.symmetry.json"
    )
    return command_line


def scalingCommandLine(
    inputExtraPath: str, outputExtraPath: str, logPath: str
) -> str:
    command_line = (
        f"{inputExtraPath}/integrated_model.expt "
        f"{inputExtraPath}/integrated_reflections.refl "
        f"output.log={logPath}/dials.scale.log "
        f"output.experiments={outputExtraPath}/scaled.expt "
        f"output.reflections={outputExtraPath}/scaled.refl "
        f"output.html={outputExtraPath}/dials.scale.html "
        f"filtering.output.scale_and_filter_results="
        f"{outputExtraPath}/scale_and_filter_results.json "
        f"cut_data.partiality_cutoff=0.4 cut_data.min_isigi=-5.0 "
        f"outlier_rejection=standard outlier_zmax=6.0 filtering.method=None "
        f"filtering.deltacchalf.max_cycles=6 "
        f"filtering.deltacchalf.max_percent_removed=10.0 "
        f"filtering.deltacchalf.mode=dataset "
        f"filtering.deltacchalf.group_size=10 "
        f"filtering.deltacchalf.stdcutoff=4.0"
    )
    return command_line


def lysoMultiScalingCommandLine(
    inputExtraList: list[str],
    outputExtraPath: str,
    logPath: str,
    mergedFile: str,
    unmergedFile: str,
    crystalName: str,
    projectName: str,
) -> str:
    temporarylist = []
    for inputExtra in inputExtraList:
        temporarylist.append(
            f"{inputExtra}/scaled.expt {inputExtra}/scaled.refl "
        )

    command_line = (
        f"{''.join(temporarylist)}"
        f"output.log={logPath}/dials.scale.log "
        f"output.experiments={outputExtraPath}/scaled.expt "
        f"output.reflections={outputExtraPath}/scaled.refl "
        f"output.html={outputExtraPath}/dials.scale.html "
        f"output.merged_mtz={mergedFile} "
        f"output.unmerged_mtz={unmergedFile} "
        f"output.crystal_name={crystalName} "
        f"output.project_name={projectName} "
        f"filtering.output.scale_and_filter_results="
        f"{outputExtraPath}/scale_and_filter_results.json "
        f"cut_data.partiality_cutoff=0.4 cut_data.min_isigi=-5.0 "
        f"scaling_options.check_consistent_indexing=True "
        f"outlier_rejection=standard outlier_zmax=6.0 filtering.method=None "
        f"filtering.deltacchalf.max_cycles=6 "
        f"filtering.deltacchalf.max_percent_removed=10.0 "
        f"filtering.deltacchalf.mode=dataset "
        f"filtering.deltacchalf.group_size=10 "
        f"filtering.deltacchalf.stdcutoff=4.0"
    )
    return command_line


def lysoMergeCommandLine(
    inputExtraPath: str,
    outputExtraPath: str,
    logPath: str,
    crystalName: str,
    projectName: str,
    mergeBins: int,
) -> str:
    command_line = (
        f"{inputExtraPath}/scaled.expt "
        f"{inputExtraPath}/scaled.refl "
        f"output.log={logPath}/dials.merge.log "
        f"output.html={outputExtraPath}/dials.merge.html "
        f"output.mtz={outputExtraPath}/merged.mtz "
        f"output.crystal_names={crystalName} "
        f"output.project_name={projectName}"
        f"assess_space_group=True "
        f"anomalous=True "
        f"truncate=True "
        f"wavelength_tolerance={1e-4} "
        f"combine_partials=True "
        f"partiality_threshold=0.4 "
        f"n_residues=200 "
        f"merging.use_internal_variance=False "
        f"merging.n_bins={mergeBins} "
        f"merging.anomalous=False"
    )
    return command_line


def exportMtzCommandLine(
    inputExtraPath: str,
    outputExtraPath: str,
    logPath: str,
    objectID,
    projectName: str,
    crystalName: str = "XTAL",
) -> str:
    command_line = (
        f"{inputExtraPath}/scaled.expt "
        f"{inputExtraPath}/scaled.refl "
        f"format=mtz mtz.hklout={outputExtraPath}/"
        f"integrated_{objectID}.mtz "
        f"output.log={logPath}/dials.export.log "
        f"mtz.combine_partials=True mtz.partiality_threshold=0.99 "
        f"mtz.min_isigi=-5.0 mtz.crystal_name={crystalName} "
        f"mtz.project_name={projectName}"
    )
    return command_line


def lysoScaleExclusionCommandLine(
    inputExtraPath: str,
    outputExtraPath: str,
    logPath: str,
    exclusionList: list[str],
) -> str:
    temporarylist = []
    for exclusion in exclusionList:
        temporarylist.append(f"exclude_images={exclusion}")

    command_line = (
        f"{inputExtraPath}/integrated_model.expt "
        f"{inputExtraPath}/integrated_reflections.refl "
        f"output.log={logPath}/dials.scale.log "
        f"output.experiments={outputExtraPath}/scaled.expt "
        f"output.reflections={outputExtraPath}/scaled.refl "
        f"output.html={outputExtraPath}/dials.scale.html "
        f"filtering.output.scale_and_filter_results="
        f"{outputExtraPath}/scale_and_filter_results.json "
        f"cut_data.partiality_cutoff=0.4 cut_data.min_isigi=-5.0 "
        f"outlier_rejection=standard outlier_zmax=6.0 filtering.method=None "
        f"filtering.deltacchalf.max_cycles=6 "
        f"filtering.deltacchalf.max_percent_removed=10.0 "
        f"filtering.deltacchalf.mode=dataset "
        f"filtering.deltacchalf.group_size=10 "
        f"filtering.deltacchalf.stdcutoff=4.0 "
        f"{' '.join(temporarylist)}"
    )
    return command_line
