""" Define constants for use in protocols and/or tests """

import dials.utils as dutils

# spot-finding algorithm choice
DISPERSION = 0
DISPERSION_EXTENDED = 1
RADIAL_PROFILE = 2

thresholdAlgorithimChoices = ["dispersion", "dispersion extended"]
thresholdDefault = DISPERSION_EXTENDED

if dutils.isMinDialsVersion("3.9.0"):
    thresholdAlgorithimChoices.append("radial profile")
    thresholdDefault = RADIAL_PROFILE

# radial profile options
NO_PREPROCESSING = 0
NARROW = 1
WIDE = 2
blurChoices = ["no preprocessing", "narrow (3×3)", "wide (5×5)"]

# scan varying methods
AUTO = 0
STATIC = 1
SCAN_VARYING = 2
UNSET = 3

# scaling outlier rejection scheme
STANDARD = 0
SIMPLE = 1

# scaling filtering method
NONE = 0
DELTA_CC_HALF = 1

# scaling filtering mode
DATASET = 0
IMAGE_GROUP = 1

# scaling dataset selection
USE_ALL = 0
DATASET_SELECTION = 1
EXCLUDE_DATASETS = 2

# report external dependencies
REMOTE = 0
LOCAL = 1
EMBED = 2

# export formats
MTZ = 0
SADABS = 1
NXS = 2
MMCIF = 3
XDS_ASCII = 4
JSON = 5
SHELX = 6
PETS = 7

exportFormatChoices = ["mtz", "sadabs", "nxs", "mmcif", "xds_ascii", "json"]

if dutils.isMinDialsVersion("3.9.0"):
    exportFormatChoices.append("shelx")
    exportFormatChoices.append("pets")
