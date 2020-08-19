''' Define constants for use in protocols and/or tests '''

# spot-finding algorithm choice
DISPERSION = 0
DISPERSION_EXTENDED = 1

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
