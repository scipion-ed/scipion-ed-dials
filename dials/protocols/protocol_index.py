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

import os
import re
from glob import glob

import pyworkflow.protocol as pwprot

from pwed.objects import DiffractionImage, SetOfDiffractionImages, DiffractionSpot, SetOfSpots, IndexedSpot, SetOfIndexedSpots
from pwed.protocols import EdProtIndexSpots
from pwed.convert import find_subranges
from dials.convert import writeJson, readRefl, writeRefl

# Create 'aliases' for indices in EnumParam choices

FFT3D = 0
FFT1D = 1
REALSPACE_GRID_SEARCH = 2
LOWRES_SPOT_MATCH = 3

FLOOD_FILL = 0
CLEAN = 1

N_SPOTS = 0
D_MIN = 1

TRIPLETS = 0
QUADS = 1

SIMPLE = 0
LOCAL = 1

AUTOMATIC = 0
SINGLE = 1
MULTIPLE = 2
HIERARCHICAL = 3

FAIL = 0
FIX = 1
REMOVE = 2

IMAGE = 0
BLOCK = 1

SIMPLE_LBFGS = 0
LBFGS_CURVS = 1
GAUSS_NEWTON = 2
LEV_MAR = 3
SPARSE_LEV_MAR = 3

REFINE_SHELLS = 0
REPREDICT_ONLY = 1
MODE_NONE = 2

FRACTION_OF_BIN_SIZE = 0
ABSOLUTE = 1

STATISTICAL = 0
STILLS = 1
CONSTANT = 2
EXTERNAL_DELTAPSI = 3

NULL = 0
AUTO = 1
MCD = 2
TUKEY = 3
SAUTER_POON = 4

LINEAR = 0
LOG = 1

DBSCAN = 0
HCLUSTER = 1

DISTANCE = 0
INCONSISTENT = 1


class DialsProtIndexSpots(EdProtIndexSpots):
    """ Protocol for indexing spots using Dials
    """
    _label = 'index'

    # -------------------------- DEFINE param functions -----------------------

    # TODO: Go through and add some more defaults and "Auto"
    def _defineParams(self, form):
        # EdProtIndexSpots._defineParams(self, form)

        # Check which parts of indexing to perform. Reindexing is probably to small to warrant
        # its own protocol.
        form.addSection(label='Processing steps')

        form.addParam('doIndex', pwprot.BooleanParam,
                      default=True, label='Do you want to index from start?')
        form.addParam('doRefineBravaisSettings', pwprot.BooleanParam,
                      default=True, label='Do you want to refine the Bravais settings?')
        form.addParam('doReindex', pwprot.BooleanParam,
                      default=True, label='Do you want to reindex after refining Bravais settings?',
                      condition='doRefineBravaisSettings')

        # The start of the actually relevant part.
        form.addSection(label='Indexing basics', condition='doIndex')

        form.addParam('inputImages', pwprot.PointerParam,
                      pointerClass='SetOfDiffractionImages',
                      label="Input diffraction images",
                      help="")

        form.addParam('inputSpots', pwprot.PointerParam,
                      pointerClass='SetOfSpots',
                      label='Input strong spots',
                      help="")

        # Help messages are copied from the DIALS documentation at
        # https://dials.github.io/documentation/programs/dials_index.html
        form.addParam('indexNproc', pwprot.IntParam,
                      label="How many processors do you want to use?",
                      default=1,
                      help="The number of processes to use.")

        form.addParam('enterSpaceGroup', pwprot.BooleanParam,
                      default=False, label='Use a known space group?')

        form.addParam('knownSpaceGroup', pwprot.StringParam,
                      condition='enterSpaceGroup', default='', label='Space group:')

        form.addParam('enterUnitCell', pwprot.BooleanParam,
                      default=False, label='Use a known unit cell?')

        form.addParam('knownUnitCell', pwprot.StringParam,
                      condition='enterUnitCell', default='', label='Unit cell:')

        form.addParam('combine_scans', pwprot.BooleanParam, label="Do you want to combine scans?",
                      default=False, expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('relative_length_tolerance', pwprot.FloatParam, default=0.1,
                      label='Relative length tolerance',
                      help="Relative tolerance for unit cell lengths in unit cell comparision.",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('absolute_angle_tolerance', pwprot.FloatParam, default=5,
                      label='Absolute angle tolerance',
                      help="Angular tolerance (in degrees) in unit cell comparison.",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('max_delta', pwprot.FloatParam, default=5,
                      label='Max delta',
                      help="Maximum allowed Le Page delta used in searching for basis "
                      "vector combinations that are consistent with the given symmetry.",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('indexMm_search_scope', pwprot.FloatParam, default=4.0,
                      help="Global radius of origin offset search.",
                      label='mm search scope',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('indexWide_search_binning', pwprot.FloatParam, default=2,
                      help="Modify the coarseness of the wide grid search for the beam centre.",
                      label='Wide search binning',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('indexMin_cell_volume', pwprot.FloatParam, default=25,
                      help="Minimum unit cell volume (in Angstrom^3).",
                      label='Min cell volume',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('indexMin_cell', pwprot.FloatParam, default=3,
                      help="Minimum length of candidate unit cell basis vectors (in Angstrom).",
                      label='Min_cell',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('indexMax_cell', pwprot.FloatParam, default=None,
                      label='Max_cell', allowsNull=True,
                      help="Maximum length of candidate unit cell basis vectors (in Angstrom).",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('misindexCheckGridScope', pwprot.IntParam, default=0,
                      help="Search scope for testing misindexing on h, k, l.",
                      label='Misindexing check grid scope',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        # New section
        form.addSection('Indexing method parameters')

        form.addParam('indexingMethod', pwprot.EnumParam,
                      label='Which indexing method do you want to use?',
                      choice=['fft3d', 'fft1d',
                              'realSpaceGridSearch', 'lowResSpotMatch'],
                      default=FFT3D)

        form.addParam('indexingOptimiseInitialBasisVectors', pwprot.BooleanParam,
                      label='Optimise initial basis vectors?', default=False,
                      help='No help currently available. Check '
                      'dials.github.io/documentation/programs/dials_index.html for recently added help text.',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('fft3dBIso', pwprot.FloatParam,
                      default=None, allowsNull=True,
                      label='b_iso for fft3d', condition='indexingMethod==FFT3D',
                      help='No help currently available. Check '
                      'dials.github.io/documentation/programs/dials_index.html for recently added help text.',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('fft3dRmsdCutoff', pwprot.FloatParam,
                      default=15, condition='indexingMethod==FFT3D',
                      label='rmsd cutoff for fft3d',
                      help='No help currently available. Check '
                      'dials.github.io/documentation/programs/dials_index.html for recently added help text.',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('fft3dPeakSearch', pwprot.EnumParam,
                      choice=['flood fill' 'clean'], default=FLOOD_FILL,
                      label='How should peak search be done?',
                      condition='indexingMethod==FFT3D',
                      help='No help currently available. Check '
                      'dials.github.io/documentation/programs/dials_index.html for recently added help text.',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('fft3dPeakVolumeCutoff', pwprot.FloatParam,
                      default=0.15, condition='indexingMethod==FFT3D',
                      label='Peak volume cutoff',
                      help='No help currently available. Check '
                      'dials.github.io/documentation/programs/dials_index.html for recently added help text.',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('fft3dReciprocSpaceGridNPoints', pwprot.IntParam,
                      default=256, condition='indexingMethod==FFT3D',
                      label='n_points',
                      help='No help currently available. Check '
                      'dials.github.io/documentation/programs/dials_index.html for recently added help text.',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('fft3dReciprocSpaceGridDMin', pwprot.FloatParam,
                      default=None, condition='indexingMethod==FFT3D',
                      label='Reciprocal space grid search d_min', allowsNull=True,
                      help='The high resolution limit in Angstrom for spots to include in  the initial indexing.',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('fft1dCharacteristicGrid', pwprot.FloatParam,
                      default=None,
                      allowsNull=True, condition='indexingMethod==FFT1D', expertLevel=pwprot.LEVEL_ADVANCED,
                      label='Characteristic grid( fft1d)',
                      help='Sampling frequency in radians. See Steller 1997.'
                      'If None, determine a grid sampling automatically using the input reflections, using at most 0.029 radians.'
                      )

        form.addParam('realSpaceGridSearchCharacteristicGrid', pwprot.FloatParam,
                      default=0.02, condition='indexingMethod==REALSPACE_GRID_SEARCH',
                      label='Characteristic grid (realspace grid search)',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('realSpaceGridSearchMaxVectors', pwprot.IntParam,
                      default=30, help='The maximum number of unique vectors to find in the grid search',
                      condition='indexingMethod==REALSPACE_GRID_SEARCH',
                      label='Maximum number of unique vectors',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('lowResSpotMatchCandidateSpotsLimitResolutionBy', pwprot.EnumParam,
                      choice=['n_spots' 'd_min'], default=N_SPOTS,
                      label='Limit resolution by', condition='indexingMethod==LOWRES_SPOT_MATCH'
                      )

        form.addParam('lowResSpotMatchCandidateSpotsDMin', pwprot.FloatParam,
                      default=15.0, label='d_min', condition='lowResSpotMatchCandidateSpotsLimitResolutionBy==D_MIN'
                      )

        form.addParam('lowResSpotMatchCandidateSpotsNSpots', pwprot.IntParam,
                      default=10, label='n_spots', condition='lowResSpotMatchCandidateSpotsLimitResolutionBy==N_SPOTS'
                      )

        form.addParam('lowResSpotMatchCandidateSpotsDStarTolerance', pwprot.FloatParam,
                      default=4.0, label='d* tolerance',
                      help="Number of sigmas from the centroid position for which to calculate d* bands",
                      condition='indexingMethod==LOWRES_SPOT_MATCH'
                      )

        form.addParam('lowResSpotMatchUseP1IndiciesAsSeeds', pwprot.BooleanParam,
                      label='Use P1 indices as seeds?', default=False,
                      condition='indexingMethod==LOWRES_SPOT_MATCH'
                      )

        form.addParam('lowResSpotMatchSearchDepth', pwprot.EnumParam,
                      choice=['triplets' 'quads'], default=TRIPLETS,
                      label='Search depth for low-res spot match?',
                      condition='indexingMethod==LOWRES_SPOT_MATCH'
                      )

        form.addParam('lowResSpotMatchBootstrapCrystal', pwprot.BooleanParam,
                      label='Bootstrap crystal for low-res spot match?', default=False,
                      condition='indexingMethod==LOWRES_SPOT_MATCH'
                      )

        form.addParam('lowResSpotMatchMaxPairs', pwprot.IntParam,
                      default=200, condition='indexingMethod==LOWRES_SPOT_MATCH',
                      label='Max pairs for low-res spot match?',
                      )

        form.addParam('lowResSpotMatchMaxTriplets', pwprot.IntParam,
                      default=600, condition='indexingMethod==LOWRES_SPOT_MATCH',
                      label='Max Triplets for low-res spot match?',
                      )

        form.addParam('lowResSpotMatchMaxQuads', pwprot.IntParam,
                      default=600, condition='indexingMethod==LOWRES_SPOT_MATCH',
                      label='Max Quads for low-res spot match?',
                      )

        # New section
        form.addSection('Max_cell_estimation',
                        expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('doFilter_ice', pwprot.BooleanParam, default=True,
                      label='Filter ice?',
                      help="Filter out reflections at typical ice ring resolutions before max_cell estimation.")

        form.addParam('doFilter_overlaps', pwprot.BooleanParam, default=True,
                      label='Filter overlaps?',
                      help="Filter out reflections with overlapping bounding boxes before max_cell estimation.")

        form.addParam('overlaps_border', pwprot.IntParam, default=0,
                      label='Overlaps border',
                      help="Optionally add a border around the bounding boxes before finding overlaps.")

        form.addParam('multiplier', pwprot.FloatParam, default=1.3,
                      help="Multiply the estimated maximum basis vector length by this value.",
                      label='Multiplier',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('step_size', pwprot.FloatParam, default=45,
                      label='Step size',
                      help="Step size, in degrees, of the blocks used to peform the max_cell  estimation.",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('nearest_neighbor_percentile', pwprot.FloatParam, default=None,
                      allowsNull=True,
                      label='Nearest neighbor percentile',
                      help="Percentile of NN histogram to use for max cell determination.",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('histogram_binning', pwprot.EnumParam,
                      label='Histogram binning',
                      choice=['linear', 'log'], default=LOG,
                      help="Choose between linear or logarithmic bins for nearest neighbour histogram analysis.")

        form.addParam('nn_per_bin', pwprot.IntParam, default=5,
                      label='Nearest neighbours per bin',
                      help="Target number of nearest neighbours per histogram bin.")

        form.addParam('max_height_fraction', pwprot.FloatParam, default=0.25,
                      label='Max height fraction',
                      expertLevel=pwprot.LEVEL_ADVANCED)

        # Section for index assignment
        form.addSection('Index assignment',
                        expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('indexAssignmentMethod', pwprot.EnumParam,
                      label='Assignment method',
                      choice=['simple' 'local'], default=SIMPLE,
                      help="Choose between simple 'global' index assignment and xds-style "
                      "'local' index assignment.")

        form.addParam('simpleHkl_tolerance', pwprot.FloatParam,
                      label='hkl tolerance',
                      condition='indexAssignmentMethod==SIMPLE', default=0.3,
                      help="Maximum allowable deviation from integer-ness for assigning a miller "
                      "index to a reciprocal lattice vector.")

        form.addParam('localEpsilon', pwprot.FloatParam, condition='indexAssignmentMethod==LOCAL',
                      label='Epsilon', default=0.05,
                      help="This corresponds to the xds parameter INDEX_ERROR=")

        form.addParam('localDelta', pwprot.IntParam, condition='indexAssignmentMethod==LOCAL',
                      label='Delta', default=8, help="This corresponds to the xds parameter INDEX_MAGNITUDE=")

        form.addParam('localLMin', pwprot.FloatParam, condition='indexAssignmentMethod==LOCAL', label='l_min',
                      default=0.8, help="This corresponds to the xds parameter INDEX_QUALITY=")

        form.addParam('localNearestNeighbours', pwprot.IntParam,
                      label='Nearest neighbours', condition='indexAssignmentMethod==LOCAL', default=20)

        # Setup some preliminary refinement
        form.addSection('preliminaryRefinement')

        form.addParam('refinementMode', pwprot.EnumParam,
                      label='Refinement mode',
                      choice=['refine_shells' 'repredict_only' 'None'], default=REFINE_SHELLS,
                      help="refine_shells: refine in increasing resolution cutoffs after indexing.\n"
                      "repredict_only: do not refine after indexing, just update spot predictions.\n"
                      "None: turns off all forms of refinement(currently only  applies to stills_indexer)",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('n_macro_cycles', pwprot.IntParam, default=5,
                      label='Number of macro cycles',
                      help="Maximum number of macro cycles of refinement, reindexing all reflections using"
                      "updated geometry at the  beginning of each cycle. Does not apply to stills.indexer=stills.",
                      condition='refinementMode!=MODE_NONE')

        form.addParam('d_min_step', pwprot.FloatParam, default=None,
                      label='d_min step', allowsNull=True,
                      help="Reduction per step in d_min for reflections to include in refinement.")

        form.addParam('d_min_start', pwprot.FloatParam,
                      default=None,
                      allowsNull=True, label='d_min start')

        form.addParam('d_min_final', pwprot.FloatParam,
                      default=None,
                      allowsNull=True, label='d_min final',
                      help="Do not ever include reflections below this value in refinement.")

        # New section
        form.addSection('multiple_lattice_search',
                        expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('cluster_analysis_search', pwprot.BooleanParam, label='Perform cluster analysis search?',
                      default=False, help="Perform cluster analysis search for multiple lattices.")

        form.addParam('recycle_unindexed_reflections_cutoff', pwprot.FloatParam,
                      default=0.1, label='Recycle unindexed reflections cutoff',
                      help="Attempt another cycle of indexing on the unindexed reflections  if "
                      "more than the fraction of input reflections are unindexed.")

        form.addParam('minimum_angular_separation', pwprot.FloatParam,
                      label='Minimum angular separation', default=5,
                      help="The minimum angular separation (in degrees) between two lattices.")

        form.addParam('max_lattices', pwprot.IntParam,
                      label='Max lattices', default=1)

        form.addParam('clusterAnalysisMethod', pwprot.EnumParam,
                      label='Cluster analysis method',
                      choice=['dbscan' 'hcluster'], default=DBSCAN)

        form.addParam('cutoff', pwprot.FloatParam,
                      label='Cutoff', default=15,
                      condition='clusterAnalysisMethod==HCLUSTER')

        form.addParam('cutoff_criterion', pwprot.EnumParam,
                      label='Cutoff criterion',
                      choice=['distance' 'inconsistent'], default=DISTANCE,
                      condition='clusterAnalysisMethod==HCLUSTER')

        form.addParam('dbscanEps', pwprot.FloatParam,
                      label='Eps', default=0.05,
                      condition='clusterAnalysisMethod==DBSCAN')

        form.addParam('dbscanMin_samples', pwprot.IntParam,
                      label='Min samples', default=30,
                      condition='clusterAnalysisMethod==DBSCAN')

        form.addParam('minClusterSize', pwprot.IntParam,
                      label='Min cluster size', default=20)

        form.addParam('intersectionUnionRatioCutoff',
                      pwprot.FloatParam, label='Intersection union ratio cutoff', default=0.4)

        # Section with options for stills
        form.addSection('stills')

        form.addParam('stillsIndexer', pwprot.EnumParam, choice=['Auto' 'stills' 'sequences'], default=AUTO,
                      help="Use the stills or sequences indexer.\n"
                      "Auto: choose based on the input imagesets (stills or sequences).",
                      expertLevel=pwprot.LEVEL_ADVANCED)

        form.addParam('ewald_proximity_resolution_cutoff', pwprot.FloatParam, default=2.0,
                      label='Ewald proximity resolution cutoff',
                      help="the acceptable volume of reciprocal space for spot prediction")

        form.addParam('stillsRefineAllCandidates', pwprot.BooleanParam,
                      label='Refine all candidates?', default=True,
                      help='If False, no attempt is made to refine the model from initial basis vector'
                      'selection. The indexing solution with the best RMSD is chosen.',
                      )

        form.addParam('stillsCandidateOutlierRejection', pwprot.BooleanParam,
                      label='Apply Sauter/ Poon (2010) outlier rejection?', default=True,
                      help='If Yes, while refining candiate basis solutions, also apply Sauter/ Poon (2010) outlier rejection',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('stillsRefineCandidatesWithKnownSymmetry', pwprot.BooleanParam,
                      label='Refine candidates with known symmetry?', default=False,
                      help='If False, when choosing the best set of candidate basis solutions, refine the candidates'
                      'in the P1  setting. If True, after indexing in P1, convert the candidates to the known'
                      'symmetry and apply the    corresponding change of basis   to the indexed reflections.',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('stillsRmsdMinPx', pwprot.FloatParam,
                      label='rmsd minimum in px', default=2,
                      help='Minimum acceptable RMSD for choosing candidate basis solutions (in pixels)'
                      )

        form.addParam('stillsEwaldProximalVolumeMax', pwprot.FloatParam,
                      label='Ewald proximal volume max', default=0.0025,
                      help='Maximum acceptable ewald proximal volume when choosing candidate basis solutions'
                      )

        # TODO: Implement isoforms. Want a group where the input can be given multiple times?

        form.addParam('stillsSetDomainSizeAngValue', pwprot.FloatParam,
                      label='Set domain size ang value', default=None,
                      allowsNull=True,
                      help='If specified, will set the domain size ang value and override the value determined from'
                      'nave refinement'
                      )

        form.addParam('stillsSetMosaicHalfDegValue', pwprot.FloatParam,
                      label='Set mosaic half deg value', default=None,
                      allowsNull=True,
                      help='If specified, will set the mosaic half degree value and override the value determined'
                      'from nave refinement'
                      )

        form.addSection('Refinement parameter configuration')

        form.addParam('refineNproc', pwprot.IntParam,
                      default=1, label='nproc',
                      help='The number of processes to use. Not all choices of refinement engine support nproc > 1.'
                      'Where multiprocessing is possible, it is helpful only in certain circumstances,'
                      'so this is not recommended for typical use.',
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        group = form.addGroup('Parametrisation')

        group.addParam('autoReductionMinNrefPerParameter', pwprot.IntParam,
                       default=5,
                       help='the smallest number of reflections per parameter for a model parameterisation below which'
                       'the parameterisation will not be made in full, but the action described below will be triggered.',
                       label='Auto reduction: Min nref per parameter',
                       )

        group.addParam('autoReductionAction', pwprot.EnumParam,
                       choice=['fail' 'fix' 'remove'], default=FAIL,
                       help='action to take if there are too few reflections across the experiments related'
                       'to a particular model parameterisation. If fail, an exception will be raised and refinement'
                       'will not proceed. If fix, refinement will continue but with the parameters relating to that'
                       'model remaining fixed at their initial values. If remove parameters relating to that model'
                       'will be fixed, and in addition all reflections related to that parameterisation will be removed.'
                       'This will therefore remove these reflections from other parameterisations of the global model'
                       'too. For example, if a crystal model could not be parameterised it will be excised completely'
                       'and not contribute to the joint refinement of the detector and beam. In the fix mode, reflections'
                       'emanating from that crystal will still form residuals and will contribute to detector and beam'
                       'refinement.',
                       label='Auto reduction action',
                       )

        group.addParam('paramScanVarying', pwprot.BooleanParam,
                       label='scan varying refinement', default=False,
                       help='Allow models that are not forced to be static to vary during the scan,'
                       'Auto will run one macrocycle with static then scan varying refinement for the crystal',
                       )

        group.addParam('composeModelPer', pwprot.EnumParam,
                       choice=['image' 'block'], default=BLOCK,
                       help='For scan-varying parameterisations, compose a new model either every image or within blocks'
                       'of a width specified in the reflections parameters. When this block width is larger than the'
                       'image width the result is faster, with a trade-off in accuracy',
                       label='Compose model per', expertLevel=pwprot.LEVEL_ADVANCED,
                       condition='paramScanVarying==True'
                       )

        group.addParam('blockWidth', pwprot.FloatParam,
                       default=1.0, label='Block width', condition='composeModelPer==BLOCK',
                       help="Width of a reflection 'block' (in degrees) determining how fine-grained the model"
                       "used for scan-varying prediction during refinement is.",
                       expertLevel=pwprot.LEVEL_ADVANCED
                       )

        group.addParam('setScanVaryingErrors', pwprot.BooleanParam,
                       label='Set scan varying errors?', default=False,
                       help='If scan-varying refinement is done, and if the estimated covariance of the model'
                       'states have been calculated by the minimiser, choose whether to return this to the models or not.'
                       'The default is not to, in order to keep the file size of the serialized model small. At the moment,'
                       'this only has an effect for crystal unit cell (B matrix) errors.',
                       condition='paramScanVarying==True'
                       )

        group.addParam('debugCentroidAnalysis', pwprot.BooleanParam,
                       label='Debug centroid analysis?', default=False,
                       help="Set True to write out a file containing the reflections used for centroid analysis for"
                       "automatic setting of the scan-varying interval width. This can then be analysed with dev.dials."
                       "plot_centroid_analysis (requires dials_scratch repository).", expertLevel=pwprot.LEVEL_ADVANCED
                       )

        group.addParam('beamFixInSpindlePlane', pwprot.BooleanParam,
                       label='Fix beam in spindle plane?', default=True,
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted to a"
                       "direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow refinement of"
                       "the unit cell volume.",
                       )

        group.addParam('beamFixOutSpindlePlane', pwprot.BooleanParam,
                       label='Fix beam out of spindle plane?', default=False,
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted to"
                       "a direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow refinement"
                       "of the unit cell volume.",
                       )

        group.addParam('beamFixWavelength', pwprot.BooleanParam,
                       label='Fix beam wavelength?', default=True,
                       help="Whether to fix beam parameters. By default, in_spindle_plane is selected, and one of the two"
                       "parameters is fixed. If a goniometer is present this leads to the beam orientation being restricted"
                       "to a direction in the initial spindle-beam plane. Wavelength is also fixed by default, to allow"
                       "refinement of the unit cell volume.",
                       )

        group.addParam('beamForceStatic', pwprot.BooleanParam,
                       label='Force static?', default=True,
                       help="Force a static parametrisation for the beam when doing scan-varying refinement",
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       condition='paramScanVarying==True'
                       )

        group.addParam('SmootherUseAbsolute', pwprot.BooleanParam,
                       label='Use number of intervals in place of interval width in degrees?', default=False,
                       help="Number of intervals between checkpoints if scan_varying refinement is requested.",
                       expertLevel=pwprot.LEVEL_ADVANCED, condition='paramScanVarying==True'
                       )

        group.addParam('SmootherIntervalWidthDegrees', pwprot.FloatParam,
                       default=36.0, help="Width of scan between checkpoints in degrees.",
                       label='Smoother: Interval width (degrees)',
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       condition='SmootherUseAbsolute==False'
                       )

        group.addParam('SmootherAbsoluteNumIntervals', pwprot.IntParam,
                       default=None,
                       allowsNull=True, help="Number of intervals between checkpoints if scan_varying refinement is requested."
                       "If set, this overrides interval_width_degrees",
                       label='Absolute number of intervals',
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       condition='SmootherUseAbsolute==True'
                       )

        group.addParam('crystalFixCell', pwprot.BooleanParam,
                       label='Crystal: Fix cell?', default=True,
                       help="Fix crystal parameters",
                       )

        group.addParam('crystalFixOrientation', pwprot.BooleanParam,
                       label='Crystal: Fix orientation?', default=True,
                       help="Fix crystal parameters",
                       )

        group.addParam('crystalUnitCellForceStatic', pwprot.BooleanParam,
                       label='Force static?', default=True,
                       help="Force a static parameterisation for the crystal unit cell when doing scan-varying refinement",
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       condition='paramScanVarying==True'
                       )

        group.addParam('crystalOrientationForceStatic', pwprot.BooleanParam,
                       label='Force static?', default=True,
                       help="Force a static parameterisation for the crystal orientation when doing scan-varying refinement",
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       condition='paramScanVarying==True'
                       )

        group.addParam('detectorPanels', pwprot.EnumParam,
                       choice=['automatic' 'single' 'multiple' 'hierarchical'], default=AUTOMATIC,
                       help="Select appropriate detector parameterisation."
                       "Both the single and multiple panel detector options treat the whole detector as a rigid body."
                       "The hierarchical parameterisation treats groups of panels as separate rigid bodies.",
                       label='Detector panel parametrisation', expertLevel=pwprot.LEVEL_ADVANCED,
                       )

        group.addParam('detectorHierarchyLevel', pwprot.IntParam,
                       default=0, help="Level of the detector hierarchy (starting from the root at 0) at which to"
                       "determine panel groups to parametrise independently",
                       label='Detector hierarchy level', expertLevel=pwprot.LEVEL_ADVANCED
                       )

        group.addParam('detectorFixPosition', pwprot.BooleanParam,
                       label='Fix detector position?', default=True,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.",
                       )

        group.addParam('detectorFixOrientation', pwprot.BooleanParam,
                       label='Fix detector orientation?', default=True,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.",
                       )

        group.addParam('detectorFixdistance', pwprot.BooleanParam,
                       label='Fix detector distance?', default=True,
                       help="Fix detector parameters. The translational parameters (position) may be set"
                       "separately to the orientation.",
                       )

        group.addParam('crystalForceStatic', pwprot.BooleanParam,
                       label='Force static?', default=True,
                       help="Force a static parameterisation for the detector when doing scan-varying refinement",
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       condition='paramScanVarying==True'
                       )

        group.addParam('goniometerFixInBeamPlane', pwprot.BooleanParam,
                       label='Fix goniometer in beam plane?', default=True,
                       help="Whether to fix goniometer parameters. By default, fix all. Alternatively the setting matrix"
                       "can be constrained to allow rotation only within the spindle-beam plane or to allow rotation only"
                       "around an axis that lies in that plane. Set to None to refine the in two orthogonal directions.",
                       condition='paramScanVarying==True'
                       )

        group.addParam('goniometerFixOutBeamPlane', pwprot.BooleanParam,
                       label='Fix goniometer outside beam plane?', default=True,
                       help="Whether to fix goniometer parameters. By default, fix all. Alternatively the setting matrix"
                       "can be constrained to allow rotation only within the spindle-beam plane or to allow rotation only"
                       "around an axis that lies in that plane. Set to None to refine the in two orthogonal directions.",
                       condition='paramScanVarying==True'
                       )

        group.addParam('goniometerForceStatic', pwprot.BooleanParam,
                       label='Force static?', default=True,
                       help="Force a static parameterisation for the goniometer when doing scan-varying refinement",
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       condition='paramScanVarying==True'
                       )

        form.addParam('treatSingleImageAsStill', pwprot.BooleanParam,
                      label='Treat a single image as a still?', default=False,
                      help="Set this to True to treat a single image scan with a non zero oscillation width as a still",
                      expertLevel=pwprot.LEVEL_ADVANCED
                      )

        form.addParam('sphericalRelpModel', pwprot.BooleanParam,
                      label='Use the spherical relp model?', default=False,
                      help="For stills refinement, set true to use the spherical relp model for prediction and gradients.",
                      expertLevel=pwprot.LEVEL_ADVANCED,
                      )
        group = form.addGroup('Refinery')

        group.addParam('refineryEngine', pwprot.EnumParam,
                       choice=['SimpleLBFGS' 'LBFGScurvs' 'GaussNewton' 'LevMar' 'SparseLevMar'], default=LEV_MAR,
                       help="The minimisation engine to use",
                       label='Minimisation engine: ',
                       )

        group.addParam('refineryMaxIterations', pwprot.IntParam,
                       default=None,
                       allowsNull=True, help="Maximum number of iterations in refinement before termination."
                       "None implies the engine supplies its own default.",
                       label='Max iterations',
                       )
        group.addParam('refineryLog', pwprot.FileParam,
                       default=None,
                       allowsNull=True, label="Minimisation engine log",
                       help="Filename for an optional log that a minimisation engine may use to write"
                       "additional information",)

        group.addParam('refineryJournalTrackStep', pwprot.BooleanParam,
                       label='Record parameter shifts?', default=False,
                       help="Record parameter shifts history in the refinement journal, if the engine supports it.",
                       )

        group.addParam('refineryJournalTrackGradient', pwprot.BooleanParam,
                       label='Record parameter gradients?', default=False,
                       help="Record parameter gradients history in the refinement journal, if the engine supports it.",
                       )

        group.addParam('refineryJournalTrackParameterCorrelation', pwprot.BooleanParam,
                       label='Record parameter parameter correlation?', default=False,
                       help="Record correlation matrix between columns of the Jacobian for each step of refinement.",
                       )

        group.addParam('refineryJournalTrackConditionNumber', pwprot.BooleanParam,
                       label='Record condition number?', default=False,
                       help="Record condition number of the Jacobian for each step of  refinement.",
                       )

        group.addParam('refineryJournalTrackOutOfSampleRmsd', pwprot.BooleanParam,
                       label='Record out of sample RMSDs?', default=False,
                       help="Record RMSDs calculated using the refined experiments with reflections not used"
                       "in refinement at each step. Only valid if a subset of input reflections was taken for refinement",
                       )

        group = form.addGroup('target')

        group.addParam('targetRmsdCutoff', pwprot.EnumParam,
                       choice=['fraction_of_bin_size' 'absolute'], default=FRACTION_OF_BIN_SIZE,
                       help="Method to choose rmsd cutoffs. This is currently either as a fraction of the discrete units of the spot positional data, i.e. (pixel width, pixel height, image thickness in phi), or a tuple of absolute values to use     as the cutoffs",
                       label='RMSD cutoff',
                       expertLevel=pwprot.LEVEL_ADVANCED
                       )

        group.addParam('targetBinSizeFraction', pwprot.FloatParam,
                       default=0.0, help="Set this to a fractional value, say 0.2, to make a cut off in the natural discrete units of positional data, viz., (pixel width, pixel height, image thickness in phi). This would then determine when the      RMSD target is achieved. Only used if rmsd_cutoff = fraction_of_bin_size.",
                       label='Target bin size fraction', condition="targetRmsdCutoff==FRACTION_OF_BIN_SIZE",
                       )

        group.addParam('targetAbsoluteCutoffs', pwprot.FloatParam,
                       default=None, allowsNull=True, label='Target absolute cutoffs',
                       help="Absolute Values for the RMSD target achieved cutoffs in X, Y and Phi. The units are (mm, mm, rad).",
                       condition="targetRmsdCutoff==ABSOLUTE",
                       )

        group.addParam('targetGradientCalculationBlockSize', pwprot.IntParam,
                       default=None, allowsNull=True,
                       help="Maximum number of reflections to use for gradient calculation. If there are more reflections than this in the manager then the minimiser must do the full calculation in blocks.",
                       label='Target gradient calculation blocksize',
                       )

        group = form.addGroup('Reflections')

        group.addParam('reflectionsPerDegree', pwprot.FloatParam,
                       default=None, allowsNull=True, expertLevel=pwprot.LEVEL_ADVANCED,
                       help="The number of centroids per degree of the sequence to use in refinement. The default (Auto) uses all reflections unless the dataset is wider than a single turn. Then the number of reflections may be reduced until a minimum of 100 per degree of the sequence is reached to speed up calculations. Set this to None to force use all of suitable reflections.",
                       label='Reflections per degree',
                       )

        group.addParam('minimumSampleSize', pwprot.IntParam,
                       default=1000, allowsNull=True,
                       help="cutoff that determines whether subsetting of the input reflection list is done",
                       label='Minimum sample size', expertLevel=pwprot.LEVEL_ADVANCED
                       )

        group.addParam('maximumSampleSize', pwprot.IntParam,
                       default=None, allowsNull=True,
                       help="The maximum number of reflections to use in refinement. Overrides reflections_per_degree if that produces a larger sample size.",
                       label='Maximum sample size', expertLevel=pwprot.LEVEL_ADVANCED
                       )

        group.addParam('randomSeed', pwprot.IntParam,
                       default=42, allowsNull=True,
                       help="Random seed to use when sampling to create a working set of reflections. May be int or None.",
                       label='Random seed', expertLevel=pwprot.LEVEL_ADVANCED
                       )

        group.addParam('closeToSpindleCutoff', pwprot.FloatParam,
                       default=0.02,
                       help="The inclusion criterion currently uses the volume of the parallelepiped formed by the spindle axis, the incident beam and the scattered beam. If this is lower than some value then the reflection is excluded from refinement. In detector space, these are the reflections located close to the rotation axis.",
                       label='Close to spindle cutoff', expertLevel=pwprot.LEVEL_ADVANCED
                       )

        group.addParam('trimScanEdges', pwprot.FloatParam,
                       default=0.0, expertLevel=pwprot.LEVEL_ADVANCED,
                       label='Trim scan edges',
                       help="",
                       )

        group.addParam('weightingStratOverride', pwprot.EnumParam,
                       choice=[
                           'statistical' 'stills' 'constant' 'external_deltapsi'],
                       label='Weighting strategy override', expertLevel=pwprot.LEVEL_ADVANCED,
                       help="selection of a strategy to override default weighting behaviour",
                       )

        group.addParam('delpsiConstant', pwprot.FloatParam,
                       default=1000000,
                       label='delpsi constant', expertLevel=pwprot.LEVEL_ADVANCED,
                       help="used by the stills strategy to choose absolute weight value for the angular distance from Ewald sphere term of the target function, whilst the X and Y parts use statistical weights",
                       )

        group.addParam('weightingStrategyConstants', pwprot.FloatParam,
                       default=1.0,
                       label='Constant weights for target function',
                       help="constant weights for three parts of the target function, whether the case is for stills or scans. The default gives unit weighting.", expertLevel=pwprot.LEVEL_ADVANCED,
                       )

        group.addParam('outlierAlgorithm', pwprot.EnumParam,
                       choice=['null' 'auto' 'mcd' 'tukey' 'sauter_poon'], default=AUTO,
                       label='Outlier rejection algorithm',
                       help="Outlier rejection algorithm. If auto is selected, the algorithm is chosen automatically",
                       )

        group.addParam('outlierMinimumNumberOfReflections', pwprot.IntParam,
                       default=20, allowsNull=True,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       label='Minimum number of reflections',
                       help="The minimum number of input observations per outlier rejection job below which all reflections in the job will be rejected as potential outliers.",
                       )

        group.addParam('outlierSeparateExperiments', pwprot.BooleanParam,
                       label='Separate experiments?', default=True,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       help="If true, outlier rejection will be performed on each experiment separately. Otherwise, the data from all experiments will be combined for outlier rejection.",
                       )

        group.addParam('outlierSeparatePanels', pwprot.BooleanParam,
                       label='Separate panels?', default=None,
                       expertLevel=pwprot.LEVEL_ADVANCED, allowsNull=True,
                       help="If true, outlier rejection will be performed separately for each panel of a multi-panel detector model. Otherwise data from across all panels will be combined for outlier rejection.",
                       )

        group.addParam('outlierSeparateBlocks', pwprot.BooleanParam,
                       label='Separate blocks?', default=True, expertLevel=pwprot.LEVEL_ADVANCED,
                       help="f true, for scans outlier rejection will be performed separately in equal-width blocks of phi, controlled by the parameter outlier.block_width.",
                       )

        group.addParam('outlierBlockWidth', pwprot.FloatParam,
                       default=None, label='Block width', expertLevel=pwprot.LEVEL_ADVANCED, allowsNull=True,
                       help="If separate_blocks, a scan will be divided into equal-sized blocks with width (in degrees) close to this value for outlier rejection. If Auto, a width of at least 18 degrees will be determined, such that each block contains enough reflections to perform outlier rejection.",
                       )

        group.addParam('outlierTukeyIqrMultiplier', pwprot.FloatParam,
                       default=1.5, condition="outlierAlgorithm==TUKEY",
                       label='iqr multiplier', expertLevel=pwprot.LEVEL_ADVANCED,
                       help="The IQR multiplier used to detect outliers. A value of 1.5 gives Tukey's rule for outlier detection",
                       )

        group.addParam('mcdAlpha', pwprot.FloatParam,
                       default=0.5, condition="outlierALgorithm==MCD",
                       label='mcd alpha', expertLevel=pwprot.LEVEL_ADVANCED,
                       help="Options for the mcd outlier rejector, which uses an algorithm based on FAST-MCD by Rousseeuw and van Driessen. See doi.org/10.1080/00401706.1999.10485670.\n"
                       "Decimal fraction controlling the size of subsets over which the covariance matrix determinant is minimised.",
                       )

        group.addParam('mcdMaxNGroups', pwprot.IntParam,
                       default=5, condition="outlierALgorithm==MCD",
                       label='max number of groups', allowsNull=True,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       help="The maximum number of groups to split the dataset into if the dataset is 'large' (more observations than twice the min_group_size).",
                       )

        group.addParam('mcdMinGroupSize', pwprot.IntParam,
                       default=300, condition="outlierALgorithm==MCD",
                       label='Min group size', allowsNull=True,
                       expertLevel=pwprot.LEVEL_ADVANCED,
                       help="The smallest sub-dataset size when splitting the dataset into a number of groups, maximally \"max number of groups\".",
                       )

        group.addParam('mcdNTrials', pwprot.IntParam,
                       default=500, allowsNull=True,
                       condition="outlierAlgorithm==MCD", expertLevel=pwprot.LEVEL_ADVANCED,
                       label='Number of trials',
                       help="The number of samples used for initial estimates to seed the search within each sub-dataset.",
                       )

        group.addParam('mcdK1', pwprot.IntParam,
                       label='k1', default=2, allowsNull=True,
                       condition="outlierAlgorithm==MCD", expertLevel=pwprot.LEVEL_ADVANCED,
                       help="The number of concentration steps to take after initial estimates.",
                       )

        group.addParam('mcdK2', pwprot.IntParam,
                       label='k2', default=2, allowsNull=True,
                       condition="outlierAlgorithm==MCD", expertLevel=pwprot.LEVEL_ADVANCED,
                       help="If the dataset is 'large', the number of concentration steps to take after applying the best subset estimates to the merged group.",
                       )

        group.addParam('mcdK3', pwprot.IntParam,
                       label='k3', default=100, allowsNull=True,
                       condition="outlierAlgorithm==MCD", expertLevel=pwprot.LEVEL_ADVANCED,
                       help="If the dataset is 'small', the number of concentration steps to take after selecting the best of the initial estimates, applied to the whole dataset.",
                       )

        group.addParam('mcdThresholdProbability', pwprot.FloatParam,
                       label='Threshold probability',
                       default=0.975,
                       condition="outlierAlgorithm==MCD", expertLevel=pwprot.LEVEL_ADVANCED,
                       help="Quantile probability from the Chi-squared distribution with number of degrees of freedom equal to the number of dimensions of the data data (e.g. 3 for X, Y and Phi residuals). Observations whose robust Mahalanobis distances are larger than the obtained quantile will be flagged as outliers.",
                       )

        group.addParam('sauterPoonPxSz', pwprot.FloatParam,
                       label='Pixel size (mm)', condition="outlierAlgorithm==SAUTER_POON",
                       default=None, expertLevel=pwprot.LEVEL_ADVANCED, allowsNull=True,
                       help="X, Y pixel size in mm. If Auto, this will be taken from the first panel of the first experiment.",
                       )

        group.addParam('sauterPoonVerbose', pwprot.BooleanParam,
                       label='Do your want verbose output?', default=False,
                       condition="outlierAlgorithm==SAUTER_POON", expertLevel=pwprot.LEVEL_ADVANCED,
                       )

        group.addParam('sauterPoonPdf', pwprot.StringParam,
                       label='Name for graph',
                       default='',
                       help="Output file name for making graphs of |dr| vs spot number and dy vs dx.",
                       )

        # Allow some options if the Bravais settings are to be refined
        form.addSection(label='Refine Bravais settings',
                        condition='doRefineBravaisSettings')

        form.addParam('refineBravNproc', pwprot.IntParam,
                      default=4, label="How many processors do you want to use?",
                      help="The number of processes to use.")

        # Allow some options that are only relevant for reindexing
        form.addSection(label='Reindex',
                        condition='doReindex')

        form.addParam('doReindexModel', pwprot.BooleanParam,
                      default=False, label="Do you want to reindex the experimental model?")

        form.addParam('doReindexReflections', pwprot.BooleanParam,
                      default=False, label="Do you want to reindex the experimental reflections?")

   # -------------------------- INSERT functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(
            'convertInputStep', self.inputImages.getObjId(), self.inputSpots.getObjId())
        if self.doIndex:
            self._insertFunctionStep('indexStep')
        if self.doRefineBravaisSettings:
            self._insertFunctionStep('refineBravaisStep')
        if self.doReindex:
            self._insertFunctionStep('reindexStep')
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, inputId):
        inputImages = self.inputImages.get()
        inputSpots = self.inputSpots.get()
        self.info("Number of images: %s" % inputImages.getSize())
        self.info("Number of spots: %s" % inputSpots.getSpots())
        writeJson(inputImages, fn=self.getModelFile())
        writeRefl(inputSpots, fn=self.getReflFile())

    def indexStep(self):
        program = 'dials.index'
        arguments = self._prepIndexCommandline(program)
        self.runJob(program, arguments)

    def refineBravaisStep(self):
        program = 'dials.refine_bravais_settings'
        arguments = self._prepBravaisCommandline(program)
        self.runJob(program, arguments)

    def indexStep(self):
        program = 'dials.reindex'
        arguments = self._prepReindexCommandline(program)
        self.runJob(program, arguments)

    def createOutputStep(self):
        reflectionData = readRefl(self.getReflFile())
        outputSet = self._createSetOfIndexedSpots()
        iSpot = IndexedSpot()
        numberOfSpots = reflectionData[2]
        reflDict = reflectionData[4]

        outputSet.setSpots(numberOfSpots)

        for i in range(0, numberOfSpots):
            iSpot.setObjId(i+1)
            iSpot.setSpotId(reflDict['id'][i])
            iSpot.setBbox(reflDict['bbox'][i])
            iSpot.setFlag(reflDict['flags'][i])
            iSpot.setIntensitySumValue(reflDict['intensity.sum.value'][i])
            iSpot.setIntensitySumVariance(
                reflDict['intensity.sum.variance'][i])
            iSpot.setNSignal(reflDict['n_signal'][i])
            iSpot.setPanel(reflDict['panel'][i])
            try:
                iSpot.setShoebox(reflDict['shoebox'][i])
            except IndexError:
                pass
            iSpot.setXyzobsPxValue(reflDict['xyzobs.px.value'][i])
            iSpot.setXyzobsPxVariance(reflDict['xyzobs.px.variance'][i])
            outputSet.append(iSpot)

        outputSet.write()

        self._defineOutputs(outputIndexedSpots=outputSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        return errors

    # -------------------------- UTILS functions ------------------------------
    def getModelFile(self):
        return self._getExtraPath('imported.expt')

    def getReflFile(self):
        return self._getExtraPath('strong.refl')

    def getChangeOfBasisOp(self):
        # TODO: extract change of basis op from result in refineBravaisStep
        change_of_basis_op = 'a,b,c'
        return change_of_basis_op

    def _prepIndexCommandline(self, program):
        "Create the command line input to run dials programs"
        # Input basic parameters
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        params = "{} output.log={} output.reflections={} {}".format(
            self.getModelFile(), logPath, self.getReflFile(), self._createScanRanges())

        # Update the command line with additional parameters

        if self.indexNproc.get() not in (None, 1):
            params += " indexing.nproc={}".format(self.indexNproc.get())

        if self.enterSpaceGroup.get():
            params += " indexing.known_symmetry.space_group={}".format(
                self.knownSpaceGroup.get())

        if self.enterUnitCell.get():
            params += " indexing.known_symmetry.unit_cell={}".format(
                self.knownUnitCell.get())

        if self.combine_scans.get():
            params += " indexing.combine_scans={}".format(
                self.combine_scans.get())

        if self.relative_length_tolerance.get() not in (None, 0.1):
            params += " relative_length_tolerance={}".format(
                self.relative_length_tolerance.get())

        if self.absolute_angle_tolerance.get() not in (None, 5):
            params += " absolute_angle_tolerance={}".format(
                self.absolute_angle_tolerance.get())

        if self.max_delta.get() not in (None, 5):
            params += " max_delta={}".format(self.max_delta.get())

        if self.indexMm_search_scope.get() not in (None, 4.0):
            params += " index.mm_search_scope={}".format(
                self.indexMm_search_scope.get())

        if self.indexWide_search_binning.get() not in (None, 2):
            params += " index.wide_search_binning={}".format(
                self.indexWide_search_binning.get())

        if self.indexMin_cell_volume.get() not in (None, 25):
            params += " index.min_cell_volume={}".format(
                self.indexMin_cell_volume.get())

        if self.indexMax_cell.get() is not None:
            params += " index.max_cell={}".format(self.indexMax_cell.get())

        if self.misindexCheckGridScope.get() not in (None, 0):
            params += " check_misindexing.grid_search_scope={}".format(
                self.misindexCheckGridScope.get())

        if self.indexingMethod.get() is FFT1D:
            params += " indexing.method=fft1d"
        elif self.indexingMethod.get() is REALSPACE_GRID_SEARCH:
            params += " indexing.method=real_space_grid_search"
        elif self.indexingMethod.get() is LOWRES_SPOT_MATCH:
            params += " indexing.method=low_res_spot_match"

        if self.indexingOptimiseInitialBasisVectors.get():
            params += " indexing.optimise_initial_basis_vectors={}".format(
                self.indexingOptimiseInitialBasisVectors.get())

        if self.fft3dBIso.get() is not None:
            params += " indexing.fft3d.b_iso={}".format(self.fft3dBIso.get())

        if self.fft3dRmsdCutoff.get() not in (None, 15):
            params += " indexing.fft3d.rmsd_cutoff={}".format(
                self.fft3dRmsdCutoff.get())

        if self.fft3dPeakSearch.get() not in (None, FLOOD_FILL):
            params += " indexing.fft3d.peak_search=clean".format(
                self.fft3dPeakSearch.get())

        if self.fft3dPeakVolumeCutoff.get() not in (None, 1.15):
            params += " indexing.fft3d.peak_volume_cutoff={}".format(
                self.fft3dPeakVolumeCutoff.get())

        if self.fft3dReciprocSpaceGridNPoints.get() not in (None, 256):
            params += " indexing.fft3d.reciprocal_space_grid.n_points={}".format(
                self.fft3dReciprocSpaceGridNPoints.get())

        if self.fft3dReciprocSpaceGridDMin.get() not in (None):
            params += " indexing.fft3d.reciprocal_space_grid.d_min={}".format(
                self.fft3dReciprocSpaceGridDMin.get())

        if self.fft1dCharacteristicGrid.get() not in (None):
            params += " indexing.fft1d.characteristic_grid={}".format(
                self.fft1dCharacteristicGrid.get())

        if self.realSpaceGridSearchCharacteristicGrid.get() not in (None, 0.02):
            params += " indexing.fft1d.real_space_grid_search.characteristic_grid={}".format(
                self.realSpaceGridSearchCharacteristicGrid.get())

        if self.realSpaceGridSearchMaxVectors.get() not in (None, 30):
            params += " indexing.fft1d.real_space_grid_search.max_vectors={}".format(
                self.realSpaceGridSearchMaxVectors.get())

        if self.lowResSpotMatchCandidateSpotsLimitResolutionBy.get() not in (None, N_SPOTS):
            params += " indexing.fft1d.low_res_spot_match.candidate_spots.limit_resolution_by=d_min"

        if self.lowResSpotMatchCandidateSpotsDMin.get() not in (None, 15.0):
            params += " indexing.fft1d.low_res_spot_match.candidate_spots.d_min={}".format(
                self.lowResSpotMatchCandidateSpotsDMin.get())

        if self.lowResSpotMatchCandidateSpotsNSpots.get() not in (None, 10):
            params += " indexing.fft1d.low_res_spot_match.candidate_spots.n_spots={}".format(
                self.lowResSpotMatchCandidateSpotsNSpots.get())

        if self.lowResSpotMatchCandidateSpotsDStarTolerance.get() not in (None, 4.0):
            params += " indexing.fft1d.low_res_spot_match.candidate_spots.d_star_tolerance={}".format(
                self.lowResSpotMatchCandidateSpotsDStarTolerance.get())

        if self.lowResSpotMatchUseP1IndiciesAsSeeds.get():
            params += " indexing.fft1d.low_res_spot_match.use_P1_indices_as_seeds={}".format(
                self.lowResSpotMatchUseP1IndiciesAsSeeds.get())

        if self.lowResSpotMatchSearchDepth.get() not in (None, TRIPLETS):
            params += " indexing.fft1d.low_res_spot_match.search_depth={}".format(
                self.lowResSpotMatchSearchDepth.get())
        if self.lowResSpotMatchBootstrapCrystal.get():
            params += " indexing.fft1d.low_res_spot_match.bootstrap_crystal={}".format(
                self.lowResSpotMatchBootstrapCrystal.get())

        if self.lowResSpotMatchMaxPairs.get() not in (None, 200):
            params += " indexing.fft1d.low_res_spot_match.max_pairs={}".format(
                self.lowResSpotMatchMaxPairs.get())

        if self.lowResSpotMatchMaxTriplets.get() not in (None, 600):
            params += " indexing.fft1d.low_res_spot_match.max_triplets={}".format(
                self.lowResSpotMatchMaxTriplets.get())

        if self.lowResSpotMatchMaxQuads.get() not in (None, 600):
            params += " indexing.fft1d.low_res_spot_match.max_quads={}".format(
                self.lowResSpotMatchMaxQuads.get())

        if self.doFilter_ice.get():
            params += " indexing.max_cell_estimation.filter_ice={}".format(
                self.doFilter_ice.get())

        if self.doFilter_overlaps.get():
            params += " indexing.max_cell_estimation.filter_overlaps={}".format(
                self.doFilter_overlaps.get())

        if self.overlaps_border.get() not in (None, 0):
            params += " indexing.max_cell_estimation.overlaps_border={}".format(
                self.overlaps_border.get())

        if self.multiplier.get() not in (None, 1.3):
            params += " indexing.max_cell_estimation.multiplier={}".format(
                self.multiplier.get())

        if self.step_size.get() not in (None, 45):
            params += " indexing.max_cell_estimation.step_size={}".format(
                self.step_size.get())

        if self.nearest_neighbor_percentile.get():
            params += " indexing.max_cell_estimation.nearest_neighbor_percentile={}".format(
                self.nearest_neighbor_percentile.get())

        if self.histogram_binning.get() not in (None, LOG):
            params += " indexing.max_cell_estimation.histogram_binning={}".format(
                self.histogram_binning.get())

        if self.nn_per_bin.get() not in (None, 5):
            params += " indexing.max_cell_estimation.nn_per_bin={}".format(
                self.nn_per_bin.get())

        if self.max_height_fraction.get() not in (None, 0.25):
            params += " indexing.max_cell_estimation.max_height_fraction={}".format(
                self.max_height_fraction.get())

        if self.indexAssignmentMethod.get() not in (None, SIMPLE):
            params += " indexing.index_assignment.method=local"

        if self.simpleHkl_tolerance.get() not in (None, 0.3):
            params += " indexing.index_assignment.simple.hkl_tolerance={}".format(
                self.simpleHkl_tolerance.get())

        if self.localEpsilon.get() not in (None, 0.05):
            params += " indexing.index_assignment.local.epsilon={}".format(
                self.localEpsilon.get())

        if self.localDelta.get() not in (None, 8):
            params += " indexing.index_assignment.local.delta={}".format(
                self.localDelta.get())

        if self.localLMin.get() not in (None, 0.8):
            params += " indexing.index_assignment.local.l_min={}".format(
                self.localLMin.get())

        if self.localNearestNeighbours.get() not in (None, 20):
            params += " indexing.index_assignment.local.nearest_neighbours={}".format(
                self.localNearestNeighbours.get())

        if self.refinementMode.get() is REPREDICT_ONLY:
            params += " indexing.refinement_protocol.mode=repredict_only"
        elif self.refinementMode.get() is MODE_NONE:
            params += " indexing.refinement_protocol.mode=None"

        if self.n_macro_cycles.get() not in (None, 5):
            params += " indexing.refinement_protocol.n_macro_cycles={}".format(
                self.n_macro_cycles.get())

        if self.d_min_step.get() is not None:
            params += " indexing.refinement_protocol.d_min_step={}".format(
                self.d_min_step.get())

        if self.d_min_start.get() is not None:
            params += " indexing.refinement_protocol.d_min_start={}".format(
                self.d_min_start.get())

        if self.d_min_final.get() is not None:
            params += " indexing.refinement_protocol.d_min_final={}".format(
                self.d_min_final.get())

        if self.cluster_analysis_search.get():
            params += " indexing.cluster_analysis_search={}".format(
                self.cluster_analysis_search.get())

        if self.recycle_unindexed_reflections_cutoff.get() not in (None, 0.1):
            params += " indexing.recycle_unindexed_reflections_cutoff={}".format(
                self.recycle_unindexed_reflections_cutoff.get())

        if self.minimum_angular_separation.get() not in (None, 5):
            params += " indexing.minimum_angular_separation={}".format(
                self.minimum_angular_separation.get())

        if self.max_lattices.get() not in (None, 1):
            params += " indexing.max_lattices={}".format(
                self.max_lattices.get())

        if self.clusterAnalysisMethod.get() not in (None, HCLUSTER):
            params += " indexing.cluster_analysis.method=hcluster".format(
                self.clusterAnalysisMethod.get())

        if self.cutoff.get() not in (None, 15):
            params += " indexing.cluster_analysis.hccluster.cutoff={}".format(
                self.cutoff.get())

        if self.cutoff_criterion.get() is INCONSISTENT:
            params += " indexing.cluster_analysis.hccluster.cutoff_criterion=inconsistent"

        if self.dbscanEps.get() not in (None, 0.05):
            params += " indexing.cluster_analysis.dbscan.eps={}".format(
                self.dbscanEps.get())

        if self.dbscanMin_samples.get() not in (None, 30):
            params += " indexing.cluster_analysis.dbscan.min_samples={}".format(
                self.dbscanMin_samples.get())

        if self.minClusterSize.get() not in (None, 20):
            params += " indexing.cluster_analysis.min_cluster_size={}".format(
                self.minClusterSize.get())

        if self.intersectionUnionRatioCutoff.get() not in (None, 0.4):
            params += " indexing.cluster_analysis.intersection_union_ratio_cutoff={}".format(
                self.intersectionUnionRatioCutoff.get())

        if self.stillsIndexer.get() is STILLS:
            params += " indexing.stills.indexer=stills"
        elif self.stillsIndexer.get() is SEQUENCES:
            params += " indexing.stills.indexer=sequences"

        if self.ewald_proximity_resolution_cutoff.get() not in (None, 2.0):
            params += " indexing.stills.ewald_proximity_resolution_cutoff={}".format(
                self.ewald_proximity_resolution_cutoff.get())

        if self.stillsRefineAllCandidates.get() not in (None, True):
            params += " indexing.stills.refine_all_candidates={}".format(
                self.stillsRefineAllCandidates.get())

        if self.stillsCandidateOutlierRejection.get() not in (None, True):
            params += " indexing.stills.candidate_outlier_rejection={}".format(
                self.stillsCandidateOutlierRejection.get())

        if self.stillsRefineCandidatesWithKnownSymmetry.get() not in (None, False):
            params += " indexing.stills.refine_candidates_with_known_symmetry={}".format(
                self.stillsRefineCandidatesWithKnownSymmetry.get())

        if self.stillsRmsdMinPx.get() not in (None, 2):
            params += " indexing.stills.rmsd_min_px={}".format(
                self.stillsRmsdMinPx.get())

        if self.stillsEwaldProximalVolumeMax.get() not in (None, 0.0025):
            params += " indexing.stills.ewald_proximal_volume_max={}".format(
                self.stillsEwaldProximalVolumeMax.get())

        if self.stillsSetDomainSizeAngValue.get() is not None:
            params += " indexing.stills.set_domain_size_ang_value={}".format(
                self.stillsSetDomainSizeAngValue.get())

        if self.stillsSetMosaicHalfDegValue.get() is not None:
            params += " indexing.stills.set_mosaic_half_deg_value={}".format(
                self.stillsSetMosaicHalfDegValue.get())

        if self.refineNproc.get() not in (None, 1):
            params += " refinement.nproc={}".format(self.refineNproc.get())

        if self.autoReductionMinNrefPerParameter.get() not in (None, 5):
            params += " refinement.parameterisation.auto_reduction.min_nref_per_parameter={}".format(
                self.autoReductionMinNrefPerParameter.get())

        if self.autoReductionAction.get() is FIX:
            params += " refinement.parameterisation.auto_reduction.action={}".format(
                'fix')
        elif self.autoReductionAction.get() is REMOVE:
            params += " refinement.parameterisation.auto_reduction.action={}".format(
                'remove')

        if self.paramScanVarying.get() not in (None, False):
            params += " refinement.parameterisation.scan_varying={}".format(
                self.paramScanVarying.get())

        if self.composeModelPer.get() is IMAGE:
            params += " refinement.parameterisation.compose_model_per={}".format(
                'image')

        if self.blockWidth.get() not in (None, 1.0):
            params += " refinement.parameterisation.block_width={}".format(
                self.blockWidth.get())

        if self.setScanVaryingErrors.get() not in (None, False):
            params += " refinement.parameterisation.set_scan_varying_errors={}".format(
                self.setScanVaryingErrors.get())

        if self.debugCentroidAnalysis.get() not in (None, False):
            params += " refinement.parameterisation.debug_centroid_analysis={}".format(
                self.debugCentroidAnalysis.get())

        # FIXME: make all beamfix one block
        beamfix = []
        if self.beamFixInSpindlePlane.get() is True:
            beamfix += "in_spindle_plane,"
        if self.beamFixOutSpindlePlane.get() is True:
            beamfix += "out_spindle_plane,"
        if self.beamFixWavelength.get() is True:
            beamfix += "wavelength"
        params += " refinement.parameterisation.beamFixInSpindlePlane={}".format(
            beamfix)

        if self.beamForceStatic.get() not in (None, True):
            params += " refinement.parameterisation.beam.force_static={}".format(
                self.beamForceStatic.get())

        if self.SmootherIntervalWidthDegrees.get() not in (None, 36.0):
            params += " refinement.parameterisation.Smoother.interval_width_degrees={}".format(
                self.SmootherIntervalWidthDegrees.get())

        if self.SmootherAbsoluteNumIntervals.get() is not None:
            params += " refinement.parameterisation.Smoother.absolute_num_intervals={}".format(
                self.SmootherAbsoluteNumIntervals.get())

        # FIXME: Combine in one line
        if self.crystalFixCell.get() not in (None, True):
            params += " refinement.parameterisation.crystal.fix.cell={}".format(
                self.crystalFixCell.get())

        if self.crystalFixOrientation.get() not in (None, True):
            params += " refinement.parameterisation.crystal.fix.orientation={}".format(
                self.crystalFixOrientation.get())

        if self.crystalUnitCellForceStatic.get() not in (None, True):
            params += " refinement.parameterisation.crystal.unit_cell.force_static={}".format(
                self.crystalUnitCellForceStatic.get())

        if self.crystalOrientationForceStatic.get() not in (None, True):
            params += " refinement.parameterisation.crystal.orientation.force_static={}".format(
                self.crystalOrientationForceStatic.get())

        if self.detectorPanels.get() is SINGLE:
            params += " refinement.parameterisation.detector.panels={}".format(
                'single')
        elif self.detectorPanels.get() is MULTIPLE:
            params += " refinement.parameterisation.detector.panels={}".format(
                'multiple')
        elif self.detectorPanels.get() is HIERARCHICAL:
            params += " refinement.parameterisation.detector.panels={}".format(
                'hierarchical')

        if self.detectorHierarchyLevel.get() not in (None, 0):
            params += " refinement.parameterisation.detector.hierarchy_level={}".format(
                self.detectorHierarchyLevel.get())

        # FIXME: Convert to one line
        if self.detectorFixPosition.get() not in (None, True):
            params += " refinement.parameterisation.detector.fix.position={}".format(
                self.detectorFixPosition.get())

        if self.detectorFixOrientation.get() not in (None, True):
            params += " refinement.parameterisation.detector.fix.orientation={}".format(
                self.detectorFixOrientation.get())

        if self.detectorFixdistance.get() not in (None, True):
            params += " refinement.parameterisation.detector.fix.distance={}".format(
                self.detectorFixdistance.get())

        if self.crystalForceStatic.get() not in (None, True):
            params += " refinement.parameterisation.crystal.force_static={}".format(
                self.crystalForceStatic.get())

        # FIXME
        if self.goniometerFixInBeamPlane.get() not in (None, True):
            params += " refinement.parameterisation.goniometer.fix.in_beam_plane={}".format(
                self.goniometerFixInBeamPlane.get())

        if self.goniometerFixOutBeamPlane.get() not in (None, True):
            params += " refinement.parameterisation.goniometer.fix.out_beam_plane={}".format(
                self.goniometerFixOutBeamPlane.get())

        if self.goniometerForceStatic.get() not in (None, True):
            params += " refinement.parameterisation.goniometer.force_static={}".format(
                self.goniometerForceStatic.get())

        if self.treatSingleImageAsStill.get() not in (None, False):
            params += " refinement.parameterisation.treat_single_image_as_still={}".format(
                self.treatSingleImageAsStill.get())

        if self.sphericalRelpModel.get() not in (None, False):
            params += " refinement.parameterisation.spherical_relp_model={}".format(
                self.sphericalRelpModel.get())

        if self.refineryEngine.get() is SIMPLE_LBFGS:
            params += " refinery.engine={}".format('SimpleLBFGS')
        elif self.refineryEngine.get() is LBFGS_CURVS:
            params += " refinery.engine={}".format('LBFGScurvs')
        elif self.refineryEngine.get() is GAUSS_NEWTON:
            params += " refinery.engine={}".format('GaussNewton')
        elif self.refineryEngine.get() is SPARSE_LEV_MAR:
            params += " refinery.engine={}".format('SparseLevMar')

        if self.refineryMaxIterations.get() is not None:
            params += " refinery.max_iterations={}".format(
                self.refineryMaxIterations.get())

        if self.refineryLog.get() is not None:
            params += " refinery.log={}".format(self.refineryLog.get())

        if self.refineryJournalTrackStep.get() not in (None, False):
            params += " refinery.Journal..track_step={}".format(
                self.refineryJournalTrackStep.get())

        if self.refineryJournalTrackGradient.get() not in (None, False):
            params += " refinery.Journal.track_gradient={}".format(
                self.refineryJournalTrackGradient.get())

        if self.refineryJournalTrackParameterCorrelation.get() not in (None, False):
            params += " refinery.Journal.track_parameter_correlation={}".format(
                self.refineryJournalTrackParameterCorrelation.get())

        if self.refineryJournalTrackConditionNumber.get() not in (None, False):
            params += " refinery.Journal.track_condition_number={}".format(
                self.refineryJournalTrackConditionNumber.get())

        if self.refineryJournalTrackOutOfSampleRmsd.get() not in (None, False):
            params += " refinery.Journal.track_out_of_sample_rmsd={}".format(
                self.refineryJournalTrackOutOfSampleRmsd.get())

        if self.targetRmsdCutoff.get() is ABSOLUTE:
            params += " target.rmsd_cutoff={}".format('absolute')

        if self.targetBinSizeFraction.get() not in (None, 0.0):
            params += " target.bin_size_fraction={}".format(
                self.targetBinSizeFraction.get())

        if self.targetAbsoluteCutoffs.get() is not None:
            params += " target.absolute_cutoffs={}".format(
                self.targetAbsoluteCutoffs.get())

        if self.targetGradientCalculationBlockSize.get() not in (None, True):
            params += " target.gradient_calculation_block_size={}".format(
                self.targetGradientCalculationBlockSize.get())

        if self.reflectionsPerDegree.get() is not None:
            params += " reflections.reflections_per_degree={}".format(
                self.reflectionsPerDegree.get())

        if self.minimumSampleSize.get() not in (None, 1000):
            params += " reflections.minimum_sample_size={}".format(
                self.minimumSampleSize.get())

        if self.maximumSampleSize.get() is not None:
            params += " reflections.maximum_sample_size={}".format(
                self.maximumSampleSize.get())

        if self.randomSeed.get() is not 42:
            params += " reflections.random_seed={}".format(
                self.randomSeed.get())

        if self.closeToSpindleCutoff.get() not in (None, 0.02):
            params += " reflections.close_to_spindle_cutoff={}".format(
                self.closeToSpindleCutoff.get())

        if self.trimScanEdges.get() not in (None, 0.0):
            params += " reflections.trim_scan_edges={}".format(
                self.trimScanEdges.get())

        if self.weightingStratOverride.get() is STATISTICAL:
            params += " reflections.weighting_strategy.override={}".format(
                'statistical')
        elif self.weightingStratOverride.get() is STILLS:
            params += " reflections.weighting_strategy.override={}".format(
                'stills')
        elif self.weightingStratOverride.get() is CONSTANT:
            params += " reflections.weighting_strategy.override={}".format(
                'constant')
        elif self.weightingStratOverride.get() is EXTERNAL_DELTAPSI:
            params += " reflections.weighting_strategy.override={}".format(
                'external_deltapsi')

        if self.delpsiConstant.get() not in (None, 1000000):
            params += " reflections.weighting_strategy.delpsi_constant={}".format(
                self.delpsiConstant.get())

        if self.weightingStrategyConstants.get() not in (None, 1.0):
            params += " reflections.weighting_strategy.constants={} {} {}".format(self.weightingStrategyConstants.get(
            ), self.weightingStrategyConstants.get(), self.weightingStrategyConstants.get())

        if self.outlierAlgorithm.get() is NULL:
            params += " reflections.outlier.algorithm={}".format('null')
        elif self.outlierAlgorithm.get() is MCD:
            params += " reflections.outlier.algorithm={}".format('mcd')
        elif self.outlierAlgorithm.get() is TUKEY:
            params += " reflections.outlier.algorithm={}".format('tukey')
        elif self.outlierAlgorithm.get() is SAUTER_POON:
            params += " reflections.outlier.algorithm={}".format('sauter_poon')

        if self.outlierMinimumNumberOfReflections.get() not in (None, 20):
            params += " reflections.outlier.minimum_number_of_reflections={}".format(
                self.outlierMinimumNumberOfReflections.get())

        if self.outlierSeparateExperiments.get() not in (None, True):
            params += " reflections.outlier.separate_experiments={}".format(
                self.outlierSeparateExperiments.get())

        if self.outlierSeparatePanels.get() is not None:
            params += " reflections.outlier.separate_panels={}".format(
                self.outlierSeparatePanels.get())

        if self.outlierSeparateBlocks.get() not in (None, True):
            params += " reflections.outlier.separate_blocks={}".format(
                self.outlierSeparateBlocks.get())

        if self.outlierBlockWidth.get() is not None:
            params += " reflections.outlier.block_width={}".format(
                self.outlierBlockWidth.get())

        if self.outlierTukeyIqrMultiplier.get() not in (None, 1.5):
            params += " reflections.outlier.tukey.iqr_multiplier={}".format(
                self.outlierTukeyIqrMultiplier.get())

        if self.mcdAlpha.get() not in (None, 0.5):
            params += " reflections.mcd.alpha={}".format(self.mcdAlpha.get())

        if self.mcdMaxNGroups.get() not in (None, 5):
            params += " reflections.mcd.max_n_groups={}".format(
                self.mcdMaxNGroups.get())

        if self.mcdMinGroupSize.get() not in (None, 300):
            params += " reflections.mcd.min_group_size={}".format(
                self.mcdMinGroupSize.get())

        if self.mcdNTrials.get() not in (None, 500):
            params += " reflections.mcd.n_trials={}".format(
                self.mcdNTrials.get())

        if self.mcdK1.get() not in (None, 2):
            params += " reflections.mcd.k1={}".format(self.mcdK1.get())

        if self.mcdK2.get() not in (None, 2):
            params += " reflections.mcd.k2={}".format(self.mcdK2.get())

        if self.mcdK3.get() not in (None, 100):
            params += " reflections.mcd.k3={}".format(self.mcdK3.get())

        if self.mcdThresholdProbability.get() not in (None, 0.975):
            params += " reflections.mcd.threshold_probability={}".format(
                self.mcdThresholdProbability.get())

        if self.sauterPoonPxSz.get() is not None:
            params += " reflections.sauter_poon.px_sz={}".format(
                self.sauterPoonPxSz.get())

        if self.sauterPoonVerbose.get() not in (None, False):
            params += " reflections.sauter_poon.verbose={}".format(
                self.sauterPoonVerbose.get())

        if self.sauterPoonPdf.get() not in (None, ''):
            params += " reflections.sauter_poon.pdf={}".format(
                self.sauterPoonPdf.get())

        return params

    def _prepBravaisCommandline(self, program):
        "Create the command line input to run dials programs"
        # Input basic parameters
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        params = "{} {} output.log={}".format(
            self.getModelFile(), self.getReflFile(), logPath)

        # Update the command line with additional parameters

        if self.refineBravNproc.get() not in (None, 4):
            params += " nproc={}".format(self.refineBravNproc.get())
        return params

    def _prepReindexCommandline(self, program):
        "Create the command line input to run dials programs"
        # Input basic parameters
        logPath = "{}/{}.log".format(self._getLogsPath(), program)
        params = "change_of_basis_op={} output.log={}".format(
            self.getChangeOfBasisOp(), logPath)

        if self.doReindexModel.get():
            params += " {}".format(self.getModelFile())

        if self.doReindexReflections.get():
            params += " {}".format(self.getReflFile())
        return params

    def _createScanRanges(self):
        # Go through the
        images = [image.getObjId() for image in self.inputImages.get()
                  if image.getIgnore() is not True]
        scanranges = find_subranges(images)
        scanrange = ' '.join('spotfinder.scan_range={},{}'.format(i, j)
                             for i, j in scanranges)
        return scanrange
