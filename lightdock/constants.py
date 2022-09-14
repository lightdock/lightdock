"""LightDock default parameters"""

# Simulation defaults
DEFAULT_NUM_SWARMS = 400
DEFAULT_NUM_GLOWWORMS = 200

"""GSO algorithm constants"""
# Maximum translation allowed at initialization (in Angstroms)
MAX_TRANSLATION = 30
# Quaternion default value for its components
MAX_ROTATION = 1.0
# Default generic GSO step (does only apply to J* functions
DEFAULT_STEP_SIZE = 0.03
# Interpolation step for translation (in %)
DEFAULT_TRANSLATION_STEP = 0.5
# Normalized SLERP step. 1 means full jump, 0 means no movement
DEFAULT_ROTATION_STEP = 0.5
# Seed for the random number generator in the GSO algorithm
GSO_SEED = 324324
# Seed for the random number generator used for calculating starting points
STARTING_POINTS_SEED = 324324

"""Normal modes"""
# Default number of normal modes to consider for receptor
DEFAULT_NMODES_REC = 10
# Default number of normal modes to consider for ligand
DEFAULT_NMODES_LIG = 10
# Default RMSD interval at which conformations will be sampled through ANM
DEFAULT_ANM_RMSD = 0.5
# Default mean and std for the Normal distribution for random ANM extents
DEFAULT_EXTENT_MU = 0.0
DEFAULT_EXTENT_SIGMA = 1.0
# Boundaries, this does not affect to modes generation
MIN_EXTENT = 0.0
MAX_EXTENT = 1.0
# Seed for the random number generator used for calculating normal modes extent
STARTING_NM_SEED = 324324
# When interpolating, which
DEFAULT_NMODES_STEP = 0.5
DEFAULT_REC_NM_FILE = "lightdock_rec.nm"
DEFAULT_LIG_NM_FILE = "lightdock_lig.nm"

"""Scoring function constants"""
# LightDock default scoring function if none is specified
DEFAULT_SCORING_FUNCTION = "fastdfire"
# Threshold for an atomic-pair to be considered a contact (in Angstroms)
DEFAULT_ATOMIC_CONTACT = 3.9
# Threshold for an atomic-pair to be considered a contact in restraints definition (in Angstroms)
DEFAULT_CONTACT_RESTRAINTS_CUTOFF = 3.9

"""Default file extensions"""
DEFAULT_REFERENCE_POINTS_EXTENSION = ".vol"
DEFAULT_LIST_EXTENSION = ".list"
DEFAULT_RMSD_EXTENSION = ".rmsd"
NUMPY_FILE_SAVE_EXTENSION = ".npy"
DEFAULT_REPRESENTATIVES_EXTENSION = ".repr"
DEFAULT_ELLIPSOID_DATA_EXTENSION = ".xyz"

"""Default file names and folders"""
# Folder which contains the initial_positions files for each swarm
DEFAULT_POSITIONS_FOLDER = "init"
# Simulation default output file
GSO_OUTPUT_FILE = "gso_%s.out"
# Folder where GSO execution for a given swarm will be stored
DEFAULT_SWARM_FOLDER = "swarm_"
# Stores simulation step information
DEFAULT_SETUP_FILE = "setup.json"
DEFAULT_PDB_STARTING_PREFIX = "starting_positions"
DEFAULT_BILD_STARTING_PREFIX = "starting_poses"
DEFAULT_STARTING_PREFIX = "initial_positions"
# A special PDB file containing all swarm centers over the receptor surface
# centered at origin of coordinates (use together with LightDock receptor parsed PDB)
SWARM_CENTERS_FILE = "swarm_centers.pdb"

RANKING_FILE = "solutions.list"
RANKING_BY_RMSD_FILE = "rank_by_rmsd.list"
RANKING_BY_LUCIFERIN_FILE = "rank_by_luciferin.list"
RANKING_BY_SCORING_FILE = "rank_by_scoring.list"
DEFAULT_LIGHTDOCK_PREFIX = "lightdock_%s"
EVALUATION_FILE = "evaluation.list"
SCORING_FILE = "scoring.list"
LIGHTDOCK_PDB_FILE = "lightdock_%s.pdb"
CLUSTER_DEFAULT_NAME = "cluster"
CLUSTER_REPRESENTATIVES_FILE = CLUSTER_DEFAULT_NAME + DEFAULT_REPRESENTATIVES_EXTENSION
# Each independent simulation generates a new file:
DEFAULT_LIGHTDOCK_INFO = "lightdock.info"
DEFAULT_MASK_FILE = "lightdock_%s_mask" + NUMPY_FILE_SAVE_EXTENSION

"""Swarm calculations"""
# Total SASA will be divided by this number
DEFAULT_SURFACE_DENSITY = 50.0
# Swarm radius for the boundary of possible translations of initial poses
DEFAULT_SWARM_RADIUS = 10.0
# If this value is 0.0, swarm distance is estimated using the ligand max diameter
DEFAULT_SWARM_DISTANCE = 0.0
# Number of swarms to keep when applying restraints
DEFAULT_SWARMS_PER_RESTRAINT = 20
# Number of sampling points per centroid in automatic swarm calculation
DEFAULT_SPHERES_PER_CENTROID = 100
# Cutoff for filtering swarms too close to the surface
SWARM_DISTANCE_TO_SURFACE_CUTOFF = 3.0
