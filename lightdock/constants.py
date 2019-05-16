"""LightDock default parameters"""

# GSO algorithm constants
MAX_TRANSLATION = 30                # In Angstroms
MAX_ROTATION = 1.0                  # Quaternion default value for its components
DEFAULT_STEP_SIZE = 0.03            # Default generic GSO step (does only apply to J* functions
DEFAULT_TRANSLATION_STEP = 0.5      # Normalized step
DEFAULT_ROTATION_STEP = 0.5         # Normalized SLERP step. 1 means full jump, 0 means no movement
GSO_SEED = 324324                   # Seed for the random number generator in the GSO algorithm
STARTING_POINTS_SEED = 324324       # Seed for the random number generator used for calculating starting points

# Normal modes
DEFAULT_NMODES_REC = 10              # Default number of normal modes to consider for receptor
DEFAULT_NMODES_LIG = 10              # Default number of normal modes to consider for ligand
DEFAULT_EXTENT_MU = 4.0
DEFAULT_EXTENT_SIGMA = 3.0
MIN_EXTENT = 0.1
MAX_EXTENT = 5.0
STARTING_NM_SEED = 324324           # Seed for the random number generator used for calculating normal modes extent
DEFAULT_NMODES_STEP = 0.5
DEFAULT_REC_NM_FILE = "lightdock_rec.nm"
DEFAULT_LIG_NM_FILE = "lightdock_lig.nm"

# Scoring function constants
DEFAULT_SCORING_FUNCTION = "fastdfire"  # Lightdock default scoring function if none is specified
DEFAULT_CONTACT_RESTRAINTS_CUTOFF = 3.9

# Default file extensions
DEFAULT_REFERENCE_POINTS_EXTENSION = ".vol"
DEFAULT_LIST_EXTENSION = ".list"
DEFAULT_RMSD_EXTENSION= ".rmsd"
NUMPY_FILE_SAVE_EXTENSION = ".npy"
DEFAULT_REPRESENTATIVES_EXTENSION = ".repr"

# Default file names and folders
DEFAULT_POSITIONS_FOLDER = "init"       # Folder which contains the initial_positions files for each cluster
GSO_OUTPUT_FILE = "gso_%s.out"
DEFAULT_SWARM_FOLDER = "swarm_"     # Folder where GSO execution for a given cluster will be stored
DEFAULT_SETUP_FILE = "setup.json"
DEFAULT_PDB_STARTING_PREFIX = "starting_positions"
DEFAULT_BILD_STARTING_PREFIX = "starting_poses"
DEFAULT_STARTING_PREFIX = "initial_positions"
CLUSTERS_CENTERS_FILE = "cluster_centers.pdb"
RANKING_FILE = "solutions.list"
RANKING_BY_RMSD_FILE = "rank_by_rmsd.list"
RANKING_BY_LUCIFERIN_FILE = "rank_by_luciferin.list"
RANKING_BY_SCORING_FILE = "rank_by_scoring.list"
DEFAULT_LIGHTDOCK_PREFIX = "lightdock_%s"
DEFAULT_ELLIPSOID_DATA_EXTENSION = ".xyz"
EVALUATION_FILE = "evaluation.list"
SCORING_FILE = "scoring.list"
LIGHTDOCK_PDB_FILE = "lightdock_%s.pdb"
CLUSTER_DEFAULT_NAME = "cluster"
CLUSTER_ANALYSIS_FILE = CLUSTER_DEFAULT_NAME + ".repr"
CLUSTER_REPRESENTATIVES_FILE = CLUSTER_DEFAULT_NAME + DEFAULT_REPRESENTATIVES_EXTENSION
DEFAULT_LIGHTDOCK_INFO = "lightdock.info"

# Surface density
MIN_SURFACE_DENSITY = 200.0

