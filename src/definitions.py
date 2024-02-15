import os


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DEFAULT_CONFIG_PATH = os.path.join(ROOT_DIR, 'config.ini')
SRC_DIR = os.path.join(ROOT_DIR, 'src')
E2P2_CLS_DIR = os.path.join(SRC_DIR, 'e2p2')
CLASSIFIERS_CLS_DIR = os.path.join(E2P2_CLS_DIR, 'classifiers')
ENSEMBLES_CLS_DIR = os.path.join(E2P2_CLS_DIR, 'ensembles')
DATA_DIR = os.path.join(ROOT_DIR, 'data')
MAPS_DIR = os.path.join(DATA_DIR, 'maps')
WEIGHTS_DIR = os.path.join(DATA_DIR, 'weights')

DEFAULT_BLAST_E_VALUE = float("1e-2")
DEFAULT_BLAST_BIT_SCORE = float("0")
DEFAULT_PRIAM_E_VALUE = float("1e-2")
DEEPEC_DIR = os.path.join(ROOT_DIR, 'deepec')
EC_TO_EF_MAPPING_PATH = os.path.join(DEEPEC_DIR, 'deepec/data/ec_to_ef.mapping')


DEFAULT_LOGGER_LEVEL = "DEBUG"
DEFAULT_LOGGER_NAME = "e2p2"
DEFAULT_OUTPUT_SUFFIX = "e2p2"
DEFAULT_LONG_OUTPUT_SUFFIX = "long"
DEFAULT_PF_OUTPUT_SUFFIX = "default.pf"
DEFAULT_ORXN_PF_OUTPUT_SUFFIX = "orxn.pf"
DEFAULT_FINAL_PF_OUTPUT_SUFFIX = "final.pf"
DEFAULT_PTOOLS_CHAR_LIMIT = 40

# Website Default
BLAST_PLUS_DOWNLOAD_LINK = "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
JAVA_8_DOWNLOAD_LINK = "https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html"
PRIAM_SEARCH_LINK = "http://priam.prabi.fr/utilities/PRIAM_search.jar"

