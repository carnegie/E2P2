from argparse import ArgumentTypeError
import os

# Directory Set Up
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(os.path.dirname(ROOT_DIR), 'data')
MAPS_DIR = os.path.join(DATA_DIR, 'maps')
WEIGHTS_DIR = os.path.join(DATA_DIR, 'weights')

# Mapping Files
EF_CLASS_MAP = os.path.join(MAPS_DIR, 'efclasses.mapping')
EC_SUPERSEDED_MAP = os.path.join(MAPS_DIR, 'pf-EC-superseded.mapping')
METACYC_RXN_MAP = os.path.join(MAPS_DIR, 'pf-metacyc-RXN-EC.mapping')
OFFICIAL_EC_METACYC_RXN_MAP = os.path.join(MAPS_DIR, 'pf-official-EC-metacyc-RXN.mapping')
TO_REMOVE_NON_SMALL_MOLECULE_METABOLISM = os.path.join(MAPS_DIR, 'pf-to-remove-non-small-molecule-metabolism.mapping')

# Weight Files
BLAST_WEIGHT = os.path.join(WEIGHTS_DIR, 'blast')
PRIAM_WEIGHT = os.path.join(WEIGHTS_DIR, 'priam')

# Default Software Commands and Arguments
# BLASTP_CMD = "blastp"
# JAVA_CMD = "java"
# BLASTP_BIN_PATH = "/usr/local/bin/"
BLASTP_NUM_THREADS = "1"

# Other Default Settings
DEFAULT_LOGGER_LEVEL = "DEBUG"
DEFAULT_LOGGER_NAME = "e2p2"
DEFAULT_ENSEMBLE_METHOD = "Maximum weight with absolute threshold"
DEFAULT_THRESHOLD = float(0.5)
DEFAULT_OUTPUT_SUFFIX = ".e2p2v4"
DEFAULT_LONG_OUTPUT_SUFFIX = ".long"
DEFAULT_PF_OUTPUT_SUFFIX = ".pf"
DEFAULT_ORXN_PF_OUTPUT_SUFFIX = ".orxn.pf"
DEFAULT_FINAL_PF_OUTPUT_SUFFIX = ".final.pf"
DEFAULT_PTOOLS_CHAR_LIMIT = 40

# Website Default
BLAST_PLUS_DOWNLOAD_LINK = "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
JAVA_8_DOWNLOAD_LINK = "https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html"
PRIAM_SEARCH_LINK = "http://priam.prabi.fr/utilities/PRIAM_search.jar"


# Define specific input types for ArgumentParser
class PathType(object):
	def __init__(self, type='file'):
		'''
		type: file, dir, parent, 'blastdb', 'blastbin', 'profiles', None, or a function returning True for valid paths
		'''

		assert type in ('file', 'dir', 'parent', 'blastdb', 'blastbin', 'profiles', None) or hasattr(type, '__call__')
		self._type = type

	def __call__(self, string):
		path = os.path.realpath(string)
		if self._type is None:
			pass
		elif self._type == 'file':
			if not os.path.isfile(path):
				raise ArgumentTypeError("Path is not a file: '%s'" % string)
		elif self._type == 'dir':
			if not os.path.isdir(path):
				raise ArgumentTypeError("Path is not a directory: '%s'" % string)
		elif self._type == 'parent':
			if not os.path.isdir(os.path.dirname(path)):
				raise ArgumentTypeError("Parent of path is not a directory: '%s'" % string)
		elif self._type == 'blastdb':
			phr_path = path + '.phr'
			pin_path = path + '.pin'
			psq_path = path + '.psq'
			if not os.path.isfile(phr_path) or not os.path.isfile(pin_path) or not os.path.isfile(psq_path):
				raise ArgumentTypeError("Cannot find blast database at path: '%s'" % string)
		elif self._type == 'blastbin':
			bin_files = os.listdir(path)
			if not os.path.isdir(path) or 'makeprofiledb' not in bin_files or \
					'rpsblast' not in bin_files or 'rpstblastn' not in bin_files:
				raise ArgumentTypeError("Path not a valid Blast+ bin folder: '%s'" % string)
		elif self._type == 'profiles':
			bin_files = os.listdir(path)
			if not os.path.isdir(path) or 'PROFILES' not in bin_files or \
					'annotation_rules.xml' not in bin_files or 'genome_rules.xml' not in bin_files:
				raise ArgumentTypeError("Path not a valid Priam profiles folder: '%s'" % string)
			else:
				profiles_path = os.path.join(path, 'PROFILES')
				annotation_rules = os.path.join(path, 'annotation_rules.xml')
				genome_rules = os.path.join(path, 'genome_rules.xml')
				if not os.path.isdir(profiles_path) or not os.path.isfile(annotation_rules) \
					or not os.path.isfile(genome_rules):
					raise ArgumentTypeError("Path not a valid Priam profiles folder: '%s'" % string)
		else:
			raise ArgumentTypeError("path not valid: '%s'" % string)
		return string
