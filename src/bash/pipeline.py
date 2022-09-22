import errno
import os
import textwrap
import time

from src.definitions import DEFAULT_DIR, DEFAULT_LOGGER_NAME, DEFAULT_LOGGER_LEVEL, DEFAULT_OUTPUT_SUFFIX
from src.lib.util import PathType, check_fasta_header, logging_helper, remove_splice_variants_from_fasta


def add_io_arguments(argument_parser):
    argument_parser.add_argument("--config", "-c", dest="config_file", type=PathType('file'),
                                 help="Path to config file", default=os.path.join(DEFAULT_DIR, "config.ini"))
    argument_parser.add_argument("--input", "-i", dest="input_file", type=PathType('file'),
                                 help="Path to input protein sequences file", required=True)
    argument_parser.add_argument("--protein_gene", "-pg", dest="protein_gene_path", type=PathType('file'),
                                 help="Provide a protein to gene map. This can be used to generate a "
                                      "splice variant removed fasta file and output the final version of e2p2.")
    argument_parser.add_argument("--remove_splice_variants", "-rm", dest="remove_splice_variants", action="store_true",
                                 help="Argument flag to remove splice variants, use it if yes.")
    argument_parser.add_argument("--output", "-o", dest="output_path", type=PathType('have_parent'),
                                 help="Path to output file. By Default would be in the same folder of the input.")
    argument_parser.add_argument("--temp_folder", "-tf", dest="temp_folder", type=PathType('dir'),
                                 help="Specify the location of the temp folder. "
                                      "By default would be in the same directory of the output.")
    argument_parser.add_argument("--log", "-l", dest="log_path", type=PathType('have_parent'),
                                 help="Specify the location of the log file. "
                                      "By default would be \"runE2P2.log\" in the temp folder.")
    verbose_message = '''Verbose level of log output. Default is 0.
            0: only step information are logged
            1: all information are logged
            '''
    argument_parser.add_argument("--verbose", "-v", dest="verbose", default="0", choices=["0", "1"],
                                 help=textwrap.dedent(verbose_message))


def add_classifier_arugments(argument_parser):
    # Argument for E2P2 classifiers
    argument_parser.add_argument("--blastp", "-b", dest="blastp_cmd",
                                 help=textwrap.dedent("Command of or path to BLAST+ \"blastp\"."))
    argument_parser.add_argument("--num_threads", "-n", dest="num_threads", type=int,
                                 help="Number of threads to run \"blastp\".")
    argument_parser.add_argument("--java", "-j", dest="java_cmd",
                                 help=textwrap.dedent("Command of or path to \"java\"."))
    argument_parser.add_argument("--priam_search", "-ps", dest="priam_search", type=PathType('file'),
                                 help=textwrap.dedent("Path to \"PRIAM_search.jar\"."))
    argument_parser.add_argument("--priam_resume", "-pr", dest="priam_resume", action='store_true',
                                 help="Whether or not to resume a found PRIAM_search.jar process.")
    argument_parser.add_argument("--blast_bin", "-bb", dest="blast_bin", type=PathType('blast_bin'),
                                 help=textwrap.dedent("Command of or path to BLAST+ bin folder."))

    # Arguments for E2P2 databases
    argument_parser.add_argument("--blast_db", "-bd", dest="blast_db", type=PathType('blast_db'),
                                 help=textwrap.dedent(
                                     "Path to rpsd blast database name.\nFor example, \"/PATH/TO/FOLDER/rpsd.fa\", "
                                     "where you can find the following files in /PATH/TO/FOLDER:rpsd.fa.phr; "
                                     "rpsd.fa.pin; rpsd.fa.psq"))
    argument_parser.add_argument("--blast_evalue", "-be", dest="blast_evalue", type=float,
                                 help=textwrap.dedent("Blastp e-value cutoff"))
    argument_parser.add_argument("--blast_weight", "-bw", dest="blast_weight", type=PathType('file'),
                                 help=textwrap.dedent("Path to weight file for the blast classifier"))
    argument_parser.add_argument("--priam_profiles", "-pp", dest="priam_profiles", type=PathType('priam_profiles'),
                                 help=textwrap.dedent(
                                     "Path to PRIAM profiles.\nFor example, \"/PATH/TO/FOLDER/profiles\", "
                                     "where you can find the following in /PATH/TO/FOLDER/profiles:\n "
                                     "files: annotation_rules.xml; genome_rules.xml\n "
                                     "folders: PROFILES: Folder contains \"LIBRARY\" folder and "
                                     "multiple \".chk\" files."))
    argument_parser.add_argument("--priam_weight", "-pw", dest="priam_weight", type=PathType('file'),
                                 help=textwrap.dedent("Path to weight file for the priam classifier"))

    # Arguments for E2P2 classifier modules
    argument_parser.add_argument("--blast_cls", "-bc", dest="blast_cls", type=PathType('file'),
                                 help=textwrap.dedent("Path to the BLAST class module."))
    argument_parser.add_argument("--priam_cls", "-pc", dest="priam_cls", type=PathType('file'),
                                 help=textwrap.dedent("Path to the PRIAM class module."))


def add_ensemble_arguments(argument_parser):
    # Arguments for E2P2 ensembles
    argument_parser.add_argument("--ensemble_cls", "-ec", dest="ensemble_cls", type=PathType('file'),
                                 help="Path to E2P2 ensemble class module.")
    argument_parser.add_argument("--threshold", "-t", dest="threshold", type=float,
                                 help="Threshold for voting results. Default is 0.5.")


def add_mapping_arguments(argument_parser):
    # Arguments for maps
    argument_parser.add_argument("--ef_map", "-ef", dest="ef_map", type=PathType('file'),
                                 help="Path to efclasses.mapping file.")
    argument_parser.add_argument("--ec_superseded", "-es", dest="ec_superseded", type=PathType('file'),
                                 help="Path to EC-superseded.mapping file.")
    argument_parser.add_argument("--rxn_ec", "-me", dest="metacyc_rxn_ec", type=PathType('file'),
                                 help="Path to metacyc-RXN-EC.mapping file.")
    argument_parser.add_argument("--official_ec_rxn", "-oer", dest="official_ec_metacyc_rxn", type=PathType('file'),
                                 help="Path to official-EC-metacyc-RXN.mapping file.")
    argument_parser.add_argument("--to_remove", "-tr", dest="to_remove_metabolism", type=PathType('file'),
                                 help="Path to to-remove-non-small-molecule-metabolism.mapping file.")


def io_helper(input_file, logger_name=DEFAULT_LOGGER_NAME, output_path=None, timestamp=str(time.time()),
              temp_folder=None, log_path=None, verbose="0"):
    check_fasta_header(input_file, logger_name)
    # Setup output paths
    if output_path is None:
        output_path = '.'.join([input_file, DEFAULT_OUTPUT_SUFFIX])
    # Setup temp folder path
    if temp_folder is None:
        temp_folder = input_file + '.temp.' + timestamp

    create_temp_folder_flag = False
    try:
        os.mkdir(temp_folder)
        create_temp_folder_flag = True
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass
    # Setup logging file path
    if log_path is None:
        log_path = os.path.join(temp_folder, '.'.join([DEFAULT_LOGGER_NAME, timestamp, 'log']))
    if verbose == "0":
        logging_level = "DEBUG"
    else:
        logging_level = "INFO"

    io_dict = {"IO": {"query": input_file, "out": temp_folder, "timestamp": timestamp}}

    return output_path, io_dict, create_temp_folder_flag, log_path, logging_level


def protein_to_gene_helper(input_file, output_path, protein_gene_path, remove_splice_variants,
                           logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
    output_folder = os.path.dirname(output_path)
    if protein_gene_path is not None and remove_splice_variants is True:
        return remove_splice_variants_from_fasta(input_file, output_folder, protein_gene_path, logger_name=logger_name)
    elif protein_gene_path is not None and remove_splice_variants is False:
        logging_helper("Protein to gene map not used to remove splice variants.", logging_level="DEBUG",
                       logger_name=logger_name)
        return input_file
    elif protein_gene_path is None and remove_splice_variants is True:
        logging_helper("Cannot remove splice variants without protein to gene map.", logging_level="WARNING",
                       logger_name=logger_name)
        return input_file
    else:
        return input_file


def blast_overwrites(blastp_cmd, blast_db, num_threads, blast_evalue, blast_weight, blast_cls, overwrites=None):
    if overwrites is None:
        overwrites = {}
    if len([i for i in [blastp_cmd, blast_db, num_threads, blast_evalue, blast_weight, blast_cls]
            if i is not None]) > 0:
        overwrites.setdefault("BLAST", {})
        if blastp_cmd is not None:
            overwrites["BLAST"].setdefault("blastp", blastp_cmd)
        if blast_db is not None:
            overwrites["BLAST"].setdefault("blastdb", blast_db)
        if num_threads is not None:
            overwrites["BLAST"].setdefault("num_threads", str(num_threads))
        if blast_evalue is not None:
            overwrites["BLAST"].setdefault("evalue_threshold", str(blast_evalue))
        if blast_weight is not None:
            overwrites["BLAST"].setdefault("weight", blast_weight)
        if blast_cls is not None:
            overwrites["BLAST"].setdefault("class", blast_cls)
    return overwrites


def priam_overwrites(java_cmd, priam_search, priam_resume, blast_bin, priam_profiles, priam_weight, priam_cls,
                     overwrites=None):
    if overwrites is None:
        overwrites = {}
    if len([i for i in [java_cmd, priam_search, priam_resume, blast_bin, priam_profiles, priam_weight, priam_cls]
            if i is not None]) > 0:
        overwrites.setdefault("PRIAM", {})
        if java_cmd is not None:
            overwrites["PRIAM"].setdefault("java_path", java_cmd)
        if priam_search is not None:
            overwrites["PRIAM"].setdefault("priam_search_path", priam_search)
        if priam_resume is not None:
            if priam_resume is True:
                overwrites["PRIAM"].setdefault("resume", 'fr')
            else:
                overwrites["PRIAM"].setdefault("resume", 'fn')
        if blast_bin is not None:
            overwrites["PRIAM"].setdefault("path_to_blast_plus_bin", blast_bin)
        if priam_profiles is not None:
            overwrites["PRIAM"].setdefault("path_to_priam_profiles", priam_profiles)
        if priam_weight is not None:
            overwrites["PRIAM"].setdefault("weight", priam_weight)
        if priam_cls is not None:
            overwrites["PRIAM"].setdefault("class", priam_cls)
    return overwrites


def ensemble_overwrites(ensemble_cls, threshold, overwrites=None):
    if overwrites is None:
        overwrites = {}
    if len([i for i in [ensemble_cls, threshold] if i is not None]) > 0:
        overwrites.setdefault("MaxWeightAbsoluteThreshold", {})
        if ensemble_cls is not None:
            overwrites["MaxWeightAbsoluteThreshold"].setdefault("class", ensemble_cls)
        if threshold is not None:
            overwrites["MaxWeightAbsoluteThreshold"].setdefault("threshold", str(threshold))
    return overwrites


def mapping_overwrite(ef_map, ec_superseded, metacyc_rxn_ec, official_ec_metacyc_rxn, to_remove_metabolism,
                      overwrites=None):
    if overwrites is None:
        overwrites = {}
    if len([i for i in [ef_map, ec_superseded, metacyc_rxn_ec, official_ec_metacyc_rxn, to_remove_metabolism]
            if i is not None]) > 0:
        overwrites.setdefault("Mapping", {})
        if ef_map is not None:
            overwrites["Mapping"].setdefault("efclasses", ef_map)
        if ec_superseded is not None:
            overwrites["Mapping"].setdefault("ec_superseded", ec_superseded)
        if metacyc_rxn_ec is not None:
            overwrites["Mapping"].setdefault("metacyc_rxn_ec", metacyc_rxn_ec)
        if official_ec_metacyc_rxn is not None:
            overwrites["Mapping"].setdefault("official_ec_metacyc_rxn", official_ec_metacyc_rxn)
        if to_remove_metabolism is not None:
            overwrites["Mapping"].setdefault("to_remove_non_small_molecule_metabolism", to_remove_metabolism)
    return overwrites


