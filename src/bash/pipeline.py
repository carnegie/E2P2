import errno
import os
import textwrap
import time

from src.definitions import DEFAULT_LOGGER_NAME, DEFAULT_OUTPUT_SUFFIX
from src.lib.process import PathType, logging_helper
from src.lib.read import check_fasta_header, remove_splice_variants_from_fasta


def add_io_arguments(argument_parser):
    """Function to add IO related arguments
    Args:
        argument_parser: argparse
    Raises:
    Returns:
    """
    argument_parser.add_argument("--input", "-i", dest="input_file", type=PathType('file'),
                                 help="Path to input protein sequences file", required=True)
    argument_parser.add_argument("--protein_gene", "-pg", dest="protein_gene_path", type=PathType('file'),
                                 help="Provide a protein to gene map. This can be used to generate a "
                                      "splice variant removed fasta file and output the final version of e2p2.")
    argument_parser.add_argument("--remove_splice_variants", "-rm", dest="remove_splice_variants", action="store_true",
                                 help="Argument flag to remove splice variants.")
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


def add_mapping_arguments(argument_parser):
    """Function to add Enzyme Function mapping related arguments
    Args:
        argument_parser: argparse
    Raises:
    Returns:
    """
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


def start_pipeline(input_file, logger_name=DEFAULT_LOGGER_NAME, output_path=None, timestamp=str(time.time()),
                   temp_folder=None, log_path=None, verbose="0"):
    """Function for setting up IO related variables
    Args:
        input_file: input file path
        logger_name: logger name
        output_path: output file path
        timestamp: time stamp
        temp_folder: path to the temp file folder
        log_path: path to the log file
        verbose: verbose level of logging
    Raises:
    Returns:
    """
    check_fasta_header(input_file, logger_name)
    # Setup output paths
    if output_path is None:
        output_path = '.'.join([input_file, DEFAULT_OUTPUT_SUFFIX])
    # Setup temp folder path
    output_folder = os.path.dirname(output_path)
    if temp_folder is None:
        input_file_name, _ = os.path.splitext(os.path.basename(input_file))
        temp_folder = os.path.join(output_folder, input_file_name + '.' + timestamp)

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
                           logger_name=DEFAULT_LOGGER_NAME):
    """Function for mapping protein IDs to gene IDs
    Args:
        input_file: input file path
        output_path: output file path
        protein_gene_path: protein to gene mapping file path
        remove_splice_variants: Boolean value to remove splice variants from input
        logger_name: logger name
    Raises:
    Returns:
    """
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

