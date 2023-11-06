import argparse
import logging.config
import os
import re
import sys

from src.bash.pipeline import *
from src.lib.classifier import run_available_classifiers
from src.lib.config import read_config_ini
from src.lib.ensemble import run_all_ensembles
from src.lib.util import LoggerConfig, logging_helper, get_all_seq_ids_from_fasta
from src.lib.write import PfFiles, write_ensemble_outputs

from src.e2p2.classifiers.blast import blast_overwrites
from src.e2p2.classifiers.priam import priam_overwrites


project_path = os.path.dirname(__file__)
sys.path.insert(0, project_path)


if __name__ == '__main__':
    name = 'e2p2.py'
    description = '''
    Runs the Ensemble Enzyme Prediction Pipeline (E2P2) on a set of input protein sequences,
    outputting enzyme functional annotations in the forms of EC numbers or MetaCyc reaction
    IDs for any predicted enzyme sequences.
    '''
    notes = '''
    - Input protein sequences should be in FASTA format.
    - Headers in the FASTA file should begin with the sequence ID followed by a space or "|".
    - Intermediate results files can be found in a temporary directory of its own subdirectory labeled with a date and 
    time stamp.
    '''
    time_stamp = str(time.time())

    parser = argparse.ArgumentParser(prog=name, description=description, formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=textwrap.dedent(notes))
    subparsers = parser.add_subparsers()
    parser_e2p2 = subparsers.add_parser('e2p2', help=textwrap.dedent("Argument to run E2P2."))
    # parser_classify = subparsers.add_parser('classify', help=textwrap.dedent("Argument to run classifications only."))
    # parser_ensemble = subparsers.add_parser('ensemble', help=textwrap.dedent("Argument to run ensembles only."))

    # Argument for IO
    add_io_arguments(parser_e2p2)
    add_classifier_arugments(parser_e2p2)
    add_ensemble_arguments(parser_e2p2)
    add_mapping_arguments(parser_e2p2)
    # Select classifiers to be used
    args = parser.parse_args()

    output_path, io_dict, create_temp_folder_flag, log_path, logging_level = \
        io_helper(args.input_file, output_path=args.output_path, temp_folder=args.temp_folder,
                  log_path=args.log_path, verbose=args.verbose, timestamp=time_stamp)

    cur_logger_config = LoggerConfig()
    if os.path.isfile(os.path.realpath(log_path)):
        cur_logger_config.add_new_logger(DEFAULT_LOGGER_NAME, log_path, logger_handler_mode='a')
    else:
        cur_logger_config.add_new_logger(DEFAULT_LOGGER_NAME, log_path)
    logging.config.dictConfig(cur_logger_config.dictConfig)
    logger = logging.getLogger(DEFAULT_LOGGER_NAME)
    if create_temp_folder_flag:
        logging_helper("Temp folder created at path %s." % io_dict["IO"]["out"], logging_level=logging_level,
                       logger_name=DEFAULT_LOGGER_NAME)
    else:
        logging_helper("Using path %s as temp folder." % io_dict["IO"]["out"], logging_level=logging_level,
                       logger_name=DEFAULT_LOGGER_NAME)
    if os.path.isfile(os.path.realpath(log_path)):
        logger.log(logging.WARNING, "Log file %s exists, will append to it..." % log_path)

    fasta_path = \
        protein_to_gene_helper(args.input_file, output_path, args.protein_gene_path, args.remove_splice_variants,
                               logger_name=DEFAULT_LOGGER_NAME)
    io_dict["IO"]["query"] = fasta_path
    all_query_ids = get_all_seq_ids_from_fasta(fasta_path)
    overwrites = \
        blast_overwrites(args.blastp_cmd, args.blast_db, args.num_threads, args.blast_evalue, args.blast_weight,
                         args.blast_cls, overwrites=None)
    overwrites = \
        priam_overwrites(args.java_cmd, args.priam_search, args.priam_resume, args.blast_bin, args.priam_profiles,
                         args.priam_weight, args.priam_cls, overwrites=overwrites)
    overwrites = ensemble_overwrites(args.ensemble_cls, args.threshold, overwrites=overwrites)
    overwrites = \
        mapping_overwrite(args.ef_map, args.ec_superseded, args.metacyc_rxn_ec, args.official_ec_metacyc_rxn,
                          args.to_remove_metabolism, overwrites=overwrites)

    classifier_sections, list_of_classifiers, ensemble_sections, list_of_ensembles, mapping_files = \
        read_config_ini(time_stamp, args.config_file, io_dict, overwrites=overwrites, logging_level="WARNING",
                        logger_name=DEFAULT_LOGGER_NAME)

    res_cls_list, skipped_classifiers = \
        run_available_classifiers(classifier_sections, list_of_classifiers, logging_level=logging_level)

    ensembles_ran, skipped_ensembles = \
        run_all_ensembles(ensemble_sections, list_of_ensembles, res_cls_list, time_stamp, logging_level,
                          DEFAULT_LOGGER_NAME)

    for ensemble_cls in ensembles_ran:
        ensemble_name = re.sub(r'[^\w\-_\. ]', '_', ensemble_cls.prediction.name)
        ensemble_classifiers = ensemble_cls.list_of_classifiers
        ensemble_output = PfFiles(ensemble_cls.get_prediction(all_query_ids))
        write_ensemble_outputs(ensemble_cls, all_query_ids, output_path, mapping_files['efclasses'],
                               mapping_files['ec_superseded'], mapping_files['metacyc_rxn_ec'],
                               mapping_files['official_ec_metacyc_rxn'],
                               mapping_files['to_remove_non_small_molecule_metabolism'],
                               prot_gene_map_path=args.protein_gene_path, logging_level=DEFAULT_LOGGER_LEVEL,
                               logger_name=DEFAULT_LOGGER_NAME)




