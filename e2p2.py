import argparse
import logging.config
import os
import re
import sys

from src.definitions import DEFAULT_CONFIG_PATH, ROOT_DIR
from src.bash.pipeline import *
from src.lib.classifier import run_available_classifiers
from src.lib.config import read_config
from src.lib.ensemble import run_all_ensembles
from src.lib.process import LoggerConfig, logging_helper, load_module_function_from_path
from src.lib.read import get_all_seq_ids_from_fasta
from src.lib.write import PfFiles, write_ensemble_outputs


project_path = os.path.dirname(__file__)
sys.path.insert(0, project_path)


def main():
    name = 'e2p2.py'
    description = '''
    Runs the Ensemble Enzyme Prediction Pipeline (E2P2) on a set of input protein sequences,
    outputting enzyme functional annotations in the forms of EC numbers or MetaCyc reaction
    IDs for any predicted enzyme sequences.
    '''
    notes = '''
    - Input protein sequences should be in FASTA format.
    - Headers in the FASTA file should begin with the sequence ID followed by a space or "|".
    - Intermediate results files can be found in a temporary directory of its own subdirectory labeled with a date and time stamp.
    '''
    time_stamp = str(int(time.time()))
    cur_logger_config = LoggerConfig()
    parser = argparse.ArgumentParser(prog=name, description=description, formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=textwrap.dedent(notes))
    add_io_arguments(parser)
    subparsers = parser.add_subparsers()
    parser_e2p2 = subparsers.add_parser('e2p2', help=textwrap.dedent("Argument to run E2P2."))

    # Config read in
    args, others = parser.parse_known_args()
    config_log_flag = False
    if args.config_ini is not None and os.path.isfile(args.config_ini):
        config_path = args.config_ini
    else:
        config_log_flag = True
        config_path = DEFAULT_CONFIG_PATH
    mapping_files, classifier_dict, ensemble_dict = read_config(config_path)
    if None in (mapping_files, classifier_dict, ensemble_dict):
        parser.print_help()
        raise SystemExit

    for cls in classifier_dict:
        cls_path = os.path.join(ROOT_DIR, classifier_dict[cls]["class"])
        cls_fn = load_module_function_from_path(cls_path, cls)
        cls_fn.add_arguments(parser_e2p2)

    for ens in ensemble_dict:
        ens_path = os.path.join(ROOT_DIR, ensemble_dict[ens]["class"])
        ens_fn = load_module_function_from_path(ens_path, ens)
        ens_fn.add_arguments(parser_e2p2)

    # Parse arguments
    args = parser.parse_args()
    output_path, io_dict, create_temp_folder_flag, log_path, logging_level = \
        start_pipeline(args.input_file, output_path=args.output_path, temp_folder=args.temp_folder,
                       log_path=args.log_path, verbose=args.verbose, timestamp=time_stamp)


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
    if config_log_flag is True:
        logging_helper("No user provided config.ini is found, attempting to use file at %s." % DEFAULT_CONFIG_PATH,
                       logging_level="INFO", logger_name=DEFAULT_LOGGER_NAME)
    fasta_path = \
        protein_to_gene_helper(args.input_file, output_path, args.protein_gene_path, args.remove_splice_variants,
                               logger_name=DEFAULT_LOGGER_NAME)
    io_dict["IO"]["query"] = fasta_path

    all_query_ids = get_all_seq_ids_from_fasta(fasta_path)

    # Overwrite config with arguments
    overwrites = {}
    for cls in classifier_dict:
        cls_path = os.path.join(ROOT_DIR, classifier_dict[cls]["class"])
        cls_fn = load_module_function_from_path(cls_path, cls)
        io_dict["IO"][cls] = cls_fn.generate_output_paths(io_dict["IO"]["query"], io_dict["IO"]["out"], cls, time_stamp)
        cls_fn.config_overwrites(args, overwrites)
    for ens in ensemble_dict:
        ens_path = os.path.join(ROOT_DIR, ensemble_dict[ens]["class"])
        ens_fn = load_module_function_from_path(ens_path, ens)
        ens_fn.config_overwrites(args, overwrites)
    _, classifier_dict, ensemble_dict = read_config(config_path, io_dict, overwrites)

    # Set up classifiers
    classifier_names = sorted(classifier_dict.keys())
    list_of_classifiers = []
    for cls in classifier_names:
        cls_path = os.path.join(ROOT_DIR, classifier_dict[cls]["class"])
        path_to_weight = classifier_dict[cls]["weight"]
        cls_fn = load_module_function_from_path(cls_path, cls)
        cls_classifier = cls_fn(time_stamp=time_stamp, path_to_weight=path_to_weight, args=args)
        cls_classifier.setup_classifier(io_dict["IO"]["query"], io_dict["IO"]["out"], classifier_dict[cls])
        list_of_classifiers.append(cls_classifier)

    # Run Classifiers
    res_cls_list, skipped_classifiers = \
        run_available_classifiers(classifier_names, list_of_classifiers, logging_level, DEFAULT_LOGGER_NAME)

    # Set up ensembles
    ensemble_names = sorted(ensemble_dict.keys())
    list_of_ensembles = []
    for ens in ensemble_names:
        ens_path = os.path.join(ROOT_DIR, ensemble_dict[ens]["class"])
        threshold = ensemble_dict[ens]["threshold"]
        ens_fn = load_module_function_from_path(ens_path, ens)
        ens_ensemble = ens_fn(res_cls_list, time_stamp, ens, threshold)
        list_of_ensembles.append(ens_ensemble)

    # Run Ensembles
    ensembles_ran, skipped_ensembles = \
        run_all_ensembles(ensemble_names, list_of_ensembles, all_query_ids, DEFAULT_LOGGER_NAME)

    for ensemble_cls in ensembles_ran:
        # ensemble_name = ensemble_cls.name
        # ensemble_classifiers = ensemble_cls.list_of_classifiers
        write_ensemble_outputs(ensemble_cls, all_query_ids, output_path,
                               os.path.join(ROOT_DIR, mapping_files['efclasses']),
                               os.path.join(ROOT_DIR, mapping_files['ec_superseded']),
                               os.path.join(ROOT_DIR, mapping_files['metacyc_rxn_ec']),
                               os.path.join(ROOT_DIR, mapping_files['official_ec_metacyc_rxn']),
                               os.path.join(ROOT_DIR, mapping_files['to_remove_non_small_molecule_metabolism']),
                               prot_gene_map_path=args.protein_gene_path, logging_level=logging_level,
                               logger_name=DEFAULT_LOGGER_NAME)


if __name__ == '__main__':
    main()



