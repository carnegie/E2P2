import os
import re
import textwrap

from src.definitions import DEFAULT_LOGGER_LEVEL, DEFAULT_LOGGER_NAME
from src.lib.classifier import FunctionClass, Classifier
from src.lib.process import logging_helper, PathType
from src.lib.read import read_groups_by_start_itr


class PRIAM(Classifier):
    def __init__(self, time_stamp, path_to_weight, name="PRIAM", e_value=float("1e-2"),
                 logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        Classifier.__init__(self, time_stamp, path_to_weight, name, logging_level, logger_name)
        # e value currently not used
        self.e_value_threshold = e_value

    def setup_classifier(self, input_path, output_path, classifier_config_dict, classifier_name="PRIAM",
                         logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Setting up PRIAM", logging_level=logging_level, logger_name=logger_name)
        self.input = input_path
        output_folder = self.generate_output_paths(self.input, output_path, classifier_name, self._time_stamp)
        self.output = os.path.join(output_folder, "PRIAM_%s" % self._time_stamp, "ANNOTATION", "sequenceECs.txt")
        [evalue_threshold, command_string] = \
            self.classifier_config_dict_helper(self._time_stamp, self.input, output_folder, classifier_config_dict,
                                               classifier_name, ["evalue_threshold", "command"],
                                               logging_level="INFO", logger_name="Log")
        try:
            self.e_value_threshold = float(evalue_threshold)
        except (TypeError, ValueError) as e:
            logging_helper("PRIAM E-value type error, using default 1e-2.", logging_level="WARNING",
                           logger_name=logger_name)
            self.e_value_threshold = float("1e-2")
        try:
            self.command = command_string.split()
        except AttributeError:
            logging_helper("PRIAM Command error, none set.", logging_level="WARNING",
                           logger_name=logger_name)
            self.command = None

    @staticmethod
    def generate_output_paths(input_path, output_path, classifier_name, time_stamp):
        input_file_name, input_file_ext = os.path.splitext(os.path.basename(input_path))
        if os.path.isfile(output_path):
            output_folder = os.path.dirname(output_path)
        elif os.path.isdir(output_path):
            output_folder = os.path.join(output_path, "PRIAM_%s_%s" % (input_file_name, time_stamp))
        else:
            input_folder = os.path.dirname(input_path)
            output_folder = os.path.join(input_folder, "PRIAM_%s_%s" % (input_file_name, time_stamp))
        # output = os.path.join(output_folder, "PRIAM_%s" % time_stamp, "ANNOTATION", "sequenceECs.txt")
        return output_folder

    def read_classifier_result(self, output_path=None, logging_level=DEFAULT_LOGGER_LEVEL,
                               logger_name=DEFAULT_LOGGER_NAME):
        if output_path is None:
            output_path = self.output
        for query_id, _, ef_cls, _, ef_e_value in self.read_priam_sequence_ec_itr(output_path, logger_name):
            try:
                ef_weight = self.weight_map[ef_cls]
            except KeyError:
                ef_weight = 0
            try:
                self.res[query_id].append(FunctionClass(ef_cls, ef_e_value, ef_weight))
            except KeyError:
                self.res.setdefault(query_id, [FunctionClass(ef_cls, ef_e_value, ef_weight)])

    @staticmethod
    def read_priam_sequence_ec_itr(path_to_sequence_ec_txt, logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Reading PRIAM sequenceEC.txt: \"" + path_to_sequence_ec_txt + "\"", logging_level="INFO",
                       logger_name=logger_name)
        try:
            with open(path_to_sequence_ec_txt, 'r') as op:
                for info in read_groups_by_start_itr(op, start=['>'], skip=['#']):
                    priam_query = [pq.strip() for pq in re.split(r'[\s|]+', info[0]) if len(pq.strip()) > 0]
                    priam_results = [pr.split('\t') for pr in info[1]]
                    for ef_class_res in priam_results:
                        try:
                            query_id, query_cls = priam_query[0], priam_query[1:]
                            ef_class = ef_class_res[0].strip()
                            ef_prob = ef_class_res[1].strip()
                            ef_e_value = ef_class_res[2].strip()
                            yield query_id.lstrip('>'), query_cls, ef_class, float(ef_prob), float(ef_e_value)
                        except (IndexError, ValueError) as e:
                            continue
        except (FileNotFoundError, TypeError) as e:
            raise e

    @staticmethod
    def add_arguments(argument_parser):
        argument_parser.add_argument("--java_path", "-j", dest="java_path",
                                     help=textwrap.dedent("Command of or path to \"java\"."))
        argument_parser.add_argument("--priam_search", "-ps", dest="priam_search", type=PathType('file'),
                                     help=textwrap.dedent("Path to \"PRIAM_search.jar\"."))
        argument_parser.add_argument("--priam_resume", "-pr", dest="priam_resume", action='store_true',
                                     help="Whether or not to resume a found PRIAM_search.jar process.")
        argument_parser.add_argument("--blast_bin", "-bb", dest="blast_bin", type=PathType('blast_bin'),
                                     help=textwrap.dedent("Command of or path to BLAST+ bin folder."))
        argument_parser.add_argument("--priam_profiles", "-pp", dest="priam_profiles", type=PathType('priam_profiles'),
                                     help=textwrap.dedent(
                                         "Path to PRIAM profiles.\nFor example, \"/PATH/TO/FOLDER/profiles\", "
                                         "where you can find the following in /PATH/TO/FOLDER/profiles:\n "
                                         "files: annotation_rules.xml; genome_rules.xml\n "
                                         "folders: PROFILES: Folder contains \"LIBRARY\" folder and "
                                         "multiple \".chk\" files."))
        argument_parser.add_argument("--priam_weight", "-pw", dest="priam_weight", type=PathType('file'),
                                     help=textwrap.dedent("Path to weight file for the priam classifier"))

    @staticmethod
    def config_overwrites(args, overwrites=None):
        priam_dest = ["java_path", "priam_search", "priam_resume", "blast_bin", "priam_profiles", "priam_weight"]
        if overwrites is None:
            overwrites = {}
        priam_overwrites = {}
        args_dict = vars(args)
        for dest in priam_dest:
            try:
                val = args_dict[dest]
                if val is not None:
                    key = dest
                    if dest in ["priam_resume", "priam_weight"]:
                        key = key.replace("priam_", "")
                        if dest == 'priam_resume':
                            if val is True:
                                val = "fr"
                            else:
                                val = "fn"
                    priam_overwrites.setdefault(key, val)
            except KeyError:
                continue
        if len(priam_overwrites) > 0:
            overwrites.setdefault("PRIAM", {})
            overwrites["PRIAM"] = priam_overwrites
        return overwrites

