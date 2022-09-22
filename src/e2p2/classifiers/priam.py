import os
import re

from src.definitions import DEFAULT_LOGGER_LEVEL, DEFAULT_LOGGER_NAME
from src.lib.classifier import FunctionClass, Classifier
from src.lib.util import file_group_by_start_itr, logging_helper


class PRIAM(Classifier):
    def __init__(self, time_stamp, path_to_weight, name="PRIAM search", logging_level=DEFAULT_LOGGER_LEVEL,
                 logger_name=DEFAULT_LOGGER_NAME):
        Classifier.__init__(self, time_stamp, path_to_weight, name, logging_level, logger_name)
        # evalue currently not used
        self.evalue_threshold = float("1e-2")

    def setup_class_w_processed_config(self, input_path, output_path, classifier_config_dict, classifier_name="PRIAM",
                                       logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Setting up PRIAM", logging_level=logging_level, logger_name=logger_name)
        self.input = input_path
        # input_file_name, input_file_ext = os.path.splitext(os.path.basename(input_path))
        if os.path.isfile(output_path):
            output_folder = os.path.dirname(output_path)
        elif os.path.isdir(output_path):
            output_folder = output_path
        else:
            output_folder = ""
        self.output = os.path.join(output_folder, "PRIAM_%s" % self._time_stamp, "ANNOTATION", "sequenceECs.txt")

        [evalue_threshold, command_string] = \
            self.retrieve_config_dict_options(self._time_stamp, self.input, output_folder, classifier_config_dict,
                                              classifier_name, ["evalue_threshold", "command"], logging_level="INFO",
                                              logger_name="Log")
        try:
            self.evalue_threshold = float(evalue_threshold)
        except (TypeError, ValueError) as e:
            logging_helper("PRIAM E-value type error, using default 1e-2.", logging_level="WARNING",
                           logger_name=logger_name)
            self.evalue_threshold = float("1e-2")
        try:
            self.command = command_string.split()
        except AttributeError:
            logging_helper("PRIAM Command error, none set.", logging_level="WARNING",
                           logger_name=logger_name)
            self.command = None

    def read_classifier_result(self, output_path=None, logging_level=DEFAULT_LOGGER_LEVEL,
                               logger_name=DEFAULT_LOGGER_NAME):
        if output_path is None:
            output_path = self.output
        for query_id, _, ef_cls, _, ef_e_value in self.read_priam_sequence_ec_itr(output_path, logger_name):
            try:
                ef_weight = self.weight[ef_cls]
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
                for info in file_group_by_start_itr(op, start=['>'], skip=['#']):
                    priam_query = [pq.strip() for pq in re.split(r'\s+|\|', info[0]) if len(pq.strip()) > 0]
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
