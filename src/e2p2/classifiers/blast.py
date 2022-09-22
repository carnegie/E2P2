import os
import random
import re

from src.definitions import DEFAULT_LOGGER_LEVEL, DEFAULT_LOGGER_NAME
from src.lib.classifier import FunctionClass, Classifier
from src.lib.util import read_delim_itr, logging_helper


class BLAST(Classifier):
    def __init__(self, time_stamp, path_to_weight, name="NCBI BLAST+ blastp", logging_level=DEFAULT_LOGGER_LEVEL,
                 logger_name=DEFAULT_LOGGER_NAME):
        Classifier.__init__(self, time_stamp, path_to_weight, name, logging_level, logger_name)
        self.evalue_threshold = float("1e-2")
        self.bitscore_threshold = float("0")

    def setup_class_w_processed_config(self, input_path, output_path, classifier_config_dict, classifier_name="BLAST",
                                       logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Setting up BLAST", logging_level=logging_level, logger_name=logger_name)
        self.input = input_path
        # input_file_name, input_file_ext = os.path.splitext(os.path.basename(input_path))
        if os.path.isfile(output_path):
            self.output = output_path
        elif os.path.isdir(output_path):
            self.output = os.path.join(output_path, '.'.join(["blast", str(self._time_stamp), "out"]))
        [evalue_threshold, bitscore_threshold, command_string] = \
            self.retrieve_config_dict_options(self._time_stamp, self.input, self.output, classifier_config_dict,
                                              classifier_name, ["evalue_threshold", "bitscore_threshold", "command"],
                                              logging_level=logging_level, logger_name=logger_name)
        try:
            self.evalue_threshold = float(evalue_threshold)
        except (TypeError, ValueError) as e:
            logging_helper("BLAST E-value type error, using default 1e-2.", logging_level="WARNING",
                           logger_name=logger_name)
            self.evalue_threshold = float("1e-2")
        try:
            self.bitscore_threshold = float(bitscore_threshold)
        except (TypeError, ValueError) as e:
            logging_helper("BLAST Bit score type error, using default 0.", logging_level="WARNING",
                           logger_name=logger_name)
            self.bitscore_threshold = float("0")
        try:
            self.command = command_string.split()
        except AttributeError:
            logging_helper("BLAST Command error, none set.", logging_level="WARNING",
                           logger_name=logger_name)
            self.command = None

    def read_classifier_result(self, output_path=None, logging_level=DEFAULT_LOGGER_LEVEL,
                               logger_name=DEFAULT_LOGGER_NAME):
        if output_path is None:
            output_path = self.output
        bit_score_dict = {}
        for query_id, _, _, hit_cls, e_value, bit_score in \
                self.blast_tab_itr(output_path, self.evalue_threshold, self.bitscore_threshold, logger_name):
            try:
                e_value_function_classes = list(self.res[query_id])
                bit_score_function_classes = list(bit_score_dict[query_id])
            except KeyError:
                e_value_function_classes = []
                bit_score_function_classes = []
            for ef_cls in hit_cls:
                try:
                    ef_weight = self.weight[ef_cls]
                except KeyError:
                    ef_weight = 0
                e_value_function_classes.append(FunctionClass(ef_cls, e_value, ef_weight))
                bit_score_function_classes.append(FunctionClass(ef_cls, bit_score, ef_weight))
            if len(e_value_function_classes) == 0 or len(bit_score_function_classes) == 0:
                self.res.setdefault(query_id, [])
                bit_score_dict.setdefault(query_id, [])
                continue
            min_e_value_function_classes = FunctionClass.min_classes(e_value_function_classes)
            min_e_value_function_classes_names = \
                list(set(FunctionClass.get_function_classes_attr_vals(min_e_value_function_classes, "name")))
            bit_score_function_classes_with_min_e_value = \
                FunctionClass.get_function_classes_by_vals(bit_score_function_classes,
                                                           min_e_value_function_classes_names, "name")
            best_bit_score_function_classes = FunctionClass.max_classes(bit_score_function_classes_with_min_e_value)
            best_bit_score_function_classes_names = \
                FunctionClass.get_function_classes_attr_vals(best_bit_score_function_classes, "name")
            best_e_value_function_classes = \
                FunctionClass.get_function_classes_by_vals(min_e_value_function_classes,
                                                           best_bit_score_function_classes_names, "name")
            try:
                self.res[query_id] = best_e_value_function_classes
                bit_score_dict[query_id] = best_bit_score_function_classes
            except KeyError:
                self.res.setdefault(query_id, best_e_value_function_classes)
                bit_score_dict.setdefault(query_id, best_bit_score_function_classes)
        for query in self.res:
            dup_removed = []
            res_of_query = self.res[query]
            function_cls_names = list(set(FunctionClass.get_function_classes_attr_vals(res_of_query, "name")))
            for func_cls_name in function_cls_names:
                all_fun_cls_of_name = FunctionClass.get_function_classes_by_vals(res_of_query, func_cls_name, "name")
                dup_removed.append(random.choice(all_fun_cls_of_name))
            self.res[query] = dup_removed

    @staticmethod
    def blast_tab_itr(path_to_blast_out, evalue_threshold=float("1e-2"), bitscore_threshold=float("0"),
                      logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Reading blast output: \"" + path_to_blast_out + "\"", logging_level="INFO",
                       logger_name=logger_name)
        try:
            with open(path_to_blast_out, 'r') as ptbo:
                for info in read_delim_itr(ptbo, val_indices=[0, 1, 10, 11]):
                    if info:
                        try:
                            e_value = float(info[1][2])
                        except ValueError:
                            e_value = float("1" + info[1][2])
                        try:
                            bit_score = float(info[1][3])
                        except ValueError:
                            bit_score = float("1" + info[1][3])
                        blast_query = [bq.strip() for bq in re.split(r'\s+|\|', info[1][0]) if len(bq.strip()) > 0]
                        blast_hits = [bh.strip() for bh in re.split(r'\s+|\|', info[1][1]) if len(bh.strip()) > 0]
                        try:
                            query_id, query_cls = blast_query[0], blast_query[1:]
                            hit_id, hit_cls = blast_hits[0], blast_hits[1:]
                            if e_value <= float(evalue_threshold) and bit_score >= bitscore_threshold:
                                yield query_id, query_cls, hit_id, hit_cls, e_value, bit_score
                        except IndexError:
                            continue
        except (FileNotFoundError, TypeError) as e:
            raise e


if __name__ == '__main__':
    pass
