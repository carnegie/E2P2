import os
import random
import re
import textwrap

from src.definitions import DEFAULT_LOGGER_LEVEL, DEFAULT_LOGGER_NAME, DEFAULT_BLAST_E_VALUE, DEFAULT_BLAST_BIT_SCORE
from src.lib.classifier import Classifier
from src.lib.function_class import FunctionClass
from src.lib.process import logging_helper, PathType
from src.lib.read import read_delim_itr

CONFIG_CLASSIFIER_NAME = "BLAST"


class BLAST(Classifier):
    def __init__(self, time_stamp, path_to_weight, name=CONFIG_CLASSIFIER_NAME, args=None,
                 logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        Classifier.__init__(self, time_stamp, path_to_weight, name, logging_level, logger_name)
        try:
            if args.blast_e_value is not None:
                self.e_value_threshold = args.blast_e_value
            else:
                self.e_value_threshold = DEFAULT_BLAST_E_VALUE
        except AttributeError:
            self.e_value_threshold = DEFAULT_BLAST_E_VALUE
        try:
            if args.blast_bit_score is not None:
                self.bit_score_threshold = args.blast_bit_score
            else:
                self.bit_score_threshold = DEFAULT_BLAST_BIT_SCORE
        except AttributeError:
            self.bit_score_threshold = DEFAULT_BLAST_BIT_SCORE

    def setup_classifier(self, input_path, output_path, classifier_config_dict, classifier_name=CONFIG_CLASSIFIER_NAME,
                         logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Setting up BLAST", logging_level=logging_level, logger_name=logger_name)
        self.input = input_path
        self.output = self.generate_output_paths(input_path, output_path, classifier_name, self._time_stamp)

        [e_value_threshold, bit_score_threshold, command_string] = \
            self.classifier_config_dict_helper(self._time_stamp, self.input, self.output, classifier_config_dict,
                                               classifier_name, ["blast_e_value", "blast_bit_score", "command"],
                                               logging_level="INFO", logger_name=logger_name)
        try:
            self.e_value_threshold = float(e_value_threshold)
        except (TypeError, ValueError) as e:
            logging_helper("BLAST E-value in config missing or type error, using default " +
                           str(DEFAULT_BLAST_E_VALUE) + ".",
                           logging_level="WARNING", logger_name=logger_name)
            self.e_value_threshold = DEFAULT_BLAST_E_VALUE
        try:
            self.bit_score_threshold = float(bit_score_threshold)
        except (TypeError, ValueError) as e:
            logging_helper("BLAST Bit score in config missing or type error, using default " +
                           str(DEFAULT_BLAST_BIT_SCORE) + ".",
                           logging_level="WARNING", logger_name=logger_name)
            self.bit_score_threshold = DEFAULT_BLAST_BIT_SCORE
        try:
            self.command = command_string.split()
        except AttributeError:
            logging_helper("BLAST Command error, none set.", logging_level="WARNING",
                           logger_name=logger_name)
            self.command = None

    @staticmethod
    def generate_output_paths(input_path, output_path, classifier_name, time_stamp):
        input_file_name, input_file_ext = os.path.splitext(os.path.basename(input_path))
        if os.path.isfile(output_path):
            output = output_path
        elif os.path.isdir(output_path):
            output = os.path.join(output_path, '.'.join(["blast", input_file_name, str(time_stamp), "out"]))
        else:
            input_folder = os.path.dirname(input_path)
            output = os.path.join(input_folder, '.'.join(["blast", input_file_name, str(time_stamp), "out"]))
        return output

    def read_classifier_result(self, output_path=None, logging_level=DEFAULT_LOGGER_LEVEL,
                               logger_name=DEFAULT_LOGGER_NAME):
        if output_path is None:
            output_path = self.output
        bit_score_dict = {}
        for query_id, _, _, hit_cls, e_value, bit_score in \
                self.blast_tab_itr(output_path, self.e_value_threshold, self.bit_score_threshold, logger_name):
            try:
                e_value_function_classes = list(self.res[query_id])
                bit_score_function_classes = list(bit_score_dict[query_id])
            except KeyError:
                e_value_function_classes = []
                bit_score_function_classes = []
            if len(hit_cls) == 0:
                # Workaround for non-enzyme hits
                hit_cls = ['#NA#']
            for ef_cls in hit_cls:
                try:
                    ef_weight = self.weight_map[ef_cls]
                except KeyError:
                    ef_weight = 0
                e_value_function_classes.append(FunctionClass(ef_cls, e_value, ef_weight))
                bit_score_function_classes.append(FunctionClass(ef_cls, bit_score, ef_weight))

            min_e_value = min(e_value_function_classes).score
            min_e_value_function_classes = (
                FunctionClass.get_function_classes_by_vals(e_value_function_classes, min_e_value, 'score'))
            min_e_value_function_classes_names = list(set([fc.name for fc in min_e_value_function_classes]))

            bit_score_function_classes_with_min_e_value = \
                FunctionClass.get_function_classes_by_vals(bit_score_function_classes,
                                                           min_e_value_function_classes_names, "name")
            max_bit_score = max(bit_score_function_classes_with_min_e_value).score
            best_bit_score_function_classes = (
                FunctionClass.get_function_classes_by_vals(bit_score_function_classes_with_min_e_value,
                                                           max_bit_score, 'score'))

            best_function_classes_names = list(set([fc.name for fc in best_bit_score_function_classes]))
            best_function_classes = (
                FunctionClass.get_function_classes_by_vals(min_e_value_function_classes,
                                                           best_function_classes_names, "name"))
            try:
                self.res[query_id] = best_function_classes
                bit_score_dict[query_id] = best_bit_score_function_classes
            except KeyError:
                self.res.setdefault(query_id, best_function_classes)
                bit_score_dict.setdefault(query_id, best_bit_score_function_classes)
        for query in self.res:
            dup_removed = []
            res_of_query = self.res[query]
            function_cls_names = list(set([fc.name for fc in res_of_query]))
            # Workaround for non-enzyme hits
            if '#NA#' in function_cls_names:
                self.res[query] = []
                continue
            for cls_name in function_cls_names:
                func_classes_w_name = FunctionClass.get_function_classes_by_vals(res_of_query, cls_name, "name")
                dup_removed.append(random.choice(func_classes_w_name))
            self.res[query] = dup_removed

    @staticmethod
    def blast_tab_itr(path_to_blast_out, e_value_threshold=DEFAULT_BLAST_E_VALUE,
                      bit_score_threshold=DEFAULT_BLAST_BIT_SCORE, logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Reading blast output: \"" + path_to_blast_out + "\"", logging_level="INFO",
                       logger_name=logger_name)
        try:
            with open(path_to_blast_out, 'r') as ptbo:
                for blast_res in read_delim_itr(ptbo, val_indices=[0, 1, 10, 11]):
                    if blast_res:
                        info = blast_res[1]
                        try:
                            e_value = float(info[2])
                        except ValueError:
                            e_value = float("1" + info[2])
                        try:
                            bit_score = float(info[3])
                        except ValueError:
                            bit_score = float("-1")
                        blast_query = [bq.strip() for bq in re.split(r'[\s|]+', info[0]) if len(bq.strip()) > 0]
                        blast_hits = [bh.strip() for bh in re.split(r'[\s|]+', info[1]) if len(bh.strip()) > 0]
                        try:
                            query_id, query_cls = blast_query[0], blast_query[1:]
                            hit_id, hit_cls = blast_hits[0], blast_hits[1:]
                            if e_value <= float(e_value_threshold) and bit_score >= bit_score_threshold:
                                yield query_id, query_cls, hit_id, hit_cls, e_value, bit_score
                        except IndexError:
                            continue
        except (FileNotFoundError, TypeError) as e:
            raise e

    @staticmethod
    def add_arguments(argument_parser):
        argument_parser.add_argument("--blastp", "-b", dest="blastp",
                                     help=textwrap.dedent("Command of or path to BLAST+ \"blastp\"."))
        argument_parser.add_argument("--num_threads", "-n", dest="num_threads", type=int,
                                     help="Number of threads to run \"blastp\".")
        argument_parser.add_argument("--blast_db", "-bd", dest="blast_db", type=PathType('blast_db'),
                                     help=textwrap.dedent(
                                         "Path to rpsd blast database name.\nFor example, \"/PATH/TO/FOLDER/rpsd.fa\", "
                                         "where you can find the following files in /PATH/TO/FOLDER:rpsd.fa.phr; "
                                         "rpsd.fa.pin; rpsd.fa.psq"))
        argument_parser.add_argument("--blast_e_value", "-be", dest="blast_e_value", type=float,
                                     default=str(DEFAULT_BLAST_E_VALUE), help=textwrap.dedent("Blastp e-value cutoff"))
        argument_parser.add_argument("--blast_weight", "-bw", dest="blast_weight", type=PathType('file'),
                                     help=textwrap.dedent("Path to weight file for the blast classifier"))

    @staticmethod
    def config_overwrites(args, overwrites=None):
        blast_dest = ["blastp", "num_threads", "blast_db", "blast_e_value", "blast_weight"]
        if overwrites is None:
            overwrites = {}
        blast_overwrites = {}
        args_dict = vars(args)
        for dest in blast_dest:
            try:
                val = args_dict[dest]
                if val is not None:
                    key = dest
                    if dest in ["blast_weight"]:
                        key = key.replace("blast_", "")
                    blast_overwrites.setdefault(key, val)
            except KeyError:
                continue
        if len(blast_overwrites) > 0:
            overwrites.setdefault(CONFIG_CLASSIFIER_NAME, {})
            overwrites[CONFIG_CLASSIFIER_NAME] = blast_overwrites
        return overwrites

