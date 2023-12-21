import os
import textwrap

from src.definitions import DEFAULT_LOGGER_LEVEL, DEFAULT_LOGGER_NAME, EC_TO_EF_MAPPING_PATH
from src.lib.classifier import Classifier
from src.lib.function_class import FunctionClass
from src.lib.process import PathType, logging_helper
from src.lib.read import read_delim_itr

CONFIG_CLASSIFIER_NAME = "DEEPEC"


class DEEPEC(Classifier):
    def __init__(self, time_stamp, path_to_weight, name=CONFIG_CLASSIFIER_NAME,
                 ec_to_ef_mapping_path=EC_TO_EF_MAPPING_PATH, logging_level=DEFAULT_LOGGER_LEVEL,
                 logger_name=DEFAULT_LOGGER_NAME):
        Classifier.__init__(self, time_stamp, path_to_weight, name, logging_level, logger_name)
        self.ec_to_ef_map = ec_to_ef_mapping_path

    def setup_classifier(self, input_path, output_path, classifier_config_dict, classifier_name=CONFIG_CLASSIFIER_NAME,
                         logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Setting up DeepEC", logging_level=logging_level, logger_name=logger_name)
        self.input = input_path
        output_folder = self.generate_output_paths(self.input, output_path, classifier_name, self._time_stamp)
        self.output = os.path.join(output_folder, "DeepEC_Result.txt")
        [command_string] = \
            self.classifier_config_dict_helper(self._time_stamp, self.input, output_folder, classifier_config_dict,
                                               classifier_name, ["command"],
                                               logging_level="INFO", logger_name="Log")
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
            output_folder = os.path.join(os.path.dirname(output_path), "DeepEC_%s_%s" % (input_file_name, time_stamp))
        elif os.path.isdir(output_path):
            output_folder = output_path
        else:
            input_folder = os.path.dirname(input_path)
            output_folder = os.path.join(input_folder, "DeepEC_%s_%s" % (input_file_name, time_stamp))
        return output_folder

    def read_classifier_result(self, output_path=None, logging_level=DEFAULT_LOGGER_LEVEL,
                               logger_name=DEFAULT_LOGGER_NAME):
        if output_path is None:
            output_path = self.output
        for query_id, ef_cls in self.read_deepec_result_itr(output_path, self.ec_to_ef_map, logger_name):
            try:
                ef_weight = self.weight_map[ef_cls]
            except KeyError:
                ef_weight = 0
            try:
                self.res[query_id].append(FunctionClass(ef_cls, 1, ef_weight))
            except KeyError:
                self.res.setdefault(query_id, [FunctionClass(ef_cls, 1, ef_weight)])

    @staticmethod
    def read_deepec_result_itr(path_to_deepec_result_txt, ec_to_ef_map=None, logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Reading DeepEC DeepEC_Result.txt: \"" + path_to_deepec_result_txt + "\"", logging_level="INFO",
                       logger_name=logger_name)
        ec_to_ef_dict = {}
        if ec_to_ef_map is not None and os.path.isfile(ec_to_ef_map):
            with open(ec_to_ef_map, 'r') as fp:
                for ec_num, ef_cls in read_delim_itr(fp):
                    mapped_efs = set()
                    for ef in ef_cls:
                        mapped_efs.update([cls.strip() for cls in ef.split('|') if cls.strip() != ""])
                    try:
                        ec_to_ef_dict[ec_num].update(mapped_efs)
                    except KeyError:
                        ec_to_ef_dict.setdefault(ec_num, mapped_efs)
        try:
            with open(path_to_deepec_result_txt, 'r') as op:
                for query_id, pred_ecs in read_delim_itr(op, skip="Query ID"):
                    if query_id != "":
                        for ec in sorted(pred_ecs):
                            if ec_to_ef_map is None or len(ec_to_ef_dict) == 0:
                                yield query_id, ec
                            else:
                                try:
                                    ef_cls = sorted(ec_to_ef_dict[ec])
                                    for ef in ef_cls:
                                        yield query_id, ef
                                except KeyError:
                                    continue
        except (FileNotFoundError, TypeError) as e:
            raise e

    @staticmethod
    def add_arguments(argument_parser):
        argument_parser.add_argument("--python_path", "-py", dest="python_path",
                                     help=textwrap.dedent("Command of or path to \"java\"."))
        argument_parser.add_argument("--deepec_path", "-dp", dest="deepec_path", type=PathType('file'),
                                     help=textwrap.dedent("Path to \"deepec.py\"."))
        argument_parser.add_argument("--ec_to_ef_mapping_path", "-ee", dest="ec_to_ef_mapping_path",
                                     type=PathType('file'), help="Path to mapping file from ECs to EFs")

    @staticmethod
    def config_overwrites(args, overwrites=None):
        priam_dest = ["python_path", "deepec_path", "ec_to_ef_mapping_path"]
        if overwrites is None:
            overwrites = {}
        deepec_overwrites = {}
        args_dict = vars(args)
        for dest in priam_dest:
            try:
                val = args_dict[dest]
                if val is not None:
                    key = dest
                    deepec_overwrites.setdefault(key, val)
            except KeyError:
                continue
        if len(deepec_overwrites) > 0:
            overwrites.setdefault(CONFIG_CLASSIFIER_NAME, {})
            overwrites[CONFIG_CLASSIFIER_NAME] = deepec_overwrites
        return overwrites
