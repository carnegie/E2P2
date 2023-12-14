import configparser
import multiprocessing
import os.path

from src.definitions import DEFAULT_LOGGER_NAME, DEFAULT_LOGGER_LEVEL
from src.lib.config import get_values_from_config_option
from src.lib.function_class import FunctionClass
from src.lib.process import RunProcess, logging_helper
from src.lib.read import read_delim_itr

_available_class_score_attr = ['weight', 'score']


class Error(Exception):
    pass


class NonCommandError(Error):
    """Error for Classifier missing "_command"
    """
    pass


class Classifier(object):
    def __init__(self, time_stamp, path_to_weight=None, name="", input_path="", output_path="",
                 logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Initializing Classifier: \"" + str(name) + "\"", logging_level=logging_level,
                       logger_name=logger_name)
        self.name = name
        self._time_stamp = time_stamp
        if path_to_weight is not None and os.path.isfile(path_to_weight):
            self.weight_map = self.read_weights(path_to_weight)
        else:
            self.weight_map = {}
        # _command: list strings of the bash call
        self.command = None
        # key: Seq ID, val: [FunctionClass, ..]
        self.res = {}
        # IO tracking
        self.input = input_path
        self.output = output_path

    def __repr__(self):
        return f'Classifier(\'{self.name}\', {self.command}, {self._time_stamp})'

    def setup_classifier(self, input_path, output_path, classifier_config_dict, classifier_name=None,
                         logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        """Placeholder function to set up a classifier using a processed config dict read from a config file
        Args:
            input_path: Input file path for the classifier
            output_path: Output file path for the classifier
            classifier_config_dict: Dictionary of the classifier read from the config
            classifier_name: Name of the classifier in the config file
            logging_level: The logging level set for this command
            logger_name: The name of the logger for this command
        Raises: AttributeError
        Returns:
        """
        if classifier_name is None:
            classifier_name = self.name
        logging_helper("Setting up " + classifier_name, logging_level=logging_level, logger_name=logger_name)
        if os.path.isfile(input_path):
            self.input = input_path
        self.output = self.generate_output_paths(input_path, output_path, classifier_name, self._time_stamp)
        [command_string] = self.classifier_config_dict_helper(self._time_stamp, self.input, self.output,
                                                              classifier_config_dict, classifier_name, ["command"],
                                                              logging_level=logging_level, logger_name=logger_name)
        try:
            self.command = command_string.split()
        except AttributeError:
            self.command = None

    @staticmethod
    def generate_output_paths(input_path, output_path, classifier_name, time_stamp):
        input_file_name, input_file_ext = os.path.splitext(os.path.basename(input_path))
        if os.path.isfile(output_path):
            output = output_path
        elif os.path.isdir(output_path):
            output = os.path.join(output_path, input_file_name, '.'.join([classifier_name, str(time_stamp)]))
        else:
            input_folder = os.path.dirname(input_path)
            output = os.path.join(input_folder, input_file_name, '.'.join([classifier_name, str(time_stamp)]))
        return output

    @staticmethod
    def classifier_config_dict_helper(timestamp, input_path, output_path, classifier_config_dict, classifier_name,
                                      option_list, logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        """Retrieve option values from a processed config dict read from a config file, include interpolation.
        Args:
            timestamp: The time stamp of the function
            input_path: Input file path for the classifier
            output_path: Output file path for the classifier
            classifier_config_dict: Dictionary of the classifier read from the config
            classifier_name: Name of the classifier in the config file
            option_list: list of config options to retrieve
            logging_level: The logging level set for this command
            logger_name: The name of the logger for this command
        Raises: AttributeError
        Returns:
        """
        config_dict = {
            "IO": {
                "query": input_path,
                classifier_name: output_path,
                "timestamp": timestamp
            },
            classifier_name: classifier_config_dict
        }
        classifier_config = configparser.ConfigParser(allow_no_value=True,
                                                      interpolation=configparser.ExtendedInterpolation())
        try:
            classifier_config.read_dict(config_dict)
            return [get_values_from_config_option(classifier_config, classifier_name, opt,
                                                  logging_level=logging_level, logger_name=logger_name)
                    for opt in option_list]
        except AttributeError:
            return None

    def read_classifier_result(self, output_path=None, logging_level=DEFAULT_LOGGER_LEVEL,
                               logger_name=DEFAULT_LOGGER_NAME):
        """Placeholder function to read classifier results from output file path.
        Args:
            output_path: Output file path for the classifier
            logging_level: The logging level set for this command
            logger_name: The name of the logger for this command
        Raises: FileNotFoundError, TypeError
        Returns:
        """
        if output_path is None:
            output_path = self.output
        try:
            with open(output_path, 'r') as op:
                for info in read_delim_itr(op, key_idx=0, val_indices=[1, 2]):
                    if info:
                        seq_id = info[0]
                        ef_class = info[1]
                        ef_score = float(info[2])
                        try:
                            ef_weight = self.weight_map[ef_class]
                        except KeyError:
                            ef_weight = float(0)
                        try:
                            self.res[seq_id].append(FunctionClass(ef_class, ef_score, ef_weight))
                        except KeyError:
                            self.res.setdefault(seq_id, [FunctionClass(ef_class, ef_score, ef_weight)])
        except (FileNotFoundError, TypeError) as e:
            raise e

    @staticmethod
    def read_weights(path_to_weight, logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        """Read in weights from file
        Args:
            path_to_weight: File path to classifier's weight
            logging_level: The logging level set for this command
            logger_name: The name of the logger for this command
        Raises: ValueError
        Returns:
        """
        logging_helper("Reading weight file: \"" + str(path_to_weight) + "\"", logging_level=logging_level,
                       logger_name=logger_name)
        weight = {}
        with open(path_to_weight, 'r') as fp:
            for function_cls, cls_weight in read_delim_itr(fp):
                try:
                    weight.setdefault(function_cls, float(cls_weight[0]))
                except ValueError:
                    weight.setdefault(function_cls, float("0.0"))
        return weight

    def get_res(self):
        return self.res


class RunClassifiers(object):
    """Object for running all the classifiers
    """
    def __init__(self, classifiers=None):
        self.queue = multiprocessing.Queue()
        self.run_process = RunProcess()
        self.workers = []
        if classifiers is None or type(classifiers) not in [list, set]:
            self.classifiers = []
        else:
            valid_classifiers = [classifier for classifier in classifiers if isinstance(classifier, Classifier)]
            self.classifiers = valid_classifiers

    def add_classifier(self, classifier):
        """Add a classifier to this object
        Args:
            classifier: A Classifier
        Raises:
        Returns:
        """
        if isinstance(classifier, Classifier):
            self.classifiers.append(classifier)

    def add_classifier_to_queue(self, classifier, logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        """Add a classifier to this object, and it's multiprocessing queue
        Args:
            classifier: A Classifier
            logging_level: The logging level set for this command
            logger_name: The name of the logger for this command
        Raises:
        Returns:
        """
        if isinstance(classifier, Classifier) and classifier.command is not None:
            self.classifiers.append(classifier)
            logging_helper("New process: " + classifier.name + ": \"" + " ".join(classifier.command) + "\"",
                           logging_level="INFO", logger_name=logger_name)
            self.run_process.add_process_to_workers(self.workers, self.queue, logging_level, logger_name,
                                                    classifier.command, classifier.name)

    def add_available_classifiers_to_queue(self, logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        """Add all classifiers of this object to its multiprocessing queue
        Args:
            logging_level: The logging level set for this command
            logger_name: The name of the logger for this command
        Raises:
        Returns:
        """
        for classifier in self.classifiers:
            logging_helper("New process: " + classifier.name + ": \"" + " ".join(classifier.command) + "\"",
                           logging_level="INFO", logger_name=logger_name)
            self.run_process.add_process_to_workers(self.workers, self.queue, logging_level, logger_name,
                                                    classifier.command, classifier.name)

    def run(self, logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        """Run all classifiers of this object
        Args:
        Raises:
        Returns:
        """
        logging_helper("Running all available processes.", logging_level=logging_level, logger_name=logger_name)
        self.run_process.run_all_worker_processes(self.workers, self.queue)

    def res(self):
        """Get list of classifier results
        Args:
        Raises:
        Returns:
            res_list: List of classifier results
        """
        res_list = []
        for classifier in self.classifiers:
            res_list.append(classifier.get_res())
        return res_list

    @staticmethod
    def add_arguments(argument_parser):
        argument_parser.add_argument('Classifier')

    @staticmethod
    def config_overwrites(args, overwrites=None):
        pass


def run_available_classifiers(classifiers_to_run, list_of_classifiers, logging_level=DEFAULT_LOGGER_LEVEL,
                              logger_name=DEFAULT_LOGGER_NAME):
    """Placeholder function to read classifier results from output file path.
    Args:
        classifiers_to_run: List of the classifier names that will be run
        list_of_classifiers: List of the classifier classes
        logging_level: The logging level set for this command
        logger_name: The name of the logger for this command
    Raises:
    Returns:
        list of classifiers that were run, list of classifiers that were skipped
    """
    run_cls = RunClassifiers()
    skipped_classifiers = []
    for idx, cls in enumerate(list_of_classifiers):
        if isinstance(cls, Classifier) and cls.command is not None and cls.name in classifiers_to_run:
            run_cls.add_classifier(cls)
        else:
            skipped_classifiers.append(classifiers_to_run[idx])
    run_cls.add_available_classifiers_to_queue(logging_level, logger_name)
    run_cls.run(logging_level, logger_name)
    for cls in run_cls.classifiers:
        cls_output = cls.output
        cls.read_classifier_result(cls_output, logging_level, logger_name)

    return run_cls.classifiers, skipped_classifiers

