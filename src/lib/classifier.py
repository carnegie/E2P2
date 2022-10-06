import configparser
import multiprocessing
import os.path

from src.definitions import DEFAULT_LOGGER_NAME, DEFAULT_LOGGER_LEVEL
from src.lib.config import get_values_from_config_option
from src.lib.util import read_delim_itr, RunProcess, logging_helper

_available_class_score_attr = ['weight', 'score']


class Error(Exception):
    pass


class FunctionClass(object):
    """Object for function classes
    """
    def __init__(self, name, score, weight=0):
        self.name = name
        self.score = score
        try:
            self.weight = float(weight)
        except TypeError:
            self.weight = float(0)

    def set_weight(self, weight_dict):
        """Set the weight for this function class
        Args:
            weight_dict: Dictionary containing weights for function classes
        Raises: KeyError
        Returns:
        """
        try:
            self.weight = weight_dict[self.name]
        except KeyError:
            self.weight = float(0)

    @staticmethod
    def lesser_than_class(function_class_1, function_class_2, attr='score'):
        """Test if attribute of function class 1 is lesser than function class 2
        Args:
            function_class_1: FunctionClass 1
            function_class_2: FunctionClass 2
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if isinstance(function_class_1, FunctionClass) and isinstance(function_class_2, FunctionClass):
            return getattr(function_class_1, attr) < getattr(function_class_2, attr)
        else:
            raise NotImplementedError

    @staticmethod
    def greater_than_class(function_class_1, function_class_2, attr='score'):
        """Test if attribute of function class 1 is greater than function class 2
        Args:
             function_class_1: FunctionClass 1
            function_class_2: FunctionClass 2
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if isinstance(function_class_1, FunctionClass) and isinstance(function_class_2, FunctionClass):
            return getattr(function_class_1, attr) > getattr(function_class_2, attr)
        else:
            raise NotImplementedError

    @staticmethod
    def lesser_and_equal_class(function_class_1, function_class_2, attr='score'):
        """Test if attribute of function class 1 is lesser than or equal to function class 2
        Args:
            function_class_1: FunctionClass 1
            function_class_2: FunctionClass 2
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if isinstance(function_class_1, FunctionClass) and isinstance(function_class_2, FunctionClass):
            return getattr(function_class_1, attr) <= getattr(function_class_2, attr)
        else:
            raise NotImplementedError

    @staticmethod
    def greater_and_equal_class(function_class_1, function_class_2, attr='score'):
        """Test if attribute of function class 1 is greater than or equal to function class 2
        Args:
            function_class_1: FunctionClass 1
            function_class_2: FunctionClass 2
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if isinstance(function_class_1, FunctionClass) and isinstance(function_class_2, FunctionClass):
            return getattr(function_class_1, attr) > getattr(function_class_2, attr)
        else:
            raise NotImplementedError

    @staticmethod
    def equal_class(function_class_1, function_class_2, attr='score'):
        """Test if attribute of function class 1 is equal to function class 2
        Args:
            function_class_1: FunctionClass 1
            function_class_2: FunctionClass 2
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if isinstance(function_class_1, FunctionClass) and isinstance(function_class_2, FunctionClass):
            return getattr(function_class_1, attr) == getattr(function_class_2, attr)
        else:
            raise NotImplementedError

    @staticmethod
    def lesser_than_threshold(function_class, threshold, attr='score'):
        """Test if attribute of function class is lesser than a threshold value
        Args:
            function_class: FunctionClass
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if isinstance(function_class, FunctionClass):
            return getattr(function_class, attr) < threshold
        else:
            raise NotImplementedError

    @staticmethod
    def greater_than_threshold(function_class, threshold, attr='score'):
        """Test if attribute of function class is greater than a threshold value
        Args:
            function_class: FunctionClass
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if isinstance(function_class, FunctionClass):
            return getattr(function_class, attr) > threshold
        else:
            raise NotImplementedError

    @staticmethod
    def lesser_and_equal_threshold(function_class, threshold, attr='score'):
        """Test if attribute of function class is lesser than or equal to a threshold value
        Args:
            function_class: FunctionClass
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if isinstance(function_class, FunctionClass):
            return getattr(function_class, attr) <= threshold
        else:
            raise NotImplementedError

    @staticmethod
    def greater_and_equal_threshold(function_class, threshold, attr='score'):
        """Test if attribute of function class is greater than or equal to a threshold value
        Args:
            function_class: FunctionClass
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if isinstance(function_class, FunctionClass):
            return getattr(function_class, attr) > threshold
        else:
            raise NotImplementedError

    @staticmethod
    def equal_threshold(function_class, threshold, attr='score'):
        """Test if attribute of function class is equal to a threshold value
        Args:
            function_class: FunctionClass
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if isinstance(function_class, FunctionClass):
            return getattr(function_class, attr) == threshold
        else:
            raise NotImplementedError

    @staticmethod
    def max_attr(list_of_function_classes, attr='score'):
        """Retrieve the maximum attribute of a list of FunctionClasses
        Args:
            list_of_function_classes: List of FunctionClasses
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            The maximum value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if False not in [isinstance(function_class, FunctionClass) for function_class in list_of_function_classes]:
            return max(getattr(function_class, attr) for function_class in list_of_function_classes)
        else:
            raise NotImplementedError

    @staticmethod
    def min_attr(list_of_function_classes, attr='score'):
        """Retrieve the minimum attribute of a list of FunctionClasses
        Args:
            list_of_function_classes: List of FunctionClasses
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            The minimum value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if False not in [isinstance(function_class, FunctionClass) for function_class in list_of_function_classes]:
            return min(getattr(function_class, attr) for function_class in list_of_function_classes)
        else:
            raise NotImplementedError

    @staticmethod
    def max_classes(list_of_function_classes, attr='score'):
        """Retrieve classes with the maximum attribute of a list of FunctionClasses
        Args:
            list_of_function_classes: List of FunctionClasses
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            List of function classes that has the maximum attribute
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if False not in [isinstance(function_class, FunctionClass) for function_class in list_of_function_classes]:
            max_attr = max(getattr(function_class, attr) for function_class in list_of_function_classes)
            return [function_class for function_class in list_of_function_classes
                    if FunctionClass.equal_threshold(function_class, max_attr)]
        else:
            raise NotImplementedError

    @staticmethod
    def min_classes(list_of_function_classes, attr='score'):
        """Retrieve classes with the minimum attribute of a list of FunctionClasses
        Args:
            list_of_function_classes: List of FunctionClasses
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            List of function classes that has the minimum attribute
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        if False not in [isinstance(function_class, FunctionClass) for function_class in list_of_function_classes]:
            min_attr = min(getattr(function_class, attr) for function_class in list_of_function_classes)
            return [function_class for function_class in list_of_function_classes
                    if FunctionClass.equal_threshold(function_class, min_attr)]
        else:
            raise NotImplementedError

    @staticmethod
    def get_function_classes_attr_vals(list_of_function_classes, attr='name'):
        """Retrieve attributes of a list of FunctionClasses
        Args:
            list_of_function_classes: List of FunctionClasses
            attr: Attribute to retrieve
        Raises: NotImplementedError
        Returns:
            List of attribute values of FunctionClasses
        """
        if attr not in _available_class_score_attr:
            attr = 'name'
        if False not in [isinstance(function_class, FunctionClass) for function_class in list_of_function_classes]:
            return [getattr(function_class, attr) for function_class in list_of_function_classes]
        else:
            raise NotImplementedError

    @staticmethod
    def get_function_classes_by_vals(list_of_function_classes, vals, attr='name'):
        """Retrieve classes that has attributes of certain values from a list of FunctionClasses
        Args:
            list_of_function_classes: List of FunctionClasses
            vals: List of values of the attribute
            attr: Attribute to retrieve
        Raises: NotImplementedError
        Returns:
            List of function classes that has the attribute
        """
        if attr not in _available_class_score_attr:
            attr = 'name'
        if False not in [isinstance(function_class, FunctionClass) for function_class in list_of_function_classes]:
            return [function_class for function_class in list_of_function_classes
                    if getattr(function_class, attr) in vals]
        else:
            raise NotImplementedError


class NonCommandError(Error):
    """Error for Classifier missing "_command"
    """
    pass


class Classifier(object):
    def __init__(self, time_stamp, path_to_weight=None, name="", logging_level=DEFAULT_LOGGER_LEVEL,
                 logger_name=DEFAULT_LOGGER_NAME):
        logging_helper("Initializing Classifier: \"" + str(name) + "\"", logging_level=logging_level,
                       logger_name=logger_name)
        self.name = name
        self._time_stamp = time_stamp
        if path_to_weight is not None and os.path.isfile(path_to_weight):
            self.weight = self.read_weights(path_to_weight)
        else:
            self.weight = {}
        # _command: list strings of the bash call
        self.command = None
        # key: Seq ID, val: [FunctionClass, ..]
        self.res = {}
        # IO tracking
        self.input = None
        self.output = None

    def setup_class_w_processed_config(self, input_path, output_path, classifier_config_dict, classifier_name=None,
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
        self.input = input_path
        input_file_name, input_file_ext = os.path.splitext(os.path.basename(input_path))
        if os.path.isfile(output_path):
            self.output = output_path
        elif os.path.isdir(output_path):
            self.output = os.path.join(output_path, input_file_name, '.'.join([classifier_name, str(self._time_stamp)]))

        [command_string] = self.retrieve_config_dict_options(self._time_stamp, self.input, self.output,
                                                             classifier_config_dict, classifier_name, ["command"],
                                                             logging_level=logging_level, logger_name=logger_name)
        try:
            self.command = command_string.split()
        except AttributeError:
            self.command = None

    @staticmethod
    def retrieve_config_dict_options(timestamp, input_path, output_path, classifier_config_dict, classifier_name,
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
                                                  logging_level=logging_level, logger_name=logger_name) for opt in option_list]
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
                            ef_weight = self.weight[ef_class]
                        except KeyError:
                            ef_weight = 0
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

    def add_command_to_queue(self, run_process, worker_list, queue, logging_level=DEFAULT_LOGGER_LEVEL,
                             logger_name=DEFAULT_LOGGER_NAME):
        """Add this classifiers command to a RunProcess worker queue
        Args:
            run_process: RunProcess class
            worker_list: List for workers
            queue: A multiprocessing Queue
            logging_level: The logging level set for this command
            logger_name: The name of the logger for this command
        Raises: ValueError
        Returns:
        """
        if self.command is not None:
            logging_helper("New process: " + self.name + ": \"" + " ".join(self.command) + "\"",
                           logging_level="INFO", logger_name=logger_name)
            run_process.add_process_to_workers(worker_list, queue, logging_level, logger_name, self.command, self.name)

    def get_classifier_result(self):
        return self.res


class RunClassifiers(object):
    """Object for running all the classifiers
    """
    def __init__(self, classifiers=None):
        self.queue = multiprocessing.Queue()
        self.run_process = RunProcess()
        self.workers = []
        if classifiers is None:
            self.classifiers = []
        else:
            self.classifiers = classifiers

    def add_classifier(self, classifier):
        """Add a classifier to this object
        Args:
            classifier: A Classifier
        Raises:
        Returns:
        """
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
        self.classifiers.append(classifier)
        classifier.add_command_to_queue(self.run_process, self.workers, self.queue, logging_level, logger_name)

    def add_all_classifiers_to_queue(self, logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        """Add all classifiers of this object to its multiprocessing queue
        Args:
            logging_level: The logging level set for this command
            logger_name: The name of the logger for this command
        Raises:
        Returns:
        """
        for classifier in self.classifiers:
            classifier.add_command_to_queue(self.run_process, self.workers, self.queue, logging_level, logger_name)

    def run_all_classifiers(self):
        """Run all classifiers of this object
        Args:
        Raises:
        Returns:
        """
        self.run_process.run_all_worker_processes(self.workers, self.queue)

    def get_classifier_results(self):
        """Get list of classifier results
        Args:
        Raises:
        Returns:
            res_list: List of classifier results
        """
        res_list = []
        for classifier in self.classifiers:
            res_list.append(classifier.get_classifier_result())
        return res_list


def run_available_classifiers(list_of_classifier_names, list_of_classifiers, logging_level=DEFAULT_LOGGER_LEVEL,
                              logger_name=DEFAULT_LOGGER_NAME):
    """Placeholder function to read classifier results from output file path.
    Args:
        list_of_classifier_names: List of the classifier names
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
        if isinstance(cls, Classifier) and cls.command is not None:
            run_cls.add_classifier(cls)
        else:
            skipped_classifiers.append(list_of_classifier_names[idx])
    run_cls.add_all_classifiers_to_queue(logging_level, logger_name)
    logging_helper("Running all available processes.", logging_level="INFO", logger_name=logger_name)
    run_cls.run_all_classifiers()
    for cls in run_cls.classifiers:
        cls_output = cls.output
        cls.read_classifier_result(cls_output, logging_level, logger_name)

    return run_cls.classifiers, skipped_classifiers
