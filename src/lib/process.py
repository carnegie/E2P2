import logging
import logging.config
import logging.handlers
import multiprocessing
import os
import subprocess
import sys
import threading
from argparse import ArgumentTypeError
from importlib import util

from src.definitions import DEFAULT_LOGGER_NAME, DEFAULT_LOGGER_LEVEL

logging_levels = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL
}


def logging_helper(log_message, logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
    """Helper function to add log messages to logger
    Args:
        log_message: Message to be logged
        logging_level: The logging level set for this command
        logger_name: The name of the logger for this command
    Raises:
    Returns:
    """
    found_logger = logging.getLogger(logger_name)
    try:
        found_logger.log(logging_levels[logging_level], log_message)
    except KeyError:
        found_logger.log(logging.DEBUG, log_message)


class LoggerConfig(object):
    """Object for generating logging config
    """
    def __init__(self, version=1):
        self.dictConfig = {
            'version': version,
            'formatters': {
                'detailed': {
                    'class': 'logging.Formatter',
                    'format': '%(asctime)s %(name)-15s %(levelname)-8s %(processName)-10s %(message)s'
                }
            },
            'handlers': {
            },
            'loggers': {
            },
            'root': {
                'level': 'DEBUG',
                'handlers': ['console']
            },
        }
        self.dictConfig['handlers'].setdefault('console', {
            'class': 'logging.StreamHandler',
            'level': 'INFO'
        })

    def add_new_logger(self, logger_name, logger_handler_filename, logger_handler_level="INFO",
                       logger_handler_mode='w'):
        """Adding a new Logger to dictConfig
        Args:
            logger_name: Name of the new Logger
            logger_handler_filename: Path to the log
            logger_handler_level: Logger lever, from [DEBUG, INFO, WARNING, ERROR, CRITICAL]
            logger_handler_mode: IO mode for logger
        Raises:
        Returns:
        """
        logger_levels = list(logging_levels.keys())
        new_logger = {
            'handlers': [logger_name]
        }
        if logger_handler_mode not in ['w', 'a', 'w+', 'a+']:
            logger_handler_mode = 'w'
        new_logger_handler = {
            'class': 'logging.FileHandler',
            'filename': logger_handler_filename,
            'mode': logger_handler_mode,
            'formatter': 'detailed',
        }
        if logger_handler_level in logger_levels:
            new_logger_handler.setdefault('level', logger_handler_level)

        self.dictConfig['handlers'].setdefault(logger_name, new_logger_handler)
        self.dictConfig['loggers'].setdefault(logger_name, new_logger)


class RunProcess(object):
    """Object for running processes
    """
    def __init__(self):
        self.mp_queue = multiprocessing.Queue()
        self.run_results = []

    @staticmethod
    def _logger_thread(mpq):
        """Separate thread for logging
        Args:
            mpq: A multiprocessing Queue
        Raises:
        Returns:
        """
        while True:
            record = mpq.get()
            if record is None:
                break
            logger = logging.getLogger(record.name)
            logger.handle(record)

    def _worker_process(self, mpq, logging_level, logger_name, cmd, process_name="Process"):
        """Worker process to run a command, puts process result in a queue
        Args:
            mpq: A multiprocessing Queue
            logging_level: The logging level set for this command
            logger_name: The name of the logger for this command
            cmd: A list for the command
            process_name: Name of the process
        Raises:
        Returns:
        """
        qh = logging.handlers.QueueHandler(mpq)
        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        root.addHandler(qh)
        try:
            command_list = []
            # Clean up command list
            for i in cmd:
                command_list += i.split()
            call_output = subprocess.check_output(command_list, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as exc:
            self.mp_queue.put(
                (' '.join(cmd), logging_level, logger_name, exc.returncode, str(exc.output.strip(), "utf-8"),
                 process_name))
        except FileNotFoundError as e:
            self.mp_queue.put(
                (' '.join(cmd), logging_level, logger_name, e.errno, e.strerror, process_name))
        else:
            self.mp_queue.put(
                (' '.join(cmd), logging_level, logger_name, 0, str(call_output.strip(), "utf-8"), process_name))

    def add_process_to_workers(self, workers, mpq, logging_level, logger_name, cmd, process_name="Process"):
        """Add a worker process to the workers
        Args:
            workers: List of workers
            mpq: A multiprocessing Queue
            logging_level: The logging level set for this command
            logger_name: The name of the logger for this command
            cmd: A string list of the command
            process_name: Name of the process
        Raises:
        Returns:
        """
        wp = multiprocessing.Process(target=self._worker_process,
                                     args=(mpq, logging_level, logger_name, cmd, process_name,))
        workers.append((wp, ' '.join(cmd), logging_level, logger_name, process_name))

    def run_all_worker_processes(self, workers, mpq):
        """add a worker process to the workers
        Args:
            workers: List of workers
            mpq: A multiprocessing Queue
        Raises:
        Returns:
        """
        for worker in workers:
            worker[0].daemon = True
            worker[0].start()
            logging_helper("Starting Process \"" + worker[1] + "\"", logging_level=worker[2],
                           logger_name=worker[3])
        lp = threading.Thread(target=self._logger_thread, args=(mpq,))
        lp.start()
        # Main process waiting for workers terminate
        try:
            for worker in workers:
                worker[0].join()
            # Finish logging
            mpq.put(None)
            lp.join()
        except KeyboardInterrupt:
            for worker in workers:
                worker[0].terminate()
            mpq.put(None)
            lp.join()
            logger = logging.getLogger()
            logger.log(logging.ERROR, "Program interrupted, exiting...")
            sys.exit(1)
        # Retrieve stdout of all queued workers
        for i in range(len(workers)):
            cmd, logging_level, logger_name, return_code, output, process_name = (self.mp_queue.get())
            self.run_results.append((cmd, return_code, output))
            if return_code != 0:
                for index, line in enumerate(output.split('\n')):
                    logging_helper("Process Error \"" + cmd + "\", stdout[" + str(index) + "]: " + line,
                                   logging_level="ERROR", logger_name=logger_name)
            else:
                for index, line in enumerate(output.split('\n')):
                    logging_helper(process_name + " Ended, stdout[" + str(index) + "]: " + line,
                                   logging_level=logging_level, logger_name=logger_name)


def load_module_function_from_path(module_path, function_name, module_name=None):
    """Load in function from a module file path
    Args:
        module_path: File path to module's '.py' file
        function_name: The name of the function to load in
        module_name: A name for the module
    Raises:
    Returns:
        fn: function to call
    """
    if module_name is None:
        module_name = os.path.basename(os.path.dirname(module_path))
    spec = util.spec_from_file_location(module_name, module_path)
    mod = util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)
    fn = mod and getattr(mod, function_name, None)
    return fn


class PathType(object):
    def __init__(self, file_type='file'):
        assert file_type in ('file', 'dir', 'have_parent', 'blast_db', 'blast_bin', 'priam_profiles', None) \
               or hasattr(file_type, '__call__')
        self._type = file_type

    def __call__(self, string):
        path = os.path.realpath(string)
        if self._type is None:
            pass
        elif self._type == 'file':
            if not os.path.isfile(path):
                raise ArgumentTypeError("Path is not a file: '%s'" % string)
        elif self._type == 'dir':
            if not os.path.isdir(path):
                raise ArgumentTypeError("Path is not a directory: '%s'" % string)
        elif self._type == 'have_parent':
            if not os.path.isdir(os.path.dirname(path)):
                raise ArgumentTypeError("Parent of path is not a directory: '%s'" % string)
        elif self._type == 'blast_db':
            phr_path = path + '.phr'
            pin_path = path + '.pin'
            psq_path = path + '.psq'
            if not os.path.isfile(phr_path) or not os.path.isfile(pin_path) or not os.path.isfile(psq_path):
                raise ArgumentTypeError("Cannot find blast database at path: '%s'" % string)
        elif self._type == 'blast_bin':
            bin_files = os.listdir(path)
            if not os.path.isdir(path) or 'makeprofiledb' not in bin_files or 'rpsblast' not in bin_files \
                    or 'rpstblastn' not in bin_files:
                raise ArgumentTypeError("Path not a valid Blast+ bin folder: '%s'" % string)
        elif self._type == 'priam_profiles':
            bin_files = os.listdir(path)
            if not os.path.isdir(path) or 'PROFILES' not in bin_files or 'annotation_rules.xml' not in bin_files \
                    or 'genome_rules.xml' not in bin_files:
                raise ArgumentTypeError("Path not a valid Priam profiles folder: '%s'" % string)
            else:
                profiles_path = os.path.join(path, 'PROFILES')
                annotation_rules = os.path.join(path, 'annotation_rules.xml')
                genome_rules = os.path.join(path, 'genome_rules.xml')
                if not os.path.isdir(profiles_path) or not os.path.isfile(annotation_rules) \
                        or not os.path.isfile(genome_rules):
                    raise ArgumentTypeError("Path not a valid Priam profiles folder: '%s'" % string)
        else:
            raise ArgumentTypeError("path not valid: '%s'" % string)
        return string

