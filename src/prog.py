"""
Name:         prog
Version:      180123
Author:       Bo Xue
Description:  The prog module implements useful objects for running processes and logging
"""

import logging
import logging.config
import logging.handlers
import multiprocessing
import subprocess
import threading

logging_levels = {
	"DEBUG": logging.DEBUG,
	"INFO": logging.INFO,
	"WARNING": logging.WARNING,
	"ERROR": logging.ERROR,
	"CRITICAL": logging.CRITICAL
}


class LoggerConfig(object):
	"""Object Class for generating logging config
		Args:
		Raises:
		Returns:
	"""

	def __init__(self):
		self.dictConfig = {
			'version': 1,
			'formatters': {
				'detailed': {
					'class': 'logging.Formatter',
					'format': '%(asctime)s %(name)-15s %(levelname)-8s %(processName)-10s %(message)s'
				}
			},
			'handlers': {
				'console': {
					'class': 'logging.StreamHandler',
					'level': 'INFO',
				},
			},
			'loggers': {
			},
			'root': {
				'level': 'DEBUG',
				'handlers': ['console']
			},
		}

	def add_new_logger(self, logger_name, logger_handler_name, logger_handler_filename, logger_handler_level=None):
		"""Adding a new Logger to dictConfig
			Args:
				logger_name: Name of the new Logger
				logger_handler_name: Name of the handler
				logger_handler_filename: Path to the log
				logger_handler_level: Logger lever, from [DEBUG, INFO, WARNING, ERROR, CRITICAL]
			Raises:
			Returns:
		"""
		logger_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
		new_logger = {
			'handlers': [logger_handler_name]
		}
		new_logger_handler = {
			'class': 'logging.FileHandler',
			'filename': logger_handler_filename,
			'mode': 'w',
			'formatter': 'detailed',
		}
		if logger_handler_level is not None and logger_handler_level in logger_levels:
			new_logger_handler.setdefault('level', logger_handler_level)
		self.dictConfig['handlers'].setdefault(logger_handler_name, new_logger_handler)
		self.dictConfig['loggers'].setdefault(logger_name, new_logger)


class RunProcess(object):
	"""Object Class for running processes
		Args:
		Raises:
		Returns:
	"""

	def __init__(self):
		self.workers = {}
		self.result_queue = multiprocessing.Queue()

	@staticmethod
	def logger_thread(mpq):
		"""Seperate thread for logging
		Args:
			mpq: a multiprocessing Queue
		Raises:
		Returns:
		"""
		while True:
			record = mpq.get()
			if record is None:
				break
			logger = logging.getLogger(record.name)
			logger.handle(record)

	def worker_process(self, mpq, logging_level, logger_name, cmd):
		"""Worker process to run a command, puts process result in a queue
		Args:
			mpq: a multiprocessing Queue
			logging_level: The logging level set for this command
			logger_name: The name of the logger for this command
			cmd: a list for the command
		Raises:
		Returns:
		"""
		qh = logging.handlers.QueueHandler(mpq)
		root = logging.getLogger()
		root.setLevel(logging.DEBUG)
		root.addHandler(qh)
		try:
			call_output = subprocess.check_output(cmd, encoding='UTF-8', stderr=subprocess.STDOUT)
		except subprocess.CalledProcessError as exc:
			self.result_queue.put((' '.join(cmd), logging_level, logger_name, exc.returncode, str(exc.output.strip())))
		else:
			self.result_queue.put((' '.join(cmd), logging_level, logger_name, 0, call_output.strip()))

	def add_process_to_workers(self, mpq, logging_level, logger_name, cmd):
		"""add a worker process to the workers
		Args:
			mpq: a multiprocessing Queue
			logging_level: The logging level set for this command
			logger_name: The name of the logger for this command
			cmd: a list for the command
		Raises:
		Returns:
		"""
		wp = multiprocessing.Process(target=self.worker_process, args=(mpq, logging_level, logger_name, cmd,))
		self.workers.setdefault(' '.join(cmd), (wp, logging_level, logger_name))

	def run_all_worker_processes(self, mpq):
		"""add a worker process to the workers
		Args:
			mpq: a multiprocessing Queue
		Raises:
		Returns:
		"""
		for worker in self.workers:
			self.workers[worker][0].start()
			logger = logging.getLogger(self.workers[worker][2])
			try:
				logger.log(logging_levels[self.workers[worker][1]], "Starting Process \"" + worker + "\"")
			except KeyError:
				logger.log(logging.DEBUG, "Starting Process \"" + worker + "\"")
		lp = threading.Thread(target=self.logger_thread, args=(mpq,))
		lp.start()
		# At this point, the main process could do some useful work of its own
		# Once it's done that, it can wait for the workers to terminate...
		for worker in self.workers:
			self.workers[worker][0].join()
		# And now tell the logging thread to finish up, too
		mpq.put(None)
		lp.join()
		# Retrieve stdout of all queued workers
		for i in range(len(self.workers)):
			cmd, logging_level, logger_name, returncode, output = (self.result_queue.get())
			logger = logging.getLogger(logger_name)
			if returncode != 0:
				for index, line in enumerate(output.split('\n')):
					logger.log(logging.ERROR, "Process Error \"" + cmd + "\", stdout[" + str(index) + "]: " + line)
			else:
				try:
					for index, line in enumerate(output.split('\n')):
						try:
							logger.log(logging_levels[logging_level],
									   "Process Ended \"" + cmd + "\", stdout[" + str(index) + "]: " + line)
						except KeyError:
							logger.log(logging.DEBUG,
									   "Process Ended \"" + cmd + "\", stdout[" + str(index) + "]: " + line)
				except KeyError:
					for index, line in enumerate(output.split('\n')):
						try:
							logger.log(logging_levels[logging_level],
									   "Process Ended \"" + cmd + "\", stdout[" + str(index) + "]: " + line)
						except KeyError:
							logger.log(logging.DEBUG,
									   "Process Ended \"" + cmd + "\", stdout[" + str(index) + "]: " + line)


if __name__ == '__main__':
	q = multiprocessing.Queue()
	cur_logger_config = LoggerConfig()
	cur_logger_config.add_new_logger("file", "file", "mplog.log")
	logging.config.dictConfig(cur_logger_config.dictConfig)

	new_process = RunProcess()
	new_process.add_process_to_workers(q, "INFO", "file", ["cp", "/Users/bxue/Documents/Carnegie/Projects.docx",
														   "/Users/bxue/Documents/Carnegie/Projects.docx_2"])
	new_process.add_process_to_workers(q, "INFO", "file", ["ls", ".."])
	new_process.add_process_to_workers(q, "DEBUG", "file", ["ls", "haha"])
	new_process.add_process_to_workers(q, "INFO", "file", ["ls"])
	new_process.run_all_worker_processes(q)
