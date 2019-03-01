import multiprocessing
import logging
import logging.config
import os
import sys


dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(dir_path, 'e2p2'))

import ensemble
import prog

if __name__ == '__main__':
	q = multiprocessing.Queue()
	cur_logger_config = prog.LoggerConfig()
	cur_logger_config.add_new_logger("file", "mplog.log")
	logging.config.dictConfig(cur_logger_config.dictConfig)

	new_process = prog.RunProcess()
	# new_process.add_process_to_workers(q, "INFO", "file", ["cp", "/Users/bxue/Documents/Carnegie/Projects.docx",
	# 													   "/Users/bxue/Documents/Carnegie/Pro jects.docx_2"], "CP")
	# new_process.add_process_to_workers(q, "INFO", "file", ["ls", ".."])
	# new_process.add_process_to_workers(q, "DEBUG", "file", ["ls", "haha"])
	# new_process.add_process_to_workers(q, "INFO", "file", ["ls"])
	new_process.add_process_to_workers(q, "INFO", "file", ["which", "blastp"])
	new_process.add_process_to_workers(q, "INFO", "file", ["which", "java"])
	new_process.add_process_to_workers(q, "INFO", "file", ["which", "haha"])
	new_process.run_all_worker_processes(q)
	print(new_process.run_results)
