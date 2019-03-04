"""
Name:			ensemble
Version:		180123
Author:			Bo Xue, Chuan Wang, Lee Chae
Description:	The ensemble module implements ensemble integration algorithms needed by
				E2P2 (Ensemble Enzyme Prediction Pipeline) to provide a final classification
				for a given query protein. Updated for Python 3

"""
import logging
import logging.config
import multiprocessing
import os
import re
import sys

# pipeline modules
import prog


class RunClassifiers(object):
	"""Object for running all classifiers
	"""

	def __init__(self, time_stamp):
		self.time_stamp = time_stamp
		self.queue = multiprocessing.Queue()
		self.classifier_processes = prog.RunProcess()
		self.output_dict = {}

	def blast_classifer(self, blastp_path, path_to_blast_db_basename, path_to_input, path_to_output, num_threads,
						logging_level, logger_name):
		"""Add subprocess of blastp
		Args:
			blastp_path: path/command for blastp
			path_to_blast_db_basename: Name for the blast database
			path_to_input: Path to the fasta input
			path_to_output: Path to the blastp output
			logging_level: The logging level set for reading weight file
			logger_name: The name of the logger for reading weight file
		Raises:
		Returns:
		"""
		logger = logging.getLogger(logger_name)
		self.output_dict.setdefault("blast", path_to_output)
		try:
			num_threads = str(int(num_threads))
		except ValueError:
			logger.log(logging.DEBUG, "num_threads input error, using \"1\" instead.")
			num_threads = "1"
		blast_cmd = [blastp_path, '-db', path_to_blast_db_basename, '-num_threads', num_threads,
					 '-query', path_to_input, '-out', path_to_output, '-outfmt', '6']
		try:
			logger.log(prog.logging_levels[logging_level], "Setup process for blastp: \"" + " ".join(blast_cmd) + "\"")
			self.classifier_processes.add_process_to_workers(self.queue, logging_level, logger_name, blast_cmd,
															 "blastp")
		except KeyError:
			logger.log(logging.DEBUG, "Setup process for blastp: \"" + " ".join(blast_cmd) + "\"")
			self.classifier_processes.add_process_to_workers(self.queue, logging_level, logger_name, blast_cmd,
															 "blastp")

	def priam_classifer(self, java_path, priam_search_path, path_to_blast_bin, path_to_priam_profiles, path_to_input,
						path_to_output_folder, logging_level, logger_name, resume=True):
		"""Add subprocess of PRIAM_search.jar
		Args:
			java_path: path/command for blastp
			priam_search_path: Name for the blast database
			path_to_blast_bin: Path to blast+ bin folder
			path_to_priam_profiles: Path to PRIAM profiles folder
			path_to_input: Path to the fasta input
			path_to_output_folder: Path to folder PRIAM_search.jar would generate output
			logging_level: The logging level set for running PRIAM_search.jar
			logger_name: The name of the logger for running PRIAM_search.jar
			resume: Whether or not to resume a preexisting job.
		Raises:
		Returns:
		"""
		logger = logging.getLogger(logger_name)
		self.output_dict.setdefault("priam",
									os.path.join(path_to_output_folder, "PRIAM_%s" % (self.time_stamp), "ANNOTATION",
												 "sequenceECs.txt"))
		if resume:
			priam_cmd = [java_path, '-Xms3072m', '-Xmx3072m', '-jar', priam_search_path, '--bd', path_to_blast_bin,
						 '--bp', '-n', self.time_stamp, '-i', path_to_input, '-p', path_to_priam_profiles, '--bh', '-o',
						 path_to_output_folder, '--fr']
		else:
			priam_cmd = [java_path, '-Xms3072m', '-Xmx3072m', '-jar', priam_search_path, '--bd', path_to_blast_bin,
						 '--bp', '-n', self.time_stamp, '-i', path_to_input, '-p', path_to_priam_profiles, '--bh', '-o',
						 path_to_output_folder, '--fn']
		try:
			logger.log(prog.logging_levels[logging_level],
					   "Setup process for PRIAM_search.jar: \"" + " ".join(priam_cmd) + "\"")
			self.classifier_processes.add_process_to_workers(self.queue, logging_level, logger_name, priam_cmd,
															 "PRIAM_search.jar")
		except KeyError:
			logger.log(logging.DEBUG, "Setup process for PRIAM_search.jar: \"" + " ".join(priam_cmd) + "\"")
			self.classifier_processes.add_process_to_workers(self.queue, logging_level, logger_name, priam_cmd,
															 "PRIAM_search.jar")

	def run_all_classifiers(self):
		"""Run all classifier subprocesses
		Args:
		Raises:
		Returns:
		"""
		self.classifier_processes.run_all_worker_processes(self.queue)


class Predictions(object):
	"""Object for generating predictions for every classifier
	"""

	def __init__(self, classifier_name):
		self.classifier_name = classifier_name
		self.seq_list = set()
		self.predictions = {}
		self.weights = {}

	def read_weight_file(self, path_to_weight, logging_level, logger_name):
		"""Read in the weight files of a classifer
		Args:
			path_to_weight: path to the weight file of a classifier
			logging_level: The logging level set for reading weight file
			logger_name: The name of the logger for reading weight file
		Raises:
		Returns:
		"""
		logger = logging.getLogger(logger_name)
		try:
			try:
				logger.log(prog.logging_levels[logging_level], "Reading weight file: \"" + path_to_weight + "\"")
			except KeyError:
				logger.log(logging.DEBUG, "Reading weight file: \"" + path_to_weight + "\"")
			with open(path_to_weight, 'r') as ptw:
				for line in ptw:
					info = line.split('\t')
					try:
						ef_name = info[0].strip()
						ef_weight = info[1].strip()
						self.weights.setdefault(ef_name, ef_weight)
					except KeyError:
						continue
		except IOError:
			logger.log(logging.ERROR, "Weight file not found: \"" + path_to_weight + "\"")
			sys.exit(1)

	@staticmethod
	def read_priam_sequence_ec_itr(sequence_ecs_fp):
		"""Iterator to read PRIAM sequenceEC.txt file
		Args:
			sequence_ecs_fp: file pointer to sequenceEC.txt
		Raises:
		Returns:
		"""
		query_id, priam_results = '', []
		for line in sequence_ecs_fp:
			line = line.strip()
			if line.startswith('>'):
				if len(query_id) > 0:
					yield query_id, priam_results
				query_id = re.sub(r'^>', '', re.split(r'\s+|\|', line)[0])
				query_id, priam_results = query_id, []
			else:
				if not line.startswith('#') and len(line) > 0:
					try:
						info = line.split('\t')
						ef_class = info[0].strip()
						try:
							e_value = float(info[2])
							priam_results.append((ef_class, e_value))
						except ValueError:
							continue
					except IndexError:
						continue
				else:
					continue
		if len(query_id) > 0:
			yield query_id, priam_results

	def generate_blast_predictions(self, path_to_blast_weight, path_to_blast_out, logging_level, logger_name):
		"""Read in the blast output and generate predictions
		Args:
			path_to_blast_weight:
			path_to_blast_out: path to the output file of blast
			logging_level: The logging level set for blast prediction
			logger_name: The name of the logger for blast prediction
		Raises:
		Returns:
		"""
		self.read_weight_file(path_to_blast_weight, logging_level, logger_name)
		logger = logging.getLogger(logger_name)
		try:
			try:
				logger.log(prog.logging_levels[logging_level], "Reading blast output: \"" + path_to_blast_out + "\"")
			except KeyError:
				logger.log(logging.DEBUG, "Reading blast output: \"" + path_to_blast_out + "\"")
			with open(path_to_blast_out, 'r') as ptbo:
				for line in ptbo:
					info = [c.strip() for c in line.strip().split('\t')]
					try:
						e_value = float(info[-2])
					except ValueError:
						e_value = float("1" + info[-2])
					query_id = re.split(r'\s+|\|', info[0])[0]
					self.seq_list.add(query_id)
					blast_hits = [bh.strip() for bh in re.split(r'\s+|\|', info[1]) if len(bh) > 0]
					# Current e-value threshold is set to 0.01
					if e_value > float("1e-2") or len(blast_hits) == 0:
						continue
					try:
						cur_e_vals = self.predictions[query_id]
						best_e_val = min(cur_e_vals.values())
						# Keep lowest e-value entries
						if e_value < best_e_val:
							for ef in blast_hits:
								try:
									cur_e_vals[ef] = e_value
								except KeyError:
									self.predictions[query_id].setdefault(ef, e_value)
							for ef in cur_e_vals:
								if ef not in blast_hits:
									self.predictions[query_id].pop(ef, None)
					except KeyError:
						self.predictions.setdefault(query_id, {})
						for ef in blast_hits:
							self.predictions[query_id].setdefault(ef, e_value)
		except IOError:
			logger.log(logging.ERROR, "Blast output not found: \"" + path_to_blast_out + "\"")
			sys.exit(1)

		# Remove non EF Classes from prediction
		for q_id in self.predictions:
			q_predictions = self.predictions[q_id]
			for p in list(q_predictions.keys()):
				if not p.startswith("EF"):
					q_predictions.pop(p, None)

	def generate_priam_predictions(self, path_to_priam_weight, path_to_sequence_ecs, logging_level, logger_name):
		"""Read in the priam output and generate predictions
		Args:
			path_to_priam_weight:
			path_to_sequence_ecs: path to the sequenceEC.txt
			logging_level: The logging level set for priam prediction
			logger_name: The name of the logger for priam prediction
		Raises:
		Returns:
		"""
		self.read_weight_file(path_to_priam_weight, logging_level, logger_name)
		logger = logging.getLogger(logger_name)
		try:
			try:
				logger.log(prog.logging_levels[logging_level],
						   "Reading sequenceECs.txt: \"" + path_to_sequence_ecs + "\"")
			except KeyError:
				logger.log(logging.DEBUG, "Reading sequenceECs.txt: \"" + path_to_sequence_ecs + "\"")
			with open(path_to_sequence_ecs) as ptse:
				for query_id, priam_results in self.read_priam_sequence_ec_itr(ptse):
					self.seq_list.add(query_id)
					if len(priam_results) > 0:
						try:
							for p_res in priam_results:
								try:
									if p_res[1] < self.predictions[query_id][p_res[0]]:
										self.predictions[query_id][p_res[0]] = p_res[1]
								except KeyError:
									self.predictions[query_id].setdefault(p_res[0], p_res[1])
						except KeyError:
							self.predictions.setdefault(query_id, {})
							for p_res in priam_results:
								cur_e_val = self.predictions[query_id].setdefault(p_res[0], p_res[1])
								if p_res[1] < cur_e_val:
									self.predictions[query_id][p_res[0]] = p_res[1]
					else:
						self.predictions.setdefault(query_id, {})

		except IOError:
			logger.log(logging.ERROR, "sequenceECs.txt not found: \"" + path_to_sequence_ecs + "\"")
			sys.exit(1)


class Ensemble(object):
	"""Object for ensemble for final prediction
	"""

	def __init__(self):
		self.classifiers = []
		self.ensemble_predictions = Predictions("Ensemble")
		self.final_predictions = {}
		self.seq_list = set()

	def add_classifier(self, classifier):
		"""Add a Classifier object to the ensemble
		Args:
			classifier: a classifier object
		Raises:
		Returns:
		"""
		self.classifiers.append(classifier)
		for seq_id in sorted(classifier.seq_list):
			self.seq_list.add(seq_id)

	def max_weight_voting(self, logging_level, logger_name):
		"""Perform maximum weight voting on all classifiers added to Ensemble object
			Args:
				logging_level: The logging level set for maximum weight voting
				logger_name: The name of the logger for maximum weight voting
			Raises:
			Returns:
		"""
		logger = logging.getLogger(logger_name)
		try:
			logger.log(prog.logging_levels[logging_level], "Performing Max Weight Voting...")
		except KeyError:
			logger.log(logging.DEBUG, "Performing Max Weight Voting...")
		for classifier in self.classifiers:
			classifier_name = classifier.classifier_name
			prediction = classifier.predictions
			weights = classifier.weights
			for seq_id in prediction:
				if len(prediction[seq_id]) <= 0:
					if self.ensemble_predictions.predictions.get(seq_id) is None:
						self.ensemble_predictions.predictions.setdefault(seq_id, {})
				else:
					for ef_class in prediction[seq_id]:
						ef_weight = float(weights.get(ef_class, 0.00))
						if self.ensemble_predictions.predictions.get(seq_id) is None:
							self.ensemble_predictions.predictions.setdefault(seq_id, {ef_class: ef_weight})
						else:
							weighted_ef_classes = self.ensemble_predictions.predictions[seq_id]
							if weighted_ef_classes.get(ef_class) is None:
								self.ensemble_predictions.predictions[seq_id].setdefault(ef_class, ef_weight)
							elif ef_weight > weighted_ef_classes[ef_class]:
								self.ensemble_predictions.predictions[seq_id][ef_class] = ef_weight

	def absolute_threshold(self, threshold, logging_level, logger_name):
		"""Perform absolute threshold on the ensemble prediction
			Args:
				threshold: the threshold to calculate absolute threshold for the prediction
				logging_level: The logging level set for maximum weight voting
				logger_name: The name of the logger for maximum weight voting
			Raises:
			Returns:
		"""
		logger = logging.getLogger(logger_name)
		try:
			logger.log(prog.logging_levels[logging_level],
					   "Performing Absolute Threshold (" + str(threshold) + ") to votes...")
		except KeyError:
			logger.log(logging.DEBUG, "Performing Absolute Threshold (" + str(threshold) + ") to votes...")
		for seq_id in self.ensemble_predictions.predictions:
			weighted_ef_classes = self.ensemble_predictions.predictions[seq_id]
			if len(weighted_ef_classes) > 0:
				max_weight = max(
					[float(v) for v in weighted_ef_classes.values() if str(v).replace('.', '', 1).isdigit()])
				t = float(max_weight - threshold)
				if t < 0.0:
					t = float(0.0)
				for ef_class in weighted_ef_classes:
					if float(weighted_ef_classes[ef_class]) > t:
						try:
							self.final_predictions[seq_id].append(ef_class)
						except KeyError:
							self.final_predictions.setdefault(seq_id, [ef_class])
			else:
				self.final_predictions.setdefault(seq_id, [])
		for seq_id in [s for s in self.seq_list if s not in self.ensemble_predictions.predictions.keys()]:
			self.final_predictions.setdefault(seq_id, [])
