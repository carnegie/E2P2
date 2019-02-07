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
import os
import re
import sys

from operator import itemgetter

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path)

import prog


class Classifier(object):
	def __init__(self, classifier_name):
		self.classifier_name = classifier_name
		self.seq_list = set()
		self.predictions = {}
		self.weights = {}

	def read_weight_file(self, path_to_weight, logging_level, logger_name):
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


class BlastClassifier(Classifier):
	def __init__(self, path_to_blast_weight, logging_level, logger_name):
		super().__init__("Blast")
		self.read_weight_file(path_to_blast_weight, logging_level, logger_name)

	def generate_blast_predictions(self, path_to_blast_out, logging_level, logger_name):
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

		# Remove non EF Classes from prediction
		for q_id in self.predictions:
			q_predictions = self.predictions[q_id]
			for p in list(q_predictions.keys()):
				if not p.startswith("EF"):
					q_predictions.pop(p, None)


class PriamClassifier(Classifier):
	def __init__(self, path_to_blast_weight, logging_level, logger_name):
		super().__init__("Priam")
		self.read_weight_file(path_to_blast_weight, logging_level, logger_name)

	@staticmethod
	def read_priam_sequence_ec_itr(sequence_ecs_fp):
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

	def generate_priam_predictions(self, path_to_sequence_ecs, logging_level, logger_name):
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


class Ensemble(object):
	def __init__(self):
		self.classifiers = []
		self.ensemble_predictions = Classifier("Ensemble")
		self.final_predictions = {}
		self.seq_list = set()

	def add_classifier(self, classifier):
		self.classifiers.append(classifier)
		for seq_id in sorted(classifier.seq_list):
			self.seq_list.add(seq_id)

	def max_weight_voting(self):
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
						ef_weight = weights.get(ef_class, 0.00)
						if self.ensemble_predictions.predictions.get(seq_id) is None:
							self.ensemble_predictions.predictions.setdefault(seq_id, {ef_class: ef_weight})
						else:
							weighted_ef_classes = self.ensemble_predictions.predictions[seq_id]
							if weighted_ef_classes.get(ef_class) is None:
								self.ensemble_predictions.predictions[seq_id].setdefault(ef_class, ef_weight)
							elif ef_weight > weighted_ef_classes[ef_class]:
								self.ensemble_predictions.predictions[seq_id][ef_class] = ef_weight

	def absolute_threshold(self, threshold):
		for seq_id in self.ensemble_predictions.predictions:
			weighted_ef_classes = self.ensemble_predictions.predictions[seq_id]
			if len(weighted_ef_classes) > 0:
				max_weight = max([float(v) for v in weighted_ef_classes.values() if v.replace('.','',1).isdigit()])
				t = float(max_weight - threshold)
				# print(seq_id, max_weight, t)
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




if __name__ == '__main__':
	cur_logger_config = prog.LoggerConfig()
	cur_logger_config.add_new_logger("file", "file", "ensemble.log")
	logging.config.dictConfig(cur_logger_config.dictConfig)

	bc = BlastClassifier("/Users/bxue/Projects/PycharmProjects/E2P2/data/weights/blast", "INFO", "file")
	bc.generate_blast_predictions(
		"/Users/bxue/Projects/GitHub/E2P2/run/Araport11_test.fasta.1549074341.11/blast.1549074341.11",
		"INFO", "file")
	# print(bc.predictions["AT1G01990.1"])
	# print(bc.weights)
	pc = PriamClassifier("/Users/bxue/Projects/PycharmProjects/E2P2/data/weights/priam", "INFO", "file")
	pc.generate_priam_predictions(
		"/Users/bxue/Projects/GitHub/E2P2/run/Araport11_test.fasta.1549074341.11/PRIAM_1549074341.11/ANNOTATION/sequenceECs.txt",
		"INFO", "file")
	# print(pc.predictions)
	# print(pc.weights)
	en = Ensemble()
	en.add_classifier(bc)
	en.add_classifier(pc)
	en.max_weight_voting()
	threshold = float(0.5)
	# print(en.ensemble_predictions.predictions)
	en.absolute_threshold(threshold)
	for seq_id in sorted(en.final_predictions):
		print(seq_id, sorted(en.final_predictions[seq_id]))
