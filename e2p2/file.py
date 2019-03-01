import datetime
import logging
import logging.config
import os
import sys

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path)

import prog


class E2P2Output(object):
	def __init__(self, final_prediction):
		self.final_prediction = final_prediction
		self.efmap = {}
		self.predictions_of_classifers = {}
		self.protein_to_gene_map = {}

	@staticmethod
	def read_pf_maps(file_path, key_idx, val_idx):
		mapping_dict = {}
		with open(file_path, 'r') as fp:
			for line in fp:
				if line.startswith('#') or len(line.strip()) == 0:
					continue
				info = line.strip().split('\t')
				try:
					mapping_dict[info[key_idx]].append(info[val_idx])
				except KeyError:
					mapping_dict.setdefault(info[key_idx], [info[val_idx]])
				except IndexError:
					continue
		return mapping_dict

	def read_efmap(self, ef_map_path, logging_level, logger_name):
		logger = logging.getLogger(logger_name)
		try:
			logger.log(prog.logging_levels[logging_level], "Loading efmap: \"" + ef_map_path + "\"")
		except KeyError:
			logger.log(logging.DEBUG, "Loading efmap: \"" + ef_map_path + "\"")
		self.efmap = self.read_pf_maps(ef_map_path, 0, 1)


	def add_predictions_of_classifer(self, classifier):
		try:
			self.predictions_of_classifers.setdefault(classifier.classifier_name, classifier)
		except AttributeError:
			pass

	def write_short_results(self, ensemble_method, output_path, logging_level, logger_name):
		cur_time = datetime.datetime.now()
		logger = logging.getLogger(logger_name)
		header = "# Result Generate time:  %s\n" \
				 "# Ensemble method used:  %s\n" % (cur_time, ensemble_method)
		with open(output_path, 'w') as op:
			op.write(header)
			try:
				for seq_id in sorted(self.final_prediction.keys()):
					predictions = self.final_prediction[seq_id]
					write_line = seq_id + '\t'
					if len(predictions) == 0:
						write_line += 'NA\n'
					else:
						write_line += '|'.join(predictions) + '\n'
					op.write(write_line)
			except AttributeError:
				logger.log(logging.ERROR, "FinalPrediction Incorrect.")
		try:
			logger.log(prog.logging_levels[logging_level], "Results written to: \"" + output_path + "\"")
		except KeyError:
			logger.log(logging.DEBUG, "Results written to: \"" + output_path + "\"")

	def write_long_results(self, ensemble_method, output_path, logging_level, logger_name):
		cur_time = datetime.datetime.now()
		logger = logging.getLogger(logger_name)
		header = "# Result Generate time:  %s\n" \
				 "# Ensemble method used:  %s\n" % (cur_time, ensemble_method)
		with open(output_path, 'w') as op:
			op.write(header)
			try:
				for seq_id in sorted(self.final_prediction.keys()):
					predictions = self.final_prediction[seq_id]
					write_line = '>' + seq_id + '\t'
					if len(predictions) == 0:
						write_line += 'NA\n'
					else:
						write_line += '|'.join(predictions) + '\n'
					op.write(write_line)
					for classifier in self.predictions_of_classifers:
						classifier_line = classifier + '\t'
						weights = self.predictions_of_classifers[classifier].weights
						try:
							classifier_predictions = set(self.predictions_of_classifers[classifier].predictions[seq_id]).intersection(set(predictions))
							weighted_predictions = []
							for result in sorted(classifier_predictions):
								try:
									weighted_predictions.append(result + ' (' + weights[result] + ')')
								except KeyError:
									continue
							if len(weighted_predictions) == 0:
								classifier_line = ''
							else:
								classifier_line += '|'.join(weighted_predictions) + '\n'
							op.write(classifier_line)
						except KeyError:
							continue
					op.write('\n')
			except AttributeError:
				logger.log(logging.ERROR, "FinalPrediction Incorrect.")
		try:
			logger.log(prog.logging_levels[logging_level], "Results (long) written to: \"" + output_path + "\"")
		except KeyError:
			logger.log(logging.DEBUG, "Results (long) written to: \"" + output_path + "\"")

	def write_pf_results(self, output_path, logging_level, logger_name):
		logger = logging.getLogger(logger_name)
		if len(self.efmap) > 0:
			with open(output_path, 'w') as op:
				for seq_id in sorted(self.final_prediction):
					predictions = self.final_prediction[seq_id]
					if len(predictions) > 0:
						op.write("ID\t%s\nNAME\t%s\nPRODUCT-TYPE\tP\n" % (seq_id, seq_id))
						for ef_class in sorted(predictions):
							try:
								mapped_ids = ["METACYC\t" + i if "RXN" in i else "EC\t" + i for i in self.efmap[ef_class]]
								op.write('\n'.join(mapped_ids) + '\n')
							except KeyError:
								logger.log(logging.ERROR, "EF Class: \"" + ef_class + "\" assigned to \"" + seq_id + "\" not found in map.")
						op.write("//\n")
			try:
				logger.log(prog.logging_levels[logging_level], "Results (pf) written to: \"" + output_path + "\"")
			except KeyError:
				logger.log(logging.DEBUG, "Results (pf) written to: \"" + output_path + "\"")
		else:
			logger.log(logging.ERROR, "EF map not found, please set up with read_efmap()...")

	def write_orxn_pf_results(self, output_path, ec_superseded, metacyc_rxn_ec, official_ec_metacyc_rxn,
							  to_remove_non_small_molecule_metabolism, logging_level, logger_name):
		logger = logging.getLogger(logger_name)

		ec_superseded_dict = self.read_pf_maps(ec_superseded, 2, 0)
		metacyc_rxn_ec_dict = self.read_pf_maps(metacyc_rxn_ec, 1, 0)
		official_ec_metacyc_rxn_dict = self.read_pf_maps(official_ec_metacyc_rxn, 0, 1)
		# metacyc_sub_reactions_dict = self.read_pf_maps(metacyc_sub_reactions, 0, 1)
		to_remove_non_small_molecule_metabolism_list = sorted(self.read_pf_maps(to_remove_non_small_molecule_metabolism, 0, 0).keys())

		if len(self.efmap) > 0:
			with open(output_path, 'w') as op:
				for seq_id in sorted(self.final_prediction):
					predictions = self.final_prediction[seq_id]
					if len(predictions) > 0:
						metacyc_ids = set()
						for ef_class in sorted(predictions):
							try:
								metacyc_ids.update([i for i in self.efmap[ef_class] if "RXN" in i])
								ec_ids = [i for i in self.efmap[ef_class] if "RXN" not in i]
								ec_ids_updated = set()
								for ec in ec_ids:
									try:
										ec_superseded_ids = [e.replace('EC-', '') for e in ec_superseded_dict["EC-" + ec]]
										ec_ids_updated.update(ec_superseded_ids)
									except KeyError:
										ec_ids_updated.add(ec)
								for ec_updated in sorted(ec_ids_updated):
									if ec_updated in official_ec_metacyc_rxn_dict:
										metacyc_ids.update(official_ec_metacyc_rxn_dict[ec_updated])
									elif ec_updated in metacyc_rxn_ec_dict:
										metacyc_ids.update(metacyc_rxn_ec_dict[ec_updated])
							except KeyError:
								logger.log(logging.ERROR, "EF Class: \"" + ef_class + "\" assigned to \"" + seq_id + "\" not found in map.")
						pf_ids = sorted(['METACYC\t' + m for m in metacyc_ids if m not in to_remove_non_small_molecule_metabolism_list])
						if len(pf_ids) > 0:
							op.write("ID\t%s\nNAME\t%s\nPRODUCT-TYPE\tP\n" % (seq_id, seq_id))
							op.write('\n'.join(pf_ids) + '\n')
							op.write("//\n")
			try:
				logger.log(prog.logging_levels[logging_level], "Results (orxn.pf) written to: \"" + output_path + "\"")
			except KeyError:
				logger.log(logging.DEBUG, "Results (orxn.pf) written to: \"" + output_path + "\"")
		else:
			logger.log(logging.ERROR, "EF map not found, please set up with read_efmap()...")
