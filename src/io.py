import logging
import logging.config
import os
import re
import sys

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path)

import prog, ensemble


class E2P2Output(object):
	def __init__(self, final_prediction):
		self.final_prediction = final_prediction
		self.efmap = {}

	def read_efmap(self, ef_map_path, logging_level, logger_name):
		logger = logging.getLogger(logger_name)
		try:
			logger.log(prog.logging_levels[logging_level], "Loading efmap: \"" + ef_map_path + "\"")
		except KeyError:
			logger.log(logging.DEBUG, "Loading efmap: \"" + ef_map_path + "\"")
		with open(ef_map_path, 'r') as emp:
			for line in emp:
				info = line.split('\t')
				try:
					ef_class = info[0].strip()
					rxn_id = info[1].strip()
					self.efmap.setdefault(ef_class, rxn_id)
				except IndexError:
					continue


if __name__ == '__main__':
	cur_logger_config = prog.LoggerConfig()
	cur_logger_config.add_new_logger("file", "file", "ensemble.log")
	logging.config.dictConfig(cur_logger_config.dictConfig)

	bc = ensemble.BlastClassifier("/Users/bxue/Projects/PycharmProjects/E2P2/data/weights/blast", "INFO", "file")
	bc.generate_blast_predictions(
		"/Users/bxue/Projects/GitHub/E2P2/run/Araport11_test.fasta.1549074341.11/blast.1549074341.11",
		"INFO", "file")
	# print(bc.predictions["AT1G01990.1"])
	# print(bc.weights)
	pc = ensemble.PriamClassifier("/Users/bxue/Projects/PycharmProjects/E2P2/data/weights/priam", "INFO", "file")
	pc.generate_priam_predictions(
		"/Users/bxue/Projects/GitHub/E2P2/run/Araport11_test.fasta.1549074341.11/PRIAM_1549074341.11/ANNOTATION/sequenceECs.txt",
		"INFO", "file")
	# print(pc.predictions)
	# print(pc.weights)
	en = ensemble.Ensemble()
	en.add_classifier(bc)
	en.add_classifier(pc)
	en.max_weight_voting("INFO", "file")
	threshold = float(0.5)
	# print(en.ensemble_predictions.predictions)
	en.absolute_threshold(threshold, "INFO", "file")
	# for seq_id in sorted(en.final_predictions):
	# 	print(seq_id, sorted(en.final_predictions[seq_id]))

	e2p2_output = E2P2Output(en.final_predictions)
	e2p2_output.read_efmap("/Users/bxue/Projects/PycharmProjects/E2P2/data/maps/fcmap", "INFO", "file")
	print(e2p2_output.efmap)


