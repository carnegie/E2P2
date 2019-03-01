import logging
import logging.config
import os
import sys


dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(dir_path, 'e2p2'))

import ensemble
import prog
import file

if __name__ == '__main__':
	cur_logger_config = prog.LoggerConfig()
	cur_logger_config.add_new_logger("file", "io.log")
	logging.config.dictConfig(cur_logger_config.dictConfig)

	bc = ensemble.BlastPredictions("/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/UpdateWeights/weights/blast", "INFO", "file")
	bc.generate_blast_predictions(
		"/Users/bxue/Documents/Carnegie/PMNProject/pmn_pipeline_test/test.blast.Araport11_test.fasta",
		"INFO", "file")

	pc = ensemble.PriamPredictions("/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/UpdateWeights/weights/priam", "INFO", "file")
	pc.generate_priam_predictions(
		"/Users/bxue/Documents/Carnegie/PMNProject/pmn_pipeline_test/test.priam.Araport11_test.fasta/PRIAM_2018/ANNOTATION/sequenceECs.txt",
		"INFO", "file")

	en = ensemble.Ensemble()
	en.add_classifier(bc)
	en.add_classifier(pc)
	en.max_weight_voting("INFO", "file")
	threshold = float(0.5)
	en.absolute_threshold(threshold, "INFO", "file")

	e2p2_output = file.E2P2Output(en.final_predictions)
	e2p2_output.add_predictions_of_classifer(bc)
	e2p2_output.add_predictions_of_classifer(pc)
	e2p2_output.read_efmap("/Users/bxue/Projects/PycharmProjects/E2P2/data/maps/efclasses.mapping", "INFO", "file")
	# e2p2_output.write_short_results("Maximum weight with absolute threshold (0.5)", "/Users/bxue/Documents/Carnegie/PMNProject/pmn_pipeline_test/test.blast.Araport11_test.fasta.e2p2v4", "INFO", "file")
	# e2p2_output.write_long_results("Maximum weight with absolute threshold (0.5)", "/Users/bxue/Documents/Carnegie/PMNProject/pmn_pipeline_test/test.blast.Araport11_test.fasta.e2p2v4.long", "INFO", "file")
	# e2p2_output.write_pf_results("/Users/bxue/Documents/Carnegie/PMNProject/pmn_pipeline_test/test.blast.Araport11_test.fasta.e2p2v4.pf", "INFO", "file")
	e2p2_output.write_orxn_pf_results("/Users/bxue/Documents/Carnegie/PMNProject/pmn_pipeline_test/test.blast.Araport11_test.fasta.e2p2v4.orxn.pf",
									  "/Users/bxue/Projects/PycharmProjects/E2P2/data/maps/pf-EC-superseded.mapping", "/Users/bxue/Projects/PycharmProjects/E2P2/data/maps/pf-metacyc-RXN-EC.mapping",
									  "/Users/bxue/Projects/PycharmProjects/E2P2/data/maps/pf-official-EC-metacyc-RXN.mapping",
									  "/Users/bxue/Projects/PycharmProjects/E2P2/data/maps/pf-to-remove-non-small-molecule-metabolism.mapping",
									  "INFO", "file")



