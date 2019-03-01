import logging.config
import os
import sys


dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(dir_path, 'e2p2'))

import ensemble
import prog

if __name__ == '__main__':
	cur_logger_config = prog.LoggerConfig()
	cur_logger_config.add_new_logger("file", "ensemble.log")
	logging.config.dictConfig(cur_logger_config.dictConfig)
	rc = ensemble.RunClassifiers("2018")
	rc.blast_classifer("/Users/bxue/Projects/GitHub/E2P2/source/blast/ncbi-blast-2.7.1+/bin/blastp",
					   "/Users/bxue/Projects/GitHub/E2P2/source/blast/db/rpsd-4.0-20190108.ef.fa",
					   "/Users/bxue/Documents/Carnegie/PMNProject/pmn_pipeline_test/Araport11_test.fasta",
					   "/Users/bxue/Documents/Carnegie/PMNProject/pmn_pipeline_test/test.blast.Araport11_test.fasta",
					   "1",
					   "INFO",
					   "file")
	rc.priam_classifer("/Users/bxue/Projects/GitHub/E2P2/source/java/jre1.8.0_192.jre/Contents/Home/bin/java",
					   "/Users/bxue/Projects/GitHub/E2P2/source/priam/PRIAM_search.jar",
					   "/Users/bxue/Projects/GitHub/E2P2/source/blast/ncbi-blast-2.7.1+/bin",
					   "/Users/bxue/Projects/GitHub/E2P2/source/priam/profiles",
					   "/Users/bxue/Documents/Carnegie/PMNProject/pmn_pipeline_test/Araport11_test.fasta",
					   "/Users/bxue/Documents/Carnegie/PMNProject/pmn_pipeline_test/test.priam.Araport11_test.fasta",
					   "INFO",
					   "file", False)
	rc.run_all_classifiers()
	print(rc.output_dict)
	bc = ensemble.BlastPredictions("/Users/bxue/Projects/PycharmProjects/E2P2/data/weights/blast", "INFO", "file")
	bc.generate_blast_predictions(
		rc.output_dict["blast"],
		"INFO", "file")
	print(bc.predictions)
	print(bc.weights)
	pc = ensemble.PriamPredictions("/Users/bxue/Projects/PycharmProjects/E2P2/data/weights/priam", "INFO", "file")
	pc.generate_priam_predictions(
		rc.output_dict["priam"],
		"INFO", "file")
	print(pc.predictions)
	print(pc.weights)
	en = ensemble.Ensemble()
	en.add_classifier(bc)
	en.add_classifier(pc)
	en.max_weight_voting("INFO", "file")
	threshold = float(0.5)
	print(en.ensemble_predictions.predictions)
	en.absolute_threshold(threshold, "INFO", "file")
	for seq_id in sorted(en.final_predictions):
		print(seq_id, sorted(en.final_predictions[seq_id]))
