#!/usr/bin/python3

import argparse
import errno
import logging.config
import multiprocessing
import os
import re
import sys
import textwrap
import time

# e2p2 modules
import definitions
import ensemble
import file
import prog


def check_commands_executable(cmd_to_test_list, logging_level, logger_name):
	q = multiprocessing.Queue()
	new_process = prog.RunProcess()
	for cmd in cmd_to_test_list:
		new_process.add_process_to_workers(q, logging_level, logger_name, ['command', '-v', cmd])

	new_process.run_all_worker_processes(q)
	executable_cmds = []
	non_executable_cmds = []
	for result in new_process.run_results:
		cmd = re.sub(r'^command\s-v', '', result[0]).strip()
		retcode = result[1]
		# output = result[2]
		if retcode == 0:
			executable_cmds.append(cmd)
		else:
			non_executable_cmds.append(cmd)
	return executable_cmds, non_executable_cmds


def get_blastp_bin_for_priam(blastp_cmd, logging_level, logger_name):
	q = multiprocessing.Queue()
	new_process = prog.RunProcess()
	new_process.add_process_to_workers(q, logging_level, logger_name, ['which', blastp_cmd])
	new_process.run_all_worker_processes(q)
	try:
		output = new_process.run_results[0][2]
		try:
			bin_folder = os.path.dirname(output.strip())
			bin_folder_files = os.listdir(bin_folder)
			if 'makeprofiledb' in bin_folder_files and 'rpsblast' in bin_folder_files and 'rpstblastn' in bin_folder_files:
				return bin_folder
			else:
				return None
		except AttributeError:
			return None
	except KeyError:
		return None


def run_e2p2_pipeline(time_stamp, input_fasta_path, blastp_cmd, blast_weight, rpsd_db_name_path, blast_output_path,
					  num_threads,
					  java_cmd, priam_search_jar, priam_weight, blast_bin_path, priam_profiles_path, priam_output_path,
					  priam_resume,
					  threshold, e2p2_short_output, e2p2_long_output, e2p2_pf_output, e2p2_orxn_pf_output,
					  e2p2_final_pf_output,
					  logging_level, protein_gene_path=None):
	logger = logging.getLogger(definitions.DEFAULT_LOGGER_NAME)
	# Set up object for running classifiers
	rc = ensemble.RunClassifiers(time_stamp)
	# Run blastp classifier
	rc.blast_classifer(blastp_cmd, rpsd_db_name_path, input_fasta_path, blast_output_path, num_threads, logging_level,
					   definitions.DEFAULT_LOGGER_NAME)
	# Run priam_search classifier
	rc.priam_classifer(java_cmd, priam_search_jar, blast_bin_path, priam_profiles_path, input_fasta_path,
					   priam_output_path, logging_level, definitions.DEFAULT_LOGGER_NAME, priam_resume)
	logger.log(logging.INFO, "Running level-0 classification processes.")
	rc.run_all_classifiers()
	try:
		logger.log(logging.INFO, "Compiling predictions.")
		# Read in prediction results
		bc = ensemble.Predictions("Blast")
		bc.generate_blast_predictions(blast_weight, rc.output_dict["blast"], logging_level,
									  definitions.DEFAULT_LOGGER_NAME)
		pc = ensemble.Predictions("Priam")
		pc.generate_priam_predictions(priam_weight, rc.output_dict["priam"], logging_level,
									  definitions.DEFAULT_LOGGER_NAME)
		# Set up object for ensemble
		en = ensemble.Ensemble()
		# Add classifier results to ensemble object
		en.add_classifier(bc)
		en.add_classifier(pc)
		logger.log(logging.INFO, "Computing ensemble predictions.")
		# Preform max weight voting on all classifier results
		en.max_weight_voting(logging_level, definitions.DEFAULT_LOGGER_NAME)
		# Preform absolute threshold on classifiers voting results
		en.absolute_threshold(threshold, logging_level, definitions.DEFAULT_LOGGER_NAME)
		# Set up object for writing E2P2 output
		e2p2 = file.E2P2files(en.final_predictions)
		logger.log(logging.INFO, "Preparing results files.")
		# Add classifer predictions to E2P2 Output for detailed results
		e2p2.add_predictions_of_classifer(bc)
		e2p2.add_predictions_of_classifer(pc)
		e2p2.read_efmap(ef_map_path, logging_level, definitions.DEFAULT_LOGGER_NAME)
		# Write Outputs
		e2p2.write_short_results(definitions.DEFAULT_ENSEMBLE_METHOD + " (" + str(threshold) + ")",
								 e2p2_short_output, logging_level, definitions.DEFAULT_LOGGER_NAME)
		e2p2.write_long_results(definitions.DEFAULT_ENSEMBLE_METHOD + " (" + str(threshold) + ")",
								e2p2_long_output, logging_level, definitions.DEFAULT_LOGGER_NAME)
		e2p2.write_pf_results(e2p2_pf_output, logging_level, definitions.DEFAULT_LOGGER_NAME)
		e2p2.write_orxn_pf_results(e2p2_orxn_pf_output,
								   definitions.EC_SUPERSEDED_MAP, definitions.METACYC_RXN_MAP,
								   definitions.OFFICIAL_EC_METACYC_RXN_MAP,
								   definitions.TO_REMOVE_NON_SMALL_MOLECULE_METABOLISM, logging_level,
								   definitions.DEFAULT_LOGGER_NAME)
		if protein_gene_path is not None:
			e2p2.write_final_pf_results(e2p2_final_pf_output,
										definitions.EC_SUPERSEDED_MAP, definitions.METACYC_RXN_MAP,
										definitions.OFFICIAL_EC_METACYC_RXN_MAP,
										definitions.TO_REMOVE_NON_SMALL_MOLECULE_METABOLISM, protein_gene_path,
										logging_level, definitions.DEFAULT_LOGGER_NAME)
		logger.log(logging.INFO, "Operation complete.")
		logger.log(logging.INFO, "Main results are in the file: %s" % e2p2_short_output)
		logger.log(logging.INFO, "Detailed results are in the file: %s" % e2p2_long_output)
		if protein_gene_path is not None:
			logger.log(logging.INFO, "To build PGDB, use .pf file: %s" % e2p2_final_pf_output)
		else:
			logger.log(logging.INFO, "To build PGDB, use .pf file: %s" % e2p2_orxn_pf_output)
	except KeyError as classifer_missing:
		logger.log(logging.ERROR, "Missing classifer result file(s): " + str(classifer_missing))
		sys.exit(1)


if __name__ == '__main__':
	name = 'run_pipeline.py'
	description = '''
	    Runs the Ensemble Enzyme Prediction Pipeline (E2P2) on a set of input protein sequences,
	    outputting enzyme functional annotations in the forms of EC numbers or MetaCyc reaction
	    IDs for any predicted enzyme sequences.
	    '''
	notes = '''
	    - Input protein sequences should be in FASTA format.
	    - Headers in the FASTA file should begin with the sequence ID followed by a space.
	    - Intermediate results files can be found in a temporary directory of its own subdirectory labeled with a
	      date and time stamp.
	'''

	parser = argparse.ArgumentParser(prog=name, description=description,
									 formatter_class=argparse.RawTextHelpFormatter,
									 epilog=textwrap.dedent(notes)
									 )
	# Argument for IO
	parser.add_argument("--input", "-i", dest="input_file", required=True, type=definitions.PathType('file'),
						help="Path to input protein sequences file")
	parser.add_argument("--output", "-o", dest="output_path", type=definitions.PathType('parent'),
						help="Path to output file. By Default would be in the same folder of the input.")
	# Argument for third party software
	parser.add_argument("--blastp", "-b", dest="blastp_cmd", required=True, help=textwrap.dedent(
		"Command of or path to BLAST+ \"blastp\".\nDownload Link:" + definitions.BLAST_PLUS_DOWNLOAD_LINK))
	parser.add_argument("--num_threads", "-n", dest="num_threads", type=int, required=False, default="1",
						help="Number of threads to run \"blastp\". Default is 1")
	parser.add_argument("--java", "-j", dest="java_cmd", required=True, help=textwrap.dedent(
		"Command of or path to \"java\".\nDownload Link:" + definitions.JAVA_8_DOWNLOAD_LINK))
	parser.add_argument("--priam_search", "-ps", dest="priam_search", required=True, type=definitions.PathType('file'),
						help=textwrap.dedent(
							"Path to \"PRIAM_search.jar\".\nDownload Link:" + definitions.PRIAM_SEARCH_LINK))
	parser.add_argument("--priam_resume", "-pr", dest="priam_resume", action='store_true',
						help="Whether or not to resume a found PRIAM_search.jar process.")
	parser.add_argument("--blast_bin", "-bb", dest="blast_bin", type=definitions.PathType('blastbin'),
						help=textwrap.dedent(
							"Command of or path to BLAST+ bin folder.\nDownload Link:" + definitions.BLAST_PLUS_DOWNLOAD_LINK))
	# Arguments for databases
	# parser.add_argument("--data", "-d", dest="data_path", type=definitions.PathType('dir'),
	# 					help=textwrap.dedent("Path to data path for E2P2 release folder"))
	parser.add_argument("--blast_weight", "-bw", dest="blast_weight", type=definitions.PathType('file'),
						help=textwrap.dedent("Path to blast weight for the blast classifier"))
	parser.add_argument("--rpsd", "-r", dest="rpsd_db", required=True, type=definitions.PathType('blastdb'),
						help=textwrap.dedent("Path to rpsd database name.\nFor example, \"/PATH/TO/FOLDER/rpsd.fa\", "
											 "where you can find the following files in /PATH/TO/FOLDER:\n"
											 "rpsd.fa.phr; rpsd.fa.pin; rpsd.fa.psq"))
	parser.add_argument("--priam_profile", "-pp", dest="priam_profile", required=True, type=definitions.PathType('profiles'),
						help=textwrap.dedent("Path to PRIAM profiles.\nFor example, \"/PATH/TO/FOLDER/profiles\", "
											 "where you can find the following in /PATH/TO/FOLDER:\n"
											 "files: annotation_rules.xml; genome_rules.xml\n"
											 "PROFILES: Folder contains \"LIBRARY\" folder and multiple \".chk\" files."))
	parser.add_argument("--priam_weight", "-pw", dest="priam_weight", type=definitions.PathType('file'),
						help=textwrap.dedent("Path to blast weight for the priam classifier"))
	parser.add_argument("--efmap", "-e", dest="ef_map", type=definitions.PathType('file'),
						help="Path to efclasses.mapping file.")
	parser.add_argument("--threshold", "-th", dest="threshold", type=float, default="0.5",
						help="Threshold for voting results. Default is 0.5.")
	parser.add_argument("--temp_folder", "-tf", dest="temp_folder", type=definitions.PathType('dir'),
						help="Specify the location of the temp folder. By default would be in the same directory of the output.")
	parser.add_argument("--log", "-l", dest="log_path", type=definitions.PathType('parent'),
						help="Specify the location of the log file. By default would be \"runE2P2.log\" in the temp folder.")
	parser.add_argument("--protein_gene", "-pg", dest="protein_gene_path", type=definitions.PathType('file'),
						help="Provide a protein to gene map. This will be used to generate a splice variant removed fasta file and output our final version of e2p2.")
	verbose_message = '''Verbose level of log output. Default is 0.
		0: only step information are logged
		1: all information are logged
		'''
	parser.add_argument("--verbose", "-v", dest="verbose", default="0", choices=["0", "1"],
						help=textwrap.dedent(verbose_message))
	args = parser.parse_args()
	timestamp = str(time.time())
	input_path_folder = os.path.dirname(args.input_file)
	input_file_name = os.path.basename(args.input_file)
	input_file_path = args.input_file
	file.E2P2files(None).check_fasta_header(input_file_path, definitions.DEFAULT_LOGGER_NAME)
	if args.output_path is None:
		output_path = input_file_path + definitions.DEFAULT_OUTPUT_SUFFIX
	else:
		output_path = args.output_path
	output_folder = os.path.dirname(output_path)
	if args.temp_folder is None:
		temp_folder = args.input_file + '.' + timestamp
	else:
		temp_folder = args.temp_folder
	create_temp_folder_flag = False
	try:
		os.mkdir(temp_folder)
		create_temp_folder_flag = True
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise
		pass
	if args.log_path is None:
		log_path = os.path.join(temp_folder, 'rune2p2.' + timestamp + '.log')
	else:
		log_path = args.log_path
	if args.verbose == "0":
		logger_handler_level = "DEBUG"
	else:
		logger_handler_level = "INFO"
	cur_logger_config = prog.LoggerConfig()
	if os.path.isfile(os.path.realpath(log_path)):
		cur_logger_config.add_new_logger(definitions.DEFAULT_LOGGER_NAME, log_path, logger_handler_mode='a')
		logging.config.dictConfig(cur_logger_config.dictConfig)
		logger = logging.getLogger(definitions.DEFAULT_LOGGER_NAME)
		if create_temp_folder_flag:
			logger.log(prog.logging_levels[logger_handler_level], "Temp folder created at path %s." % temp_folder)
		else:
			logger.log(prog.logging_levels[logger_handler_level], "Using path %s as temp folder." % temp_folder)
		logger.log(logging.WARNING, "Log file %s exists, will append to it..." % log_path)
	else:
		cur_logger_config.add_new_logger(definitions.DEFAULT_LOGGER_NAME, log_path)
		logging.config.dictConfig(cur_logger_config.dictConfig)
		logger = logging.getLogger(definitions.DEFAULT_LOGGER_NAME)
		if create_temp_folder_flag:
			logger.log(prog.logging_levels[logger_handler_level], "Temp folder created at path %s." % temp_folder)
		else:
			logger.log(prog.logging_levels[logger_handler_level], "Using path %s as temp folder." % temp_folder)

	cmd_list = [args.blastp_cmd, args.java_cmd]
	check_commands_executable(cmd_list, logger_handler_level, definitions.DEFAULT_LOGGER_NAME)

	if args.blast_bin is None:
		blast_bin_path = get_blastp_bin_for_priam(args.blastp_cmd, logger_handler_level,
												  definitions.DEFAULT_LOGGER_NAME)
		if blast_bin_path is None:
			logger.log(logging.ERROR, "Cannot verify BLAST+ bin folder from %s." % args.blastp_cmd)
			sys.exit(1)
	else:
		blast_bin_path = args.blast_bin

	if args.rpsd_db is None:
		logger.log(logging.ERROR, "RPSD blast database not found from path %s." % args.rpsd_db)
		sys.exit(1)
	if args.priam_profile is None:
		logger.log(logging.ERROR, "PRIAM profiles not found from path %s." % args.blastp_cmd)
		sys.exit(1)

	if args.ef_map is None:
		ef_map_path = definitions.EF_CLASS_MAP
	else:
		ef_map_path = args.ef_map
	if args.blast_weight is None:
		blast_weight_path = definitions.BLAST_WEIGHT
	else:
		blast_weight_path = args.blast_weight
	if args.priam_weight is None:
		priam_weight_path = definitions.PRIAM_WEIGHT
	else:
		priam_weight_path = args.priam_weight
	e2p2_long_output = output_path + definitions.DEFAULT_LONG_OUTPUT_SUFFIX
	e2p2_pf_output = output_path + definitions.DEFAULT_PF_OUTPUT_SUFFIX
	e2p2_orxn_pf_output = output_path + definitions.DEFAULT_ORXN_PF_OUTPUT_SUFFIX
	e2p2_final_pf_output = output_path + definitions.DEFAULT_FINAL_PF_OUTPUT_SUFFIX
	if os.path.isfile(output_path):
		logger.log(logging.WARNING, "Output file %s exists, will overwrite..." % output_path)
	if os.path.isfile(e2p2_long_output):
		logger.log(logging.WARNING, "Output file %s exists, will overwrite..." % e2p2_long_output)
	if os.path.isfile(e2p2_pf_output):
		logger.log(logging.WARNING, "Output file %s exists, will overwrite..." % e2p2_pf_output)
	if os.path.isfile(e2p2_orxn_pf_output):
		logger.log(logging.WARNING, "Output file %s exists, will overwrite..." % e2p2_orxn_pf_output)

	if args.protein_gene_path is not None:
		input_file_path = file.E2P2files(None).remove_splice_variants(input_file_path, output_folder,
																	  args.protein_gene_path, logger_handler_level,
																	  definitions.DEFAULT_LOGGER_NAME)
		if os.path.isfile(e2p2_final_pf_output):
			logger.log(logging.WARNING, "Output file %s exists, will overwrite..." % e2p2_final_pf_output)

	blastp_output_path = os.path.join(temp_folder, 'blast.%s.%s' % (input_file_name, timestamp))
	if os.path.isfile(blastp_output_path):
		logger.log(logging.WARNING, "Path %s for blastp result exists, will overwrite..." % blastp_output_path)
	priam_output_path = os.path.join(temp_folder, "PRIAM_" + timestamp)
	if os.path.isdir(os.path.join(temp_folder, "PRIAM_" + timestamp)):
		if args.priam_resume is True:
			logger.log(logging.WARNING, "Path %s for PRIAM result exists, will resume..." % priam_output_path)
		else:
			logger.log(logging.WARNING, "Path %s for PRIAM result exists, will overwrite..." % priam_output_path)
	run_e2p2_pipeline(timestamp, input_file_path, args.blastp_cmd, blast_weight_path, args.rpsd_db, blastp_output_path,
					  args.num_threads,
					  args.java_cmd, args.priam_search, priam_weight_path, blast_bin_path, args.priam_profile,
					  temp_folder,
					  args.priam_resume, float(args.threshold), output_path, e2p2_long_output, e2p2_pf_output,
					  e2p2_orxn_pf_output,
					  e2p2_final_pf_output, logger_handler_level, args.protein_gene_path)
	logger.log(logging.INFO, "Intermediate files are in the directory: %s" % temp_folder)
