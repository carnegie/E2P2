import logging.config
import os
import re
import sys
import time
from itertools import chain, combinations

from pipeline import definitions, prog, ensemble_test, file


print(sys.path)


def read_pf_maps(file_path, key_idx, val_idx):
    """Read Tab-delimited Files
    Args:
        file_path: Path to mapping file
        key_idx: Index for keys
        val_idx: Index for values
    Raises: KeyError, IndexError
    Returns:
        mapping_dict: Dict containing the map
    """
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


def main_test():
    ef_map_path = '/Users/bxue/Documents/Carnegie/SourceData/PMN/RPSD/release_2019-03-07/maps/efclasses.mapping'
    ef_map = read_pf_maps(ef_map_path, 0, 1)

    priam_output = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/priam/PRIAM_20201118_hold/ANNOTATION/' \
                   'sequenceECs.txt'
    # priam_e2p2 = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/UpdateWeights/classifiers/e2p2v4/' \
    #                'rpsd-4.2-20190307.ef-hold.e2p2v4.Priam'
    priam_weight = '/Users/bxue/Documents/Carnegie/SourceData/PMN/RPSD/release_2019-03-07/weights/priam'

    blast_output = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/blast/' \
                   'rpsd-4.2-20190307.ef.fasta.rpsd-4.2-20190307.ef-hold.lst.subset.fa.out'
    # blast_e2p2 = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/UpdateWeights/classifiers/e2p2v4/' \
    #              'rpsd-4.2-20190307.ef-hold.e2p2v4.BLAST'
    blast_weight = '/Users/bxue/Documents/Carnegie/SourceData/PMN/RPSD/release_2019-03-07/weights/blast'

    deepec_output = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/deepec/rpsd_split/' \
                    'deepec_output_08_23_21_v1/deepec.train_test.rpsd_true.deepec_only.hold.v1.instance.prediction.txt'
    deepec_weight = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/deepec/rpsd_split/' \
                    'deepec.weight.08_10_21.txt'
    output_folder = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/ensemble'

    metacyc_to_ec_path = '/Users/bxue/Documents/Carnegie/SourceData/PMN/RPSD/release_2019-03-07/maps/' \
                         'combined.metacyc_2_ec.mapping'
    metacyc_to_ec = read_pf_maps(metacyc_to_ec_path, 0, 1)
    input_file_path = '/Users/bxue/Documents/Carnegie/PMNProject/RPSDv4.2/train_test/rpsd-4.2-20190307.ef-hold.lst.fa'

    # input_fasta_path = os.path.basename(input_file_path)
    input_dict = {}
    with open(input_file_path, 'r') as fp:
        for header, _ in file.E2P2files(None).read_fasta(fp):
            header = header.lstrip('>')
            info = header.split('|')
            seq_id = info[0]
            efcls = set()
            for ef in info[1:]:
                if ef.strip() != '':
                    efcls.update(ef_map[ef.strip()])
            mappedecfls = set()
            for ef in sorted(efcls):
                try:
                    ef_class_list = [ec for ec in sorted(metacyc_to_ec[ef])
                                     if re.search('\d+\.\d+\.\d+\.\d+', ec)]
                    mappedecfls.update(ef_class_list)
                except KeyError:
                    mappedecfls.add(ef)
            input_dict.setdefault(seq_id, set(mappedecfls))
    log_path = os.path.join(output_folder, 'log.txt')

    logger_handler_level = "DEBUG"
    cur_logger_config = prog.LoggerConfig()
    cur_logger_config.add_new_logger(definitions.DEFAULT_LOGGER_NAME, log_path)
    logging.config.dictConfig(cur_logger_config.dictConfig)
    logger = logging.getLogger(definitions.DEFAULT_LOGGER_NAME)

    time_stamp = str(time.time())
    # rc = ensemble_test.RunClassifiers(time_stamp)

    ss = ['deepec', 'blast', 'priam']
    name_comb = list(chain(*map(lambda x: combinations(ss, x), range(0, len(ss) + 1))))
    try:
        logger.log(logging.INFO, "Compiling predictions.")
        # Read in prediction results
        dc = ensemble_test.Predictions("DeepEC")
        dc.generate_deepec_predictions(deepec_weight, deepec_output, logger_handler_level,
                                       definitions.DEFAULT_LOGGER_NAME)
        bc = ensemble_test.Predictions("Blast")
        bc.generate_blast_predictions(blast_weight, blast_output, float("1e-2"), logger_handler_level,
                                      definitions.DEFAULT_LOGGER_NAME, efmap=ef_map)
        # bc.read_e2p2_output(blast_weight, blast_e2p2, logger_handler_level, definitions.DEFAULT_LOGGER_NAME,
        #                     efmap=ef_map)

        pc = ensemble_test.Predictions("Priam")
        pc.generate_priam_predictions(priam_weight, priam_output, logger_handler_level,
                                      definitions.DEFAULT_LOGGER_NAME, efmap=ef_map)
        # pc.read_e2p2_output(priam_weight, priam_e2p2, logger_handler_level, definitions.DEFAULT_LOGGER_NAME,
        #                     efmap=ef_map)

        classifier_dict = {
            'deepec': dc,
            'blast': bc,
            'priam': pc
        }
        for cls_group in name_comb:
            classifiers_name = cls_group
            if len(cls_group) == 0:
                continue
            en = ensemble_test.Ensemble()
            for cls in cls_group:
                en.add_classifier(classifier_dict[cls])
            e2p2_short_output = os.path.join(output_folder, 'rpsd-4.2-20190307.ef-hold.' + '_'.join(classifiers_name) +
                                             definitions.DEFAULT_OUTPUT_SUFFIX)
            logger.log(logging.INFO, "Computing ensemble predictions.")
            # Preform max weight voting on all classifier results
            en.max_weight_voting(logger_handler_level, definitions.DEFAULT_LOGGER_NAME, metacyc_to_ec=metacyc_to_ec)
            # Preform absolute threshold on classifiers voting results
            en.absolute_threshold(0.5, logger_handler_level, definitions.DEFAULT_LOGGER_NAME)
            # Set up object for writing E2P2 output
            print(len(en.final_predictions))
            with open(e2p2_short_output, 'w') as op:
                op.write('Seq\tTrue\tPred\n')
                for seq in sorted(en.final_predictions):
                    # if seq in dc.predictions:
                    op.write(seq + '\t' + '|'.join(sorted(input_dict[seq])) +
                             '\t' + '|'.join(sorted(en.final_predictions[seq])) + '\n')
    except KeyError as classifer_missing:
        logger.log(logging.ERROR, "Missing classifer result file(s): " + str(classifer_missing))
        sys.exit(1)


def iterative_blast_test():
    fasta_folder = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/fasta'
    weight_folder = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/weights/blast'
    blast_folder = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/train_test/blast'
    output_folder = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/e2p2/blast'

    log_path = os.path.join(output_folder, 'log.txt')

    logger_handler_level = "DEBUG"
    cur_logger_config = prog.LoggerConfig()
    cur_logger_config.add_new_logger(definitions.DEFAULT_LOGGER_NAME, log_path)
    logging.config.dictConfig(cur_logger_config.dictConfig)
    logger = logging.getLogger(definitions.DEFAULT_LOGGER_NAME)

    for itr in range(0, 10):
        rounds = str(itr)
        blast_output = 'rpsd-4.2.test_fasta.itr_' + rounds + '.rpsd-4.2.train_fasta.itr_' + rounds + '.blastp.out'
        blast_output = os.path.join(blast_folder, blast_output)
        blast_weight = 'weight.itr_' + rounds
        blast_weight = os.path.join(weight_folder, blast_weight)
        bc = ensemble_test.Predictions("Blast")
        bc.generate_blast_predictions(blast_weight, blast_output, float("1e-2"), logger_handler_level,
                                      definitions.DEFAULT_LOGGER_NAME)

        en = ensemble_test.Ensemble()
        en.add_classifier(bc)
        e2p2_short_output = os.path.join(output_folder,
                                         'rpsd-4.2.test_fasta.itr_' + rounds + '.rpsd-4.2.train_fasta.itr_' + rounds +
                                         definitions.DEFAULT_OUTPUT_SUFFIX)
        logger.log(logging.INFO, "Computing ensemble predictions.")
        # Preform max weight voting on all classifier results
        en.max_weight_voting(logger_handler_level, definitions.DEFAULT_LOGGER_NAME)
        # Preform absolute threshold on classifiers voting results
        en.absolute_threshold(0.5, logger_handler_level, definitions.DEFAULT_LOGGER_NAME)
        # Set up object for writing E2P2 output
        print(len(en.final_predictions))

        with open(e2p2_short_output, 'w') as op:
            for seq in sorted(en.final_predictions):
                # if seq in dc.predictions:
                op.write(seq + '\t' + '|'.join(sorted(en.final_predictions[seq])) + '\n')


def iterative_priam_test():
    fasta_folder = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/fasta'
    weight_folder = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/weights/blast'
    blast_folder = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/train_test/blast'
    output_folder = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/e2p2/blast'

    log_path = os.path.join(output_folder, 'log.txt')

    logger_handler_level = "DEBUG"
    cur_logger_config = prog.LoggerConfig()
    cur_logger_config.add_new_logger(definitions.DEFAULT_LOGGER_NAME, log_path)
    logging.config.dictConfig(cur_logger_config.dictConfig)
    logger = logging.getLogger(definitions.DEFAULT_LOGGER_NAME)

    for itr in range(0, 10):
        rounds = str(itr)
        blast_output = 'rpsd-4.2.test_fasta.itr_' + rounds + '.rpsd-4.2.train_fasta.itr_' + rounds + '.blastp.out'
        blast_output = os.path.join(blast_folder, blast_output)
        blast_weight = 'weight.itr_' + rounds
        blast_weight = os.path.join(weight_folder, blast_weight)
        bc = ensemble_test.Predictions("Blast")
        bc.generate_blast_predictions(blast_weight, blast_output, float("1e-2"), logger_handler_level,
                                      definitions.DEFAULT_LOGGER_NAME)

        en = ensemble_test.Ensemble()
        en.add_classifier(bc)
        e2p2_short_output = os.path.join(output_folder,
                                         'rpsd-4.2.test_fasta.itr_' + rounds + '.rpsd-4.2.train_fasta.itr_' + rounds +
                                         definitions.DEFAULT_OUTPUT_SUFFIX)
        logger.log(logging.INFO, "Computing ensemble predictions.")
        # Preform max weight voting on all classifier results
        en.max_weight_voting(logger_handler_level, definitions.DEFAULT_LOGGER_NAME)
        # Preform absolute threshold on classifiers voting results
        en.absolute_threshold(0.5, logger_handler_level, definitions.DEFAULT_LOGGER_NAME)
        # Set up object for writing E2P2 output
        print(len(en.final_predictions))

        with open(e2p2_short_output, 'w') as op:
            for seq in sorted(en.final_predictions):
                # if seq in dc.predictions:
                op.write(seq + '\t' + '|'.join(sorted(en.final_predictions[seq])) + '\n')


def ara_fscore_filter_test():
    input_fasta_path = '/Users/bxue/Documents/Carnegie/SourceData/TAIR/Araport11_genes.201606.pep.rmspl.fasta'
    output_path = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/' \
                  'ara_fscore_filter_test/Araport11_genes.201606.pep.rmspl.fasta.blast.e2p2v4'
    protein_gene_path = '/Users/bxue/Documents/Carnegie/SourceData/TAIR/Araport11_genes.201606.pep.prot_to_gene.tsv'

    blast_weight = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/weights/avg/blast.avg.txt'
    dummy_efmap = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/weights/avg/' \
                  'blast.avg.dummy_efmaps.txt'
    blast_output = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/iterative_split/ara_fscore_filter_test/' \
                   'Araport11_genes.201606.pep.rmspl.fasta.rpsd-4.2-20190307.blastp.out'
    threshold = 0.5

    time_stamp = str(time.time())

    log_path = os.path.join(os.path.dirname(output_path), 'log.txt')

    logging_level = "DEBUG"
    cur_logger_config = prog.LoggerConfig()
    cur_logger_config.add_new_logger(definitions.DEFAULT_LOGGER_NAME, log_path)
    logging.config.dictConfig(cur_logger_config.dictConfig)
    logger = logging.getLogger(definitions.DEFAULT_LOGGER_NAME)

    # rc.output_dict["blast"] = blast_output
    logger.log(logging.INFO, "Running level-0 classification processes.")

    bc = ensemble_test.Predictions("Blast")
    bc.generate_blast_predictions(blast_weight, blast_output, float("1e-2"), logging_level, definitions.DEFAULT_LOGGER_NAME)

    print(bc.predictions)
    en = ensemble_test.Ensemble()

    en.add_classifier(bc)

    logger.log(logging.INFO, "Computing ensemble predictions.")

    en.max_weight_voting(logging_level, definitions.DEFAULT_LOGGER_NAME)
    en.absolute_threshold(float(threshold), logging_level, definitions.DEFAULT_LOGGER_NAME)

    e2p2_short_output = output_path
    e2p2_long_output = output_path + definitions.DEFAULT_LONG_OUTPUT_SUFFIX
    e2p2_pf_output = output_path + definitions.DEFAULT_PF_OUTPUT_SUFFIX
    e2p2_orxn_pf_output = output_path + definitions.DEFAULT_ORXN_PF_OUTPUT_SUFFIX
    e2p2_final_pf_output = output_path + definitions.DEFAULT_FINAL_PF_OUTPUT_SUFFIX

    e2p2 = file.E2P2files(en.final_predictions, input_fasta_path)
    logger.log(logging.INFO, "Preparing results files.")

    e2p2.add_predictions_of_classifer(bc)

    e2p2.read_efmap(dummy_efmap, logging_level, definitions.DEFAULT_LOGGER_NAME)
    e2p2.write_short_results(definitions.DEFAULT_ENSEMBLE_METHOD + " (" + str(threshold) + ")", e2p2_short_output,
                             logging_level, definitions.DEFAULT_LOGGER_NAME)
    e2p2.write_long_results(definitions.DEFAULT_ENSEMBLE_METHOD + " (" + str(threshold) + ")", e2p2_long_output,
                            logging_level, definitions.DEFAULT_LOGGER_NAME)
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



if __name__ == '__main__':
    # iterative_blast_test()
    main_test()
    # ara_fscore_filter_test()


