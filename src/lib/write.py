import re
from datetime import datetime

from src.definitions import DEFAULT_LOGGER_LEVEL, DEFAULT_LOGGER_NAME, DEFAULT_LONG_OUTPUT_SUFFIX, \
    DEFAULT_PF_OUTPUT_SUFFIX, DEFAULT_ORXN_PF_OUTPUT_SUFFIX, DEFAULT_FINAL_PF_OUTPUT_SUFFIX
from src.lib.classifier import Classifier, FunctionClass
from src.lib.ensemble import Ensemble
from src.lib.process import logging_helper
from src.lib.read import read_e2p2_maps


class PfFiles(object):
    def __init__(self, cls_to_write, input_proteins=None, logger_name=DEFAULT_LOGGER_NAME):
        """Initialize class
        Args:
            cls_to_write: Input can be an Ensemble class or a prediction dictionary
            input_proteins: List of input protein IDs
            logger_name: The name of the logger
        Raises: SystemError
        Returns:
        """
        if isinstance(cls_to_write, Ensemble):
            self.final_prediction = cls_to_write.prediction.res
        elif isinstance(cls_to_write, Classifier):
            self.final_prediction = cls_to_write.res
        elif type(cls_to_write) is dict:
            self.final_prediction = cls_to_write
        else:
            logging_helper("PfFiles initialization failure", logging_level="ERROR", logger_name=logger_name)
            raise SystemError
        if input_proteins is not None:
            for prot in input_proteins:
                if prot not in self.final_prediction:
                    self.final_prediction.setdefault(prot, [])

    def write_short_results(self, ensemble_name, output_path, logging_level=DEFAULT_LOGGER_LEVEL,
                            logger_name=DEFAULT_LOGGER_NAME):
        """Write E2P2 short version of result to output
        Args:
            ensemble_name: Name of the ensemble method
            output_path: Path to output for short version of result
            logging_level: The logging level set for write short results
            logger_name: The name of the logger for write short results
        Raises: AttributeError, NotImplementedError
        Returns:
        """
        cur_time = datetime.now()
        header = "# Result Generation time:  %s\n# Ensemble method used:  %s\n" % (cur_time, ensemble_name)
        with open(output_path, 'w') as op:
            op.write(header)
            try:
                for query in sorted(self.final_prediction.keys()):
                    predictions = self.final_prediction[query]
                    if len(predictions) == 0:
                        op.write('\t'.join([query, 'NA']) + '\n')
                    else:
                        predicted_classes = list(set([fc.name for fc in predictions]))
                        op.write('\t'.join([query, '|'.join(predicted_classes)]) + '\n')
            except (AttributeError, NotImplementedError) as e:
                logging_helper(
                    "Error when writing results: " + str(e), logging_level="ERROR", logger_name=logger_name)
        logging_helper(
            "Results written to: \"" + output_path + "\"", logging_level=logging_level, logger_name=logger_name)

    def write_long_results(self, list_of_classifiers, ensemble_name, output_path,
                           logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        """Write E2P2 long version of result to output
        Args:
            list_of_classifiers: List of classifiers used in ensemble
            ensemble_name: Name of the ensemble method
            output_path: Path to output for short version of result
            logging_level: The logging level set for write long results
            logger_name: The name of the logger for write long results
        Raises: AttributeError, NotImplementedError
        Returns:
        """
        cur_time = datetime.now()
        header = "# Result Generation time:  %s\n# Ensemble method used:  %s\n" % (cur_time, ensemble_name)
        with open(output_path, 'w') as op:
            op.write(header)
            try:
                for query in sorted(self.final_prediction.keys()):
                    predictions = self.final_prediction[query]
                    if len(predictions) == 0:
                        op.write('\t'.join(['>' + query, 'NA']) + '\n')
                    else:
                        predicted_classes = list(set([fc.name for fc in predictions]))
                        op.write('\t'.join(['>' + query, '|'.join(predicted_classes)]) + '\n')
                        for classifier in list_of_classifiers:
                            if isinstance(classifier, Classifier):
                                classifier_name = classifier.name
                                classifier_res = classifier.res
                                output_list = []
                                try:
                                    classifier_classes = classifier_res[query]
                                    for function_cls in classifier_classes:
                                        if isinstance(function_cls, FunctionClass):
                                            output_list.append('|'.join(
                                                [function_cls.name, str(function_cls.weight), str(function_cls.score)]))
                                except KeyError:
                                    continue
                                op.write('\t'.join([classifier_name + ':'] + output_list) + '\n')
            except (AttributeError, NotImplementedError) as e:
                logging_helper(
                    "Error when writing results: " + str(e), logging_level="ERROR", logger_name=logger_name)
        logging_helper(
            "Results written to: \"" + output_path + "\"", logging_level=logging_level, logger_name=logger_name)

    def write_pf_results(self, ef_map_path, output_path, logging_level=DEFAULT_LOGGER_LEVEL,
                         logger_name=DEFAULT_LOGGER_NAME):
        """Write pf file of result
        Args:
            ef_map_path: Path to EF to EC/RXN mapping file
            output_path: Path to output for short version of result
            logging_level: The logging level set for write pf results
            logger_name: The name of the logger for write pf results
        Raises: AttributeError, KeyError
        Returns:
        """
        ef_map_dict = read_e2p2_maps(ef_map_path, 0, 1)
        with open(output_path, 'w') as op:
            try:
                for query in sorted(self.final_prediction.keys()):
                    predictions = self.final_prediction[query]
                    if len(predictions) == 0:
                        continue
                    else:
                        op.write("ID\t%s\nNAME\t%s\nPRODUCT-TYPE\tP\n" % (query, query))
                        predicted_classes = list(set([fc.name for fc in predictions]))
                        for ef_class in sorted(predicted_classes):
                            try:
                                mapped_ids = ["METACYC\t" + i if "RXN" in i else "EC\t" + i for i in
                                              ef_map_dict[ef_class]]
                                op.write('\n'.join(mapped_ids) + '\n')
                            except KeyError:
                                logging_helper(
                                    "EF Class: \"" + ef_class + "\" assigned to \"" + query + "\" not found in map.",
                                    logging_level="ERROR", logger_name=logger_name)
                        op.write("//\n")
            except (AttributeError, NotImplementedError) as e:
                logging_helper(
                    "Error when writing results: " + str(e), logging_level="ERROR", logger_name=logger_name)
        logging_helper(
            "Results written to: \"" + output_path + "\"", logging_level=logging_level, logger_name=logger_name)

    @staticmethod
    def map_ec_to_rxns(list_of_ecs, ec_superseded_dict, metacyc_rxn_ec_dict, official_ec_metacyc_rxn_dict,
                       to_remove_metabolsim_list):
        """Write orxn file of result
        Args:
            list_of_ecs: EF to EC/MetaCyc RXN mapping file
            ec_superseded_dict: EC superseded mapping
            metacyc_rxn_ec_dict: MetaCyc RXN to EC mapping (EC -> RXN)
            official_ec_metacyc_rxn_dict: Official EC to MetaCyc RXN mapping
            to_remove_metabolsim_list: List of non-small molecule metabolisms
        Raises:
        Returns:
            List of official MetaCyc RXN IDs
            List of unofficial MetaCyc RXN IDs
        """
        metacyc_official = set()
        metacyc_unofficial = set()
        ec_ids_updated = set()
        for ec in sorted(list_of_ecs):
            try:
                ec_superseded_ids = [e.replace('EC-', '') for e in ec_superseded_dict["EC-" + ec]]
                ec_ids_updated.update(ec_superseded_ids)
            except KeyError:
                ec_ids_updated.add(ec)
        for ec_updated in sorted(ec_ids_updated):
            if ec_updated in official_ec_metacyc_rxn_dict:
                metacyc_official.update(official_ec_metacyc_rxn_dict[ec_updated])
            elif ec_updated in metacyc_rxn_ec_dict:
                metacyc_unofficial.update(metacyc_rxn_ec_dict[ec_updated])
        return \
            sorted(set([rxn for rxn in sorted(metacyc_official)
                        if rxn not in to_remove_metabolsim_list])), \
            sorted(set([rxn for rxn in sorted(metacyc_unofficial)
                        if rxn not in to_remove_metabolsim_list]))

    def write_orxn_results(self, ef_map_path, ec_superseded_path, metacyc_rxn_ec_path, official_ec_metacyc_rxn_path,
                           to_remove_metabolism_path, output_path, prot_gene_map_path=None,
                           logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
        """Write orxn file of result
        Args:
            ef_map_path: Path to EF to EC/RXN mapping file
            ec_superseded_path: Path to EC superseded mapping
            metacyc_rxn_ec_path: Path to MetaCyc RXN to EC mapping
            official_ec_metacyc_rxn_path: Path to official EC to MetaCyc RXN mapping
            to_remove_metabolism_path: Path to list of non-small molecule metabolisms
            output_path: Path to output for short version of result
            prot_gene_map_path: Path to protein to gene ID mapping
            logging_level: The logging level set for write orxn results
            logger_name: The name of the logger for write orxn results
        Raises: AttributeError, KeyError
        Returns:
        """
        ef_map_dict = read_e2p2_maps(ef_map_path, 0, 1)
        ec_superseded_dict = read_e2p2_maps(ec_superseded_path, 2, 0)
        metacyc_rxn_ec_dict = read_e2p2_maps(metacyc_rxn_ec_path, 1, 0)
        official_ec_metacyc_rxn_dict = read_e2p2_maps(official_ec_metacyc_rxn_path, 0, 1)
        to_remove_metabolism_list = sorted(read_e2p2_maps(to_remove_metabolism_path, 0, 0).keys())
        if prot_gene_map_path is not None:
            prot_gene_map_dict = read_e2p2_maps(prot_gene_map_path, 0, 1)
        else:
            prot_gene_map_dict = {}
        with open(output_path, 'w') as op:
            try:
                for query in sorted(self.final_prediction.keys()):
                    metacyc_ids = set()
                    metacyc_unofficial = set()
                    predictions = self.final_prediction[query]
                    if len(predictions) == 0:
                        continue
                    else:
                        predicted_classes = list(set([fc.name for fc in predictions]))
                        for ef_class in sorted(predicted_classes):
                            try:
                                metacyc_ids.update([i for i in ef_map_dict[ef_class] if "RXN" in i and
                                                    i not in to_remove_metabolism_list])
                                ec_ids = [i for i in ef_map_dict[ef_class] if "RXN" not in i and
                                          i not in to_remove_metabolism_list]
                                metacyc_from_ecs = self.map_ec_to_rxns(
                                    ec_ids, ec_superseded_dict, metacyc_rxn_ec_dict, official_ec_metacyc_rxn_dict,
                                    to_remove_metabolism_list)
                                metacyc_ids.update(metacyc_from_ecs[0])
                                metacyc_unofficial.update(metacyc_from_ecs[1])
                            except KeyError:
                                logging_helper(
                                    "EF Class: \"" + ef_class + "\" assigned to \"" + query + "\" not found in map.",
                                    logging_level="ERROR", logger_name=logger_name)
                    if len(metacyc_ids) > 0 or len(metacyc_unofficial) > 0:
                        try:
                            gene_id = prot_gene_map_dict[query][0]
                            op.write("ID\t%s\nNAME\t%s\nPRODUCT-ACCESSION\t%s\nPRODUCT-TYPE\tP\n" %
                                     (gene_id, gene_id, query))
                        except KeyError:
                            if len(prot_gene_map_dict) == 0:
                                op.write("ID\t%s\nNAME\t%s\nPRODUCT-TYPE\tP\n" % (query, query))
                            else:
                                logging_helper("No Gene found for protein: " + query, logging_level="ERROR",
                                               logger_name=logger_name)
                                raise SystemError
                        if len(metacyc_ids) > 0:
                            op.write('\n'.join(['METACYC\t' + m for m in sorted(metacyc_ids)]) + '\n')
                        if len(metacyc_unofficial) > 0:
                            op.write('\n'.join(['METACYC\t' + m + '\n#unofficial'
                                                for m in sorted(metacyc_unofficial)]) + '\n')
                        op.write("//\n")
            except (AttributeError, NotImplementedError) as e:
                logging_helper(
                    "Error when writing results: " + str(e), logging_level="ERROR", logger_name=logger_name)
        logging_helper(
            "Results written to: \"" + output_path + "\"", logging_level=logging_level, logger_name=logger_name)


def write_ensemble_outputs(ensemble_cls, all_query_ids, output_path, ef_map_path, ec_superseded_path,
                           metacyc_rxn_ec_path, official_ec_metacyc_rxn_path, to_remove_metabolism_path,
                           prot_gene_map_path=None, logging_level=DEFAULT_LOGGER_LEVEL,
                           logger_name=DEFAULT_LOGGER_NAME):
    """Write all ensemble results
    Args:
        ensemble_cls: class of the ensemble that was used
        all_query_ids: list of all query IDs
        output_path: output path to the short output file
        ef_map_path: path to efclasses.mapping
        ec_superseded_path: path to EC-superseded.mapping
        metacyc_rxn_ec_path: path to metacyc-RXN-EC.mapping
        official_ec_metacyc_rxn_path: path to official-EC-metacyc-RXN.mapping
        to_remove_metabolism_path: path to to-remove-non-small-molecule-metabolism.mapping
        prot_gene_map_path: Path to protein to gene ID mapping
        logging_level: The logging level set for write orxn results
        logger_name: The name of the logger for write orxn results
    Raises:
    Returns:
    """
    logging_helper("Writing outputs...", logging_level=logging_level, logger_name=logger_name)
    ensemble_name = re.sub(r'[^\w\-_\. ]', '_', ensemble_cls.prediction.name)
    ensemble_classifiers = ensemble_cls.list_of_classifiers

    ensemble_output = PfFiles(ensemble_cls,all_query_ids)
    short_output_path = '.'.join([output_path, ensemble_name])
    ensemble_output.write_short_results(ensemble_name, output_path, logging_level="INFO", logger_name=logger_name)

    long_output_path = '.'.join([output_path, ensemble_name, DEFAULT_LONG_OUTPUT_SUFFIX])
    ensemble_output.write_long_results(ensemble_classifiers, ensemble_name, long_output_path, logging_level="INFO",
                                       logger_name=logger_name)

    pf_output_path = '.'.join([output_path, ensemble_name, DEFAULT_PF_OUTPUT_SUFFIX])
    ensemble_output.write_pf_results(ef_map_path, pf_output_path, logging_level="INFO", logger_name=logger_name)

    orxn_output_path = '.'.join([output_path, ensemble_name, DEFAULT_ORXN_PF_OUTPUT_SUFFIX])
    ensemble_output.write_orxn_results(ef_map_path, ec_superseded_path, metacyc_rxn_ec_path,
                                       official_ec_metacyc_rxn_path, to_remove_metabolism_path, orxn_output_path,
                                       prot_gene_map_path=None, logging_level="INFO",
                                       logger_name=logger_name)
    if prot_gene_map_path is not None:
        final_output_path = '.'.join([output_path, ensemble_name, DEFAULT_FINAL_PF_OUTPUT_SUFFIX])
        ensemble_output.write_orxn_results(ef_map_path, ec_superseded_path, metacyc_rxn_ec_path,
                                           official_ec_metacyc_rxn_path, to_remove_metabolism_path, final_output_path,
                                           prot_gene_map_path=prot_gene_map_path, logging_level="INFO",
                                           logger_name=logger_name)

