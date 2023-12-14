import random

from src.definitions import DEFAULT_LOGGER_LEVEL, DEFAULT_LOGGER_NAME
from src.lib.classifier import Classifier, FunctionClass
from src.lib.process import logging_helper


class Ensemble(object):
    """Object to ensemble classifier results
    """
    def __init__(self, list_of_classifiers, time_stamp, name="", threshold=float(0), logging_level=DEFAULT_LOGGER_LEVEL,
                 logger_name=DEFAULT_LOGGER_NAME):
        """Placeholder function that retrieve and filter classifier results based on weights.
        Args:
            list_of_classifiers: A list of Classifier objects
            time_stamp: time stamp
            name: name of the ensemble
            logging_level: The logging level set for this ensemble
            logger_name: The name of the logger for this ensemble
        Raises: KeyError
        Returns:
        """
        self.list_of_classifiers = list_of_classifiers
        for idx, classifier in enumerate(list_of_classifiers):
            if not isinstance(classifier, Classifier):
                logging_helper("Exist non-Classifier object", logging_level="ERROR",
                               logger_name=logger_name)
                raise SystemError
            if len(classifier.res) > 0:
                for query in classifier.res:
                    for res_function_class in classifier.res[query]:
                        if not isinstance(res_function_class, FunctionClass):
                            logging_helper("Query\"" + str(query) + "\" results include none FunctionClass",
                                           logging_level="ERROR", logger_name=logger_name)
                            raise SystemError
        self.name = name
        self.threshold = threshold
        self.prediction = Classifier(time_stamp, None, name, logging_level, logger_name)

    def __repr__(self):
        return f'Ensemble(\'{self.name}\', {self.threshold})'

    @staticmethod
    def weighting(list_of_classifiers, weighted_res=None):
        """Placeholder function that retrieve and filter classifier results based on weights.
        Args:
            list_of_classifiers: A list of Classifier objects to be weighted
            weighted_res: A previously available result dictionary
        Raises: KeyError
        Returns:
        """
        if weighted_res is None:
            weighted_res = {}
        for classifier in list_of_classifiers:
            for query in classifier.res:
                try:
                    weighted_res[query] += list(classifier.res[query])
                except KeyError:
                    weighted_res.setdefault(query, list(classifier.res[query]))
        return weighted_res

    @staticmethod
    def voting(voted_res, threshold=float(0)):
        """Placeholder function that preforms voting on classifier results based on weights.
        Args:
            voted_res: A dictionary of results to be voted
            threshold: The threshold for the voting process
        Raises:
        Returns:
        """
        if voted_res is None:
            voted_res = {}
        for query in voted_res:
            voted_function_classes = []
            for function_class in voted_res[query]:
                if function_class.ge_threshold(threshold=threshold, attr="weight"):
                    voted_function_classes.append(function_class)
            voted_res[query] = voted_function_classes
        return voted_res

    @staticmethod
    def ensemble(ensemble_res, queries=None):
        """Function that preforms the weighting and voting process on the list of classifiers.
        Args:
            ensemble_res: A dictionary of results to be ensembled.
            queries: List of input queries
        Raises: KeyError
        Returns:
        """
        if ensemble_res is None:
            ensemble_res = {}
        for query in ensemble_res:
            res_of_query = ensemble_res[query]
            function_names_of_res = list(set([fc.name for fc in res_of_query]))
            for function_name in sorted(function_names_of_res):
                # final result: here we randomly pick one
                function_class_of_ensemble = \
                    random.choice(FunctionClass.get_function_classes_by_vals(res_of_query, function_name))
                try:
                    ensemble_res[query].append(function_class_of_ensemble)
                except KeyError:
                    ensemble_res.setdefault(query, [function_class_of_ensemble])
        if queries is not None:
            for prot in queries:
                if prot not in ensemble_res:
                    ensemble_res.setdefault(prot, [])
        return ensemble_res

    def run(self, threshold=None, previous_res=None, queries=None):
        """Function to retrieve the prediction from this ensemble class.
        Args:
            threshold: threshold for voting
            previous_res: a dictionary that contains previously ran results
            queries: List of input proteins
        Raises: KeyError
        Returns:
        """
        if threshold is None:
            threshold = self.threshold
        weighted_res = self.weighting(self.list_of_classifiers, previous_res)
        voted_res = self.voting(weighted_res, threshold)
        ensemble_res = self.ensemble(voted_res, queries)

        self.prediction.res = ensemble_res

    @staticmethod
    def add_arguments(argument_parser):
        argument_parser.add_argument('Ensemble')

    @staticmethod
    def config_overwrites(args, overwrites=None):
        pass


def run_all_ensembles(list_of_ensemble_names, list_of_ensemble_cls, queries=None, logger_name=DEFAULT_LOGGER_NAME):
    ensembles_ran = []
    skipped_ensembles = []
    for idx, ensemble_cls in enumerate(list_of_ensemble_cls):
        try:
            if isinstance(ensemble_cls, Ensemble):
                logging_helper("Performing Ensemble: %s." % list_of_ensemble_names[idx], logging_level="INFO",
                               logger_name=logger_name)
                ensemble_cls.run(queries=queries)
                ensembles_ran.append(ensemble_cls)
            else:
                skipped_ensembles.append(list_of_ensemble_names[idx])
        except TypeError:
            skipped_ensembles.append(list_of_ensemble_names[idx])
    return ensembles_ran, skipped_ensembles


