import random

from src.definitions import DEFAULT_LOGGER_LEVEL, DEFAULT_LOGGER_NAME
from src.lib.classifier import Classifier, FunctionClass
from src.lib.util import logging_helper


class Ensemble(object):
    """Object to ensemble classifier results
    """
    def __init__(self, list_of_classifiers, time_stamp, name, logging_level=DEFAULT_LOGGER_LEVEL,
                 logger_name=DEFAULT_LOGGER_NAME):
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
        self.prediction = Classifier(time_stamp, None, name, logging_level, logger_name)

    @staticmethod
    def weighting(list_of_classifiers, input_res=None):
        """Placeholder function that retrieve and filter classifier results based on weights.
        Args:
            list_of_classifiers: A list of Classifier objects
            input_res: A previously available result dictionary
        Raises: KeyError
        Returns:
        """
        if input_res is None:
            input_res = {}
        for classifier in list_of_classifiers:
            for query in classifier.res:
                try:
                    input_res[query] += list(classifier.res[query])
                except KeyError:
                    input_res.setdefault(query, list(classifier.res[query]))
        return input_res

    @staticmethod
    def voting(input_res, threshold=float(0)):
        """Placeholder function that preforms voting on classifier results based on weights.
        Args:
            input_res: A previously available result dictionary
            threshold: The threshold for the voting process
        Raises:
        Returns:
        """
        if input_res is None:
            input_res = {}
        for query in input_res:
            voted_in_pred = []
            for input_function_class in input_res[query]:
                if FunctionClass.greater_and_equal_threshold(input_function_class, threshold=threshold, attr="weight"):
                    voted_in_pred.append(input_function_class)
            input_res[query] = voted_in_pred
        return input_res

    def ensemble(self, proteins=None):
        """Function that preforms the weighting and voting process on the list of classifiers.
        Args:
            proteins: List of input proteins
        Raises: KeyError
        Returns:
        """
        pred = {}
        weighted_res = self.weighting(self.list_of_classifiers)
        voted_res = self.voting(weighted_res)
        for query in voted_res:
            voted_res_of_query = voted_res[query]
            function_names_of_voted_res = set(FunctionClass.get_function_classes_attr_vals(voted_res_of_query))
            for function_name in sorted(function_names_of_voted_res):
                # final result currently we just use 1 from all
                function_class_of_query = \
                    random.choice(FunctionClass.get_function_classes_by_vals(voted_res_of_query, function_name))
                try:
                    pred[query].append(function_class_of_query)
                except KeyError:
                    pred.setdefault(query, [function_class_of_query])
        if proteins is not None:
            for prot in proteins:
                if prot not in pred:
                    pred.setdefault(prot, [])
        setattr(self.prediction, "res", pred)

    def get_prediction(self, proteins=None):
        res = {}
        for query in getattr(self.prediction, "res"):
            res.setdefault(query, getattr(self.prediction, "res")[query])
        if proteins is not None:
            for prot in proteins:
                if prot not in getattr(self.prediction, "res"):
                    res.setdefault(prot, [])
        return res


def run_all_ensembles(list_of_ensemble_names, list_of_ensemble_fns, list_of_classifiers, time_stamp,
                      logging_level=DEFAULT_LOGGER_LEVEL, logger_name=DEFAULT_LOGGER_NAME):
    ensembles_ran = []
    skipped_ensembles = []
    for idx, ensemble_fn in enumerate(list_of_ensemble_fns):
        try:
            ensemble_cls = ensemble_fn(list_of_classifiers, time_stamp, list_of_ensemble_names[idx],
                                       logging_level, logger_name)
            if isinstance(ensemble_cls, Ensemble):
                logging_helper("Performing Ensemble: %s." % list_of_ensemble_names[idx], logging_level="INFO",
                               logger_name=logger_name)
                ensemble_cls.ensemble()
                ensembles_ran.append(ensemble_cls)
            else:
                skipped_ensembles.append(list_of_ensemble_names[idx])
        except TypeError:
            skipped_ensembles.append(list_of_ensemble_names[idx])
    return ensembles_ran, skipped_ensembles


