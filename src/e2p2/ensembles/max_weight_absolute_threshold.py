import random

from src.lib.ensemble import Ensemble
from src.lib.classifier import FunctionClass


class MaxWeightAbsoluteThreshold(Ensemble):
    @staticmethod
    def weighting(list_of_classifiers, weighted_res=None):
        if weighted_res is None:
            weighted_res = {}
        for classifier in list_of_classifiers:
            for query in classifier.res:
                query_res = list(classifier.res[query])
                try:
                    weighted_res[query] = list(weighted_res[query]) + query_res
                except KeyError:
                    weighted_res.setdefault(query, query_res)
        for query in weighted_res:
            query_res = list(weighted_res[query])
            res_function_names = list(set([fc.name for fc in query_res]))

            weighted_query_res = []
            for function_name in sorted(res_function_names):
                function_classes_w_name = (
                    FunctionClass.get_function_classes_by_vals(query_res, function_name))
                max_weight_w_name = FunctionClass.max_weight(function_classes_w_name).weight
                weighted_query_res.append(
                    random.choice(FunctionClass.get_function_classes_by_vals(function_classes_w_name, max_weight_w_name,
                                                                             "weight")))
            weighted_res[query] = weighted_query_res
        return weighted_res

    @staticmethod
    def voting(voted_res, threshold=float(0.5)):
        if voted_res is None:
            voted_res = {}
        for query in voted_res:
            query_res = list(voted_res[query])
            if len(query_res) == 0:
                continue
            max_weight = FunctionClass.max_weight(query_res).weight
            t = float(max_weight) - float(threshold)
            # Check to make sure the threshold is not negative.
            if t < 0.0:
                t = 0.0
            voted_query_res = []
            for res_function_class in query_res:
                if res_function_class.ge_threshold(t, attr='weight'):
                    voted_query_res.append(res_function_class)
            voted_res[query] = voted_query_res
        return voted_res

    @staticmethod
    def add_arguments(argument_parser):
        """Function to add E2P2 ensemble related arguments
        Args:
            argument_parser: argparse
        Raises:
        Returns:
        """
        # Arguments for E2P2 ensembles
        argument_parser.add_argument("--threshold", "-t", dest="threshold", type=float,
                                     help="Threshold for voting results. Default is 0.5.")

    @staticmethod
    def config_overwrites(args, overwrites=None):
        max_weight_abs_threshold_dest = ["threshold"]
        if overwrites is None:
            overwrites = {}
        max_weight_abs_threshold_overwrites = {}
        args_dict = vars(args)
        for dest in max_weight_abs_threshold_dest:
            try:
                val = args_dict[dest]
                if val is not None:
                    key = dest
                    max_weight_abs_threshold_overwrites.setdefault(key, val)
            except KeyError:
                continue
        if len(max_weight_abs_threshold_overwrites) > 0:
            overwrites.setdefault("MaxWeightAbsoluteThreshold", {})
            overwrites["MaxWeightAbsoluteThreshold"] = max_weight_abs_threshold_overwrites
        return overwrites
