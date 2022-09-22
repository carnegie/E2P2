import random

from src.lib.ensemble import Ensemble
from src.lib.classifier import FunctionClass


class MaxWeightAbsoluteThreshold(Ensemble):
    @staticmethod
    def weighting(list_of_classifiers, input_res=None):
        output_res = {}
        if input_res is None:
            input_res = {}
        for classifier in list_of_classifiers:
            for query in classifier.res:
                query_classifier_res = list(classifier.res[query])
                try:
                    query_all_res = list(input_res[query]) + query_classifier_res
                except KeyError:
                    query_all_res = query_classifier_res
                try:
                    output_res[query] += query_all_res
                except KeyError:
                    output_res.setdefault(query, query_all_res)
        for query in output_res:
            query_pred = []
            query_res = output_res[query]
            function_names_of_query = set(FunctionClass.get_function_classes_attr_vals(query_res))
            for function_name in sorted(function_names_of_query):
                function_classes_of_name = FunctionClass.get_function_classes_by_vals(query_res, [function_name])
                max_weight = FunctionClass.max_attr(function_classes_of_name, "weight")
                function_classes_of_max_weight = \
                    FunctionClass.get_function_classes_by_vals(function_classes_of_name, [max_weight], attr="weight")
                query_pred.append(random.choice(function_classes_of_max_weight))
            output_res[query] = query_pred
        return output_res

    @staticmethod
    def voting(input_res, threshold=float(0.5)):
        output_res = {}
        if input_res is None:
            input_res = {}
        for query in input_res:
            voted_res = []
            query_res = input_res[query]
            if len(query_res) == 0:
                output_res.setdefault(query, voted_res)
                continue
            max_weight = FunctionClass.max_attr(query_res, attr='weight')
            t = float(max_weight - threshold)
            # Check to make sure the threshold is not negative.
            if t < 0.0:
                t = 0.0
            for res_function_class in query_res:
                if FunctionClass.greater_and_equal_threshold(res_function_class, t, attr='weight'):
                    voted_res.append(res_function_class)
            output_res.setdefault(query, voted_res)
        return output_res




