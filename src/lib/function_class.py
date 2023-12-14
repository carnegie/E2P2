from functools import total_ordering

_available_class_score_attr = ['weight', 'score']


@total_ordering
class FunctionClass(object):
    """Object for function classes
    """
    def __init__(self, name, score=float(0), weight=float(0)):
        self.name = name
        try:
            self.score = float(score)
        except TypeError:
            self.score = float(0)
        try:
            self.weight = float(weight)
        except TypeError:
            self.weight = float(0)

    def __repr__(self):
        return f'FunctionClass(\'{self.name}\', {self.score}, {self.weight})'

    def retrieve_weight(self, weight_dict):
        """Set the weight for this function class
        Args:
            weight_dict: Dictionary containing weights for function classes
        Raises: KeyError
        Returns:
        """
        try:
            self.weight = weight_dict[self.name]
        except KeyError:
            self.weight = float(0)

    def __eq__(self, function_class):
        """Test if score of this class is equal to function class
        Args:
            function_class: FunctionClass
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if isinstance(function_class, FunctionClass):
            return self.score == function_class.score
        else:
            raise NotImplementedError

    def __lt__(self, function_class):
        """Test if attribute of this class is lesser than function class
        Args:
            function_class: FunctionClass
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if isinstance(function_class, FunctionClass):
            return self.score < function_class.score
        else:
            raise NotImplementedError

    def gt_threshold(self, threshold, attr='score'):
        """Test if attribute of function class is greater than a threshold value
        Args:
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        return getattr(self, attr) > threshold

    def ge_threshold(self, threshold, attr='score'):
        """Test if attribute of function class is greater than or equal to a threshold value
        Args:
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        return getattr(self, attr) >= threshold

    def eq_threshold(self, threshold, attr='score'):
        """Test if attribute of function class is equal to a threshold value
        Args:
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        return getattr(self, attr) == threshold

    def ne_threshold(self, threshold, attr='score'):
        """Test if attribute of function class is not equal to a threshold value
        Args:
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        return getattr(self, attr) != threshold

    def lt_threshold(self, threshold, attr='score'):
        """Test if attribute of function class is lesser than a threshold value
        Args:
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        return getattr(self, attr) < threshold

    def le_threshold(self, threshold, attr='score'):
        """Test if attribute of function class is lesser than or equal to a threshold value
        Args:
            threshold: Threshold for comparison
            attr: Attribute to compare
        Raises: NotImplementedError
        Returns:
            boolean value
        """
        if attr not in _available_class_score_attr:
            attr = 'score'
        return getattr(self, attr) <= threshold

    @staticmethod
    def max_weight(list_of_function_classes):
        """Retrieve class with the maximum weight of a list of FunctionClasses
        Args:
            list_of_function_classes: List of FunctionClasses
        Raises: NotImplementedError
        Returns:
            First function classes that has the maximum weight
        """
        if False not in [isinstance(function_class, FunctionClass) for function_class in list_of_function_classes]:
            top_weight = max(function_class.weight for function_class in list_of_function_classes)
            return [function_class for function_class in list_of_function_classes
                    if function_class.weight == top_weight][0]
        else:
            raise NotImplementedError

    @staticmethod
    def min_weight(list_of_function_classes):
        """Retrieve class with the minimum weight of a list of FunctionClasses
        Args:
            list_of_function_classes: List of FunctionClasses
        Raises: NotImplementedError
        Returns:
            First function classes that has the minimum weight
        """
        if False not in [isinstance(function_class, FunctionClass) for function_class in list_of_function_classes]:
            bottom_weight = min(function_class.weight for function_class in list_of_function_classes)
            return [function_class for function_class in list_of_function_classes
                    if function_class.weight == bottom_weight][0]
        else:
            raise NotImplementedError

    @staticmethod
    def get_function_classes_by_vals(list_of_function_classes, vals, attr='name'):
        """Retrieve classes that has attributes of certain values from a list of FunctionClasses
        Args:
            list_of_function_classes: List of FunctionClasses
            vals: List of values of the attribute
            attr: Attribute to retrieve
        Raises: NotImplementedError
        Returns:
            List of function classes that has the attribute
        """
        if attr not in _available_class_score_attr:
            attr = 'name'
        if type(vals) not in [list, tuple, range, 'generator', set, dict]:
            vals = [vals]
        else:
            vals = list(vals)
        if False not in [isinstance(function_class, FunctionClass) for function_class in list_of_function_classes]:
            return [function_class for function_class in list_of_function_classes
                    if getattr(function_class, attr) in vals]
        else:
            raise NotImplementedError

