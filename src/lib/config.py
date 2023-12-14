import configparser
import os
from argparse import ArgumentParser

from src.definitions import ROOT_DIR, DEFAULT_LOGGER_NAME, CONFIG_PATH, CLASSIFIERS_CLS_DIR, WEIGHTS_DIR, ENSEMBLES_CLS_DIR, MAPS_DIR
from src.lib.process import logging_helper, load_module_function_from_path


_DEFAULT_SECTIONS = ["Mapping", "Ensembles", "Classifiers"]


def get_options_from_config_section(config, section_name, name_only=True, selected_options=None,
                                    logging_level="WARNING", logger_name=DEFAULT_LOGGER_NAME):
    """Helper function to retrieve options from a configparser section
    Args:
        config: configparser
        section_name: Name of the sections to retrieve
        name_only: Boolean value for
        selected_options: specific options to retrieve
        logging_level: The logging level set for read map
        logger_name: The name of the logger for read map
    Raises:
    Returns:
        option names or options with their items
    """
    if selected_options is None:
        selected_options = set()
    elif type(selected_options) in [list, tuple, range, 'generator', set, dict]:
        selected_options = set(selected_options)
    else:
        return None
    try:
        try:
            section_options = config.items(section_name)
        except configparser.InterpolationMissingOptionError:
            section_options = config.items(section_name, raw=True)
        section_options_names = [tup[0] for tup in section_options]
        if len(selected_options) == 0:
            selected_options = set(section_options_names)
        if len(selected_options & set(section_options_names)) == 0:
            return None
        elif name_only is True:
            return sorted(selected_options & set(section_options_names))
        else:
            return {tup[0]: tup[1] for tup in section_options if tup[0] in selected_options}
    except configparser.NoSectionError:
        logging_helper("Missing '" + str(section_name) + "' section in config",
                       logging_level=logging_level, logger_name=logger_name)
        return None


def get_values_from_config_option(config, section_name, option_name, logging_level="WARNING",
                                  logger_name=DEFAULT_LOGGER_NAME):
    """Helper function to retrieve values from configparser options
    Args:
        config: configparser
        section_name: Name of the section to retrieve
        option_name: Name of the option to retrieve
        logging_level: The logging level set for read map
        logger_name: The name of the logger for read map
    Raises:
    Returns:
        option_val: values of the option
    """
    try:
        try:
            option_val = config.get(section_name, option_name)
        except configparser.InterpolationMissingOptionError:
            option_val = config.get(section_name, option_name, raw=True)
        if option_val == '':
            return None
        else:
            return option_val
    except configparser.NoSectionError:
        logging_helper("Missing '" + str(section_name) + "' section in config",
                       logging_level=logging_level, logger_name=logger_name)
        return None
    except configparser.NoOptionError:
        logging_helper("Missing '" + str(option_name) + "' in '" + str(section_name) + "' section",
                       logging_level=logging_level, logger_name=logger_name)
        return None


def config_section_to_multi_sections_helper(config, section_name, selected_options=None,
                                            logging_level="WARNING", logger_name=DEFAULT_LOGGER_NAME):
    """Helper function for a section that points to multiple sections, i.e. Ensembles & Classifiers
    Args:
        config: configparser
        section_name: Name of the section to retrieve
        selected_options: Name of the options to retrieve
        logging_level: The logging level set for read map
        logger_name: The name of the logger for read map
    Raises:
    Returns:
        sections_names: names of the referenced multiple sections
        sections_dict: a dictionary that contains options and their values
    """
    potential_sections = get_options_from_config_section(config, section_name, logging_level=logging_level,
                                                         logger_name=logger_name)
    if potential_sections is not None:
        available_sections_names = [get_values_from_config_option(config, section_name, potential_section_name,
                                                                  logging_level=logging_level, logger_name=logger_name)
                                    for potential_section_name in potential_sections]
    else:
        available_sections_names = None
    if available_sections_names is not None:
        available_sections_dicts = [get_options_from_config_section(config, section, selected_options=selected_options,
                                                                    name_only=False)
                                    for section in available_sections_names]
    else:
        available_sections_dicts = None
    return available_sections_names, available_sections_dicts


def default_path_helper(input_path, default_folder, root_dir=ROOT_DIR):
    """Get full paths using defaults
    Args:
        input_path: input file path
        default_folder: Path to default files
        root_dir: Project root
    Raises:
    Returns:
        The full path for a default file
    """
    if os.path.join(root_dir, os.path.dirname(input_path)) == default_folder:
        return os.path.join(root_dir, input_path)
    else:
        return input_path


def read_config_ini(timestamp, config_ini, io_dict, overwrites=None, logging_level="WARNING",
                    logger_name=DEFAULT_LOGGER_NAME):
    """Read in config.ini
    Args:
        timestamp:
        config_ini: configparser
        io_dict: Dictionary that represents an "IO" section for the configparser
        overwrites: A dictionary to overwrite values of the config.ini
        logging_level: The logging level set for read map
        logger_name: The name of the logger for read map
    Raises:
    Returns:
        sections_names: names of the referenced multiple sections
        sections_dict: a dictionary that contains options and their values
    """
    pipeline_config = configparser.ConfigParser(allow_no_value=True, interpolation=configparser.ExtendedInterpolation())
    pipeline_config.read_dict(io_dict)
    pipeline_config.read(config_ini)

    if overwrites is not None:
        pipeline_config.read_dict(overwrites)

    # IO
    query_path = \
        get_values_from_config_option(pipeline_config, 'IO', 'query', logging_level=logging_level, logger_name=logger_name)
    temp_path = \
        get_values_from_config_option(pipeline_config, 'IO', 'out', logging_level=logging_level, logger_name=logger_name)

    # Mapping Files
    efclasses = get_values_from_config_option(pipeline_config, 'Mapping', 'efclasses', logging_level=logging_level,
                                              logger_name=logger_name)
    efclasses = default_path_helper(efclasses, MAPS_DIR)
    ec_superseded = get_values_from_config_option(pipeline_config, 'Mapping', 'ec_superseded', logging_level=logging_level,
                                                  logger_name=logger_name)
    ec_superseded = default_path_helper(ec_superseded, MAPS_DIR)
    metacyc_rxn_ec = get_values_from_config_option(pipeline_config, 'Mapping', 'metacyc_rxn_ec', logging_level=logging_level,
                                                   logger_name=logger_name)
    metacyc_rxn_ec = default_path_helper(metacyc_rxn_ec, MAPS_DIR)
    official_ec_metacyc_rxn = get_values_from_config_option(pipeline_config, 'Mapping', 'official_ec_metacyc_rxn',
                                                            logging_level=logging_level, logger_name=logger_name)
    official_ec_metacyc_rxn = default_path_helper(official_ec_metacyc_rxn, MAPS_DIR)
    to_remove_non_small_molecule_metabolism = \
        get_values_from_config_option(pipeline_config, 'Mapping', 'to_remove_non_small_molecule_metabolism',
                                      logging_level=logging_level, logger_name=logger_name)
    to_remove_non_small_molecule_metabolism = default_path_helper(to_remove_non_small_molecule_metabolism, MAPS_DIR)
    mapping_files = {
        'efclasses': efclasses, 'ec_superseded': ec_superseded, 'metacyc_rxn_ec': metacyc_rxn_ec,
        "official_ec_metacyc_rxn": official_ec_metacyc_rxn,
        "to_remove_non_small_molecule_metabolism": to_remove_non_small_molecule_metabolism}
    # Classifiers
    classifier_sections, classifier_tuples = \
        config_section_to_multi_sections_helper(pipeline_config, "Classifiers",
                                                logging_level=logging_level, logger_name=logger_name)
    list_of_classifiers = []
    if classifier_sections is not None and classifier_tuples is not None:
        for idx, cls in enumerate(classifier_sections):
            if cls is not None:
                try:
                    cls_module_path = overwrites[cls]["class"]
                except (KeyError, TypeError) as e:
                    cls_module_path = get_values_from_config_option(pipeline_config, cls, "class")
                    if cls_module_path is not None:
                        cls_module_path = default_path_helper(cls_module_path, CLASSIFIERS_CLS_DIR)
                cls_fn = load_module_function_from_path(cls_module_path, cls)
                if cls_fn is not None:
                    cls_config_dict = classifier_tuples[idx]
                    if cls_config_dict is None:
                        list_of_classifiers.append(None)
                        continue
                    try:
                        weight_path = overwrites[cls]["weight"]
                    except (KeyError, TypeError) as e:
                        weight_path = cls_config_dict["weight"]
                        if weight_path is not None:
                            weight_path = default_path_helper(weight_path, WEIGHTS_DIR)
                    cls_class = cls_fn(timestamp, weight_path, cls)
                    cls_class.setup_classifier(query_path, temp_path, cls_config_dict, cls)
                    list_of_classifiers.append(cls_class)
                else:
                    list_of_classifiers.append(None)
            else:
                list_of_classifiers.append(None)

    # Ensembles
    ensemble_sections, ensemble_tuples = \
        config_section_to_multi_sections_helper(pipeline_config, "Ensembles",
                                                logging_level=logging_level, logger_name=logger_name)
    list_of_ensembles = []
    if ensemble_sections is not None and ensemble_tuples is not None:
        for idx, ens in enumerate(ensemble_sections):
            if ens is not None:
                try:
                    ens_module_path = overwrites[ens]["weight"]
                except (KeyError, TypeError) as e:
                    ens_module_path = get_values_from_config_option(pipeline_config, ens, "class")
                    if ens_module_path is not None:
                        ens_module_path = default_path_helper(ens_module_path, ENSEMBLES_CLS_DIR)
                ens_fn = load_module_function_from_path(ens_module_path, ens)
                list_of_ensembles.append(ens_fn)
            else:
                list_of_ensembles.append(None)
    return classifier_sections, list_of_classifiers, ensemble_sections, list_of_ensembles, mapping_files


def read_config(config_ini, io_dict=None, overwrites=None, logging_level="DEBUG", logger_name=DEFAULT_LOGGER_NAME):
    logging_helper("Processing config.ini", logging_level, logger_name)
    pipeline_config = configparser.ConfigParser(allow_no_value=True, interpolation=configparser.ExtendedInterpolation())
    if io_dict is not None and type(io_dict) is dict:
        pipeline_config.read_dict(io_dict)
    pipeline_config.read(config_ini)
    if overwrites is not None and type(overwrites) is dict:
        pipeline_config.read_dict(overwrites)

    mapping_dict = {}
    mappings_options = get_options_from_config_section(pipeline_config, "Mapping")
    if mappings_options is None:
        logging_helper("No [Mapping] section in config.ini", "ERROR", logger_name)
        raise SystemError
    for map_option in mappings_options:
        map_value = get_values_from_config_option(pipeline_config, "Mapping", map_option)
        mapping_dict.setdefault(map_option, map_value)

    classifier_dict = {}
    classifier_sections = config_section_to_multi_sections_helper(pipeline_config, "Classifiers")
    if classifier_sections is None:
        logging_helper("No [Classifiers] section in config.ini", "ERROR", logger_name)
        raise SystemError
    num_of_classifiers = len(classifier_sections[0])
    for idx in range(num_of_classifiers):
        classifier_dict.setdefault(classifier_sections[0][idx], classifier_sections[1][idx])

    ensemble_dict = {}
    ensemble_sections = config_section_to_multi_sections_helper(pipeline_config, "Ensembles")
    if ensemble_sections is None:
        logging_helper("No [Ensembles] section in config.ini", "ERROR", logger_name)
        raise SystemError
    num_of_ensembles = len(ensemble_sections[0])
    for idx in range(num_of_ensembles):
        ensemble_dict.setdefault(ensemble_sections[0][idx], ensemble_sections[1][idx])

    return mapping_dict, classifier_dict, ensemble_dict


if __name__ == '__main__':
    parser = ArgumentParser()
    config_path = CONFIG_PATH

    ensemble_dict = read_config(config_path, None)[2]
    classifier_dict = read_config(config_path, None)[1]
    print(classifier_dict)
    for cls in classifier_dict:
        cls_path = os.path.join(ROOT_DIR, classifier_dict[cls]["class"])
        cls_fn = load_module_function_from_path(cls_path, cls)
        cls_fn.add_arguments(parser)
    # parser.print_help()
    overwrites = {}
    args = parser.parse_args()
    for cls in classifier_dict:
        cls_path = os.path.join(ROOT_DIR, classifier_dict[cls]["class"])
        cls_fn = load_module_function_from_path(cls_path, cls)
        cls_fn.config_overwrites(args, overwrites)

    print("overwrite", overwrites)
    io = {"IO": {"query": "INPUT", "blast": "BLAST", "priam": "PRIAM", "timestamp": "000"}}
    classifier_dict = read_config(config_path, io, overwrites)[1]
    print(classifier_dict)

