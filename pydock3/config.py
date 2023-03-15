import logging
import itertools

import oyaml as yaml
import yamale

from pydock3.files import File
from pydock3.util import get_hexdigest_of_persistent_md5_hash_of_tuple


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Parameter(object):
    def __init__(self, name, value):
        self.name = name
        self.value = value

    @property
    def hexdigest_of_persistent_md5_hash(self):
        return get_hexdigest_of_persistent_md5_hash_of_tuple((self.name, self.value))

    def __bool__(self):
        if self.value:
            return True
        else:
            return False

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return str(self.value)

    def __eq__(self, other):
        if type(other) == type(self):
            return self.value == other.value
        return False


class ParametersConfiguration:
    def __init__(self, config_file_path, schema_file_path):
        #
        File.validate_file_exists(config_file_path)
        File.validate_file_exists(schema_file_path)

        #
        self.config_file_path = config_file_path
        self.schema_file_path = schema_file_path
        self.schema = yamale.make_schema(schema_file_path)

        # validate config based on specified schema
        data = yamale.make_data(self.config_file_path)
        try:
            yamale.validate(self.schema, data)
            logger.debug("Config validation success!")
        except ValueError as e:
            raise Exception("Config validation failed!\n%s" % str(e))

        #
        with open(self.config_file_path, "r") as f:
            self.param_dict = yaml.safe_load(f)

        #
        logger.debug(
            f"ParametersConfiguration instance initialized:\n{self.param_dict}"
        )

    @staticmethod
    def write_config_file(save_path, src_file_path, overwrite=False):
        File.validate_path(src_file_path)
        if File.file_exists(save_path):
            if overwrite:
                logger.info(f"Overwriting existing config file: {save_path}")
            else:
                logger.info(f"A config file already exists: {save_path}")
        else:
            logger.info(f"Writing config file: {save_path}")
        with open(src_file_path, "r") as infile:
            with open(save_path, "w") as outfile:
                yaml.dump(yaml.safe_load(infile), outfile)


def flatten_param_dict(d, key_prefix=""):
    new_d = {}
    for key, value in d.items():
        this_key = f"{key_prefix}{key}"
        if isinstance(value, dict):
            new_d.update(flatten_param_dict(value, f"{this_key}."))
        else:
            new_d[this_key] = value
    return new_d


def flatten_and_parameter_cast_param_dict(d, key_prefix=""):
    new_d = {}
    for key, value in d.items():
        this_key = f"{key_prefix}{key}"
        if isinstance(value, dict):
            new_d.update(flatten_and_parameter_cast_param_dict(value, f"{this_key}."))
        else:
            new_d[this_key] = Parameter(this_key, value)
    return new_d


def sort_list_of_flat_param_dicts(param_dicts):
    param_dict_hashes = []
    for p_dict in param_dicts:
        p_dict_items_interleaved_sorted_by_key_tuple = tuple(
            itertools.chain.from_iterable(
                sorted(list(zip(*list(zip(*p_dict.items())))), key=lambda x: x[0])
            )
        )
        param_dict_hashes.append(
            get_hexdigest_of_persistent_md5_hash_of_tuple(
                p_dict_items_interleaved_sorted_by_key_tuple
            )
        )
    sorted_param_dicts = [
        x
        for x, y in sorted(
            zip(param_dicts, param_dict_hashes),
            key=lambda pair: pair[1],
        )
    ]

    return sorted_param_dicts


def get_sorted_univalued_flat_parameter_cast_param_dicts_from_multivalued_param_dict(multivalued_param_dict):
    #
    keys, multivalues = zip(*sorted(flatten_param_dict(multivalued_param_dict).items(), key=lambda item: item[0]))  # sort by keys

    #
    new_multivalues = []
    for multivalue in multivalues:
        if isinstance(multivalue, list):
            new_multivalues.append(sorted(multivalue))  # is multivalue, sort it
        else:
            new_multivalues.append([multivalue])  # is univalue, so cast as multivalue
    multivalues = new_multivalues

    #
    univalued_flat_parameter_cast_param_dicts = []
    for univalues_combination in itertools.product(*multivalues):
        univalued_flat_parameter_cast_param_dict = {}
        for i, univalue in enumerate(univalues_combination):
            key = keys[i]
            univalued_flat_parameter_cast_param_dict[key] = Parameter(
                name=key, value=univalue
            )
        univalued_flat_parameter_cast_param_dicts.append(
            univalued_flat_parameter_cast_param_dict
        )

    return sort_list_of_flat_param_dicts(univalued_flat_parameter_cast_param_dicts)
