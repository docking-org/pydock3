import logging

import oyaml as yaml
import yamale

from pydock3.files import File


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Parameter(object):
    def __init__(self, name, value):
        self.name = name
        self.value = value

    def __bool__(self):
        if self.value:
            return True
        else:
            return False

    def __str__(self):
        return str(self.value)

    def __eq__(self, other):
        if type(other) == type(self):
            return self.value == other.value
        return False

    def __hash__(self):
        return hash((self.name, self.value))


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
