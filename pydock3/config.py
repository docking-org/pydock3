import logging

import pandas as pd
import yamale
import oyaml as yaml

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
        self.config_file_path = config_file_path

        #
        self.schema_file_path = schema_file_path
        self.schema = yamale.make_schema(schema_file_path)

        #
        data = yamale.make_data(self.config_file_path)
        try:
            yamale.validate(self.schema, data)
            logger.debug('Config validation success!')
        except ValueError as e:
            raise Exception('Config validation failed!\n%s' % str(e))

        #
        param_dict, = pd.json_normalize(data[0][0]).to_dict('records')  # TODO: add validation
        self.param_dict = {key: Parameter(name=key, value=value) for key, value in param_dict.items()}

        #
        logger.debug(f"Parameters:\n{self.param_dict}")

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
        with open(src_file_path, 'r') as infile:
            with open(save_path, "w") as outfile:
                yaml.dump(yaml.safe_load(infile), outfile)
