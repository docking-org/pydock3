import os
import logging

from ucsfdock.util import ParametersConfiguration
from ucsfdock.blastermaster import __file__ as BLASTERMASTER_INIT_FILE_PATH
BLASTERMASTER_CONFIG_SCHEMA_FILE_PATH = os.path.join(os.path.dirname(BLASTERMASTER_INIT_FILE_PATH), "blastermaster_config_schema.yaml")


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class BlastermasterParametersConfiguration(ParametersConfiguration):

    def __init__(self, config_file_path):
        super().__init__(config_file_path=config_file_path, schema_file_path=BLASTERMASTER_CONFIG_SCHEMA_FILE_PATH)
