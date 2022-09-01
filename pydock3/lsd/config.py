import os
import logging

from pydock3.config import ParametersConfiguration
from pydock3.lsd import __file__ as LSD_INIT_FILE_PATH

#
LSD_CONFIG_SCHEMA_FILE_PATH = os.path.join(os.path.dirname(LSD_INIT_FILE_PATH),
                                               "lsd_config_schema.yaml")


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class LSDParametersConfiguration(ParametersConfiguration):

    def __init__(self, config_file_path):
        super().__init__(config_file_path=config_file_path, schema_file_path=LSD_CONFIG_SCHEMA_FILE_PATH)
