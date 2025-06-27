import os
import logging

from pydock3.config import ParametersConfiguration
from pydock3.decoy_gen import __file__ as DECOY_GEN_INIT_FILE_PATH

#
DECOY_GEN_CONFIG_SCHEMA_FILE_PATH = os.path.join(
    os.path.dirname(DECOY_GEN_INIT_FILE_PATH), "decoy_gen_config_schema.yaml"
)


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class DecoyGenParametersConfiguration(ParametersConfiguration):
    def __init__(self, config_file_path):
        super().__init__(
            config_file_path=config_file_path,
            schema_file_path=DECOY_GEN_CONFIG_SCHEMA_FILE_PATH,
        )