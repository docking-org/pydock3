import os
import logging

from pydock3.config import Parameter, ParametersConfiguration
from pydock3.dockopt import __file__ as DOCKOPT_INIT_FILE_PATH


#
DOCKOPT_CONFIG_SCHEMA_FILE_PATH = os.path.join(
    os.path.dirname(DOCKOPT_INIT_FILE_PATH), "dockopt_config_schema.yaml"
)


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class DockoptParametersConfiguration(ParametersConfiguration):
    def __init__(self, config_file_path):
        super().__init__(
            config_file_path=config_file_path,
            schema_file_path=DOCKOPT_CONFIG_SCHEMA_FILE_PATH,
        )
