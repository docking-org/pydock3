import os
import itertools
import logging
from copy import deepcopy

from pydock3.config import Parameter, ParametersConfiguration
from pydock3.dockopt import __file__ as DOCKOPT_INIT_FILE_PATH

#
DOCKOPT_CONFIG_SCHEMA_FILE_PATH = os.path.join(os.path.dirname(DOCKOPT_INIT_FILE_PATH),
                                                  "dockopt_config_schema.yaml")


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class DockoptParametersConfiguration(ParametersConfiguration):

    def __init__(self, config_file_path):
        super().__init__(config_file_path=config_file_path, schema_file_path=DOCKOPT_CONFIG_SCHEMA_FILE_PATH)

        # if value type is list then we treat key as optimizable (i.e. it is a list of possible values to try)
        self.optimizable_parameter_keys = [key for key, param in self.param_dict.items() if isinstance(param.value, list)]

        self.job_param_dicts = self.get_param_dicts_for_different_optimization_configurations(self.param_dict)

    def get_param_dicts_for_different_optimization_configurations(self, opt_param_dict):
        param_dicts = []
        opt_param_values_list = [opt_param_dict[opt_param_key].value for opt_param_key in self.optimizable_parameter_keys]
        for opt_param_value_combination_tuple in itertools.product(*opt_param_values_list):
            param_dict = deepcopy(opt_param_dict)
            for i, opt_param_value in enumerate(opt_param_value_combination_tuple):
                opt_param_key = self.optimizable_parameter_keys[i]
                param_dict[opt_param_key] = Parameter(name=opt_param_key, value=opt_param_value)
            param_dicts.append(param_dict)

        return param_dicts
