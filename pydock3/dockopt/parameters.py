from typing import List, Tuple, Union
from functools import reduce
from operator import getitem
from copy import deepcopy

import pandas as pd

from pydock3.config import flatten_param_dict
from pydock3.jobs import DOCK3_EXECUTABLE_PATH


class ParametersManager(object):
    def __init__(self, parameters_dict):
        self._parameters_dict = parameters_dict

    @property
    def parameters_dict(self):
        return self._parameters_dict

    @property
    def flattened_parameters_dict(self):
        return flatten_param_dict(self._parameters_dict)


class DockoptComponentParametersManager(ParametersManager):
    def __init__(self, parameters_dict, last_component_completed=None):
        #
        if last_component_completed is not None:
            for row_index, row in last_component_completed.load_results_dataframe().head(last_component_completed.top_n).iterrows():
                nested_target_keys_and_value_tuples = self._load_nested_target_keys_and_value_tuples_from_dataframe_row(row, identifier_prefix='parameters.', include_prefix=True)
                for nested_target_keys, value in nested_target_keys_and_value_tuples:
                    parameters_dict = self._get_parameters_dict_with_next_step_reference_value_replaced(parameters_dict, nested_target_keys, new_ref=value, old_ref='^')

        #
        parameters_dict = self._get_parameters_dict_with_next_step_numerical_operators_applied(parameters_dict)

        #
        if isinstance(parameters_dict["parameters"]["dock_executable_path"], list):
            new_dock_executable_path_value = []
            for dock_executable_path in parameters_dict["parameters"]["dock_executable_path"]:
                if dock_executable_path is None:
                    new_dock_executable_path_value.append(DOCK3_EXECUTABLE_PATH)
                else:
                    new_dock_executable_path_value.append(dock_executable_path)
        else:
            if parameters_dict["parameters"]["dock_executable_path"] is None:
                new_dock_executable_path_value = DOCK3_EXECUTABLE_PATH
            else:
                new_dock_executable_path_value = parameters_dict["parameters"]["dock_executable_path"]
        parameters_dict["parameters"]["dock_executable_path"] = new_dock_executable_path_value

        #
        super().__init__(parameters_dict)

    @staticmethod
    def _get_parameters_dict_with_next_step_reference_value_replaced(parameters_dict: dict, nested_target_keys: List[str], new_ref: float, old_ref: str = '^') -> dict:
        """Takes a set of parameters, finds the next nested step to be run, and, if it
        contains numerical operators, replaces the `reference_value` of the `target_key`
        with the specified float `new_ref` if `reference_value` matches the string
        `old_ref`."""

        def get_nested_dict_item(dic, nested_keys):
            """Get item in nested dictionary"""
            return reduce(getitem, nested_keys, dic)

        def set_nested_dict_item(dic, nested_keys, value):
            """Set item in nested dictionary"""
            reduce(getitem, nested_keys[:-1], dic)[nested_keys[-1]] = value
            return dic

        def traverse(obj):
            if isinstance(obj, dict):
                try:
                    nested_target = get_nested_dict_item(obj, nested_target_keys)
                except KeyError:
                    for key, value in obj.items():
                        obj[key] = traverse(value)
                    return obj
                if isinstance(nested_target, dict):
                    if 'reference_value' in nested_target and 'arguments' in nested_target and 'operator' in nested_target:  # numerical operator detected
                        # replace old ref with new ref
                        if nested_target['reference_value'] == old_ref:
                            obj = set_nested_dict_item(obj, nested_target_keys + ['reference_value'], new_ref)
                else:
                    if nested_target == old_ref:
                        obj = set_nested_dict_item(obj, nested_target_keys, new_ref)
                return obj
            elif isinstance(obj, list):  # obj is sequence
                obj[0] = traverse(obj[0])  # only change next step to be run, which will be found in the first element
                return obj
            else:
                return obj

        return traverse(deepcopy(parameters_dict))

    @staticmethod
    def _load_nested_target_keys_and_value_tuples_from_dataframe_row(row: pd.Series, identifier_prefix: str = 'parameters.', include_prefix: bool = False) -> List[Tuple[List[str], Union[float, str]]]:
        """Loads the parameters in a dataframe row according to the column names."""

        dic = row.to_dict()
        nested_target_keys_and_value_tuples = [(key.split('.'), value) for key, value in dic.items() if key.startswith(identifier_prefix)]

        if not include_prefix:
            nested_target_keys_and_value_tuples = [(x[0][1:], x[1]) for x in nested_target_keys_and_value_tuples]

        return nested_target_keys_and_value_tuples

    @staticmethod
    def _get_parameters_dict_with_next_step_numerical_operators_applied(parameters_dict: dict) -> dict:
        """Takes a set of parameters, finds the next nested step to be run, and, if it
        contains numerical operators, applies them."""

        def traverse(obj):
            if isinstance(obj, dict):
                if 'reference_value' in obj and 'arguments' in obj and 'operator' in obj:  # numerical operator detected
                    # apply operators
                    if obj['operator'] == '+':
                        obj = [float(obj['reference_value']) + float(x) for x in obj['arguments']]
                    elif obj['operator'] == '-':
                        obj = [float(obj['reference_value']) - float(x) for x in obj['arguments']]
                    elif obj['operator'] == '*':
                        obj = [float(obj['reference_value']) * float(x) for x in obj['arguments']]
                    elif obj['operator'] == '/':
                        obj = [float(obj['reference_value']) / float(x) for x in obj['arguments']]
                    else:
                        raise ValueError(
                            f"Witnessed operator `{obj['operator']}`. Only the following numerical operators are supported: `+`, `-`, `*`, `/`")
                else:
                    for key, value in obj.items():
                        obj[key] = traverse(value)
                return obj
            elif isinstance(obj, list):  # obj is sequence
                obj[0] = traverse(obj[0])  # only change next step to be run, which will be found in the first element
                return obj
            else:
                return obj

        return traverse(deepcopy(parameters_dict))
