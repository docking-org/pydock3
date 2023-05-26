from typing import TypeVar, Callable, Iterable, Hashable, Any
from typing_extensions import ParamSpec
from dataclasses import fields
import subprocess
import logging
import os
import sys
import traceback
import hashlib
import inspect
from functools import reduce
from operator import getitem


#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#
logging_formatter = logging.Formatter(
    "%(asctime)s;%(levelname)s;%(message)s", "%Y-%m-%d %H:%M:%S"
)

#
T = TypeVar('T')  # generic type for type hinting
P = ParamSpec('P')  # generic parameters spec for type hinting


def validate_variable_type(var, allowed_instance_types):
    if not isinstance(var, allowed_instance_types):
        raise Exception(
            f"Variable '{var}' must be an instance of one of allowed_instance_types={allowed_instance_types}. Type witnessed: {type(var)}"
        )


def get_hexdigest_of_persistent_md5_hash_of_tuple(t: tuple) -> str:
    m = hashlib.md5()
    for s in t:
        m.update(str(s).encode())
    return m.hexdigest()


def get_dataclass_as_dict(data_cls):
    return {field.name: getattr(data_cls, field.name) for field in fields(data_cls)}


def system_call(command_str, cwd=os.getcwd(), timeout_seconds=None, env_vars_dict=None):
    logger.debug(
        f"Running system call.\nCurrent working directory: {cwd}\n Command:\n{command_str}"
    )
    proc = subprocess.run(
        command_str,
        cwd=cwd,
        shell=True,
        text=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout_seconds,
        env=env_vars_dict,
    )
    logger.debug(
        f"System call returned: {proc}\n\nstdout:{proc.stdout}\n\nstderr:{proc.stderr}\n"
    )

    return proc


class Script(object):
    """Base class for all classes intended to serve as scripts in the package."""

    def __init__(self):
        pass


def get_logger_for_script(log_file_path=None, debug=False):
    # get highest-level logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    #
    sh = logging.StreamHandler(sys.stdout)
    if debug:
        sh.setLevel(logging.DEBUG)
    else:
        sh.setLevel(logging.INFO)
    sh.setFormatter(logging_formatter)
    logger.addHandler(sh)

    #
    if log_file_path is not None:
        fh = logging.FileHandler(log_file_path)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging_formatter)
        logger.addHandler(fh)

    return logger


class CleanExit(object):
    def __init__(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        if exc_type is not None:
            logger.exception(f"Shutting down job due to exception: {exc_type}")
            #logger.exception(traceback.format_exc())
            logger.debug(str(exc_value))
            return True
        return exc_type is None


def filter_kwargs_for_callable(d: dict, callable: Callable[P, T]) -> dict:
    return {k: v for k, v in d.items() if k in inspect.getfullargspec(callable).args}


def get_nested_dict_item(dic: dict, nested_keys: Iterable[Hashable]) -> Any:
    """Get item in nested dictionary"""

    return reduce(getitem, nested_keys, dic)


def set_nested_dict_item(dic: dict, nested_keys: Iterable[Hashable], value: Any) -> dict:
    """Set item in nested dictionary"""

    reduce(getitem, nested_keys[:-1], dic)[nested_keys[-1]] = value
    return dic


def find_key_values_in_dict(nested_dict, key):
    result = []

    def traverse_dict(d):
        for k, v in d.items():
            if k == key:
                result.append(v)
            elif isinstance(v, dict):
                traverse_dict(v)

    traverse_dict(nested_dict)
    return result


def get_ordinal(n: int) -> str:
    """Get ordinal number (e.g. 1st, 2nd, 3rd, 4th, etc.)"""
    return "%d%s" % (
        n,
        "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10:: 4],
    )


def sort_list_by_another_list(list_to_be_sorted: list, list_to_sort_by: list) -> list:
    """Sort one list by the sort order of another list"""
    return tuple(zip(*sorted(zip(list_to_be_sorted, list_to_sort_by), key=lambda x: x[1])))[0]
