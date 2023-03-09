from typing import TypeVar, Callable
from typing_extensions import ParamSpec
from dataclasses import fields
import subprocess
import logging
import os
import sys
import traceback
import hashlib
import inspect


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


def get_logger_for_script(log_file_path, debug=False):
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
            logger.info(f"Shutting down job due to exception: {exc_type}")
            logger.info(traceback.format_exc())
            logger.debug(str(exc_value))
            return True
        return exc_type is None


def filter_kwargs_for_callable(d: dict, callable: Callable[P, T]) -> dict:
    return {k: v for k, v in d.items() if k in inspect.getfullargspec(callable).args}
