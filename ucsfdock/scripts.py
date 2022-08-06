import logging

import fire

from ucsfdock.util import get_logger_for_script
from ucsfdock.blastermaster.blastermaster import Blastermaster
from ucsfdock.dockmaster.dockmaster import Dockmaster
from ucsfdock.files import SDIFile

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


SCRIPT_CLASSES = [
    Blastermaster,
    Dockmaster,
    SDIFile,
]
SCRIPT_CLASSES_DICT = {cls.__name__.lower(): cls for cls in SCRIPT_CLASSES}


def get_script_class(script_class_name, *args, **kwargs):
    logger = get_logger_for_script("ucsfdock.log", debug=False)

    if script_class_name not in SCRIPT_CLASSES_DICT:
        logger.error(f"script_class_name must be one of:\n{sorted(list(SCRIPT_CLASSES_DICT.keys()))}")
        return

    return SCRIPT_CLASSES_DICT[script_class_name](*args, **kwargs)


def main():
    fire.Fire(get_script_class)


if __name__ == '__main__':
    main()
