import logging

import fire

from pydock3.util import get_logger_for_script
from pydock3.blastermaster.blastermaster import Blastermaster
from pydock3.retrodock.retrodock import Retrodock
from pydock3.dockopt.dockopt import Dockopt
from pydock3.top_poses import TopPoses
from pydock3.files import SDIFile

#
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


SCRIPT_CLASSES = [  # script classes are those that inherit from the Script class defined in this module
    Blastermaster,
    # Retrodock,
    Dockopt,
    # TopPoses,
]
NON_SCRIPT_CLASSES_TO_TREAT_AS_SCRIPTS = [  # non-script classes can also be used as scripts through fire; such classes whose functions may be desirable to use as scripts should be included here
    # SDIFile,
]

SCRIPT_CLASSES_DICT = {
    **{cls.__name__.lower(): cls for cls in SCRIPT_CLASSES},
    **{cls.__name__.lower(): cls for cls in NON_SCRIPT_CLASSES_TO_TREAT_AS_SCRIPTS},
}


def get_script_class(script_class_name, *args, **kwargs):
    logger = get_logger_for_script(debug=False)

    if script_class_name not in SCRIPT_CLASSES_DICT:
        logger.error(
            f"script_class_name must be one of:\n{sorted(list(SCRIPT_CLASSES_DICT.keys()))}"
        )
        return

    return SCRIPT_CLASSES_DICT[script_class_name](*args, **kwargs)


def main():
    fire.Fire(get_script_class)


if __name__ == "__main__":
    main()
