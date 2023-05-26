import logging
from typing import TYPE_CHECKING, Union, List, Tuple, Dict, Any, Optional

import fire

from pydock3.util import get_logger_for_script

if TYPE_CHECKING:
    from pydock3.util import Script


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

SCRIPT_CLASSES = [  # script classes are those that inherit from the Script class defined in this module
    "Blastermaster",
    "Retrodock",
    "Dockopt",
    # "TopPoses",
]

NON_SCRIPT_CLASSES_TO_TREAT_AS_SCRIPTS = [  # non-script classes can also be used as scripts through fire; such classes whose functions may be desirable to use as scripts should be included here
    "SDIFile",
]

SCRIPT_CLASSES_DICT = {
    **{cls.__name__.lower(): cls for cls in SCRIPT_CLASSES},
    **{cls.__name__.lower(): cls for cls in NON_SCRIPT_CLASSES_TO_TREAT_AS_SCRIPTS},
}

def get_script_class(script_class_name, *args, **kwargs) -> Union[Script, None]:
    logger = get_logger_for_script(debug=False)

    if script_class_name not in SCRIPT_CLASSES_DICT:
        logger.error(
            f"script_class_name must be one of:\n{sorted(list(SCRIPT_CLASSES_DICT.keys()))}"
        )
        return

    if script_class_name in SCRIPT_CLASSES_DICT:
        if script_class_name == "blastermaster":
            from pydock3.blastermaster.blastermaster import Blastermaster as cls
        elif script_class_name == "retrodock":
            from pydock3.retrodock.retrodock import Retrodock as cls
        elif script_class_name == "dockopt":
            from pydock3.dockopt.dockopt import Dockopt as cls
        elif script_class_name == "SDIFile":
            from pydock3.files import SDIFile as cls
        else:
            raise NotImplementedError

    return cls(*args, **kwargs)


def main():
    fire.Fire(get_script_class)


if __name__ == "__main__":
    main()
