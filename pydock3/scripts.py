import logging
from typing import TYPE_CHECKING, Union, List, Tuple, Dict, Any, Optional

import fire

from pydock3.util import get_logger_for_script, Script


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

SCRIPT_CLASSES = [  # script classes are those that inherit from the Script class defined in this module
    "Blastermaster",
    "Retrodock",
    "Dockopt",
    "DecoyGen",
    "Configure",
    # "TopPoses",
]

NON_SCRIPT_CLASSES_TO_TREAT_AS_SCRIPTS = [  # non-script classes can also be used as scripts through fire; such classes whose functions may be desirable to use as scripts should be included here
    "SDIFile",
]

# Standalone function-based commands
FUNCTION_COMMANDS = [
    "protonate",
    "configure",
]

SCRIPT_CLASSES_DICT = {
    **{cls.lower(): cls for cls in SCRIPT_CLASSES},
    **{cls.lower(): cls for cls in NON_SCRIPT_CLASSES_TO_TREAT_AS_SCRIPTS},
    **{cmd.lower(): cmd for cmd in FUNCTION_COMMANDS},
}

def get_script_class(script_class_name, *args, **kwargs) -> Union[Script, None]:
    logger = get_logger_for_script(debug=True)

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
        elif script_class_name == "decoygen":
            from pydock3.decoy_gen.decoy_gen import DecoyGen as cls
        elif script_class_name == "configure":
            from pydock3.protonation.configure import Configure as cls
        elif script_class_name == "SDIFile":
            from pydock3.files import SDIFile as cls
        elif script_class_name == "protonate":
            return protonate_command
        else:
            raise NotImplementedError

    return cls(*args, **kwargs)

def protonate_command(input_file: str, output_file: str = None, ph: float = 7.4, 
                     score_cutoff: float = 10.0) -> None:
    """
    Generate protomers for SMILES in input file
    
    Args:
        input_file: Path to input file with SMILES (format: "SMILES name" per line)
        output_file: Path to output file (default: input_file + "_protomers.txt")
        ph: pH for protonation (default: 7.4)
        score_cutoff: Minimum combined score (default: 10.0)
    """
    from pydock3.protonation import protonate_smiles_file, ProtonationError
    
    if output_file is None:
        output_file = input_file.rsplit('.', 1)[0] + "_protomers.txt"
    
    try:
        print(f"Generating protomers for {input_file}")
        print(f"Parameters: pH={ph}, score_cutoff={score_cutoff}")
        print(f"Output will be written to: {output_file}")
        
        protonate_smiles_file(input_file, output_file, ph=ph, score_cutoff=score_cutoff)
        print(f"Success! Protomers written to {output_file}")
        
    except ProtonationError as e:
        print(f"Protonation failed: {e}")
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1
    
    return 0


def main():
    fire.Fire(get_script_class)


if __name__ == "__main__":
    main()
