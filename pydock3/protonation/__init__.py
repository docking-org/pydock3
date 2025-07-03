"""
PyDock3 Protonation Module

Provides functionality for generating protomers and tautomers using ChemAxon tools,
and filtering molecules based on protonation state matching.
"""

from .core import generate_protomers, ProtonationError, check_chemaxon_tools, protonate_smiles_file
from .config import configure_chemaxon, show_current_config, is_chemaxon_configured

__all__ = ['generate_protomers', 'ProtonationError', 'check_chemaxon_tools', 
           'protonate_smiles_file', 'configure_chemaxon', 'show_current_config', 'is_chemaxon_configured']