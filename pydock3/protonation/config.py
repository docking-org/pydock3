"""
Configuration management for PyDock3 protonation module
"""

import os
import json
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any

# Save config with the PyDock installation, not in user home directory
PYDOCK_DIR = Path(__file__).parent.parent  # pydock3/ directory
CONFIG_FILE = PYDOCK_DIR / "protonation_config.json"

def get_config() -> Dict[str, Any]:
    """Load protonation configuration from file"""
    if CONFIG_FILE.exists():
        try:
            with open(CONFIG_FILE, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    return {}

def save_config(config: Dict[str, Any]) -> None:
    """Save protonation configuration to file"""
    try:
        with open(CONFIG_FILE, 'w') as f:
            json.dump(config, f, indent=2)
    except IOError as e:
        raise RuntimeError(f"Failed to save configuration: {e}")

def get_chemaxon_paths() -> tuple[Optional[str], Optional[str]]:
    """
    Get paths to cxcalc and molconvert executables
    
    Returns:
        (cxcalc_path, molconvert_path) or (None, None) if not configured
    """
    config = get_config()
    jchem_path = config.get('jchem_path')
    
    if not jchem_path:
        return None, None
    
    jchem_bin = Path(jchem_path) / "bin"
    cxcalc_path = jchem_bin / "cxcalc"
    molconvert_path = jchem_bin / "molconvert"
    
    # Check if files exist and are executable
    if cxcalc_path.exists() and molconvert_path.exists():
        return str(cxcalc_path), str(molconvert_path)
    
    return None, None

def setup_chemaxon_environment() -> bool:
    """
    Setup ChemAxon environment variables
    
    Returns:
        True if setup successful, False otherwise
    """
    config = get_config()
    license_path = config.get('license_path')
    
    if license_path:
        os.environ['CHEMAXON_LICENSE_URL'] = license_path
        return True
    
    return False

def validate_chemaxon_installation(jchem_path: str, license_path: str) -> tuple[bool, str]:
    """
    Validate ChemAxon installation and license
    
    Args:
        jchem_path: Path to JChem installation directory
        license_path: Path to license file or license URL
        
    Returns:
        (is_valid, error_message)
    """
    jchem_path = Path(jchem_path)
    
    # Check if JChem directory exists
    if not jchem_path.exists():
        return False, f"JChem directory does not exist: {jchem_path}"
    
    # Check for executables
    cxcalc_path = jchem_path / "bin" / "cxcalc"
    molconvert_path = jchem_path / "bin" / "molconvert"
    
    if not cxcalc_path.exists():
        return False, f"cxcalc not found at: {cxcalc_path}"
    
    if not molconvert_path.exists():
        return False, f"molconvert not found at: {molconvert_path}"
    
    # Check if executables are runnable by testing with license
    old_license = os.environ.get('CHEMAXON_LICENSE_URL')
    try:
        os.environ['CHEMAXON_LICENSE_URL'] = license_path
        
        # Test cxcalc
        result = subprocess.run([str(cxcalc_path), "--help"], 
                              capture_output=True, timeout=10)
        if result.returncode != 0:
            return False, f"cxcalc failed to run (check license): {result.stderr.decode()}"
        
        # Test molconvert (use -h for older versions)
        result = subprocess.run([str(molconvert_path), "-h"], 
                              capture_output=True, timeout=10)
        if result.returncode != 0:
            return False, f"molconvert failed to run (check license): {result.stderr.decode()}"
        
        return True, "ChemAxon tools validated successfully"
        
    except subprocess.TimeoutExpired:
        return False, "ChemAxon tools timed out (possible license issue)"
    except FileNotFoundError:
        return False, "ChemAxon executables not found or not executable"
    except Exception as e:
        return False, f"Validation failed: {e}"
    finally:
        # Restore original license
        if old_license:
            os.environ['CHEMAXON_LICENSE_URL'] = old_license
        elif 'CHEMAXON_LICENSE_URL' in os.environ:
            del os.environ['CHEMAXON_LICENSE_URL']

def configure_chemaxon(jchem_path: str, license_path: str) -> bool:
    """
    Configure ChemAxon tools for PyDock3
    
    Args:
        jchem_path: Path to JChem installation directory
        license_path: Path to license file or license URL
        
    Returns:
        True if configuration successful, False otherwise
    """
    print("Validating ChemAxon installation...")
    
    is_valid, message = validate_chemaxon_installation(jchem_path, license_path)
    
    if not is_valid:
        print(f"❌ Validation failed: {message}")
        return False
    
    print(f"✅ {message}")
    
    # Save configuration
    config = {
        'jchem_path': str(Path(jchem_path).absolute()),
        'license_path': license_path
    }
    
    try:
        save_config(config)
        print(f"✅ Configuration saved to {CONFIG_FILE}")
        return True
    except Exception as e:
        print(f"❌ Failed to save configuration: {e}")
        return False

def show_current_config() -> None:
    """Display current configuration"""
    config = get_config()
    
    if not config:
        print("❌ No ChemAxon configuration found")
        print(f"Run 'python -m pydock3.scripts configure' to set up ChemAxon tools")
        return
    
    print("📋 Current ChemAxon Configuration:")
    print(f"   JChem Path: {config.get('jchem_path', 'Not set')}")
    print(f"   License Path: {config.get('license_path', 'Not set')}")
    
    # Test if tools are working
    cxcalc_path, molconvert_path = get_chemaxon_paths()
    if cxcalc_path and molconvert_path:
        setup_chemaxon_environment()
        try:
            result1 = subprocess.run([cxcalc_path, "--help"], 
                                   capture_output=True, timeout=5)
            result2 = subprocess.run([molconvert_path, "-h"], 
                                   capture_output=True, timeout=5)
            if result1.returncode == 0 and result2.returncode == 0:
                print("   Status: ✅ Tools are working")
            else:
                print("   Status: ❌ Tools not responding (check license)")
        except Exception:
            print("   Status: ❌ Tools not accessible")
    else:
        print("   Status: ❌ Tools not found")

def is_chemaxon_configured() -> bool:
    """Check if ChemAxon tools are properly configured"""
    cxcalc_path, molconvert_path = get_chemaxon_paths()
    license_configured = setup_chemaxon_environment()
    return cxcalc_path is not None and molconvert_path is not None and license_configured