"""
Configuration class for PyDock3 protonation tools
"""

import logging
from typing import Optional

from pydock3.util import Script

logger = logging.getLogger(__name__)

class Configure(Script):
    """Configuration management for PyDock3 protonation tools"""
    
    def __init__(self):
        super().__init__()
    
    def run(self) -> None:
        """
        Run all configuration steps for PyDock3
        
        This will configure all external software dependencies needed for PyDock3.
        Currently includes: ChemAxon tools for protonation filtering.
        """
        print("🔧 PyDock3 Complete Configuration")
        print("=" * 50)
        print()
        print("This will configure all external software for PyDock3.")
        print()
        
        # List of all configuration steps
        config_steps = [
            ("ChemAxon Tools", "Required for protonation filtering", self._configure_chemaxon_interactive),
            # Future: Add more software here
            # ("OpenEye Tools", "Required for 3D generation", self._configure_openeye_interactive),
            # ("Schrodinger Tools", "Required for ligand preparation", self._configure_schrodinger_interactive),
        ]
        
        print("Configuration steps:")
        for i, (name, description, _) in enumerate(config_steps, 1):
            print(f"  {i}. {name} - {description}")
        print()
        
        # Run each configuration step
        for i, (name, description, config_func) in enumerate(config_steps, 1):
            print(f"Step {i}/{len(config_steps)}: Configuring {name}")
            print("-" * 40)
            
            try:
                config_func()
                print(f"✅ {name} configuration completed")
            except KeyboardInterrupt:
                print(f"\n❌ {name} configuration cancelled by user")
                print("You can configure individual components later with:")
                print(f"  python -m pydock3.scripts configure chemaxon")
                return
            except Exception as e:
                print(f"❌ {name} configuration failed: {e}")
                print("You can try configuring individual components later with:")
                print(f"  python -m pydock3.scripts configure chemaxon")
            
            print()
        
        # Final status check
        print("🎯 Configuration Summary")
        print("=" * 30)
        self.show()
        
        print("\n🎉 PyDock3 configuration completed!")
        print("You can now use all configured features.")
    
    def _configure_chemaxon_interactive(self) -> None:
        """Interactive ChemAxon configuration"""
        from pydock3.protonation.config import configure_chemaxon
        
        print("ChemAxon tools are required for protonation filtering in decoy generation.")
        print()
        
        jchem_path = input("Enter path to JChem installation directory: ").strip()
        if not jchem_path:
            raise ValueError("JChem path is required")
        
        license_path = input("Enter path to ChemAxon license file or license URL: ").strip()
        if not license_path:
            raise ValueError("License path is required")
        
        print("\nValidating ChemAxon installation...")
        success = configure_chemaxon(jchem_path, license_path)
        
        if not success:
            raise RuntimeError("ChemAxon validation failed")
    
    def chemaxon(self, jchem_path: Optional[str] = None, license_path: Optional[str] = None) -> None:
        """
        Configure ChemAxon tools for PyDock3 protonation
        
        Args:
            jchem_path: Path to JChem installation directory
            license_path: Path to license file or license URL
        """
        from pydock3.protonation.config import configure_chemaxon
        
        if not jchem_path or not license_path:
            print("🔧 PyDock3 ChemAxon Configuration")
            print("=" * 40)
            print()
            
            if not jchem_path:
                jchem_path = input("Enter path to JChem installation directory: ").strip()
            
            if not license_path:
                license_path = input("Enter path to ChemAxon license file or license URL: ").strip()
            
            if not jchem_path or not license_path:
                print("❌ Both JChem path and license path are required")
                return
        
        print()
        print("Configuring ChemAxon tools...")
        
        success = configure_chemaxon(jchem_path, license_path)
        
        if success:
            print()
            print("🎉 ChemAxon configuration completed successfully!")
            print("You can now use protonation filtering in PyDock3.")
        else:
            print()
            print("❌ ChemAxon configuration failed.")
            print("Please check the paths and license file.")
    
    def show(self) -> None:
        """Show current ChemAxon configuration"""
        from pydock3.protonation.config import show_current_config
        show_current_config()
    
    def test(self) -> None:
        """Test ChemAxon tools configuration"""
        from pydock3.protonation import check_chemaxon_tools, is_chemaxon_configured
        
        print("🔧 Testing ChemAxon Configuration")
        print("=" * 40)
        
        configured = is_chemaxon_configured()
        tools_work = check_chemaxon_tools()
        
        print(f"Configuration status: {'✅ Configured' if configured else '❌ Not configured'}")
        print(f"Tools working: {'✅ Working' if tools_work else '❌ Not working'}")
        
        if configured and tools_work:
            print("\n🎉 ChemAxon tools are ready to use!")
        else:
            print(f"\n❌ ChemAxon tools are not ready.")
            if not configured:
                print("Run: python -m pydock3.scripts configure chemaxon")
            else:
                print("Check your ChemAxon installation and license.")