"""
Dynamic Environment Manager for Roshambo Shape Similarity Components.
This module handles dynamic configuration of environment variables and conda environments
for flexible Roshambo deployment across different user setups.
"""

import os
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, Optional, List, Tuple
import logging


class RoshamboDynamicEnvironmentManager:
    """
    Manages dynamic environment configuration for Roshambo components.
    Handles different user setups, conda installations, and RDKit build locations.
    """

    def __init__(self,
                 rdbase_path: str,
                 conda_env_name: str,
                 conda_base_path: Optional[str] = None,
                 cuda_home_path: str = "/usr/local/cuda",
                 auto_setup_env: bool = True,
                 debug: bool = False):
        """
        Initialize the dynamic environment manager.

        Args:
            rdbase_path: REQUIRED - Path to RDKit build directory
            conda_env_name: REQUIRED - Name of the conda environment with Roshambo
            conda_base_path: OPTIONAL - Custom path to conda installation (auto-discovered if not provided)
            cuda_home_path: OPTIONAL - Path to CUDA installation
            auto_setup_env: Whether to automatically set up environment variables
            debug: Enable debug logging
        """
        if not rdbase_path:
            raise ValueError("rdbase_path is required - specify the path to your RDKit build directory")
        if not conda_env_name:
            raise ValueError("conda_env_name is required - specify the name of your conda environment with Roshambo")

        self.rdbase_path = rdbase_path
        self.conda_env_name = conda_env_name
        self.conda_base_path = conda_base_path
        self.cuda_home_path = cuda_home_path
        self.auto_setup_env = auto_setup_env
        self.debug = debug

        # Set up logging
        self.logger = logging.getLogger(__name__)
        if debug:
            logging.basicConfig(level=logging.DEBUG)

        # Discovered paths
        self._conda_executable = None
        self._conda_env_path = None
        self._validated_rdbase = None

        # Environment validation results
        self._validation_results = {}

    def discover_conda_installation(self) -> Tuple[Optional[str], Optional[str]]:
        """
        Discover conda installation and environment paths.

        Returns:
            Tuple of (conda_executable_path, conda_env_path)
        """
        conda_executable = None
        conda_env_path = None

        # Try different conda locations
        conda_candidates = [
            "conda",  # In PATH
            "/opt/conda/bin/conda",
            "/home/*/miniconda3/bin/conda",
            "/home/*/anaconda3/bin/conda",
        ]

        if self.conda_base_path:
            conda_candidates.insert(0, os.path.join(self.conda_base_path, "bin", "conda"))

        for candidate in conda_candidates:
            try:
                if "*" in candidate:
                    # Handle wildcard paths
                    import glob
                    matches = glob.glob(candidate)
                    if matches:
                        candidate = matches[0]

                result = subprocess.run([candidate, "--version"],
                                      capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    conda_executable = candidate
                    if self.debug:
                        self.logger.debug(f"Found conda at: {conda_executable}")
                    break
            except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
                continue

        if conda_executable:
            # Get conda environment path
            try:
                result = subprocess.run([conda_executable, "env", "list"],
                                      capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    for line in result.stdout.split('\n'):
                        if self.conda_env_name in line and '*' not in line:
                            parts = line.split()
                            if len(parts) >= 2:
                                conda_env_path = parts[-1]
                                break
            except (subprocess.TimeoutExpired, subprocess.SubprocessError):
                pass

        self._conda_executable = conda_executable
        self._conda_env_path = conda_env_path

        if self.debug:
            self.logger.debug(f"Conda executable: {conda_executable}")
            self.logger.debug(f"Conda env path: {conda_env_path}")

        return conda_executable, conda_env_path

    def discover_rdbase_path(self) -> Optional[str]:
        """
        Validate the provided RDKit build directory (RDBASE) path.

        Returns:
            Path to RDKit build directory or None if not valid
        """
        if self.validate_rdbase_path(self.rdbase_path):
            self._validated_rdbase = self.rdbase_path
            if self.debug:
                self.logger.debug(f"Validated RDBASE at: {self.rdbase_path}")
            return self.rdbase_path
        else:
            if self.debug:
                self.logger.debug(f"Invalid RDBASE path: {self.rdbase_path}")
            return None

    def validate_rdbase_path(self, rdbase_path: str) -> bool:
        """
        Validate that the RDBASE path contains required RDKit components.

        Args:
            rdbase_path: Path to validate

        Returns:
            True if valid RDBASE path
        """
        if not rdbase_path or not os.path.exists(rdbase_path):
            return False

        required_components = [
            "lib",
            "Code",
            "Data",
        ]

        # Check for Python module (can be in different locations)
        python_module_locations = [
            "rdkit",
            "Code/rdkit",
            "lib/python*/site-packages/rdkit",
        ]

        # Check required directories
        for component in required_components:
            component_path = os.path.join(rdbase_path, component)
            if not os.path.exists(component_path):
                if self.debug:
                    self.logger.debug(f"Missing required component: {component_path}")
                return False

        # Check for Python module
        python_module_found = False
        for location in python_module_locations:
            if "*" in location:
                import glob
                matches = glob.glob(os.path.join(rdbase_path, location))
                if matches:
                    python_module_found = True
                    break
            else:
                module_path = os.path.join(rdbase_path, location)
                if os.path.exists(module_path):
                    python_module_found = True
                    break

        if not python_module_found:
            if self.debug:
                self.logger.debug("RDKit Python module not found in any expected location")
            return False

        return True

    def validate_environment(self) -> Dict[str, bool]:
        """
        Validate the complete environment setup.

        Returns:
            Dictionary with validation results for each component
        """
        results = {}

        # Validate conda
        conda_exe, conda_env = self.discover_conda_installation()
        results['conda_executable'] = conda_exe is not None
        results['conda_environment'] = conda_env is not None

        # Validate RDBASE
        rdbase = self.discover_rdbase_path()
        results['rdbase_path'] = rdbase is not None

        # Validate CUDA
        results['cuda_available'] = os.path.exists(self.cuda_home_path)

        self._validation_results = results
        return results

    def generate_environment_script(self, output_path: Optional[str] = None) -> str:
        """
        Generate a shell script to set up the environment.

        Args:
            output_path: Path to save the script (optional)

        Returns:
            Path to the generated script
        """
        if output_path is None:
            fd, output_path = tempfile.mkstemp(suffix='.sh', prefix='roshambo_env_')
            os.close(fd)

        # Discover paths
        conda_exe, conda_env = self.discover_conda_installation()
        rdbase = self.discover_rdbase_path()

        script_content = f"""#!/bin/bash
# Auto-generated Roshambo environment setup script
# Generated by RoshamboDynamicEnvironmentManager

# Error handling
set -e

# Initialize conda
"""

        # Add conda initialization
        if conda_exe:
            conda_base = os.path.dirname(os.path.dirname(conda_exe))
            script_content += f"""
# Initialize conda from discovered installation
if [ -f "{conda_base}/etc/profile.d/conda.sh" ]; then
    source {conda_base}/etc/profile.d/conda.sh
elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source /opt/conda/etc/profile.d/conda.sh
fi

# Activate environment
conda activate {self.conda_env_name}
"""
        else:
            script_content += f"""
# Fallback conda activation
conda activate {self.conda_env_name} 2>/dev/null || {{
    echo "Error: Could not activate conda environment '{self.conda_env_name}'" >&2
    exit 1
}}
"""

        # Add environment variables
        if rdbase:
            script_content += f"""
# Set RDKit environment variables
export RDBASE="{rdbase}"
export RDKIT_LIB_DIR="$RDBASE/lib"
export RDKIT_INCLUDE_DIR="$RDBASE/Code"
export RDKIT_DATA_DIR="$RDBASE/Data"

# Clear and set PYTHONPATH
unset PYTHONPATH
export PYTHONPATH="$RDBASE:$CONDA_PREFIX/lib/python*/site-packages"

# Set library path
export LD_LIBRARY_PATH="$RDBASE/lib:$CONDA_PREFIX/lib:{self.cuda_home_path}/lib64:$LD_LIBRARY_PATH"
"""

        # Add CUDA configuration
        script_content += f"""
# Set CUDA environment
export CUDA_HOME="{self.cuda_home_path}"
export PATH="$CUDA_HOME/bin:$PATH"

# Verify critical imports (silent check)
python -c "
import sys
try:
    from rdkit import Chem
    import roshambo
except ImportError as e:
    print(f'Environment setup error: {{e}}', file=sys.stderr)
    sys.exit(1)
" 2>/dev/null || {{
    echo "Error: Failed to import required packages" >&2
    exit 1
}}
"""

        # Write script
        with open(output_path, 'w') as f:
            f.write(script_content)

        # Make executable
        os.chmod(output_path, 0o755)

        if self.debug:
            self.logger.debug(f"Generated environment script: {output_path}")

        return output_path

    def setup_environment_variables(self) -> Dict[str, str]:
        """
        Set up environment variables in the current process.

        Returns:
            Dictionary of environment variables that were set
        """
        env_vars = {}

        # Discover paths
        conda_exe, conda_env = self.discover_conda_installation()
        rdbase = self.discover_rdbase_path()

        if rdbase:
            env_vars['RDBASE'] = rdbase
            env_vars['RDKIT_LIB_DIR'] = os.path.join(rdbase, 'lib')
            env_vars['RDKIT_INCLUDE_DIR'] = os.path.join(rdbase, 'Code')
            env_vars['RDKIT_DATA_DIR'] = os.path.join(rdbase, 'Data')

            # Set PYTHONPATH
            current_pythonpath = os.environ.get('PYTHONPATH', '')
            new_pythonpath = rdbase
            if conda_env:
                new_pythonpath += f":{conda_env}/lib/python*/site-packages"
            env_vars['PYTHONPATH'] = new_pythonpath

            # Set LD_LIBRARY_PATH
            current_ld_path = os.environ.get('LD_LIBRARY_PATH', '')
            new_ld_path = f"{rdbase}/lib"
            if conda_env:
                new_ld_path += f":{conda_env}/lib"
            new_ld_path += f":{self.cuda_home_path}/lib64"
            if current_ld_path:
                new_ld_path += f":{current_ld_path}"
            env_vars['LD_LIBRARY_PATH'] = new_ld_path

        # Set CUDA variables
        env_vars['CUDA_HOME'] = self.cuda_home_path

        # Apply environment variables
        if self.auto_setup_env:
            for key, value in env_vars.items():
                os.environ[key] = value
                if self.debug:
                    self.logger.debug(f"Set {key}={value}")

        return env_vars

    def get_environment_summary(self) -> Dict[str, any]:
        """
        Get a summary of the current environment configuration.

        Returns:
            Dictionary with environment information
        """
        conda_exe, conda_env = self.discover_conda_installation()
        rdbase = self.discover_rdbase_path()
        validation = self.validate_environment()

        return {
            'conda_executable': conda_exe,
            'conda_environment_path': conda_env,
            'conda_environment_name': self.conda_env_name,
            'rdbase_path': rdbase,
            'cuda_home': self.cuda_home_path,
            'validation_results': validation,
            'auto_setup_enabled': self.auto_setup_env,
            'debug_enabled': self.debug
        }
