from typing import List
import os
import numpy as np
from pathlib import Path

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.enums.roshambo_specific_parameters_enum import RoshamboSpecificParametersEnum
from reinvent_scoring.scoring.score_components.roshambo.rdkit_conformer_generator import RoshamboConformerGenerator
from reinvent_scoring.scoring.score_components.roshambo.dynamic_environment_manager import RoshamboDynamicEnvironmentManager


class RoshamboShapeSimilarity(BaseScoreComponent):
    """
    Scoring component that uses Roshambo for GPU-accelerated molecular shape comparison.
    This component requires the Roshambo package to be installed.
    """

    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

        # Initialize parameter enum
        self.param_enum = RoshamboSpecificParametersEnum()

        # Debug mode
        self.debug = self.parameters.specific_parameters.get(
            self.param_enum.DEBUG, False
        )

        # Dynamic environment configuration - REQUIRED PARAMETERS
        self.rdbase_path = self.parameters.specific_parameters.get(self.param_enum.RDBASE_PATH)
        self.conda_env_name = self.parameters.specific_parameters.get(self.param_enum.CONDA_ENV_NAME)

        # Validate required parameters
        if not self.rdbase_path:
            error_msg = (
                "'rdbase_path' is required for Roshambo component. "
                "Please specify the path to your RDKit build directory in the configuration. "
                "Example: 'rdbase_path': '/home/user/Desktop/roshambo/rdkit'"
            )
            if self.debug:
                print(f"❌ Configuration Error: {error_msg}")
                print(f"Available parameters: {list(self.parameters.specific_parameters.keys())}")
            raise ValueError(error_msg)
        if not self.conda_env_name:
            error_msg = (
                "'conda_env_name' is required for Roshambo component. "
                "Please specify the name of your conda environment with Roshambo installed. "
                "Example: 'conda_env_name': 'roshambo'"
            )
            if self.debug:
                print(f"❌ Configuration Error: {error_msg}")
                print(f"Available parameters: {list(self.parameters.specific_parameters.keys())}")
            raise ValueError(error_msg)

        # Optional parameters with defaults
        self.conda_base_path = self.parameters.specific_parameters.get(self.param_enum.CONDA_BASE_PATH, None)
        self.cuda_home_path = self.parameters.specific_parameters.get(self.param_enum.CUDA_HOME_PATH, "/usr/local/cuda")
        self.auto_setup_env = self.parameters.specific_parameters.get(self.param_enum.AUTO_SETUP_ENV, True)

        # Debug output for configuration
        if self.debug:
            print("🔧 Roshambo Configuration:")
            print(f"  rdbase_path: {self.rdbase_path}")
            print(f"  conda_env_name: {self.conda_env_name}")
            print(f"  conda_base_path: {self.conda_base_path}")
            print(f"  cuda_home_path: {self.cuda_home_path}")
            print(f"  auto_setup_env: {self.auto_setup_env}")

        # Legacy environment path support (for backward compatibility)
        self.environment_path = self.parameters.specific_parameters.get(
            self.param_enum.ENVIRONMENT_PATH, ""
        )

        # Initialize dynamic environment manager
        self.env_manager = RoshamboDynamicEnvironmentManager(
            rdbase_path=self.rdbase_path,
            conda_env_name=self.conda_env_name,
            conda_base_path=self.conda_base_path,
            cuda_home_path=self.cuda_home_path,
            auto_setup_env=self.auto_setup_env,
            debug=self.debug
        )

        # Validate and setup environment
        self._setup_environment()

        # Import Roshambo here to avoid dependency issues if not installed
        try:
            if self.environment_path:
                # If environment path is provided, we'll use subprocess to call roshambo
                self.use_subprocess = True
            else:
                # Try direct import first
                from roshambo.api import get_similarity_scores
                self.get_similarity_scores = get_similarity_scores
                self.use_subprocess = False
        except ImportError:
            # If direct import fails, check if we have dynamic environment setup
            if self.conda_env_name and not self.environment_path:
                if self.debug:
                    print(f"Roshambo not available in current environment, will use conda environment: {self.conda_env_name}")
                # Use subprocess execution with dynamic environment
                self.use_subprocess = True
            else:
                raise ImportError("Roshambo package is required for this component. "
                                 "Please install it with: pip install git+https://github.com/molecularinformatics/roshambo.git "
                                 "or provide an environment_path to a conda environment with roshambo installed.")

        # Basic parameters using the enum
        self.shape_weight = self.parameters.specific_parameters.get(self.param_enum.SHAPE_WEIGHT, 0.5)
        self.color_weight = self.parameters.specific_parameters.get(self.param_enum.COLOR_WEIGHT, 0.5)
        self.reference_file = self.parameters.specific_parameters.get(self.param_enum.REFERENCE_FILE, "")

        # Roshambo specific parameters
        self.gpu_id = self.parameters.specific_parameters.get(self.param_enum.GPU_ID, 0)
        self.n_confs = self.parameters.specific_parameters.get(self.param_enum.N_CONFS, 0)  # 0 means use existing conformers
        self.ignore_hs = self.parameters.specific_parameters.get(self.param_enum.IGNORE_HS, True)
        self.use_carbon_radii = self.parameters.specific_parameters.get(self.param_enum.USE_CARBON_RADII, True)

        # Overlay saving parameters
        self.save_overlays = self.parameters.specific_parameters.get(self.param_enum.SAVE_OVERLAYS, False)
        self.overlays_dir = self.parameters.specific_parameters.get(self.param_enum.OVERLAYS_DIR, "roshambo_overlays")
        self.overlay_prefix = self.parameters.specific_parameters.get(self.param_enum.OVERLAY_PREFIX, "overlay")

        # Multiple reference files for PROTAC design
        self.warhead1_reference = self.parameters.specific_parameters.get(self.param_enum.WARHEAD1_REFERENCE, "")
        self.warhead2_reference = self.parameters.specific_parameters.get(self.param_enum.WARHEAD2_REFERENCE, "")



        # RDKit conformer generation parameters
        self.use_rdkit_conformers = self.n_confs > 0
        if self.use_rdkit_conformers:
            self.conformer_generator = RoshamboConformerGenerator(
                n_confs=self.n_confs,
                method=self.parameters.specific_parameters.get(self.param_enum.RDKIT_METHOD, "ETKDGv3"),
                random_seed=self.parameters.specific_parameters.get(self.param_enum.RDKIT_RANDOM_SEED, 42),
                ff=self.parameters.specific_parameters.get(self.param_enum.RDKIT_FF, "MMFF94s"),
                add_hs=self.parameters.specific_parameters.get(self.param_enum.RDKIT_ADD_HS, True),
                opt_confs=self.parameters.specific_parameters.get(self.param_enum.RDKIT_OPT_CONFS, True),
                num_threads=self.parameters.specific_parameters.get(self.param_enum.RDKIT_NUM_THREADS, 1)
            )

        # Create overlays directory if saving is enabled
        if self.save_overlays:
            Path(self.overlays_dir).mkdir(parents=True, exist_ok=True)

        # Validate reference file
        if not self.reference_file and not self.warhead1_reference and not self.warhead2_reference:
            raise ValueError("At least one reference file must be provided")

        for ref_file in [self.reference_file, self.warhead1_reference, self.warhead2_reference]:
            if ref_file and not os.path.exists(ref_file):
                raise FileNotFoundError(f"Reference file not found: {ref_file}")

    def _setup_environment(self):
        """Set up and validate the Roshambo environment."""
        try:
            # Validate environment
            validation_results = self.env_manager.validate_environment()

            if self.debug:
                print("Environment validation results:")
                for component, status in validation_results.items():
                    status_str = "✅" if status else "❌"
                    print(f"  {component}: {status_str}")

            # Check for critical failures
            critical_failures = []
            if not validation_results.get('conda_environment', False):
                critical_failures.append(f"Conda environment '{self.conda_env_name}' not found")

            if not validation_results.get('rdbase_path', False):
                if self.rdbase_path:
                    critical_failures.append(f"Specified RDBASE path not valid: {self.rdbase_path}")
                else:
                    critical_failures.append("No valid RDBASE path found")

            if critical_failures:
                error_msg = "Environment setup failed:\n" + "\n".join(f"  - {failure}" for failure in critical_failures)

                # Provide helpful suggestions
                env_summary = self.env_manager.get_environment_summary()
                error_msg += "\n\nEnvironment Summary:"
                error_msg += f"\n  Conda executable: {env_summary.get('conda_executable', 'Not found')}"
                error_msg += f"\n  Conda environment: {env_summary.get('conda_environment_path', 'Not found')}"
                error_msg += f"\n  RDBASE path: {env_summary.get('rdbase_path', 'Not found')}"

                error_msg += "\n\nSuggestions:"
                error_msg += f"\n  1. Ensure conda environment '{self.conda_env_name}' exists"
                error_msg += "\n  2. Specify 'rdbase_path' in configuration if RDKit is in a custom location"
                error_msg += "\n  3. Set 'conda_base_path' if conda is installed in a non-standard location"

                raise EnvironmentError(error_msg)

            # Set up environment variables if auto_setup is enabled
            if self.auto_setup_env:
                env_vars = self.env_manager.setup_environment_variables()
                if self.debug:
                    print("Set environment variables:")
                    for key, value in env_vars.items():
                        print(f"  {key}={value[:100]}{'...' if len(value) > 100 else ''}")

            # Generate environment script for subprocess execution if needed
            # Always generate if we might need subprocess execution
            if self.environment_path or not self.auto_setup_env or self.conda_env_name:
                self.dynamic_env_script = self.env_manager.generate_environment_script()
                if self.debug:
                    print(f"Generated dynamic environment script: {self.dynamic_env_script}")
            else:
                self.dynamic_env_script = None

        except Exception as e:
            error_msg = f"Failed to set up Roshambo environment: {e}"
            if self.debug:
                import traceback
                error_msg += f"\n\nFull traceback:\n{traceback.format_exc()}"
            raise EnvironmentError(error_msg)

    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        """Calculate shape similarity scores for a list of molecules."""
        if self.debug:
            print(f"🧮 Roshambo calculate_score called with {len(molecules)} molecules, step={step}")
            print(f"  Component initialized with rdbase_path: {self.rdbase_path}")
            print(f"  Component initialized with conda_env_name: {self.conda_env_name}")

        # Convert input to SMILES strings
        smiles_list = []
        for mol in molecules:
            if isinstance(mol, str):
                smiles_list.append(mol)
            else:
                try:
                    from rdkit import Chem
                    smiles = Chem.MolToSmiles(mol)
                    smiles_list.append(smiles)
                except:
                    smiles_list.append("")  # Empty string will result in zero score

        # Calculate scores
        scores = self._calculate_shape_scores(smiles_list, step)

        # Create and return component summary
        score_summary = ComponentSummary(total_score=scores, parameters=self.parameters)
        return score_summary

    def calculate_score_for_step(self, molecules: List, step=-1) -> ComponentSummary:
        """Calculate shape similarity scores for a specific step (used in reinforcement learning)."""
        return self.calculate_score(molecules, step)

    def _get_similarity_scores_with_env(self, **kwargs):
        """
        Execute Roshambo using subprocess with conda environment.
        Simplified version that focuses on CSV output.
        """
        import json
        import tempfile
        import subprocess

        # Create temporary files for input/output
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as input_file:
            json.dump(kwargs, input_file)
            input_path = input_file.name

        if self.debug:
            print(f"Roshambo subprocess parameters: {kwargs}")

        output_path = input_path.replace('.json', '_output.json')

        try:
            # Create simplified Python script
            script_content = f'''
import json
import sys
import os
import pandas as pd

# Set up environment variables
os.environ["RDBASE"] = "{self.rdbase_path}"
os.environ["RDKIT_LIB_DIR"] = "{self.rdbase_path}/lib"
os.environ["RDKIT_INCLUDE_DIR"] = "{self.rdbase_path}/Code"
os.environ["RDKIT_DATA_DIR"] = "{self.rdbase_path}/Data"
os.environ["CUDA_HOME"] = "{self.cuda_home_path}"

# Set up PYTHONPATH
current_pythonpath = os.environ.get("PYTHONPATH", "")
new_pythonpath = "{self.rdbase_path}"
if current_pythonpath:
    new_pythonpath += ":" + current_pythonpath
os.environ["PYTHONPATH"] = new_pythonpath

# Set up LD_LIBRARY_PATH
current_ld_path = os.environ.get("LD_LIBRARY_PATH", "")
new_ld_path = "{self.rdbase_path}/lib"
if current_ld_path:
    new_ld_path += ":" + current_ld_path
os.environ["LD_LIBRARY_PATH"] = new_ld_path

# Add RDKit to Python path
sys.path.insert(0, "{self.rdbase_path}")

try:
    from roshambo.api import get_similarity_scores

    # Load input parameters
    with open("{input_path}", "r") as f:
        kwargs = json.load(f)

    # Force write_to_file=True
    kwargs['write_to_file'] = True

    # Change to working directory if specified
    working_dir = kwargs.get("working_dir")
    if working_dir and os.path.exists(working_dir):
        os.chdir(working_dir)

    # Execute roshambo (returns None but writes CSV)
    get_similarity_scores(**kwargs)

    # Read CSV file and convert to records
    csv_files = [f for f in os.listdir('.') if f.endswith('.csv')]
    if csv_files:
        df = pd.read_csv(csv_files[0], sep='\\t')
        result_dict = df.to_dict('records')
    else:
        result_dict = None

    # Save result
    with open("{output_path}", "w") as f:
        json.dump(result_dict, f)

except Exception as e:
    import traceback
    error_info = {{
        "error": str(e),
        "traceback": traceback.format_exc()
    }}
    with open("{output_path}", "w") as f:
        json.dump(error_info, f)
    sys.exit(1)
'''

            # Create temporary script file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as script_file:
                script_file.write(script_content)
                script_path = script_file.name

            # Execute script
            if self.environment_path:
                full_command = f"{self.environment_path} python {script_path}"
            else:
                full_command = f"conda run -n {self.conda_env_name} python {script_path}"

            if self.debug:
                print(f"Executing: {full_command}")

            result = subprocess.run(full_command, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                error_msg = f"Roshambo subprocess failed: {result.stderr}"
                if self.debug:
                    print(f"Command output: {result.stdout}")
                    print(f"Command error: {result.stderr}")
                raise RuntimeError(error_msg)

            # Load result
            output = None
            try:
                if os.path.exists(output_path):
                    with open(output_path, 'r') as f:
                        content = f.read().strip()
                        if content:
                            output = json.loads(content)
            except Exception as e:
                if self.debug:
                    print(f"Error reading output file: {e}")

            if output and isinstance(output, dict) and 'error' in output:
                raise RuntimeError(f"Roshambo error: {output['error']}")

            # Convert list of records back to DataFrame-like object
            if output and isinstance(output, list):
                class DataFrameProxy:
                    def __init__(self, records):
                        self.records = records

                    def iterrows(self):
                        for i, record in enumerate(self.records):
                            yield i, record

                return DataFrameProxy(output)

            return None

        finally:
            # Cleanup temporary files
            for path in [input_path, output_path, script_path]:
                try:
                    os.unlink(path)
                except:
                    pass

    def _calculate_shape_scores(self, smiles_list: List[str], step: int) -> np.array:
        """Calculate shape similarity scores using Roshambo."""
        if self.debug:
            print(f"Processing {len(smiles_list)} SMILES: {smiles_list[:5]}...")  # Show first 5

        if not smiles_list:
            return np.array([], dtype=np.float32)

        # Create a temporary directory
        temp_dir = os.path.join(self.overlays_dir if self.save_overlays else ".", "temp")
        Path(temp_dir).mkdir(parents=True, exist_ok=True)

        # Prepare input file
        if self.use_rdkit_conformers:
            # Use RDKit to generate conformers and create SDF file
            input_file = self.conformer_generator.create_temp_sdf_from_smiles(smiles_list, temp_dir)
        else:
            # Create SMILES file for Roshambo to handle conformer generation
            input_file = os.path.join(temp_dir, f"query_molecules_{step}.smi")
            with open(input_file, "w") as f:
                for i, smi in enumerate(smiles_list):
                    if smi:  # Skip empty SMILES
                        f.write(f"{smi} mol_{i}\n")

        if self.debug:
            print(f"Created input file: {input_file}")
            if os.path.exists(input_file):
                with open(input_file, 'r') as f:
                    content = f.read()
                    print(f"Input file content (first 500 chars): {content[:500]}")
            else:
                print(f"ERROR: Input file {input_file} was not created!")

        # Initialize scores array
        scores = [0.0] * len(smiles_list)

        # Process each reference file
        ref_files = [f for f in [self.reference_file, self.warhead1_reference, self.warhead2_reference] if f]

        if self.debug:
            print(f"Using reference files: {ref_files}")
            for ref_file in ref_files:
                if os.path.exists(ref_file):
                    print(f"  ✓ {ref_file} exists")
                    # Check if the file has content
                    try:
                        with open(ref_file, 'r') as f:
                            content = f.read()
                            if content.strip():
                                print(f"    File size: {len(content)} chars, first 100 chars: {content[:100]}")
                            else:
                                print(f"    WARNING: File is empty!")
                    except Exception as e:
                        print(f"    ERROR reading file: {e}")
                else:
                    print(f"  ✗ {ref_file} does not exist")

        for ref_file in ref_files:
            try:
                # Set up working directory for this reference
                ref_name = os.path.basename(ref_file).split('.')[0]
                working_dir = os.path.join(self.overlays_dir if self.save_overlays else temp_dir,
                                          f"{ref_name}_step{step}")

                # Run Roshambo
                # When using RDKit conformers, set n_confs=0 to avoid double conformer generation
                roshambo_n_confs = 0 if self.use_rdkit_conformers else self.n_confs

                # Copy files to working directory and use relative paths (like in working example)
                if working_dir:
                    os.makedirs(working_dir, exist_ok=True)

                    # Convert and copy reference file to working directory (roshambo-compatible format)
                    ref_basename = os.path.basename(ref_file)
                    local_ref_file = os.path.join(working_dir, ref_basename)
                    self._convert_sdf_for_roshambo(ref_file, local_ref_file)

                    # Convert and copy dataset file to working directory
                    dataset_basename = os.path.basename(input_file)
                    local_dataset_file = os.path.join(working_dir, dataset_basename)
                    self._convert_sdf_for_roshambo(input_file, local_dataset_file)

                    # Use relative paths from working directory
                    rel_ref_file = ref_basename
                    rel_dataset_file = dataset_basename
                    rel_working_dir = working_dir
                else:
                    # Fallback to absolute paths
                    rel_ref_file = ref_file
                    rel_dataset_file = input_file
                    rel_working_dir = None

                if self.debug:
                    print(f"Calling roshambo with ref_file={rel_ref_file}, dataset_file={rel_dataset_file}, n_confs={roshambo_n_confs}")
                    print(f"Working directory: {rel_working_dir}")

                # Execute roshambo (always with write_to_file=True due to API bug)
                if self.use_subprocess:
                    mol_scores = self._run_roshambo_subprocess(rel_ref_file, rel_dataset_file, roshambo_n_confs, rel_working_dir)
                else:
                    mol_scores = self._run_roshambo_direct(rel_ref_file, rel_dataset_file, roshambo_n_confs, rel_working_dir)

                # Update scores with the best score for each molecule
                for idx, score in mol_scores.items():
                    if 0 <= idx < len(scores):
                        scores[idx] = max(scores[idx], score)

                if self.debug:
                    print(f"Scores from {ref_file}: {mol_scores}")

            except Exception as e:
                error_msg = f"Error calculating Roshambo similarity with {ref_file}: {e}"
                print(f"❌ {error_msg}")

                if self.debug:
                    import traceback
                    print(f"Full traceback: {traceback.format_exc()}")
                    print(f"Configuration used:")
                    print(f"  rdbase_path: {self.rdbase_path}")
                    print(f"  conda_env_name: {self.conda_env_name}")
                    print(f"  use_subprocess: {self.use_subprocess}")
                    print(f"  reference file exists: {os.path.exists(ref_file)}")

                # If roshambo fails, assign small non-zero scores to ensure diversity filter gets data
                print(f"⚠️  Roshambo failed for {ref_file}, assigning fallback scores")
                for i in range(len(smiles_list)):
                    if smiles_list[i]:  # Only assign scores to valid SMILES
                        scores[i] = max(scores[i], 0.001)  # Small non-zero score

        # Clean up temporary files if not saving
        if not self.save_overlays:
            try:
                import shutil
                shutil.rmtree(temp_dir)
            except:
                pass

        # Ensure we have a score for each input molecule
        while len(scores) < len(smiles_list):
            scores.append(0.0)

        if self.debug:
            print(f"Final scores for {len(smiles_list)} molecules: {scores}")
            non_zero_count = sum(1 for s in scores if s > 0)
            print(f"Non-zero scores: {non_zero_count}/{len(scores)}")

        return np.array(scores, dtype=np.float32)

    def _run_roshambo_direct(self, ref_file: str, dataset_file: str, n_confs: int, working_dir: str) -> dict:
        """Run roshambo directly and extract scores from CSV."""
        try:
            # Execute roshambo API (always returns None due to bug, but writes CSV)
            self.get_similarity_scores(
                ref_file=ref_file,
                dataset_file=dataset_file,  # Fixed parameter name
                ignore_hs=self.ignore_hs,
                n_confs=n_confs,
                use_carbon_radii=self.use_carbon_radii,
                color=self.color_weight > 0,
                sort_by="ComboTanimoto",
                write_to_file=True,  # Always write to file due to API bug
                gpu_id=self.gpu_id,
                working_dir=working_dir
            )

            # Read scores from CSV file
            return self._read_scores_from_csv(working_dir)

        except Exception as e:
            if self.debug:
                print(f"Error running roshambo directly: {e}")
            return {}

    def _run_roshambo_subprocess(self, ref_file: str, dataset_file: str, n_confs: int, working_dir: str) -> dict:
        """Run roshambo via subprocess and extract scores from CSV."""
        try:
            # Use the existing subprocess method
            results = self._get_similarity_scores_with_env(
                ref_file=ref_file,
                dataset_file=dataset_file,  # Fixed parameter name
                ignore_hs=self.ignore_hs,
                n_confs=n_confs,
                use_carbon_radii=self.use_carbon_radii,
                color=self.color_weight > 0,
                sort_by="ComboTanimoto",
                write_to_file=True,
                gpu_id=self.gpu_id,
                working_dir=working_dir
            )

            # If subprocess returned parsed results, use them
            if results and hasattr(results, 'iterrows'):
                return self._extract_scores_from_dataframe(results)
            else:
                # Otherwise read from CSV file
                return self._read_scores_from_csv(working_dir)

        except Exception as e:
            if self.debug:
                print(f"Error running roshambo via subprocess: {e}")
            return {}

    def _read_scores_from_csv(self, working_dir: str) -> dict:
        """Read and parse scores from roshambo CSV output."""
        mol_scores = {}

        if not working_dir:
            return mol_scores

        try:
            import glob
            import pandas as pd

            csv_files = glob.glob(os.path.join(working_dir, "*.csv"))

            if self.debug:
                print(f"Looking for CSV files in {working_dir}: {csv_files}")

            for csv_file in csv_files:
                try:
                    df = pd.read_csv(csv_file, sep='\t')

                    if self.debug:
                        print(f"Reading CSV {csv_file}: {df.shape[0]} rows, {df.shape[1]} columns")
                        print(f"Columns: {df.columns.tolist()}")

                    mol_scores.update(self._extract_scores_from_dataframe(df))

                except Exception as e:
                    if self.debug:
                        print(f"Error reading CSV file {csv_file}: {e}")
                    continue

        except Exception as e:
            if self.debug:
                print(f"Error looking for CSV files: {e}")

        return mol_scores

    def _extract_scores_from_dataframe(self, df) -> dict:
        """Extract molecule scores from DataFrame or DataFrame-like object."""
        mol_scores = {}

        try:
            for _, row in df.iterrows():
                name = row.get("Molecule", "")
                if name and name.startswith("mol_"):
                    try:
                        # Extract molecule index from name like "mol_0_0" -> 0
                        parts = name.split("_")
                        if len(parts) >= 2:
                            idx = int(parts[1])
                            shape_score = float(row.get("ShapeTanimoto", 0.0))
                            color_score = float(row.get("ColorTanimoto", 0.0))

                            # Calculate weighted combination
                            combo_score = (self.shape_weight * shape_score +
                                         self.color_weight * color_score) / (self.shape_weight + self.color_weight)

                            mol_scores[idx] = max(mol_scores.get(idx, 0.0), combo_score)

                            if self.debug:
                                print(f"  Molecule {idx}: shape={shape_score:.3f}, color={color_score:.3f}, combo={combo_score:.3f}")

                    except Exception as e:
                        if self.debug:
                            print(f"  Error processing row for {name}: {e}")
                        continue

        except Exception as e:
            if self.debug:
                print(f"Error extracting scores from DataFrame: {e}")

        return mol_scores

    def _convert_sdf_for_roshambo(self, input_sdf: str, output_sdf: str):
        """
        Convert SDF file to a format that roshambo can read properly.

        The issue is that roshambo's prepare_mols function incorrectly tries to parse
        some SDF files as SMILES. This function creates a clean SDF file that roshambo
        can read correctly.
        """
        try:
            from rdkit import Chem

            if self.debug:
                print(f"Converting SDF for roshambo: {input_sdf} -> {output_sdf}")

            # Read molecules from input SDF
            suppl = Chem.SDMolSupplier(input_sdf)
            molecules = []

            for i, mol in enumerate(suppl):
                if mol is not None:
                    molecules.append(mol)
                    if self.debug:
                        print(f"  Read molecule {i}: {mol.GetNumAtoms()} atoms")

            if not molecules:
                if self.debug:
                    print(f"  No valid molecules found in {input_sdf}")
                # Copy original file as fallback
                import shutil
                shutil.copy2(input_sdf, output_sdf)
                return

            # Write molecules to output SDF with clean format
            writer = Chem.SDWriter(output_sdf)

            for i, mol in enumerate(molecules):
                # Set a clean molecule name
                mol.SetProp("_Name", f"mol_{i}")

                # Ensure molecule has 3D coordinates
                if mol.GetNumConformers() == 0:
                    # Generate 3D coordinates if missing
                    from rdkit.Chem import AllChem
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol)
                    AllChem.MMFFOptimizeMolecule(mol)

                writer.write(mol)

            writer.close()

            if self.debug:
                print(f"  Successfully converted {len(molecules)} molecules")

        except Exception as e:
            if self.debug:
                print(f"  Error converting SDF: {e}")
            # Fallback: copy original file
            import shutil
            shutil.copy2(input_sdf, output_sdf)


