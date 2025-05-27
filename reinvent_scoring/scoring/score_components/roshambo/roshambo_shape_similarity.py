"""
Simplified Roshambo Shape Similarity Component for Reinvent Scoring.
This component uses the Roshambo Flask API for GPU-accelerated molecular shape comparison.
"""

import os
import pandas as pd
import numpy as np
import requests
from pathlib import Path
from typing import List

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.enums.roshambo_specific_parameters_enum import RoshamboSpecificParametersEnum
from reinvent_scoring.scoring.score_components.roshambo.rdkit_conformer_generator import RoshamboConformerGenerator


class RoshamboShapeSimilarity(BaseScoreComponent):
    """
    Simplified Roshambo Shape Similarity Component.
    Uses Flask API for GPU-accelerated molecular shape comparison.
    """

    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

        # Initialize parameter enum
        self.param_enum = RoshamboSpecificParametersEnum()

        # Basic parameters
        self.reference_file = self.parameters.specific_parameters.get(self.param_enum.REFERENCE_FILE, "")
        self.shape_weight = self.parameters.specific_parameters.get(self.param_enum.SHAPE_WEIGHT, 0.5)
        self.color_weight = self.parameters.specific_parameters.get(self.param_enum.COLOR_WEIGHT, 0.5)
        self.n_confs = self.parameters.specific_parameters.get(self.param_enum.N_CONFS, 0)
        self.ignore_hs = self.parameters.specific_parameters.get(self.param_enum.IGNORE_HS, True)
        self.use_carbon_radii = self.parameters.specific_parameters.get(self.param_enum.USE_CARBON_RADII, True)
        self.gpu_id = self.parameters.specific_parameters.get(self.param_enum.GPU_ID, 0)

        # Flask API configuration
        self.roshambo_api_url = self.parameters.specific_parameters.get("roshambo_api_url", "http://127.0.0.1:5000/")

        # Overlay saving parameters
        self.save_overlays = self.parameters.specific_parameters.get(self.param_enum.SAVE_OVERLAYS, True)
        self.overlays_dir = self.parameters.specific_parameters.get(self.param_enum.OVERLAYS_DIR, "roshambo_overlays")

        # Debug mode
        self.debug = self.parameters.specific_parameters.get("debug", False)

        # RDKit conformer generation
        self.use_rdkit_conformers = self.n_confs > 0
        if self.use_rdkit_conformers:
            self.conformer_generator = RoshamboConformerGenerator(
                n_confs=self.n_confs,
                method=self.parameters.specific_parameters.get("rdkit_method", "ETKDGv3"),
                random_seed=self.parameters.specific_parameters.get("rdkit_random_seed", 42),
                add_hs=self.parameters.specific_parameters.get("rdkit_add_hs", True),
                opt_confs=self.parameters.specific_parameters.get("rdkit_opt_confs", True),
                num_threads=self.parameters.specific_parameters.get("rdkit_num_threads", 1)
            )

        # Create overlays directory
        if self.save_overlays:
            Path(self.overlays_dir).mkdir(parents=True, exist_ok=True)

        # Validate reference file
        if not self.reference_file:
            raise ValueError("Reference file must be provided")
        if not os.path.exists(self.reference_file):
            raise FileNotFoundError(f"Reference file not found: {self.reference_file}")

        # Initialize epoch counter
        self.current_epoch = 0

    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        """Calculate shape similarity scores for a list of molecules."""
        if self.debug:
            print(f"🧮 Roshambo calculate_score called with {len(molecules)} molecules, step={step}")

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

    def _calculate_shape_scores(self, smiles_list: List[str], step: int) -> np.array:
        """Calculate shape similarity scores using Roshambo Flask API."""
        if self.debug:
            print(f"Processing {len(smiles_list)} SMILES for step {step}")

        if not smiles_list:
            return np.array([], dtype=np.float32)

        # Create epoch folder
        epoch_folder = os.path.join(self.overlays_dir, f"epoch_{step}")
        Path(epoch_folder).mkdir(parents=True, exist_ok=True)

        # Prepare dataset file
        if self.use_rdkit_conformers:
            # Use RDKit to generate conformers and create SDF file
            dataset_file = self.conformer_generator.create_temp_sdf_from_smiles(smiles_list, epoch_folder)
        else:
            # Create SMILES file for Roshambo to handle conformer generation
            dataset_file = os.path.join(epoch_folder, f"dataset_{step}.smi")
            with open(dataset_file, "w") as f:
                for i, smi in enumerate(smiles_list):
                    if smi:  # Skip empty SMILES
                        f.write(f"{smi} mol_{i}\n")

        if self.debug:
            print(f"Created dataset file: {dataset_file}")

        # Call Roshambo Flask API
        try:
            scores = self._call_roshambo_api(self.reference_file, dataset_file, epoch_folder)
            return np.array(scores, dtype=np.float32)
        except Exception as e:
            if self.debug:
                print(f"Error calling Roshambo API: {e}")
            return np.zeros(len(smiles_list), dtype=np.float32)

    def _call_roshambo_api(self, reference_file: str, dataset_file: str, epoch_folder: str) -> List[float]:
        """Call Roshambo Flask API and return scores."""
        try:
            # Prepare API request
            api_data = {
                "reference_file": reference_file,
                "dataset_file": dataset_file,
                "ignore_hs": self.ignore_hs,
                "n_confs": 0,  # Use existing conformers
                "use_carbon_radii": self.use_carbon_radii,
                "color": self.color_weight > 0,
                "sort_by": "ComboTanimoto",
                "write_to_file": True,
                "gpu_id": self.gpu_id,
                "working_dir": epoch_folder
            }

            if self.debug:
                print(f"Calling Roshambo API at {self.roshambo_api_url}/similarity")
                print(f"API data: {api_data}")

            # Call Flask API
            response = requests.post(
                f"{self.roshambo_api_url}/similarity",
                json=api_data,
                timeout=300  # 5 minute timeout
            )

            if response.status_code == 200:
                result = response.json()
                if result.get("success"):
                    # Read CSV file and extract scores
                    csv_file = os.path.join(epoch_folder, "roshambo.csv")
                    if os.path.exists(csv_file):
                        return self._extract_scores_from_csv(csv_file)
                    else:
                        if self.debug:
                            print(f"CSV file not found: {csv_file}")
                        return []
                else:
                    if self.debug:
                        print(f"API returned error: {result.get('error', 'Unknown error')}")
                    return []
            else:
                if self.debug:
                    print(f"API request failed with status {response.status_code}: {response.text}")
                return []

        except Exception as e:
            if self.debug:
                print(f"Error calling Roshambo API: {e}")
            return []

    def _extract_scores_from_csv(self, csv_file: str) -> List[float]:
        """Extract scores from Roshambo CSV output."""
        scores = []
        try:
            df = pd.read_csv(csv_file, sep='\t')

            if self.debug:
                print(f"Reading CSV {csv_file}: {df.shape[0]} rows")
                print(f"Columns: {df.columns.tolist()}")

            # Create a mapping from molecule names to scores
            mol_scores = {}
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

                    except Exception as e:
                        if self.debug:
                            print(f"Error processing row for {name}: {e}")
                        continue

            # Convert to ordered list
            max_idx = max(mol_scores.keys()) if mol_scores else -1
            scores = [mol_scores.get(i, 0.0) for i in range(max_idx + 1)]

        except Exception as e:
            if self.debug:
                print(f"Error reading CSV file {csv_file}: {e}")

        return scores