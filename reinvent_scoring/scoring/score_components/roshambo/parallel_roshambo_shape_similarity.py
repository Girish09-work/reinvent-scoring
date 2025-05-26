from typing import List
import os
import numpy as np
from pathlib import Path
import multiprocessing
from multiprocessing import Pool

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.enums.roshambo_specific_parameters_enum import RoshamboSpecificParametersEnum


class ParallelRoshamboShapeSimilarity(RoshamboShapeSimilarity):
    """
    Parallel version of the Roshambo shape similarity component.
    This distributes molecules across multiple GPUs for faster processing.
    """

    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)

        # Parallel processing parameters using the enum
        self.max_gpus = self.parameters.specific_parameters.get(self.param_enum.MAX_GPUS, 1)
        self.gpu_ids = self.parameters.specific_parameters.get(self.param_enum.GPU_IDS, [0])

        # If gpu_ids not specified, use sequential IDs up to max_gpus
        if not self.gpu_ids:
            self.gpu_ids = list(range(self.max_gpus))

        # Limit to available GPUs
        try:
            import torch
            num_available_gpus = torch.cuda.device_count()
            self.gpu_ids = [gpu_id for gpu_id in self.gpu_ids if gpu_id < num_available_gpus]
            if self.debug:
                print(f"Available GPUs: {num_available_gpus}, Using GPUs: {self.gpu_ids}")
        except:
            # If torch not available or no CUDA, just use the provided IDs
            if self.debug:
                print("PyTorch not available or no CUDA support detected")
            pass

        if not self.gpu_ids:
            self.gpu_ids = [0]  # Fallback to GPU 0

        # Number of processes to use (one per GPU)
        self.num_processes = min(len(self.gpu_ids), multiprocessing.cpu_count())

    def _calculate_shape_scores(self, smiles_list: List[str], step: int) -> np.array:
        """Calculate shape similarity scores using Roshambo in parallel across multiple GPUs."""
        if not smiles_list:
            return np.array([], dtype=np.float32)

        # Split molecules into chunks for parallel processing
        chunk_size = max(1, len(smiles_list) // self.num_processes)
        chunks = [smiles_list[i:i + chunk_size] for i in range(0, len(smiles_list), chunk_size)]

        # Create a temporary directory for each process
        temp_dir = os.path.join(self.overlays_dir if self.save_overlays else ".", "temp")
        Path(temp_dir).mkdir(parents=True, exist_ok=True)

        # Prepare input for each process
        inputs = []
        for i, chunk in enumerate(chunks):
            gpu_id = self.gpu_ids[i % len(self.gpu_ids)]

            # Create a temporary file with the SMILES for this chunk
            chunk_file = os.path.join(temp_dir, f"query_molecules_{step}_chunk{i}.smi")
            with open(chunk_file, "w") as f:
                for j, smi in enumerate(chunk):
                    if smi:  # Skip empty SMILES
                        f.write(f"{smi} mol_{i}_{j}\n")

            # Process each reference file
            ref_files = [f for f in [self.reference_file, self.warhead1_reference, self.warhead2_reference] if f]

            for ref_file in ref_files:
                ref_name = os.path.basename(ref_file).split('.')[0]
                working_dir = os.path.join(self.overlays_dir if self.save_overlays else temp_dir,
                                          f"{ref_name}_step{step}_chunk{i}")

                inputs.append({
                    "ref_file": ref_file,
                    "smiles_file": chunk_file,
                    "gpu_id": gpu_id,
                    "working_dir": working_dir,
                    "chunk_index": i,
                    "chunk_size": len(chunk),
                    "original_indices": list(range(i * chunk_size, min(i * chunk_size + len(chunk), len(smiles_list))))
                })

        # Process chunks in parallel
        all_scores = [0.0] * len(smiles_list)  # Initialize with zeros

        if self.num_processes > 1:
            with Pool(processes=self.num_processes) as pool:
                results = pool.map(self._process_chunk, inputs)

            # Combine results
            for chunk_scores, original_indices in results:
                for idx, score in zip(original_indices, chunk_scores):
                    if 0 <= idx < len(all_scores):
                        all_scores[idx] = max(all_scores[idx], score)
        else:
            # Process sequentially if only one process
            for input_data in inputs:
                chunk_scores, original_indices = self._process_chunk(input_data)
                for idx, score in zip(original_indices, chunk_scores):
                    if 0 <= idx < len(all_scores):
                        all_scores[idx] = max(all_scores[idx], score)

        # Clean up temporary files if not saving
        if not self.save_overlays:
            try:
                import shutil
                shutil.rmtree(temp_dir)
            except:
                pass

        return np.array(all_scores, dtype=np.float32)

    def _process_chunk(self, input_data):
        """Process a chunk of molecules on a specific GPU."""
        ref_file = input_data["ref_file"]
        smiles_file = input_data["smiles_file"]
        gpu_id = input_data["gpu_id"]
        working_dir = input_data["working_dir"]
        original_indices = input_data["original_indices"]

        chunk_scores = [0.0] * len(original_indices)

        try:
            # Run Roshambo on this chunk
            results = self.get_similarity_scores(
                ref_file=ref_file,
                dataset_files_pattern=smiles_file,
                ignore_hs=self.ignore_hs,
                n_confs=self.n_confs,
                use_carbon_radii=self.use_carbon_radii,
                color=self.color_weight > 0,
                sort_by="ComboTanimoto",
                write_to_file=self.save_overlays,
                gpu_id=gpu_id,
                working_dir=working_dir if self.save_overlays else None
            )

            # Extract scores
            if results and hasattr(results, "scores") and len(results.scores) > 0:
                # Map results back to original molecules
                for result in results.scores:
                    name = result.get("Name", "")
                    if name.startswith("mol_"):
                        try:
                            parts = name.split("_")
                            if len(parts) >= 3:
                                mol_idx = int(parts[2])

                                shape_score = result.get("ShapeTanimoto", 0.0)
                                color_score = result.get("ColorTanimoto", 0.0)
                                combo_score = (self.shape_weight * shape_score +
                                              self.color_weight * color_score) / (self.shape_weight + self.color_weight)

                                if 0 <= mol_idx < len(chunk_scores):
                                    chunk_scores[mol_idx] = max(chunk_scores[mol_idx], combo_score)
                        except:
                            continue

        except Exception as e:
            print(f"Error processing chunk on GPU {gpu_id} with {ref_file}: {e}")

        return chunk_scores, original_indices
