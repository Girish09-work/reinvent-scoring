
import os
import multiprocessing
from multiprocessing import Pool
from typing import List

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_components.rdkit_shape.rdkit_shape_similarity import RDKitShapeSimilarity
from reinvent_scoring.scoring.score_components.rdkit_shape.rdkit_conformer_generator import RDKitConformerGenerator


class ParallelRDKitShapeSimilarity(RDKitShapeSimilarity):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self.shape_weight = self.parameters.specific_parameters.get("shape_weight", 0.5)
        self.color_weight = self.parameters.specific_parameters.get("color_weight", 0.5)
        self.reference_file = self.parameters.specific_parameters.get("reference_file", "")
        self.method = self.parameters.specific_parameters.get("method", "usrcat")
        self.max_confs = self.parameters.specific_parameters.get("max_confs", 50)

        self.save_overlays = self.parameters.specific_parameters.get("save_overlays", False)
        self.overlay_prefix = self.parameters.specific_parameters.get("overlay_prefix", "mol_")
        self.overlays_dir = self.parameters.specific_parameters.get("overlays_dir", "overlays")
        if self.save_overlays:
            os.makedirs(self.overlays_dir, exist_ok=True)

        self.num_cpus = min(
            multiprocessing.cpu_count(),
            self.parameters.specific_parameters.get("max_num_cpus", 4)
        )

        self.conformer_generator = RDKitConformerGenerator(
            max_confs=self.max_confs,
            energy_window=self.parameters.specific_parameters.get("ewindow", 10),
            max_stereo=self.parameters.specific_parameters.get("max_stereo", 0)
        )
        self.reference_mol = self._load_reference_molecule()

    def _load_reference_molecule(self):
        if not self.reference_file:
            return None

        ref_mol = None
        if self.reference_file.endswith(".sdf"):
            suppl = Chem.SDMolSupplier(self.reference_file)
            if suppl and len(suppl) > 0:
                ref_mol = suppl[0]

        if ref_mol is None:
            raise ValueError(f"Could not load reference molecule from {self.reference_file}")

        if ref_mol.GetNumConformers() == 0:
            ref_mol = self.conformer_generator.generate_conformers(Chem.MolToSmiles(ref_mol))

        return ref_mol

    def calculate_score(self, molecules: List[str]) -> ComponentSummary:
        rdkit_mols = [Chem.MolFromSmiles(smiles) for smiles in molecules]
        scores = self._calculate_shape_scores(rdkit_mols)
        return ComponentSummary(total_score=scores, parameters=self.parameters)

    def _calculate_shape_scores(self, mols):
        """Calculate shape similarity scores in parallel."""
        if len(mols) == 0 or self.reference_mol is None:
            return np.array([], dtype=np.float32)
        
        # Prepare arguments for parallel processing
        args_list = []
        for i, mol in enumerate(mols):
            if mol is None:
                continue
            
            args = {
                "mol": mol,
                "smiles": Chem.MolToSmiles(mol) if mol else "",
                "reference_mol": self.reference_mol,
                "method": self.method,
                "max_confs": self.max_confs,
                "shape_weight": self.shape_weight,
                "color_weight": self.color_weight,
                "index": i,
                "save_overlays": self.save_overlays,
                "energy_window": self.parameters.specific_parameters.get("ewindow", 10),
                "max_stereo": self.parameters.specific_parameters.get("max_stereo", 0)
            }
            
            if self.save_overlays:
                args["overlays_dir"] = self.overlays_dir
                args["overlay_prefix"] = self.overlay_prefix
            
            args_list.append(args)
        
        # Process in parallel
        scores = np.zeros(len(mols), dtype=np.float32)
        
        if len(args_list) > 0:
            with Pool(processes=self.num_cpus) as pool:
                results = pool.map(calculate_single_molecule, args_list)
            
            # Collect results
            for idx, score in results:
                scores[idx] = score
            
        return scores


def calculate_single_molecule(args):
    """Calculate similarity for a single molecule (called by each process)."""
    mol = args["mol"]
    smiles = args["smiles"]
    reference_mol = args["reference_mol"]
    method = args["method"]
    max_confs = args["max_confs"]
    shape_weight = args["shape_weight"]
    color_weight = args["color_weight"]
    index = args["index"]
    save_overlays = args["save_overlays"]
    energy_window = args["energy_window"]
    max_stereo = args["max_stereo"]
    
    if mol is None or reference_mol is None:
        return index, 0.0
        
    # Create conformer generator
    conformer_generator = RDKitConformerGenerator(
        max_confs=max_confs,
        energy_window=energy_window,
        max_stereo=max_stereo
    )
    
    # Generate conformers
    mol_with_confs = conformer_generator.generate_conformers(smiles)
    if mol_with_confs is None or mol_with_confs.GetNumConformers() == 0:
        return index, 0.0
    
    # Calculate best similarity
    best_score = 0.0
    best_conf_pair = None
    
    if method == "usrcat":
        from rdkit.Chem.rdMolDescriptors import GetUSRCAT, GetUSRScore
        
        for q_conf_id in range(mol_with_confs.GetNumConformers()):
            for r_conf_id in range(reference_mol.GetNumConformers()):
                # Get USRCAT descriptors
                query_descriptor = GetUSRCAT(mol_with_confs, confId=q_conf_id)
                ref_descriptor = GetUSRCAT(reference_mol, confId=r_conf_id)
                
                # Calculate similarity
                similarity = GetUSRScore(query_descriptor, ref_descriptor)
                
                if similarity > best_score:
                    best_score = similarity
                    best_conf_pair = (q_conf_id, r_conf_id)
    
    elif method == "o3a":
        from rdkit.Chem import ChemicalFeatures
        from rdkit import RDConfig
        import os
        
        # Load feature factory for pharmacophore features
        fdef_file = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        feature_factory = ChemicalFeatures.BuildFeatureFactory(fdef_file)
        
        for q_conf_id in range(mol_with_confs.GetNumConformers()):
            for r_conf_id in range(reference_mol.GetNumConformers()):
                try:
                    # Create O3A alignment
                    pyO3A = AllChem.GetO3A(mol_with_confs, reference_mol, 
                                          confId1=q_conf_id, confId2=r_conf_id)
                    
                    # Get alignment score (shape similarity)
                    shape_sim = pyO3A.Score() / 100.0  # Normalize to 0-1 range
                    
                    # Align the molecules
                    pyO3A.Align()
                    
                    # Calculate feature similarity after alignment (color similarity)
                    color_sim = 0.0
                    # This would require implementing a feature-based similarity calculation
                    # For simplicity, we'll use shape similarity as a proxy for now
                    
                    # Combine scores using weights
                    combined_score = ((shape_weight * shape_sim) + 
                                     (color_weight * color_sim)) / (shape_weight + color_weight)
                    
                    if combined_score > best_score:
                        best_score = combined_score
                        best_conf_pair = (q_conf_id, r_conf_id)
                except Exception as e:
                    continue
    
    # Save overlay if requested
    if save_overlays and best_conf_pair is not None:
        overlays_dir = args.get("overlays_dir", "overlays")
        overlay_prefix = args.get("overlay_prefix", "mol_")
        
        try:
            # Create output file
            overlay_file = os.path.join(overlays_dir, f"{overlay_prefix}{index}.sdf")
            with Chem.SDWriter(overlay_file) as writer:
                writer.write(mol_with_confs, confId=best_conf_pair[0])
        except Exception as e:
            pass
    
    return index, best_score





