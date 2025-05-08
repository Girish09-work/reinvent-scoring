from typing import List
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_components.rdkit_shape.rdkit_conformer_generator import RDKitConformerGenerator


class RDKitShapeSimilarity(BaseScoreComponent):
    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self.shape_weight = self.parameters.specific_parameters.get("shape_weight", 0.5)
        self.color_weight = self.parameters.specific_parameters.get("color_weight", 0.5)
        self.reference_file = self.parameters.specific_parameters.get("reference_file", "")
        self.method = self.parameters.specific_parameters.get("method", "usrcat")
        self.max_confs = self.parameters.specific_parameters.get("max_confs", 50)
        
        # Initialize conformer generator
        self.conformer_generator = RDKitConformerGenerator(
            max_confs=self.max_confs,
            energy_window=self.parameters.specific_parameters.get("ewindow", 10),
            max_stereo=self.parameters.specific_parameters.get("max_stereo", 0)
        )
        
        # Load reference molecule
        self.reference_mol = self._load_reference_molecule()
        
    def _load_reference_molecule(self):
        """Load the reference molecule from file and generate conformers."""
        if not self.reference_file:
            return None
        
        # Load reference molecule
        ref_mol = None
        if self.reference_file.endswith('.sdf'):
            suppl = Chem.SDMolSupplier(self.reference_file)
            if suppl and len(suppl) > 0:
                ref_mol = suppl[0]
    
        # Generate conformers for reference molecule if needed
        if ref_mol and ref_mol.GetNumConformers() == 0:
            ref_smiles = Chem.MolToSmiles(ref_mol)
            ref_mol = self.conformer_generator.generate_conformers(ref_smiles)
        
        return ref_mol
        
    def calculate_score(self, molecules: List) -> ComponentSummary:
        """Calculate shape similarity scores for a list of molecules."""
        # Convert SMILES to RDKit molecules
        rdkit_mols = [Chem.MolFromSmiles(smiles) for smiles in molecules]
        
        # Calculate scores
        scores = self._calculate_shape_scores(rdkit_mols)
        
        # Create and return component summary
        score_summary = ComponentSummary(total_score=scores, parameters=self.parameters)
        return score_summary
        
    def _calculate_shape_scores(self, mols):
        """Calculate shape similarity scores for a list of RDKit molecules."""
        scores = []
        
        for mol in mols:
            if mol is None or self.reference_mol is None:
                scores.append(0.0)
                continue
                
            # Generate conformers using the conformer generator
            mol_with_confs = self.conformer_generator.generate_conformers(Chem.MolToSmiles(mol))
            if mol_with_confs is None or mol_with_confs.GetNumConformers() == 0:
                scores.append(0.0)
                continue
            
            # Calculate best similarity score across all conformer pairs
            best_score = self._calculate_best_similarity(mol_with_confs)
            scores.append(best_score)
            
        return np.array(scores, dtype=np.float32)
        
    def _calculate_best_similarity(self, mol):
        """Calculate the best similarity score across all conformer pairs."""
        best_score = 0.0
        
        if self.method == "usrcat":
            best_score = self._calculate_usrcat_similarity(mol)
        elif self.method == "o3a":
            best_score = self._calculate_o3a_similarity(mol)
            
        return best_score
        
    def _calculate_usrcat_similarity(self, mol):
        """Calculate similarity using USRCAT method."""
        from rdkit.Chem.rdMolDescriptors import GetUSRCAT, GetUSRScore
        
        best_score = 0.0
        
        for q_conf_id in range(mol.GetNumConformers()):
            for r_conf_id in range(self.reference_mol.GetNumConformers()):
                # Get USRCAT descriptors
                query_descriptor = GetUSRCAT(mol, confId=q_conf_id)
                ref_descriptor = GetUSRCAT(self.reference_mol, confId=r_conf_id)
                
                # Calculate similarity
                similarity = GetUSRScore(query_descriptor, ref_descriptor)
                
                if similarity > best_score:
                    best_score = similarity
                    
        return best_score
        
    def _calculate_o3a_similarity(self, mol):
        """Calculate similarity using O3A method."""
        best_score = 0.0
        
        for q_conf_id in range(mol.GetNumConformers()):
            for r_conf_id in range(self.reference_mol.GetNumConformers()):
                # Create O3A alignment
                pyO3A = AllChem.GetO3A(mol, self.reference_mol, 
                                      confId1=q_conf_id, confId2=r_conf_id)
                
                # Get alignment score (shape similarity)
                shape_sim = pyO3A.Score() / 100.0  # Normalize to 0-1 range
                
                if shape_sim > best_score:
                    best_score = shape_sim
                    
        return best_score

