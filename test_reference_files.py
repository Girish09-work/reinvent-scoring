#!/usr/bin/env python3
"""
Test script to validate reference SDF files for roshambo.
"""

import os
import sys

def test_sdf_file(filepath):
    """Test if an SDF file is valid and can be read by RDKit."""
    print(f"\n=== Testing {filepath} ===")
    
    if not os.path.exists(filepath):
        print(f"❌ File does not exist: {filepath}")
        return False
    
    # Check file size
    file_size = os.path.getsize(filepath)
    print(f"📁 File size: {file_size} bytes")
    
    if file_size == 0:
        print("❌ File is empty")
        return False
    
    # Read first few lines
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            print(f"📄 Total lines: {len(lines)}")
            print(f"📄 First 5 lines:")
            for i, line in enumerate(lines[:5]):
                print(f"  {i+1}: {line.rstrip()}")
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return False
    
    # Test with RDKit
    try:
        from rdkit import Chem
        suppl = Chem.SDMolSupplier(filepath)
        mol_count = 0
        valid_mols = 0
        
        for mol in suppl:
            mol_count += 1
            if mol is not None:
                valid_mols += 1
                if valid_mols == 1:  # Show info for first valid molecule
                    print(f"🧪 First molecule:")
                    print(f"   Atoms: {mol.GetNumAtoms()}")
                    print(f"   Bonds: {mol.GetNumBonds()}")
                    print(f"   SMILES: {Chem.MolToSmiles(mol)}")
        
        print(f"🧪 Total molecules in file: {mol_count}")
        print(f"🧪 Valid molecules: {valid_mols}")
        
        if valid_mols > 0:
            print("✅ File is valid and contains readable molecules")
            return True
        else:
            print("❌ File contains no valid molecules")
            return False
            
    except Exception as e:
        print(f"❌ Error reading with RDKit: {e}")
        return False

def main():
    """Test all reference files."""
    reference_files = [
        "/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BRD9/protac.sdf",
        "/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BRD9/warhead.sdf",
        "/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BTK/BTK_sel.sdf"
    ]
    
    print("🔍 Testing reference SDF files for roshambo compatibility...")
    
    all_valid = True
    for ref_file in reference_files:
        is_valid = test_sdf_file(ref_file)
        if not is_valid:
            all_valid = False
    
    print(f"\n{'='*50}")
    if all_valid:
        print("✅ All reference files are valid!")
    else:
        print("❌ Some reference files have issues!")
    
    return all_valid

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
