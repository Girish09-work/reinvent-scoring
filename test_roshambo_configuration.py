#!/usr/bin/env python3
"""
Test script to diagnose Roshambo configuration issues.
This script helps identify why the Roshambo component is returning zero scores.
"""

import os
import sys
import json
from pathlib import Path

def test_roshambo_component():
    """Test the Roshambo component with a minimal configuration."""
    print("🔧 Testing Roshambo Component Configuration")
    print("=" * 60)
    
    # Add the reinvent_scoring module to the path
    current_dir = Path(__file__).parent
    sys.path.insert(0, str(current_dir))
    
    try:
        from reinvent_scoring.scoring.component_parameters import ComponentParameters
        from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity
        print("✅ Successfully imported Roshambo component")
    except ImportError as e:
        print(f"❌ Failed to import Roshambo component: {e}")
        return False
    
    # Test configuration - UPDATE THESE PATHS FOR YOUR SYSTEM
    test_config = {
        "component_type": "roshambo_shape_similarity",
        "name": "Test Roshambo",
        "weight": 1.0,
        "specific_parameters": {
            # REQUIRED: Update these paths for your system
            "rdbase_path": "/home/protacinvent/Desktop/roshambo/rdkit",
            "conda_env_name": "roshambo",
            
            # Optional: Update if needed
            "conda_base_path": "/home/protacinvent/.conda",
            "cuda_home_path": "/usr/local/cuda",
            "auto_setup_env": True,
            
            # Reference file - UPDATE THIS PATH
            "reference_file": "/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BTK/BTK_sel.sdf",
            
            # Roshambo parameters
            "shape_weight": 0.6,
            "color_weight": 0.4,
            "n_confs": 0,  # Use existing conformers
            "ignore_hs": True,
            "use_carbon_radii": True,
            "gpu_id": 0,
            
            # Enable debug mode
            "debug": True
        }
    }
    
    print("\n📋 Test Configuration:")
    print(f"  rdbase_path: {test_config['specific_parameters']['rdbase_path']}")
    print(f"  conda_env_name: {test_config['specific_parameters']['conda_env_name']}")
    print(f"  reference_file: {test_config['specific_parameters']['reference_file']}")
    
    # Check if reference file exists
    ref_file = test_config['specific_parameters']['reference_file']
    if os.path.exists(ref_file):
        print(f"  ✅ Reference file exists: {ref_file}")
    else:
        print(f"  ❌ Reference file not found: {ref_file}")
        print("     Please update the 'reference_file' path in this script")
        return False
    
    # Test component initialization
    print("\n🔧 Testing Component Initialization...")
    try:
        params = ComponentParameters(
            component_type=test_config["component_type"],
            name=test_config["name"],
            weight=test_config["weight"],
            specific_parameters=test_config["specific_parameters"]
        )
        
        component = RoshamboShapeSimilarity(params)
        print("✅ Component initialized successfully")
        
    except Exception as e:
        print(f"❌ Component initialization failed: {e}")
        print("\nThis error indicates that required parameters are missing or invalid.")
        print("Please check:")
        print("1. rdbase_path points to a valid RDKit build directory")
        print("2. conda_env_name refers to an existing conda environment")
        print("3. reference_file exists and is accessible")
        return False
    
    # Test scoring with simple molecules
    print("\n🧮 Testing Scoring...")
    test_smiles = [
        "CCO",  # ethanol
        "CCC",  # propane
        "c1ccccc1",  # benzene
    ]
    
    try:
        result = component.calculate_score(test_smiles, step=0)
        scores = result.total_score
        
        print(f"✅ Scoring completed")
        print(f"  Input molecules: {len(test_smiles)}")
        print(f"  Output scores: {scores}")
        print(f"  Non-zero scores: {sum(1 for s in scores if s > 0)}/{len(scores)}")
        
        if all(s == 0.0 for s in scores):
            print("⚠️  All scores are zero - this indicates a problem with Roshambo execution")
            print("   Check the debug output above for error messages")
            return False
        else:
            print("✅ Got non-zero scores - Roshambo is working correctly")
            return True
            
    except Exception as e:
        print(f"❌ Scoring failed: {e}")
        return False

def main():
    """Main function to run the test."""
    print("🧪 Roshambo Configuration Test")
    print("This script tests the Roshambo component configuration")
    print("and helps diagnose why scores might be returning zero.\n")
    
    success = test_roshambo_component()
    
    print("\n" + "=" * 60)
    if success:
        print("✅ Test completed successfully!")
        print("The Roshambo component is configured correctly.")
    else:
        print("❌ Test failed!")
        print("Please check the error messages above and:")
        print("1. Update the paths in this script for your system")
        print("2. Ensure the conda environment 'roshambo' exists")
        print("3. Verify the RDKit build path is correct")
        print("4. Check that the reference file exists")
        print("\nFor more help, see the configuration examples in:")
        print("- examples/roshambo_configuration_examples.json")
        print("- ROSHAMBO_MANDATORY_PARAMETERS.md")

if __name__ == "__main__":
    main()
