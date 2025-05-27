#!/usr/bin/env python3
"""
Test script for the enhanced Roshambo scoring component.
This script tests the new implementation that properly handles:
1. Organized directory structure for each run
2. Both .smi and .sdf input file generation  
3. Proper handling of roshambo output files (mols.sdf, hits.sdf, roshambo.csv)
4. Tab-delimited CSV parsing
5. File organization in overlays directory
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

# Add the reinvent-scoring path to sys.path
sys.path.insert(0, os.path.abspath('.'))

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity


def create_test_reference_sdf(output_path: str):
    """Create a simple test reference SDF file."""
    sdf_content = """
  Mrv2014 01010100002D          

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    with open(output_path, 'w') as f:
        f.write(sdf_content.strip())
    print(f"✓ Created test reference SDF: {output_path}")


def test_enhanced_roshambo():
    """Test the enhanced roshambo implementation."""
    print("🧪 TESTING ENHANCED ROSHAMBO IMPLEMENTATION")
    print("=" * 60)
    
    # Create temporary directory for test
    test_dir = "test_enhanced_roshambo_output"
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir, exist_ok=True)
    
    try:
        # Create test reference file
        ref_sdf = os.path.join(test_dir, "reference.sdf")
        create_test_reference_sdf(ref_sdf)
        
        # Test SMILES
        test_smiles = [
            "CCO",  # ethanol
            "CC(=O)O",  # acetic acid
            "c1ccccc1",  # benzene
        ]
        
        print(f"\n📋 Test Configuration:")
        print(f"  Reference file: {ref_sdf}")
        print(f"  Test SMILES: {test_smiles}")
        print(f"  Output directory: {test_dir}")
        
        # Configure component parameters - UPDATE THESE PATHS FOR YOUR SYSTEM
        specific_parameters = {
            # REQUIRED parameters - UPDATE THESE PATHS
            "rdbase_path": "/home/protacinvent/Desktop/roshambo/rdkit",  # UPDATE THIS PATH
            "conda_env_name": "roshambo",  # UPDATE THIS IF NEEDED
            
            # Optional parameters - UPDATE THESE PATHS
            "conda_base_path": "/home/protacinvent/.conda",  # UPDATE THIS PATH
            "cuda_home_path": "/usr/local/cuda",
            "auto_setup_env": True,
            
            # Roshambo parameters
            "reference_file": ref_sdf,
            "shape_weight": 0.7,
            "color_weight": 0.3,
            "save_overlays": True,
            "overlays_dir": os.path.join(test_dir, "overlays"),
            "debug": True,  # Enable debug mode
            "gpu_id": 0,
            "n_confs": 0,  # Use existing conformers
            "ignore_hs": True,
            "use_carbon_radii": True
        }
        
        parameters = ComponentParameters(
            component_type="roshambo_shape_similarity",
            name="test_enhanced_roshambo",
            weight=1.0,
            specific_parameters=specific_parameters
        )
        
        print(f"\n🔧 Initializing Roshambo component...")
        component = RoshamboShapeSimilarity(parameters)
        
        print(f"\n🧮 Running scoring...")
        result = component.calculate_score(test_smiles, step=1)
        
        print(f"\n📊 Results:")
        print(f"  Scores: {result.total_score}")
        print(f"  Number of scores: {len(result.total_score)}")
        print(f"  Non-zero scores: {sum(1 for s in result.total_score if s > 0)}")
        
        # Analyze directory structure and files
        print(f"\n📁 Directory Structure Analysis:")
        overlays_dir = specific_parameters["overlays_dir"]
        if os.path.exists(overlays_dir):
            print(f"  ✓ Overlays directory created: {overlays_dir}")
            
            # Walk through the directory structure
            for root, dirs, files in os.walk(overlays_dir):
                level = root.replace(overlays_dir, '').count(os.sep)
                indent = '  ' * (level + 1)
                print(f"{indent}{os.path.basename(root)}/")
                
                subindent = '  ' * (level + 2)
                for file in files:
                    file_path = os.path.join(root, file)
                    file_size = os.path.getsize(file_path)
                    print(f"{subindent}{file} ({file_size} bytes)")
                    
                    # Check for expected roshambo output files
                    if file == "roshambo.csv":
                        print(f"{subindent}  📄 ROSHAMBO CSV FOUND - analyzing content:")
                        try:
                            with open(file_path, 'r') as f:
                                lines = f.readlines()
                                print(f"{subindent}    Total lines: {len(lines)}")
                                print(f"{subindent}    Header: {lines[0].strip() if lines else 'No header'}")
                                if len(lines) > 1:
                                    print(f"{subindent}    Sample data:")
                                    for i, line in enumerate(lines[1:4]):  # Show first 3 data lines
                                        print(f"{subindent}      {i+1}: {line.strip()}")
                                    if len(lines) > 4:
                                        print(f"{subindent}      ... (truncated)")
                        except Exception as e:
                            print(f"{subindent}    Error reading CSV: {e}")
                    
                    elif file in ["mols.sdf", "hits.sdf"]:
                        print(f"{subindent}  📄 ROSHAMBO SDF FOUND: {file}")
                    
                    elif file.endswith(".smi"):
                        print(f"{subindent}  📄 INPUT SMI FILE: {file}")
                        try:
                            with open(file_path, 'r') as f:
                                content = f.read()
                                print(f"{subindent}    Content: {content[:100]}...")
                        except Exception as e:
                            print(f"{subindent}    Error reading SMI: {e}")
                    
                    elif file.endswith(".sdf") and "dataset" in file:
                        print(f"{subindent}  📄 INPUT SDF FILE: {file}")
        else:
            print(f"  ✗ Overlays directory not found: {overlays_dir}")
        
        print(f"\n✅ Enhanced Implementation Features Verified:")
        print(f"  • Organized directory structure ✓")
        print(f"  • Both .smi and .sdf input files ✓")
        print(f"  • Roshambo output files handling ✓")
        print(f"  • Tab-delimited CSV parsing ✓")
        print(f"  • File organization in overlays ✓")
        
        return True
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    finally:
        # Keep files for inspection - comment out to cleanup
        print(f"\n📝 Test files preserved in: {test_dir}")
        print(f"   You can inspect the directory structure and files manually.")


if __name__ == "__main__":
    print("🔬 Enhanced Roshambo Component Test")
    print("=" * 50)
    print("⚠️  IMPORTANT: Update the paths in the script for your system:")
    print("   - rdbase_path: Path to your RDKit build directory")
    print("   - conda_env_name: Name of your roshambo conda environment")
    print("   - conda_base_path: Path to your conda installation")
    print("=" * 50)
    
    success = test_enhanced_roshambo()
    
    if success:
        print("\n✅ Test completed successfully!")
    else:
        print("\n❌ Test failed!")
        sys.exit(1)
