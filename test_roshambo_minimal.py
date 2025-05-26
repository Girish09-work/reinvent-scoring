#!/usr/bin/env python3
"""
Minimal test to reproduce the roshambo error.
"""

import os
import tempfile
from pathlib import Path

def create_test_sdf():
    """Create a simple test SDF file with a few molecules."""
    # Create a temporary SDF file with simple molecules
    temp_dir = "/tmp/roshambo_test"
    Path(temp_dir).mkdir(exist_ok=True)
    
    sdf_file = os.path.join(temp_dir, "test_molecules.sdf")
    
    # Simple SDF content with a few molecules
    sdf_content = """mol_0
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    2.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  6  1  1  0
M  END
$$$$
mol_1
     RDKit          3D

  8  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  6  7  1  0
  7  8  1  0
  8  1  1  0
M  END
$$$$
"""
    
    with open(sdf_file, 'w') as f:
        f.write(sdf_content)
    
    print(f"Created test SDF file: {sdf_file}")
    return sdf_file

def test_roshambo_call():
    """Test the exact roshambo call that's failing."""
    
    # Create test dataset
    dataset_file = create_test_sdf()
    
    # Use the first reference file that we know is valid
    ref_file = "/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BRD9/protac.sdf"
    
    print(f"Testing roshambo with:")
    print(f"  ref_file: {ref_file}")
    print(f"  dataset_file: {dataset_file}")
    
    try:
        from roshambo.api import get_similarity_scores
        
        print("Calling roshambo.api.get_similarity_scores...")
        
        results = get_similarity_scores(
            ref_file=ref_file,
            dataset_files_pattern=dataset_file,
            ignore_hs=True,
            n_confs=0,  # Don't generate conformers since we have SDF
            use_carbon_radii=True,
            color=True,
            sort_by="ComboTanimoto",
            write_to_file=False,  # Don't write files for this test
            gpu_id=0,
            working_dir=None
        )
        
        print("✅ Roshambo call succeeded!")
        print(f"Results type: {type(results)}")
        if hasattr(results, 'scores'):
            print(f"Number of scores: {len(results.scores)}")
        
        return True
        
    except Exception as e:
        print(f"❌ Roshambo call failed: {e}")
        import traceback
        print("Full traceback:")
        traceback.print_exc()
        return False

def main():
    """Run the test."""
    print("🧪 Testing minimal roshambo call...")
    
    # Test with simple SDF
    success = test_roshambo_call()
    
    if success:
        print("\n✅ Test passed! The issue might be with the specific dataset being generated.")
    else:
        print("\n❌ Test failed! There's a fundamental issue with roshambo setup.")
    
    return success

if __name__ == "__main__":
    import sys
    success = main()
    sys.exit(0 if success else 1)
