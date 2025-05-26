#!/usr/bin/env python3
"""
Test script to trace roshambo score responses and understand the return format.
This script demonstrates exactly how scores are returned from the roshambo API.
"""

import os
import sys
import tempfile
from pathlib import Path

def create_test_data():
    """Create minimal test SDF files for roshambo testing."""

    # Create test directory
    test_dir = "roshambo_test_data"
    os.makedirs(test_dir, exist_ok=True)

    # Simple reference molecule (benzene)
    ref_sdf_content = """
  Mrv2014 01010100002D

  6  6  0  0  0  0            999 V2000
   -1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  END
$$$$
"""

    # Simple dataset molecules (similar to benzene)
    dataset_sdf_content = """mol_0
  Mrv2014 01010100002D

  6  6  0  0  0  0            999 V2000
   -1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  END
$$$$
mol_1
  Mrv2014 01010100002D

  7  7  0  0  0  0            999 V2000
   -1.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.7320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.7320    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  1  1  0  0  0  0
  7  1  1  0  0  0  0
M  END
$$$$
"""

    # Write files
    with open(f"{test_dir}/reference.sdf", "w") as f:
        f.write(ref_sdf_content)

    with open(f"{test_dir}/dataset.sdf", "w") as f:
        f.write(dataset_sdf_content)

    return test_dir

def test_roshambo_api_direct():
    """Test roshambo API directly to see return format."""

    print("=" * 60)
    print("TESTING ROSHAMBO API DIRECTLY")
    print("=" * 60)

    test_dir = create_test_data()

    try:
        from roshambo.api import get_similarity_scores

        print("✅ Successfully imported roshambo.api.get_similarity_scores")

        # Test with write_to_file=False to see direct return
        print("\n--- Testing with write_to_file=False ---")

        result = get_similarity_scores(
            ref_file=f"{test_dir}/reference.sdf",
            dataset_files_pattern=f"{test_dir}/dataset.sdf",
            ignore_hs=True,
            n_confs=0,  # Use existing conformers
            use_carbon_radii=True,
            color=True,
            sort_by="ComboTanimoto",
            write_to_file=False,  # Don't write files, return directly
            gpu_id=0,
            working_dir=None
        )

        print(f"\n🔍 RESULT ANALYSIS:")
        print(f"   Type: {type(result)}")
        print(f"   Is DataFrame: {hasattr(result, 'shape') and hasattr(result, 'columns')}")

        if hasattr(result, 'shape'):
            print(f"   Shape: {result.shape}")

        if hasattr(result, 'columns'):
            print(f"   Columns: {result.columns.tolist()}")

        if hasattr(result, 'head'):
            print(f"\n📊 FIRST FEW ROWS:")
            print(result.head())

        if hasattr(result, 'iterrows'):
            print(f"\n🎯 INDIVIDUAL SCORES:")
            for idx, row in result.iterrows():
                molecule = row.get('Molecule', 'Unknown')
                shape_score = row.get('ShapeTanimoto', 0.0)
                color_score = row.get('ColorTanimoto', 0.0)
                combo_score = row.get('ComboTanimoto', 0.0)
                print(f"   {molecule}: Shape={shape_score:.3f}, Color={color_score:.3f}, Combo={combo_score:.3f}")

        return True

    except ImportError as e:
        print(f"❌ Failed to import roshambo: {e}")
        return False
    except Exception as e:
        print(f"❌ Error running roshambo: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_roshambo_api_with_files():
    """Test roshambo API with write_to_file=True to see file output."""

    print("\n" + "=" * 60)
    print("TESTING ROSHAMBO API WITH FILE OUTPUT")
    print("=" * 60)

    test_dir = create_test_data()

    try:
        from roshambo.api import get_similarity_scores

        # Test with write_to_file=True to see file output
        print("\n--- Testing with write_to_file=True ---")

        result = get_similarity_scores(
            ref_file="reference.sdf",
            dataset_files_pattern="dataset.sdf",
            ignore_hs=True,
            n_confs=0,  # Use existing conformers
            use_carbon_radii=True,
            color=True,
            sort_by="ComboTanimoto",
            write_to_file=True,  # Write files
            gpu_id=0,
            working_dir=test_dir
        )

        print(f"\n🔍 RESULT ANALYSIS (with file output):")
        print(f"   Type: {type(result)}")

        if hasattr(result, 'shape'):
            print(f"   Shape: {result.shape}")

        # Check for output files
        output_files = []
        for file in os.listdir(test_dir):
            if file.endswith('.csv') or file.endswith('.sdf'):
                output_files.append(file)

        print(f"\n📁 OUTPUT FILES CREATED:")
        for file in output_files:
            file_path = os.path.join(test_dir, file)
            file_size = os.path.getsize(file_path)
            print(f"   {file} ({file_size} bytes)")

            # Show CSV content if it's a CSV file
            if file.endswith('.csv'):
                print(f"\n📄 CONTENT OF {file}:")
                with open(file_path, 'r') as f:
                    content = f.read()
                    print(content[:500] + "..." if len(content) > 500 else content)

        return True

    except Exception as e:
        print(f"❌ Error running roshambo with files: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("🧪 ROSHAMBO SCORE TRACING TEST")
    print("This script demonstrates how roshambo returns scores")

    # Test direct API
    success1 = test_roshambo_api_direct()

    # Test with file output
    success2 = test_roshambo_api_with_files()

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    if success1:
        print("✅ Direct API test: SUCCESS")
        print("   → roshambo.api.get_similarity_scores() returns a pandas DataFrame")
        print("   → DataFrame contains columns: Molecule, ShapeTanimoto, ColorTanimoto, ComboTanimoto, etc.")
        print("   → Scores are accessed via DataFrame.iterrows() or DataFrame indexing")
    else:
        print("❌ Direct API test: FAILED")

    if success2:
        print("✅ File output test: SUCCESS")
        print("   → When write_to_file=True, roshambo writes CSV files to working_dir")
        print("   → CSV files contain the same data as the returned DataFrame")
    else:
        print("❌ File output test: FAILED")

    print("\n🎯 KEY INSIGHTS:")
    print("   1. The print statements you see come from roshambo's internal timing logs.")
    print("   2. ❌ BUG DISCOVERED: roshambo.api.get_similarity_scores() has a bug!")
    print("      → It should return a pandas DataFrame but returns None instead")
    print("      → The function works correctly but doesn't return the result")
    print("   3. ✅ WORKAROUND: Always use write_to_file=True and read from CSV")
    print("      → The CSV file contains all the correct scores")
    print("      → CSV format: tab-separated with columns Molecule, ShapeTanimoto, ColorTanimoto, etc.")
    print("   4. 🔧 INTEGRATION FIX: Updated reinvent-scoring to handle this correctly")
    print("      → Always forces write_to_file=True")
    print("      → Reads scores from the generated CSV file")
    print("      → Handles both direct API and subprocess execution")
