#!/usr/bin/env python3
"""
Test script to debug Roshambo CSV file handling.
This script tests the updated roshambo scoring component with enhanced debugging.
"""

import os
import sys
import tempfile
from pathlib import Path

def create_test_data():
    """Create test SDF files for roshambo testing."""
    print("🧪 Creating test data...")

    # Create test directory
    test_dir = "test_roshambo_csv"
    Path(test_dir).mkdir(exist_ok=True)

    # Create reference SDF file
    ref_sdf = os.path.join(test_dir, "reference.sdf")
    with open(ref_sdf, "w") as f:
        f.write("""reference_mol
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.9517    0.7811   -0.6622 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2847    1.3329   -0.3121 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2365    0.5518    0.3512 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0365   -0.8137    0.6637 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2002   -1.3654    0.3135 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1519   -0.5844   -0.3499 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
$$$$
""")

    # Create dataset SDF file with multiple molecules
    dataset_sdf = os.path.join(test_dir, "dataset.sdf")
    with open(dataset_sdf, "w") as f:
        # Molecule 0 - similar to reference
        f.write("""mol_0
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.9517    0.7811   -0.6622 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2847    1.3329   -0.3121 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2365    0.5518    0.3512 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0365   -0.8137    0.6637 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2002   -1.3654    0.3135 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1519   -0.5844   -0.3499 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
$$$$
""")

        # Molecule 1 - different structure
        f.write("""mol_1
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
M  END
$$$$
""")

        # Molecule 2 - another different structure
        f.write("""mol_2
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
""")

    print(f"✅ Created test files:")
    print(f"   Reference: {ref_sdf}")
    print(f"   Dataset: {dataset_sdf}")

    return test_dir, ref_sdf, dataset_sdf

def test_roshambo_component():
    """Test the roshambo component with CSV debugging."""
    print("\n🧮 Testing Roshambo Component...")

    try:
        # Import the component
        from reinvent_scoring.scoring.component_parameters import ComponentParameters
        from reinvent_scoring.scoring.score_components.roshambo.roshambo_shape_similarity import RoshamboShapeSimilarity

        # Create test data
        test_dir, ref_sdf, dataset_sdf = create_test_data()

        # Configure component parameters
        specific_parameters = {
            "rdbase_path": "/home/protacinvent/Desktop/roshambo/rdkit",  # Update this path as needed
            "conda_env_name": "roshambo",
            "conda_base_path": "/home/protacinvent/.conda",  # Update this path as needed
            "reference_file": ref_sdf,
            "shape_weight": 0.7,
            "color_weight": 0.3,
            "save_overlays": True,
            "overlays_dir": "test_roshambo_overlays",
            "debug": True,  # Enable debug mode
            "gpu_id": 0,
            "n_confs": 0,  # Use existing conformers
            "ignore_hs": True,
            "use_carbon_radii": True
        }

        parameters = ComponentParameters(
            component_type="roshambo_shape_similarity",
            name="test_roshambo",
            weight=1.0,
            specific_parameters=specific_parameters
        )

        # Create component
        print("🔧 Creating Roshambo component...")
        component = RoshamboShapeSimilarity(parameters)

        # Test molecules (SMILES)
        test_smiles = [
            "c1ccccc1",  # benzene - should be similar to reference
            "CCCC",      # butane - different
            "CCO"        # ethanol - different
        ]

        print(f"\n🧪 Testing with {len(test_smiles)} molecules:")
        for i, smi in enumerate(test_smiles):
            print(f"   {i}: {smi}")

        # Calculate scores
        print("\n🚀 Calculating scores...")
        result = component.calculate_score(test_smiles, step=0)

        print(f"\n📊 Results:")
        print(f"   Scores: {result.total_score}")
        print(f"   Score type: {type(result.total_score)}")
        print(f"   Score shape: {result.total_score.shape if hasattr(result.total_score, 'shape') else 'N/A'}")

        # Check if overlays directory was created and contains files
        overlays_dir = specific_parameters["overlays_dir"]
        if os.path.exists(overlays_dir):
            print(f"\n📁 Overlays directory contents:")
            for root, _, files in os.walk(overlays_dir):
                level = root.replace(overlays_dir, '').count(os.sep)
                indent = ' ' * 2 * level
                print(f"   {indent}{os.path.basename(root)}/")
                subindent = ' ' * 2 * (level + 1)
                for file in files:
                    file_path = os.path.join(root, file)
                    file_size = os.path.getsize(file_path) if os.path.exists(file_path) else 0
                    print(f"   {subindent}{file} ({file_size} bytes)")

                    # If it's a CSV file, show its contents
                    if file.endswith('.csv'):
                        print(f"   {subindent}CSV contents:")
                        try:
                            with open(file_path, 'r') as f:
                                lines = f.readlines()[:10]  # Show first 10 lines
                                for line in lines:
                                    print(f"   {subindent}  {line.strip()}")
                                if len(f.readlines()) > 10:
                                    print(f"   {subindent}  ... (truncated)")
                        except Exception as e:
                            print(f"   {subindent}  Error reading file: {e}")

        return True

    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("🔬 Roshambo CSV Debug Test")
    print("=" * 50)

    success = test_roshambo_component()

    if success:
        print("\n✅ Test completed successfully!")
    else:
        print("\n❌ Test failed!")
        sys.exit(1)
