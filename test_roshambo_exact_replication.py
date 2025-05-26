#!/usr/bin/env python3
"""
Exact replication of working roshambo setup to identify the issue.
"""

import os
import tempfile
import shutil
from pathlib import Path

def setup_test_environment():
    """Set up the exact environment that works."""
    # Create test directory structure
    test_dir = "/tmp/roshambo_exact_test"
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    
    Path(test_dir).mkdir(parents=True, exist_ok=True)
    
    # Copy reference file to test directory
    ref_source = "/home/protacinvent/Desktop/Getting Started/protac-invent/Protac-invent/data/protac/BRD9/protac.sdf"
    ref_dest = os.path.join(test_dir, "query.sdf")
    shutil.copy2(ref_source, ref_dest)
    
    # Create a simple dataset SDF (using first molecule from reference)
    dataset_dest = os.path.join(test_dir, "dataset.sdf")
    
    # Extract first molecule from reference to create dataset
    try:
        from rdkit import Chem
        suppl = Chem.SDMolSupplier(ref_source)
        mol = next(suppl)
        if mol:
            writer = Chem.SDWriter(dataset_dest)
            writer.write(mol)
            writer.close()
            print(f"✅ Created dataset SDF: {dataset_dest}")
        else:
            print("❌ Could not read molecule from reference")
            return None
    except Exception as e:
        print(f"❌ Error creating dataset: {e}")
        return None
    
    return test_dir

def test_direct_roshambo():
    """Test roshambo directly with the exact working parameters."""
    
    test_dir = setup_test_environment()
    if not test_dir:
        return False
    
    # Change to test directory (like in working example)
    original_cwd = os.getcwd()
    os.chdir(test_dir)
    
    try:
        print(f"🔍 Testing roshambo from directory: {os.getcwd()}")
        print(f"Files in directory: {os.listdir('.')}")
        
        # Import roshambo
        from roshambo.api import get_similarity_scores
        
        # Call with exact working parameters
        print("📞 Calling roshambo with working parameters...")
        
        results = get_similarity_scores(
            ref_file="query.sdf",
            dataset_files_pattern="dataset.sdf",
            ignore_hs=True,
            n_confs=0,
            use_carbon_radii=True,
            color=True,
            sort_by="ComboTanimoto",
            write_to_file=True,
            gpu_id=0,
            working_dir="inpdata",
        )
        
        print("✅ Direct roshambo call succeeded!")
        print(f"Results type: {type(results)}")
        print(f"Results: {results}")
        
        # Check for output files
        if os.path.exists("inpdata"):
            files = os.listdir("inpdata")
            print(f"Output files created: {files}")
            
            # Look for CSV files
            csv_files = [f for f in files if f.endswith('.csv')]
            if csv_files:
                print(f"CSV files found: {csv_files}")
                
                # Read first CSV
                csv_path = os.path.join("inpdata", csv_files[0])
                try:
                    import pandas as pd
                    df = pd.read_csv(csv_path)
                    print(f"CSV content:\n{df}")
                except Exception as e:
                    print(f"Error reading CSV: {e}")
        
        return True
        
    except Exception as e:
        print(f"❌ Direct roshambo call failed: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        os.chdir(original_cwd)

def test_subprocess_roshambo():
    """Test roshambo via subprocess (like our integration does)."""
    
    test_dir = setup_test_environment()
    if not test_dir:
        return False
    
    # Create the exact script that our integration creates
    script_content = f'''
import os
import sys
import json

# Set up environment variables (exact copy from our integration)
os.environ["RDBASE"] = "/home/protacinvent/Desktop/roshambo/rdkit"
os.environ["RDKIT_LIB_DIR"] = "/home/protacinvent/Desktop/roshambo/rdkit/lib"
os.environ["RDKIT_INCLUDE_DIR"] = "/home/protacinvent/Desktop/roshambo/rdkit/Code"
os.environ["RDKIT_DATA_DIR"] = "/home/protacinvent/Desktop/roshambo/rdkit/Data"
os.environ["CUDA_HOME"] = "/usr/local/cuda"

# Set up PYTHONPATH
current_pythonpath = os.environ.get("PYTHONPATH", "")
new_pythonpath = "/home/protacinvent/Desktop/roshambo/rdkit"
if current_pythonpath:
    new_pythonpath += ":" + current_pythonpath
os.environ["PYTHONPATH"] = new_pythonpath

# Set up LD_LIBRARY_PATH
current_ld_path = os.environ.get("LD_LIBRARY_PATH", "")
new_ld_path = "/home/protacinvent/Desktop/roshambo/rdkit/lib"
if current_ld_path:
    new_ld_path += ":" + current_ld_path
os.environ["LD_LIBRARY_PATH"] = new_ld_path

# Add RDKit to Python path
sys.path.insert(0, "/home/protacinvent/Desktop/roshambo/rdkit")

# Change to test directory
os.chdir("{test_dir}")

try:
    from roshambo.api import get_similarity_scores
    
    print("Environment setup complete, calling roshambo...")
    print(f"Current directory: {{os.getcwd()}}")
    print(f"Files available: {{os.listdir('.')}}")
    
    result = get_similarity_scores(
        ref_file="query.sdf",
        dataset_files_pattern="dataset.sdf",
        ignore_hs=True,
        n_confs=0,
        use_carbon_radii=True,
        color=True,
        sort_by="ComboTanimoto",
        write_to_file=True,
        gpu_id=0,
        working_dir="inpdata"
    )
    
    print("SUCCESS: Subprocess roshambo call worked!")
    print(f"Result: {{result}}")
    
    # Check output files
    if os.path.exists("inpdata"):
        files = os.listdir("inpdata")
        print(f"Output files: {{files}}")
    
except Exception as e:
    print(f"ERROR: {{e}}")
    import traceback
    traceback.print_exc()
'''
    
    script_path = os.path.join(test_dir, "test_subprocess.py")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Execute via subprocess (like our integration)
    import subprocess
    
    print("🔍 Testing roshambo via subprocess...")
    
    # Try the exact command our integration uses
    commands_to_try = [
        f"conda run -n roshambo python {script_path}",
        f"/home/protacinvent/.conda/envs/roshambo/bin/python {script_path}",
        f"bash -c 'source /home/protacinvent/Desktop/Getting\\ Started/protac-invent/Protac-invent/setup_roshambo_env.sh && python {script_path}'"
    ]
    
    for i, cmd in enumerate(commands_to_try):
        print(f"\n📞 Trying command {i+1}: {cmd}")
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=test_dir)
            
            print(f"Return code: {result.returncode}")
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")
            
            if result.returncode == 0:
                print(f"✅ Command {i+1} succeeded!")
                return True
            else:
                print(f"❌ Command {i+1} failed")
                
        except Exception as e:
            print(f"❌ Command {i+1} exception: {e}")
    
    return False

def main():
    """Run comprehensive tests."""
    print("🧪 Comprehensive Roshambo Testing")
    print("=" * 50)
    
    # Test 1: Direct call
    print("\n🔬 TEST 1: Direct roshambo call")
    direct_success = test_direct_roshambo()
    
    # Test 2: Subprocess call
    print("\n🔬 TEST 2: Subprocess roshambo call")
    subprocess_success = test_subprocess_roshambo()
    
    # Summary
    print("\n" + "=" * 50)
    print("📊 SUMMARY:")
    print(f"Direct call: {'✅ SUCCESS' if direct_success else '❌ FAILED'}")
    print(f"Subprocess call: {'✅ SUCCESS' if subprocess_success else '❌ FAILED'}")
    
    if direct_success and not subprocess_success:
        print("\n💡 DIAGNOSIS: Environment setup issue in subprocess")
        print("   - Roshambo works directly")
        print("   - Subprocess environment is not set up correctly")
        print("   - Need to fix environment variable setup")
    elif not direct_success:
        print("\n💡 DIAGNOSIS: Fundamental roshambo issue")
        print("   - Check roshambo installation")
        print("   - Check file permissions")
        print("   - Check working directory setup")
    else:
        print("\n🎉 Both tests passed - integration should work!")

if __name__ == "__main__":
    main()
