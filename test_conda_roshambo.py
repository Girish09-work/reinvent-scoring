#!/usr/bin/env python3
"""
Test script to verify conda path and roshambo availability.
"""

import os
import subprocess
import sys

def test_conda_paths():
    """Test different conda paths to find the correct one."""
    print("🔍 Testing conda paths...")
    
    possible_paths = [
        "/home/protacinvent/.conda/bin/conda",
        "/opt/conda/bin/conda", 
        "/usr/local/bin/conda",
        "conda"  # System PATH
    ]
    
    for conda_path in possible_paths:
        try:
            print(f"\n   Testing: {conda_path}")
            result = subprocess.run([conda_path, "--version"], capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                print(f"   ✅ Success: {result.stdout.strip()}")
                
                # Test listing environments
                env_result = subprocess.run([conda_path, "env", "list"], capture_output=True, text=True, timeout=10)
                if env_result.returncode == 0:
                    print(f"   📋 Available environments:")
                    for line in env_result.stdout.split('\n'):
                        if 'roshambo' in line.lower():
                            print(f"      🎯 {line.strip()}")
                        elif line.strip() and not line.startswith('#'):
                            print(f"         {line.strip()}")
                
                return conda_path
            else:
                print(f"   ❌ Failed: {result.stderr.strip()}")
        except FileNotFoundError:
            print(f"   ❌ Not found: {conda_path}")
        except subprocess.TimeoutExpired:
            print(f"   ❌ Timeout: {conda_path}")
        except Exception as e:
            print(f"   ❌ Error: {e}")
    
    return None

def test_roshambo_in_current_env():
    """Test if roshambo is available in the current environment."""
    print("\n🧪 Testing roshambo in current environment...")
    
    try:
        from roshambo.api import get_similarity_scores
        print("   ✅ roshambo.api.get_similarity_scores imported successfully")
        return True
    except ImportError as e:
        print(f"   ❌ Import failed: {e}")
        return False

def test_roshambo_in_conda_env(conda_path, env_name="roshambo"):
    """Test if roshambo is available in a specific conda environment."""
    print(f"\n🧪 Testing roshambo in conda environment '{env_name}'...")
    
    if not conda_path:
        print("   ❌ No valid conda path found")
        return False
    
    try:
        # Test if the environment exists
        env_result = subprocess.run([conda_path, "env", "list"], capture_output=True, text=True, timeout=10)
        if env_name not in env_result.stdout:
            print(f"   ❌ Environment '{env_name}' not found")
            return False
        
        print(f"   ✅ Environment '{env_name}' exists")
        
        # Test importing roshambo in that environment
        test_script = f'''
import sys
try:
    from roshambo.api import get_similarity_scores
    print("SUCCESS: roshambo imported")
    print(f"Python path: {{sys.executable}}")
    print(f"Python version: {{sys.version}}")
except ImportError as e:
    print(f"FAILED: {{e}}")
    sys.exit(1)
'''
        
        result = subprocess.run([conda_path, "run", "-n", env_name, "python", "-c", test_script], 
                              capture_output=True, text=True, timeout=30)
        
        print(f"   Return code: {result.returncode}")
        print(f"   Output: {result.stdout.strip()}")
        if result.stderr:
            print(f"   Error: {result.stderr.strip()}")
        
        return result.returncode == 0
        
    except Exception as e:
        print(f"   ❌ Error testing conda environment: {e}")
        return False

def test_environment_variables():
    """Test current environment variables."""
    print("\n🔧 Current environment variables:")
    
    env_vars = [
        "RDBASE", "RDKIT_LIB_DIR", "RDKIT_INCLUDE_DIR", "RDKIT_DATA_DIR",
        "PYTHONPATH", "LD_LIBRARY_PATH", "CUDA_HOME", "PATH"
    ]
    
    for var in env_vars:
        value = os.environ.get(var, "Not set")
        if len(value) > 100:
            value = value[:100] + "..."
        print(f"   {var}: {value}")

def main():
    print("🔬 Conda and Roshambo Environment Test")
    print("=" * 50)
    
    # Test conda paths
    conda_path = test_conda_paths()
    
    # Test roshambo in current environment
    current_env_works = test_roshambo_in_current_env()
    
    # Test roshambo in conda environment
    conda_env_works = test_roshambo_in_conda_env(conda_path)
    
    # Show environment variables
    test_environment_variables()
    
    print("\n📊 Summary:")
    print(f"   Valid conda path: {conda_path or 'None found'}")
    print(f"   Roshambo in current env: {'✅' if current_env_works else '❌'}")
    print(f"   Roshambo in conda env: {'✅' if conda_env_works else '❌'}")
    
    if conda_path and conda_env_works:
        print("\n✅ Recommended configuration:")
        print(f"   conda_base_path: {os.path.dirname(os.path.dirname(conda_path))}")
        print(f"   conda_env_name: roshambo")
    elif current_env_works:
        print("\n✅ Use direct execution (no subprocess needed)")
    else:
        print("\n❌ Roshambo setup issues detected")
        print("   Please check:")
        print("   1. Conda installation and PATH")
        print("   2. Roshambo environment exists")
        print("   3. Roshambo is installed in the environment")

if __name__ == "__main__":
    main()
