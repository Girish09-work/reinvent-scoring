#!/bin/bash

# =============================================================================
# Simple Roshambo Test Script
# =============================================================================
# This script tests if Roshambo is working correctly in your environment
# Usage: bash test_roshambo_simple.sh
# =============================================================================

echo "Testing Roshambo Environment..."
echo "================================"

# Since you're already in the roshambo environment, let's test directly
echo "Current environment: $CONDA_DEFAULT_ENV"
echo "Python executable: $(which python)"
echo ""

# Test basic imports
echo "1. Testing Python imports:"
python -c "
import sys
print('Python version:', sys.version)
print('Python path entries:')
for i, path in enumerate(sys.path[:5]):  # Show first 5 entries
    print(f'  {i}: {path}')
if len(sys.path) > 5:
    print(f'  ... and {len(sys.path)-5} more entries')
print()

# Test RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    print('✅ RDKit imported successfully')
    
    # Test basic RDKit functionality
    mol = Chem.MolFromSmiles('CCO')
    if mol:
        print('✅ RDKit can create molecules from SMILES')
    else:
        print('❌ RDKit molecule creation failed')
        
except ImportError as e:
    print('❌ RDKit import failed:', e)

# Test Roshambo
try:
    import roshambo
    print('✅ Roshambo imported successfully')
    
    # Try to import the API
    from roshambo.api import get_similarity_scores
    print('✅ Roshambo API imported successfully')
    
except ImportError as e:
    print('❌ Roshambo import failed:', e)
    print('   Install with: pip install git+https://github.com/molecularinformatics/roshambo.git')

# Test PyTorch/CUDA
try:
    import torch
    print('✅ PyTorch imported successfully')
    print(f'   PyTorch version: {torch.__version__}')
    
    if torch.cuda.is_available():
        print(f'✅ CUDA available - {torch.cuda.device_count()} GPU(s)')
        for i in range(torch.cuda.device_count()):
            print(f'   GPU {i}: {torch.cuda.get_device_name(i)}')
    else:
        print('⚠  CUDA not available')
        
except ImportError as e:
    print('❌ PyTorch import failed:', e)

# Test NumPy
try:
    import numpy as np
    print('✅ NumPy imported successfully')
    print(f'   NumPy version: {np.__version__}')
except ImportError as e:
    print('❌ NumPy import failed:', e)
"

echo ""
echo "2. Testing environment variables:"
echo "RDBASE: $RDBASE"
echo "PYTHONPATH: $PYTHONPATH"
echo "CUDA_HOME: $CUDA_HOME"
echo "LD_LIBRARY_PATH (first 100 chars): ${LD_LIBRARY_PATH:0:100}..."

echo ""
echo "3. Testing file access:"
RDBASE="/home/protacinvent/Desktop/roshambo/rdkit"
if [ -f "$RDBASE/rdkit/__init__.py" ]; then
    echo "✅ RDKit Python module found"
elif [ -f "$RDBASE/Code/rdkit/__init__.py" ]; then
    echo "✅ RDKit Python module found in Code directory"
else
    echo "❌ RDKit Python module not found"
    echo "   Looking for: $RDBASE/rdkit/__init__.py"
fi

echo ""
echo "4. Testing a simple Roshambo operation:"
python -c "
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import tempfile
    import os
    
    # Create a simple test molecule
    mol = Chem.MolFromSmiles('CCO')
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    
    # Save to temporary SDF
    temp_file = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False)
    writer = Chem.SDWriter(temp_file.name)
    writer.write(mol)
    writer.close()
    
    print('✅ Created test SDF file:', temp_file.name)
    
    # Try to import Roshambo and test basic functionality
    try:
        from roshambo.api import get_similarity_scores
        print('✅ Roshambo API ready for use')
        print('   Note: Actual similarity calculation requires reference files')
    except Exception as e:
        print('❌ Roshambo API test failed:', e)
    
    # Cleanup
    os.unlink(temp_file.name)
    
except Exception as e:
    print('❌ Test molecule creation failed:', e)
"

echo ""
echo "================================"
echo "Test completed!"
echo ""
echo "If all tests passed (✅), your environment is ready for Roshambo."
echo "If any tests failed (❌), check the error messages above."
