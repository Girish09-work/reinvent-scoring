#!/bin/bash

# =============================================================================
# Roshambo Environment Debugging Script
# =============================================================================
# This script helps diagnose environment setup issues
# Usage: bash debug_roshambo_env.sh
# =============================================================================

echo "=============================================================================
Roshambo Environment Debugging
============================================================================="

# Check conda installation
echo ""
echo "1. Checking Conda Installation:"
which conda
conda --version

# Check if roshambo environment exists
echo ""
echo "2. Checking Roshambo Environment:"
conda env list | grep roshambo || echo "❌ Roshambo environment not found"

# Activate environment and check packages
echo ""
echo "3. Checking Installed Packages:"
# Try different conda initialization paths
if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source /opt/conda/etc/profile.d/conda.sh
elif [ -f "/home/protacinvent/miniconda3/etc/profile.d/conda.sh" ]; then
    source /home/protacinvent/miniconda3/etc/profile.d/conda.sh
elif [ -f "/home/protacinvent/.conda/etc/profile.d/conda.sh" ]; then
    source /home/protacinvent/.conda/etc/profile.d/conda.sh
else
    echo "⚠ Could not find conda.sh, trying direct activation"
fi

conda activate roshambo 2>/dev/null && {
    echo "✅ Successfully activated roshambo environment"
    echo "Python version: $(python --version)"
    echo "Pip packages:"
    pip list | grep -E "(rdkit|roshambo|torch|numpy)" || echo "❌ Key packages not found"
} || {
    echo "❌ Failed to activate roshambo environment"
}

# Check directory structure
echo ""
echo "4. Checking Directory Structure:"
RDBASE="/home/protacinvent/Desktop/roshambo/rdkit"
echo "RDBASE: $RDBASE"
[ -d "$RDBASE" ] && echo "✅ RDKit directory exists" || echo "❌ RDKit directory missing"
[ -d "$RDBASE/lib" ] && echo "✅ RDKit lib directory exists" || echo "❌ RDKit lib directory missing"
[ -d "$RDBASE/Code" ] && echo "✅ RDKit Code directory exists" || echo "❌ RDKit Code directory missing"
[ -d "$RDBASE/Data" ] && echo "✅ RDKit Data directory exists" || echo "❌ RDKit Data directory missing"

# Check CUDA
echo ""
echo "5. Checking CUDA:"
which nvcc && nvcc --version || echo "❌ CUDA not found in PATH"
[ -d "/usr/local/cuda" ] && echo "✅ CUDA directory exists" || echo "❌ CUDA directory missing"

# Check current environment variables
echo ""
echo "6. Current Environment Variables:"
echo "CONDA_PREFIX: $CONDA_PREFIX"
echo "PYTHONPATH: $PYTHONPATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "CUDA_HOME: $CUDA_HOME"

# Test Python imports
echo ""
echo "7. Testing Python Imports:"
# Re-initialize conda if needed
if [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source /opt/conda/etc/profile.d/conda.sh
elif [ -f "/home/protacinvent/miniconda3/etc/profile.d/conda.sh" ]; then
    source /home/protacinvent/miniconda3/etc/profile.d/conda.sh
fi

conda activate roshambo 2>/dev/null && {
    python -c "
import sys
print('Python executable:', sys.executable)
print('Python path:')
for p in sys.path:
    print('  ', p)

print()
print('Testing imports:')
try:
    import rdkit
    print('✅ rdkit imported')
    from rdkit import Chem
    print('✅ rdkit.Chem imported')
    print('   RDKit version:', rdkit.__version__)
except Exception as e:
    print('❌ RDKit import failed:', e)

try:
    import roshambo
    print('✅ roshambo imported')
    from roshambo.api import get_similarity_scores
    print('✅ roshambo.api imported')
except Exception as e:
    print('❌ Roshambo import failed:', e)

try:
    import torch
    print('✅ torch imported')
    print('   CUDA available:', torch.cuda.is_available())
    if torch.cuda.is_available():
        print('   GPU count:', torch.cuda.device_count())
except Exception as e:
    print('❌ PyTorch import failed:', e)

try:
    import numpy
    print('✅ numpy imported')
except Exception as e:
    print('❌ NumPy import failed:', e)
"
} || {
    echo "❌ Cannot activate environment for testing"
}

echo ""
echo "=============================================================================
Debugging Complete
============================================================================="

# Provide recommendations
echo ""
echo "Recommendations:"
echo "1. If roshambo environment doesn't exist:"
echo "   conda create -n roshambo python=3.8"
echo "   conda activate roshambo"
echo "   pip install rdkit torch numpy"
echo "   pip install git+https://github.com/molecularinformatics/roshambo.git"
echo ""
echo "2. If RDKit directories are missing:"
echo "   Check if RDKit was built correctly in /home/protacinvent/Desktop/roshambo/rdkit"
echo ""
echo "3. If CUDA is not available:"
echo "   Check CUDA installation and PATH configuration"
echo ""
echo "4. To fix PYTHONPATH issues:"
echo "   unset PYTHONPATH before setting it in the script"
