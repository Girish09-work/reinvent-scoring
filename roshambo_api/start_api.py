#!/usr/bin/env python3
"""
Startup script for Roshambo Flask API.
This script sets up the environment and starts the Flask API server.
"""

import os
import sys
from pathlib import Path

def setup_environment():
    """Set up the environment for Roshambo."""
    print("🔧 Setting up Roshambo environment...")

    # Create necessary directories
    inpdata_dir = Path("inpdata")
    inpdata_dir.mkdir(exist_ok=True)
    print(f"✅ Created directory: {inpdata_dir}")

    # Check if roshambo is available
    try:
        import roshambo
        print("✅ Roshambo package is available")
        return True
    except ImportError:
        print("❌ Roshambo package not found")
        print("Please install roshambo in your conda environment:")
        print("  conda activate your_roshambo_env")
        print("  pip install git+https://github.com/rashatwi/roshambo.git")
        return False

def check_gpu():
    """Check GPU availability."""
    try:
        import torch
        if torch.cuda.is_available():
            gpu_count = torch.cuda.device_count()
            print(f"✅ GPU available: {gpu_count} device(s)")
            for i in range(gpu_count):
                gpu_name = torch.cuda.get_device_name(i)
                print(f"  GPU {i}: {gpu_name}")
            return True
        else:
            print("⚠️  No GPU available, will use CPU")
            return False
    except ImportError:
        print("⚠️  PyTorch not available, cannot check GPU")
        return False

def start_api(host='127.0.0.1', port=5000, debug=True):
    """Start the Flask API server."""
    print(f"🚀 Starting Roshambo Flask API on {host}:{port}")

    # Import and run the Flask app
    from app import app
    app.run(host=host, port=port, debug=debug)

def main():
    """Main function."""
    print("🧬 Roshambo Flask API Startup")
    print("=" * 40)

    # Setup environment
    if not setup_environment():
        sys.exit(1)

    # Check GPU
    check_gpu()

    # Start API
    try:
        start_api()
    except KeyboardInterrupt:
        print("\n👋 Shutting down Roshambo API")
    except Exception as e:
        print(f"❌ Error starting API: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
