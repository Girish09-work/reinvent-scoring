"""
Roshambo Flask API for Reinvent Scoring.
This API provides a RESTful interface to the Roshambo shape similarity engine.
"""

import os
import shutil
import tempfile
from pathlib import Path
from flask import Flask, request, jsonify
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)

# Configuration
ROSHAMBO_WORKING_DIR = os.path.join(os.path.dirname(__file__), "inpdata")
Path(ROSHAMBO_WORKING_DIR).mkdir(parents=True, exist_ok=True)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({"status": "healthy", "service": "roshambo-api"})

@app.route('/similarity', methods=['POST'])
def calculate_similarity():
    """
    Calculate molecular shape similarity using Roshambo.

    Expected JSON payload:
    {
        "reference_file": "path/to/reference.sdf",
        "dataset_file": "path/to/dataset.sdf",
        "ignore_hs": true,
        "n_confs": 0,
        "use_carbon_radii": true,
        "color": true,
        "sort_by": "ComboTanimoto",
        "write_to_file": true,
        "gpu_id": 0,
        "working_dir": "path/to/working/dir"
    }
    """
    try:
        data = request.get_json()

        if not data:
            return jsonify({"success": False, "error": "No JSON data provided"}), 400

        # Extract parameters
        reference_file = data.get("reference_file")
        dataset_file = data.get("dataset_file")
        working_dir = data.get("working_dir", ROSHAMBO_WORKING_DIR)

        # Validate required parameters
        if not reference_file:
            return jsonify({"success": False, "error": "reference_file is required"}), 400
        if not dataset_file:
            return jsonify({"success": False, "error": "dataset_file is required"}), 400

        # Handle relative paths - make them relative to working directory
        if not os.path.isabs(reference_file):
            reference_file = os.path.join(working_dir, reference_file)
        if not os.path.isabs(dataset_file):
            dataset_file = os.path.join(working_dir, dataset_file)

        logger.info(f"Looking for reference file: {reference_file}")
        logger.info(f"Looking for dataset file: {dataset_file}")
        logger.info(f"Working directory: {working_dir}")

        # Validate files exist
        if not os.path.exists(reference_file):
            return jsonify({"success": False, "error": f"Reference file not found: {reference_file}"}), 400
        if not os.path.exists(dataset_file):
            return jsonify({"success": False, "error": f"Dataset file not found: {dataset_file}"}), 400

        # Create working directory
        Path(working_dir).mkdir(parents=True, exist_ok=True)

        # Copy files to working directory (only if they're not already there)
        ref_basename = os.path.basename(reference_file)
        dataset_basename = os.path.basename(dataset_file)
        local_ref_file = os.path.join(working_dir, ref_basename)
        local_dataset_file = os.path.join(working_dir, dataset_basename)

        # Only copy if source and destination are different
        if os.path.abspath(reference_file) != os.path.abspath(local_ref_file):
            shutil.copy2(reference_file, local_ref_file)
            logger.info(f"Copied reference file to: {local_ref_file}")
        else:
            logger.info(f"Reference file already in working directory: {local_ref_file}")

        if os.path.abspath(dataset_file) != os.path.abspath(local_dataset_file):
            shutil.copy2(dataset_file, local_dataset_file)
            logger.info(f"Copied dataset file to: {local_dataset_file}")
        else:
            logger.info(f"Dataset file already in working directory: {local_dataset_file}")

        logger.info(f"Processing similarity calculation in {working_dir}")
        logger.info(f"Reference: {ref_basename}, Dataset: {dataset_basename}")

        # Prepare roshambo parameters
        roshambo_params = {
            "ref_file": ref_basename,
            "dataset_files_pattern": dataset_basename,
            "ignore_hs": data.get("ignore_hs", True),
            "n_confs": data.get("n_confs", 0),
            "use_carbon_radii": data.get("use_carbon_radii", True),
            "color": data.get("color", True),
            "sort_by": data.get("sort_by", "ComboTanimoto"),
            "write_to_file": data.get("write_to_file", True),
            "gpu_id": data.get("gpu_id", 0),
            "working_dir": working_dir
        }

        # Call roshambo
        result = call_roshambo_api(**roshambo_params)

        if result["success"]:
            # Move output files to working directory if they're not already there
            move_roshambo_outputs(working_dir)

            return jsonify({
                "success": True,
                "message": "Similarity calculation completed",
                "working_dir": working_dir,
                "output_files": {
                    "csv": os.path.join(working_dir, "roshambo.csv"),
                    "mols_sdf": os.path.join(working_dir, "mols.sdf"),
                    "hits_sdf": os.path.join(working_dir, "hits.sdf")
                }
            })
        else:
            return jsonify(result), 500

    except Exception as e:
        logger.error(f"Error in similarity calculation: {e}")
        return jsonify({"success": False, "error": str(e)}), 500

def call_roshambo_api(**kwargs):
    """Call the roshambo API with the given parameters."""
    try:
        # Import roshambo here to handle import errors gracefully
        from roshambo.api import get_similarity_scores

        # Change to working directory
        original_cwd = os.getcwd()
        working_dir = kwargs.get("working_dir", ROSHAMBO_WORKING_DIR)

        try:
            os.chdir(working_dir)
            logger.info(f"Changed to working directory: {working_dir}")
            logger.info(f"Calling roshambo with parameters: {kwargs}")

            # Call roshambo (returns None but writes files)
            get_similarity_scores(**kwargs)

            # Check if output files were created
            expected_files = ["roshambo.csv", "mols.sdf", "hits.sdf"]
            created_files = []

            for filename in expected_files:
                if os.path.exists(filename):
                    created_files.append(filename)
                    logger.info(f"Created file: {filename}")

            if created_files:
                return {
                    "success": True,
                    "message": f"Roshambo completed successfully. Created files: {created_files}",
                    "created_files": created_files
                }
            else:
                return {
                    "success": False,
                    "error": "Roshambo completed but no output files were created"
                }

        finally:
            os.chdir(original_cwd)

    except ImportError as e:
        logger.error(f"Roshambo import error: {e}")
        return {
            "success": False,
            "error": f"Roshambo not available: {e}"
        }
    except Exception as e:
        logger.error(f"Roshambo execution error: {e}")
        return {
            "success": False,
            "error": f"Roshambo execution failed: {e}"
        }

def move_roshambo_outputs(target_dir):
    """Move roshambo output files to the target directory."""
    try:
        output_files = ["mols.sdf", "hits.sdf", "roshambo.csv"]

        for filename in output_files:
            # Check if file exists in current directory
            if os.path.exists(filename):
                target_path = os.path.join(target_dir, filename)
                if not os.path.exists(target_path):
                    shutil.move(filename, target_path)
                    logger.info(f"Moved {filename} to {target_path}")

    except Exception as e:
        logger.warning(f"Error moving output files: {e}")

if __name__ == '__main__':
    # Create inpdata directory
    Path(ROSHAMBO_WORKING_DIR).mkdir(parents=True, exist_ok=True)

    # Run the Flask app (localhost only)
    app.run(host='127.0.0.1', port=5000, debug=True)
