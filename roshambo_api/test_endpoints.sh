#!/bin/bash

# Test script for Roshambo Flask API endpoints

echo "🧪 Testing Roshambo Flask API Endpoints"
echo "========================================"

# Get the current directory
CURRENT_DIR=$(pwd)
API_DIR="$CURRENT_DIR/roshambo_api"

echo "API Directory: $API_DIR"
echo ""

# Test 1: Health Check
echo "1. Testing Health Check Endpoint"
echo "--------------------------------"
echo "Command: curl -X GET http://localhost:5000/health"
echo ""

curl -X GET http://localhost:5000/health
echo ""
echo ""

# Test 2: Similarity Calculation with example files
echo "2. Testing Similarity Calculation Endpoint"
echo "------------------------------------------"

# Create the JSON payload
cat > /tmp/roshambo_test.json << EOF
{
  "reference_file": "$API_DIR/inpdata/query.sdf",
  "dataset_file": "$API_DIR/inpdata/dataset.sdf",
  "ignore_hs": true,
  "n_confs": 0,
  "use_carbon_radii": true,
  "color": true,
  "sort_by": "ComboTanimoto",
  "write_to_file": true,
  "gpu_id": 0,
  "working_dir": "$API_DIR/inpdata"
}
EOF

echo "JSON Payload:"
cat /tmp/roshambo_test.json
echo ""
echo ""

echo "Command: curl -X POST http://localhost:5000/similarity -H 'Content-Type: application/json' -d @/tmp/roshambo_test.json"
echo ""

curl -X POST http://localhost:5000/similarity \
  -H "Content-Type: application/json" \
  -d @/tmp/roshambo_test.json

echo ""
echo ""

# Test 3: Check if output files were created
echo "3. Checking Output Files"
echo "------------------------"

if [ -f "$API_DIR/inpdata/roshambo.csv" ]; then
    echo "✅ roshambo.csv created"
    echo "First few lines:"
    head -5 "$API_DIR/inpdata/roshambo.csv"
else
    echo "❌ roshambo.csv not found"
fi

echo ""

if [ -f "$API_DIR/inpdata/mols.sdf" ]; then
    echo "✅ mols.sdf created"
    echo "File size: $(wc -l < "$API_DIR/inpdata/mols.sdf") lines"
else
    echo "❌ mols.sdf not found"
fi

echo ""

if [ -f "$API_DIR/inpdata/hits.sdf" ]; then
    echo "✅ hits.sdf created"
    echo "File size: $(wc -l < "$API_DIR/inpdata/hits.sdf") lines"
else
    echo "❌ hits.sdf not found"
fi

echo ""
echo "🏁 Test completed!"

# Cleanup
rm -f /tmp/roshambo_test.json
