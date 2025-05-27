#!/usr/bin/env python3
"""
Test script to verify CSV parsing functionality with actual Roshambo CSV data.
"""

import pandas as pd
import tempfile
import os

def test_csv_parsing():
    """Test CSV parsing with the actual data format you provided."""
    
    # Create a test CSV file with your actual data format
    csv_content = """Molecule	OriginalName	ComboTanimoto	ShapeTanimoto	ColorTanimoto	FitTverskyCombo	FitTversky	FitColorTversky	RefTverskyCombo	RefTversky	RefColorTversky	Overlap
mol_126_0	mol_126	0.35	0.291	0.059	0.82	0.678	0.142	0.43	0.338	0.092	1269.983
mol_2200_0	mol_2200	0.332	0.242	0.091	0.718	0.535	0.184	0.458	0.306	0.151	1156.798
mol_1524_0	mol_1524	0.323	0.266	0.057	0.813	0.664	0.149	0.391	0.307	0.084	1154.463
mol_939_0	mol_939	0.321	0.252	0.069	0.822	0.645	0.178	0.394	0.293	0.101	1100.4
mol_67_0	mol_67	0.312	0.241	0.071	0.668	0.512	0.156	0.429	0.313	0.116	1184.284
mol_1483_0	mol_1483	0.311	0.223	0.088	0.922	0.662	0.259	0.369	0.252	0.117	941.963
mol_2294_0	mol_2294	0.307	0.25	0.057	0.707	0.568	0.139	0.398	0.309	0.089	1164.302
mol_1379_0	mol_1379	0.296	0.221	0.075	0.679	0.504	0.175	0.399	0.282	0.117	1065.392
mol_302_0	mol_302	0.296	0.219	0.077	0.816	0.605	0.211	0.363	0.256	0.108	957.823
mol_1063_0	mol_1063	0.294	0.244	0.05	0.705	0.593	0.112	0.376	0.292	0.084	1099.79
mol_1125_0	mol_1125	0.294	0.221	0.073	0.703	0.547	0.156	0.39	0.27	0.12	1017.247"""

    # Write to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_content)
        csv_file = f.name

    try:
        print("🧪 Testing CSV parsing with actual Roshambo data format")
        print("=" * 60)
        
        # Test reading with different separators
        print("\n📄 Testing CSV reading with different separators:")
        
        for sep_name, sep in [("tab", '\t'), ("comma", ',')]:
            try:
                df = pd.read_csv(csv_file, sep=sep)
                print(f"   {sep_name}: ✅ Success - {df.shape[0]} rows, {df.shape[1]} columns")
                if len(df.columns) > 1:
                    print(f"      Columns: {df.columns.tolist()}")
                    break
            except Exception as e:
                print(f"   {sep_name}: ❌ Failed - {e}")
        
        if 'df' not in locals() or df.empty:
            print("❌ Could not read CSV file with any separator")
            return False
        
        print(f"\n📊 CSV Data Analysis:")
        print(f"   Shape: {df.shape}")
        print(f"   Columns: {df.columns.tolist()}")
        print(f"   Sample data:")
        print(df.head(3).to_string(index=False))
        
        # Test the parsing logic from the component
        print(f"\n🔍 Testing molecule name parsing:")
        
        shape_weight = 0.7
        color_weight = 0.3
        
        mol_scores = {}
        
        for _, row in df.iterrows():
            # Try different column names for molecule identifier
            name = row.get("Molecule", "") or row.get("Name", "") or row.get("OriginalName", "")
            
            print(f"   Processing row with name: '{name}'")
            
            if name and ("mol_" in str(name)):
                try:
                    # Handle different name formats: mol_126_0 -> 126, mol_126 -> 126
                    name_str = str(name)
                    if name_str.startswith("mol_"):
                        parts = name_str.split("_")
                        if len(parts) >= 2:
                            idx = int(parts[1])  # Extract the molecule index
                            print(f"      Extracted index: {idx}")
                        else:
                            print(f"      ❌ Could not extract index from {name_str}")
                            continue
                    else:
                        print(f"      ❌ Name doesn't start with 'mol_': {name_str}")
                        continue
                        
                    shape_score = float(row.get("ShapeTanimoto", 0.0))
                    color_score = float(row.get("ColorTanimoto", 0.0))
                    
                    # Calculate combined score
                    if shape_weight + color_weight > 0:
                        combo_score = ((shape_weight * shape_score) + 
                                     (color_weight * color_score)) / (shape_weight + color_weight)
                    else:
                        combo_score = 0.0
                        
                    mol_scores[idx] = max(mol_scores.get(idx, 0.0), combo_score)
                    
                    print(f"      ✅ Molecule {idx}: shape={shape_score:.3f}, color={color_score:.3f}, combo={combo_score:.3f}")
                        
                except (ValueError, IndexError) as e:
                    print(f"      ❌ Error processing row with name '{name}': {e}")
        
        print(f"\n📊 Final Results:")
        print(f"   Processed {len(mol_scores)} molecules")
        print(f"   Score mapping: {dict(sorted(mol_scores.items()))}")
        
        # Test creating a scores array
        if mol_scores:
            max_idx = max(mol_scores.keys())
            scores_array = [0.0] * (max_idx + 1)
            
            for idx, score in mol_scores.items():
                scores_array[idx] = score
            
            print(f"   Scores array (first 10): {scores_array[:10]}")
            print(f"   Non-zero scores: {[(i, s) for i, s in enumerate(scores_array) if s > 0]}")
        
        return True
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        # Clean up
        try:
            os.unlink(csv_file)
        except:
            pass

if __name__ == "__main__":
    print("🔬 CSV Parsing Test")
    print("Testing with actual Roshambo CSV data format")
    
    success = test_csv_parsing()
    
    if success:
        print("\n✅ CSV parsing test completed successfully!")
    else:
        print("\n❌ CSV parsing test failed!")
