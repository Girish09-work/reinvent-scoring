#!/usr/bin/env python3
"""
Test script to verify that the scaffold CSV generation fix works correctly.

This script tests that the base score component properly passes the step parameter
to scoring components during reinforcement learning, which is required for the
diversity filter system to generate scaffold_memory.csv files.
"""

import os
import sys
import tempfile
import numpy as np
from typing import List

# Add the current directory to Python path to import reinvent_scoring
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from reinvent_scoring.scoring.component_parameters import ComponentParameters
from reinvent_scoring.scoring.score_summary import ComponentSummary
from reinvent_scoring.scoring.score_components.base_score_component import BaseScoreComponent


class TestScoreComponent(BaseScoreComponent):
    """Test scoring component to verify step parameter passing."""

    def __init__(self, parameters: ComponentParameters):
        super().__init__(parameters)
        self.last_step = None
        self.call_count = 0

    def calculate_score(self, molecules: List, step=-1) -> ComponentSummary:
        """Test implementation that tracks the step parameter."""
        self.last_step = step
        self.call_count += 1

        # Return dummy scores
        scores = np.array([0.5] * len(molecules), dtype=np.float32)
        return ComponentSummary(total_score=scores, parameters=self.parameters)


def test_step_parameter_passing():
    """Test that calculate_score_for_step properly passes the step parameter."""
    print("🧪 Testing step parameter passing...")

    # Create test component
    params = ComponentParameters(
        component_type="test_component",
        name="Test Component",
        weight=1.0
    )
    component = TestScoreComponent(params)

    # Test molecules (dummy SMILES)
    test_molecules = ["CCO", "CCC", "CCCC"]
    test_step = 42

    # Call calculate_score_for_step (this is what the system calls during RL)
    result = component.calculate_score_for_step(test_molecules, step=test_step)

    # Verify the step was passed correctly
    assert component.last_step == test_step, f"Expected step {test_step}, got {component.last_step}"
    assert component.call_count == 1, f"Expected 1 call, got {component.call_count}"
    assert len(result.total_score) == len(test_molecules), "Score array length mismatch"

    print(f"✅ Step parameter correctly passed: {component.last_step}")
    print(f"✅ Component called {component.call_count} time(s)")
    print(f"✅ Returned {len(result.total_score)} scores for {len(test_molecules)} molecules")

    return True


def test_old_vs_new_behavior():
    """Test the difference between old (broken) and new (fixed) behavior."""
    print("\n🧪 Testing old vs new behavior...")

    # Create test component
    params = ComponentParameters(
        component_type="test_component",
        name="Test Component",
        weight=1.0
    )
    component = TestScoreComponent(params)

    # Test molecules
    test_molecules = ["CCO", "CCC"]
    test_step = 123

    # Test the NEW behavior (after fix)
    print("  Testing NEW behavior (after fix):")
    component.last_step = None
    result = component.calculate_score_for_step(test_molecules, step=test_step)
    print(f"    calculate_score_for_step({test_step}) -> calculate_score received step: {component.last_step}")

    # Verify the fix works
    assert component.last_step == test_step, f"NEW behavior failed: Expected {test_step}, got {component.last_step}"
    print(f"    ✅ Step parameter correctly passed: {test_step}")

    # Simulate OLD behavior (what would happen before the fix)
    print("  Simulating OLD behavior (before fix):")
    component.last_step = None
    # This simulates the old broken base class method: return self.calculate_score(molecules)
    old_result = component.calculate_score(test_molecules)  # No step parameter!
    print(f"    calculate_score() without step -> received step: {component.last_step}")

    # Show the difference
    assert component.last_step == -1, f"OLD behavior simulation failed: Expected -1, got {component.last_step}"
    print(f"    ❌ Step parameter was lost (default -1 used)")

    print("  📊 Comparison:")
    print(f"    OLD: calculate_score_for_step(step={test_step}) -> calculate_score() -> step=-1")
    print(f"    NEW: calculate_score_for_step(step={test_step}) -> calculate_score(step={test_step}) -> step={test_step}")

    return True


def test_base_class_behavior():
    """Test the base class calculate_score_for_step method directly."""
    print("\n🧪 Testing base class behavior...")

    # Create test component
    params = ComponentParameters(
        component_type="test_component",
        name="Test Component",
        weight=1.0
    )
    component = TestScoreComponent(params)

    # Test with different step values
    test_cases = [
        (["CCO"], -1),  # Default step
        (["CCC", "CCCC"], 0),  # Step 0
        (["CCCCC"], 100),  # Step 100
    ]

    for molecules, step in test_cases:
        component.last_step = None  # Reset
        result = component.calculate_score_for_step(molecules, step=step)

        assert component.last_step == step, f"Step {step}: Expected {step}, got {component.last_step}"
        print(f"✅ Step {step}: Correctly passed to calculate_score")

    return True


def main():
    """Run all tests."""
    print("🔧 Testing scaffold CSV generation fix...")
    print("=" * 60)

    try:
        # Test 1: Basic step parameter passing
        test_step_parameter_passing()

        # Test 2: Old vs new behavior comparison
        test_old_vs_new_behavior()

        # Test 3: Base class behavior
        test_base_class_behavior()

        print("\n" + "=" * 60)
        print("🎉 ALL TESTS PASSED!")
        print("\n📋 Summary:")
        print("  ✅ Base score component properly passes step parameter")
        print("  ✅ calculate_score_for_step works correctly")
        print("  ✅ Scoring components will now integrate with diversity filter")
        print("  ✅ scaffold_memory.csv should be generated during reinforcement learning")

        print("\n💡 Next steps:")
        print("  1. Run reinforcement learning with rdkit_shape or roshambo components")
        print("  2. Check that scaffold_memory.csv is generated in results directory")
        print("  3. Verify CSV contains component scores and scaffold information")

        return True

    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
