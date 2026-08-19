#include "MacrocycleGenerator.h"
#include <iostream>
#include <iomanip>

using namespace RDDepict;

int main() {
  std::cout << "Testing MacrocycleGenerator Phase 2: Constraints\n";
  std::cout << "=================================================\n\n";

  // Test 1: FIXED constraint - 3-atom fusion pattern (RLR)
  {
    std::cout << "Test 1: 14-member ring with FIXED constraint (RLR at position 0)\n";
    MacrocycleGenerator gen(14);

    // Add constraint: positions 0, 1, 2 must be R, L, R (fused ring pattern)
    TurnConstraint constraint;
    constraint.position = 0;
    constraint.type = ConstraintType::FIXED;
    constraint.pattern = {1, -1, 1};  // R, L, R
    constraint.reason = "fused_ring_3_at_0";
    gen.addConstraint(constraint);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";

    if (solved) {
      std::cout << "  Closure error: " << std::fixed << std::setprecision(6)
                << gen.getClosureError() << " Å\n";

      const auto& turns = gen.getTurns();
      std::cout << "  Turn sequence: ";
      for (int t : turns) {
        std::cout << (t == 1 ? "R" : "L");
      }
      std::cout << "\n";

      // Verify constraint was applied
      if (turns[0] == 1 && turns[1] == -1 && turns[2] == 1) {
        std::cout << "  ✓ Constraint verified: positions 0,1,2 = RLR\n";
      } else {
        std::cout << "  ✗ Constraint VIOLATED!\n";
      }
    }
    std::cout << "\n";
  }

  // Test 2: SAME constraint - simulates cis double bond
  {
    std::cout << "Test 2: 12-member ring with SAME constraint at position 3\n";
    MacrocycleGenerator gen(12);

    // Add constraint: positions 3 and 4 must be same (cis double bond)
    TurnConstraint constraint;
    constraint.position = 3;
    constraint.type = ConstraintType::SAME;
    constraint.reason = "cis_double_bond_at_3";
    gen.addConstraint(constraint);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";

    if (solved) {
      std::cout << "  Closure error: " << std::fixed << std::setprecision(6)
                << gen.getClosureError() << " Å\n";

      const auto& turns = gen.getTurns();
      std::cout << "  Turn sequence: ";
      for (int t : turns) {
        std::cout << (t == 1 ? "R" : "L");
      }
      std::cout << "\n";

      // Verify constraint was applied
      if (turns[3] == turns[4]) {
        std::cout << "  ✓ Constraint verified: positions 3,4 are same ("
                  << (turns[3] == 1 ? "RR" : "LL") << ")\n";
      } else {
        std::cout << "  ✗ Constraint VIOLATED!\n";
      }
    }
    std::cout << "\n";
  }

  // Test 3: OPPOSITE constraint - simulates trans double bond
  {
    std::cout << "Test 3: 14-member ring with OPPOSITE constraint at position 5\n";
    MacrocycleGenerator gen(14);

    // Add constraint: positions 5 and 6 must be opposite (trans double bond)
    TurnConstraint constraint;
    constraint.position = 5;
    constraint.type = ConstraintType::OPPOSITE;
    constraint.reason = "trans_double_bond_at_5";
    gen.addConstraint(constraint);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";

    if (solved) {
      std::cout << "  Closure error: " << std::fixed << std::setprecision(6)
                << gen.getClosureError() << " Å\n";

      const auto& turns = gen.getTurns();
      std::cout << "  Turn sequence: ";
      for (int t : turns) {
        std::cout << (t == 1 ? "R" : "L");
      }
      std::cout << "\n";

      // Verify constraint was applied
      if (turns[5] != turns[6]) {
        std::cout << "  ✓ Constraint verified: positions 5,6 are opposite ("
                  << (turns[5] == 1 ? "RL" : "LR") << ")\n";
      } else {
        std::cout << "  ✗ Constraint VIOLATED!\n";
      }
    }
    std::cout << "\n";
  }

  // Test 4: Multiple constraints
  {
    std::cout << "Test 4: 16-member ring with multiple constraints\n";
    MacrocycleGenerator gen(16);

    // Add SAME constraint at position 2
    TurnConstraint c1;
    c1.position = 2;
    c1.type = ConstraintType::SAME;
    c1.reason = "cis_double_bond_at_2";
    gen.addConstraint(c1);

    // Add OPPOSITE constraint at position 7
    TurnConstraint c2;
    c2.position = 7;
    c2.type = ConstraintType::OPPOSITE;
    c2.reason = "trans_double_bond_at_7";
    gen.addConstraint(c2);

    // Add FIXED constraint at position 12
    TurnConstraint c3;
    c3.position = 12;
    c3.type = ConstraintType::FIXED;
    c3.pattern = {1, -1, 1};
    c3.reason = "fused_ring_3_at_12";
    gen.addConstraint(c3);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";

    if (solved) {
      std::cout << "  Closure error: " << std::fixed << std::setprecision(6)
                << gen.getClosureError() << " Å\n";

      const auto& turns = gen.getTurns();
      std::cout << "  Turn sequence: ";
      for (int t : turns) {
        std::cout << (t == 1 ? "R" : "L");
      }
      std::cout << "\n";

      // Verify all constraints
      bool allValid = true;
      if (turns[2] != turns[3]) {
        std::cout << "  ✗ Constraint 1 VIOLATED (positions 2,3)\n";
        allValid = false;
      }
      if (turns[7] == turns[8]) {
        std::cout << "  ✗ Constraint 2 VIOLATED (positions 7,8)\n";
        allValid = false;
      }
      if (turns[12] != 1 || turns[13] != -1 || turns[14] != 1) {
        std::cout << "  ✗ Constraint 3 VIOLATED (positions 12,13,14)\n";
        allValid = false;
      }
      if (allValid) {
        std::cout << "  ✓ All constraints verified!\n";
      }
    }
    std::cout << "\n";
  }

  // Test 5: Conflicting constraints (should fail)
  {
    std::cout << "Test 5: Conflicting constraints (should fail)\n";
    MacrocycleGenerator gen(12);

    // Add constraint: position 3 must be R
    TurnConstraint c1;
    c1.position = 3;
    c1.type = ConstraintType::FIXED;
    c1.pattern = {1};  // R
    c1.reason = "force_R_at_3";
    gen.addConstraint(c1);

    // Add constraint: position 3 must be L (conflicts!)
    TurnConstraint c2;
    c2.position = 3;
    c2.type = ConstraintType::FIXED;
    c2.pattern = {-1};  // L
    c2.reason = "force_L_at_3";
    gen.addConstraint(c2);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";
    std::cout << "  (Expected: FAILED - conflicting constraints)\n";
    std::cout << "\n";
  }

  // Test 6: Over-constrained (too many fixed turns)
  {
    std::cout << "Test 6: Over-constrained system (should fail)\n";
    MacrocycleGenerator gen(12);

    // Force all 12 positions to be R
    // But we need R-L = 6, so we need 9R and 3L, not 12R
    for (size_t i = 0; i < 12; ++i) {
      TurnConstraint c;
      c.position = i;
      c.type = ConstraintType::FIXED;
      c.pattern = {1};  // R
      c.reason = "force_all_R";
      gen.addConstraint(c);
    }

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";
    std::cout << "  (Expected: FAILED - angular closure cannot be satisfied)\n";
    std::cout << "\n";
  }

  std::cout << "All constraint tests completed.\n";
  return 0;
}
