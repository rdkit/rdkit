#include "MacrocycleGenerator.h"
#include <iostream>
#include <iomanip>

using namespace RDDepict;

int main() {
  std::cout << "Testing Self-Crossing Detection with Hexagonal Coordinates\n";
  std::cout << "===========================================================\n\n";

  // Test 1: Known non-crossing pattern for 14-member ring
  {
    std::cout << "Test 1: 14-member ring (should find non-crossing solution)\n";
    MacrocycleGenerator gen(14);

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

      // Calculate hex coordinates and verify no duplicates
      auto coords = gen.calculateHexCoords(turns);
      std::set<HexCoord> unique_coords;
      for (size_t i = 0; i < coords.size() - 1; ++i) {
        unique_coords.insert(coords[i]);
      }

      std::cout << "  Unique positions: " << unique_coords.size()
                << " (expected: " << (coords.size() - 1) << ")\n";

      if (unique_coords.size() == coords.size() - 1) {
        std::cout << "  ✓ No self-crossing detected\n";
      } else {
        std::cout << "  ✗ Self-crossing found!\n";
      }
    }
    std::cout << "\n";
  }

  // Test 2: 12-member ring
  {
    std::cout << "Test 2: 12-member ring (should find non-crossing solution)\n";
    MacrocycleGenerator gen(12);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";

    if (solved) {
      std::cout << "  Closure error: " << std::fixed << std::setprecision(6)
                << gen.getClosureError() << " Å\n";

      // Check for self-crossing using calculateMinDistance (returns 0 if crossing)
      const auto& turns = gen.getTurns();
      int minDist = gen.calculateMinDistance(turns);
      bool hasCrossing = (minDist == 0);
      std::cout << "  Self-crossing: " << (hasCrossing ? "YES (bad)" : "NO (good)") << "\n";
      std::cout << "  Min distance: " << minDist << " hex units\n";

      if (!hasCrossing) {
        std::cout << "  ✓ No self-crossing\n";
      } else {
        std::cout << "  ✗ Self-crossing detected!\n";
      }
    }
    std::cout << "\n";
  }

  // Test 3: Manually test a crossing pattern
  {
    std::cout << "Test 3: Manual crossing pattern test\n";

    // Create a pattern that should cross itself
    // For example: RRRRRR... (all right turns) would spiral inward
    std::vector<int> crossing_pattern(12, 1);  // All R turns

    MacrocycleGenerator gen(12);
    int minDist = gen.calculateMinDistance(crossing_pattern);
    bool crosses = (minDist == 0);

    std::cout << "  Pattern: ";
    for (int t : crossing_pattern) {
      std::cout << (t == 1 ? "R" : "L");
    }
    std::cout << "\n";
    std::cout << "  Self-crossing: " << (crosses ? "YES" : "NO") << "\n";
    std::cout << "  Min distance: " << minDist << " hex units\n";
    std::cout << "  " << (crosses ? "✓ Correctly detected crossing"
                                  : "✗ Should have detected crossing") << "\n";
    std::cout << "\n";
  }

  // Test 4: Print hex coordinates for visual inspection
  {
    std::cout << "Test 4: Hex coordinate visualization for simple pattern\n";

    std::vector<int> pattern = {1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1};  // RRLLRRLLRRLL

    MacrocycleGenerator gen(12);
    auto coords = gen.calculateHexCoords(pattern);

    std::cout << "  Pattern: ";
    for (int t : pattern) {
      std::cout << (t == 1 ? "R" : "L");
    }
    std::cout << "\n";

    std::cout << "  Hex coordinates:\n";
    for (size_t i = 0; i < coords.size(); ++i) {
      std::cout << "    [" << i << "] ("
                << std::setw(3) << coords[i].x << ", "
                << std::setw(3) << coords[i].y << ", "
                << std::setw(3) << coords[i].z << ")\n";
    }

    int minDist = gen.calculateMinDistance(pattern);
    bool crosses = (minDist == 0);
    std::cout << "  Self-crossing: " << (crosses ? "YES" : "NO") << "\n";
    std::cout << "  Min distance: " << minDist << " hex units\n";
    std::cout << "\n";
  }

  // Test 5: Test with constraints
  {
    std::cout << "Test 5: 14-member ring with SAME constraint\n";
    MacrocycleGenerator gen(14);

    TurnConstraint c1;
    c1.position = 3;
    c1.type = ConstraintType::SAME;
    c1.reason = "test_constraint";
    gen.addConstraint(c1);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";

    if (solved) {
      const auto& turns = gen.getTurns();
      int minDist = gen.calculateMinDistance(turns);
      bool crosses = (minDist == 0);
      std::cout << "  Self-crossing: " << (crosses ? "YES (bad)" : "NO (good)") << "\n";
      std::cout << "  Min distance: " << minDist << " hex units\n";
      std::cout << "  Closure error: " << std::fixed << std::setprecision(6)
                << gen.getClosureError() << " Å\n";
    }
    std::cout << "\n";
  }

  std::cout << "All self-crossing tests completed.\n";
  return 0;
}
