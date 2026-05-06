#include "MacrocycleGenerator.h"
#include <iostream>
#include <iomanip>

int main() {
  std::cout << "Testing MacrocycleGenerator Phase 1 Implementation\n";
  std::cout << "==================================================\n\n";

  // Test 1: Simple 14-member ring
  {
    std::cout << "Test 1: 14-member macrocycle\n";
    RDDepict::MacrocycleGenerator gen(14);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";

    if (solved) {
      std::cout << "  Closure error: " << std::fixed << std::setprecision(6)
                << gen.getClosureError() << " Å\n";

      const auto& turns = gen.getTurns();
      int numR = 0, numL = 0;
      for (int t : turns) {
        if (t == 1) numR++;
        else if (t == -1) numL++;
      }
      std::cout << "  Turn counts: R=" << numR << ", L=" << numL
                << ", (R-L)=" << (numR - numL) << "\n";

      std::cout << "  Turn sequence: ";
      for (int t : turns) {
        std::cout << (t == 1 ? "R" : "L");
      }
      std::cout << "\n";

      auto coords = gen.generateCoordinates();
      std::cout << "  Generated " << coords.size() << " coordinates\n";

      std::cout << "  First few coordinates:\n";
      for (size_t i = 0; i < std::min(size_t(5), coords.size()); ++i) {
        std::cout << "    [" << i << "] (" << coords[i].x << ", "
                  << coords[i].y << ")\n";
      }
    }
    std::cout << "\n";
  }

  // Test 2: 12-member ring
  {
    std::cout << "Test 2: 12-member macrocycle\n";
    RDDepict::MacrocycleGenerator gen(12);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";

    if (solved) {
      std::cout << "  Closure error: " << std::fixed << std::setprecision(6)
                << gen.getClosureError() << " Å\n";
    }
    std::cout << "\n";
  }

  // Test 3: 10-member ring
  {
    std::cout << "Test 3: 10-member macrocycle\n";
    RDDepict::MacrocycleGenerator gen(10);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";

    if (solved) {
      std::cout << "  Closure error: " << std::fixed << std::setprecision(6)
                << gen.getClosureError() << " Å\n";
    }
    std::cout << "\n";
  }

  // Test 4: Ring too small (should fail)
  {
    std::cout << "Test 4: 5-member ring (too small, should fail)\n";
    RDDepict::MacrocycleGenerator gen(5);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";
    std::cout << "  (Expected: FAILED - ring too small for angular closure)\n";
    std::cout << "\n";
  }

  // Test 5: Ring with odd parity (should succeed with ~1.5 Å gap)
  {
    std::cout << "Test 5: 13-member ring (odd parity, imperfect closure)\n";
    RDDepict::MacrocycleGenerator gen(13);

    bool solved = gen.solve();
    std::cout << "  Solve result: " << (solved ? "SUCCESS" : "FAILED") << "\n";

    if (solved) {
      double error = gen.getClosureError();
      std::cout << "  Closure error: " << std::fixed << std::setprecision(6)
                << error << " Å\n";
      std::cout << "  (Expected: ~1.5 Å - one bond length gap for odd rings)\n";
    }
    std::cout << "\n";
  }

  std::cout << "All tests completed.\n";
  return 0;
}
