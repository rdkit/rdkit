//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MacrocycleGenerator.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <set>

namespace RDDepict {

MacrocycleGenerator::MacrocycleGenerator(size_t ringSize, double bondLength,
                                         bool useJacobianRefinement)
    : d_ringSize(ringSize),
      d_bondLength(bondLength),
      d_closureError(std::numeric_limits<double>::max()),
      d_solved(false),
      d_useJacobianRefinement(useJacobianRefinement) {
  d_turns.resize(ringSize, 0);
}

void MacrocycleGenerator::addConstraint(const TurnConstraint &constraint) {
  d_constraints.push_back(constraint);
}

bool MacrocycleGenerator::applyConstraints(std::vector<int> &turns) const {
  for (const auto &constraint : d_constraints) {
    size_t pos = constraint.position;

    if (constraint.type == ConstraintType::FIXED) {
      std::cerr << "Applying FIXED constraint at position " << pos
                << ": pattern = ";
      for (int turn : constraint.pattern) {
        std::cerr << (turn == 1 ? "R " : turn == -1 ? "L " : "0 ");
      }
      std::cerr << std::endl;
      std::cerr << "Before applying FIXED constraint: ";
      for (size_t i = 0; i < turns.size(); ++i) {
        std::cerr << "Atom " << i << ": "
                  << (turns[i] == 1    ? "R "
                      : turns[i] == -1 ? "L "
                                       : "0 ")
                  << std::endl;
      }
      // Apply fixed pattern
      for (size_t i = 0; i < constraint.pattern.size(); ++i) {
        auto idx = (pos + i) % d_ringSize;
        if (turns[idx] == 0) {
          turns[idx] = constraint.pattern[i];
        } else if (turns[idx] != constraint.pattern[i]) {
          // Conflict: position already set to different value
          std::cerr << "Constraints conflict - no solution possible"
                    << std::endl;
          return false;
        }
      }
    } else if (constraint.type == ConstraintType::SAME) {
      // Two consecutive turns must be same
      size_t pos1 = pos % d_ringSize;
      size_t pos2 = (pos + 1) % d_ringSize;

      if (turns[pos1] != 0 && turns[pos2] == 0) {
        turns[pos2] = turns[pos1];
      } else if (turns[pos1] == 0 && turns[pos2] != 0) {
        turns[pos1] = turns[pos2];
      } else if (turns[pos1] != 0 && turns[pos2] != 0 &&
                 turns[pos1] != turns[pos2]) {
        // Conflict: already set to different values
        std::cerr << "Constraints 2222 conflict - no solution possible"
                  << std::endl;
        return false;
      }
    } else if (constraint.type == ConstraintType::OPPOSITE) {
      // Two consecutive turns must be opposite
      size_t pos1 = pos % d_ringSize;
      size_t pos2 = (pos + 1) % d_ringSize;

      if (turns[pos1] != 0 && turns[pos2] == 0) {
        turns[pos2] = -turns[pos1];
      } else if (turns[pos1] == 0 && turns[pos2] != 0) {
        turns[pos1] = -turns[pos2];
      } else if (turns[pos1] != 0 && turns[pos2] != 0 &&
                 turns[pos1] == turns[pos2]) {
        // Conflict: already set to same value
        return false;
      }
    }
  }

  return true;
}

bool MacrocycleGenerator::validateConstraints(
    const std::vector<int> &turns) const {
  for (const auto &constraint : d_constraints) {
    size_t pos = constraint.position;
    std::cerr << "validating constraints, sequence is:" << std::endl;
    for (size_t i = 0; i < turns.size(); ++i) {
      std::cerr << "Atom " << i << ": "
                << (turns[i] == 1    ? "R "
                    : turns[i] == -1 ? "L "
                                     : "0 ")
                << std::endl;
    }

    if (constraint.type == ConstraintType::FIXED) {
      // Check fixed pattern matches
      for (size_t i = 0; i < constraint.pattern.size(); ++i) {
        auto idx = (pos + i) % d_ringSize;
        if (turns[idx] != constraint.pattern[i]) {
          return false;
        }
      }
    } else if (constraint.type == ConstraintType::SAME) {
      // Check two consecutive turns are same
      size_t pos1 = pos % d_ringSize;
      size_t pos2 = (pos + 1) % d_ringSize;
      if (turns[pos1] != turns[pos2]) {
        return false;
      }
    } else if (constraint.type == ConstraintType::OPPOSITE) {
      // Check two consecutive turns are opposite
      size_t pos1 = pos % d_ringSize;
      size_t pos2 = (pos + 1) % d_ringSize;
      if (turns[pos1] == turns[pos2]) {
        return false;
      }
    }
  }

  return true;
}

bool MacrocycleGenerator::solve() {
  // Initialize all turns as undecided
  d_turns.assign(d_ringSize, 0);

  // Apply structural constraints
  if (!applyConstraints(d_turns)) {
    std::cerr << "Constraints conflict1 - no solution possible" << std::endl;
    return false;  // Constraints conflict
  }

  // Count already-decided turns
  int decidedRight = 0;
  int decidedLeft = 0;
  int freeCount = 0;

  for (int turn : d_turns) {
    if (turn == 1) {
      ++decidedRight;
    } else if (turn == -1) {
      ++decidedLeft;
    } else {
      ++freeCount;
    }
  }

  // Angular closure constraint: count(R) - count(L) = 6
  // For even-numbered rings: exact angular closure with R - L = 6
  // For odd-numbered rings: use R - L = 5, then distribute extra 60° across all
  // angles
  bool isOddRing = (d_ringSize % 2) == 1;
  int targetDiff = isOddRing ? 5 : 6;

  // We need: (decidedRight + freeRight) - (decidedLeft + freeLeft) = targetDiff
  // Where: freeRight + freeLeft = freeCount
  // Solving: freeRight = (targetDiff - decidedRight + decidedLeft + freeCount)
  // / 2

  int totalRight = (d_ringSize + targetDiff) / 2;
  int totalLeft = (d_ringSize - targetDiff) / 2;

  // Check parity
  if ((d_ringSize + targetDiff) % 2 != 0) {
    return false;  // Parity mismatch - no integer solution
  }

  // Check if solution is possible
  if (totalRight < 0 || totalLeft < 0) {
    return false;  // Ring too small for angular closure
  }

  int freeRight = totalRight - decidedRight;
  int freeLeft = totalLeft - decidedLeft;

  if (freeRight < 0 || freeLeft < 0) {
    return false;  // Constraints over-constrain the solution
  }

  if (freeRight + freeLeft != freeCount) {
    return false;  // Sanity check failed
  }

  // Find optimal turn sequence for free variables
  d_solved = findOptimalTurnSequence(freeRight, freeLeft);
  return d_solved;
}

bool MacrocycleGenerator::findOptimalTurnSequence(int numRight, int numLeft) {
  std::cerr << "Finding optimal turn sequence: numRight = " << numRight
            << ", numLeft = " << numLeft << std::endl;
  // Find positions of free variables (currently 0)
  std::vector<size_t> freePositions;
  for (size_t i = 0; i < d_ringSize; ++i) {
    if (d_turns[i] == 0) {
      freePositions.push_back(i);
    }
  }

  // If no free variables, we're done (constraints fully determined the
  // solution)
  if (freePositions.empty()) {
    d_closureError = calculateClosureError(d_turns);
    return true;
  }

  // If all positions are free (no constraints), use old logic
  if (freePositions.size() == d_ringSize) {
    std::cerr << "No constraints - enumerating all combinations" << std::endl;
    std::vector<int> candidate(d_ringSize, -1);  // Initialize all to L
    std::vector<size_t> rightPositions(numRight);

    // Initialize first combination: [0, 1, 2, ..., numRight-1]
    for (int i = 0; i < numRight; ++i) {
      rightPositions[i] = i;
      candidate[i] = 1;  // Set to R
    }

    double bestError = std::numeric_limits<double>::max();
    std::vector<int> bestTurns;

    // Evaluate first combination
    if (validateConstraints(candidate) && !hasSelfCrossing(candidate)) {
      double error = calculateClosureError(candidate);
      if (error < bestError) {
        bestError = error;
        bestTurns = candidate;
      }
    }
    std::cerr << "Testing combinations: " << std::endl;

    // Enumerate remaining combinations
    while (nextCombination(rightPositions, d_ringSize)) {
      std::cerr << "Testing combination: ";
      for (size_t pos : rightPositions) {
        std::cerr << pos << " ";
      }
      // Reset candidate to all L
      std::fill(candidate.begin(), candidate.end(), -1);

      // Set R positions
      for (size_t pos : rightPositions) {
        candidate[pos] = 1;
      }

      // Validate constraints and check for self-crossing before evaluating
      if (validateConstraints(candidate)) {
        if (hasSelfCrossing(candidate)) {
          std::cerr << "Self-crossing detected - rejecting" << std::endl;
        } else {
          std::cerr << "Valid combination. Calculating error..." << std::endl;
          // Calculate error
          double error = calculateClosureError(candidate);
          std::cerr << "Closure error: " << error << std::endl;
          if (error < bestError) {
            bestError = error;
            bestTurns = candidate;
          }
        }
      }
    }

    // Store best solution
    if (bestTurns.empty()) {
      return false;  // No valid solution found that satisfies constraints
    }

    d_turns = bestTurns;
    d_closureError = bestError;

    return true;
  }

  std::cerr
      << "Constraints present - enumerating combinations for free positions only"
      << std::endl;

  // Mixed case: some positions fixed, some free
  // Enumerate ways to assign R/L to free positions
  std::vector<size_t> rightPositions(numRight);
  size_t numFree = freePositions.size();

  // Initialize first combination: choose first numRight free positions for R
  for (int i = 0; i < numRight; ++i) {
    rightPositions[i] = i;
  }

  double bestError = std::numeric_limits<double>::max();
  std::vector<int> bestTurns;

  // Helper: create candidate from current combination
  auto createCandidate = [&]() {
    std::vector<int> candidate = d_turns;  // Start with constrained values
    std::cerr << "  Creating candidate from constraints: ";
    for (size_t i = 0; i < d_turns.size(); ++i) {
      std::cerr << (d_turns[i] == 1 ? "R" : d_turns[i] == -1 ? "L" : "0");
    }
    std::cerr << std::endl;

    // Set free positions to L
    for (size_t freeIdx : freePositions) {
      candidate[freeIdx] = -1;
    }
    // Set selected positions to R
    std::cerr << "  Setting R at free position indices: ";
    for (size_t idx : rightPositions) {
      std::cerr << idx << " (pos " << freePositions[idx] << ") ";
      candidate[freePositions[idx]] = 1;
    }
    std::cerr << std::endl;
    std::cerr << "  Resulting candidate: ";
    for (size_t i = 0; i < candidate.size(); ++i) {
      std::cerr << (candidate[i] == 1 ? "R" : candidate[i] == -1 ? "L" : "0");
    }
    std::cerr << std::endl;
    return candidate;
  };

  // Evaluate first combination
  std::cerr << "Evaluating first combination:" << std::endl;
  std::vector<int> candidate = createCandidate();
  if (validateConstraints(candidate)) {
    std::cerr << "  Validation passed. Checking for self-crossing..."
              << std::endl;
    if (hasSelfCrossing(candidate)) {
      std::cerr << "  Self-crossing detected - rejecting" << std::endl;
    } else {
      std::cerr << "  No self-crossing. Calculating error..." << std::endl;
      double error = calculateClosureError(candidate);
      std::cerr << "  Closure error: " << error << std::endl;
      if (error < bestError) {
        bestError = error;
        bestTurns = candidate;
      }
    }
  } else {
    std::cerr << "  Validation FAILED for first combination" << std::endl;
  }

  // Enumerate remaining combinations
  int combinationCount = 1;
  while (nextCombination(rightPositions, numFree)) {
    combinationCount++;
    std::cerr << "\nEvaluating combination " << combinationCount << ":"
              << std::endl;
    candidate = createCandidate();
    if (validateConstraints(candidate)) {
      std::cerr << "  Validation passed. Checking for self-crossing..."
                << std::endl;
      if (hasSelfCrossing(candidate)) {
        std::cerr << "  Self-crossing detected - rejecting" << std::endl;
      } else {
        std::cerr << "  No self-crossing. Calculating error..." << std::endl;
        double error = calculateClosureError(candidate);
        std::cerr << "  Closure error: " << error << std::endl;
        if (error < bestError) {
          bestError = error;
          bestTurns = candidate;
        }
      }
    } else {
      std::cerr << "  Validation FAILED" << std::endl;
    }
  }

  std::cerr << "\nTested " << combinationCount << " combinations total"
            << std::endl;

  // Store best solution
  if (bestTurns.empty()) {
    std::cerr << "ERROR: No valid solution found that satisfies constraints!"
              << std::endl;
    return false;  // No valid solution found that satisfies constraints
  }

  std::cerr << "Found best solution with error " << bestError << ": ";
  for (size_t i = 0; i < bestTurns.size(); ++i) {
    std::cerr << (bestTurns[i] == 1 ? "R" : bestTurns[i] == -1 ? "L" : "0");
  }
  std::cerr << std::endl;

  d_turns = bestTurns;
  d_closureError = bestError;

  return true;
}

bool MacrocycleGenerator::nextCombination(std::vector<size_t> &positions,
                                          size_t n) {
  // Generate next lexicographic combination
  // positions is sorted in increasing order

  int k = positions.size();
  if (k == 0) return false;

  // Find rightmost position that can be incremented
  int i = k - 1;
  while (i >= 0 && positions[i] == n - k + i) {
    --i;
  }

  if (i < 0) {
    return false;  // No more combinations
  }

  // Increment this position
  positions[i]++;

  // Set all positions to the right to consecutive values
  for (int j = i + 1; j < k; ++j) {
    positions[j] = positions[j - 1] + 1;
  }

  return true;
}

double MacrocycleGenerator::calculateClosureError(
    const std::vector<int> &turns) const {
  // Calculate final position after following the turn sequence

  double x = 0.0;
  double y = 0.0;
  double angle = 0.0;  // Start facing east (0 radians)

  for (size_t i = 0; i < turns.size(); ++i) {
    // Turn at atom i: ±60° = ±π/3
    angle += turns[i] * (M_PI / 3.0);

    // Advance in new direction
    x += d_bondLength * std::cos(angle);
    y += d_bondLength * std::sin(angle);
  }

  // Calculate distance from final position to origin
  return std::sqrt(x * x + y * y);
}

std::vector<HexCoord> MacrocycleGenerator::calculateHexCoords(
    const std::vector<int> &turns) const {
  // Direction vectors in hexagonal cube coordinates
  // 6 directions at 0°, 60°, 120°, 180°, 240°, 300°
  static const HexCoord directions[6] = {
      HexCoord(+1, -1, 0),  // 0° (East)
      HexCoord(+1, 0, -1),  // 60° (Northeast)
      HexCoord(0, +1, -1),  // 120° (Northwest)
      HexCoord(-1, +1, 0),  // 180° (West)
      HexCoord(-1, 0, +1),  // 240° (Southwest)
      HexCoord(0, -1, +1)   // 300° (Southeast)
  };

  std::vector<HexCoord> coords;
  coords.reserve(turns.size() + 1);

  // Start at origin
  HexCoord pos(0, 0, 0);
  coords.push_back(pos);

  // Start facing direction 0 (East)
  int direction = 0;

  for (size_t i = 0; i < turns.size(); ++i) {
    // Apply turn at atom i (affects direction from i to i+1)
    if (turns[i] == 1) {
      direction = (direction + 1) % 6;  // Right turn
    } else if (turns[i] == -1) {
      direction = (direction + 5) % 6;  // Left turn
    }

    // Move one step in current direction
    const HexCoord &dir = directions[direction];
    pos.x += dir.x;
    pos.y += dir.y;
    pos.z += dir.z;
    coords.push_back(pos);
  }

  return coords;
}

bool MacrocycleGenerator::hasSelfCrossing(const std::vector<int> &turns) const {
  // Calculate hex coordinates for all nodes
  auto coords = calculateHexCoords(turns);

  // Check for duplicate positions (excluding first and last which should match)
  std::set<HexCoord> visited;

  for (size_t i = 0; i < coords.size() - 1; ++i) {
    if (!visited.insert(coords[i]).second) {
      // Duplicate found - self-crossing detected
      return true;
    }
  }

  return false;
}

std::vector<RDGeom::Point2D> MacrocycleGenerator::generateCoordinates() const {
  if (!d_solved) {
    throw std::runtime_error(
        "MacrocycleGenerator: must call solve() before generateCoordinates()");
  }

  std::vector<RDGeom::Point2D> coords;
  coords.reserve(d_ringSize + 1);  // +1 for dummy atom

  bool isOddRing = (d_ringSize % 2) == 1;

  if (isOddRing) {
    // Odd rings: Generate coordinates with extra angle per atom, then use dummy
    double extraPerAtom = (M_PI / 3.0) / d_ringSize;

    RDGeom::Point2D pos(0, 0);
    coords.push_back(pos);  // coords[0]

    double direction = 0.0;

    for (size_t i = 0; i < d_ringSize; ++i) {
      // Apply turn at atom i with extra angle
      if (d_turns[i] == 1) {
        direction += (M_PI / 3.0) + extraPerAtom;  // R turn with extra
      } else if (d_turns[i] == -1) {
        direction += -(M_PI / 3.0) + extraPerAtom;  // L turn with extra
      }

      // Move in current direction to reach next atom
      pos.x += d_bondLength * std::cos(direction);
      pos.y += d_bondLength * std::sin(direction);
      coords.push_back(pos);  // Add atoms 1 to N (where N is the dummy)
    }
    for (size_t i = 0; i < coords.size(); ++i) {
      std::cerr << "ODD: Atom " << i << ": (" << coords[i].x << ", "
                << coords[i].y << ")" << std::endl;
    }

    // Now coords has N+1 elements: [0, 1, 2, ..., N-1, N_dummy]
  } else {
    // Even rings: use standard hexagonal grid, then add dummy
    auto hexCoords = calculateHexCoords(d_turns);

    // Convert hex coordinates to 2D Cartesian coordinates
    const double SQRT3_2 = std::sqrt(3.0) / 2.0;

    for (size_t i = 0; i < d_ringSize; ++i) {
      const auto &hex = hexCoords[i];
      double x = d_bondLength * (hex.x + hex.z / 2.0);
      double y = d_bondLength * (-SQRT3_2 * hex.z);
      coords.emplace_back(x, y);
    }
    for (size_t i = 0; i < coords.size(); ++i) {
      std::cerr << "EVEN: Atom " << i << ": (" << coords[i].x << ", "
                << coords[i].y << ")" << std::endl;
    }
  }
  // Analytical refinement: adjust angles with Jacobian (if enabled)
  if (d_useJacobianRefinement) {
    adjustAnglesForClosure(coords, isOddRing);
  }

  // Remove the dummy atom (for both even and odd rings)
  if (coords.size() == d_ringSize + 1) {
    coords.pop_back();
  }

  // Center the coordinates (translate centroid to origin)
  RDGeom::Point2D centroid(0.0, 0.0);
  for (const auto &coord : coords) {
    centroid.x += coord.x;
    centroid.y += coord.y;
  }
  centroid.x /= d_ringSize;
  centroid.y /= d_ringSize;

  for (auto &coord : coords) {
    coord.x -= centroid.x;
    coord.y -= centroid.y;
  }

  return coords;
}

void MacrocycleGenerator::adjustAnglesForClosure(
    std::vector<RDGeom::Point2D> &coords, bool isOddRing) const {
  // Iterative Jacobian pseudo-inverse with lever-arm based derivatives
  // coords has N+1 elements: [0, 1, ..., N-1, N_dummy]

  size_t N = d_ringSize;
  double initialGap = (coords[0] - coords[N]).length();

  std::cerr << "DEBUG Angle Adjustment (Iterative Jacobian): Initial gap = "
            << initialGap << " Å" << std::endl;

  if (initialGap < 1e-6) {
    std::cerr << "  Gap already closed" << std::endl;
    return;
  }

  // Extra angle per atom for odd rings
  double extraPerAtom = isOddRing ? (M_PI / 3.0) / N : 0.0;

  // Cumulative angle adjustments
  std::vector<double> angleAdjustments(N, 0.0);

  // Iterative refinement with damping
  const int maxIters = 3;
  const double damping = 0.8;  // Apply 80% of computed correction

  for (int iter = 0; iter < maxIters; ++iter) {
    // Rebuild coordinates with current adjustments
    RDGeom::Point2D pos(0, 0);
    coords[0] = pos;

    double direction = 0.0;

    for (size_t i = 0; i < N; ++i) {
      if (d_turns[i] == 1) {
        direction += (M_PI / 3.0) + extraPerAtom + angleAdjustments[i];
      } else if (d_turns[i] == -1) {
        direction += -(M_PI / 3.0) + extraPerAtom + angleAdjustments[i];
      }

      pos.x += d_bondLength * std::cos(direction);
      pos.y += d_bondLength * std::sin(direction);
      coords[i + 1] = pos;
    }

    // Measure gap
    RDGeom::Point2D gap = coords[0] - coords[N];
    double gapMag = gap.length();

    std::cerr << "  Iter " << iter << ": gap = " << gapMag << " Å";

    if (gapMag < 0.01) {
      std::cerr << " - Converged!" << std::endl;
      break;
    }

    // Build Jacobian: ∂gap/∂θ_k = perpendicular(p_N - p_k)
    std::vector<RDGeom::Point2D> jacobian(N);

    for (size_t k = 0; k < N; ++k) {
      RDGeom::Point2D r = coords[N] - coords[k];
      jacobian[k].x = r.y;   // ∂gap_x/∂θ_k
      jacobian[k].y = -r.x;  // ∂gap_y/∂θ_k
    }

    // Compute J·J^T (2×2)
    double A00 = 0.0, A01 = 0.0, A11 = 0.0;
    for (size_t k = 0; k < N; ++k) {
      A00 += jacobian[k].x * jacobian[k].x;
      A01 += jacobian[k].x * jacobian[k].y;
      A11 += jacobian[k].y * jacobian[k].y;
    }

    // Invert 2×2 matrix
    double det = A00 * A11 - A01 * A01;

    if (std::abs(det) < 1e-12) {
      std::cerr << " - Singular Jacobian!" << std::endl;
      break;
    }

    double invA00 = A11 / det;
    double invA01 = -A01 / det;
    double invA11 = A00 / det;

    // Compute J^T · (J·J^T)^{-1} · gap
    double temp_x = invA00 * gap.x + invA01 * gap.y;
    double temp_y = invA01 * gap.x + invA11 * gap.y;

    // Compute Δθ = -J^T · temp
    std::vector<double> deltaTheta(N);
    double maxDelta = 0.0;

    for (size_t k = 0; k < N; ++k) {
      deltaTheta[k] = -(jacobian[k].x * temp_x + jacobian[k].y * temp_y);
      maxDelta = std::max(maxDelta, std::abs(deltaTheta[k]));
    }

    std::cerr << ", max Δθ = " << (maxDelta * 180.0 / M_PI) << "°" << std::endl;

    // Apply damped adjustments
    for (size_t k = 0; k < N; ++k) {
      angleAdjustments[k] += damping * deltaTheta[k];
    }
  }

  // Final gap after adjustments
  RDGeom::Point2D finalGap = coords[0] - coords[N];
  double finalGapMag = finalGap.length();
  std::cerr << "  Final gap: " << finalGapMag << " Å" << std::endl;

  // Debug: Print turn sequence vs actual angles
  std::cerr << "  Turn sequence: ";
  for (size_t i = 0; i < N; ++i) {
    if (d_turns[i] == 1)
      std::cerr << "R";
    else if (d_turns[i] == -1)
      std::cerr << "L";
    else
      std::cerr << ".";
  }
  std::cerr << std::endl;

  // Measure actual angles at positions 0 to N-1
  std::cerr << "  Actual exterior angles: ";
  for (size_t i = 0; i < N; ++i) {
    size_t prev = (i == 0) ? N : (i - 1);
    size_t next = i + 1;

    RDGeom::Point2D v1 = coords[i] - coords[prev];
    RDGeom::Point2D v2 = coords[next] - coords[i];

    double angle1 = std::atan2(v1.y, v1.x);
    double angle2 = std::atan2(v2.y, v2.x);
    double exteriorTurn = angle2 - angle1;

    while (exteriorTurn > M_PI) exteriorTurn -= 2 * M_PI;
    while (exteriorTurn < -M_PI) exteriorTurn += 2 * M_PI;

    double exteriorDeg = exteriorTurn * 180 / M_PI;
    std::cerr << (exteriorDeg > 0 ? "R" : "L") << "(" << std::abs(exteriorDeg)
              << ") ";
  }
  std::cerr << std::endl;

  // Check for inversions
  std::cerr << "  Checking angles after adjustment:" << std::endl;
  int inverted = 0;
  int largeDev = 0;
  for (size_t i = 0; i < N; ++i) {
    size_t prev = (i == 0) ? N : (i - 1);
    size_t next = i + 1;

    RDGeom::Point2D v1 = coords[i] - coords[prev];
    RDGeom::Point2D v2 = coords[next] - coords[i];

    double angle1 = std::atan2(v1.y, v1.x);
    double angle2 = std::atan2(v2.y, v2.x);
    double exteriorTurn = angle2 - angle1;

    while (exteriorTurn > M_PI) exteriorTurn -= 2 * M_PI;
    while (exteriorTurn < -M_PI) exteriorTurn += 2 * M_PI;

    double exteriorDeg = exteriorTurn * 180 / M_PI;

    size_t turnIdx = i;

    double expectedExterior = 0.0;
    if (d_turns[turnIdx] == 1) {
      expectedExterior = 60.0 + (isOddRing ? 4.0 : 0.0);
    } else if (d_turns[turnIdx] == -1) {
      expectedExterior = -60.0 + (isOddRing ? 4.0 : 0.0);
    }

    double deviation = std::abs(exteriorDeg - expectedExterior);

    bool isRTurn = (d_turns[turnIdx] == 1);
    bool hasCorrectSign =
        (isRTurn && exteriorTurn > 0) || (!isRTurn && exteriorTurn < 0);

    if (!hasCorrectSign && d_turns[turnIdx] != 0) {
      inverted++;
      std::cerr << "    Pos " << i << " (turn[" << turnIdx
                << "]=" << (isRTurn ? "R" : "L") << ") INVERTED! Got "
                << exteriorDeg << "°" << std::endl;
    } else if (deviation > 10.0 && d_turns[turnIdx] != 0) {
      largeDev++;
      std::cerr << "    Pos " << i << " (turn[" << turnIdx
                << "]=" << (isRTurn ? "R" : "L") << "): " << exteriorDeg
                << "° (expected " << expectedExterior << "°, dev=" << deviation
                << "°)" << std::endl;
    }
  }

  if (inverted == 0 && largeDev == 0) {
    std::cerr << "  ✓ All angles within tolerance" << std::endl;
  } else {
    std::cerr << "  INVERSIONS: " << inverted << "/" << N << std::endl;
    std::cerr << "  LARGE DEVIATIONS (>10°): " << largeDev << "/" << N
              << std::endl;
  }
}

}  // namespace RDDepict
