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
#include "DepictUtils.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <set>
#include <map>
#include <unordered_set>

namespace RDDepict {

MacrocycleGenerator::MacrocycleGenerator(size_t ringSize, double bondLength)
    : d_ringSize(ringSize),
      d_bondLength(bondLength),
      d_innerTurnSign(0),
      d_closureError(std::numeric_limits<double>::max()),
      d_solved(false) {
  d_turns.resize(ringSize, 0);
}

void MacrocycleGenerator::addConstraint(const TurnConstraint &constraint) {
  d_constraints.push_back(constraint);
}

void MacrocycleGenerator::addAngleConstraint(
    const AngleConstraint &constraint) {
  // Check if position already has a constraint
  for (auto &existing : d_angleConstraints) {
    if (existing.position == constraint.position) {
      // Keep the constraint with smallest absolute angle (most constrained)
      if (std::abs(constraint.targetAngle) < std::abs(existing.targetAngle)) {
        existing = constraint;
      }
      return;
    }
  }
  d_angleConstraints.push_back(constraint);
}

void MacrocycleGenerator::setSubstituentInfo(
    const std::map<size_t, int> &substituentSizes, int innerTurnSign) {
  d_substituentSizes = substituentSizes;
  d_innerTurnSign = innerTurnSign;
}

bool MacrocycleGenerator::hasAngleConstraint(size_t position) const {
  for (const auto &constraint : d_angleConstraints) {
    if (constraint.position == position) {
      return true;
    }
  }
  return false;
}

double MacrocycleGenerator::getConstraintAngle(size_t position) const {
  for (const auto &constraint : d_angleConstraints) {
    if (constraint.position == position) {
      return constraint.targetAngle;
    }
  }
  return 0.0;  // No constraint
}

bool MacrocycleGenerator::applyConstraints(std::vector<int> &turns) const {
  for (const auto &constraint : d_constraints) {
    size_t pos = constraint.position;

    if (constraint.type == ConstraintType::FIXED) {
      // Apply fixed pattern
      for (size_t i = 0; i < constraint.pattern.size(); ++i) {
        auto idx = (pos + i) % d_ringSize;
        if (turns[idx] == 0) {
          turns[idx] = constraint.pattern[i];
        } else if (turns[idx] != constraint.pattern[i]) {
          // Conflict: position already set to different value
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
    return false;  // Constraints conflict
  }

  bool isOddRing = (d_ringSize % 2) == 1;

  // Try different target differences with relaxation
  // Normal: R-L = 6 (even) or 5 (or 7) (odd)
  std::vector<int> targetDiffs;

  if (isOddRing) {
    targetDiffs = {5, 7, 3};
  } else {
    targetDiffs = {6, 4};
  }

  for (int targetDiff : targetDiffs) {
    if (trySolveWithTargetDiff(targetDiff)) {
      return true;
    }
  }

  return false;
}

bool MacrocycleGenerator::trySolveWithTargetDiff(int targetDiff) {
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

size_t MacrocycleGenerator::getNumFreePositions() const {
  return std::count(d_turns.begin(), d_turns.end(), 0);
}

MacrocycleGenerator::FastHeuristicResult MacrocycleGenerator::tryFastHeuristic(
    std::vector<size_t> &freePositions, int &numRight, int &numLeft) {
  size_t numFree = freePositions.size();

  // Strategy: Use substituent info to reduce degrees of freedom
  // Positions with big substituents (size > 1) should be outer turns
  // This is both chemically sensible and reduces the search space
  std::vector<int> candidate = d_turns;  // Start with constrained values

  int outerTurn = -d_innerTurnSign;  // Outer turn is opposite of inner
  int assignedRight = 0;
  int assignedLeft = 0;

  // Step 1: Assign outer turns for positions with big substituents (size > 1)
  int substituentAssignments = 0;
  for (size_t pos : freePositions) {
    auto it = d_substituentSizes.find(pos);
    if (it != d_substituentSizes.end() && it->second > 1) {
      // Big substituent - force outer turn
      candidate[pos] = outerTurn;
      if (outerTurn == 1) {
        assignedRight++;
      } else {
        assignedLeft++;
      }
      substituentAssignments++;
    }
  }

  // Step 2: Check if we have capacity for these assignments
  int remainingRight = numRight - assignedRight;
  int remainingLeft = numLeft - assignedLeft;

  if (remainingRight < 0 || remainingLeft < 0) {
    return FastHeuristicResult::FAILURE;
  }

  // Step 3: If still too many free positions, use alternating for most, keep
  // few truly free
  std::vector<size_t> stillFree;
  for (size_t pos : freePositions) {
    if (candidate[pos] == 0) {
      stillFree.push_back(pos);
    }
  }

  const size_t MAX_ENUMERATION_SIZE = 10;  // Only enumerate over ~10 positions

  if (stillFree.size() > MAX_ENUMERATION_SIZE) {
    // How many to fix vs keep free
    size_t numToFix = stillFree.size() - MAX_ENUMERATION_SIZE;

    // Use alternating pattern for first numToFix positions
    int currentTurn = (remainingRight >= remainingLeft) ? 1 : -1;
    size_t fixed = 0;

    for (size_t pos : stillFree) {
      if (fixed >= numToFix) break;  // Leave the rest free

      // Assign with alternating pattern
      if (currentTurn == 1 && remainingRight > 0) {
        candidate[pos] = 1;
        remainingRight--;
        currentTurn = -1;
      } else if (currentTurn == -1 && remainingLeft > 0) {
        candidate[pos] = -1;
        remainingLeft--;
        currentTurn = 1;
      } else if (remainingRight > 0) {
        candidate[pos] = 1;
        remainingRight--;
      } else if (remainingLeft > 0) {
        candidate[pos] = -1;
        remainingLeft--;
      }
      fixed++;
    }
    // Update d_turns with the fixed positions, then fall through to normal
    // enumeration
    d_turns = candidate;

    // Recalculate free positions for enumeration
    freePositions.clear();
    for (size_t i = 0; i < d_ringSize; ++i) {
      if (d_turns[i] == 0) {
        freePositions.push_back(i);
      }
    }
    // Update numRight and numLeft to reflect the remaining assignments needed
    numRight = remainingRight;
    numLeft = remainingLeft;

    // Fall through to enumeration with reduced free positions
    return FastHeuristicResult::CONTINUE;

  } else {
    // Few enough positions - just assign them all with alternating
    int currentTurn = (remainingRight >= remainingLeft) ? 1 : -1;

    for (size_t pos : stillFree) {
      // Assign based on what we still need
      if (currentTurn == 1 && remainingRight > 0) {
        candidate[pos] = 1;
        remainingRight--;
        currentTurn = -1;  // Try to alternate
      } else if (currentTurn == -1 && remainingLeft > 0) {
        candidate[pos] = -1;
        remainingLeft--;
        currentTurn = 1;  // Try to alternate
      } else if (remainingRight > 0) {
        candidate[pos] = 1;
        remainingRight--;
      } else if (remainingLeft > 0) {
        candidate[pos] = -1;
        remainingLeft--;
      }
    }

    // Validate and return (all positions assigned)
    if (validateConstraints(candidate)) {
      int minDistance = calculateMinDistance(candidate);
      if (minDistance > 0) {
        d_turns = candidate;
        d_closureError = calculateClosureError(candidate);
        return FastHeuristicResult::SUCCESS;
      }
    }
    // All positions assigned but failed - can't enumerate further
    return FastHeuristicResult::FAILURE;
  }
}

bool MacrocycleGenerator::findOptimalTurnSequence(int numRight, int numLeft) {
  // Find positions of free variables
  std::vector<size_t> freePositions;
  for (size_t i = 0; i < d_ringSize; ++i) {
    if (d_turns[i] == 0) {
      freePositions.push_back(i);
    }
  }

  const size_t ENUMERATION_THRESHOLD = 15;

  if (freePositions.size() > ENUMERATION_THRESHOLD) {
    // System is too complex for exhaustive search, try a fast heuristic to
    // reduce the search space
    FastHeuristicResult heuristicResult =
        tryFastHeuristic(freePositions, numRight, numLeft);
    if (heuristicResult == FastHeuristicResult::SUCCESS) {
      return true;
    } else if (heuristicResult == FastHeuristicResult::FAILURE) {
      return false;
    }
    // Otherwise CONTINUE with enumeration, the problem has been reduced and can
    // now be tackled with a exhaustive search
  }

  // If no free variables, we're done (constraints fully determined the
  // solution)
  if (freePositions.empty()) {
    d_closureError = calculateClosureError(d_turns);
    return true;
  }

  // Enumerate ways to assign R/L to free positions
  std::vector<size_t> rightPositions(numRight);
  size_t numFree = freePositions.size();

  // Initialize first combination: choose first numRight free positions for R
  for (int i = 0; i < numRight; ++i) {
    rightPositions[i] = i;
  }

  double bestScore = std::numeric_limits<double>::max();
  std::vector<int> bestTurns;

  // Helper: create candidate from current combination
  auto createCandidate = [&]() {
    std::vector<int> candidate = d_turns;  // Start with constrained values
    // Set free positions to L
    for (size_t freeIdx : freePositions) {
      candidate[freeIdx] = -1;
    }
    // Set selected positions to R
    for (size_t idx : rightPositions) {
      candidate[freePositions[idx]] = 1;
    }
    return candidate;
  };

  // Weight for minimum distance in scoring
  const double MIN_DISTANCE_WEIGHT = 0.5;

  // Helper: evaluate a candidate and update best solution if better
  auto evaluateCandidate = [&](const std::vector<int> &candidate) {
    if (validateConstraints(candidate)) {
      int minDistance = calculateMinDistance(candidate);
      // Skip if self-crossing detected
      if (minDistance > 0) {
        double closureError = calculateClosureError(candidate);
        double substituentPenalty = calculateSubstituentPenalty(candidate);
        double score = closureError - MIN_DISTANCE_WEIGHT * minDistance +
                       substituentPenalty;
        if (score < bestScore) {
          bestScore = score;
          bestTurns = candidate;
        }
      }
    }
  };

  // Enumerate all combinations
  std::vector<int> candidate;
  do {
    candidate = createCandidate();
    evaluateCandidate(candidate);
  } while (nextCombination(rightPositions, numFree));

  // Store best solution
  if (bestTurns.empty()) {
    return false;  // No valid solution found that satisfies constraints
  }

  d_turns = bestTurns;
  d_closureError = calculateClosureError(bestTurns);
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

int MacrocycleGenerator::calculateMinDistance(
    const std::vector<int> &turns) const {
  // Calculate hex coordinates for all atoms
  auto coords = calculateHexCoords(turns);

  // Find minimum distance between non-adjacent atoms
  // Distance in hex cube coordinates: (|dx| + |dy| + |dz|) / 2
  // Also detect clashes (distance = 0, self-crossing)
  int minDist = std::numeric_limits<int>::max();

  for (size_t i = 0; i < coords.size() - 1; ++i) {  // -1 to exclude dummy
    for (size_t j = i + 1; j < coords.size() - 1; ++j) {
      // Skip adjacent atoms (bonded atoms always have distance 1)
      if (j == i + 1 || (i == 0 && j == coords.size() - 2)) {
        continue;
      }

      int dx = std::abs(coords[i].x - coords[j].x);
      int dy = std::abs(coords[i].y - coords[j].y);
      int dz = std::abs(coords[i].z - coords[j].z);
      int dist = (dx + dy + dz) / 2;

      // Early exit on clash (self-crossing)
      if (dist == 0) {
        return 0;
      }

      if (dist < minDist) {
        minDist = dist;
      }
    }
  }

  return minDist;
}

double MacrocycleGenerator::calculateSubstituentPenalty(
    const std::vector<int> &turns) const {
  if (d_substituentSizes.empty() || d_innerTurnSign == 0) {
    return 0.0;  // No substituents or inner turn not set
  }

  // Calculate penalty for substituents on inner turns
  // Penalty is weighted by substituent size
  double penalty = 0.0;
  const double PENALTY_PER_ATOM =
      1.0;  // Penalty per substituent atom on inner turn

  for (const auto &[pos, size] : d_substituentSizes) {
    if (pos < turns.size()) {
      int turnAtPos = turns[pos];
      bool isInner = (turnAtPos == d_innerTurnSign);
      if (isInner) {
        penalty += size * PENALTY_PER_ATOM;
      }
    }
  }
  return penalty;
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
    // Count constrained angles
    size_t numConstrained = d_angleConstraints.size();
    size_t numFree = d_ringSize - numConstrained;

    // For odd rings, distribute extra angle only across free angles
    double extraPerAtom = 0.0;
    if (numFree > 0) {
      extraPerAtom = (M_PI / 3.0) / numFree;
    }

    RDGeom::Point2D pos(0, 0);
    coords.push_back(pos);  // coords[0]

    double direction = 0.0;

    for (size_t i = 0; i < d_ringSize; ++i) {
      double turnAngle;

      // Check for angle constraint at this position
      if (hasAngleConstraint(i)) {
        turnAngle =
            getConstraintAngle(i);  // Use exact ideal angle, NO extraPerAtom
      } else {
        // Normal turn computation with extraPerAtom (only for free angles)
        if (d_turns[i] == 1) {
          turnAngle = (M_PI / 3.0) + extraPerAtom;  // R turn with extra
        } else if (d_turns[i] == -1) {
          turnAngle = -(M_PI / 3.0) + extraPerAtom;  // L turn with extra
        } else {
          turnAngle = 0.0;  // No turn
        }
      }

      direction += turnAngle;

      // Move in current direction to reach next atom
      pos.x += d_bondLength * std::cos(direction);
      pos.y += d_bondLength * std::sin(direction);
      coords.push_back(pos);  // Add atoms 1 to N (where N is the dummy)
    }

    // Now coords has N+1 elements: [0, 1, 2, ..., N-1, N_dummy]
  } else {
    // Even rings: use standard hexagonal grid, then add dummy
    auto hexCoords = calculateHexCoords(d_turns);

    // Convert hex coordinates to 2D Cartesian coordinates
    const double SQRT3_2 = std::sqrt(3.0) / 2.0;

    // Convert first N coordinates (atoms 0 to N-1)
    for (size_t i = 0; i < d_ringSize; ++i) {
      const auto &hex = hexCoords[i];
      double x = d_bondLength * (hex.x + hex.z / 2.0);
      double y = d_bondLength * (-SQRT3_2 * hex.z);
      coords.emplace_back(x, y);
    }

    // Add dummy atom (position N) - should be close to atom 0 for even rings
    const auto &dummyHex = hexCoords[d_ringSize];
    double dummyX = d_bondLength * (dummyHex.x + dummyHex.z / 2.0);
    double dummyY = d_bondLength * (-SQRT3_2 * dummyHex.z);
    coords.emplace_back(dummyX, dummyY);

    // Now coords has N+1 elements: [0, 1, 2, ..., N-1, N_dummy]
  }
  // Analytical refinement: adjust angles with Jacobian
  adjustAnglesForClosure(coords, isOddRing);

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

  if (initialGap < 1e-6) {
    return;
  }

  // Extra angle per atom for odd rings
  // Count constrained angles and distribute extraPerAtom only across free
  // angles
  size_t numConstrained = d_angleConstraints.size();
  size_t numFree = N - numConstrained;
  double extraPerAtom = 0.0;
  if (isOddRing && numFree > 0) {
    extraPerAtom = (M_PI / 3.0) / numFree;
  }

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
      double turnAngle;

      // Check for angle constraint at this position
      if (hasAngleConstraint(i)) {
        turnAngle = getConstraintAngle(i);  // Use exact ideal angle
      } else {
        // Normal turn computation with extraPerAtom and adjustments
        if (d_turns[i] == 1) {
          turnAngle = (M_PI / 3.0) + extraPerAtom + angleAdjustments[i];
        } else if (d_turns[i] == -1) {
          turnAngle = -(M_PI / 3.0) + extraPerAtom + angleAdjustments[i];
        } else {
          turnAngle = 0.0;
        }
      }

      direction += turnAngle;

      pos.x += d_bondLength * std::cos(direction);
      pos.y += d_bondLength * std::sin(direction);
      coords[i + 1] = pos;
    }

    // Measure gap
    RDGeom::Point2D gap = coords[0] - coords[N];
    double gapMag = gap.length();

    if (gapMag < 0.01) {
      break;
    }

    // Build Jacobian: ∂gap/∂θ_k = perpendicular(p_N - p_k)
    std::vector<RDGeom::Point2D> jacobian(N);

    for (size_t k = 0; k < N; ++k) {
      // Skip constrained angles - zero out their Jacobian row
      if (hasAngleConstraint(k)) {
        jacobian[k].x = 0.0;
        jacobian[k].y = 0.0;
        continue;
      }

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

    for (size_t k = 0; k < N; ++k) {
      deltaTheta[k] = -(jacobian[k].x * temp_x + jacobian[k].y * temp_y);
    }

    // Apply damped adjustments
    for (size_t k = 0; k < N; ++k) {
      // Skip constrained angles
      if (hasAngleConstraint(k)) {
        continue;
      }
      angleAdjustments[k] += damping * deltaTheta[k];
    }
  }

  // Final gap measurement
  double finalGap = (coords[0] - coords[N]).length();
}

// ============================================================================
// Helper Functions for Angle Constraint Detection and Refinement
// ============================================================================

double computeIdealAngle(int ringSize) {
  return 2.0 * M_PI / ringSize;  // 360/n degrees
}

// Find positions in the macrocycle ring that correspond to atoms shared with
// the given ring
std::vector<size_t> findSharedPositions(const std::vector<int> &macrocycleRing,
                                        const std::vector<int> &ring) {
  // Build set of ring atoms for O(1) lookup instead of O(N) linear search
  std::unordered_set<int> ringAtoms(ring.begin(), ring.end());

  std::vector<size_t> sharedPositions;
  for (size_t i = 0; i < macrocycleRing.size(); ++i) {
    int atom = macrocycleRing[i];
    if (ringAtoms.count(atom)) {
      sharedPositions.push_back(i);
    }
  }
  return sharedPositions;
}

// Verify that shared positions are contiguous (allowing wraparound) and reorder
// them so they appear in the order they are encountered in the macrocycle.
// E.g., if shared positions are [0, 1, 14] in a 15-atom macrocycle, this is
// contiguous with wraparound and will be reordered to [14, 0, 1] so that
// constraint position calculation works correctly.
bool verifyAndReorderSharedPositions(std::vector<size_t> &sharedPositions,
                                     size_t macrocycleSize) {
  if (sharedPositions.size() <= 1) {
    return true;  // Single position is always contiguous
  }

  bool contiguous = true;
  size_t gapIdx = 0;  // Track where wraparound gap is (if any)

  // Check if all consecutive pairs are adjacent
  for (size_t i = 1; i < sharedPositions.size(); ++i) {
    size_t expected = (sharedPositions[i - 1] + 1) % macrocycleSize;
    if (sharedPositions[i] != expected) {
      contiguous = false;
      break;
    }
  }

  // If not contiguous in forward direction, check for wraparound
  // E.g., positions [0, 1, 14] or [13, 14, 0, 1]
  if (!contiguous) {
    // Find the gap: where consecutive positions have a large jump
    gapIdx = 0;
    for (size_t i = 1; i < sharedPositions.size(); ++i) {
      if (sharedPositions[i] != sharedPositions[i - 1] + 1) {
        gapIdx = i;
        break;
      }
    }

    // If there's a gap, check if it's the wraparound
    if (gapIdx > 0) {
      // Check if positions after gap are contiguous
      bool wrapContiguous = true;

      for (size_t i = gapIdx + 1; i < sharedPositions.size(); ++i) {
        if (sharedPositions[i] != sharedPositions[i - 1] + 1) {
          wrapContiguous = false;
          break;
        }
      }

      // Check wraparound: last position + 1 should equal first position (mod
      // ring size)
      if (wrapContiguous) {
        size_t lastPos = sharedPositions[sharedPositions.size() - 1];
        size_t firstPos = sharedPositions[0];
        if ((lastPos + 1) % macrocycleSize == firstPos) {
          contiguous = true;
        }
      }
    }
  }

  if (!contiguous) {
    return false;
  }

  // If wraparound was detected, reorder positions to start from the wraparound
  // point E.g., [0, 1, 14] → [14, 0, 1] so constraint position calculation
  // works correctly
  if (gapIdx > 0) {
    // Rotate the vector: move elements from gapIdx to end to the front
    std::vector<size_t> reordered;
    for (size_t i = gapIdx; i < sharedPositions.size(); ++i) {
      reordered.push_back(sharedPositions[i]);
    }
    for (size_t i = 0; i < gapIdx; ++i) {
      reordered.push_back(sharedPositions[i]);
    }
    sharedPositions = reordered;
  }
  return true;
}

void trackEndpoint(
    std::map<size_t, std::vector<EndpointInfo>> &endpointPositions, size_t pos,
    int ringSize, const std::vector<int> &ringAtoms) {
  EndpointInfo info;
  info.ringSize = ringSize;
  info.ringAtoms = ringAtoms;
  endpointPositions[pos].push_back(info);
}

AngleConstraint computeSharedEndpointConstraint(
    const std::vector<EndpointInfo> &endpoints, size_t pos) {
  AngleConstraint constraint;
  constraint.position = SIZE_MAX;  // Invalid by default

  if (endpoints.size() != 2) {
    return constraint;  // Only valid for exactly 2 endpoints
  }
  int ringSize1 = endpoints[0].ringSize;
  int ringSize2 = endpoints[1].ringSize;
  const auto &ring1Atoms = endpoints[0].ringAtoms;
  const auto &ring2Atoms = endpoints[1].ringAtoms;

  // Count shared atoms between the two small rings
  std::set<int> ring1Set(ring1Atoms.begin(), ring1Atoms.end());
  size_t sharedCount = 0;
  for (int atom : ring2Atoms) {
    if (ring1Set.count(atom)) {
      ++sharedCount;
    }
  }

  // Get the ideal turn magnitude for each ring
  double turn1 = computeIdealAngle(ringSize1);
  double turn2 = computeIdealAngle(ringSize2);

  double combinedAngle;
  if (sharedCount == 1) {
    // Spiro fusion: the two rings are spiro-connected at this macrocycle
    // position The total angle is half of the remaining space, so the
    // connection is symmetrical

    double macrocycleInternal = (turn1 + turn2) / 2.0;
    combinedAngle = M_PI - macrocycleInternal;
  } else if (sharedCount == 2) {
    // Edge fusion (2 shared atoms): the two rings share an edge
    // At a shared endpoint, three rings meet: ring1, ring2, and the macrocycle
    // The sum of their internal angles must be 2π (360°)
    // Internal angle = π - |turn_angle|
    // So: (π - |turn1|) + (π - |turn2|) + macrocycle_internal = 2π
    // Therefore: macrocycle_internal = 2π - (π - |turn1|) - (π - |turn2|)
    //                                 = |turn1| + |turn2|
    // And macrocycle_turn = π - macrocycle_internal = π - (|turn1| + |turn2|)
    double macrocycleInternal = turn1 + turn2;
    combinedAngle = M_PI - macrocycleInternal;
  }

  constraint.position = pos;
  constraint.targetAngle = combinedAngle;
  return constraint;
}

// add constraints for fused rings. Each angle formed by 3 atoms shared between
// the macrocycle and a fused ring can be constrained to the ideal angle for
// that small ring (e.g. 120° for a 6-membered ring, 108° for 5-membered ring
// etc. When an atom is fused between the macrocycle and two other rings, its
// target angle becomes 360* - the sum of the two small ring angles, so that
// both can remain regular polygons)
std::vector<AngleConstraint> identifyAngleConstraintsForFusedRings(
    const std::vector<int> &macrocycleRing,
    const std::vector<std::vector<int>> &allRings) {
  std::vector<AngleConstraint> angleConstraints;
  std::map<size_t, std::vector<EndpointInfo>> endpointPositions;

  // Process each ring to find fusions with the macrocycle
  for (const auto &ring : allRings) {
    if (ring.size() > MACROCYCLE_SIZE_THRESHOLD) {
      continue;
    }
    // Find shared atoms between this ring and the macrocycle
    std::vector<size_t> sharedPositions =
        findSharedPositions(macrocycleRing, ring);

    // Need at least 1 shared atom (spiro or fusion)
    // Skip if more than 4 shared atoms (unusual/complex fusion)
    if (sharedPositions.size() < 1 || sharedPositions.size() > 4) {
      continue;
    }

    // Verify shared atoms are contiguous and reorder if needed
    if (!verifyAndReorderSharedPositions(sharedPositions,
                                         macrocycleRing.size())) {
      continue;
    }
    size_t numShared = sharedPositions.size();

    // Track endpoints for this fused ring
    size_t firstPos = sharedPositions[0];
    size_t lastPos = sharedPositions[numShared - 1];

    trackEndpoint(endpointPositions, firstPos, ring.size(), ring);
    trackEndpoint(endpointPositions, lastPos, ring.size(), ring);

    // Add constraints for internal positions (only when numShared >= 3)
    // Loop runs 0 times when numShared == 2 (since 1 < 2-1 is false)
    for (size_t i = 1; i < numShared - 1; ++i) {
      size_t pos = sharedPositions[i];
      int turnSign = -1;  // L turn
      double idealAngle = computeIdealAngle(ring.size()) * turnSign;

      AngleConstraint angleConstraint;
      angleConstraint.position = pos;
      angleConstraint.targetAngle = idealAngle;
      angleConstraints.push_back(angleConstraint);
    }
  }
  // Process shared endpoints: positions that are endpoints of exactly 2 small
  // rings
  for (const auto &[pos, endpoints] : endpointPositions) {
    if (endpoints.size() == 2) {
      AngleConstraint sharedConstraint =
          computeSharedEndpointConstraint(endpoints, pos);
      if (sharedConstraint.position != SIZE_MAX) {  // Valid constraint
        angleConstraints.push_back(sharedConstraint);
      }
    }
  }
  return angleConstraints;
}

void refineWithJacobian(std::vector<RDGeom::Point2D> &coords,
                        std::vector<double> &angles,
                        const std::set<size_t> &constrainedPositions,
                        double bondLength, bool isOddRing) {
  // Iterative Jacobian pseudo-inverse refinement to close the gap
  // coords has N+1 elements: [0, 1, ..., N-1, N_dummy]

  size_t N = coords.size() - 1;  // Last element is dummy atom
  double initialGap = (coords[0] - coords[N]).length();

  const double gapThreshold = 0.01;  // Stop if gap < 0.01 Å

  if (initialGap < gapThreshold) {
    return;  // Already closed
  }

  // Cumulative angle adjustments
  std::vector<double> angleAdjustments(N, 0.0);

  // Iterative refinement with damping
  const int maxIters = 4;
  const double damping = 0.8;  // Apply 80% of computed correction

  for (int iter = 0; iter < maxIters; ++iter) {
    // Rebuild coordinates with current angle adjustments
    RDGeom::Point2D pos(0, 0);
    coords[0] = pos;

    double direction = 0.0;

    for (size_t i = 0; i < N; ++i) {
      double turnAngle;

      // Check if this position is constrained
      if (constrainedPositions.count(i) > 0) {
        turnAngle = angles[i];  // Use exact constrained angle
      } else {
        turnAngle = angles[i] + angleAdjustments[i];  // Apply adjustment
      }

      direction += turnAngle;

      pos.x += bondLength * std::cos(direction);
      pos.y += bondLength * std::sin(direction);
      coords[i + 1] = pos;
    }

    // Measure gap
    RDGeom::Point2D gap = coords[0] - coords[N];
    double gapMag = gap.length();

    if (gapMag < gapThreshold) {
      break;
    }

    // Build Jacobian: ∂gap/∂θ_k = perpendicular(p_N - p_k)
    std::vector<RDGeom::Point2D> jacobian(N);

    for (size_t k = 0; k < N; ++k) {
      // Skip constrained angles - zero out their Jacobian row
      if (constrainedPositions.count(k) > 0) {
        jacobian[k].x = 0.0;
        jacobian[k].y = 0.0;
        continue;
      }

      RDGeom::Point2D r = coords[N] - coords[k];
      jacobian[k].x = r.y;   // ∂gap_x/∂θ_k (perpendicular)
      jacobian[k].y = -r.x;  // ∂gap_y/∂θ_k
    }

    // Compute J·J^T (2×2 matrix)
    double A00 = 0.0, A01 = 0.0, A11 = 0.0;
    for (size_t k = 0; k < N; ++k) {
      A00 += jacobian[k].x * jacobian[k].x;
      A01 += jacobian[k].x * jacobian[k].y;
      A11 += jacobian[k].y * jacobian[k].y;
    }

    // Invert 2×2 matrix
    double det = A00 * A11 - A01 * A01;

    if (std::abs(det) < 1e-12) {
      break;  // Singular matrix, cannot proceed
    }

    double invA00 = A11 / det;
    double invA01 = -A01 / det;
    double invA11 = A00 / det;

    // Compute J^T · (J·J^T)^{-1} · gap
    double temp_x = invA00 * gap.x + invA01 * gap.y;
    double temp_y = invA01 * gap.x + invA11 * gap.y;

    // Compute Δθ = -J^T · temp
    std::vector<double> deltaTheta(N);

    for (size_t k = 0; k < N; ++k) {
      deltaTheta[k] = -(jacobian[k].x * temp_x + jacobian[k].y * temp_y);
    }

    // Apply damped adjustments (only to non-constrained angles)
    for (size_t k = 0; k < N; ++k) {
      if (constrainedPositions.count(k) == 0) {
        angleAdjustments[k] += damping * deltaTheta[k];
      }
    }
  }

  // Update angles with final adjustments
  for (size_t i = 0; i < N; ++i) {
    if (constrainedPositions.count(i) == 0) {
      angles[i] += angleAdjustments[i];
    }
  }
}

std::vector<RDGeom::Point2D> refineMacrocycleWithAngleConstraints(
    const std::vector<RDGeom::Point2D> &templateCoords,
    const std::vector<AngleConstraint> &angleConstraints, double bondLength) {
  size_t N = templateCoords.size();

  // Extract current turn angles from template coordinates
  std::vector<double> currentAngles(N);

  // Initialize with the direction of the last bond (from N-1 to 0)
  // This ensures the first turn angle is calculated correctly for a closed
  // polygon
  size_t prev = N - 1;
  double dx_initial = templateCoords[0].x - templateCoords[prev].x;
  double dy_initial = templateCoords[0].y - templateCoords[prev].y;
  double currentDirection = std::atan2(dy_initial, dx_initial);

  for (size_t i = 0; i < N; ++i) {
    size_t next = (i + 1) % N;

    // Compute bond vector from current atom to next atom
    double dx = templateCoords[next].x - templateCoords[i].x;
    double dy = templateCoords[next].y - templateCoords[i].y;
    double newDirection = std::atan2(dy, dx);

    // Turn angle is change in direction
    double turnAngle = newDirection - currentDirection;

    // Normalize to [-π, π]
    while (turnAngle > M_PI) turnAngle -= 2.0 * M_PI;
    while (turnAngle < -M_PI) turnAngle += 2.0 * M_PI;

    currentAngles[i] = turnAngle;
    currentDirection = newDirection;
  }

  // Apply angle constraints and compute angular difference
  //    Adjust constraint signs based on template
  std::vector<double> adjustedAngles = currentAngles;
  double constraintDelta = 0.0;
  std::set<size_t> constrainedPositions;
  for (const auto &constraint : angleConstraints) {
    size_t pos = constraint.position;
    double oldAngle = currentAngles[pos];

    // preserve the sign of the template's existing angle, only adjust the
    // magnitude
    double templateSign = (oldAngle >= 0) ? 1.0 : -1.0;
    double newAngle = templateSign * std::abs(constraint.targetAngle);

    adjustedAngles[pos] = newAngle;
    constraintDelta += (newAngle - oldAngle);
    constrainedPositions.insert(pos);
  }

  // Distribute angular difference across non-constrained angles to preserve
  // total angular sum
  size_t numFree = N - constrainedPositions.size();
  if (numFree > 0) {
    double adjustmentPerFree = -constraintDelta / numFree;
    for (size_t i = 0; i < N; ++i) {
      if (constrainedPositions.count(i) == 0) {
        adjustedAngles[i] += adjustmentPerFree;
      }
    }
  }

  // Generate coordinates from adjusted angles
  std::vector<RDGeom::Point2D> coords(N + 1);  // N atoms + 1 dummy
  coords[0] = RDGeom::Point2D(0, 0);

  double direction = 0.0;
  RDGeom::Point2D pos(0, 0);

  for (size_t i = 0; i < N; ++i) {
    direction += adjustedAngles[i];
    pos.x += bondLength * std::cos(direction);
    pos.y += bondLength * std::sin(direction);
    coords[i + 1] = pos;
  }

  // Close gap using Jacobian minimization
  refineWithJacobian(coords, adjustedAngles, constrainedPositions, bondLength,
                     N % 2 == 1);

  // Remove dummy atom and return
  coords.pop_back();
  return coords;
}

void maybeReflectSymmetricFusedRings(const RDKit::ROMol &mol,
                                     const RDKit::INT_VECT &macrocycleRing,
                                     const RDKit::VECT_INT_VECT &fusedRings,
                                     RDGeom::INT_POINT2D_MAP &eatoms) {
  // Convert to std::vector<int> for helper functions
  std::vector<int> macrocycleRingVec(macrocycleRing.begin(),
                                     macrocycleRing.end());

  // Check each fused ring to see if it's a 4 or 6-membered ring fused at axial
  // positions
  for (const auto &ring : fusedRings) {
    if (ring.size() != 4 && ring.size() != 6) {
      continue;
    }

    // find shared positions in macrocycle
    std::vector<int> ringVec(ring.begin(), ring.end());
    std::vector<size_t> sharedPositions =
        findSharedPositions(macrocycleRingVec, ringVec);

    // Verify and reorder shared positions
    if (!verifyAndReorderSharedPositions(sharedPositions,
                                         macrocycleRing.size())) {
      continue;  // Not contiguous, skip
    }

    size_t numShared = sharedPositions.size();

    // Check if we have axial fusion
    // For 6-membered: 4 shared, 2 non-shared
    // For 4-membered: 3 shared, 1 non-shared
    if ((ring.size() == 6 && numShared == 4) ||
        (ring.size() == 4 && numShared == 3)) {
      // Axial atoms are at first and last shared positions
      int axial1 = macrocycleRing[sharedPositions.front()];
      int axial2 = macrocycleRing[sharedPositions.back()];

      // Middle shared atoms are all except first and last
      std::vector<int> middleShared;
      for (size_t i = 1; i < numShared - 1; ++i) {
        middleShared.push_back(macrocycleRing[sharedPositions[i]]);
      }

      std::unordered_set<int> sharedSet;
      for (auto pos : sharedPositions) {
        sharedSet.insert(macrocycleRing[pos]);
      }
      std::vector<int> nonSharedAtoms;
      for (int atom : ring) {
        if (sharedSet.find(atom) == sharedSet.end()) {
          nonSharedAtoms.push_back(atom);
        }
      }

      // Count substituents on middle shared atoms (degree - 2)
      int middleSubstituents = 0;
      for (int atom : middleShared) {
        const auto *ratom = mol.getAtomWithIdx(atom);
        middleSubstituents += ratom->getDegree() - 2;
      }

      // Count substituents on non-shared atoms
      int nonSharedSubstituents = 0;
      for (int atom : nonSharedAtoms) {
        const auto *ratom = mol.getAtomWithIdx(atom);
        nonSharedSubstituents += ratom->getDegree() - 2;
      }

      // If middle shared atoms are more substituted, reflect them
      if (middleSubstituents > nonSharedSubstituents) {
        // Check that axial atoms have coordinates
        if (eatoms.find(axial1) != eatoms.end() &&
            eatoms.find(axial2) != eatoms.end()) {
          const auto &pin1 = eatoms[axial1];
          const auto &pin2 = eatoms[axial2];

          // Reflect middle shared atoms. Non shared atoms are not placed yet so
          // there's no need to flip them, they will fall in the right place
          // when the small ring is embedded later
          for (int atom : middleShared) {
            if (eatoms.find(atom) != eatoms.end()) {
              eatoms[atom] = reflectPoint(eatoms[atom], pin1, pin2);
            }
          }
        }
      }
    }
  }
}

// Helper: Add turn constraints for fused rings
static void addFusedRingTurnConstraints(
    MacrocycleGenerator &generator, const std::vector<int> &macrocycleRingVec,
    const RDKit::VECT_INT_VECT &allRings) {
  // Add turn constraints for fused rings
  // Pattern: first turn = R, middle turns = L, last turn = R
  // Examples: 2-atom fusion → RR, 3-atom → RLR, 4-atom → RLLR
  for (const auto &ring : allRings) {
    // Skip macrocycles
    if (ring.size() > MACROCYCLE_SIZE_THRESHOLD) {
      continue;
    }

    // Find shared atoms between this ring and the macrocycle
    std::vector<size_t> sharedPositions = findSharedPositions(
        macrocycleRingVec, std::vector<int>(ring.begin(), ring.end()));

    // Need at least 1 shared atom (spiro or fusion)
    // Skip if more than 4 shared atoms (unusual/complex fusion)
    if (sharedPositions.size() < 1 || sharedPositions.size() > 4) {
      continue;
    }

    // Verify shared atoms are contiguous and reorder if needed
    if (!verifyAndReorderSharedPositions(sharedPositions,
                                         macrocycleRingVec.size())) {
      continue;
    }

    // Create turn pattern based on number of shared atoms:
    // Pattern length = N (number of shared atoms)
    // First and last turn: R (convex, prevalent)
    // Middle turns: L (concave)
    // 1 shared atom: R
    // 2 shared atoms: RR
    // 3 shared atoms: RLR
    // 4 shared atoms: RLLR
    std::vector<int> pattern;
    size_t numShared = sharedPositions.size();

    if (numShared == 1) {
      pattern.push_back(1);  // R
    } else if (numShared == 2) {
      // RR pattern for 2-atom fusion
      pattern.push_back(1);  // R
      pattern.push_back(1);  // R
    } else {
      // First turn: R (endpoint)
      pattern.push_back(1);
      // Middle turns: L
      for (size_t i = 1; i < numShared - 1; ++i) {
        pattern.push_back(-1);  // L
      }
      // Last turn: R (endpoint)
      pattern.push_back(1);
    }

    // Add FIXED turn constraint
    TurnConstraint constraint;
    constraint.position = sharedPositions[0];
    constraint.type = ConstraintType::FIXED;
    constraint.pattern = pattern;
    constraint.reason = "fused_ring_" + std::to_string(ring.size()) + "_at_" +
                        std::to_string(sharedPositions[0]);

    generator.addConstraint(constraint);
  }
}

DoubleBondStereoAtoms getDoubleBondStereoAtoms(const RDKit::Bond *bond,
                                               const RDKit::ROMol *mol) {
  DoubleBondStereoAtoms result;

  // Get the double bond atoms
  result.atom1 = bond->getBeginAtomIdx();
  result.atom2 = bond->getEndAtomIdx();

  // Get the stereo-defining atoms (based on CIP priority)
  const auto &neighbors = bond->getStereoAtoms();
  if (neighbors.size() != 2) {
    return result;  // invalid
  }

  // Figure out which stereo atom is connected to which double bond atom
  bool neighbor0ConnectedToAtom1 =
      mol->getBondBetweenAtoms(result.atom1, neighbors[0]) != nullptr;
  bool neighbor0ConnectedToAtom2 =
      mol->getBondBetweenAtoms(result.atom2, neighbors[0]) != nullptr;

  if (neighbor0ConnectedToAtom1 && !neighbor0ConnectedToAtom2) {
    result.atom1Neighbor1 = neighbors[0];
    result.atom2Neighbor1 = neighbors[1];
  } else if (neighbor0ConnectedToAtom2 && !neighbor0ConnectedToAtom1) {
    result.atom1Neighbor1 = neighbors[1];
    result.atom2Neighbor1 = neighbors[0];
  } else {
    return result;  // invalid - ambiguous connectivity
  }

  // Find other neighbors for each atom (if degree > 2)
  if (mol->getAtomWithIdx(result.atom1)->getDegree() > 2) {
    for (auto neighbor :
         mol->atomNeighbors(mol->getAtomWithIdx(result.atom1))) {
      int nidx = neighbor->getIdx();
      if (nidx != result.atom1Neighbor1 && nidx != result.atom2) {
        result.atom1Neighbor2 = nidx;
        break;
      }
    }
  }

  if (mol->getAtomWithIdx(result.atom2)->getDegree() > 2) {
    for (auto neighbor :
         mol->atomNeighbors(mol->getAtomWithIdx(result.atom2))) {
      int nidx = neighbor->getIdx();
      if (nidx != result.atom2Neighbor1 && nidx != result.atom1) {
        result.atom2Neighbor2 = nidx;
        break;
      }
    }
  }

  result.valid = true;
  return result;
}

// Helper: Add stereochemistry constraints for double bonds
static void addStereochemistryConstraints(
    MacrocycleGenerator &generator, const RDKit::ROMol *mol,
    const RDKit::INT_VECT &macrocycleRing) {
  // Add constraints for double bonds with E/Z stereochemistry
  for (size_t i = 0; i < macrocycleRing.size(); ++i) {
    size_t nextIdx = (i + 1) % macrocycleRing.size();
    int atom1 = macrocycleRing[i];
    int atom2 = macrocycleRing[nextIdx];

    const RDKit::Bond *bond = mol->getBondBetweenAtoms(atom1, atom2);
    if (!bond || bond->getBondType() != RDKit::Bond::DOUBLE) {
      continue;
    }

    auto stereoType = bond->getStereo();
    if (stereoType == RDKit::Bond::STEREOANY ||
        stereoType == RDKit::Bond::STEREONONE) {
      continue;
    }

    // Extract all atoms around the double bond
    auto stereoAtoms = getDoubleBondStereoAtoms(bond, mol);
    if (!stereoAtoms.valid) {
      continue;
    }

    // Determine the macrocycle ring neighbors
    size_t prevIdx = (i - 1 + macrocycleRing.size()) % macrocycleRing.size();
    size_t nextNextIdx = (i + 2) % macrocycleRing.size();
    int macrocycle_atom1_neighbor = macrocycleRing[prevIdx];
    int macrocycle_atom2_neighbor = macrocycleRing[nextNextIdx];

    // Check if stereo-defining atoms match the macrocycle neighbors
    // If they don't, we need to swap the interpretation
    bool swapStereo = false;

    // Check if atom1_neighbor1 is the macrocycle neighbor
    // If not, use atom1_neighbor2 (if it exists)
    if (stereoAtoms.atom1Neighbor1 != macrocycle_atom1_neighbor) {
      if (stereoAtoms.atom1Neighbor2 == macrocycle_atom1_neighbor) {
        swapStereo = !swapStereo;
      } else {
        continue;
      }
    }

    // Check if atom2_neighbor1 is the macrocycle neighbor
    // If not, use atom2_neighbor2 (if it exists)
    if (stereoAtoms.atom2Neighbor1 != macrocycle_atom2_neighbor) {
      if (stereoAtoms.atom2Neighbor2 == macrocycle_atom2_neighbor) {
        swapStereo = !swapStereo;
      } else {
        continue;
      }
    }

    // Determine expected stereochemistry
    bool expected_cis = (stereoType == RDKit::Bond::STEREOZ ||
                         stereoType == RDKit::Bond::STEREOCIS);

    // Apply swap if needed
    if (swapStereo) {
      expected_cis = !expected_cis;
    }

    // Add constraint based on cis/trans relationship
    // For a double bond at positions i → i+1 in the ring:
    // The turns at positions i and i+1 determine the stereochemistry
    // Cis: turns should be SAME
    // Trans: turns should be OPPOSITE
    // Note: OPPOSITE/SAME constraint at position X constrains positions X and
    // X+1

    TurnConstraint constraint;
    constraint.position = i;  // This will constrain turns at i and i+1

    if (expected_cis) {
      constraint.type = ConstraintType::SAME;
      constraint.reason = "cis_double_bond_at_" + std::to_string(i);
    } else {
      constraint.type = ConstraintType::OPPOSITE;
      constraint.reason = "trans_double_bond_at_" + std::to_string(i);
    }

    generator.addConstraint(constraint);
  }
}

// Shared implementation for identifying triple bond constraints
std::vector<AngleConstraint> identifyAngleConstraintsForTripleBonds(
    const RDKit::ROMol *mol, const std::vector<int> &macrocycleRing) {
  std::vector<AngleConstraint> angleConstraints;

  if (!mol) {
    return angleConstraints;
  }

  // Iterate through macrocycle bonds to find triple bonds
  for (size_t i = 0; i < macrocycleRing.size(); ++i) {
    size_t nextIdx = (i + 1) % macrocycleRing.size();
    int atom1 = macrocycleRing[i];
    int atom2 = macrocycleRing[nextIdx];

    const auto *bond = mol->getBondBetweenAtoms(atom1, atom2);
    if (!bond || bond->getBondType() != RDKit::Bond::TRIPLE) {
      continue;
    }

    // Add 180° angle constraint for both atoms of the triple bond
    AngleConstraint constraint1;
    constraint1.position = i;
    constraint1.targetAngle = 0;  // flat
    angleConstraints.push_back(constraint1);

    AngleConstraint constraint2;
    constraint2.position = nextIdx;
    constraint2.targetAngle = 0;
    angleConstraints.push_back(constraint2);
  }

  return angleConstraints;
}

// Helper: Add angle constraints for triple bonds to generator
static void addTripleBondAngleConstraints(
    MacrocycleGenerator &generator, const RDKit::ROMol *mol,
    const RDKit::INT_VECT &macrocycleRing) {
  std::vector<int> ringVec(macrocycleRing.begin(), macrocycleRing.end());
  auto constraints = identifyAngleConstraintsForTripleBonds(mol, ringVec);
  for (const auto &constraint : constraints) {
    generator.addAngleConstraint(constraint);
  }
}

std::vector<RDGeom::Point2D> generateMacrocycleCoordinates(
    const RDKit::ROMol *mol, const RDKit::INT_VECT &macrocycleRing,
    const RDKit::VECT_INT_VECT &allRings,
    const std::map<size_t, int> &substituentSizesByPosition,
    double bondLength) {
  // Only attempt for macrocycles (size > MACROCYCLE_SIZE_THRESHOLD)
  if (macrocycleRing.size() <= MACROCYCLE_SIZE_THRESHOLD) {
    return {};
  }

  // Create MacrocycleGenerator
  MacrocycleGenerator generator(macrocycleRing.size(), bondLength);

  // Convert allRings to std::vector<std::vector<int>> for helper function
  std::vector<std::vector<int>> allRingsVec;
  for (const auto &ring : allRings) {
    std::vector<int> ringVec(ring.begin(), ring.end());
    allRingsVec.push_back(ringVec);
  }
  std::vector<int> macrocycleRingVec(macrocycleRing.begin(),
                                     macrocycleRing.end());

  // Identify all angle constraints for fused small rings
  auto angleConstraints =
      identifyAngleConstraintsForFusedRings(macrocycleRingVec, allRingsVec);
  for (const auto &angleConstraint : angleConstraints) {
    generator.addAngleConstraint(angleConstraint);
  }

  // Use pre-computed substituent sizes by position
  generator.setSubstituentInfo(substituentSizesByPosition,
                               -1);  // -1 = L turns point inward

  // Add constraints for fused rings and stereochemistry
  addFusedRingTurnConstraints(generator, macrocycleRingVec, allRings);
  addStereochemistryConstraints(generator, mol, macrocycleRing);
  addTripleBondAngleConstraints(generator, mol, macrocycleRing);

  // Get number of free positions before solving (for threshold scaling)
  size_t numFreePositions = generator.getNumFreePositions();

  // Solve for optimal turn sequence
  if (!generator.solve()) {
    return {};  // No solution found
  }

  // Check closure quality
  double closureError = generator.getClosureError();
  double MAX_CLOSURE_ERROR = 6.0;  // base threshold in Angstroms

  // Scale threshold based on number of free turn positions (not total ring
  // size) Fast heuristic is used when free positions > 15, creating approximate
  // patterns
  if (numFreePositions > 15) {
    // Allow proportionally larger initial closure errors for fast heuristic
    // Jacobian refinement will close the gap
    MAX_CLOSURE_ERROR = numFreePositions * 1.5;  // ~1.5 Å per free position
  }

  if (closureError > MAX_CLOSURE_ERROR) {
    return {};  // Closure error too large
  }

  // Generate and return 2D coordinates
  return generator.generateCoordinates();
}

}  // namespace RDDepict
