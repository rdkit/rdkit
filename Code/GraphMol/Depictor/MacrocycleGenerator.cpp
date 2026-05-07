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
#include <map>

namespace RDDepict {

MacrocycleGenerator::MacrocycleGenerator(size_t ringSize, double bondLength,
                                         bool useJacobianRefinement)
    : d_ringSize(ringSize),
      d_bondLength(bondLength),
      d_innerTurnSign(0),
      d_closureError(std::numeric_limits<double>::max()),
      d_solved(false),
      d_useJacobianRefinement(useJacobianRefinement) {
  d_turns.resize(ringSize, 0);
}

void MacrocycleGenerator::addConstraint(const TurnConstraint &constraint) {
  d_constraints.push_back(constraint);
}

void MacrocycleGenerator::addAngleConstraint(const AngleConstraint &constraint) {
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

void MacrocycleGenerator::setSubstituentInfo(const std::map<size_t, int> &substituentSizes, int innerTurnSign) {
  d_substituentSizes = substituentSizes;
  d_innerTurnSign = innerTurnSign;

  // Debug: Show what substituent info was set
  std::cerr << "SUBST_INFO_SET: innerTurnSign=" << innerTurnSign << ", substituents={";
  for (const auto &[pos, size] : substituentSizes) {
    std::cerr << " [" << pos << "]=" << size;
  }
  std::cerr << " }" << std::endl;
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

void MacrocycleGenerator::debugPrintAngleConstraints(
    const std::vector<RDGeom::Point2D> &coords) const {
  std::cerr << "\n[DEBUG] Final Angles (macrocycle size=" << d_ringSize
            << ", coords.size()=" << coords.size() << "):" << std::endl;

  // Print ALL angles
  size_t N = d_ringSize;
  for (size_t pos = 0; pos < N; ++pos) {
    const auto &p_prev = coords[(pos + N - 1) % N];
    const auto &p_curr = coords[pos];
    const auto &p_next = coords[(pos + 1) % N];

    double v_before_x = p_curr.x - p_prev.x;
    double v_before_y = p_curr.y - p_prev.y;
    double v_after_x = p_next.x - p_curr.x;
    double v_after_y = p_next.y - p_curr.y;

    double dir_before = std::atan2(v_before_y, v_before_x);
    double dir_after = std::atan2(v_after_y, v_after_x);
    double angle = dir_after - dir_before;

    while (angle > M_PI) angle -= 2 * M_PI;
    while (angle < -M_PI) angle += 2 * M_PI;

    double angleDeg = angle * 180.0 / M_PI;

    // Check if this position has a constraint
    bool isConstrained = hasAngleConstraint(pos);
    double targetDeg = 0.0;
    if (isConstrained) {
      targetDeg = getConstraintAngle(pos) * 180.0 / M_PI;
    }

    std::cerr << "  Pos " << pos << ": " << angleDeg << "°";
    if (isConstrained) {
      double dev = std::abs(angle - getConstraintAngle(pos)) * 180.0 / M_PI;
      std::cerr << " [CONSTRAINED target=" << targetDeg << "° dev=" << dev << "°]";
    }
    std::cerr << std::endl;
  }

  if (d_angleConstraints.empty()) {
    return;  // No constraints, no output
  }

  std::cerr << "\n[DEBUG] Angle Constraints: " << d_angleConstraints.size()
            << " constraints found" << std::endl;

  for (const auto &constraint : d_angleConstraints) {
    size_t pos = constraint.position;
    double targetAngle = constraint.targetAngle;
    double targetAngleDeg = targetAngle * 180.0 / M_PI;

    // Compute actual angle at this position
    // The turn at position i is the rotation from direction (i-1→i) to (i→i+1)
    size_t N = d_ringSize;
    if (coords.size() < N) {
      continue;
    }

    // Turn at position 'pos' is the rotation AFTER placing atom pos
    // It affects the direction from atom pos to atom pos+1
    // To measure it, we need: direction_before_turn → direction_after_turn
    // direction_before = from (pos-1) to pos
    // direction_after = from pos to (pos+1)
    const auto &p_prev = coords[(pos + N - 1) % N];  // pos-1 (wrap around)
    const auto &p_curr = coords[pos];
    const auto &p_next = coords[(pos + 1) % N];

    // Compute direction vectors
    double v_before_x = p_curr.x - p_prev.x;
    double v_before_y = p_curr.y - p_prev.y;
    double v_after_x = p_next.x - p_curr.x;
    double v_after_y = p_next.y - p_curr.y;

    // Compute directions as angles
    double dir_before = std::atan2(v_before_y, v_before_x);
    double dir_after = std::atan2(v_after_y, v_after_x);

    // Turn angle is the change in direction
    double actualAngle = dir_after - dir_before;

    // Normalize to [-π, π]
    while (actualAngle > M_PI) actualAngle -= 2 * M_PI;
    while (actualAngle < -M_PI) actualAngle += 2 * M_PI;

    double actualAngleDeg = actualAngle * 180.0 / M_PI;
    double deviation = std::abs(actualAngle - targetAngle) * 180.0 / M_PI;

    std::cerr << "  Position " << pos << ": target=" << targetAngleDeg
              << "°, actual=" << actualAngleDeg << "°, deviation=" << deviation
              << "°" << std::endl;
  }

  std::cerr << std::endl;
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
  // Normal: R-L = 6 (even) or 5 (odd)
  // First relaxation: ±2 from normal
  // Second relaxation (odd only): ±4 from normal
  std::vector<int> targetDiffs;

  if (isOddRing) {
    targetDiffs = {5, 3, 7, 1, 9};  // normal, then ±2, then ±4
  } else {
    targetDiffs = {6, 4, 8, 2, 10};  // normal, then ±2, then ±4
  }

  for (int targetDiff : targetDiffs) {
    std::cerr << "SOLVE_ATTEMPT: Trying targetDiff(R-L)=" << targetDiff << std::endl;
    if (trySolveWithTargetDiff(targetDiff)) {
      std::cerr << "SOLVE_SUCCESS: Found solution with targetDiff=" << targetDiff << std::endl;
      return true;
    }
  }

  std::cerr << "SOLVE_FAILED: All targetDiff values failed" << std::endl;
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

  bool isOddRing = (d_ringSize % 2) == 1;

  std::cerr << "CONSTRAINT_INFO: ringSize=" << d_ringSize << ", isOddRing=" << isOddRing
            << ", targetDiff(R-L)=" << targetDiff << std::endl;
  std::cerr << "CONSTRAINT_INFO: decidedRight=" << decidedRight
            << ", decidedLeft=" << decidedLeft << ", freeCount=" << freeCount << std::endl;

  // Show angle constraints
  if (!d_angleConstraints.empty()) {
    std::cerr << "ANGLE_CONSTRAINTS: " << d_angleConstraints.size() << " constraints set:" << std::endl;
    for (const auto &c : d_angleConstraints) {
      int turn_sign = (c.targetAngle > 0) ? 1 : -1;
      std::cerr << "  pos=" << c.position << ", targetAngle=" << (c.targetAngle * 180.0 / M_PI)
                << "°, turn=" << (turn_sign > 0 ? "R" : "L") << std::endl;
    }
  }

  // Show which positions are already decided (constrained to R or L)
  std::cerr << "DECIDED_TURNS: ";
  for (size_t i = 0; i < d_turns.size(); ++i) {
    if (d_turns[i] != 0) {
      std::cerr << " [" << i << "]=" << (d_turns[i] == 1 ? "R" : "L");
    }
  }
  std::cerr << std::endl;

  // We need: (decidedRight + freeRight) - (decidedLeft + freeLeft) = targetDiff
  // Where: freeRight + freeLeft = freeCount
  // Solving: freeRight = (targetDiff - decidedRight + decidedLeft + freeCount)
  // / 2

  int totalRight = (d_ringSize + targetDiff) / 2;
  int totalLeft = (d_ringSize - targetDiff) / 2;

  std::cerr << "SOLVE_CALC: totalRight=" << totalRight << ", totalLeft=" << totalLeft << std::endl;

  // Check parity
  if ((d_ringSize + targetDiff) % 2 != 0) {
    std::cerr << "SOLVE_FAILED: Parity mismatch" << std::endl;
    return false;  // Parity mismatch - no integer solution
  }

  // Check if solution is possible
  if (totalRight < 0 || totalLeft < 0) {
    std::cerr << "SOLVE_FAILED: Negative total (totalRight=" << totalRight << ", totalLeft=" << totalLeft << ")" << std::endl;
    return false;  // Ring too small for angular closure
  }

  int freeRight = totalRight - decidedRight;
  int freeLeft = totalLeft - decidedLeft;

  std::cerr << "SOLVE_CALC: freeRight=" << freeRight << ", freeLeft=" << freeLeft << std::endl;

  if (freeRight < 0 || freeLeft < 0) {
    std::cerr << "SOLVE_FAILED: Negative free (freeRight=" << freeRight << ", freeLeft=" << freeLeft << ")" << std::endl;
    return false;  // Constraints over-constrain the solution
  }

  if (freeRight + freeLeft != freeCount) {
    std::cerr << "SOLVE_FAILED: Sanity check (freeRight+freeLeft=" << (freeRight+freeLeft) << " != freeCount=" << freeCount << ")" << std::endl;
    return false;  // Sanity check failed
  }

  // Find optimal turn sequence for free variables
  std::cerr << "SOLVE: Calling findOptimalTurnSequence(freeRight=" << freeRight << ", freeLeft=" << freeLeft << ")" << std::endl;
  d_solved = findOptimalTurnSequence(freeRight, freeLeft);
  if (!d_solved) {
    std::cerr << "SOLVE_FAILED: findOptimalTurnSequence returned false" << std::endl;
  } else {
    std::cerr << "SOLVE_SUCCESS: Found valid turn sequence" << std::endl;
  }
  return d_solved;
}

size_t MacrocycleGenerator::getNumFreePositions() const {
  size_t count = 0;
  for (size_t i = 0; i < d_ringSize; ++i) {
    if (d_turns[i] == 0) {
      count++;
    }
  }
  return count;
}

bool MacrocycleGenerator::findOptimalTurnSequence(int numRight, int numLeft) {
  // Find positions of free variables (currently 0)
  std::vector<size_t> freePositions;
  for (size_t i = 0; i < d_ringSize; ++i) {
    if (d_turns[i] == 0) {
      freePositions.push_back(i);
    }
  }

  // Fast heuristic for large rings: use alternating pattern instead of enumeration
  // When there are many free positions (>15), enumeration becomes computationally infeasible
  const size_t ENUMERATION_THRESHOLD = 15;

  if (freePositions.size() > ENUMERATION_THRESHOLD) {
    size_t numFree = freePositions.size();
    std::cerr << "FAST_HEURISTIC: Using substituent-based pattern for " << numFree << " free positions" << std::endl;

    // Strategy: Use substituent info to reduce degrees of freedom
    // Positions with big substituents (size > 1) should be outer turns
    // This is both chemically sensible and reduces the search space
    std::vector<int> candidate = d_turns;  // Start with constrained values

    int outerTurn = -d_innerTurnSign;  // Outer turn is opposite of inner
    int assignedRight = 0;
    int assignedLeft = 0;

    std::cerr << "FAST_HEURISTIC: innerTurnSign=" << d_innerTurnSign
              << ", outerTurn=" << (outerTurn == 1 ? "R" : "L") << std::endl;

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
        std::cerr << "  pos[" << pos << "] subst_size=" << it->second
                  << " -> OUTER (" << (outerTurn == 1 ? "R" : "L") << ")" << std::endl;
      }
    }

    // Step 2: Check if we have capacity for these assignments
    int remainingRight = numRight - assignedRight;
    int remainingLeft = numLeft - assignedLeft;

    std::cerr << "FAST_HEURISTIC: After " << substituentAssignments << " substituent assignments: "
              << "assignedR=" << assignedRight << " assignedL=" << assignedLeft
              << ", remainingR=" << remainingRight << " remainingL=" << remainingLeft << std::endl;

    if (remainingRight < 0 || remainingLeft < 0) {
      std::cerr << "FAST_HEURISTIC: Substituent assignment exceeds required counts" << std::endl;
      return false;
    }

    // Step 3: If still too many free positions, use alternating for most, keep few truly free
    std::vector<size_t> stillFree;
    for (size_t pos : freePositions) {
      if (candidate[pos] == 0) {
        stillFree.push_back(pos);
      }
    }

    const size_t MAX_ENUMERATION_SIZE = 10;  // Only enumerate over ~10 positions

    if (stillFree.size() > MAX_ENUMERATION_SIZE) {
      std::cerr << "FAST_HEURISTIC: Still have " << stillFree.size()
                << " free positions after substituents, fixing most with alternating" << std::endl;

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

      std::cerr << "FAST_HEURISTIC: Fixed " << fixed << " positions with alternating, "
                << (stillFree.size() - fixed) << " still free for enumeration" << std::endl;

      // Update d_turns with the fixed positions, then fall through to normal enumeration
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
      std::cerr << "FAST_HEURISTIC: Reduced free positions from " << numFree
                << " to " << freePositions.size()
                << ", need " << numRight << " R and " << numLeft << " L"
                << ", continuing with enumeration" << std::endl;

      // Don't return - fall through to enumeration code below

    } else {
      // Few enough positions - just assign them all with alternating
      std::cerr << "FAST_HEURISTIC: Only " << stillFree.size()
                << " positions remain, assigning all with alternating" << std::endl;

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

          std::cerr << "FAST_HEURISTIC: Success (closureError=" << d_closureError
                    << ", minDist=" << minDistance << ")" << std::endl;
          return true;
        } else {
          std::cerr << "FAST_HEURISTIC: Self-crossing detected (minDist=0)" << std::endl;
        }
      } else {
        std::cerr << "FAST_HEURISTIC: Constraint validation failed" << std::endl;
      }

      // All positions assigned but failed - can't enumerate further
      std::cerr << "FAST_HEURISTIC: Failed after assigning all " << numFree
                << " free positions" << std::endl;
      return false;
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
    std::vector<int> candidate(d_ringSize, -1);  // Initialize all to L
    std::vector<size_t> rightPositions(numRight);

    // Initialize first combination: [0, 1, 2, ..., numRight-1]
    for (int i = 0; i < numRight; ++i) {
      rightPositions[i] = i;
      candidate[i] = 1;  // Set to R
    }

    double bestScore = std::numeric_limits<double>::max();
    std::vector<int> bestTurns;

    // Weight for minimum distance in scoring (balances closure error vs roundness)
    const double MIN_DISTANCE_WEIGHT = 0.5;

    // Evaluate first combination
    if (validateConstraints(candidate)) {
      int minDistance = calculateMinDistance(candidate);
      // Skip if self-crossing detected (minDistance == 0)
      if (minDistance > 0) {
        double closureError = calculateClosureError(candidate);
        double substituentPenalty = calculateSubstituentPenalty(candidate);
        // Combined score: closure error - roundness bonus + substituent penalty
        double score = closureError - MIN_DISTANCE_WEIGHT * minDistance + substituentPenalty;
        if (score < bestScore) {
          bestScore = score;
          bestTurns = candidate;
        }
      }
    }

    // Enumerate remaining combinations
    while (nextCombination(rightPositions, d_ringSize)) {
      // Reset candidate to all L
      std::fill(candidate.begin(), candidate.end(), -1);

      // Set R positions
      for (size_t pos : rightPositions) {
        candidate[pos] = 1;
      }

      // Validate constraints and calculate min distance (includes clash detection)
      if (validateConstraints(candidate)) {
        int minDistance = calculateMinDistance(candidate);
        // Skip if self-crossing detected
        if (minDistance > 0) {
          double closureError = calculateClosureError(candidate);
          double substituentPenalty = calculateSubstituentPenalty(candidate);
          double score = closureError - MIN_DISTANCE_WEIGHT * minDistance + substituentPenalty;
          if (score < bestScore) {
            bestScore = score;
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
    d_closureError = calculateClosureError(bestTurns);

    // Debug: Show scoring components for winning solution
    {
      int minDistance = calculateMinDistance(bestTurns);
      double substituentPenalty = calculateSubstituentPenalty(bestTurns);
      double minDistWeighted = MIN_DISTANCE_WEIGHT * minDistance;
      double totalScore = d_closureError - minDistWeighted + substituentPenalty;
      std::cerr << "MACROCYCLE_SCORE (fully-free): closureError=" << d_closureError
                << ", minDist=" << minDistance << " (-0.5*" << minDistance << "=" << -minDistWeighted << ")"
                << ", substPenalty=" << substituentPenalty
                << ", totalScore=" << totalScore << std::endl;

      // Show turn sequence for debugging
      std::cerr << "WINNING_TURNS: [";
      int countR = 0, countL = 0;
      for (size_t i = 0; i < bestTurns.size(); ++i) {
        if (i > 0) std::cerr << ", ";
        std::cerr << bestTurns[i];
        if (bestTurns[i] == 1) countR++;
        else if (bestTurns[i] == -1) countL++;
      }
      std::cerr << "]" << std::endl;
      std::cerr << "WINNING_TURNS: R=" << countR << ", L=" << countL << ", R-L=" << (countR - countL) << std::endl;
    }

    return true;
  }

  // Mixed case: some positions fixed, some free
  // Enumerate ways to assign R/L to free positions
  std::vector<size_t> rightPositions(numRight);
  size_t numFree = freePositions.size();

  // Initialize first combination: choose first numRight free positions for R
  for (int i = 0; i < numRight; ++i) {
    rightPositions[i] = i;
  }

  double bestScore = std::numeric_limits<double>::max();
  std::vector<int> bestTurns;

  // Weight for minimum distance in scoring
  const double MIN_DISTANCE_WEIGHT = 0.5;

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

  // Evaluate first combination
  int combinationCount = 0;
  int rejectedConstraints = 0, rejectedMinDist = 0;
  std::vector<int> candidate = createCandidate();
  combinationCount++;
  if (validateConstraints(candidate)) {
    int minDistance = calculateMinDistance(candidate);
    // Skip if self-crossing detected
    if (minDistance > 0) {
      double closureError = calculateClosureError(candidate);
      double substituentPenalty = calculateSubstituentPenalty(candidate);
      double score = closureError - MIN_DISTANCE_WEIGHT * minDistance + substituentPenalty;
      if (score < bestScore) {
        bestScore = score;
        bestTurns = candidate;
      }
    } else {
      rejectedMinDist++;
    }
  } else {
    rejectedConstraints++;
  }

  // Enumerate remaining combinations
  while (nextCombination(rightPositions, numFree)) {
    candidate = createCandidate();
    combinationCount++;
    if (validateConstraints(candidate)) {
      int minDistance = calculateMinDistance(candidate);
      // Skip if self-crossing detected
      if (minDistance > 0) {
        double closureError = calculateClosureError(candidate);
        double substituentPenalty = calculateSubstituentPenalty(candidate);
        double score = closureError - MIN_DISTANCE_WEIGHT * minDistance + substituentPenalty;
        if (score < bestScore) {
          bestScore = score;
          bestTurns = candidate;
        }
      } else {
        rejectedMinDist++;
      }
    } else {
      rejectedConstraints++;
    }
  }

  std::cerr << "ENUMERATION: Checked " << combinationCount << " combinations: "
            << "rejectedConstraints=" << rejectedConstraints
            << ", rejectedMinDist=" << rejectedMinDist
            << ", accepted=" << (combinationCount - rejectedConstraints - rejectedMinDist) << std::endl;

  // Store best solution
  if (bestTurns.empty()) {
    std::cerr << "ENUMERATION_FAILED: No valid combinations found" << std::endl;
    return false;  // No valid solution found that satisfies constraints
  }

  d_turns = bestTurns;
  d_closureError = calculateClosureError(bestTurns);

  // Debug: Test alternative sequence suggested by user (flip first turn)
  if (bestTurns.size() > 0) {
    std::vector<int> altTurns = bestTurns;
    altTurns[0] = -altTurns[0];  // Flip first turn

    // Check if still valid with constraints
    if (validateConstraints(altTurns)) {
      int altR = 0, altL = 0;
      for (int t : altTurns) {
        if (t == 1) altR++;
        else if (t == -1) altL++;
      }
      double altClosure = calculateClosureError(altTurns);
      int altMinDist = calculateMinDistance(altTurns);
      double altSubstPenalty = calculateSubstituentPenalty(altTurns);
      double altScore = altClosure - MIN_DISTANCE_WEIGHT * altMinDist + altSubstPenalty;

      std::cerr << "ALTERNATIVE_SEQUENCE (flip pos 0): R=" << altR << ", L=" << altL
                << ", R-L=" << (altR - altL)
                << ", closureError=" << altClosure
                << ", minDist=" << altMinDist
                << ", substPenalty=" << altSubstPenalty
                << ", score=" << altScore << std::endl;
    } else {
      std::cerr << "ALTERNATIVE_SEQUENCE (flip pos 0): INVALID (violates constraints)" << std::endl;
    }
  }

  // Debug: Show scoring components for winning solution
  {
    int minDistance = calculateMinDistance(bestTurns);
    double substituentPenalty = calculateSubstituentPenalty(bestTurns);
    double minDistWeighted = MIN_DISTANCE_WEIGHT * minDistance;
    double totalScore = d_closureError - minDistWeighted + substituentPenalty;
    std::cerr << "MACROCYCLE_SCORE (partial-constrained): closureError=" << d_closureError
              << ", minDist=" << minDistance << " (-0.5*" << minDistance << "=" << -minDistWeighted << ")"
              << ", substPenalty=" << substituentPenalty
              << ", totalScore=" << totalScore << std::endl;

    // Show turn sequence for debugging
    std::cerr << "WINNING_TURNS: [";
    int countR = 0, countL = 0;
    for (size_t i = 0; i < bestTurns.size(); ++i) {
      if (i > 0) std::cerr << ", ";
      std::cerr << bestTurns[i];
      if (bestTurns[i] == 1) countR++;
      else if (bestTurns[i] == -1) countL++;
    }
    std::cerr << "]" << std::endl;
    std::cerr << "WINNING_TURNS: R=" << countR << ", L=" << countL << ", R-L=" << (countR - countL) << std::endl;
  }

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
    std::cerr << "SUBST_PENALTY: NO_PENALTY (empty=" << d_substituentSizes.empty()
              << ", innerSign=" << d_innerTurnSign << ")" << std::endl;
    return 0.0;  // No substituents or inner turn not set
  }

  // Calculate penalty for substituents on inner turns
  // Penalty is weighted by substituent size
  double penalty = 0.0;
  const double PENALTY_PER_ATOM = 1.0;  // Penalty per substituent atom on inner turn

  std::cerr << "SUBST_PENALTY_CHECK: innerTurnSign=" << d_innerTurnSign << std::endl;
  for (const auto &[pos, size] : d_substituentSizes) {
    if (pos < turns.size()) {
      int turnAtPos = turns[pos];
      bool isInner = (turnAtPos == d_innerTurnSign);
      std::cerr << "  pos=" << pos << ", size=" << size << ", turn=" << turnAtPos
                << (isInner ? " [INNER - PENALTY]" : " [OUTER - OK]") << std::endl;
      if (isInner) {
        penalty += size * PENALTY_PER_ATOM;
      }
    }
  }
  std::cerr << "SUBST_PENALTY: total=" << penalty << std::endl;

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
        turnAngle = getConstraintAngle(i);  // Use exact ideal angle, NO extraPerAtom
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

  std::cerr << "JACOBIAN_START: initialGap=" << initialGap << std::endl;

  if (initialGap < 1e-6) {
    std::cerr << "JACOBIAN_SKIP: gap already closed" << std::endl;
    return;
  }

  // Extra angle per atom for odd rings
  // Count constrained angles and distribute extraPerAtom only across free angles
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

    std::cerr << "JACOBIAN_ITER[" << iter << "]: gapMag=" << gapMag << std::endl;

    if (gapMag < 0.01) {
      std::cerr << "JACOBIAN_CONVERGED: gap < 0.01" << std::endl;
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
  std::cerr << "JACOBIAN_END: finalGap=" << finalGap
            << ", improvement=" << (initialGap - finalGap) << std::endl;
}

// ============================================================================
// Helper Functions for Angle Constraint Detection and Refinement
// ============================================================================

double computeIdealAngle(int ringSize, int turnSign) {
  double turnAngle = 2.0 * M_PI / ringSize;  // 360/n degrees
  return turnSign * turnAngle;
}

std::vector<size_t> findSharedPositions(
    const std::vector<int> &macrocycleRing,
    const std::vector<int> &ring) {
  std::vector<size_t> sharedPositions;
  for (size_t i = 0; i < macrocycleRing.size(); ++i) {
    int atom = macrocycleRing[i];
    if (std::find(ring.begin(), ring.end(), atom) != ring.end()) {
      sharedPositions.push_back(i);
    }
  }
  return sharedPositions;
}

bool verifyAndReorderSharedPositions(
    std::vector<size_t> &sharedPositions,
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

    // If there's exactly one gap, check if it's the wraparound
    if (gapIdx > 0) {
      // Check if positions after gap are contiguous
      bool wrapContiguous = true;

      for (size_t i = gapIdx + 1; i < sharedPositions.size(); ++i) {
        if (sharedPositions[i] != sharedPositions[i - 1] + 1) {
          wrapContiguous = false;
          break;
        }
      }

      // Check wraparound: last position + 1 should equal first position (mod ring size)
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

  // If wraparound was detected, reorder positions to start from the wraparound point
  // E.g., [0, 1, 14] → [14, 0, 1] so constraint position calculation works correctly
  if (gapIdx > 0 && contiguous) {
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
    std::map<size_t, std::vector<EndpointInfo>> &endpointPositions,
    size_t pos,
    int ringSize,
    size_t adjacentInternalPos,
    bool isFirst) {
  EndpointInfo info;
  info.ringSize = ringSize;
  info.adjacentInternalPos = adjacentInternalPos;
  info.isFirst = isFirst;
  endpointPositions[pos].push_back(info);
}

AngleConstraint computeSharedEndpointConstraint(
    const std::vector<EndpointInfo> &endpoints,
    size_t pos) {
  AngleConstraint constraint;
  constraint.position = SIZE_MAX;  // Invalid by default

  if (endpoints.size() != 2) {
    return constraint;  // Only valid for exactly 2 endpoints
  }

  size_t adj1 = endpoints[0].adjacentInternalPos;
  size_t adj2 = endpoints[1].adjacentInternalPos;
  int ringSize1 = endpoints[0].ringSize;
  int ringSize2 = endpoints[1].ringSize;

  // At a shared endpoint, three rings meet: ring1, ring2, and the macrocycle
  // The sum of their internal angles must be 2π (360°)
  // Internal angle = π - |turn_angle|
  // So: (π - |turn1|) + (π - |turn2|) + macrocycle_internal = 2π
  // Therefore: macrocycle_internal = 2π - (π - |turn1|) - (π - |turn2|)
  //                                 = |turn1| + |turn2|
  // And macrocycle_turn = π - macrocycle_internal = π - (|turn1| + |turn2|)

  bool ring1HasInternal = (adj1 != pos);
  bool ring2HasInternal = (adj2 != pos);

  // Get the ideal turn magnitude for each ring
  double turn1 = std::abs(computeIdealAngle(ringSize1, ring1HasInternal ? -1 : 1));
  double turn2 = std::abs(computeIdealAngle(ringSize2, ring2HasInternal ? -1 : 1));

  // Macrocycle internal angle = |turn1| + |turn2|
  double macrocycleInternal = turn1 + turn2;
  // Macrocycle turn angle = π - macrocycle_internal
  double combinedAngle = M_PI - macrocycleInternal;

  constraint.position = pos;
  constraint.targetAngle = combinedAngle;

  return constraint;
}

std::vector<AngleConstraint> identifyAngleConstraintsForFusedRings(
    const std::vector<int> &macrocycleRing,
    const std::vector<std::vector<int>> &allRings) {
  std::vector<AngleConstraint> angleConstraints;
  std::map<size_t, std::vector<EndpointInfo>> endpointPositions;

  std::cerr << ">>> identifyAngleConstraintsForFusedRings: macrocycle_size="
            << macrocycleRing.size() << ", allRings.size=" << allRings.size()
            << std::endl;

  // Process each ring to find fusions with the macrocycle
  int ringIdx = 0;
  for (const auto &ring : allRings) {
    std::cerr << "  Processing ring[" << ringIdx << "]: size=" << ring.size();

    // Skip large rings (only handle small rings 3-7)
    if (ring.size() < 3 || ring.size() > 7) {
      std::cerr << " [SKIP - not small (3-7)]" << std::endl;
      ringIdx++;
      continue;
    }

    // Find shared atoms between this ring and the macrocycle
    std::vector<size_t> sharedPositions = findSharedPositions(macrocycleRing, ring);

    std::cerr << ", shared=" << sharedPositions.size() << " [";
    for (size_t i = 0; i < sharedPositions.size(); ++i) {
      if (i > 0) std::cerr << ", ";
      std::cerr << sharedPositions[i];
    }
    std::cerr << "]";

    // Need at least 1 shared atom (spiro or fusion)
    // Skip if more than 4 shared atoms (unusual/complex fusion)
    if (sharedPositions.size() < 1 || sharedPositions.size() > 4) {
      std::cerr << " [SKIP - invalid shared count]" << std::endl;
      ringIdx++;
      continue;
    }

    // Verify shared atoms are contiguous and reorder if needed
    if (!verifyAndReorderSharedPositions(sharedPositions, macrocycleRing.size())) {
      std::cerr << " [SKIP - not contiguous]" << std::endl;
      ringIdx++;
      continue;
    }

    bool isSmallRing = (ring.size() >= 3 && ring.size() <= 7);
    size_t numShared = sharedPositions.size();

    std::cerr << ", isSmall=" << isSmallRing << std::endl;

    // Handle different fusion patterns
    if (numShared == 2) {
      // RR pattern for 2-atom fusion
      // Track both endpoints for small rings (needed to detect shared endpoints)
      if (isSmallRing) {
        // First endpoint
        size_t firstPos = sharedPositions[0];
        trackEndpoint(endpointPositions, firstPos, ring.size(), firstPos, true);

        // Last endpoint
        size_t lastPos = sharedPositions[1];
        trackEndpoint(endpointPositions, lastPos, ring.size(), lastPos, false);
      }
    } else if (numShared >= 3) {
      // RLR or RLLR pattern
      // Track first endpoint for all small rings
      if (isSmallRing) {
        size_t firstPos = sharedPositions[0];
        size_t adjacentInternalPos = sharedPositions[1];
        trackEndpoint(endpointPositions, firstPos, ring.size(), adjacentInternalPos, true);
      }

      // Add constraints for internal positions (skip first and last)
      std::cerr << "    Adding constraints for internal positions (numShared="
                << numShared << "):" << std::endl;
      for (size_t i = 1; i < numShared - 1; ++i) {
        if (isSmallRing) {
          size_t pos = sharedPositions[i];
          int turnSign = -1;  // L turn
          double idealAngle = computeIdealAngle(ring.size(), turnSign);

          std::cerr << "      pos[" << pos << "]: targetAngle="
                    << (idealAngle * 180.0 / M_PI) << "° (ring size="
                    << ring.size() << ")" << std::endl;

          AngleConstraint angleConstraint;
          angleConstraint.position = pos;
          angleConstraint.targetAngle = idealAngle;
          angleConstraints.push_back(angleConstraint);
        }
      }

      // Track last endpoint for all small rings
      if (isSmallRing) {
        size_t lastPos = sharedPositions[numShared - 1];
        size_t adjacentInternalPos = sharedPositions[numShared - 2];
        trackEndpoint(endpointPositions, lastPos, ring.size(), adjacentInternalPos, false);
      }
    }

    ringIdx++;
  }

  std::cerr << ">>> Processing shared endpoints (tracked "
            << endpointPositions.size() << " positions):" << std::endl;

  // Process shared endpoints: positions that are endpoints of exactly 2 small rings
  for (const auto &[pos, endpoints] : endpointPositions) {
    std::cerr << "  pos[" << pos << "]: " << endpoints.size()
              << " endpoint(s)";
    if (endpoints.size() == 2) {
      AngleConstraint sharedConstraint = computeSharedEndpointConstraint(endpoints, pos);
      if (sharedConstraint.position != SIZE_MAX) {  // Valid constraint
        std::cerr << " -> adding shared endpoint constraint, targetAngle="
                  << (sharedConstraint.targetAngle * 180.0 / M_PI) << "°";
        angleConstraints.push_back(sharedConstraint);
      } else {
        std::cerr << " -> constraint invalid";
      }
    }
    std::cerr << std::endl;
  }

  std::cerr << ">>> Total angle constraints: " << angleConstraints.size()
            << std::endl;

  return angleConstraints;
}

void refineWithJacobian(
    std::vector<RDGeom::Point2D> &coords,
    std::vector<double> &angles,
    const std::set<size_t> &constrainedPositions,
    double bondLength,
    bool isOddRing) {
  // Iterative Jacobian pseudo-inverse refinement to close the gap
  // coords has N+1 elements: [0, 1, ..., N-1, N_dummy]

  size_t N = coords.size() - 1;  // Last element is dummy atom
  double initialGap = (coords[0] - coords[N]).length();

  if (initialGap < 1e-6) {
    return;  // Already closed
  }

  // Cumulative angle adjustments
  std::vector<double> angleAdjustments(N, 0.0);

  // Iterative refinement with damping
  const int maxIters = 3;
  const double damping = 0.8;  // Apply 80% of computed correction
  const double gapThreshold = 0.01;  // Stop if gap < 0.01 Å

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
    const std::vector<AngleConstraint> &angleConstraints,
    double bondLength) {
  size_t N = templateCoords.size();

  // 1. Extract current turn angles from template coordinates
  std::vector<double> currentAngles(N);

  // Initialize with the direction of the last bond (from N-1 to 0)
  // This ensures the first turn angle is calculated correctly for a closed polygon
  size_t prev = N - 1;
  double dx_initial = templateCoords[0].x - templateCoords[prev].x;
  double dy_initial = templateCoords[0].y - templateCoords[prev].y;
  double currentDirection = std::atan2(dy_initial, dx_initial);

  std::cerr << ">>> EXTRACTING ANGLES FROM TEMPLATE (N=" << N << ")" << std::endl;

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

  std::cerr << ">>> EXTRACTED ANGLES (radians and degrees):" << std::endl;
  for (size_t i = 0; i < N; ++i) {
    double degrees = currentAngles[i] * 180.0 / M_PI;
    std::cerr << "  pos[" << i << "]: " << currentAngles[i] << " rad = "
              << degrees << "°" << std::endl;
  }

  // 2. Determine template chirality by counting angle signs
  //    Most prevalent sign = internal/convex turns (R or L depending on direction)
  int positiveCount = 0;
  int negativeCount = 0;
  for (double angle : currentAngles) {
    if (angle > 0) positiveCount++;
    else if (angle < 0) negativeCount++;
  }

  // If majority is negative, template goes CCW (R turns are negative)
  // If majority is positive, template goes CW (L turns are positive)
  bool templateIsCCW = (negativeCount > positiveCount);
  // Constraints are generated assuming CCW with turnSign=-1 for L turns
  // For CCW template: L turns are positive, so flip constraint sign
  // For CW template: L turns are negative, so keep constraint sign
  int signMultiplier = templateIsCCW ? -1 : 1;

  // 3. Apply angle constraints and compute angular difference
  //    Adjust constraint signs based on template chirality
  std::vector<double> adjustedAngles = currentAngles;
  double constraintDelta = 0.0;
  std::set<size_t> constrainedPositions;

  std::cerr << ">>> Template chirality: " << (templateIsCCW ? "CCW" : "CW")
            << ", signMultiplier=" << signMultiplier << std::endl;

  for (const auto &constraint : angleConstraints) {
    size_t pos = constraint.position;
    double oldAngle = adjustedAngles[pos];

    // CRITICAL FIX: Preserve the sign of the template's existing angle
    // Only adjust the magnitude, don't flip the turn direction
    // Template angle sign tells us if this is R (+) or L (-) in this specific template
    double templateSign = (currentAngles[pos] >= 0) ? 1.0 : -1.0;
    double newAngle = templateSign * std::abs(constraint.targetAngle);

    std::cerr << "  Constraint at pos[" << pos << "]: "
              << (oldAngle * 180.0 / M_PI) << "° → "
              << (newAngle * 180.0 / M_PI) << "° (target magnitude was "
              << (std::abs(constraint.targetAngle) * 180.0 / M_PI)
              << "°, preserved template sign)" << std::endl;

    adjustedAngles[pos] = newAngle;
    constraintDelta += (newAngle - oldAngle);
    constrainedPositions.insert(pos);
  }

  std::cerr << ">>> Total constraintDelta = " << (constraintDelta * 180.0 / M_PI)
            << "°" << std::endl;

  // 4. Distribute angular difference across non-constrained angles
  //    to preserve total angular sum
  size_t numFree = N - constrainedPositions.size();
  if (numFree > 0) {
    double adjustmentPerFree = -constraintDelta / numFree;

    std::cerr << ">>> Distributing " << (-constraintDelta * 180.0 / M_PI)
              << "° across " << numFree << " free positions = "
              << (adjustmentPerFree * 180.0 / M_PI) << "° per position"
              << std::endl;

    for (size_t i = 0; i < N; ++i) {
      if (constrainedPositions.count(i) == 0) {
        adjustedAngles[i] += adjustmentPerFree;
      }
    }
  }

  std::cerr << ">>> ADJUSTED ANGLES (after constraints):" << std::endl;
  for (size_t i = 0; i < N; ++i) {
    double degrees = adjustedAngles[i] * 180.0 / M_PI;
    std::cerr << "  pos[" << i << "]: " << adjustedAngles[i] << " rad = "
              << degrees << "°";
    if (constrainedPositions.count(i)) {
      std::cerr << " [CONSTRAINED]";
    }
    std::cerr << std::endl;
  }

  // 5. Generate coordinates from adjusted angles
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

  // 6. Close gap using Jacobian minimization
  refineWithJacobian(coords, adjustedAngles, constrainedPositions,
                     bondLength, N % 2 == 1);

  // 7. Remove dummy atom and return
  coords.pop_back();
  return coords;
}

}  // namespace RDDepict
