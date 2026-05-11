//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MACROCYCLE_GENERATOR_H
#define RD_MACROCYCLE_GENERATOR_H

#include <vector>
#include <string>
#include <Geometry/point.h>

namespace RDKit {
class ROMol;
}

namespace RDDepict {

//! Hexagonal grid coordinates for tracking macrocycle positions
/*!
  Uses cube coordinates (x, y, z) where x + y + z = 0
  Each bond moves one unit in one of 6 directions (0°, 60°, 120°, 180°, 240°,
  300°) This makes collision detection trivial - just check if two nodes occupy
  same position
*/
struct HexCoord {
  int x, y, z;

  HexCoord() : x(0), y(0), z(0) {}
  HexCoord(int x_, int y_, int z_) : x(x_), y(y_), z(z_) {}

  bool operator==(const HexCoord &other) const {
    return x == other.x && y == other.y && z == other.z;
  }

  bool operator!=(const HexCoord &other) const { return !(*this == other); }

  // For use in std::set or std::map
  bool operator<(const HexCoord &other) const {
    if (x != other.x) return x < other.x;
    if (y != other.y) return y < other.y;
    return z < other.z;
  }
};

//! Type of constraint on turn sequence
enum class ConstraintType {
  SAME,  //!< Two consecutive turns must be same (RR or LL) - e.g., cis double
         //!< bond
  OPPOSITE,  //!< Two consecutive turns must be opposite (RL or LR) - e.g.,
             //!< trans double bond
  FIXED      //!< Fixed pattern of specific turns - e.g., fused rings
};

//! Constraint on turn sequence at a specific position
struct TurnConstraint {
  size_t position;           //!< Starting position in macrocycle (0-indexed)
  ConstraintType type;       //!< Type of constraint
  std::vector<int> pattern;  //!< For FIXED: specific turns (+1=R, -1=L)
                             //!< For SAME/OPPOSITE: empty (applied to position
                             //!< and position+1)
  std::string
      reason;  //!< Description for debugging (e.g., "cis_double_bond_at_5")
};

//! Constraint on angle at a specific position (for regular polygon small rings)
struct AngleConstraint {
  size_t position;     //!< Position in macrocycle (0-indexed)
  double targetAngle;  //!< Target angle in radians
};

//! Generate 2D coordinates for macrocycles using turn-based encoding
/*!
  This class generates macrocycle coordinates by encoding them as a sequence
  of "turn right" (R) or "turn left" (L) decisions. Each step advances by
  bond length and rotates ±60°.

  The algorithm satisfies two closure constraints:
  1. Angular closure: count(R) - count(L) = 6 (ensures 360° rotation)
  2. Positional closure: final position ≈ start position (minimize error)

  Structural constraints (Phase 2):
  - Cis double bonds: two consecutive turns must be same (SAME constraint)
  - Trans double bonds: two consecutive turns must be opposite (OPPOSITE
  constraint)
  - Fused rings: fixed turn patterns (FIXED constraint)
*/
class MacrocycleGenerator {
 public:
  //! Constructor
  /*!
    \param ringSize: Number of atoms in the macrocycle
    \param bondLength: Length of each bond (default 1.5 Å)
    \param useJacobianRefinement: Whether to use Jacobian angle adjustment to
    close gaps (default true)
  */
  MacrocycleGenerator(size_t ringSize, double bondLength = 1.5,
                      bool useJacobianRefinement = true);

  //! Add a structural constraint
  /*!
    \param constraint: Constraint to add to the turn sequence
  */
  void addConstraint(const TurnConstraint &constraint);

  //! Add an angle constraint
  /*!
    \param constraint: Angle constraint for a small ring fusion
  */
  void addAngleConstraint(const AngleConstraint &constraint);

  //! Set substituent sizes at each position (for penalty scoring)
  /*!
    \param substituentSizes: Map from macrocycle position to total substituent
    size \param innerTurnSign: Which turn direction points inward (+1 for R, -1
    for L)
  */
  void setSubstituentInfo(const std::map<size_t, int> &substituentSizes,
                          int innerTurnSign);

  //! Check if a position has an angle constraint
  /*!
    \param position: Position in macrocycle to check
    \return true if position has an angle constraint
  */
  bool hasAngleConstraint(size_t position) const;

  //! Get the constraint angle at a position
  /*!
    \param position: Position in macrocycle
    \return target angle in radians (0 if no constraint)
  */
  double getConstraintAngle(size_t position) const;

  //! Solve for optimal turn sequence
  /*!
    Finds a turn sequence that satisfies angular closure and minimizes
    positional closure error.

    \return true if a valid solution was found
  */
  bool solve();

  //! Generate 2D coordinates from the solved turn sequence
  /*!
    Must call solve() first and ensure it returned true.

    \return vector of 2D coordinates for each atom
  */
  std::vector<RDGeom::Point2D> generateCoordinates() const;

  //! Get the closure error of the current solution
  /*!
    Distance between final position and start position.

    \return closure error in Angstroms
  */
  double getClosureError() const { return d_closureError; }

  //! Get the turn sequence (for debugging/analysis)
  /*!
    \return vector of turns: +1 = R (right), -1 = L (left)
  */
  const std::vector<int> &getTurns() const { return d_turns; }

  //! Calculate hexagonal coordinates for debugging/testing
  /*!
    \param turns: Turn sequence to convert to hex coordinates
    \return vector of hexagonal coordinates for each node
  */
  std::vector<HexCoord> calculateHexCoords(const std::vector<int> &turns) const;

  //! Calculate minimum distance between non-adjacent atoms in hex coordinates
  /*!
    Combines clash detection and minimum distance calculation for efficiency.
    Calculates hex coordinates once and checks all pairs.

    \param turns: Turn sequence to evaluate
    \return minimum distance (in hex units) between any two non-adjacent atoms,
            or 0 if self-crossing detected (clash)
  */
  int calculateMinDistance(const std::vector<int> &turns) const;

  //! Analytically refine coordinates to fix closure gap and geometry
  /*!
    Applies three analytical steps to improve coordinate quality:
    1. Distribute closure gap linearly across all atoms
    2. Adjust bond lengths to exactly bondLength
    3. Distribute angular error evenly across all angles

    Non-iterative - runs in a single pass (~0.1-0.5ms).

    \param coords: Input coordinates (modified in-place)
    \return Final closure error after refinement (should be < 0.02 Å)
  */
  double refineCoordinatesAnalytical(
      std::vector<RDGeom::Point2D> &coords) const;

  //! Adjust angles analytically to close the gap (O(N) least-squares solution)
  /*!
    Uses Jacobian pseudo-inverse to compute minimal angle adjustments that close
    the gap between dummy atom and first atom. Preserves geometry better than
    linear gap distribution.

    \param coords: Input coordinates with dummy atom (N+1 elements)
    \param isOddRing: Whether this is an odd-numbered ring
  */
  void adjustAnglesForClosure(std::vector<RDGeom::Point2D> &coords,
                              bool isOddRing) const;

  //! Get number of free (undecided) turn positions
  /*!
    \return number of positions with turn = 0 (not yet decided)
  */
  size_t getNumFreePositions() const;

 private:
  //! Step 1: Distribute closure gap linearly
  void distributeClosureGap(std::vector<RDGeom::Point2D> &coords) const;

  //! Get current closure gap magnitude
  double getClosureGap(const std::vector<RDGeom::Point2D> &coords) const;
  //! Calculate positional closure error for a given turn sequence
  double calculateClosureError(const std::vector<int> &turns) const;

  //! Calculate penalty for substituents on inner turns
  /*!
    Penalizes substituents placed on inner turns (pointing toward center).
    Penalty is weighted by substituent size.

    \param turns: Turn sequence to evaluate
    \return penalty value (higher = worse)
  */
  double calculateSubstituentPenalty(const std::vector<int> &turns) const;

  //! Apply constraints to turn sequence, return false if constraints conflict
  bool applyConstraints(std::vector<int> &turns) const;

  //! Validate that a turn sequence satisfies all constraints
  bool validateConstraints(const std::vector<int> &turns) const;

  //! Try to solve with a specific R-L target difference
  bool trySolveWithTargetDiff(int targetDiff);

  //! Enumerate combinations and find optimal solution
  bool findOptimalTurnSequence(int numRight, int numLeft);

  //! Helper: generate next combination
  bool nextCombination(std::vector<size_t> &positions, size_t n);

  size_t d_ringSize;         //!< Number of atoms in the ring
  double d_bondLength;       //!< Length of each bond
  std::vector<int> d_turns;  //!< Turn sequence: +1 = R, -1 = L, 0 = undecided
  std::vector<TurnConstraint> d_constraints;  //!< Structural constraints
  std::vector<AngleConstraint>
      d_angleConstraints;  //!< Angle constraints for small rings
  std::map<size_t, int>
      d_substituentSizes;  //!< Map: position -> total substituent size
  int d_innerTurnSign;     //!< Which turn direction points inward (+1 for R, -1
                           //!< for L)
  double d_closureError;   //!< Positional closure error
  bool d_solved;           //!< Whether solve() has been called successfully
  bool d_useJacobianRefinement;  //!< Whether to use Jacobian angle adjustment
};

//! Endpoint information for shared vertex constraint computation
struct EndpointInfo {
  int ringSize;                //!< Size of the small ring
  size_t adjacentInternalPos;  //!< Adjacent internal shared position (for turn
                               //!< direction)
  bool isFirst;  //!< Whether this is the first endpoint of the pattern
};

// ============================================================================
// Helper Functions for Angle Constraint Detection and Refinement
// ============================================================================

//! Compute ideal turn angle for a regular polygon small ring
/*!
  Turn angle (exterior angle) = 360° / n
  For square (n=4): 90°, pentagon (n=5): 72°, hexagon (n=6): 60°

  \param ringSize: Number of atoms in the ring
  \param turnSign: +1 for R turn, -1 for L turn
  \return Turn angle in radians
*/
double computeIdealAngle(int ringSize, int turnSign);

//! Find positions where a small ring shares atoms with the macrocycle
/*!
  \param macrocycleRing: Indices of atoms in the macrocycle ring
  \param ring: Indices of atoms in a small ring
  \return Vector of positions (indices into macrocycleRing) where atoms are
  shared
*/
std::vector<size_t> findSharedPositions(const std::vector<int> &macrocycleRing,
                                        const std::vector<int> &ring);

//! Verify shared positions are contiguous and reorder if needed
/*!
  Checks that shared positions form a contiguous sequence (allowing wraparound).
  Reorders positions to start from the lowest index if valid.

  \param sharedPositions: Vector of shared positions (modified in-place if
  valid) \param macrocycleSize: Size of the macrocycle ring \return true if
  positions are contiguous, false otherwise
*/
bool verifyAndReorderSharedPositions(std::vector<size_t> &sharedPositions,
                                     size_t macrocycleSize);

//! Track endpoint information for shared vertex constraint computation
/*!
  Records information about an endpoint of a shared pattern (RR, RLR, RLLR).

  \param endpointPositions: Map from position to vector of endpoint infos
  \param pos: Position of the endpoint in the macrocycle
  \param ringSize: Size of the small ring
  \param adjacentInternalPos: Adjacent internal shared position (for turn
  direction)
  \param isFirst: Whether this is the first endpoint of the pattern
*/
void trackEndpoint(
    std::map<size_t, std::vector<EndpointInfo>> &endpointPositions, size_t pos,
    int ringSize, size_t adjacentInternalPos, bool isFirst);

//! Compute angle constraint for a shared endpoint (where two small rings meet)
/*!
  Uses the formula: macrocycle_turn = π - (|turn1| + |turn2|)
  where turn1 and turn2 are the exterior angles of the two small rings.

  \param endpoints: Vector of exactly 2 EndpointInfo structures
  \param pos: Position of the shared endpoint
  \return AngleConstraint with computed target angle (position = SIZE_MAX if
  invalid)
*/
AngleConstraint computeSharedEndpointConstraint(
    const std::vector<EndpointInfo> &endpoints, size_t pos);

//! Identify all angle constraints for fused small rings in a macrocycle
/*!
  Analyzes all small rings  fused to the macrocycle and
  generates angle constraints for:
  - ideal polygon angles for internal positions of fused rings
  - Shared endpoints where two fused small rings meet

  \param macrocycleRing: Indices of atoms in the macrocycle ring
  \param allRings: All rings in the molecule
  \return Vector of angle constraints to apply
*/
std::vector<AngleConstraint> identifyAngleConstraintsForFusedRings(
    const std::vector<int> &macrocycleRing,
    const std::vector<std::vector<int>> &allRings);

//! Refine macrocycle coordinates with angle constraints
/*!
  Takes existing coordinates, applies angle constraints, and uses Jacobian
  minimization to close the gap. Preserves total angular sum by distributing
  the constraint-induced difference across non-constrained angles.

  \param templateCoords: Existing 2D coordinates
  \param angleConstraints: Angle constraints to apply
  \param bondLength: Target bond length (default 1.5)
  \return Refined coordinates with constraints applied
*/
std::vector<RDGeom::Point2D> refineMacrocycleWithAngleConstraints(
    const std::vector<RDGeom::Point2D> &templateCoords,
    const std::vector<AngleConstraint> &angleConstraints,
    double bondLength = 1.5);

//! Close gap using Jacobian pseudo-inverse minimization
/*!
  Applies Jacobian-based angle adjustments to close the gap between the
  dummy atom (last coordinate) and the first atom. Preserves constrained angles.

  \param coords: Coordinates with dummy atom (N+1 elements, modified in-place)
  \param angles: Turn angles (modified in-place)
  \param constrainedPositions: Set of positions that should not be adjusted
  \param bondLength: Target bond length
  \param isOddRing: Whether this is an odd-numbered ring
*/
void refineWithJacobian(std::vector<RDGeom::Point2D> &coords,
                        std::vector<double> &angles,
                        const std::set<size_t> &constrainedPositions,
                        double bondLength, bool isOddRing);

//! Flip certain fused rings if it improves the geometry
/*!
  For macrocycles with 4 or 6-membered rings fused at axial positions (1,3 or
  1,4), checks if the middle shared atoms are more substituted than the
  non-shared atoms. If so, reflects the ring across the axial axis to ensure
  substitutions point outward, improving geometry and reducing clashes.

  \param mol: The molecule
  \param macrocycleRing: The macrocycle ring atoms
  \param fusedRings: All fused rings in the system
  \param eatoms: Embedded atom coordinates (modified in-place)
*/
void maybeReflectSymmetricFusedRings(const RDKit::ROMol &mol,
                                     const RDKit::INT_VECT &macrocycleRing,
                                     const RDKit::VECT_INT_VECT &fusedRings,
                                     RDGeom::INT_POINT2D_MAP &eatoms);

}  // namespace RDDepict

#endif
