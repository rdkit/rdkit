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
#include "DepictUtils.h"

namespace RDKit {
class ROMol;
class Bond;
}  // namespace RDKit

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

//! Struct to hold substituent size information
struct SubstituentInfo {
  std::map<size_t, int>
      sizesByPosition;  //!< Map: macrocycle position -> total substituent size
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
  */
  MacrocycleGenerator(size_t ringSize, double bondLength = 1.5);

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
    size
  */
  void setSubstituentInfo(const std::map<size_t, int> &substituentSizes);

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

  //! Print all constraints for debugging
  void printConstraints() const;

  //! Get number of constraints
  size_t getNumConstraints() const;

 private:
  //! Simplify the system by adding constraints for large rings
  /*!
    Reduces enumeration space by fixing positions and adding constraints.

    Two mechanisms are used:
    1. Direct assignment: Big substituents are fixed to outer turns (chemically
    sensible) 2. Constraint addition: OPPOSITE constraints create dependencies
    between positions
       - One position in each pair is "independent" (enumerated)
       - The other is "dependent" (computed as opposite of its pair)

    \param freePositions [in/out]: Initially all free positions, returns only
    independent ones
  */
  void simplifySystem(std::vector<size_t> &freePositions);

  //! Phase 1: Add FIXED constraints for big substituents to outer turns
  /*!
    Instead of directly modifying turns, adds FIXED constraints to d_constraints
    that will be applied by applyConstraints(). This provides symmetry with
    other constraint types.

    \param freePositions: Positions currently free
    \return Number of FIXED constraints added
  */
  size_t addBigSubstituentConstraints(const std::vector<size_t> &freePositions);

  //! Phase 2: Add OPPOSITE constraints to reduce enumeration space
  /*!
    Scans free positions for consecutive pairs in the ring and adds OPPOSITE
    constraints between them.

    \param freePositions: Positions currently free
    \param constraintsNeeded: How many OPPOSITE constraints to try to add
    \return Number of OPPOSITE constraints actually added
  */
  size_t addOppositeConstraintsToReduce(
      const std::vector<size_t> &freePositions, size_t constraintsNeeded);

  //! Build map of dependent positions via OPPOSITE constraints
  /*!
    Analyzes OPPOSITE constraints to determine which positions are dependent
    (their value is determined by another position's value).

    \param freePositions: Positions still free
    \return Map from dependent position to its controlling independent position
  */
  std::map<size_t, size_t> buildDependencyMap(
      const std::vector<size_t> &freePositions) const;

  //! Collect all positions where turn value is 0 (undecided)
  /*!
    \param turns: Turn sequence to scan
    \return Vector of positions with turn value 0
  */
  std::vector<size_t> collectFreePositions(const std::vector<int> &turns) const;
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

  //! Maximum number of independent positions for enumeration
  //! Systems with more free positions will be simplified via constraints
  static constexpr size_t MAX_ENUMERATION_POSITIONS = 15;

  size_t d_ringSize;         //!< Number of atoms in the ring
  double d_bondLength;       //!< Length of each bond
  std::vector<int> d_turns;  //!< Turn sequence: +1 = R, -1 = L, 0 = undecided
  std::vector<TurnConstraint> d_constraints;  //!< Structural constraints
  std::vector<AngleConstraint>
      d_angleConstraints;  //!< Angle constraints for small rings
  std::map<size_t, int>
      d_substituentSizes;  //!< Map: position -> total substituent size
  const int d_innerTurnSign =
      1;                  //!< Which turn direction points inward (+1 for R)
  double d_closureError;  //!< Positional closure error
  bool d_solved;          //!< Whether solve() has been called successfully
};

//! Endpoint information for shared vertex constraint computation
struct EndpointInfo {
  int ringSize;                //!< Size of the small ring
  std::vector<int> ringAtoms;  //!< Atoms in the small ring (to check overlap)
};

//! Atoms around a double bond with stereochemistry
struct DoubleBondStereoAtoms {
  int atom1;           //!< First atom of the double bond
  int atom2;           //!< Second atom of the double bond
  int atom1Neighbor1;  //!< Stereo-controlling neighbor of atom1 (from
                       //!< getStereoAtoms)
  int atom2Neighbor1;  //!< Stereo-controlling neighbor of atom2 (from
                       //!< getStereoAtoms)
  int atom1Neighbor2;  //!< Other neighbor of atom1 (degree > 2), or -1
  int atom2Neighbor2;  //!< Other neighbor of atom2 (degree > 2), or -1
  bool swappedStereo;  //!< Whether stereo interpretation needs to be swapped
  bool valid;          //!< Whether all required atoms were found

  DoubleBondStereoAtoms()
      : atom1(-1),
        atom2(-1),
        atom1Neighbor1(-1),
        atom2Neighbor1(-1),
        atom1Neighbor2(-1),
        atom2Neighbor2(-1),
        swappedStereo(false),
        valid(false) {}
};

//! Information about a small ring fused to a macrocycle
struct FusedRingInfo {
  RDKit::INT_VECT ringAtoms;    //!< All atoms in the fused ring
  RDKit::INT_VECT sharedAtoms;  //!< Atoms shared with macrocycle, in the order
                                //!< they appear in the macrocycle
  size_t ringSize;

  //! For validation: positions of shared atoms in macrocycle sequence
  std::vector<size_t> macrocyclePositions;
};

// ============================================================================
// Helper Functions for Angle Constraint Detection and Refinement
// ============================================================================

//! Extract atoms around a double bond with stereochemistry
/*!
  Identifies all atoms involved in double bond stereochemistry:
  - The two double bond atoms
  - Their stereo-controlling neighbors (from getStereoAtoms)
  - Additional neighbors if degree > 2

  \param bond: The double bond with stereochemistry
  \param mol: The molecule containing the bond
  \return DoubleBondStereoAtoms struct with all atoms identified
*/
DoubleBondStereoAtoms getDoubleBondStereoAtoms(const RDKit::Bond *bond,
                                               const RDKit::ROMol *mol);

//! Compute ideal turn angle for a regular polygon small ring
/*!
  Turn angle (exterior angle) = 360° / n
  For square (n=4): 90°, pentagon (n=5): 72°, hexagon (n=6): 60°

  \param ringSize: Number of atoms in the ring
  \return Turn angle in radians
*/
double computeIdealAngle(int ringSize);

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

//! Generate turn constraints for fused rings
/*!
  Generates ready-to-use TurnConstraint objects for fused ring geometry.
  Pattern: first=external(R), middle=internal(L), last=external(R)

  For small rings: returns 1 constraint with full pattern
  For macrocycles: returns 2 constraints (first and last positions only)

  \param sharedPositions: Positions in macrocycle where fusion occurs (must be
  contiguous) \param ringSize: Size of the fused ring (for labeling) \param
  isMacrocycle: True if the fused ring is a macrocycle \return Vector of
  TurnConstraint objects ready to add to MacrocycleGenerator

  Examples:
  - Small ring, 2 shared atoms: 1 constraint with pattern {R, R}
  - Small ring, 3 shared atoms: 1 constraint with pattern {R, L, R}
  - Macrocycle, 3 shared atoms: 2 constraints (first pos: {R}, last pos: {R})
*/
std::vector<TurnConstraint> generateFusionConstraints(
    const std::vector<size_t> &sharedPositions, size_t ringSize,
    bool isMacrocycle);

//! Track endpoint information for shared vertex constraint computation
/*!
  Records information about an endpoint of a shared pattern (RR, RLR, RLLR).

  \param endpointPositions: Map from position to vector of endpoint infos
  \param pos: Position of the endpoint in the macrocycle
  \param ringSize: Size of the small ring
*/
void trackEndpoint(
    std::map<size_t, std::vector<EndpointInfo>> &endpointPositions, size_t pos,
    int ringSize, const std::vector<int> &ringAtoms);

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

//! Identify angle constraints for triple bonds in macrocycle
/*!
  \param mol: The molecule
  \param macrocycleRing: Vector of atom indices in the macrocycle
  \return Vector of angle constraints for triple bonds (180° angles)
*/
RDKIT_DEPICTOR_EXPORT
std::vector<AngleConstraint> identifyAngleConstraintsForTripleBonds(
    const RDKit::ROMol *mol, const std::vector<int> &macrocycleRing);

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
    double bondLength = RDDepict::BOND_LEN);

//! Close gap using Jacobian pseudo-inverse minimization
/*!
  Applies Jacobian-based angle adjustments to close the gap between the
  dummy atom (last coordinate) and the first atom. Preserves constrained angles.

  \param coords: Coordinates with dummy atom (N+1 elements, modified in-place)
  \param angles: Turn angles (modified in-place)
  \param constrainedPositions: Set of positions that should not be adjusted
  \param bondLength: Target bond length
*/
void refineWithJacobian(std::vector<RDGeom::Point2D> &coords,
                        std::vector<double> &angles,
                        const std::set<size_t> &constrainedPositions,
                        double bondLength);

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
                                     RDGeom::INT_POINT2D_MAP &eatoms,
                                     int currentRingIndex = -1);

//! Refine template-matched macrocycle coordinates with angle constraints
/*!
  For template-matched macrocycles with fused small rings or triple bonds,
  applies angle constraint refinement to ensure ideal geometry.

  \param mol: The molecule
  \param macrocycleRing: Atom indices in the macrocycle
  \param allRings: All rings in the molecule (macrocycle will be filtered out)
  \param coords: Coordinate map (modified in-place)
*/
void maybeRefineTemplateMatchedMacrocycle(const RDKit::ROMol *mol,
                                          const RDKit::INT_VECT &macrocycleRing,
                                          const RDKit::VECT_INT_VECT &allRings,
                                          RDGeom::INT_POINT2D_MAP &coords);

//! Generate de-novo 2D coordinates for a macrocycle using turn-based encoding
/*!
  Uses MacrocycleGenerator to create coordinates that satisfy geometric
  constraints from fused rings, stereochemistry, and substituents.

  \param mol: The molecule
  \param macrocycleRing: Indices of atoms in the macrocycle ring
  \param allRings: All rings in the molecule
  \param substituentSizesByPosition: Map of macrocycle position -> total
  substituent size \param bondLength: Target bond length (default 1.5) \return
  Coordinates for macrocycle atoms (empty if generation failed)
*/
std::vector<RDGeom::Point2D> generateMacrocycleCoordinates(
    const RDKit::ROMol *mol, const RDKit::INT_VECT &macrocycleRing,
    const RDKit::VECT_INT_VECT &allRings,
    const std::map<size_t, int> &substituentSizesByPosition,
    double bondLength = RDDepict::BOND_LEN, int currentRingIndex = -1,
    const RDGeom::INT_POINT2D_MAP *existingCoords = nullptr);

//! Match macrocycle to a template and extract coordinates
/*!
  Attempts to match the macrocycle to a template from the coordinate template
  library. If a match is found, the coordinates are extracted and refined

  \param mol: The molecule
  \param macrocycleRing: Atom indices in the macrocycle
  \param allRings: All rings in the molecule
  \param substituentSizesByPosition: Map of position -> total substituent size
  \param coords: Output coordinate map (populated on success)
  \return true if a template match was found, false otherwise
*/
bool matchToTemplateMacrocycle(
    const RDKit::ROMol *mol, const RDKit::INT_VECT &macrocycleRing,
    const RDKit::VECT_INT_VECT &allRings,
    const std::map<size_t, int> &substituentSizesByPosition,
    RDGeom::INT_POINT2D_MAP &coords, int currentRingIndex = -1);

//! Compute substituent sizes for all positions in a macrocycle
/*!
  \param mol: The molecule
  \param macrocycleRing: Atom indices in the macrocycle
  \param ringAtoms: Bitset of all atoms in rings
  \return SubstituentInfo with size map by position
*/
SubstituentInfo computeSubstituentInfo(
    const RDKit::ROMol *mol, const RDKit::INT_VECT &macrocycleRing,
    const boost::dynamic_bitset<> &ringAtoms);

//! Check if stereochemistry of template matches molecule
/*!
  \param mol: The molecule being matched
  \param templateMol: The template molecule
  \param match: The match vector from substructure matching
  \return true if stereochemistry is compatible
*/
bool checkStereoChemistry(const RDKit::ROMol &mol,
                          const RDKit::ROMol &templateMol,
                          const RDKit::MatchVectType &match);

//! Check if a macrocycle ring vector should be reversed for fusion
/*!
  Determines if the macrocycle atom list should be reversed based on the
  geometric layout (CW vs CCW) of an already-embedded fused macrocycle. This
  ensures that when coordinates are generated for the second macrocycle, which
  is always done CW, the shared atoms will be traversed in the correct order

  \param macrocycleRing: Atom indices in the macrocycle to check
  \param allRings: All rings in the molecule
  \param currentRingIndex: Index of the current macrocycle in allRings
  \param existingCoords: Coordinate map of already-embedded atoms
  \return true if the ring should be reversed, false otherwise
*/
RDKIT_DEPICTOR_EXPORT bool shouldReverseMacrocycleForFusion(
    const RDKit::INT_VECT &macrocycleRing, const RDKit::VECT_INT_VECT &allRings,
    int currentRingIndex, const RDGeom::INT_POINT2D_MAP *existingCoords);

}  // namespace RDDepict

#endif
