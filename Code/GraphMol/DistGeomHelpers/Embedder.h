//
//  Copyright (C) 2004-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RD_EMBEDDER_H_GUARD
#define RD_EMBEDDER_H_GUARD

#include <map>
#include <Geometry/point.h>
#include <GraphMol/ROMol.h>
#include <boost/shared_ptr.hpp>
#include <DistGeom/BoundsMatrix.h>

namespace RDKit {
namespace DGeomHelpers {

//! Parameter object for controlling embedding
/*!
  numConfs       Number of conformations to be generated
  numThreads     Sets the number of threads to use (more than one thread
                 will only be used if the RDKit was build with multithread
                 support) If set to zero, the max supported by the system will
                 be used.
  maxIterations  Max. number of times the embedding will be tried if
                 coordinates are not obtained successfully. The default
                 value is 10x the number of atoms.
  randomSeed     provides a seed for the random number generator (so that
                 the same coordinates can be obtained for a
                 molecule on multiple runs) If -1, the
                 RNG will not be seeded.
  clearConfs     Clear all existing conformations on the molecule
  useRandomCoords  Start the embedding from random coordinates instead of
                   using eigenvalues of the distance matrix.
  boxSizeMult    Determines the size of the box that is used for
                 random coordinates. If this is a positive number, the
                 side length will equal the largest element of the distance
                 matrix times \c boxSizeMult. If this is a negative number,
                 the side length will equal \c -boxSizeMult (i.e. independent
                 of the elements of the distance matrix).
  randNegEig     Picks coordinates at random when a embedding process produces
                 negative eigenvalues
  numZeroFail    Fail embedding if we find this many or more zero eigenvalues
                 (within a tolerance)
  pruneRmsThresh Retain only the conformations out of 'numConfs' after
                 embedding that are at least this far apart from each other.
                 RMSD is computed on the heavy atoms.
                 Prunining is greedy; i.e. the first embedded conformation is
                 retained and from then on only those that are at least
                 \c pruneRmsThresh away from already
                 retained conformations are kept. The pruning is done
                 after embedding and bounds violation minimization.
                 No pruning by default.
  coordMap       a map of int to Point3D, between atom IDs and their locations
                 their locations.  If this container is provided, the
                 coordinates are used to set distance constraints on the
                 embedding. The resulting conformer(s) should have distances
                 between the specified atoms that reproduce those between the
                 points in \c coordMap. Because the embedding produces a
                 molecule in an arbitrary reference frame, an alignment step
                 is required to actually reproduce the provided coordinates.
  optimizerForceTol set the tolerance on forces in the DGeom optimizer
                    (this shouldn't normally be altered in client code).
  ignoreSmoothingFailures  try to embed the molecule even if triangle bounds
                           smoothing fails
  enforceChirality  enforce the correct chirality if chiral centers are present
  useExpTorsionAnglePrefs  impose experimental torsion-angle preferences
  useBasicKnowledge  impose "basic knowledge" terms such as flat
                     aromatic rings, ketones, etc.
  ETversion      version of the experimental torsion-angle preferences
  verbose        print output of experimental torsion-angle preferences
  basinThresh    set the basin threshold for the DGeom force field,
                 (this shouldn't normally be altered in client code).
  onlyHeavyAtomsForRMS  only use the heavy atoms when doing RMS filtering
  boundsMat	custom bound matrix to specify upper and lower bounds of atom
  pairs embedFragmentsSeparately	embed each fragment of molecule in turn
  useSmallRingTorsions	optional torsions to improve small ring conformer
  sampling

  useMacrocycleTorsions	optional torsions to improve macrocycle conformer
  sampling useMacrocycle14config  If 1-4 distances bound heuristics for
  macrocycles is used

  CPCI	custom columbic interactions between atom pairs
*/
struct RDKIT_DISTGEOMHELPERS_EXPORT EmbedParameters {
  unsigned int maxIterations{0};
  int numThreads{1};
  int randomSeed{-1};
  bool clearConfs{true};
  bool useRandomCoords{false};
  double boxSizeMult{2.0};
  bool randNegEig{true};
  unsigned int numZeroFail{1};
  const std::map<int, RDGeom::Point3D> *coordMap{nullptr};
  double optimizerForceTol{1e-3};
  bool ignoreSmoothingFailures{false};
  bool enforceChirality{true};
  bool useExpTorsionAnglePrefs{false};
  bool useBasicKnowledge{false};
  bool verbose{false};
  double basinThresh{5.0};
  double pruneRmsThresh{-1.0};
  bool onlyHeavyAtomsForRMS{false};
  unsigned int ETversion{1};
  boost::shared_ptr<const DistGeom::BoundsMatrix> boundsMat;
  bool embedFragmentsSeparately{true};
  bool useSmallRingTorsions{false};
  bool useMacrocycleTorsions{false};
  bool useMacrocycle14config{false};
  std::shared_ptr<std::map<std::pair<unsigned int, unsigned int>, double>> CPCI;
  EmbedParameters()
      : 
        boundsMat(nullptr),
        
        CPCI(nullptr){};
  EmbedParameters(
      unsigned int maxIterations, int numThreads, int randomSeed,
      bool clearConfs, bool useRandomCoords, double boxSizeMult,
      bool randNegEig, unsigned int numZeroFail,
      const std::map<int, RDGeom::Point3D> *coordMap, double optimizerForceTol,
      bool ignoreSmoothingFailures, bool enforceChirality,
      bool useExpTorsionAnglePrefs, bool useBasicKnowledge, bool verbose,
      double basinThresh, double pruneRmsThresh, bool onlyHeavyAtomsForRMS,
      unsigned int ETversion = 1,
      const DistGeom::BoundsMatrix *boundsMat = nullptr,
      bool embedFragmentsSeparately = true, bool useSmallRingTorsions = false,
      bool useMacrocycleTorsions = false, bool useMacrocycle14config = false,
      std::shared_ptr<std::map<std::pair<unsigned int, unsigned int>, double>>
          CPCI = nullptr)
      : maxIterations(maxIterations),
        numThreads(numThreads),
        randomSeed(randomSeed),
        clearConfs(clearConfs),
        useRandomCoords(useRandomCoords),
        boxSizeMult(boxSizeMult),
        randNegEig(randNegEig),
        numZeroFail(numZeroFail),
        coordMap(coordMap),
        optimizerForceTol(optimizerForceTol),
        ignoreSmoothingFailures(ignoreSmoothingFailures),
        enforceChirality(enforceChirality),
        useExpTorsionAnglePrefs(useExpTorsionAnglePrefs),
        useBasicKnowledge(useBasicKnowledge),
        verbose(verbose),
        basinThresh(basinThresh),
        pruneRmsThresh(pruneRmsThresh),
        onlyHeavyAtomsForRMS(onlyHeavyAtomsForRMS),
        ETversion(ETversion),
        boundsMat(boundsMat),
        embedFragmentsSeparately(embedFragmentsSeparately),
        useSmallRingTorsions(useSmallRingTorsions),
        useMacrocycleTorsions(useMacrocycleTorsions),
        useMacrocycle14config(useMacrocycle14config),
        CPCI(CPCI){};
};

//*! Embed multiple conformations for a molecule
RDKIT_DISTGEOMHELPERS_EXPORT void EmbedMultipleConfs(
    ROMol &mol, INT_VECT &res, unsigned int numConfs,
    const EmbedParameters &params);
inline INT_VECT EmbedMultipleConfs(ROMol &mol, unsigned int numConfs,
                                   const EmbedParameters &params) {
  INT_VECT res;
  EmbedMultipleConfs(mol, res, numConfs, params);
  return res;
}

//! Compute an embedding (in 3D) for the specified molecule using Distance
// Geometry
inline int EmbedMolecule(ROMol &mol, const EmbedParameters &params) {
  INT_VECT confIds;
  EmbedMultipleConfs(mol, confIds, 1, params);

  int res;
  if (confIds.size()) {
    res = confIds[0];
  } else {
    res = -1;
  }
  return res;
}

//! Compute an embedding (in 3D) for the specified molecule using Distance
// Geometry
/*!
  The following operations are performed (in order) here:
   -# Build a distance bounds matrix based on the topology, including 1-5
      distances but not VDW scaling
   -# Triangle smooth this bounds matrix
   -# If step 2 fails - repeat step 1, this time without 1-5 bounds and with vdW
      scaling, and repeat step 2
   -# Pick a distance matrix at random using the bounds matrix
   -# Compute initial coordinates from the distance matrix
   -# Repeat steps 3 and 4 until maxIterations is reached or embedding is
  successful
   -# Adjust initial coordinates by minimizing a Distance Violation error
  function
   **NOTE**: if the molecule has multiple fragments, they will be embedded
  separately,
     this means that they will likely occupy the same region of space.
  \param mol            Molecule of interest
  \param maxIterations  Max. number of times the embedding will be tried if
                        coordinates are not obtained successfully. The default
                        value is 10x the number of atoms.
  \param seed           provides a seed for the random number generator (so that
                        the same coordinates can be obtained for a molecule on
                        multiple runs). If negative, the RNG will not be seeded.
  \param clearConfs     Clear all existing conformations on the molecule
  \param useRandomCoords  Start the embedding from random coordinates instead of
                          using eigenvalues of the distance matrix.
  \param boxSizeMult    Determines the size of the box that is used for
                        random coordinates. If this is a positive number, the
                        side length will equal the largest element of the
                        distance matrix times \c boxSizeMult. If this is a
                        negative number, the side length will equal
                        \c -boxSizeMult (i.e. independent of the elements of the
                        distance matrix).
  \param randNegEig     Picks coordinates at random when a embedding process
                        produces negative eigenvalues
  \param numZeroFail    Fail embedding if we find this many or more zero
                        eigenvalues (within a tolerance)
  \param coordMap  a map of int to Point3D, between atom IDs and their locations
                   their locations.  If this container is provided, the
                   coordinates are used to set distance constraints on the
                   embedding. The resulting conformer(s) should have distances
                   between the specified atoms that reproduce those between the
                   points in \c coordMap. Because the embedding produces a
                   molecule in an arbitrary reference frame, an alignment step
                   is required to actually reproduce the provided coordinates.
  \param optimizerForceTol set the tolerance on forces in the distgeom optimizer
                           (this shouldn't normally be altered in client code).
  \param ignoreSmoothingFailures  try to embed the molecule even if triangle
                                  bounds smoothing fails
  \param enforceChirality  enforce the correct chirality if chiral centers are
                           present
  \param useExpTorsionAnglePrefs  impose experimental torsion-angle preferences
  \param useBasicKnowledge        impose "basic knowledge" terms such as flat
                                  aromatic rings, ketones, etc.
  \param verbose        print output of experimental torsion-angle preferences
  \param basinThresh    set the basin threshold for the DGeom force field,
                        (this shouldn't normally be altered in client code).
  \param onlyHeavyAtomsForRMS  only use the heavy atoms when doing RMS filtering
  \param ETversion	version of torsion preferences to use
  \param useSmallRingTorsions	optional torsions to improve small ring
  conformer sampling

  \param useMacrocycleTorsions	optional torsions to improve macrocycle
  conformer sampling \param useMacrocycle14config  If 1-4 distances bound
  heuristics for macrocycles is used \return ID of the conformations added to
  the molecule, -1 if the emdedding failed
*/
inline int EmbedMolecule(
    ROMol &mol, unsigned int maxIterations = 0, int seed = -1,
    bool clearConfs = true, bool useRandomCoords = false,
    double boxSizeMult = 2.0, bool randNegEig = true,
    unsigned int numZeroFail = 1,
    const std::map<int, RDGeom::Point3D> *coordMap = nullptr,
    double optimizerForceTol = 1e-3, bool ignoreSmoothingFailures = false,
    bool enforceChirality = true, bool useExpTorsionAnglePrefs = false,
    bool useBasicKnowledge = false, bool verbose = false,
    double basinThresh = 5.0, bool onlyHeavyAtomsForRMS = false,
    unsigned int ETversion = 1, bool useSmallRingTorsions = false,
    bool useMacrocycleTorsions = false, bool useMacrocycle14config = false) {
  EmbedParameters params(
      maxIterations, 1, seed, clearConfs, useRandomCoords, boxSizeMult,
      randNegEig, numZeroFail, coordMap, optimizerForceTol,
      ignoreSmoothingFailures, enforceChirality, useExpTorsionAnglePrefs,
      useBasicKnowledge, verbose, basinThresh, -1.0, onlyHeavyAtomsForRMS,
      ETversion, nullptr, true, useSmallRingTorsions, useMacrocycleTorsions,
      useMacrocycle14config);
  return EmbedMolecule(mol, params);
};

//*! Embed multiple conformations for a molecule
/*!
  This is kind of equivalent to calling EmbedMolecule multiple times - just that
  the bounds
  matrix is computed only once from the topology
   **NOTE**: if the molecule has multiple fragments, they will be embedded
  separately,
     this means that they will likely occupy the same region of space.
  \param mol            Molecule of interest
  \param res            Used to return the resulting conformer ids
  \param numConfs       Number of conformations to be generated
  \param numThreads     Sets the number of threads to use (more than one thread
                        will only be used if the RDKit was build with
  multithread
                        support). If set to zero, the max supported by the
  system
                        will be used.
  \param maxIterations  Max. number of times the embedding will be tried if
                        coordinates are not obtained successfully. The default
                        value is 10x the number of atoms.
  \param seed           provides a seed for the random number generator (so that
                        the same coordinates can be obtained for a molecule on
                        multiple runs). If negative, the RNG will not be seeded.
  \param clearConfs     Clear all existing conformations on the molecule
  \param useRandomCoords  Start the embedding from random coordinates instead of
                          using eigenvalues of the distance matrix.
  \param boxSizeMult    Determines the size of the box that is used for
                        random coordinates. If this is a positive number, the
                        side length will equal the largest element of the
                        distance matrix times \c boxSizeMult. If this is a
                        negative number, the side length will equal
                        \c -boxSizeMult (i.e. independent of the elements of the
                        distance matrix).
  \param randNegEig     Picks coordinates at random when a embedding process
                        produces negative eigenvalues
  \param numZeroFail    Fail embedding if we find this many or more zero
                        eigenvalues (within a tolerance)
  \param pruneRmsThresh Retain only the conformations out of 'numConfs' after
                        embedding that are at least this far apart from each
                        other. RMSD is computed on the heavy atoms.
                        Pruning is greedy; i.e. the first embedded conformation
                        is retained and from then on only those that are at
  least
                        pruneRmsThresh away from already retained conformations
                        are kept. The pruning is done after embedding and
                        bounds violation minimization. No pruning by default.
  \param coordMap  a map of int to Point3D, between atom IDs and their locations
                   their locations.  If this container is provided, the
                   coordinates are used to set distance constraints on the
                   embedding. The resulting conformer(s) should have distances
                   between the specified atoms that reproduce those between the
                   points in \c coordMap. Because the embedding produces a
                   molecule in an arbitrary reference frame, an alignment step
                   is required to actually reproduce the provided coordinates.
  \param optimizerForceTol set the tolerance on forces in the DGeom optimizer
                           (this shouldn't normally be altered in client code).
  \param ignoreSmoothingFailures  try to embed the molecule even if triangle
                                  bounds smoothing fails
  \param enforceChirality  enforce the correct chirality if chiral centers are
                           present
  \param useExpTorsionAnglePrefs  impose experimental torsion-angle preferences
  \param useBasicKnowledge        impose "basic knowledge" terms such as flat
                                  aromatic rings, ketones, etc.
  \param verbose        print output of experimental torsion-angle preferences
  \param basinThresh    set the basin threshold for the DGeom force field,
                        (this shouldn't normally be altered in client code).
  \param onlyHeavyAtomsForRMS  only use the heavy atoms when doing RMS filtering
  \param ETversion	version of torsion preferences to use
  \param useSmallRingTorsions	optional torsions to improve small ring
  conformer sampling

  \param useMacrocycleTorsions	optional torsions to improve macrocycle
  conformer sampling \param useMacrocycle14config  If 1-4 distances bound
  heuristics for macrocycles is used

*/
inline void EmbedMultipleConfs(
    ROMol &mol, INT_VECT &res, unsigned int numConfs = 10, int numThreads = 1,
    unsigned int maxIterations = 30, int seed = -1, bool clearConfs = true,
    bool useRandomCoords = false, double boxSizeMult = 2.0,
    bool randNegEig = true, unsigned int numZeroFail = 1,
    double pruneRmsThresh = -1.0,
    const std::map<int, RDGeom::Point3D> *coordMap = nullptr,
    double optimizerForceTol = 1e-3, bool ignoreSmoothingFailures = false,
    bool enforceChirality = true, bool useExpTorsionAnglePrefs = false,
    bool useBasicKnowledge = false, bool verbose = false,
    double basinThresh = 5.0, bool onlyHeavyAtomsForRMS = false,
    unsigned int ETversion = 1, bool useSmallRingTorsions = false,
    bool useMacrocycleTorsions = false, bool useMacrocycle14config = false) {
  EmbedParameters params(
      maxIterations, numThreads, seed, clearConfs, useRandomCoords, boxSizeMult,
      randNegEig, numZeroFail, coordMap, optimizerForceTol,
      ignoreSmoothingFailures, enforceChirality, useExpTorsionAnglePrefs,
      useBasicKnowledge, verbose, basinThresh, pruneRmsThresh,
      onlyHeavyAtomsForRMS, ETversion, nullptr, true, useSmallRingTorsions,
      useMacrocycleTorsions, useMacrocycle14config);
  EmbedMultipleConfs(mol, res, numConfs, params);
};
//! \overload
inline INT_VECT EmbedMultipleConfs(
    ROMol &mol, unsigned int numConfs = 10, unsigned int maxIterations = 30,
    int seed = -1, bool clearConfs = true, bool useRandomCoords = false,
    double boxSizeMult = 2.0, bool randNegEig = true,
    unsigned int numZeroFail = 1, double pruneRmsThresh = -1.0,
    const std::map<int, RDGeom::Point3D> *coordMap = nullptr,
    double optimizerForceTol = 1e-3, bool ignoreSmoothingFailures = false,
    bool enforceChirality = true, bool useExpTorsionAnglePrefs = false,
    bool useBasicKnowledge = false, bool verbose = false,
    double basinThresh = 5.0, bool onlyHeavyAtomsForRMS = false,
    unsigned int ETversion = 1, bool useSmallRingTorsions = false,
    bool useMacrocycleTorsions = false, bool useMacrocycle14config = false) {
  EmbedParameters params(
      maxIterations, 1, seed, clearConfs, useRandomCoords, boxSizeMult,
      randNegEig, numZeroFail, coordMap, optimizerForceTol,
      ignoreSmoothingFailures, enforceChirality, useExpTorsionAnglePrefs,
      useBasicKnowledge, verbose, basinThresh, pruneRmsThresh,
      onlyHeavyAtomsForRMS, ETversion, nullptr, true, useSmallRingTorsions,
      useMacrocycleTorsions, useMacrocycle14config);
  INT_VECT res;
  EmbedMultipleConfs(mol, res, numConfs, params);
  return res;
};

//! Parameters corresponding to Sereina Riniker's KDG approach
RDKIT_DISTGEOMHELPERS_EXPORT extern const EmbedParameters KDG;
//! Parameters corresponding to Sereina Riniker's ETDG approach
RDKIT_DISTGEOMHELPERS_EXPORT extern const EmbedParameters ETDG;
//! Parameters corresponding to Sereina Riniker's ETKDG approach
RDKIT_DISTGEOMHELPERS_EXPORT extern const EmbedParameters ETKDG;
//! Parameters corresponding to Sereina Riniker's ETKDG approach - version 2
RDKIT_DISTGEOMHELPERS_EXPORT extern const EmbedParameters ETKDGv2;
//! Parameters corresponding improved ETKDG by Wang, Witek, Landrum and Riniker
//! (10.1021/acs.jcim.0c00025) - the macrocycle part
RDKIT_DISTGEOMHELPERS_EXPORT extern const EmbedParameters ETKDGv3;
//! Parameters corresponding improved ETKDG by Wang, Witek, Landrum and Riniker
//! (10.1021/acs.jcim.0c00025) - the small ring part
RDKIT_DISTGEOMHELPERS_EXPORT extern const EmbedParameters srETKDGv3;
}  // namespace DGeomHelpers
}  // namespace RDKit

#endif
