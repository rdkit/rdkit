//
//  Copyright (C) 2018-2022 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RD_RDFINGERPRINTGEN_H_2018_07
#define RD_RDFINGERPRINTGEN_H_2018_07

#include <GraphMol/Fingerprints/FingerprintGenerator.h>

namespace RDKit {
namespace RDKitFP {

class RDKIT_FINGERPRINTS_EXPORT RDKitFPArguments : public FingerprintArguments {
 public:
  unsigned int d_minPath;
  unsigned int d_maxPath;
  bool df_useHs;
  bool df_branchedPaths;
  bool df_useBondOrder;

  std::string infoString() const override;

  /**
   \brief Construct a new RDKitFPArguments object

   \param minPath the minimum path length (in bonds) to be included
   \param maxPath the maximum path length (in bonds) to be included
   \param useHs toggles inclusion of Hs in paths (if the molecule has
   explicit Hs)
   \param branchedPaths toggles generation of branched subgraphs, not just
   linear paths
   \param useBondOrder toggles inclusion of bond orders in the path hashes
   \param countSimulation         if set, use count simulation while
   generating the fingerprint
   \param countBounds  boundaries for count simulation, corresponding bit will
   be set if the count is higher than the number provided for that spot
   \param fpSize size of the generated fingerprint, does not affect the sparse
   versions
   \param numBitsPerFeature controls the number of bits that are set for each
   path/subgraph found

   */
  RDKitFPArguments(unsigned int minPath, unsigned int maxPath, bool useHs,
                   bool branchedPaths, bool useBondOrder, bool countSimulation,
                   const std::vector<std::uint32_t> countBounds,
                   std::uint32_t fpSize, std::uint32_t numBitsPerFeature);
};

class RDKIT_FINGERPRINTS_EXPORT RDKitFPAtomInvGenerator
    : public AtomInvariantsGenerator {
 public:
  std::vector<std::uint32_t> *getAtomInvariants(
      const ROMol &mol) const override;

  std::string infoString() const override;
  RDKitFPAtomInvGenerator *clone() const override;
};

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT RDKitFPAtomEnv
    : public AtomEnvironment<OutputType> {
  const OutputType d_bitId;
  const boost::dynamic_bitset<> d_atomsInPath;
  const INT_VECT d_bondPath;

 public:
  OutputType getBitId(
      FingerprintArguments *arguments,                   // unused
      const std::vector<std::uint32_t> *atomInvariants,  // unused
      const std::vector<std::uint32_t> *bondInvariants,  // unused
      AdditionalOutput *additionalOutput,                // unused
      bool hashResults = false,                          // unused
      const std::uint64_t fpSize = 0                     // unused
  ) const override;
  void updateAdditionalOutput(AdditionalOutput *output,
                              size_t bitId) const override;

  /**
  \brief Construct a new RDKitFPAtomEnv object

  \param bitId bitId generated for this environment
  \param atomsInPath holds atoms in this environment to set additional output
  \param bondPath the bond path defining the environment

  */
  RDKitFPAtomEnv(const OutputType bitId, boost::dynamic_bitset<> atomsInPath,
                 INT_VECT bondPath)
      : d_bitId(bitId),
        d_atomsInPath(std::move(atomsInPath)),
        d_bondPath(std::move(bondPath)) {}
};

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT RDKitFPEnvGenerator
    : public AtomEnvironmentGenerator<OutputType> {
 public:
  std::vector<AtomEnvironment<OutputType> *> getEnvironments(
      const ROMol &mol, FingerprintArguments *arguments,
      const std::vector<std::uint32_t> *fromAtoms,
      const std::vector<std::uint32_t> *ignoreAtoms, int confId,
      const AdditionalOutput *additionalOutput,
      const std::vector<std::uint32_t> *atomInvariants,
      const std::vector<std::uint32_t> *bondInvariants,
      bool hashResults = false) const override;

  std::string infoString() const override;
  OutputType getResultSize() const override;

};  // namespace RDKitFP

/**
 \brief Get a RDKit fingerprint generator with given parameters

 \tparam OutputType determines the size of the bitIds and the result, can be 32
 or 64 bit unsigned integer
 \param minPath the minimum path length (in bonds) to be included
 \param maxPath the maximum path length (in bonds) to be included
 \param useHs toggles inclusion of Hs in paths (if the molecule has
 explicit Hs)
 \param branchedPaths toggles generation of branched subgraphs, not just
 linear paths
 \param useBondOrder toggles inclusion of bond orders in the path hashes
 \param atomInvariantsGenerator custom atom invariants generator to use
 \param countSimulation         if set, use count simulation while
 generating the fingerprint
 \param countBounds  boundaries for count simulation, corresponding bit will be
 set if the count is higher than the number provided for that spot
 \param fpSize size of the generated fingerprint, does not affect the sparse
 versions
 \param numBitsPerFeature controls the number of bits that are set for each
 path/subgraph found
 \param ownsAtomInvGen  if set atom invariants generator is destroyed with the
 fingerprint generator

 /return FingerprintGenerator<OutputType>* that generates RDKit fingerprints

 This generator supports the following \c AdditionalOutput types:
  - \c atomToBits : which bits each atom is involved in
  - \c atomCounts : how many bits each atom sets
  - \c bitPaths : map from bitId to vectors of bond indices for the individual
 subgraphs

 */
template <typename OutputType>
RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<OutputType> *getRDKitFPGenerator(
    unsigned int minPath = 1, unsigned int maxPath = 7, bool useHs = true,
    bool branchedPaths = true, bool useBondOrder = true,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    bool countSimulation = false,
    const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    std::uint32_t fpSize = 2048, std::uint32_t numBitsPerFeature = 2,
    bool ownsAtomInvGen = false);

}  // namespace RDKitFP
}  // namespace RDKit

#endif
