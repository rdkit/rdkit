//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
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

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT RDKitFPArguments
    : public FingerprintArguments<OutputType> {
 public:
  const unsigned int d_minPath;
  const unsigned int d_maxPath;
  const bool df_useHs;
  const bool df_branchedPaths;
  const bool df_useBondOrder;

  OutputType getResultSize() const;

  std::string infoString() const;

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
  std::vector<std::uint32_t> *getAtomInvariants(const ROMol &mol) const;

  std::string infoString() const;
  RDKitFPAtomInvGenerator *clone() const;
};

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT RDKitFPAtomEnv
    : public AtomEnvironment<OutputType> {
  const OutputType d_bitId;
  const boost::dynamic_bitset<> d_atomsInPath;

 public:
  OutputType getBitId(FingerprintArguments<OutputType> *arguments,
                      const std::vector<std::uint32_t> *atomInvariants,
                      const std::vector<std::uint32_t> *bondInvariants,
                      const AdditionalOutput *additionalOutput,
                      bool hashResults = false) const;

  /**
  \brief Construct a new RDKitFPAtomEnv object

  \param bitId bitId generated for this environment
  \param atomsInPath holds atoms in this environment to set additional output

  */
  RDKitFPAtomEnv(const OutputType bitId, boost::dynamic_bitset<> atomsInPath);
};

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT RDKitFPEnvGenerator
    : public AtomEnvironmentGenerator<OutputType> {
 public:
  std::vector<AtomEnvironment<OutputType> *> getEnvironments(
      const ROMol &mol, FingerprintArguments<OutputType> *arguments,
      const std::vector<std::uint32_t> *fromAtoms,
      const std::vector<std::uint32_t> *ignoreAtoms, int confId,
      const AdditionalOutput *additionalOutput,
      const std::vector<std::uint32_t> *atomInvariants,
      const std::vector<std::uint32_t> *bondInvariants,
      bool hashResults = false) const;

  std::string infoString() const;
};

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

 /return FingerprintGenerator<OutputType>* that generated RDKit fingerprints
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