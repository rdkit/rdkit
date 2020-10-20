//
//  Copyright (C) 2018 Boran Adas, Google Summer of Code
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <RDGeneral/hash/hash.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/foreach.hpp>

#include <GraphMol/Fingerprints/FingerprintUtil.h>

namespace RDKit {
namespace MorganFingerprint {

using namespace MorganFingerprints;

MorganAtomInvGenerator::MorganAtomInvGenerator(const bool includeRingMembership)
    : df_includeRingMembership(includeRingMembership) {}

std::vector<std::uint32_t> *MorganAtomInvGenerator::getAtomInvariants(
    const ROMol &mol) const {
  unsigned int nAtoms = mol.getNumAtoms();
  std::vector<std::uint32_t> *atomInvariants =
      new std::vector<std::uint32_t>(nAtoms);
  getConnectivityInvariants(mol, *atomInvariants, df_includeRingMembership);
  return atomInvariants;
}

std::string MorganAtomInvGenerator::infoString() const {
  return "MorganInvariantGenerator includeRingMembership=" +
         std::to_string(df_includeRingMembership);
}

MorganAtomInvGenerator *MorganAtomInvGenerator::clone() const {
  return new MorganAtomInvGenerator(df_includeRingMembership);
}

MorganFeatureAtomInvGenerator::MorganFeatureAtomInvGenerator(
    std::vector<const ROMol *> *patterns) {
  dp_patterns = patterns;
}

std::string MorganFeatureAtomInvGenerator::infoString() const {
  return "MorganFeatureInvariantGenerator";
}

MorganFeatureAtomInvGenerator *MorganFeatureAtomInvGenerator::clone() const {
  return new MorganFeatureAtomInvGenerator(dp_patterns);
}

std::vector<std::uint32_t> *MorganFeatureAtomInvGenerator::getAtomInvariants(
    const ROMol &mol) const {
  unsigned int nAtoms = mol.getNumAtoms();
  std::vector<std::uint32_t> *result = new std::vector<std::uint32_t>(nAtoms);

  getFeatureInvariants(mol, *result, dp_patterns);
  return result;
}

MorganBondInvGenerator::MorganBondInvGenerator(const bool useBondTypes,
                                               const bool useChirality)
    : df_useBondTypes(useBondTypes), df_useChirality(useChirality) {}

std::vector<std::uint32_t> *MorganBondInvGenerator::getBondInvariants(
    const ROMol &mol) const {
  std::vector<std::uint32_t> *result =
      new std::vector<std::uint32_t>(mol.getNumBonds());
  for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    Bond const *bond = mol.getBondWithIdx(i);
    int32_t bondInvariant = 1;
    if (df_useBondTypes) {
      if (!df_useChirality || bond->getBondType() != Bond::DOUBLE ||
          bond->getStereo() == Bond::STEREONONE) {
        bondInvariant = static_cast<int32_t>(bond->getBondType());
      } else {
        const int32_t stereoOffset = 100;
        const int32_t bondTypeOffset = 10;
        bondInvariant =
            stereoOffset +
            bondTypeOffset * static_cast<int32_t>(bond->getBondType()) +
            static_cast<int32_t>(bond->getStereo());
      }
    }
    (*result)[bond->getIdx()] = static_cast<int32_t>(bondInvariant);
  }
  return result;
}

std::string MorganBondInvGenerator::infoString() const {
  return "MorganInvariantGenerator useBondTypes=" +
         std::to_string(df_useBondTypes) +
         " useChirality=" + std::to_string(df_useChirality);
}

MorganBondInvGenerator *MorganBondInvGenerator::clone() const {
  return new MorganBondInvGenerator(df_useBondTypes, df_useChirality);
}

template <typename OutputType>
OutputType MorganArguments<OutputType>::getResultSize() const {
  return std::numeric_limits<OutputType>::max();
}

template <typename OutputType>
MorganArguments<OutputType>::MorganArguments(
    const unsigned int radius, const bool countSimulation,
    const bool includeChirality, const bool onlyNonzeroInvariants,
    const std::vector<std::uint32_t> countBounds, const std::uint32_t fpSize)
    : FingerprintArguments<OutputType>(countSimulation, countBounds, fpSize),
      df_includeChirality(includeChirality),
      df_onlyNonzeroInvariants(onlyNonzeroInvariants),
      d_radius(radius) {}

template <typename OutputType>
std::string MorganArguments<OutputType>::infoString() const {
  return "MorganArguments includeChirality=" +
         std::to_string(df_includeChirality) +
         " onlyNonzeroInvariants=" + std::to_string(df_onlyNonzeroInvariants) +
         " radius=" + std::to_string(d_radius);
}

template <typename OutputType>
OutputType MorganAtomEnv<OutputType>::getBitId(
    FingerprintArguments<OutputType> *, // arguments
    const std::vector<std::uint32_t> *, // atomInvariants
    const std::vector<std::uint32_t> *, // bondInvariants
    const AdditionalOutput *additionalOutput,
    const bool // hashResults
) const {
  if (additionalOutput) {
    // todo: set additional outputs
  }

  return d_code;
}

template <typename OutputType>
MorganAtomEnv<OutputType>::MorganAtomEnv(const std::uint32_t code,
                                         const unsigned int atomId,
                                         const unsigned int layer)
    : d_code(code), d_atomId(atomId), d_layer(layer) {}

template <typename OutputType>
std::vector<AtomEnvironment<OutputType> *>
MorganEnvGenerator<OutputType>::getEnvironments(
    const ROMol &mol, FingerprintArguments<OutputType> *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *, //ignoreAtoms
    const int, // confId
    const AdditionalOutput *, // additionalOutput
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *bondInvariants,
    const bool // hashResults
) const {
  PRECONDITION(atomInvariants && (atomInvariants->size() >= mol.getNumAtoms()),
               "bad atom invariants size");
  PRECONDITION(bondInvariants && (bondInvariants->size() >= mol.getNumBonds()),
               "bad bond invariants size");
  unsigned int nAtoms = mol.getNumAtoms();
  std::vector<AtomEnvironment<OutputType> *> result =
      std::vector<AtomEnvironment<OutputType> *>();
  auto *morganArguments =
      dynamic_cast<MorganArguments<OutputType> *>(arguments);

  std::vector<OutputType> currentInvariants(atomInvariants->size());
  std::copy(atomInvariants->begin(), atomInvariants->end(),
            currentInvariants.begin());

  boost::dynamic_bitset<> includeAtoms(nAtoms);
  if (fromAtoms) {
    BOOST_FOREACH (uint32_t idx, *fromAtoms) { includeAtoms.set(idx, 1); }
  } else {
    includeAtoms.set();
  }

  boost::dynamic_bitset<> chiralAtoms(nAtoms);

  // these are the neighborhoods that have already been added
  // to the fingerprint
  std::vector<boost::dynamic_bitset<>> neighborhoods;
  // these are the environments around each atom:
  std::vector<boost::dynamic_bitset<>> atomNeighborhoods(
      nAtoms, boost::dynamic_bitset<>(mol.getNumBonds()));
  boost::dynamic_bitset<> deadAtoms(nAtoms);

  // if df_onlyNonzeroInvariants is set order the atoms to make sure atoms with
  // zero invariants are processed last so that in case of duplicate
  // environments atoms with non-zero invariants are used
  std::vector<unsigned int> atomOrder(nAtoms);
  if (morganArguments->df_onlyNonzeroInvariants) {
    std::vector<std::pair<int32_t, uint32_t>> ordering;
    for (unsigned int i = 0; i < nAtoms; ++i) {
      if (!currentInvariants[i]) {
        ordering.emplace_back(1, i);
      } else {
        ordering.emplace_back(0, i);
      }
    }
    std::sort(ordering.begin(), ordering.end());
    for (unsigned int i = 0; i < nAtoms; ++i) {
      atomOrder[i] = ordering[i].second;
    }
  } else {
    for (unsigned int i = 0; i < nAtoms; ++i) {
      atomOrder[i] = i;
    }
  }

  // add the round 0 invariants to the result
  for (unsigned int i = 0; i < nAtoms; ++i) {
    if (includeAtoms[i]) {
      if (!morganArguments->df_onlyNonzeroInvariants || currentInvariants[i]) {
        result.push_back(
            new MorganAtomEnv<OutputType>(currentInvariants[i], i, 0));
      }
    }
  }

  // now do our subsequent rounds:
  for (unsigned int layer = 0; layer < morganArguments->d_radius; ++layer) {
    // will hold bit ids calculated this round to be used as invariants next
    // round
    std::vector<OutputType> nextLayerInvariants(nAtoms);

    // holds atoms in the environment (neighborhood) for the current layer for
    // each atom, starts with the immediate neighbors of atoms and expands with
    // every iteration
    std::vector<boost::dynamic_bitset<>> roundAtomNeighborhoods =
        atomNeighborhoods;
    std::vector<AccumTuple> allNeighborhoodsThisRound;
    BOOST_FOREACH (unsigned int atomIdx, atomOrder) {
      // skip atoms which will not generate unique environments (neighborhoods)
      // anymore
      if (!deadAtoms[atomIdx]) {
        const Atom *tAtom = mol.getAtomWithIdx(atomIdx);
        if (!tAtom->getDegree()) {
          deadAtoms.set(atomIdx, 1);
          continue;
        }

        // will hold up to date invariants of neighboring atoms with bond types,
        // these invariants hold information from atoms around radius as big as
        // current layer around the current atom
        std::vector<std::pair<int32_t, uint32_t>> neighborhoodInvariants;
        // add up to date invariants of neighbors
        for(auto *bond : tAtom->bonds()) {
          roundAtomNeighborhoods[atomIdx][bond->getIdx()] = 1;

          unsigned int oIdx = bond->getOtherAtomIdx(atomIdx);
          roundAtomNeighborhoods[atomIdx] |= atomNeighborhoods[oIdx];

          auto bt = static_cast<int32_t>((*bondInvariants)[bond->getIdx()]);
          neighborhoodInvariants.push_back(
              std::make_pair(bt, currentInvariants[oIdx]));
        }

        // sort the neighbor list:
        std::sort(neighborhoodInvariants.begin(), neighborhoodInvariants.end());
        // and now calculate the new invariant and test if the atom is newly
        // "chiral"
        std::uint32_t invar = layer;
        gboost::hash_combine(invar, currentInvariants[atomIdx]);
        bool looksChiral = (tAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
        for (std::vector<std::pair<int32_t, uint32_t>>::const_iterator it =
                 neighborhoodInvariants.begin();
             it != neighborhoodInvariants.end(); ++it) {
          // add the contribution to the new invariant:
          gboost::hash_combine(invar, *it);

          // update our "chirality":
          if (morganArguments->df_includeChirality && looksChiral &&
              chiralAtoms[atomIdx]) {
            if (it->first != static_cast<int32_t>(Bond::SINGLE)) {
              looksChiral = false;
            } else if (it != neighborhoodInvariants.begin() &&
                       it->second == (it - 1)->second) {
              looksChiral = false;
            }
          }
        }

        if (morganArguments->df_includeChirality && looksChiral) {
          chiralAtoms[atomIdx] = 1;
          // add an extra value to the invariant to reflect chirality:
          std::string cip = "";
          tAtom->getPropIfPresent(common_properties::_CIPCode, cip);
          if (cip == "R") {
            gboost::hash_combine(invar, 3);
          } else if (cip == "S") {
            gboost::hash_combine(invar, 2);
          } else {
            gboost::hash_combine(invar, 1);
          }
        }

        // this rounds bit id will be next rounds atom invariant, so we save it
        // here
        nextLayerInvariants[atomIdx] = static_cast<OutputType>(invar);

        // store the environment that generated this bit id along with the bit
        // id and the atom id
        allNeighborhoodsThisRound.push_back(
            boost::make_tuple(roundAtomNeighborhoods[atomIdx],
                              static_cast<OutputType>(invar), atomIdx));
        if (std::find(neighborhoods.begin(), neighborhoods.end(),
                      roundAtomNeighborhoods[atomIdx]) != neighborhoods.end()) {
          // we have seen this exact environment before, this atom
          // is now out of consideration:
          deadAtoms[atomIdx] = 1;
        }
      }
    }

    std::sort(allNeighborhoodsThisRound.begin(),
              allNeighborhoodsThisRound.end());
    for (std::vector<AccumTuple>::const_iterator iter =
             allNeighborhoodsThisRound.begin();
         iter != allNeighborhoodsThisRound.end(); ++iter) {
      // if we haven't seen this exact environment before, add it to the result
      if (std::find(neighborhoods.begin(), neighborhoods.end(),
                    iter->get<0>()) == neighborhoods.end()) {
        if (!morganArguments->df_onlyNonzeroInvariants ||
            (*atomInvariants)[iter->get<2>()]) {
          if (includeAtoms[iter->get<2>()]) {
            result.push_back(new MorganAtomEnv<OutputType>(
                iter->get<1>(), iter->get<2>(), layer + 1));
            neighborhoods.push_back(iter->get<0>());
          }
        }
      } else {
        // we have seen this exact environment before, this atom
        // is now out of consideration:
        deadAtoms[iter->get<2>()] = 1;
      }
    }

    // the invariants from this round become the next round invariants:
    std::copy(nextLayerInvariants.begin(), nextLayerInvariants.end(),
              currentInvariants.begin());

    // this rounds calculated neighbors will be next rounds initial neighbors,
    // so the radius can grow every iteration
    atomNeighborhoods = roundAtomNeighborhoods;
  }

  return result;
}

template <typename OutputType>
std::string MorganEnvGenerator<OutputType>::infoString() const {
  return "MorganEnvironmentGenerator";
}

template <typename OutputType>
FingerprintGenerator<OutputType> *getMorganGenerator(
    const unsigned int radius, const bool countSimulation,
    const bool includeChirality, const bool useBondTypes,
    const bool onlyNonzeroInvariants,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator,
    const std::uint32_t fpSize, const std::vector<std::uint32_t> countBounds,
    const bool ownsAtomInvGen, const bool // ownsBondInvGen
) {
  AtomEnvironmentGenerator<OutputType> *morganEnvGenerator =
      new MorganEnvGenerator<OutputType>();
  FingerprintArguments<OutputType> *morganArguments =
      new MorganArguments<OutputType>(radius, countSimulation, includeChirality,
                                      onlyNonzeroInvariants, countBounds,
                                      fpSize);

  bool ownsAtomInvGenerator = ownsAtomInvGen;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator = new MorganAtomInvGenerator();
    ownsAtomInvGenerator = true;
  }

  bool ownsBondInvGenerator = false;
  if (!bondInvariantsGenerator) {
    bondInvariantsGenerator =
        new MorganBondInvGenerator(useBondTypes, includeChirality);
    ownsBondInvGenerator = true;
  }

  return new FingerprintGenerator<OutputType>(
      morganEnvGenerator, morganArguments, atomInvariantsGenerator,
      bondInvariantsGenerator, ownsAtomInvGenerator, ownsBondInvGenerator);
}

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint32_t> *getMorganGenerator(
    const unsigned int radius, const bool countSimulation,
    const bool includeChirality, const bool useBondTypes,
    const bool onlyNonzeroInvariants,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator,
    const std::uint32_t fpSize, const std::vector<std::uint32_t> countBounds,
    const bool ownsAtomInvGen, const bool ownsBondInvGen);

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t> *getMorganGenerator(
    const unsigned int radius, const bool countSimulation,
    const bool includeChirality, const bool useBondTypes,
    const bool onlyNonzeroInvariants,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator,
    const std::uint32_t fpSize, const std::vector<std::uint32_t> countBounds,
    const bool ownsAtomInvGen, const bool ownsBondInvGen);

}  // namespace MorganFingerprint
}  // namespace RDKit
