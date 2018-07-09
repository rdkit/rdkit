#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <cstdint>
#include <RDGeneral/hash/hash.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

namespace RDKit {
namespace MorganFingerprint {

typedef boost::tuple<boost::dynamic_bitset<>, uint32_t, unsigned int>
    AccumTuple;

MorganAtomInvGenerator::MorganAtomInvGenerator(const bool includeRingMembership)
    : df_includeRingMembership(includeRingMembership) {}

std::vector<std::uint32_t> *MorganAtomInvGenerator::getAtomInvariants(
    const ROMol &mol) const {
  unsigned int nAtoms = mol.getNumAtoms();
  std::vector<std::uint32_t> *atomInvariants =
      new std::vector<std::uint32_t>(nAtoms);
  gboost::hash<std::vector<uint32_t>> vectHasher;
  for (unsigned int i = 0; i < nAtoms; ++i) {
    Atom const *atom = mol.getAtomWithIdx(i);
    std::vector<uint32_t> components;
    components.push_back(atom->getAtomicNum());
    components.push_back(atom->getTotalDegree());
    components.push_back(atom->getTotalNumHs());
    components.push_back(atom->getFormalCharge());
    int deltaMass = static_cast<int>(
        atom->getMass() -
        PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
    components.push_back(deltaMass);

    if (df_includeRingMembership &&
        atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx())) {
      components.push_back(1);
    }
    (*atomInvariants)[i] = vectHasher(components);
  }
  return atomInvariants;
}

std::string MorganAtomInvGenerator::infoString() const {
  return "MorganInvariantGenerator";
}

std::uint64_t MorganArguments::getResultSize() const {
  return std::numeric_limits<uint32_t>::max();
}

MorganArguments::MorganArguments(const unsigned int radius,
                                 const bool countSimulation,
                                 const bool includeChirality,
                                 const bool useBondTypes,
                                 const bool onlyNonzeroInvariants,
                                 const std::vector<std::uint32_t> countBounds,
                                 const std::uint32_t foldedSize)
    : FingerprintArguments(countSimulation, countBounds, foldedSize),
      df_includeChirality(includeChirality),
      df_useBondTypes(useBondTypes),
      df_onlyNonzeroInvariants(onlyNonzeroInvariants),
      d_radius(radius) {}

std::string MorganArguments::infoString() const {
  return "MorganArguments includeChirality=" +
         std::to_string(df_includeChirality) +
         " useBondTypes=" + std::to_string(df_useBondTypes) +
         " onlyNonzeroInvariants=" + std::to_string(df_onlyNonzeroInvariants) +
         " radius=" + std::to_string(d_radius);
}

std::uint32_t MorganAtomEnv::getBitId(
    FingerprintArguments *arguments,
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *bondInvariants,
    const AdditionalOutput *additionalOutput) const {
  if (additionalOutput) {
    // todo: set additional outputs
  }

  return d_code;
}

MorganAtomEnv::MorganAtomEnv(const std::uint32_t code,
                             const unsigned int atomId,
                             const unsigned int layer)
    : d_code(code), d_atomId(atomId), d_layer(layer) {}

std::vector<AtomEnvironment *> MorganEnvGenerator::getEnvironments(
    const ROMol &mol, FingerprintArguments *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *bondInvariants) const {
  PRECONDITION(atomInvariants && (atomInvariants->size() >= mol.getNumAtoms()),
               "bad atom invariants size");
  unsigned int nAtoms = mol.getNumAtoms();
  std::vector<AtomEnvironment *> result = std::vector<AtomEnvironment *>();
  MorganArguments *morganArguments = dynamic_cast<MorganArguments *>(arguments);
  std::vector<uint32_t> *currentAtomInvariants =
      new std::vector<uint32_t>(nAtoms);
  std::vector<uint32_t> *nextAtomInvariants = new std::vector<uint32_t>(nAtoms);
  std::copy(atomInvariants->begin(), atomInvariants->end(),
            currentAtomInvariants->begin());
  std::copy(atomInvariants->begin(), atomInvariants->end(),
            nextAtomInvariants->begin());

  boost::dynamic_bitset<> *chiralAtomsP = new boost::dynamic_bitset<>(nAtoms);

  unsigned int *lastLayer = new unsigned int(0);

  std::vector<uint32_t> currentInvariants(nAtoms);
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

  std::vector<unsigned int> atomOrder(nAtoms);
  if (morganArguments->df_onlyNonzeroInvariants) {
    std::vector<std::pair<int32_t, uint32_t>> ordering;
    for (unsigned int i = 0; i < nAtoms; ++i) {
      if (!currentInvariants[i])
        ordering.push_back(std::make_pair(1, i));
      else
        ordering.push_back(std::make_pair(0, i));
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
        result.push_back(new MorganAtomEnv(currentInvariants[i], i, 0));
      }
    }
  }

  // now do our subsequent rounds:
  for (unsigned int layer = 0; layer < morganArguments->d_radius; ++layer) {
    std::vector<uint32_t> roundInvariants(nAtoms);
    std::vector<boost::dynamic_bitset<>> roundAtomNeighborhoods =
        atomNeighborhoods;
    std::vector<AccumTuple> neighborhoodsThisRound;

    BOOST_FOREACH (unsigned int atomIdx, atomOrder) {
      if (!deadAtoms[atomIdx]) {
        const Atom *tAtom = mol.getAtomWithIdx(atomIdx);
        if (!tAtom->getDegree()) {
          deadAtoms.set(atomIdx, 1);
          continue;
        }
        std::vector<std::pair<int32_t, uint32_t>> nbrs;
        ROMol::OEDGE_ITER beg, end;
        boost::tie(beg, end) = mol.getAtomBonds(tAtom);
        while (beg != end) {
          const Bond *bond = mol[*beg];
          roundAtomNeighborhoods[atomIdx][bond->getIdx()] = 1;

          unsigned int oIdx = bond->getOtherAtomIdx(atomIdx);
          roundAtomNeighborhoods[atomIdx] |= atomNeighborhoods[oIdx];

          int32_t bt = 1;
          if (morganArguments->df_useBondTypes) {
            if (!morganArguments->df_includeChirality ||
                bond->getBondType() != Bond::DOUBLE ||
                bond->getStereo() == Bond::STEREONONE) {
              bt = static_cast<int32_t>(bond->getBondType());
            } else {
              const int32_t stereoOffset = 100;
              const int32_t bondTypeOffset = 10;
              bt = stereoOffset +
                   bondTypeOffset * static_cast<int32_t>(bond->getBondType()) +
                   static_cast<int32_t>(bond->getStereo());
            }
          }
          nbrs.push_back(std::make_pair(bt, currentInvariants[oIdx]));

          ++beg;
        }

        // sort the neighbor list:
        std::sort(nbrs.begin(), nbrs.end());
        // and now calculate the new invariant and test if the atom is newly
        // "chiral"
        boost::uint32_t invar = layer;
        gboost::hash_combine(invar, currentInvariants[atomIdx]);
        bool looksChiral = (tAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
        for (std::vector<std::pair<int32_t, uint32_t>>::const_iterator it =
                 nbrs.begin();
             it != nbrs.end(); ++it) {
          // add the contribution to the new invariant:
          gboost::hash_combine(invar, *it);

          // update our "chirality":
          if (morganArguments->df_includeChirality && looksChiral &&
              chiralAtoms[atomIdx]) {
            if (it->first != static_cast<int32_t>(Bond::SINGLE)) {
              looksChiral = false;
            } else if (it != nbrs.begin() && it->second == (it - 1)->second) {
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
        roundInvariants[atomIdx] = static_cast<uint32_t>(invar);
        neighborhoodsThisRound.push_back(
            boost::make_tuple(roundAtomNeighborhoods[atomIdx],
                              static_cast<uint32_t>(invar), atomIdx));
        if (std::find(neighborhoods.begin(), neighborhoods.end(),
                      roundAtomNeighborhoods[atomIdx]) != neighborhoods.end()) {
          // we have seen this exact environment before, this atom
          // is now out of consideration:
          deadAtoms[atomIdx] = 1;
        }
      }
    }

    std::vector<MorganAtomEnv *> *deadEnvThisRound =
        new std::vector<MorganAtomEnv *>();

    std::sort(neighborhoodsThisRound.begin(), neighborhoodsThisRound.end());
    for (std::vector<AccumTuple>::const_iterator iter =
             neighborhoodsThisRound.begin();
         iter != neighborhoodsThisRound.end(); ++iter) {
      // if we haven't seen this exact environment before, add it to the result
      if (std::find(neighborhoods.begin(), neighborhoods.end(),
                    iter->get<0>()) == neighborhoods.end()) {
        if (!morganArguments->df_onlyNonzeroInvariants ||
            (*atomInvariants)[iter->get<2>()]) {
          if (includeAtoms[iter->get<2>()]) {
            result.push_back(
                new MorganAtomEnv(iter->get<1>(), iter->get<2>(), layer + 1));
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
    std::copy(roundInvariants.begin(), roundInvariants.end(),
              currentInvariants.begin());

    atomNeighborhoods = roundAtomNeighborhoods;
  }

  return result;
}

std::string MorganEnvGenerator::infoString() const {
  return "MorganEnvironmentGenerator";
}

FingerprintGenerator *getMorganGenerator(
    const unsigned int radius, const bool countSimulation,
    const bool includeChirality, const bool useBondTypes,
    const bool onlyNonzeroInvariants,
    const std::vector<std::uint32_t> countBounds,
    const std::uint32_t foldedSize,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator) {
  AtomEnvironmentGenerator *morganEnvGenerator = new MorganEnvGenerator();
  FingerprintArguments *morganArguments = new MorganArguments(
      radius, countSimulation, includeChirality, useBondTypes,
      onlyNonzeroInvariants, countBounds, foldedSize);

  bool ownsAtomInvGenerator = false;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator = new MorganAtomInvGenerator();
    ownsAtomInvGenerator = true;
  }

  return new FingerprintGenerator(
      morganEnvGenerator, morganArguments, atomInvariantsGenerator,
      bondInvariantsGenerator, ownsAtomInvGenerator, false);
}

}  // namespace MorganFingerprint
}  // namespace RDKit