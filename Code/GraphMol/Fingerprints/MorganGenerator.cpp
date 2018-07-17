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
//#include <algorithm>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace MorganFingerprint {

class ss_matcher {
 public:
  ss_matcher(){};
  ss_matcher(const std::string &pattern) {
    RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    TEST_ASSERT(p);
    m_matcher.reset(p);
  };

  // const RDKit::ROMOL_SPTR &getMatcher() const { return m_matcher; };
  const RDKit::ROMol *getMatcher() const { return m_matcher.get(); };

 private:
  RDKit::ROMOL_SPTR m_matcher;
};

// Definitions for feature points adapted from:
// Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)
const char *smartsPatterns[6] = {
    "[$([N;!H0;v3,v4&+1]),\
$([O,S;H1;+0]),\
n&H1&+0]",                                                  // Donor
    "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
$([O,S;H0;v2]),\
$([O,S;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),\
n&H0&+0,\
$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",                    // Acceptor
    "[a]",                                                  // Aromatic
    "[F,Cl,Br,I]",                                          // Halogen
    "[#7;+,\
$([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",  // Basic
    "[$([C,S](=[O,S,P])-[O;H1,-1])]"                        // Acidic
};
std::vector<std::string> defaultFeatureSmarts(smartsPatterns,
                                              smartsPatterns + 6);
typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
                         boost::flyweights::no_tracking>
    pattern_flyweight;

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
  return "MorganInvariantGenerator includeRingMembership=" +
         std::to_string(df_includeRingMembership);
}

MorganFeatureAtomInvGenerator::MorganFeatureAtomInvGenerator(
    std::vector<const ROMol *> *patterns) {
  if (!patterns) {
    dp_patterns = new std::vector<const ROMol *>();
    dp_patterns->reserve(defaultFeatureSmarts.size());
    for (std::vector<std::string>::const_iterator smaIt =
             defaultFeatureSmarts.begin();
         smaIt != defaultFeatureSmarts.end(); ++smaIt) {
      const ROMol *matcher = pattern_flyweight(*smaIt).get().getMatcher();
      CHECK_INVARIANT(matcher, "bad smarts");
      dp_patterns->push_back(matcher);
    }
    df_ownsData = true;
  } else {
    dp_patterns = patterns;
    df_ownsData = false;
  }
}

MorganFeatureAtomInvGenerator::~MorganFeatureAtomInvGenerator() {
  if (df_ownsData) {
    delete dp_patterns;
  }
}

std::string MorganFeatureAtomInvGenerator::infoString() const {
  return "MorganFeatureInvariantGenerator";
}

std::vector<std::uint32_t> *MorganFeatureAtomInvGenerator::getAtomInvariants(
    const ROMol &mol) const {
  unsigned int nAtoms = mol.getNumAtoms();
  std::vector<std::uint32_t> *result = new std::vector<std::uint32_t>(nAtoms);

  std::fill(result->begin(), result->end(), 0);
  for (unsigned int i = 0; i < dp_patterns->size(); ++i) {
    std::uint32_t mask = 1 << i;
    std::vector<MatchVectType> matchVect;
    // to maintain thread safety, we have to copy the pattern
    // molecules:
    SubstructMatch(mol, ROMol(*(*dp_patterns)[i], true), matchVect);
    for (std::vector<MatchVectType>::const_iterator mvIt = matchVect.begin();
         mvIt != matchVect.end(); ++mvIt) {
      for (const auto &mIt : *mvIt) {
        (*result)[mIt.second] |= mask;
      }
    }
  }
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

std::uint64_t MorganArguments::getResultSize() const {
  return std::numeric_limits<uint32_t>::max();
}

MorganArguments::MorganArguments(const unsigned int radius,
                                 const bool countSimulation,
                                 const bool includeChirality,
                                 const bool onlyNonzeroInvariants,
                                 const std::vector<std::uint32_t> countBounds,
                                 const std::uint32_t foldedSize)
    : FingerprintArguments(countSimulation, countBounds, foldedSize),
      df_includeChirality(includeChirality),
      df_onlyNonzeroInvariants(onlyNonzeroInvariants),
      d_radius(radius) {}

std::string MorganArguments::infoString() const {
  return "MorganArguments includeChirality=" +
         std::to_string(df_includeChirality) +
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

          int32_t bt = static_cast<int32_t>((*bondInvariants)[bond->getIdx()]);

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
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator,
    const std::uint32_t foldedSize,
    const std::vector<std::uint32_t> countBounds) {
  AtomEnvironmentGenerator *morganEnvGenerator = new MorganEnvGenerator();
  FingerprintArguments *morganArguments =
      new MorganArguments(radius, countSimulation, includeChirality,
                          onlyNonzeroInvariants, countBounds, foldedSize);

  bool ownsAtomInvGenerator = false;
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

  return new FingerprintGenerator(
      morganEnvGenerator, morganArguments, atomInvariantsGenerator,
      bondInvariantsGenerator, ownsAtomInvGenerator, ownsBondInvGenerator);
}

}  // namespace MorganFingerprint
}  // namespace RDKit