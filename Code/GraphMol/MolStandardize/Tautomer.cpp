//
//  Copyright (C) 2018-2021 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Tautomer.h"
#include "Fragment.h"
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogUtils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/new_canon.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_set>

#include <utility>

// #define VERBOSE_ENUMERATION 1

#ifdef VERBOSE_ENUMERATION
#include <GraphMol/SmilesParse/SmartsWrite.h>
#endif

namespace RDKit {

namespace MolStandardize {

namespace {

// Count the number of bonds with non-trivial stereo (stereo != STEREONONE).
// Used to break ties when multiple chemically-equivalent tautomers exist
// in the state-key map but differ in preserved stereo information.
unsigned int countBondStereo(const ROMol &mol) {
  unsigned int count = 0;
  for (const auto bond : mol.bonds()) {
    if (bond->getStereo() != Bond::STEREONONE) {
      ++count;
    }
  }
  return count;
}

// Count the number of atoms with non-trivial chiral tag (CHI_UNSPECIFIED).
unsigned int countAtomStereo(const ROMol &mol) {
  unsigned int count = 0;
  for (const auto atom : mol.atoms()) {
    if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
      ++count;
    }
  }
  return count;
}

inline std::uint64_t bondKey(unsigned int b, unsigned int e) {
  const auto lo = std::min(b, e);
  const auto hi = std::max(b, e);
  return (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint64_t>(hi);
}

void copyStereoGroupsByIndices(const ROMol &src, RWMol &dst,
                               bool preferBondIdxLookup) {
  // Copy StereoGroups from src to dst, assuming src and dst represent the same
  // molecular graph and that atom indices refer to the same atoms.
  //
  // Bonds in StereoGroups are stored as Bond* pointers. To translate src-group
  // bonds to dst-group bonds we have two options:
  // - If dst was constructed with the same bond insertion order as src, then
  //   bond indices match and we can do O(1) lookup via getBondWithIdx().
  // - Otherwise we must find the dst bond by its endpoints. We avoid repeated
  //   getBondBetweenAtoms() searches by building a one-time (minIdx,maxIdx)
  //   -> Bond* map for dst.
  //
  // preferBondIdxLookup selects the first fast-path. Even when true, we still
  // fall back to endpoint lookup if the index-based lookup fails.
  //
  // Why this helper lives here (and not in GraphMol/StereoGroup.*):
  // this relies on a tautomer-canonicalization-specific invariant --
  // (identical // atom indices between src/dst)
  const auto &srcGroups = src.getStereoGroups();
  if (srcGroups.empty()) {
    return;
  }

  // Precompute a fast bond lookup by (minAtomIdx,maxAtomIdx) once.
  // This avoids repeated adjacency-list searches in getBondBetweenAtoms().
  std::unordered_map<std::uint64_t, Bond *> bondMap;
  bondMap.reserve(dst.getNumBonds());
  for (auto *b : dst.bonds()) {
    bondMap.emplace(bondKey(b->getBeginAtomIdx(), b->getEndAtomIdx()), b);
  }

  std::vector<StereoGroup> newGroups;
  newGroups.reserve(srcGroups.size());
  for (const auto &group : srcGroups) {
    std::vector<Atom *> atoms;
    atoms.reserve(group.getAtoms().size());
    for (const auto atom : group.getAtoms()) {
      atoms.push_back(dst.getAtomWithIdx(atom->getIdx()));
    }

    std::vector<Bond *> bonds;
    bonds.reserve(group.getBonds().size());
    for (const auto bond : group.getBonds()) {
      Bond *dstBond = nullptr;
      if (preferBondIdxLookup && bond->getIdx() < dst.getNumBonds()) {
        dstBond = dst.getBondWithIdx(bond->getIdx());
      }
      if (!dstBond) {
        const auto it = bondMap.find(
            bondKey(bond->getBeginAtomIdx(), bond->getEndAtomIdx()));
        if (it != bondMap.end()) {
          dstBond = it->second;
        }
      }
      if (dstBond) {
        bonds.push_back(dstBond);
      }
    }

    StereoGroup newGroup(group.getGroupType(), std::move(atoms), std::move(bonds),
                         group.getReadId());
    newGroup.setWriteId(group.getWriteId());
    newGroups.push_back(std::move(newGroup));
  }
  dst.setStereoGroups(std::move(newGroups));
}

// Canonical ordering for state key computation.
// By iterating atoms and bonds in canonical rank order, the state key
// becomes independent of the input atom numbering, making BFS traversal
// deterministic regardless of atom order.
struct CanonicalKeyOrder {
  std::vector<unsigned int> atomOrder;  // atomOrder[canonPos] = atomIdx
  std::vector<unsigned int> atomRank;   // atomRank[atomIdx] = canonPos
  std::vector<unsigned int> bondOrder;  // bondOrder[canonPos] = bondIdx
  std::vector<unsigned int> bondRank;   // bondRank[bondIdx] = canonPos
};

CanonicalKeyOrder computeCanonicalKeyOrder(const ROMol &mol) {
  CanonicalKeyOrder order;
  unsigned int numAtoms = mol.getNumAtoms();
  unsigned int numBonds = mol.getNumBonds();

  // Compute canonical atom ranks
  order.atomRank.resize(numAtoms);
  Canon::rankMolAtoms(mol, order.atomRank);

  // Build inverse: atomOrder[rank] = atomIdx
  order.atomOrder.resize(numAtoms);
  for (unsigned int i = 0; i < numAtoms; ++i) {
    order.atomOrder[order.atomRank[i]] = i;
  }

  // Compute canonical bond order: sort bonds by
  // (min(rank[begin], rank[end]), max(rank[begin], rank[end]))
  std::vector<std::pair<std::pair<unsigned int, unsigned int>, unsigned int>>
      bondPairs;
  bondPairs.reserve(numBonds);
  for (const auto bond : mol.bonds()) {
    unsigned int r1 = order.atomRank[bond->getBeginAtomIdx()];
    unsigned int r2 = order.atomRank[bond->getEndAtomIdx()];
    if (r1 > r2) {
      std::swap(r1, r2);
    }
    bondPairs.push_back({{r1, r2}, bond->getIdx()});
  }
  std::sort(bondPairs.begin(), bondPairs.end());

  order.bondOrder.resize(numBonds);
  order.bondRank.resize(numBonds);
  for (unsigned int pos = 0; pos < numBonds; ++pos) {
    order.bondOrder[pos] = bondPairs[pos].second;
    order.bondRank[bondPairs[pos].second] = pos;
  }

  return order;
}

char encodeBondType(Bond::BondType bt) {
  switch (bt) {
    case Bond::SINGLE:
      return '1';
    case Bond::DOUBLE:
      return '2';
    case Bond::TRIPLE:
      return '3';
    case Bond::AROMATIC:
      return '4';
    default:
      return '0';
  }
}

// Generate a cheap tautomer state key for deduplication.
// Since all tautomers share the same molecular graph, we only need to
// encode what differs: H counts, formal charges, bond orders, and aromaticity.
// This is O(n) vs O(n log n) for canonical SMILES.
//
// Atoms and bonds are iterated in canonical rank order (from `order`),
// making the key independent of input atom numbering.
//
// `collapseAromaticBondOrder` should be false for the normal/stateful keys
// (the defaults used for init/product/reinserted map keys). We only set it
// to true in the pre-sanitize, kekulized fast-path dedup checks (`kmolKey`
// and `perturbedKey`) to avoid representation-dependent aromatic alternation
// creating false "new" states before sanitization/aromaticity perception.
std::string getTautomerStateKey(const ROMol &mol,
                                const CanonicalKeyOrder &order,
                                bool collapseAromaticBondOrder = false) {
  std::string key;
  unsigned int numAtoms = mol.getNumAtoms();
  unsigned int numBonds = mol.getNumBonds();
  // Reserve approximate space: 3 chars per atom + 1 char per bond
  key.reserve(numAtoms * 3 + numBonds + 1);

  // Encode atom state: H count, formal charge, and aromaticity
  for (unsigned int pos = 0; pos < numAtoms; ++pos) {
    const auto atom = mol.getAtomWithIdx(order.atomOrder[pos]);
    // Total H count (0-9 should cover most cases)
    unsigned int totalH = atom->getTotalNumHs();
    key += static_cast<char>('0' + std::min(totalH, 9u));
    // Formal charge (-4 to +4 mapped to '0'-'8')
    int charge = atom->getFormalCharge();
    key += static_cast<char>('4' + std::max(-4, std::min(4, charge)));
    key += atom->getIsAromatic() ? 'a' : 'A';
  }

  // Separator between atoms and bonds
  key += '|';

  for (unsigned int pos = 0; pos < numBonds; ++pos) {
    const auto bond = mol.getBondWithIdx(order.bondOrder[pos]);
    if (collapseAromaticBondOrder && bond->getBeginAtom()->getIsAromatic() &&
        bond->getEndAtom()->getIsAromatic()) {
      key += '4';
    } else {
      key += encodeBondType(bond->getBondType());
    }
  }

  return key;
}

// Compute what the state key WOULD be after applying a tautomer transform,
// without actually modifying the molecule. This allows us to check for
// duplicates before making an expensive molecule copy.
//
// Returns the perturbed key. The baseKey should be from the source molecule
// (kmol). The match, transform describe what would change.
std::string getPerturbedStateKey(const std::string &baseKey, const ROMol &mol,
                                 const MatchVectType &match,
                                 const TautomerTransform &transform,
                                 const CanonicalKeyOrder &order,
                                 bool collapseAromaticBondOrder = false) {
  // Key format: [H-charge-arom per atom] | [bondtype per bond]
  // Each atom takes 3 chars, separator is 1 char, each bond takes 1 char
  const unsigned int numAtoms = mol.getNumAtoms();
  const size_t atomSectionEnd = numAtoms * 3;  // position of '|'

  std::string key = baseKey;

  // The position of atom/bond `idx` in the key is determined by its
  // canonical rank, not its raw index.
  auto atomKeyPos = [&order](unsigned int idx) -> size_t {
    return static_cast<size_t>(order.atomRank[idx]) * 3;
  };
  auto bondKeyPos = [&order, atomSectionEnd](unsigned int idx) -> size_t {
    return atomSectionEnd + 1 +
           static_cast<size_t>(order.bondRank[idx]);
  };

  // Modify H counts for first and last atoms
  int firstIdx = match.front().second;
  int lastIdx = match.back().second;

  // First atom: H count decreases by 1
  size_t firstHPos = atomKeyPos(firstIdx);
  if (key[firstHPos] > '0') {
    key[firstHPos] = static_cast<char>(key[firstHPos] - 1);
  }

  // Last atom: H count increases by 1
  size_t lastHPos = atomKeyPos(lastIdx);
  if (key[lastHPos] < '9') {
    key[lastHPos] = static_cast<char>(key[lastHPos] + 1);
  }

  // Modify formal charges if specified
  if (!transform.Charges.empty()) {
    unsigned int ci = 0;
    for (const auto &pair : match) {
      int chargeAdj = transform.Charges[ci++];
      if (chargeAdj != 0) {
        size_t chargePos = atomKeyPos(pair.second) + 1;
        int newCharge = (key[chargePos] - '4') + chargeAdj;
        key[chargePos] =
            static_cast<char>('4' + std::max(-4, std::min(4, newCharge)));
      }
    }
  }

  // Modify bond orders
  for (size_t i = 0; i < transform.Mol->getNumBonds(); ++i) {
    const auto tbond = transform.Mol->getBondWithIdx(i);
    const Bond *bond = mol.getBondBetweenAtoms(
        match[tbond->getBeginAtomIdx()].second,
        match[tbond->getEndAtomIdx()].second);
    if (!bond) {
      continue;
    }

    // Bond position in key: uses canonical bond rank when ordering is active
    size_t bondPos = bondKeyPos(bond->getIdx());

    if (collapseAromaticBondOrder && bond->getBeginAtom()->getIsAromatic() &&
        bond->getEndAtom()->getIsAromatic()) {
      key[bondPos] = '4';
      continue;
    }

    if (!transform.BondTypes.empty()) {
      key[bondPos] = encodeBondType(transform.BondTypes[i]);
    } else {
      // Flip SINGLE <-> DOUBLE
      if (key[bondPos] == '1') {
        key[bondPos] = '2';
      } else if (key[bondPos] == '2') {
        key[bondPos] = '1';
      }
    }
  }

  return key;
}

}  // namespace

namespace TautomerScoringFunctions {
int scoreRings(const ROMol &mol) {
  int score = 0;
  auto ringInfo = mol.getRingInfo();
  // For aromatic ring scoring, we only need SSSR (not symmetrized SSSR).
  // Aromatic rings are always essential to the ring system and will be in SSSR.
  // Avoiding symmetrizeSSSR and the molecule copy saves significant time
  // since scoreRings is called for each tautomer during canonicalization.
  if (!ringInfo->isSssrOrBetter()) {
    MolOps::findSSSR(mol);
    // ringInfo pointer remains valid - just refresh from the mol
    ringInfo = mol.getRingInfo();
  }
  boost::dynamic_bitset<> isArom(mol.getNumBonds());
  boost::dynamic_bitset<> bothCarbon(mol.getNumBonds());
  for (const auto &bnd : mol.bonds()) {
    if (bnd->getIsAromatic()) {
      isArom.set(bnd->getIdx());
      if (bnd->getBeginAtom()->getAtomicNum() == 6 &&
          bnd->getEndAtom()->getAtomicNum() == 6) {
        bothCarbon.set(bnd->getIdx());
      }
    }
  }
  for (const auto &bring : ringInfo->bondRings()) {
    bool allC = true;
    bool allAromatic = true;
    for (const auto bidx : bring) {
      if (!isArom[bidx]) {
        allAromatic = false;
        break;
      }
      if (!bothCarbon[bidx]) {
        allC = false;
      }
    }
    if (allAromatic) {
      score += 100;
      if (allC) {
        score += 150;
      }
    }
  }
  return score;
};

SubstructTerm::SubstructTerm(std::string aname, std::string asmarts, int ascore,
                             std::vector<int> reqElements,
                             std::string connSmarts)
    : name(std::move(aname)),
      smarts(std::move(asmarts)),
      score(ascore),
      requiredElements(std::move(reqElements)),
      connectivitySmarts(std::move(connSmarts)) {
  std::unique_ptr<ROMol> pattern(SmartsToMol(smarts));
  if (pattern) {
    matcher = std::move(*pattern);
  }
  // Initialize connectivity matcher if provided
  if (!connectivitySmarts.empty()) {
    std::unique_ptr<ROMol> connPattern(SmartsToMol(connectivitySmarts));
    if (connPattern) {
      connectivityMatcher = std::move(*connPattern);
    }
  }
}

const std::vector<SubstructTerm> &getDefaultTautomerScoreSubstructs() {
  // Each term specifies:
  //   - name, SMARTS, score
  //   - requiredElements: atomic numbers that must be present (empty = no filter)
  //   - connectivitySmarts: bond-order-agnostic pattern for pre-screening
  // Since tautomerization only moves H and changes bond orders (never creates/
  // destroys heavy-atom bonds), we can skip patterns whose connectivity
  // prerequisites aren't met by the input molecule.
  static std::vector<SubstructTerm> substructureTerms{
      {"benzoquinone", "[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]", 25, {6}, "[#6]1(~[#6]~[#6]~[#6](~[#6]~[#6]~1)~[N,S,O])~[N,S,O]"},
      {"oxim", "[#6]=[N][OH]", 4, {6, 7, 8}, "[#6]~[#7]~[#8]"},
      {"C=O", "[#6]=,:[#8]", 2, {6, 8}, "[#6]~[#8]"},
      {"N=O", "[#7]=,:[#8]", 2, {7, 8}, "[#7]~[#8]"},
      {"P=O", "[#15]=,:[#8]", 2, {15, 8}, "[#15]~[#8]"},
      {"C=hetero", "[C]=[!#1;!#6]", 1, {6}, "[C]~[!#1;!#6]"},
      {"C(=hetero)-hetero", "[C](=[!#1;!#6])[!#1;!#6]", 2, {6}, "[C](~[!#1;!#6])~[!#1;!#6]"},
      {"aromatic C = exocyclic N", "[c]=!@[N]", -1, {6, 7}, "[c]~[N]"},
      {"methyl", "[CX4H3]", 1, {6}, ""},
      {"guanidine terminal=N", "[#7]C(=[NR0])[#7H0]", 1, {6, 7}, "[#7]~[#6]~[#7]"},
      {"guanidine endocyclic=N", "[#7;R][#6;R]([N])=[#7;R]", 2, {6, 7}, "[#7]~[#6]~[#7]"},
      {"aci-nitro", "[#6]=[N+]([O-])[OH]", -4, {6, 7, 8}, "[#6]~[#7]~[#8]"}};
  return substructureTerms;
}

std::vector<size_t> getRelevantSubstructTermIndices(
    const ROMol &mol, const std::vector<SubstructTerm> &terms) {
  // Prepare relevant SubstructTerms for this molecule by filtering in two
  // stages:
  //   1. Element check: skip terms requiring elements not in the molecule
  //   2. Connectivity check: skip terms whose bond-order-agnostic pattern
  //      doesn't match (since tautomerization doesn't create/destroy bonds)

  std::unordered_set<int> presentElements;
  for (const auto atom : mol.atoms()) {
    presentElements.insert(atom->getAtomicNum());
  }

  std::vector<size_t> relevantIndices;
  relevantIndices.reserve(terms.size());

  SubstructMatchParameters params;
  params.maxMatches = 1;

  for (size_t i = 0; i < terms.size(); ++i) {
    const auto &term = terms[i];

    bool hasAllElements = true;
    for (int elem : term.requiredElements) {
      if (presentElements.find(elem) == presentElements.end()) {
        hasAllElements = false;
        break;
      }
    }
    if (!hasAllElements) {
      continue;
    }
    if (term.connectivityMatcher.getNumAtoms() > 0) {
      auto matches = SubstructMatch(mol, term.connectivityMatcher, params);
      if (matches.empty()) {
        continue;
      }
    }
    relevantIndices.push_back(i);
  }

  return relevantIndices;
}

// ---------------------------------------------------------------------------
// Specialized matchers for simple scoring patterns.
// These avoid VF2 graph matching overhead by directly iterating atoms/bonds.
// Each function is equivalent to the corresponding SMARTS pattern but runs
// in O(n) time without the graph isomorphism overhead.
// ---------------------------------------------------------------------------

// Count bonds between two element types with double or aromatic bond order.
// Matches: [#elem1]=,:[#elem2]
inline unsigned int countDoubleOrAromaticBonds(const ROMol &mol, int elem1,
                                               int elem2) {
  unsigned int count = 0;
  for (const auto bond : mol.bonds()) {
    const auto bt = bond->getBondType();
    if (bt != Bond::DOUBLE && bt != Bond::AROMATIC) {
      continue;
    }
    const int a1 = bond->getBeginAtom()->getAtomicNum();
    const int a2 = bond->getEndAtom()->getAtomicNum();
    if ((a1 == elem1 && a2 == elem2) || (a1 == elem2 && a2 == elem1)) {
      ++count;
    }
  }
  return count;
}

// Count methyl groups: [CX4H3] - sp3 carbon with exactly 3 total H
inline unsigned int countMethyls(const ROMol &mol) {
  unsigned int count = 0;
  for (const auto atom : mol.atoms()) {
    if (atom->getAtomicNum() != 6) {
      continue;
    }
    // X4 means total degree 4 (including implicit H)
    if (atom->getTotalDegree() != 4) {
      continue;
    }
    // H3 means exactly 3 total hydrogens
    if (atom->getTotalNumHs() != 3) {
      continue;
    }
    ++count;
  }
  return count;
}

// Count [C]=[!#1;!#6] - aliphatic carbon double-bonded to heteroatom
inline unsigned int countCarbonDoubleHetero(const ROMol &mol) {
  unsigned int count = 0;
  for (const auto bond : mol.bonds()) {
    if (bond->getBondType() != Bond::DOUBLE) {
      continue;
    }
    const auto *a1 = bond->getBeginAtom();
    const auto *a2 = bond->getEndAtom();
    // [C] is aliphatic carbon (not aromatic)
    // Check C double-bonded to heteroatom (not H, not C)
    auto isAliphaticCarbonToHetero = [](const Atom *c, const Atom *het) {
      return c->getAtomicNum() == 6 && !c->getIsAromatic() &&
             het->getAtomicNum() != 1 && het->getAtomicNum() != 6;
    };
    if (isAliphaticCarbonToHetero(a1, a2) ||
        isAliphaticCarbonToHetero(a2, a1)) {
      ++count;
    }
  }
  return count;
}

// Count [c]=!@[N] - aromatic carbon double-bonded to nitrogen (exocyclic)
inline unsigned int countAromaticCarbonExocyclicN(const ROMol &mol) {
  unsigned int count = 0;
  const auto *ringInfo = mol.getRingInfo();
  for (const auto bond : mol.bonds()) {
    if (bond->getBondType() != Bond::DOUBLE) {
      continue;
    }
    // !@ means bond not in ring
    if (ringInfo->numBondRings(bond->getIdx()) > 0) {
      continue;
    }
    const auto *a1 = bond->getBeginAtom();
    const auto *a2 = bond->getEndAtom();
    // [c] is aromatic carbon, [N] is any nitrogen
    auto isAromaticCarbonToN = [](const Atom *c, const Atom *n) {
      return c->getAtomicNum() == 6 && c->getIsAromatic() &&
             n->getAtomicNum() == 7;
    };
    if (isAromaticCarbonToN(a1, a2) || isAromaticCarbonToN(a2, a1)) {
      ++count;
    }
  }
  return count;
}

// Pattern indices in getDefaultTautomerScoreSubstructs().
// These must match the order of patterns in that function.
// Used for O(1) dispatch to specialized matchers.
namespace {
enum PatternIdx {
  kBenzoquinone = 0,
  kOxim = 1,
  kCarbonylO = 2,
  kNO = 3,
  kPO = 4,
  kCHetero = 5,
  kCHeteroHetero = 6,
  kAromaticCN = 7,
  kMethyl = 8,
  kGuanidineTerminal = 9,
  kGuanidineEndocyclic = 10,
  kAciNitro = 11
};
}  // namespace

int scoreSubstructsFiltered(const ROMol &mol,
                            const std::vector<SubstructTerm> &terms,
                            const std::vector<size_t> &relevantIndices) {
  int score = 0;
  SubstructMatchParameters params;

  for (const size_t idx : relevantIndices) {
    const auto &term = terms[idx];
    if (!term.matcher.getNumAtoms()) {
      continue;
    }

    unsigned int nMatches = 0;

    // Use specialized matchers for simple patterns, VF2 for complex ones
    switch (idx) {
      case kCarbonylO:
        nMatches = countDoubleOrAromaticBonds(mol, 6, 8);  // C=O
        break;
      case kNO:
        nMatches = countDoubleOrAromaticBonds(mol, 7, 8);  // N=O
        break;
      case kPO:
        nMatches = countDoubleOrAromaticBonds(mol, 15, 8);  // P=O
        break;
      case kCHetero:
        nMatches = countCarbonDoubleHetero(mol);
        break;
      case kAromaticCN:
        nMatches = countAromaticCarbonExocyclicN(mol);
        break;
      case kMethyl:
        nMatches = countMethyls(mol);
        break;
      default:
        // Fall back to VF2 for complex patterns (benzoquinone, oxim,
        // C(=hetero)-hetero, guanidine, aci-nitro)
        nMatches = SubstructMatchCount(mol, term.matcher, params);
        break;
    }

    score += static_cast<int>(nMatches) * term.score;
  }
  return score;
}

int scoreSubstructs(const ROMol &mol,
                    const std::vector<SubstructTerm> &substructureTerms) {
  int score = 0;
  SubstructMatchParameters params;
  for (const auto &term : substructureTerms) {
    if (!term.matcher.getNumAtoms()) {
      BOOST_LOG(rdErrorLog) << " matcher for term " << term.name
                            << " is invalid, ignoring it." << std::endl;
      continue;
    }
    const auto nMatches = SubstructMatchCount(mol, term.matcher, params);
    score += static_cast<int>(nMatches) * term.score;
  }
  return score;
}

int scoreHeteroHs(const ROMol &mol) {
  int score = 0;
  for (const auto &at : mol.atoms()) {
    int anum = at->getAtomicNum();
    if (anum == 15 || anum == 16 || anum == 34 || anum == 52) {
      score -= at->getTotalNumHs();
    }
  }
  return score;
}
}  // namespace TautomerScoringFunctions

TautomerEnumerator::TautomerEnumerator(const CleanupParameters &params)
    : d_maxTautomers(params.maxTautomers),
      d_maxTransforms(params.maxTransforms),
      d_removeSp3Stereo(params.tautomerRemoveSp3Stereo),
      d_removeBondStereo(params.tautomerRemoveBondStereo),
      d_removeIsotopicHs(params.tautomerRemoveIsotopicHs),
      d_reassignStereo(params.tautomerReassignStereo) {
  std::unique_ptr<TautomerCatalogParams> tautParams;
  if (params.tautomerTransformData.empty()) {
    tautParams.reset(new TautomerCatalogParams(params.tautomerTransforms));
  } else {
    tautParams.reset(new TautomerCatalogParams(params.tautomerTransformData));
  }
  dp_catalog.reset(new TautomerCatalog(tautParams.get()));
}

bool TautomerEnumerator::setTautomerStereoAndIsoHs(
    const ROMol &mol, ROMol &taut, const TautomerEnumeratorResult &res) const {
  bool modified = false;
  // Iterate only the atoms/bonds actually modified by transforms.
  for (auto atomIdx = res.d_modifiedAtoms.find_first();
       atomIdx != boost::dynamic_bitset<>::npos;
       atomIdx = res.d_modifiedAtoms.find_next(atomIdx)) {
    const auto atom = mol.getAtomWithIdx(static_cast<unsigned int>(atomIdx));
    auto tautAtom = taut.getAtomWithIdx(atomIdx);
    // clear chiral tag on sp2 atoms (also sp3 if d_removeSp3Stereo is true)
    if (tautAtom->getHybridization() == Atom::SP2 || d_removeSp3Stereo) {
      modified |= (tautAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
      tautAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
      if (tautAtom->hasProp(common_properties::_CIPCode)) {
        tautAtom->clearProp(common_properties::_CIPCode);
      }
    } else {
      modified |= (tautAtom->getChiralTag() != atom->getChiralTag());
      tautAtom->setChiralTag(atom->getChiralTag());
      if (atom->hasProp(common_properties::_CIPCode)) {
        tautAtom->setProp(
            common_properties::_CIPCode,
            atom->getProp<std::string>(common_properties::_CIPCode));
      }
    }
    // remove isotopic Hs if present (and if d_removeIsotopicHs is true)
    if (tautAtom->hasProp(common_properties::_isotopicHs) &&
        (d_removeIsotopicHs || !tautAtom->getTotalNumHs())) {
      tautAtom->clearProp(common_properties::_isotopicHs);
    }
  }
  // remove stereochemistry on bonds that are part of a tautomeric path
  // Build bond lookup caches to avoid O(n) getBondWithIdx calls.
  // getBondWithIdx iterates through all bonds to find the one with the given
  // index, which is expensive when called multiple times in a loop.
  std::vector<const Bond *> molBonds;
  std::vector<Bond *> tautBonds;
  if (res.d_modifiedBonds.any()) {
    const auto numBonds = mol.getNumBonds();
    molBonds.resize(numBonds);
    tautBonds.resize(numBonds);
    for (auto bond : mol.bonds()) {
      molBonds[bond->getIdx()] = bond;
    }
    for (auto bond : taut.bonds()) {
      tautBonds[bond->getIdx()] = bond;
    }
  }
  for (auto bondIdx = res.d_modifiedBonds.find_first();
       bondIdx != boost::dynamic_bitset<>::npos;
       bondIdx = res.d_modifiedBonds.find_next(bondIdx)) {
    const auto bond = molBonds[bondIdx];
    std::vector<unsigned int> bondsToClearDirs;
    if (bond->getBondType() == Bond::DOUBLE &&
        bond->getStereo() > Bond::STEREOANY) {
      // look around the beginning and end atoms and check for bonds with
      // direction set
      for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
        for (const auto &nbri :
             boost::make_iterator_range(mol.getAtomBonds(atom))) {
          const auto &obnd = mol[nbri];
          if (obnd->getBondDir() == Bond::ENDDOWNRIGHT ||
              obnd->getBondDir() == Bond::ENDUPRIGHT) {
            bondsToClearDirs.push_back(obnd->getIdx());
          }
        }
      }
    }
    auto tautBond = tautBonds[bondIdx];
    if (tautBond->getBondType() != Bond::DOUBLE || d_removeBondStereo) {
      modified |= (tautBond->getStereo() != Bond::STEREONONE);
      tautBond->setStereo(Bond::STEREONONE);
      tautBond->getStereoAtoms().clear();
      for (auto bi : bondsToClearDirs) {
        tautBonds[bi]->setBondDir(Bond::NONE);
      }
    } else {
      const INT_VECT &sa = bond->getStereoAtoms();
      modified |= (tautBond->getStereo() != bond->getStereo() ||
                   sa.size() != tautBond->getStereoAtoms().size());
      if (sa.size() == 2) {
        tautBond->setStereoAtoms(sa.front(), sa.back());
      }
      tautBond->setStereo(bond->getStereo());
      for (auto bi : bondsToClearDirs) {
        tautBonds[bi]->setBondDir(molBonds[bi]->getBondDir());
      }
    }
  }
  if (d_reassignStereo) {
    static const bool cleanIt = true;
    static const bool force = true;
    MolOps::assignStereochemistry(taut, cleanIt, force);
  } else {
    taut.setProp(common_properties::_StereochemDone, 1);
  }
  return modified;
}

std::vector<ROMOL_SPTR> TautomerEnumerator::enumerate(
    const ROMol &mol, boost::dynamic_bitset<> *modifiedAtoms,
    boost::dynamic_bitset<> *modifiedBonds) const {
  TautomerEnumeratorResult tresult = enumerate(mol);
  if (modifiedAtoms) {
    *modifiedAtoms = tresult.modifiedAtoms();
  }
  if (modifiedBonds) {
    *modifiedBonds = tresult.modifiedBonds();
  }
  return tresult.tautomers();
}

TautomerEnumeratorResult TautomerEnumerator::enumerate(const ROMol &mol) const {
#ifdef VERBOSE_ENUMERATION
  std::cout << "**********************************" << std::endl;
#endif
  PRECONDITION(dp_catalog, "no catalog!");
  const TautomerCatalogParams *tautparams = dp_catalog->getCatalogParams();
  PRECONDITION(tautparams, "");

  TautomerEnumeratorResult res;

  // Additional fast path: track state keys computed *before* the expensive
  // sanitize steps (Kekulize/setAromaticity/...) so we can skip redundant
  // candidates early.
  std::unordered_set<std::string> preSanitizeStateKeys;
  preSanitizeStateKeys.reserve(d_maxTautomers * 2);

  const std::vector<TautomerTransform> &transforms =
      tautparams->getTransforms();

  // Enumerate all possible tautomers and return them as a vector.
  // taut is a copy of the input molecule
  ROMOL_SPTR taut(new ROMol(mol));
  // do whatever sanitization bits are required
  if (taut->needsUpdatePropertyCache()) {
    taut->updatePropertyCache(false);
  }
  if (!taut->getRingInfo()->isSymmSssr()) {
    MolOps::symmetrizeSSSR(*taut);
  }

  // Compute a canonical atom/bond ordering for state key computation.
  // This makes the BFS traversal deterministic regardless of input atom
  // numbering, ensuring enumerate() discovers the same tautomer set for
  // the same molecule given in any atom order.
  CanonicalKeyOrder canonOrder = computeCanonicalKeyOrder(*taut);

  // Kekulized form will be created lazily when needed for transform matching
  // NOTE: we key the internal map by a cheap "state key" rather than SMILES.
  // This is an experiment to avoid per-tautomer canonical SMILES generation
  // during enumeration. See perf/tautomer-canonicalization.md for notes.
  std::string initKey = getTautomerStateKey(*taut, canonOrder);
  res.d_tautomers = {{initKey, Tautomer(taut, 0, 0)}};
  res.d_modifiedAtoms.resize(mol.getNumAtoms());
  res.d_modifiedBonds.resize(mol.getNumBonds());

  // Keep running counts of modified atoms/bonds.
  // `boost::dynamic_bitset<>::count()` is O(n) in the number of blocks, and we
  // were previously calling it once per new tautomer, which is avoidable.
  size_t numModifiedAtoms = 0;
  size_t numModifiedBonds = 0;
  const auto markAtomModified = [&res, &numModifiedAtoms](unsigned int idx) {
    if (!res.d_modifiedAtoms.test(idx)) {
      res.d_modifiedAtoms.set(idx);
      ++numModifiedAtoms;
    }
  };
  const auto markBondModified = [&res, &numModifiedBonds](unsigned int idx) {
    if (!res.d_modifiedBonds.test(idx)) {
      res.d_modifiedBonds.set(idx);
      ++numModifiedBonds;
    }
  };

  bool completed = false;
  bool bailOut = false;
  unsigned int nTransforms = 0;
  static const std::array<const char *, 4> statusMsg{
      "completed", "max tautomers reached", "max transforms reached",
      "canceled"};

  while (!completed && !bailOut) {
    // std::map automatically sorts res.d_tautomers into alphabetical order
    // (state key)
    for (auto &smilesTautomerPair : res.d_tautomers) {
#ifdef VERBOSE_ENUMERATION
      std::cout << "Current tautomers: " << std::endl;
      for (const auto &smilesTautomerPair : res.d_tautomers) {
        std::cout << smilesTautomerPair.first << " done "
                  << smilesTautomerPair.second.d_done << std::endl;
      }
#endif
      if (smilesTautomerPair.second.d_done) {
#ifdef VERBOSE_ENUMERATION
        std::cout << "Skipping " << smilesTautomerPair.first
                  << " as already done" << std::endl;
#endif
        continue;
      }
#ifdef VERBOSE_ENUMERATION
      std::cout << "Looking at tautomer: " << smilesTautomerPair.first
                << std::endl;
#endif
      // tautomer not yet done
      for (const auto &transform : transforms) {
        if (bailOut) {
          break;
        }
        // kmol is the kekulized version of the tautomer (created lazily)
        const auto &kmol = smilesTautomerPair.second.getKekulized();
        std::vector<MatchVectType> matches;
        unsigned int matched = SubstructMatch(*kmol, *(transform.Mol), matches);

        if (!matched) {
          continue;
        }
        ++nTransforms;
#ifdef VERBOSE_ENUMERATION
        std::string name;
        (transform.Mol)->getProp(common_properties::_Name, name);
        SmilesWriteParams smilesWriteParams;
        smilesWriteParams.allBondsExplicit = true;
        std::cout << "kmol for " << smilesTautomerPair.first << " : "
                  << MolToSmiles(*kmol, smilesWriteParams) << std::endl;
        std::cout << "transform mol: " << MolToSmarts(*(transform.Mol))
                  << std::endl;

        std::cout << "Matched: " << name << std::endl;
#endif
        // Compute the base state key from kmol once, before the match loop.
        // This will be used to compute perturbed keys for each match.
        // Use aromatic-collapsed keying only in this pre-sanitize fast path
        // to avoid representation-dependent aromatic alternation duplicates.
        std::string kmolKey = getTautomerStateKey(*kmol, canonOrder, true);

        // loop over transform matches
        for (const auto &match : matches) {
          if (nTransforms >= d_maxTransforms) {
            res.d_status = TautomerEnumeratorStatus::MaxTransformsReached;
            bailOut = true;
          } else if (res.d_tautomers.size() >= d_maxTautomers) {
            res.d_status = TautomerEnumeratorStatus::MaxTautomersReached;
            bailOut = true;
          } else if (d_callback.get() && !(*d_callback)(mol, res)) {
            res.d_status = TautomerEnumeratorStatus::Canceled;
            bailOut = true;
          }
          if (bailOut) {
            break;
          }

          // Compute what the state key WOULD be after applying this transform,
          // WITHOUT copying the molecule. This allows us to skip duplicate
          // tautomers before the expensive molecule copy.
          std::string perturbedKey =
              // Keep aromatic collapsing aligned with kmolKey in this
              // pre-sanitize duplicate check path only.
              getPerturbedStateKey(kmolKey, *kmol, match, transform,
                        canonOrder, true);
          bool isNewPerturbed = preSanitizeStateKeys.insert(perturbedKey).second;
          if (!isNewPerturbed) {
#ifdef VERBOSE_ENUMERATION
            std::cout << "Previous tautomer state seen again (pre-copy check)"
                      << std::endl;
#endif
            continue;
          }

          // This is a potentially new tautomer - now create the copy
          RWMOL_SPTR product(new RWMol(*kmol, true));
          // Remove a hydrogen from the first matched atom and add one to the
          // last
          int firstIdx = match.front().second;
          int lastIdx = match.back().second;
          Atom *first = product->getAtomWithIdx(firstIdx);
          Atom *last = product->getAtomWithIdx(lastIdx);
          markAtomModified(static_cast<unsigned int>(firstIdx));
          markAtomModified(static_cast<unsigned int>(lastIdx));
          first->setNumExplicitHs(
              std::max(0, static_cast<int>(first->getTotalNumHs()) - 1));
          last->setNumExplicitHs(last->getTotalNumHs() + 1);
          // Remove any implicit hydrogens from the first and last atoms
          // now we have set the count explicitly
          first->setNoImplicit(true);
          last->setNoImplicit(true);

          // Adjust bond orders
          unsigned int bi = 0;
          for (size_t i = 0; i < transform.Mol->getNumBonds(); ++i) {
            const auto tbond = transform.Mol->getBondWithIdx(i);
            Bond *bond = product->getBondBetweenAtoms(
                match[tbond->getBeginAtomIdx()].second,
                match[tbond->getEndAtomIdx()].second);
            ASSERT_INVARIANT(bond, "required bond not found");
            // check if bonds is specified in tautomer.in file
            if (!transform.BondTypes.empty()) {
              bond->setBondType(transform.BondTypes[bi]);
              ++bi;
            } else {
              Bond::BondType bondtype = bond->getBondType();
#ifdef VERBOSE_ENUMERATION
              std::cout << "Bond as double: " << bond->getBondTypeAsDouble()
                        << std::endl;
              std::cout << bondtype << std::endl;
#endif
              if (bondtype == Bond::SINGLE) {
                bond->setBondType(Bond::DOUBLE);
#ifdef VERBOSE_ENUMERATION
                std::cout << "Set bond to double" << std::endl;
#endif
              }
              if (bondtype == Bond::DOUBLE) {
                bond->setBondType(Bond::SINGLE);
#ifdef VERBOSE_ENUMERATION
                std::cout << "Set bond to single" << std::endl;
#endif
              }
            }
            markBondModified(bond->getIdx());
          }
          // TODO adjust charges
          if (!transform.Charges.empty()) {
            unsigned int ci = 0;
            for (const auto &pair : match) {
              Atom *atom = product->getAtomWithIdx(pair.second);
              atom->setFormalCharge(atom->getFormalCharge() +
                                    transform.Charges[ci++]);
            }
          }
#ifdef VERBOSE_ENUMERATION
          {
            SmilesWriteParams smilesWriteParams;
            smilesWriteParams.allBondsExplicit = true;
            std::cout << "pre-sanitize: "
                      << MolToSmiles(*product, smilesWriteParams) << std::endl;
          }
#endif

          try {
            // We only change bond orders/H counts/charges; the molecular graph
            // (and therefore ring topology) is unchanged.
            // `sanitizeMol()` always calls `clearComputedProps()` which resets
            // ring info and forces ring-finding for each generated tautomer.
            // Avoid that by clearing computed props without touching rings,
            // then running the specific sanitize steps we need.
            product->clearComputedProps(false);
            product->updatePropertyCache(false);

            // Note: We already checked for duplicates using the perturbed key
            // BEFORE copying the molecule. No need to check again here.

            MolOps::Kekulize(*product);
            MolOps::setAromaticity(*product);
            MolOps::setConjugation(*product);
            MolOps::setHybridization(*product);
            MolOps::adjustHs(*product);
          } catch (const KekulizeException &) {
            continue;
          }
#ifdef VERBOSE_ENUMERATION
          SmilesWriteParams smilesWriteParams;
          smilesWriteParams.allBondsExplicit = true;
          std::cout << "pre-setTautomerStereo: "
                    << MolToSmiles(*product, smilesWriteParams) << std::endl;
#endif
          setTautomerStereoAndIsoHs(mol, *product, res);

          // Quick duplicate check using cheap state key (also the map key)
          std::string stateKey = getTautomerStateKey(*product, canonOrder);
          bool stateExists = res.d_tautomers.find(stateKey) != res.d_tautomers.end();
          if (stateExists) {
#ifdef VERBOSE_ENUMERATION
            std::cout << "Previous tautomer state seen again (cheap check)"
                      << std::endl;
#endif
            continue;
          }

#ifdef VERBOSE_ENUMERATION
          (transform.Mol)->getProp(common_properties::_Name, name);
          std::cout << "Applied rule: " << name << " to "
                    << smilesTautomerPair.first << std::endl;
#endif
          // in addition to the above transformations, sanitization may modify
          // bonds, e.g. Cc1nc2ccccc2[nH]1
          // Use parallel iteration to avoid O(n) getBondWithIdx lookups
          {
            auto productBondIt = product->bonds().begin();
            for (const auto molBond : mol.bonds()) {
              const auto productBond = *productBondIt;
              ++productBondIt;
              auto i = molBond->getIdx();
              if (molBond->getBondType() != productBond->getBondType() &&
                  !res.d_modifiedBonds.test(i)) {
#ifdef VERBOSE_ENUMERATION
                std::cout << "Sanitization has modified bond " << i
                          << std::endl;
#endif
                markBondModified(static_cast<unsigned int>(i));
              }
            }
          }
          // Kekulized form will be created lazily when needed
#ifdef VERBOSE_ENUMERATION
          std::cout << "New tautomer added with state key: " << stateKey
                    << std::endl;
#endif
          // BOOST_LOG(rdInfoLog)
          //     << "Tautomer transform "
          //     <<
          //     transform.Mol->getProp<std::string>(common_properties::_Name)
          //     << " produced tautomer " << tsmiles << std::endl;
          res.d_tautomers[stateKey] =
              Tautomer(std::move(product), numModifiedAtoms, numModifiedBonds);
        }
      }
      smilesTautomerPair.second.d_done = true;
    }
    completed = true;
    size_t maxNumModifiedAtoms = numModifiedAtoms;
    size_t maxNumModifiedBonds = numModifiedBonds;
    for (auto it = res.d_tautomers.begin(); it != res.d_tautomers.end();) {
      auto &taut = it->second;
      if (!taut.d_done) {
        completed = false;
      }
      if ((taut.d_numModifiedAtoms < maxNumModifiedAtoms ||
           taut.d_numModifiedBonds < maxNumModifiedBonds) &&
          setTautomerStereoAndIsoHs(mol, *taut.tautomer, res)) {
        Tautomer tautStored = std::move(taut);
        it = res.d_tautomers.erase(it);
        tautStored.d_numModifiedAtoms = maxNumModifiedAtoms;
        tautStored.d_numModifiedBonds = maxNumModifiedBonds;
        auto insertRes = res.d_tautomers.insert(std::make_pair(
            getTautomerStateKey(*tautStored.tautomer, canonOrder),
            std::move(tautStored)));
        if (insertRes.second) {
          it = insertRes.first;
        }
      } else {
        ++it;
      }
    }
    if (bailOut && res.d_tautomers.size() < d_maxTautomers &&
        res.d_status == TautomerEnumeratorStatus::MaxTautomersReached) {
      res.d_status = TautomerEnumeratorStatus::Completed;
      bailOut = false;
    }
  }  // while
  res.fillTautomersItVec();
  if (!completed) {
    BOOST_LOG(rdWarningLog)
        << "Tautomer enumeration stopped at " << res.d_tautomers.size()
        << " tautomers: " << statusMsg.at(static_cast<size_t>(res.d_status))
        << std::endl;
  }

  return res;
}

TautomerEnumeratorResult TautomerEnumeratorResult::collapsedToSmilesKeys()
    const {
  TautomerEnumeratorResult out;
  out.d_status = d_status;
  out.d_modifiedAtoms = d_modifiedAtoms;
  out.d_modifiedBonds = d_modifiedBonds;

  for (const auto &kv : d_tautomers) {
    if (!kv.second.tautomer) {
      continue;
    }
    std::string smi = MolToSmiles(*kv.second.tautomer, true);
    // Deduplicate by SMILES. When multiple state-key entries share the same
    // canonical SMILES (e.g. because tautomer transforms reversed the bond
    // directions within a tautomeric path), prefer the one that preserves
    // more stereo information from the original molecule.
    auto it = out.d_tautomers.find(smi);
    if (it == out.d_tautomers.end()) {
      out.d_tautomers.emplace(std::move(smi), kv.second);
    } else {
      // Replace if the new tautomer has more stereo info
      unsigned int existingStereo =
          countBondStereo(*it->second.tautomer) +
          countAtomStereo(*it->second.tautomer);
      unsigned int newStereo =
          countBondStereo(*kv.second.tautomer) +
          countAtomStereo(*kv.second.tautomer);
      if (newStereo > existingStereo) {
        it->second = kv.second;
      }
    }
  }
  out.fillTautomersItVec();
  return out;
}

// pickCanonical non-templated overload that avoids recomputing SMILES
ROMol *TautomerEnumerator::pickCanonical(
    const TautomerEnumeratorResult &tautRes,
    boost::function<int(const ROMol &mol)> scoreFunc) const {
  ROMOL_SPTR bestMol;
  if (tautRes.d_tautomers.size() == 1) {
    bestMol = tautRes.d_tautomers.begin()->second.tautomer;
  } else {
    // Calculate score for each tautomer
    int bestScore = std::numeric_limits<int>::min();
    std::string bestSmiles;
    bool bestSmilesInitialized = false;
    for (const auto &t : tautRes.d_tautomers) {
      auto score = scoreFunc(*t.second.tautomer);
#ifdef VERBOSE_ENUMERATION
      std::cerr << "  " << t.first << " " << score << std::endl;
#endif
      if (score > bestScore) {
        bestScore = score;
        bestMol = t.second.tautomer;
        bestSmilesInitialized = false;
      } else if (score == bestScore) {
        // Tie-break by canonical SMILES (computed lazily only on ties).
        if (!bestSmilesInitialized && bestMol) {
          bestSmiles = MolToSmiles(*bestMol, true);
          bestSmilesInitialized = true;
        }
        auto curSmiles = MolToSmiles(*t.second.tautomer, true);
        if (!bestSmilesInitialized || curSmiles < bestSmiles) {
          bestSmiles = std::move(curSmiles);
          bestSmilesInitialized = true;
          bestMol = t.second.tautomer;
        } else if (curSmiles == bestSmiles) {
          // Same score and same SMILES — these are chemically equivalent
          // tautomers that differ only in which bonds along the tautomeric
          // path are single vs double (state-key dedup can produce these).
          // Prefer the one that preserves more stereo information, since
          // stereo on single bonds is lost during setTautomerStereoAndIsoHs.
          unsigned int curStereo = countBondStereo(*t.second.tautomer) +
                                   countAtomStereo(*t.second.tautomer);
          unsigned int bestStereo = countBondStereo(*bestMol) +
                                    countAtomStereo(*bestMol);
          if (curStereo > bestStereo) {
            bestMol = t.second.tautomer;
          }
        }
      }
    }
  }
  ROMol *res = new ROMol(*bestMol);
  static const bool cleanIt = true;
  static const bool force = true;
  MolOps::assignStereochemistry(*res, cleanIt, force);

  return res;
}

// static
NormalizedMolOrder TautomerEnumerator::normalizeAtomBondOrder(
    const ROMol &mol) {
  NormalizedMolOrder result;

  // Compute canonical atom output order via MolToSmiles. This resolves
  // symmetric-atom tie-breaking deterministically via canonical DFS.
  // The SMILES string itself is discarded — we only need the ordering.
  SmilesWriteParams smiPs;
  MolToSmiles(mol, smiPs);
  std::vector<unsigned int> smilesOrder;
  mol.getProp(common_properties::_smilesAtomOutputOrder, smilesOrder);
  CHECK_INVARIANT(smilesOrder.size() == mol.getNumAtoms(),
                  "MolToSmiles _smilesAtomOutputOrder size mismatch");

  // renumberAtoms(mol, order): newMol.atom[i] = mol.atom[order[i]]
  std::unique_ptr<ROMol> renumbered(MolOps::renumberAtoms(mol, smilesOrder));

  // Rebuild the molecule with bonds in sorted order so that neighbor
  // traversal (which depends on bond storage order) is deterministic.
  RWMol sortedMol;
  for (unsigned int i = 0; i < renumbered->getNumAtoms(); ++i) {
    sortedMol.addAtom(new Atom(*renumbered->getAtomWithIdx(i)), true, true);
  }
  // Sort bonds by (min_endpoint, max_endpoint) to determine addition
  // order. Preserve original begin/end atom indices on each bond to
  // maintain double-bond stereo semantics.
  std::vector<std::pair<std::pair<unsigned int, unsigned int>, unsigned int>>
      bondKeys;
  bondKeys.reserve(renumbered->getNumBonds());
  for (const auto bond : renumbered->bonds()) {
    unsigned int b = bond->getBeginAtomIdx();
    unsigned int e = bond->getEndAtomIdx();
    unsigned int lo = std::min(b, e);
    unsigned int hi = std::max(b, e);
    bondKeys.push_back({{lo, hi}, bond->getIdx()});
  }
  std::sort(bondKeys.begin(), bondKeys.end());
  for (const auto &bk : bondKeys) {
    const Bond *origBond = renumbered->getBondWithIdx(bk.second);
    Bond *newBond = new Bond(*origBond);
    // Keep original begin/end — do NOT normalize to (min, max)
    sortedMol.addBond(newBond, true);
  }

  // StereoGroups are not molecule properties and are not copied by default.
  // Atom indices are preserved by the rebuild, but bond indices are likely NOT
  // preserved because we add bonds in sorted order.
  copyStereoGroupsByIndices(*renumbered, sortedMol, false);
  // Note: chirality tags (CW/CCW) are NOT adjusted here even though
  // bond sorting changes neighbor iteration order. The enumeration
  // does not depend on chirality, and denormalization restores the
  // original bond storage order so chirality tags retain their meaning.

  // Build normPerm: smilesOrder[newIdx] = origIdx, so
  // normPerm[origIdx] = newIdx
  result.normPerm.resize(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    result.normPerm[smilesOrder[i]] = i;
  }
  result.mol.reset(new ROMol(sortedMol));
  return result;
}

// static
ROMol *TautomerEnumerator::denormalizeAtomBondOrder(
    const ROMol &canonical, const ROMol &origMol,
    const std::vector<unsigned int> &normPerm) {
  // normPerm[origIdx] = normalizedIdx → use as the permutation directly:
  // renumberAtoms(mol, perm) puts old-atom perm[i] at new position i.
  std::unique_ptr<ROMol> reordered{
      MolOps::renumberAtoms(canonical, normPerm)};

  RWMol fixedMol;
  for (unsigned int i = 0; i < reordered->getNumAtoms(); ++i) {
    fixedMol.addAtom(new Atom(*reordered->getAtomWithIdx(i)), true, true);
  }
  for (const auto origBond : origMol.bonds()) {
    unsigned int b = origBond->getBeginAtomIdx();
    unsigned int e = origBond->getEndAtomIdx();
    const Bond *canBond = reordered->getBondBetweenAtoms(b, e);
    if (canBond) {
      Bond *newBond = new Bond(*canBond);
      // Use original begin/end to preserve neighbor iteration order
      newBond->setBeginAtomIdx(b);
      newBond->setEndAtomIdx(e);
      fixedMol.addBond(newBond, true);
    }
  }

  // Here we add bonds in origMol order, so bond indices match and we can
  // (usually) do O(1) bond lookup by index.
  copyStereoGroupsByIndices(origMol, fixedMol, true);
  auto *res = new ROMol(fixedMol);

  // Recompute CIP codes: bond storage order changed, which affects
  // neighbor iteration order around chiral atoms.
  static const bool cleanIt = true;
  static const bool force = true;
  MolOps::assignStereochemistry(*res, cleanIt, force);

  return res;
}

TautomerEnumeratorResult
TautomerEnumeratorResult::denormalizedToOriginalOrder(
    const ROMol &origMol,
    const std::vector<unsigned int> &normPerm) const {
  TautomerEnumeratorResult out;
  out.d_status = d_status;
  out.d_modifiedAtoms = d_modifiedAtoms;
  out.d_modifiedBonds = d_modifiedBonds;

  for (const auto &kv : d_tautomers) {
    if (!kv.second.tautomer) {
      continue;
    }
    std::unique_ptr<ROMol> denorm{
        TautomerEnumerator::denormalizeAtomBondOrder(
            *kv.second.tautomer, origMol, normPerm)};
    ROMOL_SPTR denormSptr(denorm.release());
    // Re-key by the denormalized SMILES since the state-key is
    // meaningless after atom reordering.
    std::string smi = MolToSmiles(*denormSptr, true);
    auto it = out.d_tautomers.find(smi);
    if (it == out.d_tautomers.end()) {
      out.d_tautomers.emplace(std::move(smi), Tautomer(denormSptr, 0, 0));
    }
  }
  out.fillTautomersItVec();
  return out;
}

ROMol *TautomerEnumerator::canonicalize(
    const ROMol &mol, boost::function<int(const ROMol &mol)> scoreFunc) const {
  auto thisCopy = TautomerEnumerator(*this);
  thisCopy.setReassignStereo(false);

  auto normOrder = TautomerEnumerator::normalizeAtomBondOrder(mol);
  const ROMol &molForEnumeration = *normOrder.mol;

  auto res = thisCopy.enumerate(molForEnumeration);
  if (res.empty()) {
    BOOST_LOG(rdWarningLog)
        << "no tautomers found, returning input molecule" << std::endl;
    return new ROMol(mol);
  }
  // When no custom scorer provided, use optimized scoring that pre-filters
  // SubstructTerm patterns once for the input molecule rather than evaluating
  // all 12 substructure matches per tautomer. This is safe because
  // tautomerization only moves H and changes bond orders, never creates or
  // destroys heavy-atom bonds.
  if (!scoreFunc) {
    scoreFunc = TautomerScoringFunctions::makeOptimizedScorer(molForEnumeration);
  }
  std::unique_ptr<ROMol> canonical{pickCanonical(res, scoreFunc)};

  // Remap canonical tautomer back to original atom order.
  ROMol *result =
      TautomerEnumerator::denormalizeAtomBondOrder(*canonical, mol, normOrder.normPerm);

  // quickCopy during enumeration doesn't copy molecule properties or
  // conformers.  Restore both from the original molecule so that
  // downstream code (e.g. InChI generation) that relies on 2D/3D
  // coordinates or mol-level properties works correctly.
  result->updateProps(mol);
  for (auto confIt = mol.beginConformers(); confIt != mol.endConformers();
       ++confIt) {
    result->addConformer(new Conformer(**confIt), true);
  }
  return result;
}

void TautomerEnumerator::canonicalizeInPlace(
    RWMol &mol, boost::function<int(const ROMol &mol)> scoreFunc) const {
  auto thisCopy = TautomerEnumerator(*this);
  thisCopy.setReassignStereo(false);

  auto normOrder = TautomerEnumerator::normalizeAtomBondOrder(mol);
  const ROMol &molForEnumeration = *normOrder.mol;

  auto res = thisCopy.enumerate(molForEnumeration);
  if (res.empty()) {
    BOOST_LOG(rdWarningLog)
        << "no tautomers found, molecule unchanged" << std::endl;
    return;
  }
  // When no custom scorer provided, use optimized scoring that pre-filters
  // SubstructTerm patterns once for the input molecule.
  if (!scoreFunc) {
    scoreFunc = TautomerScoringFunctions::makeOptimizedScorer(molForEnumeration);
  }
  std::unique_ptr<ROMol> canonical{pickCanonical(res, scoreFunc)};

  // Remap canonical tautomer back to original atom order.
  std::unique_ptr<ROMol> tmp{
      TautomerEnumerator::denormalizeAtomBondOrder(*canonical, mol, normOrder.normPerm)};

  TEST_ASSERT(tmp->getNumAtoms() == mol.getNumAtoms());
  TEST_ASSERT(tmp->getNumBonds() == mol.getNumBonds());
  // now copy the info from the canonical tautomer over to the input molecule
  // iterate both molecules' atoms and bonds in parallel - they have matching indices
  {
    auto molAtomIt = mol.atoms().begin();
    for (const auto tmpAtom : tmp->atoms()) {
      auto atom = *molAtomIt;
      ++molAtomIt;
      TEST_ASSERT(tmpAtom->getAtomicNum() == atom->getAtomicNum());
      atom->setFormalCharge(tmpAtom->getFormalCharge());
      atom->setNoImplicit(tmpAtom->getNoImplicit());
      atom->setIsAromatic(tmpAtom->getIsAromatic());
      atom->setNumExplicitHs(tmpAtom->getNumExplicitHs());
      atom->setNumRadicalElectrons(tmpAtom->getNumRadicalElectrons());
      atom->setChiralTag(tmpAtom->getChiralTag());
    }
  }
  {
    auto molBondIt = mol.bonds().begin();
    for (const auto tmpBond : tmp->bonds()) {
      auto bond = *molBondIt;
      ++molBondIt;
      TEST_ASSERT(tmpBond->getBeginAtomIdx() == bond->getBeginAtomIdx());
      TEST_ASSERT(tmpBond->getEndAtomIdx() == bond->getEndAtomIdx());
      bond->setBondType(tmpBond->getBondType());
      bond->setBondDir(tmpBond->getBondDir());
      bond->setIsAromatic(tmpBond->getIsAromatic());
      bond->setIsConjugated(tmpBond->getIsConjugated());
      if (tmpBond->getStereoAtoms().size() == 2) {
        bond->setStereoAtoms(tmpBond->getStereoAtoms()[0],
                             tmpBond->getStereoAtoms()[1]);
      }
      bond->setStereo(tmpBond->getStereo());
    }
  }
  mol.updatePropertyCache(false);
}

}  // namespace MolStandardize
}  // namespace RDKit
