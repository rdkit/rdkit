//
//  Copyright (C) 2011-2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/types.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <vector>
#include <string>

namespace {
class ss_matcher {
 public:
  ss_matcher(const std::string &pattern) : m_pattern(pattern) {
    m_needCopies = (pattern.find_first_of("$") != std::string::npos);
    RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    m_matcher = p;
    POSTCONDITION(m_matcher, "no matcher");
  };
  const RDKit::ROMol *getMatcher() const { return m_matcher; };
  unsigned int countMatches(const RDKit::ROMol &mol) const {
    PRECONDITION(m_matcher, "no matcher");
    std::vector<RDKit::MatchVectType> matches;
    // This is an ugly one. Recursive queries aren't thread safe.
    // Unfortunately we have to take a performance hit here in order
    // to guarantee thread safety
    if (m_needCopies) {
      const RDKit::ROMol nm(*(m_matcher), true);
      RDKit::SubstructMatch(mol, nm, matches);
    } else {
      const RDKit::ROMol &nm = *m_matcher;
      RDKit::SubstructMatch(mol, nm, matches);
    }
    return matches.size();
  }
  ~ss_matcher() { delete m_matcher; };

 private:
  ss_matcher() : m_pattern(""){};
  std::string m_pattern;
  bool m_needCopies{false};
  const RDKit::ROMol *m_matcher{nullptr};
};
}  // namespace

typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
                         boost::flyweights::no_tracking>
    pattern_flyweight;
#define SMARTSCOUNTFUNC(nm, pattern, vers)         \
  const std::string nm##Version = vers;            \
  unsigned int calc##nm(const RDKit::ROMol &mol) { \
    pattern_flyweight m(pattern);                  \
    return m.get().countMatches(mol);              \
  }                                                \
  extern int no_such_variable

namespace RDKit {
namespace Descriptors {

unsigned int calcLipinskiHBA(const ROMol &mol) {
  unsigned int res = 0;
  for (ROMol::ConstAtomIterator iter = mol.beginAtoms(); iter != mol.endAtoms();
       ++iter) {
    if ((*iter)->getAtomicNum() == 7 || (*iter)->getAtomicNum() == 8) {
      ++res;
    }
  }
  return res;
}

unsigned int calcLipinskiHBD(const ROMol &mol) {
  unsigned int res = 0;
  for (ROMol::ConstAtomIterator iter = mol.beginAtoms(); iter != mol.endAtoms();
       ++iter) {
    if (((*iter)->getAtomicNum() == 7 || (*iter)->getAtomicNum() == 8)) {
      res += (*iter)->getTotalNumHs(true);
    }
  }
  return res;
}

namespace {
#ifdef RDK_USE_STRICT_ROTOR_DEFINITION
const NumRotatableBondsOptions DefaultStrictDefinition = Strict;
#else
const NumRotatableBondsOptions DefaultStrictDefinition = NonStrict;
#endif
}  // namespace

const std::string NumRotatableBondsVersion = "3.1.0";
unsigned int calcNumRotatableBonds(const ROMol &mol,
                                   NumRotatableBondsOptions strict) {
  if (strict == Default) {
    strict = DefaultStrictDefinition;
  }

  if (strict == NonStrict) {
    std::string pattern = "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]";
    pattern_flyweight m(pattern);
    return m.get().countMatches(mol);
  } else if (strict == Strict) {
    std::string strict_pattern =
        "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])("
        "[CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]="
        "[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$"
        "(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])(["
        "CH3])[CH3])]";
    pattern_flyweight m(strict_pattern);
    return m.get().countMatches(mol);
  } else {
    // Major changes in definition relative to the original GPS calculator:
    //   Bonds linking ring systems:
    //     - Single bonds between aliphatic ring Cs are always rotatable. This
    //     means that the
    //       central bond in CC1CCCC(C)C1-C1C(C)CCCC1C is now considered
    //       rotatable; it was not
    //       before.
    //     - Heteroatoms in the linked rings no longer affect whether or not the
    //     linking bond
    //       is rotatable
    //     - the linking bond in systems like Cc1cccc(C)c1-c1c(C)cccc1 is now
    //     considered
    //       non-rotatable
    pattern_flyweight rotBonds_matcher("[!$([D1&!#1])]-,:;!@[!$([D1&!#1])]");
    pattern_flyweight nonRingAmides_matcher("[C&!R](=O)NC");
    pattern_flyweight symRings_matcher(
        "[a;r6;$(a(-,:;!@[a;r6])(a[!#1])a[!#1])]-,:;!@[a;r6;$(a(-,:;!@[a;r6])("
        "a[!#1])a)]");
    pattern_flyweight terminalTripleBonds_matcher("C#[#6,#7]");

    std::vector<MatchVectType> matches;

    // initialize to the number of bonds matching the base pattern:
    int res = rotBonds_matcher.get().countMatches(mol);
    if (!res) {
      return 0;
    }

    // remove symmetrical rings:
    res -= symRings_matcher.get().countMatches(mol);
    if (res < 0) {
      res = 0;
    }

    // remove triple bonds
    res -= terminalTripleBonds_matcher.get().countMatches(mol);
    if (res < 0) {
      res = 0;
    }

    // removing amides is more complex
    boost::dynamic_bitset<> atomsSeen(mol.getNumAtoms());
    SubstructMatch(mol, *(nonRingAmides_matcher.get().getMatcher()), matches);
    for (const auto &iv : matches) {
      bool distinct = true;
      for (const auto &mIt : iv) {
        if (atomsSeen[mIt.second]) {
          distinct = false;
        }
        atomsSeen.set(mIt.second);
      }
      if (distinct && res > 0) {
        --res;
      }
    }

    if (res < 0) {
      res = 0;
    }
    return static_cast<unsigned int>(res);
  }
}

unsigned int calcNumRotatableBonds(const ROMol &mol, bool strict) {
  return calcNumRotatableBonds(mol, (strict) ? Strict : NonStrict);
}

// SMARTSCOUNTFUNC(NumHBD,
// "[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]","2.0.1" ) ;
SMARTSCOUNTFUNC(NumHBD, "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]",
                "2.0.1");
SMARTSCOUNTFUNC(NumHBA,
                "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$("
                "[N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]",
                "2.0.1");
SMARTSCOUNTFUNC(NumHeteroatoms, "[!#6;!#1]", "1.0.1");
SMARTSCOUNTFUNC(NumAmideBonds, "C(=[O;!R])N", "1.0.0");

const std::string NumRingsVersion = "1.0.1";
unsigned int calcNumRings(const ROMol &mol) {
  return mol.getRingInfo()->numRings();
}

const std::string FractionCSP3Version = "1.0.0";
double calcFractionCSP3(const ROMol &mol) {
  unsigned int nCSP3 = 0;
  unsigned int nC = 0;
  ROMol::VERTEX_ITER atBegin, atEnd;
  boost::tie(atBegin, atEnd) = mol.getVertices();
  while (atBegin != atEnd) {
    const Atom *at = mol[*atBegin];
    if (at->getAtomicNum() == 6) {
      ++nC;
      if (at->getTotalDegree() == 4) {
        ++nCSP3;
      }
    }
    ++atBegin;
  }
  if (!nC) {
    return 0;
  }
  return static_cast<double>(nCSP3) / nC;
}

const std::string NumHeterocyclesVersion = "1.0.0";
unsigned int calcNumHeterocycles(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto &iv : mol.getRingInfo()->atomRings()) {
    for (auto i : iv) {
      if (mol.getAtomWithIdx(i)->getAtomicNum() != 6) {
        ++res;
        break;
      }
    }
  }
  return res;
}
const std::string NumAromaticRingsVersion = "1.0.0";
unsigned int calcNumAromaticRings(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto &iv : mol.getRingInfo()->bondRings()) {
    ++res;
    for (auto i : iv) {
      if (!mol.getBondWithIdx(i)->getIsAromatic()) {
        --res;
        break;
      }
    }
  }
  return res;
}
const std::string NumSaturatedRingsVersion = "1.0.0";
unsigned int calcNumSaturatedRings(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto &iv : mol.getRingInfo()->bondRings()) {
    ++res;
    for (int i : iv) {
      if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
          mol.getBondWithIdx(i)->getIsAromatic()) {
        --res;
        break;
      }
    }
  }
  return res;
}
const std::string NumAliphaticRingsVersion = "1.0.0";
unsigned int calcNumAliphaticRings(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto &iv : mol.getRingInfo()->bondRings()) {
    for (auto i : iv) {
      if (!mol.getBondWithIdx(i)->getIsAromatic()) {
        ++res;
        break;
      }
    }
  }
  return res;
}
const std::string NumAromaticHeterocyclesVersion = "1.0.0";
unsigned int calcNumAromaticHeterocycles(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto &iv : mol.getRingInfo()->bondRings()) {
    bool countIt = false;
    for (auto i : iv) {
      if (!mol.getBondWithIdx(i)->getIsAromatic()) {
        countIt = false;
        break;
      }
      // we're checking each atom twice, which is kind of doofy, but this
      // function is hopefully not going to be a big time sink.
      if (!countIt &&
          (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
        countIt = true;
      }
    }
    if (countIt) {
      ++res;
    }
  }
  return res;
}
const std::string NumAromaticCarbocyclesVersion = "1.0.0";
unsigned int calcNumAromaticCarbocycles(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto &iv : mol.getRingInfo()->bondRings()) {
    bool countIt = true;
    for (auto i : iv) {
      if (!mol.getBondWithIdx(i)->getIsAromatic()) {
        countIt = false;
        break;
      }
      // we're checking each atom twice, which is kind of doofy, but this
      // function is hopefully not going to be a big time sync.
      if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
          mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
        countIt = false;
        break;
      }
    }
    if (countIt) {
      ++res;
    }
  }
  return res;
}
const std::string NumAliphaticHeterocyclesVersion = "1.0.0";
unsigned int calcNumAliphaticHeterocycles(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto &iv : mol.getRingInfo()->bondRings()) {
    bool hasAliph = false;
    bool hasHetero = false;
    for (auto i : iv) {
      if (!mol.getBondWithIdx(i)->getIsAromatic()) {
        hasAliph = true;
      }
      // we're checking each atom twice, which is kind of doofy, but this
      // function is hopefully not going to be a big time sink.
      if (!hasHetero &&
          (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
        hasHetero = true;
      }
    }
    if (hasHetero && hasAliph) {
      ++res;
    }
  }
  return res;
}
const std::string NumAliphaticCarbocyclesVersion = "1.0.0";
unsigned int calcNumAliphaticCarbocycles(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto &iv : mol.getRingInfo()->bondRings()) {
    bool hasAliph = false;
    bool hasHetero = false;
    for (auto i : iv) {
      if (!mol.getBondWithIdx(i)->getIsAromatic()) {
        hasAliph = true;
      }
      // we're checking each atom twice, which is kind of doofy, but this
      // function is hopefully not going to be a big time sync.
      if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
          mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
        hasHetero = true;
        break;
      }
    }
    if (hasAliph && !hasHetero) {
      ++res;
    }
  }
  return res;
}
const std::string NumSaturatedHeterocyclesVersion = "1.0.0";
unsigned int calcNumSaturatedHeterocycles(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto &iv : mol.getRingInfo()->bondRings()) {
    bool countIt = false;
    for (auto i : iv) {
      if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
          mol.getBondWithIdx(i)->getIsAromatic()) {
        countIt = false;
        break;
      }
      // we're checking each atom twice, which is kind of doofy, but this
      // function is hopefully not going to be a big time sync.
      if (!countIt &&
          (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
        countIt = true;
      }
    }
    if (countIt) {
      ++res;
    }
  }
  return res;
}
const std::string NumSaturatedCarbocyclesVersion = "1.0.0";
unsigned int calcNumSaturatedCarbocycles(const ROMol &mol) {
  unsigned int res = 0;
  for (const auto &iv : mol.getRingInfo()->bondRings()) {
    bool countIt = true;
    for (auto i : iv) {
      if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
          mol.getBondWithIdx(i)->getIsAromatic()) {
        countIt = false;
        break;
      }
      // we're checking each atom twice, which is kind of doofy, but this
      // function is hopefully not going to be a big time sync.
      if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
          mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
        countIt = false;
        break;
      }
    }
    if (countIt) {
      ++res;
    }
  }
  return res;
}

const std::string NumSpiroAtomsVersion = "1.0.0";
unsigned int calcNumSpiroAtoms(const ROMol &mol,
                               std::vector<unsigned int> *atoms) {
  if (!mol.getRingInfo() || !mol.getRingInfo()->isSssrOrBetter()) {
    MolOps::findSSSR(mol);
  }
  const RingInfo *rInfo = mol.getRingInfo();
  std::vector<unsigned int> lAtoms;
  if (!atoms) {
    atoms = &lAtoms;
  }

  for (unsigned int i = 0; i < rInfo->atomRings().size(); ++i) {
    const INT_VECT &ri = rInfo->atomRings()[i];
    for (unsigned int j = i + 1; j < rInfo->atomRings().size(); ++j) {
      const INT_VECT &rj = rInfo->atomRings()[j];
      // EFF: using intersect here does more work and memory allocation than is
      // required
      INT_VECT inter;
      Intersect(ri, rj, inter);
      if (inter.size() == 1) {
        if (std::find(atoms->begin(), atoms->end(), inter[0]) == atoms->end()) {
          atoms->push_back(inter[0]);
        }
      }
    }
  }
  return atoms->size();
}

const std::string NumBridgeheadAtomsVersion = "2.0.0";
unsigned int calcNumBridgeheadAtoms(const ROMol &mol,
                                    std::vector<unsigned int> *atoms) {
  if (!mol.getRingInfo() || !mol.getRingInfo()->isSssrOrBetter()) {
    MolOps::findSSSR(mol);
  }
  const RingInfo *rInfo = mol.getRingInfo();
  std::vector<unsigned int> lAtoms;
  if (!atoms) {
    atoms = &lAtoms;
  }

  for (unsigned int i = 0; i < rInfo->bondRings().size(); ++i) {
    const INT_VECT &ri = rInfo->bondRings()[i];
    for (unsigned int j = i + 1; j < rInfo->bondRings().size(); ++j) {
      const INT_VECT &rj = rInfo->bondRings()[j];
      // EFF: using intersect here does more work and memory allocation than is
      // required
      INT_VECT inter;
      Intersect(ri, rj, inter);
      if (inter.size() > 1) {
        INT_VECT atomCounts(mol.getNumAtoms(), 0);
        for (auto ii : inter) {
          atomCounts[mol.getBondWithIdx(ii)->getBeginAtomIdx()] += 1;
          atomCounts[mol.getBondWithIdx(ii)->getEndAtomIdx()] += 1;
        }
        for (unsigned int ti = 0; ti < atomCounts.size(); ++ti) {
          if (atomCounts[ti] == 1) {
            if (std::find(atoms->begin(), atoms->end(), ti) == atoms->end()) {
              atoms->push_back(ti);
            }
          }
        }
      }
    }
  }
  return atoms->size();
}

namespace {
bool hasStereoAssigned(const ROMol &mol) {
  return mol.hasProp(common_properties::_StereochemDone);
}
}  // namespace
const std::string NumAtomStereoCentersVersion = "1.0.1";
unsigned int numAtomStereoCenters(const ROMol &mol) {
  std::unique_ptr<ROMol> tmol;
  const ROMol *mptr = &mol;
  if (!hasStereoAssigned(mol)) {
    tmol.reset(new ROMol(mol));
    constexpr bool cleanIt = true;
    constexpr bool force = true;
    constexpr bool flagPossible = true;
    MolOps::assignStereochemistry(*tmol, cleanIt, force, flagPossible);
    mptr = tmol.get();
  }

  unsigned int res = 0;
  for (const auto &atom : mptr->atoms()) {
    if (atom->hasProp(common_properties::_ChiralityPossible)) {
      ++res;
    }
  }
  return res;
}

const std::string NumUnspecifiedAtomStereoCentersVersion = "1.0.1";
unsigned int numUnspecifiedAtomStereoCenters(const ROMol &mol) {
  std::unique_ptr<ROMol> tmol;
  const ROMol *mptr = &mol;
  if (!hasStereoAssigned(mol)) {
    tmol.reset(new ROMol(mol));
    constexpr bool cleanIt = true;
    constexpr bool force = true;
    constexpr bool flagPossible = true;
    MolOps::assignStereochemistry(*tmol, cleanIt, force, flagPossible);
    mptr = tmol.get();
  }

  unsigned int res = 0;
  for (const auto &atom : mptr->atoms()) {
    if (atom->hasProp(common_properties::_ChiralityPossible) &&
        atom->getChiralTag() == Atom::CHI_UNSPECIFIED) {
      ++res;
    }
  }
  return res;
}

}  // end of namespace Descriptors
}  // end of namespace RDKit
