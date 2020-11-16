//
//  Copyright (C) 2013-2018 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MolFragmenter.h"
#include "ChemTransforms.h"
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/StreamOps.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <boost/functional/hash.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <sstream>
#include <map>

namespace RDKit {
namespace MolFragmenter {
std::size_t hash_value(const FragmenterBondType &fbt) {
  size_t res = boost::hash<int>()((int)fbt.bondType);
  boost::hash_combine(res, fbt.atom1Label);
  boost::hash_combine(res, fbt.atom2Label);
  return res;
}
bool operator==(const FragmenterBondType &v1, const FragmenterBondType &v2) {
  return (v1.atom1Label == v2.atom1Label) && (v1.atom2Label == v2.atom2Label) &&
         (v1.bondType == v2.bondType);
}
void constructFragmenterAtomTypes(
    std::istream *inStream, std::map<unsigned int, std::string> &defs,
    const std::string &comment, bool validate,
    std::map<unsigned int, ROMOL_SPTR> *environs) {
  PRECONDITION(inStream, "no stream");
  defs.clear();
  unsigned int line = 0;
  while (!inStream->eof() && !inStream->fail()) {
    ++line;
    std::string tempStr = getLine(inStream);
    if (tempStr == "" || tempStr.find(comment) == 0) {
      continue;
    }
    std::vector<std::string> tokens;
    boost::split(tokens, tempStr, boost::is_any_of(" \t"),
                 boost::token_compress_on);
    if (tokens.size() < 2) {
      BOOST_LOG(rdWarningLog)
          << "line " << line << " is too short" << std::endl;
      continue;
    }
    auto idx = boost::lexical_cast<unsigned int>(tokens[0]);
    if (defs.find(idx) != defs.end()) {
      BOOST_LOG(rdWarningLog)
          << "definition #" << idx
          << " encountered more than once. Using the first occurrence."
          << std::endl;
      continue;
    }
    if (validate || environs) {
      ROMol *p = SmartsToMol(tokens[1]);
      if (!p) {
        BOOST_LOG(rdWarningLog) << "cannot convert SMARTS " << tokens[1]
                                << " to molecule at line " << line << std::endl;
        continue;
      }
      if (!environs) {
        delete p;
      } else {
        (*environs)[idx] = ROMOL_SPTR(p);
      }
    }
    defs[idx] = tokens[1];
  }
}
void constructFragmenterAtomTypes(
    const std::string &str, std::map<unsigned int, std::string> &defs,
    const std::string &comment, bool validate,
    std::map<unsigned int, ROMOL_SPTR> *environs) {
  std::stringstream istr(str);
  constructFragmenterAtomTypes(&istr, defs, comment, validate, environs);
}
void constructBRICSAtomTypes(std::map<unsigned int, std::string> &defs,
                             std::map<unsigned int, ROMOL_SPTR> *environs) {
  /*
     After some discussion, the L2 definitions ("N.pl3" in the original
     paper) have been removed and incorporated into a (almost) general
     purpose amine definition in L5 ("N.sp3" in the paper).

     The problem is one of consistency.
        Based on the original definitions you should get the following
        fragmentations:
          C1CCCCC1NC(=O)C -> C1CCCCC1N[2*].[1*]C(=O)C
          c1ccccc1NC(=O)C -> c1ccccc1[16*].[2*]N[2*].[1*]C(=O)C
        This difference just didn't make sense to us. By switching to
        the unified definition we end up with:
          C1CCCCC1NC(=O)C -> C1CCCCC1[15*].[5*]N[5*].[1*]C(=O)C
          c1ccccc1NC(=O)C -> c1ccccc1[16*].[5*]N[5*].[1*]C(=O)C
  */
  const std::string BRICSdefs =
      "1 [C;D3]([#0,#6,#7,#8])(=O)\n\
3 [O;D2]-;!@[#0,#6,#1]\n\
5 [N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]\n\
9 [n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]\n\
10 [N;R;$(N(@C(=O))@[C,N,O,S])]\n\
11 [S;D2](-;!@[#0,#6])\n\
12 [S;D4]([#6,#0])(=O)(=O)\n\
6 [C;D3;!R](=O)-;!@[#0,#6,#7,#8]\n\
13 [C;$(C(-;@[C,N,O,S])-;@[N,O,S])]\n\
14 [c;$(c(:[c,n,o,s]):[n,o,s])]\n\
15 [C;$(C(-;@C)-;@C)]\n\
4 [C;!D1;!$(C=*)]-;!@[#6]\n\
7 [C;D2,D3]-[#6]\n\
8 [C;!R;!D1;!$(C!-*)]\n\
16 [c;$(c(:c):c)]";
  constructFragmenterAtomTypes(BRICSdefs, defs, "//", true, environs);
}

void constructFragmenterBondTypes(
    std::istream *inStream,
    const std::map<unsigned int, std::string> &atomTypes,
    std::vector<FragmenterBondType> &defs, const std::string &comment,
    bool validate, bool labelByConnector) {
  PRECONDITION(inStream, "no stream");
  defs.clear();
  defs.resize(0);

  unsigned int line = 0;
  while (!inStream->eof() && !inStream->fail()) {
    ++line;
    std::string tempStr = getLine(inStream);
    if (tempStr == "" || tempStr.find(comment) == 0) {
      continue;
    }
    std::vector<std::string> tokens;
    boost::split(tokens, tempStr, boost::is_any_of(" \t"),
                 boost::token_compress_on);
    if (tokens.size() < 3) {
      BOOST_LOG(rdWarningLog)
          << "line " << line << " is too short" << std::endl;
      continue;
    }
    auto idx1 = boost::lexical_cast<unsigned int>(tokens[0]);
    if (atomTypes.find(idx1) == atomTypes.end()) {
      BOOST_LOG(rdWarningLog)
          << "atom type #" << idx1 << " not recognized." << std::endl;
      continue;
    }
    auto idx2 = boost::lexical_cast<unsigned int>(tokens[1]);
    if (atomTypes.find(idx2) == atomTypes.end()) {
      BOOST_LOG(rdWarningLog)
          << "atom type #" << idx2 << " not recognized." << std::endl;
      continue;
    }
    std::string sma1 = atomTypes.find(idx1)->second;
    std::string sma2 = atomTypes.find(idx2)->second;
    std::string smarts = "[$(" + sma1 + ")]" + tokens[2] + "[$(" + sma2 + ")]";
    ROMol *p = SmartsToMol(smarts);
    if (validate) {
      if (!p) {
        BOOST_LOG(rdWarningLog) << "cannot convert SMARTS " << smarts
                                << " to molecule at line " << line << std::endl;
        continue;
      }
    }
    FragmenterBondType fbt;
    fbt.atom1Type = idx1;
    fbt.atom2Type = idx2;
    if (labelByConnector) {
      fbt.atom1Label = idx1;
      fbt.atom2Label = idx2;
    } else {
      fbt.atom1Label = idx2;
      fbt.atom2Label = idx1;
    }
    if (p) {
      // for the purposes of replacing the bond, we'll use just the first
      // character to set the bond type (if we recognize it):
      switch (tokens[2][0]) {
        case '-':
          fbt.bondType = Bond::SINGLE;
          break;
        case '=':
          fbt.bondType = Bond::DOUBLE;
          break;
        case '#':
          fbt.bondType = Bond::TRIPLE;
          break;
        case ':':
          fbt.bondType = Bond::AROMATIC;
          break;
        default:
          fbt.bondType = p->getBondWithIdx(0)->getBondType();
      }
      fbt.query = ROMOL_SPTR(p);
    } else {
      fbt.bondType = Bond::UNSPECIFIED;
      fbt.query = ROMOL_SPTR();
    }
    defs.push_back(fbt);
  }
}

void constructFragmenterBondTypes(
    const std::string &str,
    const std::map<unsigned int, std::string> &atomTypes,
    std::vector<FragmenterBondType> &defs, const std::string &comment,
    bool validate, bool labelByConnector) {
  std::stringstream istr(str);
  constructFragmenterBondTypes(&istr, atomTypes, defs, comment, validate,
                               labelByConnector);
}
void constructBRICSBondTypes(std::vector<FragmenterBondType> &defs) {
  const std::string BRICSdefs =
      "// L1\n\
1 3 -;!@\n\
1 5 -;!@\n\
1 10 -;!@\n\
// L3 \n\
3 4 -;!@\n\
3 13 -;!@\n\
3 14 -;!@\n\
3 15 -;!@\n\
3 16 -;!@\n\
// L4\n\
4 5 -;!@\n\
4 11 -;!@\n\
// L5\n\
5 12 -;!@\n\
5 14 -;!@\n\
5 16 -;!@\n\
5 13 -;!@\n\
5 15 -;!@\n\
// L6\n\
6 13 -;!@\n\
6 14 -;!@\n\
6 15 -;!@\n\
6 16 -;!@\n\
// L7\n\
7 7 =;!@\n\
// L8\n\
8 9 -;!@\n\
8 10 -;!@\n\
8 13 -;!@\n\
8 14 -;!@\n\
8 15 -;!@\n\
8 16 -;!@\n\
// L9\n\
9 13 -;!@ // not in original paper\n\
9 14 -;!@ // not in original paper\n\
9 15 -;!@\n\
9 16 -;!@\n\
// L10\n\
10 13 -;!@\n\
10 14 -;!@\n\
10 15 -;!@\n\
10 16 -;!@\n\
// L11\n\
11 13 -;!@\n\
11 14 -;!@\n\
11 15 -;!@\n\
11 16 -;!@\n\
// L12\n\
// none left\n\
// L13\n\
13 14 -;!@\n\
13 15 -;!@\n\
13 16 -;!@\n\
// L14\n\
14 14 -;!@ // not in original paper\n\
14 15 -;!@\n\
14 16 -;!@\n\
// L15\n\
15 16 -;!@\n\
// L16\n\
16 16 -;!@ // not in original paper";
  std::map<unsigned int, std::string> atTypes;
  constructBRICSAtomTypes(atTypes);
  constructFragmenterBondTypes(BRICSdefs, atTypes, defs, "//", true, false);
}

namespace {
boost::uint64_t nextBitCombo(boost::uint64_t v) {
  // code from:
  // http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
  boost::uint64_t t = (v | (v - 1)) + 1;
  return t | ((((t & -t) / (v & -v)) >> 1) - 1);
}
}  // namespace

void fragmentOnSomeBonds(
    const ROMol &mol, const std::vector<unsigned int> &bondIndices,
    std::vector<ROMOL_SPTR> &resMols, unsigned int maxToCut, bool addDummies,
    const std::vector<std::pair<unsigned int, unsigned int>> *dummyLabels,
    const std::vector<Bond::BondType> *bondTypes,
    std::vector<std::vector<unsigned int>> *nCutsPerAtom) {
  PRECONDITION((!dummyLabels || dummyLabels->size() == bondIndices.size()),
               "bad dummyLabel vector");
  PRECONDITION((!bondTypes || bondTypes->size() == bondIndices.size()),
               "bad bondType vector");
  if (bondIndices.size() > 63) {
    throw ValueErrorException("currently can only fragment on up to 63 bonds");
  }
  if (!maxToCut || !mol.getNumAtoms() || !bondIndices.size()) {
    return;
  }

  boost::uint64_t state = (0x1L << maxToCut) - 1;
  boost::uint64_t stop = 0x1L << bondIndices.size();
  std::vector<unsigned int> fragmentHere(maxToCut);
  std::vector<std::pair<unsigned int, unsigned int>> *dummyLabelsHere = nullptr;
  if (dummyLabels) {
    dummyLabelsHere =
        new std::vector<std::pair<unsigned int, unsigned int>>(maxToCut);
  }
  std::vector<Bond::BondType> *bondTypesHere = nullptr;
  if (bondTypes) {
    bondTypesHere = new std::vector<Bond::BondType>(maxToCut);
  }
  while (state < stop) {
    unsigned int nSeen = 0;
    for (unsigned int i = 0; i < bondIndices.size() && nSeen < maxToCut; ++i) {
      if (state & (0x1L << i)) {
        fragmentHere[nSeen] = bondIndices[i];
        if (dummyLabelsHere) {
          (*dummyLabelsHere)[nSeen] = (*dummyLabels)[i];
        }
        if (bondTypesHere) {
          (*bondTypesHere)[nSeen] = (*bondTypes)[i];
        }
        ++nSeen;
      }
    }
    std::vector<unsigned int> *lCutsPerAtom = nullptr;
    if (nCutsPerAtom) {
      nCutsPerAtom->push_back(std::vector<unsigned int>(mol.getNumAtoms()));
      lCutsPerAtom = &(nCutsPerAtom->back());
    }
    ROMol *nm = fragmentOnBonds(mol, fragmentHere, addDummies, dummyLabelsHere,
                                bondTypesHere, lCutsPerAtom);
    resMols.emplace_back(nm);

    state = nextBitCombo(state);
  }
  delete dummyLabelsHere;
  delete bondTypesHere;
}

namespace {
void checkChiralityPostMove(const ROMol &mol, const Atom *oAt, Atom *nAt,
                            const Bond *bond) {
  static const std::string newBondOrder = "_newBondOrder";
  INT_LIST newOrder;
  INT_LIST incomingOrder;

  const int check_bond_index = static_cast<int>(bond->getIdx());
  // since we may call this function more than once, we need to keep track of
  // whether or not we've already been called and what the new atom order is.
  // we do this with a property.
  // this was github #1734
  if (nAt->getPropIfPresent(newBondOrder, incomingOrder)) {
    for (int bidx : incomingOrder) {
      if (bidx != check_bond_index) {
        newOrder.push_back(bidx);
      }
    }
  } else {
    ROMol::OEDGE_ITER beg, end;
    boost::tie(beg, end) = mol.getAtomBonds(oAt);
    while (beg != end) {
      const Bond *obond = mol[*beg];
      ++beg;
      if (obond == bond) {
        continue;
      }
      newOrder.push_back(obond->getIdx());
    }
  }
  newOrder.push_back(bond->getIdx());
  nAt->setProp(newBondOrder, newOrder, true);
  unsigned int nSwaps = oAt->getPerturbationOrder(newOrder);
  // std::copy(newOrder.begin(), newOrder.end(),
  //           std::ostream_iterator<int>(std::cerr, ", "));
  // std::cerr << std::endl;
  // std::cerr<<"ccpm: "<<oAt->getIdx()<<"->"<<nAt->getIdx()<<" bond:
  // "<<bond->getIdx()<<" swaps: "<<nSwaps<<std::endl;
  nAt->setChiralTag(oAt->getChiralTag());
  if (nSwaps % 2) {
    nAt->invertChirality();
  }
}
}  // namespace

ROMol *fragmentOnBonds(
    const ROMol &mol, const std::vector<unsigned int> &bondIndices,
    bool addDummies,
    const std::vector<std::pair<unsigned int, unsigned int>> *dummyLabels,
    const std::vector<Bond::BondType> *bondTypes,
    std::vector<unsigned int> *nCutsPerAtom) {
  PRECONDITION((!dummyLabels || dummyLabels->size() == bondIndices.size()),
               "bad dummyLabel vector");
  PRECONDITION((!bondTypes || bondTypes->size() == bondIndices.size()),
               "bad bondType vector");
  PRECONDITION((!nCutsPerAtom || nCutsPerAtom->size() == mol.getNumAtoms()),
               "bad nCutsPerAtom vector");
  if (nCutsPerAtom) {
    BOOST_FOREACH (unsigned int &nCuts, *nCutsPerAtom) { nCuts = 0; }
  }
  auto *res = new RWMol(mol);
  if (!mol.getNumAtoms()) {
    return res;
  }

  std::vector<Bond *> bondsToRemove;
  bondsToRemove.reserve(bondIndices.size());
  BOOST_FOREACH (unsigned int bondIdx, bondIndices) {
    bondsToRemove.push_back(res->getBondWithIdx(bondIdx));
  }
  for (unsigned int i = 0; i < bondsToRemove.size(); ++i) {
    const Bond *bond = bondsToRemove[i];
    unsigned int bidx = bond->getBeginAtomIdx();
    unsigned int eidx = bond->getEndAtomIdx();
    Bond::BondType bT = bond->getBondType();
    Bond::BondDir bD = bond->getBondDir();
    unsigned int bondidx;
    res->removeBond(bidx, eidx);
    if (nCutsPerAtom) {
      (*nCutsPerAtom)[bidx] += 1;
      (*nCutsPerAtom)[eidx] += 1;
    }
    if (addDummies) {
      Atom *at1, *at2;
      at1 = new Atom(0);
      at2 = new Atom(0);
      if (dummyLabels) {
        at1->setIsotope((*dummyLabels)[i].first);
        at2->setIsotope((*dummyLabels)[i].second);
      } else {
        at1->setIsotope(bidx);
        at2->setIsotope(eidx);
      }
      unsigned int idx1 = res->addAtom(at1, false, true);
      if (bondTypes) {
        bT = (*bondTypes)[i];
      }
      bondidx = res->addBond(at1->getIdx(), eidx, bT) - 1;
      // the dummy replaces the original start atom, so the
      // direction will be ok as long as it's one of the
      // states associated with double bond stereo
      if (bD == Bond::ENDDOWNRIGHT || bD == Bond::ENDUPRIGHT) {
        res->getBondWithIdx(bondidx)->setBondDir(bD);
      }

      unsigned int idx2 = res->addAtom(at2, false, true);
      bondidx = res->addBond(bidx, at2->getIdx(), bT) - 1;
      // this bond starts at the same atom, so its direction should always be
      // correct:
      res->getBondWithIdx(bondidx)->setBondDir(bD);

      // figure out if we need to change the stereo tags on the atoms:
      if (mol.getAtomWithIdx(bidx)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW ||
          mol.getAtomWithIdx(bidx)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW) {
        checkChiralityPostMove(mol, mol.getAtomWithIdx(bidx),
                               res->getAtomWithIdx(bidx),
                               mol.getBondBetweenAtoms(bidx, eidx));
      }
      if (mol.getAtomWithIdx(eidx)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CCW ||
          mol.getAtomWithIdx(eidx)->getChiralTag() ==
              Atom::CHI_TETRAHEDRAL_CW) {
        checkChiralityPostMove(mol, mol.getAtomWithIdx(eidx),
                               res->getAtomWithIdx(eidx),
                               mol.getBondBetweenAtoms(bidx, eidx));
      }

      for (auto confIt = res->beginConformers(); confIt != res->endConformers();
           ++confIt) {
        Conformer *conf = (*confIt).get();
        conf->setAtomPos(idx1, conf->getAtomPos(bidx));
        conf->setAtomPos(idx2, conf->getAtomPos(eidx));
      }
    } else {
      // was github issue 429
      Atom *tatom = res->getAtomWithIdx(bidx);
      if (tatom->getIsAromatic() && tatom->getAtomicNum() != 6) {
        tatom->setNumExplicitHs(tatom->getNumExplicitHs() + 1);
      }
      tatom = res->getAtomWithIdx(eidx);
      if (tatom->getIsAromatic() && tatom->getAtomicNum() != 6) {
        tatom->setNumExplicitHs(tatom->getNumExplicitHs() + 1);
      }
    }
  }
  res->clearComputedProps();
  return static_cast<ROMol *>(res);
}

ROMol *fragmentOnBonds(const ROMol &mol,
                       const std::vector<FragmenterBondType> &bondPatterns,
                       const std::map<unsigned int, ROMOL_SPTR> *atomEnvirons,
                       std::vector<unsigned int> *nCutsPerAtom) {
  PRECONDITION((!nCutsPerAtom || nCutsPerAtom->size() == mol.getNumAtoms()),
               "bad nCutsPerAtom vector");
  std::vector<unsigned int> bondIndices;
  std::vector<std::pair<unsigned int, unsigned int>> dummyLabels;
  std::vector<Bond::BondType> bondTypes;

  std::map<unsigned int, bool> environsMatch;
  if (atomEnvirons) {
    for (const auto &atomEnviron : *atomEnvirons) {
      MatchVectType mv;
      environsMatch[atomEnviron.first] =
          SubstructMatch(mol, *(atomEnviron.second), mv);
    }
  }

  boost::dynamic_bitset<> bondsUsed(mol.getNumBonds(), 0);
  // the bond definitions are organized (more or less) general -> specific, so
  // loop
  // over them backwards
  BOOST_REVERSE_FOREACH(const FragmenterBondType &fbt, bondPatterns) {
    if (fbt.query->getNumAtoms() != 2 || fbt.query->getNumBonds() != 1) {
      BOOST_LOG(rdErrorLog)
          << "fragmentation queries must have 2 atoms and 1 bond" << std::endl;
      continue;
    }
    if (atomEnvirons &&
        (!environsMatch[fbt.atom1Type] || !environsMatch[fbt.atom2Type])) {
      continue;
    }
    // std::cerr<<"  >>> "<<fbt.atom1Label<<" "<<fbt.atom2Label<<std::endl;
    std::vector<MatchVectType> bondMatches;
    SubstructMatch(mol, *fbt.query.get(), bondMatches);
    BOOST_FOREACH (const MatchVectType &mv, bondMatches) {
      const Bond *bond = mol.getBondBetweenAtoms(mv[0].second, mv[1].second);
      // std::cerr<<"          "<<bond->getIdx()<<std::endl;
      TEST_ASSERT(bond);
      if (bondsUsed[bond->getIdx()]) {
        // BOOST_LOG(rdWarningLog)<<"bond #"<<bond->getIdx()<<" matched multiple
        // times in decomposition. Later matches ignored."<<std::endl;
        continue;
      }
      bondsUsed.set(bond->getIdx());
      bondIndices.push_back(bond->getIdx());
      if (bond->getBeginAtomIdx() == static_cast<unsigned int>(mv[0].second)) {
        dummyLabels.emplace_back(fbt.atom1Label, fbt.atom2Label);
      } else {
        dummyLabels.emplace_back(fbt.atom2Label, fbt.atom1Label);
      }
      bondTypes.push_back(fbt.bondType);
    }
  }
  return fragmentOnBonds(mol, bondIndices, true, &dummyLabels, &bondTypes,
                         nCutsPerAtom);
}

boost::flyweight<std::vector<FragmenterBondType>,
                 boost::flyweights::no_tracking>
    bondPatterns;
boost::flyweight<std::map<unsigned int, ROMOL_SPTR>,
                 boost::flyweights::no_tracking>
    atomEnvs;
ROMol *fragmentOnBRICSBonds(const ROMol &mol) {
  if (bondPatterns.get().size() == 0) {
    std::map<unsigned int, std::string> adefs;
    std::map<unsigned int, ROMOL_SPTR> aenvs;
    constructBRICSAtomTypes(adefs, &aenvs);
    atomEnvs = aenvs;
    std::vector<FragmenterBondType> tbondPatterns;
    constructBRICSBondTypes(tbondPatterns);
    bondPatterns = tbondPatterns;
  }
  return fragmentOnBonds(mol, bondPatterns, &(atomEnvs.get()));
}
} // End of MolFragmenter

namespace {
// Get the atom label - this might be useful as a util class
unsigned int get_label(const Atom *a, const MolzipParams &p) {
  unsigned int idx = 0; // 0 is no label?
  switch(p.label) {
  case AtomLabel::AtomMapNumber:
    return a->getAtomMapNum();
  case AtomLabel::Isotope:
    return a->getIsotope();
  case AtomLabel::AtomType: {
    idx = std::distance(p.atomSymbols.begin(),
                        std::find(p.atomSymbols.begin(), p.atomSymbols.end(),
                                  a->getSymbol()));
    if(idx == p.atomSymbols.size()) {
      idx = 0; // not found since iterator == end()
    }
    break;
  }
  }
  return idx;
}
 
// Map the atoms based no the current label settings
//  label_id -> atom
// no two atoms int the same molecule can have the same label
std::map<unsigned int, const Atom*> get_mappings(const ROMol &a, const MolzipParams &p,
                                                 int mol) {
  std::map<unsigned int, const Atom*> mappings;
  for(auto *atom : a.atoms()) {
      atom->setProp<int>("mol", mol);
    unsigned int idx = get_label(atom, p);
    if(idx) {
        PRECONDITION(mappings.find(idx) == mappings.end(),
             "duplicate label in molecule, can't molzip");
        mappings[idx] = atom;
    }
  }
  return mappings;
}

// Return the connected atom
//  n.b. There can be only one connection from a mapped atom
const Atom *get_other_atom(const Atom *a) {
  auto &m = a->getOwningMol();
  if(m.getAtomDegree(a) != 1)
    return nullptr;
  
  for(const auto &nbrIdx : boost::make_iterator_range(m.getAtomNeighbors(a))) {
    return m[nbrIdx];
  }
  return nullptr; // can't get here
}

// Are two vectors simple rotations of each other
bool is_rotation(std::vector<unsigned int> orders1, const std::vector<unsigned int> &orders2) {
  assert(orders1.size() == orders2.size());//, "Order vectors aren't the same length");
  for(size_t i=0; i<orders1.size(); ++i) {
    std::rotate(orders1.begin(), orders1.begin()+1, orders1.end());
    if (orders1 == orders2)
      return true;
  }
  return false;
}

const std::string CHIRAL_TAG = "chiral_tag";
const std::string BOND_ME = "BOND_ME";
const std::string DELETE_ME = "DELETE_ME";

std::vector<unsigned int> get_marked_chiral_orders(const Atom *atom) {
    std::vector<unsigned int> orders;
    const std::string tag = atom->getProp<std::string>(CHIRAL_TAG);
    auto &mol = atom->getOwningMol();
    for(const auto &nbrIdx : boost::make_iterator_range(mol.getAtomNeighbors(atom))) {
        orders.push_back(mol[nbrIdx]->getProp<unsigned int>(tag));
    }
    return orders;
}

// mark the original orders of chirality around the atom with a tag unique to
//  the original atom
void mark_original_chiral_orders(const Atom *atom) {
  std::string tag = "chirality_around_atom_" + std::to_string(atom->getIdx());
  atom->setProp<std::string>(CHIRAL_TAG, tag);
  int order = 0;
  auto &mol = atom->getOwningMol();
  for(const auto &nbrIdx : boost::make_iterator_range(mol.getAtomNeighbors(atom))) {
    mol[nbrIdx]->setProp<unsigned int>(tag, order++);
  }
}
  
//  Save the original atom chiral orderings
struct ChiralBookmark {
  Atom *atom;
  std::vector<unsigned int> orders;
  ChiralBookmark() : atom(nullptr), orders() {
    
  }

  ChiralBookmark(Atom *atom) : atom(atom), orders(get_marked_chiral_orders(atom)) {
  }

  ChiralBookmark(const ChiralBookmark &rhs) : atom(rhs.atom), orders(rhs.orders) {
  }

  ChiralBookmark& operator=(const ChiralBookmark&rhs) {
        if(this != &rhs) {
            atom = rhs.atom;
            orders = rhs.orders;
        }
            return *this;
  }
  void restore_chirality() {
    // if our ordering is not a rotation of the orignial orders
    //  (with the dummy now a real atom) our chirality must have flipped.
    auto orders2 = get_marked_chiral_orders(atom);
    if (!is_rotation(orders, orders2))
      atom->invertChirality();
  }
};


// copy the chiral order from an atom being replaced (i.e. dummy -> C)
//  and return the bookmarked original chiralities prior to any reordering
void copy_chirality(Atom *chiral_atom, Atom *future_bond_atom,
			std::map<const Atom *, ChiralBookmark> &remap_chiral,
			std::vector<const Atom*> &deletions) {
  std::string tag;
  PRECONDITION(chiral_atom && chiral_atom->getPropIfPresent(CHIRAL_TAG, tag), "Atom is not chiral");
  
  remap_chiral[chiral_atom] = ChiralBookmark(chiral_atom);
  auto &m = chiral_atom->getOwningMol();
  for(const auto &nbrIdx : boost::make_iterator_range(m.getAtomNeighbors(chiral_atom))) {
    Atom *oatom = m[nbrIdx];
      int count = 0;
    if (std::find(deletions.begin(), deletions.end(), oatom) != deletions.end()) {
      auto order = oatom->getProp<unsigned int>(tag);
      future_bond_atom->setProp<unsigned int>(tag, order);
      count++;
    }
      if (count > 1) {
          BOOST_LOG(rdWarningLog) <<  "Can't handle more than one deletion around a chiral atom, chirality will be suspect." << std::endl;
      }
  }
}
}

std::unique_ptr<ROMol> molzip(
        const ROMol &a, const ROMol &b, const MolzipParams &params) {
  auto mapsa = get_mappings(a, params, 1);
  auto mapsb = get_mappings(b, params, 2);

  // mark the ations we need to perform to combine the molecules, basically we need
  //  to combine the mols, and bond the mapped atoms connected to the dummies and delete the dummies.
  //  i.e. C[*:1] + N[*:1] ==> CN  
  int count = 0;
  for(auto &kv : mapsa) {
    auto it = mapsb.find(kv.first);
    if(it != mapsb.end()) {
      auto atom_a = kv.second;
      auto anchor_a = get_other_atom(atom_a);
      auto atom_b = it->second;
      auto anchor_b = get_other_atom(atom_b);
      PRECONDITION(anchor_a, "Molzip Labelled atom in molecule a must only have one connection");
      PRECONDITION(anchor_b, "Molzip Labelled atom in molecule b must only have one connection");
      if(params.preserveChirality) {
        if(anchor_a->getChiralTag()) {
          mark_original_chiral_orders(anchor_a);
        }
        if(anchor_b->getChiralTag()) {
          mark_original_chiral_orders(anchor_b);
        }
      }

      count++;
      anchor_a->setProp<int>(BOND_ME, count);
      anchor_b->setProp<int>(BOND_ME, count);
      atom_a->setProp<int>(DELETE_ME, count);
      atom_b->setProp<int>(DELETE_ME, count);
    }
  }

  auto actions = dynamic_cast<RWMol*>(combineMols(a, b));
  assert (actions);

  if(count) {
    std::map<unsigned int, std::vector<Atom*>> bonds;
    std::vector<const Atom*> deletions;
    for(auto *atom : actions->atoms()) {
      if(atom->hasProp(BOND_ME)) {
	     bonds[atom->getProp<int>(BOND_ME)].push_back(atom);
      } else if ( atom->hasProp(DELETE_ME) ) {
          deletions.push_back(atom);
      }
    }
    
    std::map<const Atom *, ChiralBookmark> remap_chiral;
    // add the bonds and mark original chirality where appropriate
    for(auto &kv : bonds) {
      assert(kv.second.size() == 2);
      // PRECONODITION?
      Atom *a = kv.second[0];
      Atom *b = kv.second[1];
      if (a->hasProp(CHIRAL_TAG)) {
          copy_chirality(a, b, remap_chiral, deletions);
      }
      if(b->hasProp(CHIRAL_TAG)) {
          copy_chirality(b, a, remap_chiral, deletions);
      }
      actions->addBond(a->getIdx(), b->getIdx(), Bond::SINGLE);
    }

    for(auto &kv : deletions) {
      actions->removeAtom(kv->getIdx());
    }

    for(auto &kv : remap_chiral) {
      kv.second.restore_chirality();
    }
    actions->updatePropertyCache();
  }

  ROMol *m = dynamic_cast<ROMol*>(actions);
  return std::unique_ptr<ROMol>(m);
}
  
}  // end of namespace RDKit
