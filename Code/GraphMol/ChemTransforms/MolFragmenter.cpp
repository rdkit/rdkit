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

#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/utils.h>
#include "ChemTransforms.h"

#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/flyweight.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <boost/functional/hash.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/tokenizer.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <algorithm>
#include <cstdint>
#include <map>
#include <optional>
#include <sstream>
#include <vector>

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

std::vector<std::pair<Bond *, std::vector<int>>> getNbrBondStereo(
    RWMol &mol, const Bond *bnd) {
  PRECONDITION(bnd, "null bond");
  // loop over neighboring double bonds and remove their stereo atom
  std::vector<std::pair<Bond *, std::vector<int>>> res;
  const auto bgn = bnd->getBeginAtom();
  const auto end = bnd->getEndAtom();
  for (const auto *atom : {bgn, end}) {
    ROMol::OEDGE_ITER a1, a2;
    for (boost::tie(a1, a2) = mol.getAtomBonds(atom); a1 != a2; ++a1) {
      Bond *obnd = mol[*a1];
      if (obnd->getIdx() != bnd->getIdx() && !obnd->getStereoAtoms().empty()) {
        res.emplace_back(obnd, obnd->getStereoAtoms());
      }
    }
  }
  return res;
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
    for (auto &nCuts : *nCutsPerAtom) {
      nCuts = 0;
    }
  }
  auto *res = new RWMol(mol);
  if (!mol.getNumAtoms()) {
    return res;
  }

  std::vector<Bond *> bondsToRemove;
  bondsToRemove.reserve(bondIndices.size());
  for (auto bondIdx : bondIndices) {
    bondsToRemove.push_back(res->getBondWithIdx(bondIdx));
  }
  for (unsigned int i = 0; i < bondsToRemove.size(); ++i) {
    const Bond *bond = bondsToRemove[i];
    unsigned int bidx = bond->getBeginAtomIdx();
    unsigned int eidx = bond->getEndAtomIdx();
    Bond::BondType bT = bond->getBondType();
    Bond::BondDir bD = bond->getBondDir();
    unsigned int bondidx;
    auto nbr_bond_stereo = getNbrBondStereo(*res, bond);
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

      // restore stereo atoms
      for (auto &stereo_atoms : nbr_bond_stereo) {
        std::replace(stereo_atoms.second.begin(), stereo_atoms.second.end(),
                     bidx, idx1);
        std::replace(stereo_atoms.second.begin(), stereo_atoms.second.end(),
                     eidx, idx2);
        stereo_atoms.first->getStereoAtoms().swap(stereo_atoms.second);
      }

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
      // was github issues 429, 6034
      for (auto idx : {bidx, eidx}) {
        if (auto tatom = res->getAtomWithIdx(idx);
            tatom->getNoImplicit() ||
            (tatom->getIsAromatic() && tatom->getAtomicNum() != 6)) {
          tatom->setNumExplicitHs(tatom->getNumExplicitHs() + 1);
        } else {
          tatom->updatePropertyCache(false);
        }
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
  // loop over them backwards
  for (const auto &fbt : boost::adaptors::reverse(bondPatterns)) {
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
    for (const auto &mv : bondMatches) {
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
}  // namespace MolFragmenter

namespace {
const unsigned int NOLABEL = std::numeric_limits<unsigned int>::max();
// Get the atom label - this might be useful as a util class
unsigned int get_label(const Atom *a, const MolzipParams &p) {
  PRECONDITION(a, "bad atom in MolZip::get_label")
  unsigned int idx = NOLABEL;
  switch (p.label) {
    case MolzipLabel::AtomMapNumber:
      if (a->getAtomicNum() == 0) {
        auto mapno = a->getAtomMapNum();
        return mapno ? mapno : NOLABEL;
      }
      break;

    case MolzipLabel::Isotope:
      if (a->getAtomicNum() == 0) {
        auto iso = a->getIsotope();
        return iso ? iso : NOLABEL;
      }
      break;

    case MolzipLabel::AtomType:
      idx = std::distance(p.atomSymbols.begin(),
                          std::find(p.atomSymbols.begin(), p.atomSymbols.end(),
                                    a->getSymbol()));
      if (idx == p.atomSymbols.size()) {
        idx = NOLABEL;
      }
      break;

    case MolzipLabel::FragmentOnBonds:
      // shouldn't ever get here
      CHECK_INVARIANT(
          0, "FragmentOnBonds is not an atom label, it is an atom index");
      break;

    case MolzipLabel::AtomProperty:
      a->getPropIfPresent<unsigned int>(p.atomProperty, idx);
      break;

    default:
      CHECK_INVARIANT(0, "bogus MolZipLabel value in MolZip::get_label");
  }
  return idx;
}

// Return the connected atom
//  n.b. There can be only one connection from a mapped atom
Atom *get_other_atom(Atom *a) {
  PRECONDITION(a, "null atom in MolZip::get_other_atom");
  auto &m = a->getOwningMol();
  if (m.getAtomDegree(a) != 1) {
    return nullptr;
  }

  return m[*m.getAtomNeighbors(a).first];
}

int num_swaps_to_interconvert(std::vector<unsigned int> &orders) {
  int nswaps = 0;
  std::vector<bool> seen(orders.size());
  for (size_t i = 0; i < orders.size(); ++i) {
    if (!seen[i]) {
      auto j = i;
      while (orders[j] != i) {
        j = orders[j];
        CHECK_INVARIANT(
            j < orders.size(),
            "molzip: bond index outside of number of bonds for atom")
        seen[j] = true;
        nswaps++;
      }
    }
  }
  return nswaps;
}

// Simple bookkeeping class to bond attachments and handle stereo
struct ZipBond {
  Atom *a;        // atom being bonded
  Atom *a_dummy;  // Labelled atom, i.e. [*:1]-C  will bond the C to something
  Atom *b;        // atom being bonded
  Atom *b_dummy;  // Labelled atom, i.e. [*:1]-O will bond the O to something

  // Backup the original chirality mark_chirality must be called first
  //  as it checks the datastructure for validity;
  void mark_chirality() const {
    PRECONDITION(a, "Must have a begin atom to bond");
    PRECONDITION(b, "Must have an end atom to bond");
    PRECONDITION(a_dummy, "Must have a begin dummy atom");
    PRECONDITION(b_dummy, "Must have an end dummy atom");

    mark(a, a_dummy, b);
    mark(b, b_dummy, a);
  }

  // bond a<->b for now only use single bonds
  //  XXX FIX ME take the highest bond order.
  bool bond(RWMol &newmol, const MolzipParams &params) const {
    if (!a || !b || !a_dummy || !b_dummy) {
      BOOST_LOG(rdWarningLog)
          << "Incomplete atom labelling, cannot make bond" << std::endl;
      return false;
    }

    // Fragment on bonds allows multiple links to the same atom
    // i.e. C.[1C].[1C]
    //  otherwise throw an invariant error
    CHECK_INVARIANT(
        params.label == MolzipLabel::FragmentOnBonds ||
            !a->getOwningMol().getBondBetweenAtoms(a->getIdx(), b->getIdx()),
        "molzip: zipped Bond already exists, perhaps labels are duplicated");

    if (!a->getOwningMol().getBondBetweenAtoms(a->getIdx(), b->getIdx())) {
      CHECK_INVARIANT(&a->getOwningMol() == &newmol,
                      "Owning mol is not the combined molecule!!");
      auto bnd = newmol.getBondBetweenAtoms(a->getIdx(), a_dummy->getIdx());
      CHECK_INVARIANT(bnd != nullptr,
                      "molzip: begin atom and specified dummy atom connection "
                      "are not bonded.")
      auto bond_type_a = bnd->getBondType();
      auto bond_dir_a = bnd->getBondDir();
      auto a_is_start = bnd->getBeginAtom() == a;

      bnd = newmol.getBondBetweenAtoms(b->getIdx(), b_dummy->getIdx());
      CHECK_INVARIANT(bnd != nullptr,
                      "molzip: end atom and specified dummy connection atom "
                      "are not bonded.")
      auto bond_type_b = bnd->getBondType();

      auto bond_dir_b = bnd->getBondDir();
      auto b_is_start = bnd->getBeginAtom() == b;

      unsigned int bnd_idx = 0;
      // Fusion bond-dir logic table
      // a-* b-* => a-b
      //  < = wedge

      //  a<* b-* => a<b
      //  a>* b-* => a>b
      //  a-* b>* => a<b
      //  a-* b<* => a>b
      Bond::BondDir bond_dir{Bond::BondDir::NONE};
      auto start = a;
      auto end = b;
      if (bond_dir_a != Bond::BondDir::NONE &&
          bond_dir_b != Bond::BondDir::NONE) {
        // are we consistent between the two bond orders check for the case of
        // fragment on bonds where a<* and b>* or a>* and b<* when < is either a
        // hash or wedge bond but not both.
        bool consistent_directions = false;
        if (bond_dir_a == bond_dir_b) {
          if ((a_is_start != b_is_start)) {
            consistent_directions = true;
          }
        }
        if (!consistent_directions) {
          BOOST_LOG(rdWarningLog)
              << "inconsistent bond directions when merging fragments, ignoring..."
              << std::endl;
          bond_dir_a = bond_dir_b = Bond::BondDir::NONE;
        } else {
          bond_dir_b = Bond::BondDir::NONE;
        }
      }

      if (bond_dir_a != Bond::BondDir::NONE) {
        if (!a_is_start) {
          start = b;
          end = a;
        }
        bond_dir = bond_dir_a;
      } else if (bond_dir_b != Bond::BondDir::NONE) {
        if (b_is_start) {
          start = b;
          end = a;
        }
        bond_dir = bond_dir_b;
      }

      if (bond_type_a != Bond::BondType::SINGLE) {
        bnd_idx = newmol.addBond(start, end, bond_type_a);
      } else if (bond_type_b != Bond::BondType::SINGLE) {
        bnd_idx = newmol.addBond(start, end, bond_type_b);
      } else {
        bnd_idx = newmol.addBond(start, end, Bond::BondType::SINGLE);
      }

      newmol.getBondWithIdx(bnd_idx - 1)->setBondDir(bond_dir);
    }
    a_dummy->setProp("__molzip_used", true);
    b_dummy->setProp("__molzip_used", true);

    return true;
  }

  // Restore the marked chirality (mark_chirality must be called first)
  void restore_chirality(std::set<Atom *> &already_checked) const {
    PRECONDITION(a, "Must have a begin atom to bond");
    PRECONDITION(b, "Must have an end atom to bond");
    PRECONDITION(a_dummy, "Must have a begin dummy atom");
    PRECONDITION(b_dummy, "Must have an end dummy atom");
    if (already_checked.find(a) == already_checked.end()) {
      restore(a);
      already_checked.insert(a);
    }
    if (already_checked.find(b) == already_checked.end()) {
      restore(b);
      already_checked.insert(b);
    }

    // now do bond stereo
    std::string mark = "__molzip_bond_stereo_mark";
    for (auto *bond : a->getOwningMol().bonds()) {
      if (bond->hasProp(mark)) {
        std::vector<int> atoms;
        for (auto *atom : bond->getProp<std::vector<Atom *>>(mark)) {
          atoms.push_back(rdcast<int>(atom->getIdx()));
        }
        bond->getStereoAtoms().swap(atoms);
        bond->setStereo(
            bond->getProp<Bond::BondStereo>("__molzip_bond_stereo"));
      }
    }
  }

 private:
  // Mark the original order of the nbr atoms including the dummy
  //  The goal is to copy the dummy chiral order over to the
  //  atom being bonded
  void mark(Atom *chiral_atom, Atom *dummy_atom, Atom *new_atom) const {
    if (chiral_atom->getChiralTag()) {
      std::string mark =
          "__molzip_mark_" + std::to_string(chiral_atom->getIdx());
      chiral_atom->setProp("__molzip_chiral_mark", mark);
      int order = 0;
      auto &m = chiral_atom->getOwningMol();
      for (auto nbrIdx :
           boost::make_iterator_range(m.getAtomNeighbors(chiral_atom))) {
        m[nbrIdx]->setProp(mark, order);
        ++order;
      }
      new_atom->setProp(mark, dummy_atom->getProp<int>(mark));
    }

    // check bond stereo
    auto &m = chiral_atom->getOwningMol();
    for (auto nbrIdx :
         boost::make_iterator_range(m.getAtomNeighbors(chiral_atom))) {
      auto bond = m.getBondBetweenAtoms(chiral_atom->getIdx(), nbrIdx);
      if (bond->getStereo()) {
        std::string mark = "__molzip_bond_stereo_mark";
        std::vector<Atom *> atoms;
        bool has_dummy = false;
        for (auto idx : bond->getStereoAtoms()) {
          if (static_cast<unsigned>(idx) == dummy_atom->getIdx()) {
            atoms.push_back(new_atom);
            has_dummy = true;
          } else {
            atoms.push_back(m.getAtomWithIdx(idx));
          }
        }
        if (has_dummy) {
          bond->setProp(mark, atoms);
          bond->setProp<Bond::BondStereo>("__molzip_bond_stereo",
                                          bond->getStereo());
        }
      }
    }
  }

  // Restore the atom's chirality by comparing the original order
  //  to the current
  void restore(Atom *chiral_atom) const {
    if (!chiral_atom->getChiralTag()) {
      return;
    }
    std::string mark =
        chiral_atom->getProp<std::string>("__molzip_chiral_mark");
    // std::vector<unsigned int> orders1;
    std::vector<unsigned int> orders2;
    auto &m = chiral_atom->getOwningMol();
    for (auto nbrIdx :
         boost::make_iterator_range(m.getAtomNeighbors(chiral_atom))) {
      orders2.push_back(m[nbrIdx]->getProp<int>(mark));
    }
    if (num_swaps_to_interconvert(orders2) % 2 == 1) {
      chiral_atom->invertChirality();
    }
  }
};
}  // namespace

static const std::string indexPropName("__zipIndex");

std::unique_ptr<ROMol> molzip(
    const ROMol &a, const ROMol &b, const MolzipParams &params,
    std::optional<std::map<int, int>> &attachmentMapping) {
  if (attachmentMapping) {
    attachmentMapping->clear();
  }
  std::unique_ptr<RWMol> newmol;
  if (b.getNumAtoms()) {
    newmol.reset(static_cast<RWMol *>(combineMols(a, b)));
  } else {
    newmol.reset(new RWMol(a));
  }

  std::map<unsigned int, ZipBond> mappings;
  std::map<Atom *, std::vector<const ZipBond *>> mappings_by_atom;
  std::vector<Atom *> deletions;
  if (params.label == MolzipLabel::FragmentOnBonds) {
    for (auto *atom : newmol->atoms()) {
      if (atom->getAtomicNum() == 0) {
        auto molno = atom->getIsotope();
        auto attached_atom = get_other_atom(atom);
        auto &bond = mappings[molno];
        bond.a = attached_atom;
        bond.a_dummy = atom;
        bond.b = newmol->getAtomWithIdx(molno);
        for (auto nbrIdx :
             boost::make_iterator_range(newmol->getAtomNeighbors(bond.b))) {
          auto *nbr = (*newmol)[nbrIdx];
          if (nbr->getAtomicNum() == 0 &&
              nbr->getIsotope() == attached_atom->getIdx()) {
            bond.b_dummy = nbr;
            break;
          }
        }
        if (!bond.b_dummy) {
          BOOST_LOG(rdErrorLog)
              << "Cannot find atom to bond using FragmentOnBond labelling"
              << std::endl;
          return std::unique_ptr<ROMol>();
        }
        mappings_by_atom[atom].push_back(&bond);
        deletions.push_back(atom);
        if (attachmentMapping) {
          if (int otherIndex, dummyIndex;
              atom->getPropIfPresent(indexPropName, dummyIndex) &&
              bond.b->getPropIfPresent(indexPropName, otherIndex)) {
            (*attachmentMapping)[dummyIndex] = otherIndex;
          }
        }
      }
    }
  } else {
    for (auto *atom : newmol->atoms()) {
      auto molno = get_label(atom, params);
      if (molno != NOLABEL) {
        auto attached_atom = get_other_atom(atom);
        if (mappings.find(molno) == mappings.end()) {
          auto &bond = mappings[molno];
          CHECK_INVARIANT(
              !bond.a,
              "molzip: bond info already setup for bgn atom with label:" +
                  std::to_string(molno));
          bond.a = attached_atom;
          bond.a_dummy = atom;
        } else {
          auto &bond = mappings[molno];
          CHECK_INVARIANT(
              bond.a,
              "molzip: bond info not properly setup for bgn atom with label:" +
                  std::to_string(molno));
          CHECK_INVARIANT(
              !bond.b,
              "molzip: bond info already exists for end atom with label:" +
                  std::to_string(molno));
          bond.b = attached_atom;
          bond.b_dummy = atom;
          mappings_by_atom[bond.a].push_back(&bond);
          if (attachmentMapping) {
            if (int otherIndex, dummyIndex;
                bond.a_dummy->getPropIfPresent(indexPropName, dummyIndex) &&
                bond.b->getPropIfPresent(indexPropName, otherIndex)) {
              (*attachmentMapping)[dummyIndex] = otherIndex;
            }
            if (int otherIndex, dummyIndex;
                bond.b_dummy->getPropIfPresent(indexPropName, dummyIndex) &&
                bond.a->getPropIfPresent(indexPropName, otherIndex)) {
              (*attachmentMapping)[dummyIndex] = otherIndex;
            }
          }
        }
        deletions.push_back(atom);
      }
    }
  }
  for (auto &kv : mappings_by_atom) {
    for (auto &bond : kv.second) {
      bond->mark_chirality();
    }
  }
  for (auto &kv : mappings) {
    kv.second.bond(*newmol, params);
  }
  newmol->beginBatchEdit();
  for (auto &atom : deletions) {
    if (atom->hasProp("__molzip_used")) {
      newmol->removeAtom(atom);
    }
  }
  newmol->commitBatchEdit();

  std::set<Atom *> already_checked;
  for (auto &kv : mappings_by_atom) {
    for (auto &bond : kv.second) {
      bond->restore_chirality(already_checked);
    }
  }

  // remove all molzip tags
  for (auto *atom : newmol->atoms()) {
    auto propnames = atom->getPropList();
    for (auto &prop : propnames) {
      if (prop.find("__molzip") == 0) {
        atom->clearProp(prop);
      }
    }
  }
  for (auto *bond : newmol->bonds()) {
    auto propnames = bond->getPropList();
    for (auto &prop : propnames) {
      if (prop.find("__molzip") == 0) {
        bond->clearProp(prop);
      }
    }
  }
  newmol->updatePropertyCache(params.enforceValenceRules);
  newmol->setProp(common_properties::_StereochemDone, true);
  return newmol;
}

RDKIT_CHEMTRANSFORMS_EXPORT std::unique_ptr<ROMol> molzip(
    const ROMol &a, const ROMol &b, const MolzipParams &params) {
  std::optional<std::map<int, int>> opt(std::nullopt);
  return molzip(a, b, params, opt);
}

std::unique_ptr<ROMol> molzip(const ROMol &a, const MolzipParams &params) {
  const static ROMol b;
  return molzip(a, b, params);
}

std::unique_ptr<ROMol> molzip(std::vector<ROMOL_SPTR> &decomposition,
                              const MolzipParams &params) {
  if (params.generateCoordinates) {
    int index = 0;
    for (const auto &mol : decomposition) {
      for (const auto atom : mol->atoms()) {
        atom->setProp(indexPropName, index++);
      }
    }
  }

  const auto combinedMol = std::accumulate(
      decomposition.begin() + 1, decomposition.end(), decomposition[0],
      [](const auto &combined, const auto &mol) {
        auto c = combineMols(*combined, *mol);
        ROMOL_SPTR ptr(c);
        return ptr;
      });

  const static ROMol b;
  std::optional attachmentMappingOption = std::map<int, int>();
  auto zippedMol = molzip(*combinedMol, b, params, attachmentMappingOption);

  if (params.generateCoordinates && zippedMol->getNumAtoms() > 0) {
    const auto confId = RDDepict::compute2DCoords(*zippedMol);
    const auto zippedConf = zippedMol->getConformer(confId);
    auto attachmentMapping = *attachmentMappingOption;
    for (auto &mol : decomposition) {
      const auto newConf = new Conformer(mol->getNumAtoms());
      newConf->set3D(false);
      for (const auto atom : mol->atoms()) {
        int zippedIndex = atom->getProp<int>(indexPropName);
        atom->clearProp(indexPropName);
        if (const auto attachment = attachmentMapping.find(zippedIndex);
            attachment != attachmentMapping.end()) {
          zippedIndex = (*attachment).second;
        }
        auto zipppedAtoms = zippedMol->atoms();
        auto zippedAtom = std::find_if(
            zipppedAtoms.begin(), zipppedAtoms.end(),
            [zippedIndex](const Atom *zippedAtom) {
              const auto index = zippedAtom->getProp<int>(indexPropName);
              return index == zippedIndex;
            });

        newConf->setAtomPos(atom->getIdx(),
                            zippedConf.getAtomPos((*zippedAtom)->getIdx()));
      }
      mol->addConformer(newConf, true);
    }
    for (const auto atom : zippedMol->atoms()) {
      atom->clearProp(indexPropName);
    }

    return zippedMol;
  }

  return zippedMol;
}

}  // end of namespace RDKit
