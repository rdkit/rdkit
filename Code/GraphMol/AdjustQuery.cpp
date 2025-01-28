//
//  Copyright (c) 2015-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/GraphMol.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <vector>
#include <algorithm>

namespace RDKit {
namespace {
bool isMapped(const Atom *atom) {
  return atom->hasProp(common_properties::molAtomMapNumber);
}
}  // namespace

namespace MolOps {

namespace {
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

unsigned int parseWhichString(const std::string &txt) {
  unsigned int res = MolOps::ADJUST_IGNORENONE;
  boost::char_separator<char> sep("|");
  tokenizer tokens(txt, sep);
  for (const auto &token : tokens) {
    if (token == "IGNORENONE") {
      res |= MolOps::ADJUST_IGNORENONE;
    } else if (token == "IGNORERINGS") {
      res |= MolOps::ADJUST_IGNORERINGS;
    } else if (token == "IGNORECHAINS") {
      res |= MolOps::ADJUST_IGNORECHAINS;
    } else if (token == "IGNOREDUMMIES") {
      res |= MolOps::ADJUST_IGNOREDUMMIES;
    } else if (token == "IGNORENONDUMMIES") {
      res |= MolOps::ADJUST_IGNORENONDUMMIES;
    } else if (token == "IGNOREALL") {
      res |= MolOps::ADJUST_IGNOREALL;
    } else {
      std::string msg = "Unknown flag value: '" + token + "'. Flags ignored.";
      throw ValueErrorException(msg);
    }
  }
  return res;
}
constexpr const char *conjugatedOrAromatic = "_conjugatedOrAromatic";
void adjustConjugatedFiveRings(RWMol &mol) {
  /*
   The idea here is to allow conjugated five-rings to match either aromatic or
   aliphatic rings

   five-rings which contain at least 3 conjugated bonds have all of their
   non-query bonds replaced with a SINGLE|DOUBLE|AROMATIC query.

   */

  std::vector<Bond::BondType> bondTypesToModify = {
      Bond::BondType::SINGLE, Bond::BondType::DOUBLE, Bond::BondType::AROMATIC};
  if (!mol.getRingInfo()->isSymmSssr()) {
    MolOps::symmetrizeSSSR(mol);
  }
  for (auto ring : mol.getRingInfo()->bondRings()) {
    // only consider 5-rings with at least 3 conjugated bonds
    if (ring.size() != 5) {
      continue;
    }
    unsigned int nconj = 0;
    for (auto bi : ring) {
      const auto bond = mol.getBondWithIdx(bi);
      if (bond->getIsConjugated()) {
        ++nconj;
        if (nconj >= 3) {
          break;
        }
      }
    }
    if (nconj < 3) {
      continue;
    }
    // now make the adjustments
    QueryBond qb;
    qb.setQuery(makeSingleOrDoubleOrAromaticBondQuery());
    for (auto bi : ring) {
      const auto bond = mol.getBondWithIdx(bi);
      bond->getBeginAtom()->setProp(conjugatedOrAromatic, 1, true);
      bond->getEndAtom()->setProp(conjugatedOrAromatic, 1, true);
      if (std::find(bondTypesToModify.begin(), bondTypesToModify.end(),
                    bond->getBondType()) != bondTypesToModify.end()) {
        if (bond->hasQuery()) {
          BOOST_LOG(rdWarningLog)
              << "adjustConjugatedFiveRings: replacing a bond "
                 "that already has a query"
              << std::endl;
        }
        mol.replaceBond(bi, &qb);
      }
    }
  }
}

bool isAromaticOrConjugated(const Atom &atom) {
  return atom.getIsAromatic() || atom.hasProp(conjugatedOrAromatic);
}
void adjustSingleBondsFromAromaticAtoms(RWMol &mol, bool toDegreeOneNeighbors,
                                        bool betweenAromaticAtoms) {
  /*
  The idea here is to allow single bonds coming from aromatic atoms to match
  aromatic bonds under particular circumstances. The conditions are:

  1. toDegreeOneNeighbors: [D1]-[a] -> [D1]-,:[a]
  2. betweenAromaticAtoms: [a]-[a] -> [a]-,:[a]

  */
  if (!toDegreeOneNeighbors && !betweenAromaticAtoms) {
    return;
  }
  QueryBond qb;
  qb.setQuery(makeSingleOrAromaticBondQuery());
  if (!mol.getRingInfo()->isSymmSssr()) {
    MolOps::symmetrizeSSSR(mol);
  }
  for (auto bond : mol.bonds()) {
    const auto bAt = bond->getBeginAtom();
    const auto eAt = bond->getEndAtom();
    if (!bond->hasQuery() && bond->getBondType() == Bond::BondType::SINGLE) {
      auto bAtIsAromatic = isAromaticOrConjugated(*bAt);
      auto eAtIsAromatic = isAromaticOrConjugated(*eAt);

      if (toDegreeOneNeighbors && (bAtIsAromatic ^ eAtIsAromatic)) {
        if ((bAtIsAromatic && eAt->getDegree() == 1) ||
            (eAtIsAromatic && bAt->getDegree() == 1)) {
          mol.replaceBond(bond->getIdx(), &qb);
        }
      } else if (betweenAromaticAtoms && bAtIsAromatic && eAtIsAromatic) {
        mol.replaceBond(bond->getIdx(), &qb);
      }
    }
  }
}

void setMDLAromaticity(RWMol &mol) {
  /*
  The idea here is to make aromatic 5-rings that contain an "A" atom in the CTAB
  match both aromatic and aliphatic rings.
  Schematically, this converts the ring from:
     ["A"]1:c:c:c:c:1
  to:
     ["A"]1-,:c=,:c-,:c=,:c-,:1
  Note that "A" is an A atom from a CTAB, not a SMARTS aliphatic query

  */

  // it would be simpler to use the substructure matcher for this, but we can't
  // use SubstructMatch in the core GraphMol lib
  if (!mol.getRingInfo()->isSymmSssr()) {
    MolOps::symmetrizeSSSR(mol);
  }
  for (auto ring : mol.getRingInfo()->atomRings()) {
    if (ring.size() != 5) {
      continue;
    }
    bool keepIt = true;
    size_t dummy = ring.size() + 1;
    for (size_t i = 0; i < ring.size(); ++i) {
      auto ai = ring[i];
      const auto atom = mol.getAtomWithIdx(ai);
      if (!atom->getIsAromatic()) {
        // we only do fully aromatic rings:
        keepIt = false;
        break;
      } else if (atom->getAtomicNum() == 0 && atom->hasQuery() &&
                 atom->getQuery()->getTypeLabel() == "A") {
        if (dummy >= ring.size()) {
          dummy = i;
        } else {
          // second dummy encountered, we won't do this ring.
          keepIt = false;
          break;
        }
      } else if (atom->getAtomicNum() != 6) {
        // we only do rings consisting solely of C and *
        keepIt = false;
        break;
      }
      // we can't handle rings that have query bonds already:
      auto oidx = ring[4];
      if (i > 0) {
        oidx = ring[i - 1];
      }
      auto bond = mol.getBondBetweenAtoms(ring[i], oidx);
      ASSERT_INVARIANT(bond, "expected bond not found");
      if (bond->hasQuery()) {
        keepIt = false;
        break;
      }
    }
    if (keepIt && dummy < ring.size()) {
      // we think about the 5-ring in three layers:
      //   layer 0: the dummy
      //   layer 1: the two atoms connected to the dummy
      //   layer 2: the two atoms not connected to the dummy
      auto l0 = ring[dummy];
      std::vector<int> l1;
      std::vector<int> l2;
      for (auto ai : ring) {
        if (ai == l0) {
          continue;
        } else if (mol.getBondBetweenAtoms(ai, l0)) {
          l1.push_back(ai);
        } else {
          l2.push_back(ai);
        }
      }
      ASSERT_INVARIANT(l1.size() == 2, "bad layer 1 size");
      ASSERT_INVARIANT(l2.size() == 2, "bad layer 2 size");

      QueryBond qbSingleAromatic;
      {
        BOND_OR_QUERY *q = new BOND_OR_QUERY;
        q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
            makeBondOrderEqualsQuery(Bond::SINGLE)));
        q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
            makeBondOrderEqualsQuery(Bond::AROMATIC)));
        q->setDescription("BondOr");
        qbSingleAromatic.setQuery(q);
      }
      QueryBond qbDoubleAromatic;
      {
        BOND_OR_QUERY *q = new BOND_OR_QUERY;
        q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
            makeBondOrderEqualsQuery(Bond::DOUBLE)));
        q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
            makeBondOrderEqualsQuery(Bond::AROMATIC)));
        q->setDescription("BondOr");
        qbDoubleAromatic.setQuery(q);
      }
      for (auto ai : l1) {
        // l0 - l1 bonds:
        auto bond = mol.getBondBetweenAtoms(ai, l0);
        ASSERT_INVARIANT(bond, "expected l0-l1 bond not found");
        mol.replaceBond(bond->getIdx(), &qbSingleAromatic);
        // l1 - l2 bonds:
        bond = mol.getBondBetweenAtoms(ai, l2[0]);
        if (!bond) {
          bond = mol.getBondBetweenAtoms(ai, l2[1]);
        }
        ASSERT_INVARIANT(bond, "expected l1-l2 bond not found");
        mol.replaceBond(bond->getIdx(), &qbDoubleAromatic);
      }
      // l2 - l2 bond:
      auto bond = mol.getBondBetweenAtoms(l2[0], l2[1]);
      ASSERT_INVARIANT(bond, "expected l2-l2 bond not found");
      mol.replaceBond(bond->getIdx(), &qbSingleAromatic);
    }
  }
}
}  // namespace
void parseAdjustQueryParametersFromJSON(MolOps::AdjustQueryParameters &p,
                                        const std::string &json) {
  PRECONDITION(!json.empty(), "empty JSON provided");
  std::istringstream ss;
  ss.str(json);

  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  p.adjustDegree = pt.get("adjustDegree", p.adjustDegree);
  p.adjustRingCount = pt.get("adjustRingCount", p.adjustRingCount);
  p.makeDummiesQueries = pt.get("makeDummiesQueries", p.makeDummiesQueries);
  p.aromatizeIfPossible = pt.get("aromatizeIfPossible", p.aromatizeIfPossible);
  p.makeBondsGeneric = pt.get("makeBondsGeneric", p.makeBondsGeneric);
  p.makeAtomsGeneric = pt.get("makeAtomsGeneric", p.makeAtomsGeneric);
  p.adjustHeavyDegree = pt.get("adjustHeavyDegree", p.adjustHeavyDegree);
  p.adjustRingChain = pt.get("adjustRingChain", p.adjustRingChain);
  p.useStereoCareForBonds =
      pt.get("useStereoCareForBonds", p.useStereoCareForBonds);
  p.adjustConjugatedFiveRings =
      pt.get("adjustConjugatedFiveRings", p.adjustConjugatedFiveRings);
  p.setMDLFiveRingAromaticity =
      pt.get("setMDLFiveRingAromaticity", p.setMDLFiveRingAromaticity);
  p.adjustSingleBondsToDegreeOneNeighbors =
      pt.get("adjustSingleBondsToDegreeOneNeighbors",
             p.adjustSingleBondsToDegreeOneNeighbors);
  p.adjustSingleBondsBetweenAromaticAtoms =
      pt.get("adjustSingleBondsBetweenAromaticAtoms",
             p.adjustSingleBondsBetweenAromaticAtoms);

  std::string which;
  which = boost::to_upper_copy<std::string>(pt.get("adjustDegreeFlags", ""));
  if (!which.empty()) {
    p.adjustDegreeFlags = parseWhichString(which);
  }
  which =
      boost::to_upper_copy<std::string>(pt.get("adjustHeavyDegreeFlags", ""));
  if (!which.empty()) {
    p.adjustHeavyDegreeFlags = parseWhichString(which);
  }
  which = boost::to_upper_copy<std::string>(pt.get("adjustRingCountFlags", ""));
  if (!which.empty()) {
    p.adjustRingCountFlags = parseWhichString(which);
  }
  which =
      boost::to_upper_copy<std::string>(pt.get("makeBondsGenericFlags", ""));
  if (!which.empty()) {
    p.makeBondsGenericFlags = parseWhichString(which);
  }
  which =
      boost::to_upper_copy<std::string>(pt.get("makeAtomsGenericFlags", ""));
  if (!which.empty()) {
    p.makeAtomsGenericFlags = parseWhichString(which);
  }
  which = boost::to_upper_copy<std::string>(pt.get("adjustRingChainFlags", ""));
  if (!which.empty()) {
    p.adjustRingChainFlags = parseWhichString(which);
  }
}  // namespace MolOps

ROMol *adjustQueryProperties(const ROMol &mol,
                             const AdjustQueryParameters *params) {
  auto *res = new RWMol(mol);
  try {
    adjustQueryProperties(*res, params);
  } catch (MolSanitizeException &se) {
    delete res;
    throw se;
  }
  return static_cast<ROMol *>(res);
}

void adjustQueryProperties(RWMol &mol, const AdjustQueryParameters *inParams) {
  AdjustQueryParameters params;
  if (inParams) {
    params = *inParams;
  }
  const auto ringInfo = mol.getRingInfo();

  if (params.aromatizeIfPossible) {
    unsigned int failed;
    sanitizeMol(mol, failed, SANITIZE_SYMMRINGS | SANITIZE_SETAROMATICITY);
  } else {
    if (!ringInfo->isSymmSssr()) {
      MolOps::symmetrizeSSSR(mol);
    }
  }
  QueryAtom qaTmpl;
  QueryBond qbTmpl;

  if (params.makeAtomsGeneric) {
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      if (!((params.makeAtomsGenericFlags & ADJUST_IGNORECHAINS) &&
            !ringInfo->numAtomRings(i)) &&
          !((params.makeAtomsGenericFlags & ADJUST_IGNORERINGS) &&
            ringInfo->numAtomRings(i)) &&
          !((params.makeAtomsGenericFlags & ADJUST_IGNOREMAPPED) &&
            isMapped(mol.getAtomWithIdx(i)))) {
        qaTmpl.setQuery(makeAtomNullQuery());
        constexpr bool updateLabel = false;
        constexpr bool preserveProps = true;
        mol.replaceAtom(i, &qaTmpl, updateLabel, preserveProps);
      }
    }
  }  // end of makeAtomsGeneric
  if (params.makeBondsGeneric) {
    for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
      if (!((params.makeBondsGenericFlags & ADJUST_IGNORECHAINS) &&
            !ringInfo->numBondRings(i)) &&
          !((params.makeBondsGenericFlags & ADJUST_IGNORERINGS) &&
            ringInfo->numBondRings(i))) {
        qbTmpl.setQuery(makeBondNullQuery());
        constexpr bool preserveProps = true;
        mol.replaceBond(i, &qbTmpl, preserveProps);
      }
    }
  }  // end of makeBondsGeneric
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    auto *at = mol.getAtomWithIdx(i);
    // pull properties we need from the atom here, once we
    // create a query atom they may no longer be valid.
    auto nRings = ringInfo->numAtomRings(i);
    auto atomicNum = at->getAtomicNum();
    if (params.makeDummiesQueries && atomicNum == 0 && !at->hasQuery() &&
        !at->getIsotope()) {
      qaTmpl.setQuery(makeAtomNullQuery());
      constexpr bool updateLabel = false;
      constexpr bool preserveProps = true;
      mol.replaceAtom(i, &qaTmpl, updateLabel, preserveProps);
      at = mol.getAtomWithIdx(i);
    }  // end of makeDummiesQueries
    if (params.adjustDegree &&
        !((params.adjustDegreeFlags & ADJUST_IGNORECHAINS) && !nRings) &&
        !((params.adjustDegreeFlags & ADJUST_IGNORERINGS) && nRings) &&
        !((params.adjustDegreeFlags & ADJUST_IGNOREDUMMIES) && !atomicNum) &&
        !((params.adjustDegreeFlags & ADJUST_IGNORENONDUMMIES) && atomicNum) &&
        !((params.adjustDegreeFlags & ADJUST_IGNOREMAPPED) && isMapped(at))) {
      QueryAtom *qa;
      if (!at->hasQuery()) {
        QueryAtom atQueryAtom(*at);
        constexpr bool updateLabel = false;
        constexpr bool preserveProps = true;
        mol.replaceAtom(i, &atQueryAtom, updateLabel, preserveProps);
        qa = static_cast<QueryAtom *>(mol.getAtomWithIdx(i));
        at = static_cast<Atom *>(qa);
      } else {
        qa = static_cast<QueryAtom *>(at);
      }
      qa->expandQuery(makeAtomExplicitDegreeQuery(qa->getDegree()));
    }  // end of adjust degree
    if (params.adjustHeavyDegree &&
        !((params.adjustHeavyDegreeFlags & ADJUST_IGNORECHAINS) && !nRings) &&
        !((params.adjustHeavyDegreeFlags & ADJUST_IGNORERINGS) && nRings) &&
        !((params.adjustHeavyDegreeFlags & ADJUST_IGNOREDUMMIES) &&
          !atomicNum) &&
        !((params.adjustHeavyDegreeFlags & ADJUST_IGNORENONDUMMIES) &&
          atomicNum) &&
        !((params.adjustHeavyDegreeFlags & ADJUST_IGNOREMAPPED) &&
          isMapped(at))) {
      QueryAtom *qa;
      if (!at->hasQuery()) {
        QueryAtom atQueryAtom(*at);
        constexpr bool updateLabel = false;
        constexpr bool preserveProps = true;
        mol.replaceAtom(i, &atQueryAtom, updateLabel, preserveProps);
        qa = static_cast<QueryAtom *>(mol.getAtomWithIdx(i));
        at = static_cast<Atom *>(qa);
      } else {
        qa = static_cast<QueryAtom *>(at);
      }
      qa->expandQuery(makeAtomHeavyAtomDegreeQuery(qa->getTotalDegree() -
                                                   qa->getTotalNumHs(true)));
    }  // end of adjust heavy degree
    if (params.adjustRingCount &&
        !((params.adjustRingCountFlags & ADJUST_IGNORECHAINS) && !nRings) &&
        !((params.adjustRingCountFlags & ADJUST_IGNORERINGS) && nRings) &&
        !((params.adjustRingCountFlags & ADJUST_IGNOREDUMMIES) && !atomicNum) &&
        !((params.adjustRingCountFlags & ADJUST_IGNORENONDUMMIES) &&
          atomicNum) &&
        !((params.adjustRingCountFlags & ADJUST_IGNOREMAPPED) &&
          isMapped(at))) {
      QueryAtom *qa;
      if (!at->hasQuery()) {
        QueryAtom atQueryAtom(*at);
        constexpr bool updateLabel = false;
        constexpr bool preserveProps = true;
        mol.replaceAtom(i, &atQueryAtom, updateLabel, preserveProps);
        qa = static_cast<QueryAtom *>(mol.getAtomWithIdx(i));
        at = static_cast<Atom *>(qa);
      } else {
        qa = static_cast<QueryAtom *>(at);
      }
      qa->expandQuery(makeAtomInNRingsQuery(nRings));
    }  // end of adjust ring count
    if (params.adjustRingChain &&
        !((params.adjustRingChainFlags & ADJUST_IGNORECHAINS) && !nRings) &&
        !((params.adjustRingChainFlags & ADJUST_IGNORERINGS) && nRings) &&
        !((params.adjustRingChainFlags & ADJUST_IGNOREDUMMIES) && !atomicNum) &&
        !((params.adjustRingChainFlags & ADJUST_IGNORENONDUMMIES) &&
          atomicNum) &&
        !((params.adjustRingChainFlags & ADJUST_IGNOREMAPPED) &&
          isMapped(at))) {
      QueryAtom *qa;
      if (!at->hasQuery()) {
        QueryAtom atQueryAtom(*at);
        constexpr bool updateLabel = false;
        constexpr bool preserveProps = true;
        mol.replaceAtom(i, &atQueryAtom, updateLabel, preserveProps);
        qa = static_cast<QueryAtom *>(mol.getAtomWithIdx(i));
        at = static_cast<Atom *>(qa);
      } else {
        qa = static_cast<QueryAtom *>(at);
      }
      ATOM_EQUALS_QUERY *nq = makeAtomInRingQuery();
      if (!nRings) {
        nq->setNegation(true);
      }
      qa->expandQuery(nq);
    }  // end of adjust ring chain
  }  // end of loop over atoms
  if (params.useStereoCareForBonds) {
    for (auto bnd : mol.bonds()) {
      if (bnd->getBondType() == Bond::BondType::DOUBLE) {
        if (bnd->getStereo() > Bond::BondStereo::STEREOANY) {
          bool preserve = false;
          int val = 0;
          // is stereoCare set on the bond or both atoms?
          if (bnd->getPropIfPresent(common_properties::molStereoCare, val) &&
              val) {
            preserve = true;
          }
          if (!preserve) {
            bnd->setStereo(Bond::BondStereo::STEREONONE);
          }
        }
      }
    }
  }
  if (params.setMDLFiveRingAromaticity) {
    setMDLAromaticity(mol);
  }
  if (params.adjustConjugatedFiveRings) {
    adjustConjugatedFiveRings(mol);
  }
  if (params.adjustSingleBondsToDegreeOneNeighbors ||
      params.adjustSingleBondsBetweenAromaticAtoms) {
    adjustSingleBondsFromAromaticAtoms(
        mol, params.adjustSingleBondsToDegreeOneNeighbors,
        params.adjustSingleBondsBetweenAromaticAtoms);
  }
  if (params.setMDLFiveRingAromaticity) {
    for (auto atom : mol.atoms()) {
      atom->clearProp(conjugatedOrAromatic);
    }
  }
}
}  // namespace MolOps
}  // namespace RDKit
