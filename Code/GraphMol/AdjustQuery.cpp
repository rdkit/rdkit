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

void adjustConjugatedFiveRings(RWMol &mol) {
  std::vector<Bond::BondType> bondTypesToModify = {
      Bond::BondType::SINGLE, Bond::BondType::DOUBLE, Bond::BondType::AROMATIC};
  if (!mol.getRingInfo()->isInitialized()) {
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
    // FIX: replace this with the SingleOrDoubleOrAromatic query once the
    // tautomer PR is merged
    QueryBond qb;
    qb.setQuery(makeBondOrderEqualsQuery(Bond::BondType::SINGLE));
    qb.expandQuery(makeBondOrderEqualsQuery(Bond::BondType::DOUBLE),
                   Queries::COMPOSITE_OR);
    qb.expandQuery(makeBondOrderEqualsQuery(Bond::BondType::AROMATIC),
                   Queries::COMPOSITE_OR);
    for (auto bi : ring) {
      const auto bond = mol.getBondWithIdx(bi);
      if (!bond->hasQuery() &&
          std::find(bondTypesToModify.begin(), bondTypesToModify.end(),
                    bond->getBondType()) != bondTypesToModify.end()) {
        mol.replaceBond(bi, &qb);
      }
    }
  }
}
void adjustSingleBondsFromAromaticAtoms(RWMol &mol, bool toDegreeOneNeighbors,
                                        bool betweenAromaticAtoms) {
  QueryBond qb;
  qb.setQuery(makeBondOrderEqualsQuery(Bond::BondType::SINGLE));
  qb.expandQuery(makeBondOrderEqualsQuery(Bond::BondType::AROMATIC),
                 Queries::COMPOSITE_OR);
  if (!mol.getRingInfo()->isInitialized()) {
    MolOps::symmetrizeSSSR(mol);
  }
  for (auto bond : mol.bonds()) {
    const auto bAt = bond->getBeginAtom();
    const auto eAt = bond->getEndAtom();
    if (!bond->hasQuery() && bond->getBondType() == Bond::BondType::SINGLE &&
        (bAt->getIsAromatic() || eAt->getIsAromatic()) &&
        !mol.getRingInfo()->numBondRings(bond->getIdx())) {
      if (toDegreeOneNeighbors &&
          (bAt->getIsAromatic() ^ eAt->getIsAromatic())) {
        if ((bAt->getIsAromatic() && eAt->getDegree() == 1) ||
            (eAt->getIsAromatic() && bAt->getDegree() == 1)) {
          mol.replaceBond(bond->getIdx(), &qb);
        }
      } else if (betweenAromaticAtoms && bAt->getIsAromatic() &&
                 eAt->getIsAromatic()) {
        mol.replaceBond(bond->getIdx(), &qb);
      }
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
    p.adjustRingCountFlags = parseWhichString(which);
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
  const RingInfo *ringInfo = mol.getRingInfo();

  if (params.aromatizeIfPossible) {
    unsigned int failed;
    sanitizeMol(mol, failed, SANITIZE_SYMMRINGS | SANITIZE_SETAROMATICITY);
  } else {
    if (!ringInfo->isInitialized()) {
      MolOps::symmetrizeSSSR(mol);
    }
  }

  if (params.makeAtomsGeneric) {
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      if (!((params.makeAtomsGenericFlags & ADJUST_IGNORECHAINS) &&
            !ringInfo->numAtomRings(i)) &&
          !((params.makeAtomsGenericFlags & ADJUST_IGNORERINGS) &&
            ringInfo->numAtomRings(i)) &&
          !((params.adjustDegreeFlags & ADJUST_IGNOREMAPPED) &&
            isMapped(mol.getAtomWithIdx(i)))) {
        auto *qa = new QueryAtom();
        qa->setQuery(makeAtomNullQuery());
        const bool updateLabel = false;
        const bool preserveProps = true;
        mol.replaceAtom(i, qa, updateLabel, preserveProps);
        delete qa;
      }
    }
  }  // end of makeAtomsGeneric
  if (params.makeBondsGeneric) {
    for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
      if (!((params.makeBondsGenericFlags & ADJUST_IGNORECHAINS) &&
            !ringInfo->numBondRings(i)) &&
          !((params.makeBondsGenericFlags & ADJUST_IGNORERINGS) &&
            ringInfo->numBondRings(i))) {
        auto *qb = new QueryBond();
        qb->setQuery(makeBondNullQuery());
        const bool preserveProps = true;
        mol.replaceBond(i, qb, preserveProps);
        delete qb;
      }
    }
  }  // end of makeBondsGeneric
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    Atom *at = mol.getAtomWithIdx(i);
    // pull properties we need from the atom here, once we
    // create a query atom they may no longer be valid.
    unsigned int nRings = ringInfo->numAtomRings(i);
    int atomicNum = at->getAtomicNum();
    if (params.makeDummiesQueries && atomicNum == 0 && !at->hasQuery() &&
        !at->getIsotope()) {
      auto *qa = new QueryAtom();
      qa->setQuery(makeAtomNullQuery());
      const bool updateLabel = false;
      const bool preserveProps = true;
      mol.replaceAtom(i, qa, updateLabel, preserveProps);
      delete qa;
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
        qa = new QueryAtom(*at);
        const bool updateLabel = false;
        const bool preserveProps = true;
        mol.replaceAtom(i, qa, updateLabel, preserveProps);
        delete qa;
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
        qa = new QueryAtom(*at);
        const bool updateLabel = false;
        const bool preserveProps = true;
        mol.replaceAtom(i, qa, updateLabel, preserveProps);
        delete qa;
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
        qa = new QueryAtom(*at);
        const bool updateLabel = false;
        const bool preserveProps = true;
        mol.replaceAtom(i, qa, updateLabel, preserveProps);
        delete qa;
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
        qa = new QueryAtom(*at);
        const bool updateLabel = false;
        const bool preserveProps = true;
        mol.replaceAtom(i, qa, updateLabel, preserveProps);
        delete qa;
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
  }    // end of loop over atoms
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
  if (params.adjustConjugatedFiveRings) {
    adjustConjugatedFiveRings(mol);
  }
  if (params.adjustSingleBondsToDegreeOneNeighbors ||
      params.adjustSingleBondsBetweenAromaticAtoms) {
    adjustSingleBondsFromAromaticAtoms(
        mol, params.adjustSingleBondsToDegreeOneNeighbors,
        params.adjustSingleBondsBetweenAromaticAtoms);
  }
}
}  // namespace MolOps
}  // namespace RDKit
