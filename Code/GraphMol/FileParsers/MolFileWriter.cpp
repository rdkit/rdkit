//
//  Copyright (C) 2003-2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  23/12/2013:
//     V3000 mol block writer contributed by Jan Holst Jensen
//
#include "FileParsers.h"
#include "FileParserUtils.h"
#include "MolSGroupWriting.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/Atropisomers.h>
#include <RDGeneral/Ranking.h>
#include <RDGeneral/LocaleSwitcher.h>

#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>

#include <boost/format.hpp>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/GenericGroups/GenericGroups.h>

#include <boost/algorithm/string.hpp>

using namespace RDKit::SGroupWriting;

namespace RDKit {

//*************************************
//
// Every effort has been made to adhere to MDL's standard
// for mol files
//
//*************************************

namespace {

int getQueryBondTopology(const Bond *bond) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(bond->hasQuery(), "no query");
  int res = 0;
  Bond::QUERYBOND_QUERY *qry = bond->getQuery();
  // start by catching combined bond order + bond topology queries

  if (qry->getDescription() == "BondAnd" && !qry->getNegation() &&
      qry->endChildren() - qry->beginChildren() == 2) {
    auto child1 = qry->beginChildren();
    auto child2 = child1 + 1;
    if (((*child1)->getDescription() == "BondInRing") !=
        ((*child2)->getDescription() == "BondInRing")) {
      if ((*child1)->getDescription() != "BondInRing") {
        std::swap(child1, child2);
      }
      qry = child1->get();
    }
  }
  if (qry->getDescription() == "BondInRing") {
    if (qry->getNegation()) {
      res = 2;
    } else {
      res = 1;
    }
  }
  return res;
}

// returns 0 if there's a basic bond-order query
int getQueryBondSymbol(const Bond *bond) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(bond->hasQuery(), "no query");
  int res = 8;

  Bond::QUERYBOND_QUERY *qry = bond->getQuery();
  if (qry->getDescription() == "BondOrder" || getQueryBondTopology(bond)) {
    // trap the simple bond-order query
    res = 0;
  } else {
    // start by catching combined bond order + bond topology queries
    if (qry->getDescription() == "BondAnd" && !qry->getNegation() &&
        qry->endChildren() - qry->beginChildren() == 2) {
      auto child1 = qry->beginChildren();
      auto child2 = child1 + 1;
      if ((*child2)->getDescription() == "BondInRing") {
        qry = child1->get();
      } else if ((*child1)->getDescription() == "BondInRing") {
        qry = child2->get();
      }
    }
    if (qry->getDescription() == "BondOr" && !qry->getNegation()) {
      if (qry->endChildren() - qry->beginChildren() == 2) {
        auto child1 = qry->beginChildren();
        auto child2 = child1 + 1;
        if ((*child1)->getDescription() == "BondOrder" &&
            !(*child1)->getNegation() &&
            (*child2)->getDescription() == "BondOrder" &&
            !(*child2)->getNegation()) {
          // ok, it's a bond query we have a chance of dealing with
          int t1 = static_cast<BOND_EQUALS_QUERY *>(child1->get())->getVal();
          int t2 = static_cast<BOND_EQUALS_QUERY *>(child2->get())->getVal();
          if (t1 > t2) {
            std::swap(t1, t2);
          }
          if (t1 == Bond::SINGLE && t2 == Bond::DOUBLE) {
            res = 5;
          } else if (t1 == Bond::SINGLE && t2 == Bond::AROMATIC) {
            res = 6;
          } else if (t1 == Bond::DOUBLE && t2 == Bond::AROMATIC) {
            res = 7;
          }
        }
      }
    } else if (qry->getDescription() == "SingleOrAromaticBond" &&
               !qry->getNegation()) {
      res = 6;
    } else if (qry->getDescription() == "SingleOrDoubleBond" &&
               !qry->getNegation()) {
      res = 5;
    } else if (qry->getDescription() == "DoubleOrAromaticBond" &&
               !qry->getNegation()) {
      res = 7;
    }
  }
  return res;
}

bool isAtomRGroup(const Atom &atom) {
  return atom.getAtomicNum() == 0 &&
         atom.hasProp(common_properties::_MolFileRLabel);
}
}  // namespace

const std::string GetMolFileChargeInfo(const RWMol &mol) {
  std::stringstream res;
  std::stringstream chgss;
  std::stringstream radss;
  std::stringstream massdiffss;
  unsigned int nChgs = 0;
  unsigned int nRads = 0;
  unsigned int nMassDiffs = 0;
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    const Atom *atom = *atomIt;
    if (atom->getFormalCharge() != 0) {
      ++nChgs;
      chgss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
                   atom->getFormalCharge();
      if (nChgs == 8) {
        res << boost::format("M  CHG%3d") % nChgs << chgss.str() << "\n";
        chgss.str("");
        nChgs = 0;
      }
    }
    unsigned int nRadEs = atom->getNumRadicalElectrons();
    if (nRadEs != 0 && atom->getTotalDegree() != 0) {
      ++nRads;
      if (nRadEs % 2) {
        nRadEs = 2;
      } else {
        nRadEs = 3;  // we use triplets, not singlets:
      }
      radss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % nRadEs;
      if (nRads == 8) {
        res << boost::format("M  RAD%3d") % nRads << radss.str() << "\n";
        radss.str("");
        nRads = 0;
      }
    }
    if (!isAtomRGroup(*atom)) {
      int isotope = atom->getIsotope();
      if (isotope != 0) {
        ++nMassDiffs;
        massdiffss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
                          isotope;
        if (nMassDiffs == 8) {
          res << boost::format("M  ISO%3d") % nMassDiffs << massdiffss.str()
              << "\n";
          massdiffss.str("");
          nMassDiffs = 0;
        }
      }
    }
  }
  if (nChgs) {
    res << boost::format("M  CHG%3d") % nChgs << chgss.str() << "\n";
  }
  if (nRads) {
    res << boost::format("M  RAD%3d") % nRads << radss.str() << "\n";
  }
  if (nMassDiffs) {
    res << boost::format("M  ISO%3d") % nMassDiffs << massdiffss.str() << "\n";
  }
  return res.str();
}

bool hasComplexQuery(const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  bool res = false;
  if (atom->hasQuery()) {
    res = true;
    // counter examples:
    //  1) atomic number
    //  2) the smarts parser inserts AtomAnd queries
    //     for "C" or "c":
    //
    std::string descr = atom->getQuery()->getDescription();
    if (descr == "AtomAtomicNum") {
      res = false;
    } else if (descr == "AtomAnd") {
      if ((*atom->getQuery()->beginChildren())->getDescription() ==
          "AtomAtomicNum") {
        res = false;
      }
    }
  }
  return res;
}

const std::string GetMolFileQueryInfo(
    const RWMol &mol, const boost::dynamic_bitset<> &queryListAtoms) {
  std::stringstream ss;
  boost::dynamic_bitset<> listQs(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    if (isAtomListQuery(atom) && !queryListAtoms[atom->getIdx()]) {
      listQs.set(atom->getIdx());
    }
  }
  for (const auto atom : mol.atoms()) {
    bool wrote_query = false;
    if (!listQs[atom->getIdx()] && !queryListAtoms[atom->getIdx()] &&
        hasComplexQuery(atom)) {
      std::string sma =
          SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(atom));
      ss << "V  " << std::setw(3) << atom->getIdx() + 1 << " " << sma << "\n";
      wrote_query = true;
    }
    std::string molFileValue;
    if (!wrote_query &&
        atom->getPropIfPresent(common_properties::molFileValue, molFileValue)) {
      ss << "V  " << std::setw(3) << atom->getIdx() + 1 << " " << molFileValue
         << "\n";
    }
  }
  for (const auto atom : mol.atoms()) {
    if (listQs[atom->getIdx()]) {
      INT_VECT vals;
      getAtomListQueryVals(atom->getQuery(), vals);
      ss << "M  ALS " << std::setw(3) << atom->getIdx() + 1 << " ";
      ss << std::setw(2) << vals.size();
      if (atom->getQuery()->getNegation()) {
        ss << " T ";
      } else {
        ss << " F ";
      }
      for (auto val : vals) {
        ss << std::setw(4) << std::left
           << (PeriodicTable::getTable()->getElementSymbol(val));
      }
      ss << "\n";
    }
  }
  return ss.str();
}

const std::string GetMolFileRGroupInfo(const RWMol &mol) {
  std::stringstream ss;
  unsigned int nEntries = 0;
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    unsigned int lbl;
    if ((*atomIt)->getPropIfPresent(common_properties::_MolFileRLabel, lbl)) {
      ss << " " << std::setw(3) << (*atomIt)->getIdx() + 1 << " "
         << std::setw(3) << lbl;
      ++nEntries;
    }
  }
  std::stringstream ss2;
  if (nEntries) {
    ss2 << "M  RGP" << std::setw(3) << nEntries << ss.str() << "\n";
  }
  return ss2.str();
}

const std::string GetMolFileAliasInfo(const RWMol &mol) {
  std::stringstream ss;
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    std::string lbl;
    if ((*atomIt)->getPropIfPresent(common_properties::molFileAlias, lbl)) {
      if (!lbl.empty()) {
        ss << "A  " << std::setw(3) << (*atomIt)->getIdx() + 1 << "\n"
           << lbl << "\n";
      }
    }
  }
  return ss.str();
}

const std::string GetMolFilePXAInfo(const RWMol &mol) {
  std::string res;
  for (const auto atom : mol.atoms()) {
    if (atom->hasProp("_MolFile_PXA")) {
      res +=
          boost::str(boost::format("M  PXA % 3d%s\n") % (atom->getIdx() + 1) %
                     atom->getProp<std::string>("_MolFile_PXA"));
    }
  }
  return res;
}
const std::string GetMolFileZBOInfo(const RWMol &mol) {
  std::stringstream res;
  std::stringstream ss;
  unsigned int nEntries = 0;
  boost::dynamic_bitset<> atomsAffected(mol.getNumAtoms(), 0);
  for (ROMol::ConstBondIterator bondIt = mol.beginBonds();
       bondIt != mol.endBonds(); ++bondIt) {
    if ((*bondIt)->getBondType() == Bond::ZERO) {
      ++nEntries;
      ss << " " << std::setw(3) << (*bondIt)->getIdx() + 1 << " "
         << std::setw(3) << 0;
      if (nEntries == 8) {
        res << "M  ZBO" << std::setw(3) << nEntries << ss.str() << "\n";
        nEntries = 0;
        ss.str("");
      }
      atomsAffected[(*bondIt)->getBeginAtomIdx()] = 1;
      atomsAffected[(*bondIt)->getEndAtomIdx()] = 1;
    }
  }
  if (nEntries) {
    res << "M  ZBO" << std::setw(3) << nEntries << ss.str() << "\n";
  }
  if (atomsAffected.count()) {
    std::stringstream hydss;
    unsigned int nhyd = 0;
    std::stringstream zchss;
    unsigned int nzch = 0;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      if (!atomsAffected[i]) {
        continue;
      }
      const Atom *atom = mol.getAtomWithIdx(i);
      nhyd++;
      hydss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
                   atom->getTotalNumHs();
      if (nhyd == 8) {
        res << boost::format("M  HYD%3d") % nhyd << hydss.str() << "\n";
        hydss.str("");
        nhyd = 0;
      }
      if (atom->getFormalCharge()) {
        nzch++;
        zchss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
                     atom->getFormalCharge();
        if (nzch == 8) {
          res << boost::format("M  ZCH%3d") % nzch << zchss.str() << "\n";
          zchss.str("");
          nzch = 0;
        }
      }
    }
    if (nhyd) {
      res << boost::format("M  HYD%3d") % nhyd << hydss.str() << "\n";
    }
    if (nzch) {
      res << boost::format("M  ZCH%3d") % nzch << zchss.str() << "\n";
    }
  }
  return res.str();
}

const std::string AtomGetMolFileSymbol(
    const Atom *atom, bool padWithSpaces,
    boost::dynamic_bitset<> &queryListAtoms) {
  PRECONDITION(atom, "");

  std::string res;
  if (atom->hasProp(common_properties::_MolFileRLabel)) {
    res = "R#";
    //    } else if(!atom->hasQuery() && atom->getAtomicNum()){
  } else if (atom->getAtomicNum()) {
    res = atom->getSymbol();
  } else {
    if (!atom->hasProp(common_properties::dummyLabel)) {
      if (atom->hasQuery() &&
          (atom->getQuery()->getTypeLabel() == "A" ||
           (atom->getQuery()->getNegation() &&
            atom->getQuery()->getDescription() == "AtomAtomicNum" &&
            static_cast<ATOM_EQUALS_QUERY *>(atom->getQuery())->getVal() ==
                1))) {
        res = "A";
        queryListAtoms.set(atom->getIdx());
      } else if (atom->hasQuery() &&
                 (atom->getQuery()->getTypeLabel() == "Q" ||
                  (atom->getQuery()->getNegation() &&
                   atom->getQuery()->getDescription() == "AtomOr" &&
                   atom->getQuery()->endChildren() -
                           atom->getQuery()->beginChildren() ==
                       2 &&
                   (*atom->getQuery()->beginChildren())->getDescription() ==
                       "AtomAtomicNum" &&
                   static_cast<ATOM_EQUALS_QUERY *>(
                       (*atom->getQuery()->beginChildren()).get())
                           ->getVal() == 6 &&
                   (*++(atom->getQuery()->beginChildren()))->getDescription() ==
                       "AtomAtomicNum" &&
                   static_cast<ATOM_EQUALS_QUERY *>(
                       (*++(atom->getQuery()->beginChildren())).get())
                           ->getVal() == 1))) {
        res = "Q";
        queryListAtoms.set(atom->getIdx());
      } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "X") {
        res = "X";
        queryListAtoms.set(atom->getIdx());
      } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "M") {
        res = "M";
        queryListAtoms.set(atom->getIdx());
      } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "AH") {
        res = "AH";
        queryListAtoms.set(atom->getIdx());
      } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "QH") {
        res = "QH";
        queryListAtoms.set(atom->getIdx());
      } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "XH") {
        res = "XH";
        queryListAtoms.set(atom->getIdx());
      } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "MH") {
        res = "MH";
        queryListAtoms.set(atom->getIdx());
      } else if (hasComplexQuery(atom)) {
        if (isAtomListQuery(atom)) {
          res = "L";
        } else {
          res = "*";
        }
      } else {
        res = "R";
      }
    } else {
      std::string symb;
      atom->getProp(common_properties::dummyLabel, symb);
      if (symb == "*") {
        res = "R";
      } else if (symb == "X") {
        res = "R";
      } else if (symb == "Xa") {
        res = "R1";
      } else if (symb == "Xb") {
        res = "R2";
      } else if (symb == "Xc") {
        res = "R3";
      } else if (symb == "Xd") {
        res = "R4";
      } else if (symb == "Xf") {
        res = "R5";
      } else if (symb == "Xg") {
        res = "R6";
      } else if (symb == "Xh") {
        res = "R7";
      } else if (symb == "Xi") {
        res = "R8";
      } else if (symb == "Xj") {
        res = "R9";
      } else {
        res = symb;
      }
    }
  }
  // pad the end with spaces
  if (padWithSpaces) {
    while (res.size() < 3) {
      res += " ";
    }
  }
  return res;
}

namespace {
unsigned int getAtomParityFlag(const Atom *atom, const Conformer *conf) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(conf, "bad conformer");
  if (!conf->is3D() ||
      !(atom->getDegree() >= 3 && atom->getTotalDegree() == 4)) {
    return 0;
  }

  const ROMol &mol = atom->getOwningMol();
  RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
  std::vector<std::pair<unsigned int, RDGeom::Point3D>> vs;
  ROMol::ADJ_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
  while (nbrIdx != endNbrs) {
    const Atom *at = mol.getAtomWithIdx(*nbrIdx);
    unsigned int idx = at->getIdx();
    RDGeom::Point3D v = conf->getAtomPos(idx);
    v -= pos;
    if (at->getAtomicNum() == 1) {
      idx += mol.getNumAtoms();
    }
    vs.emplace_back(idx, v);
    ++nbrIdx;
  }
  std::sort(vs.begin(), vs.end(), Rankers::pairLess);
  double vol;
  if (vs.size() == 4) {
    vol = vs[0].second.crossProduct(vs[1].second).dotProduct(vs[3].second);
  } else {
    vol = -vs[0].second.crossProduct(vs[1].second).dotProduct(vs[2].second);
  }
  if (vol < 0) {
    return 2;
  } else if (vol > 0) {
    return 1;
  }
  return 0;
}
}  // namespace

bool hasNonDefaultValence(const Atom *atom) {
  if (atom->getNumRadicalElectrons() != 0) {
    return true;
  }
  // for queries and atoms which don't have computed properties, the answer is
  // always no:
  if (atom->hasQuery() || atom->needsUpdatePropertyCache()) {
    return false;
  }

  if (atom->getAtomicNum() == 1 ||
      SmilesWrite ::inOrganicSubset(atom->getAtomicNum())) {
    // for the ones we "know", we may have to specify the valence if it's
    // not the default value
    return atom->getNoImplicit() &&
           (atom->getExplicitValence() !=
            PeriodicTable::getTable()->getDefaultValence(atom->getAtomicNum()));
  }
  return true;
}

void GetMolFileAtomProperties(const Atom *atom, const Conformer *conf,
                              int &totValence, int &atomMapNumber,
                              unsigned int &parityFlag, double &x, double &y,
                              double &z) {
  PRECONDITION(atom, "");
  totValence = 0;
  atomMapNumber = 0;
  parityFlag = 0;
  x = y = z = 0.0;

  if (!atom->getPropIfPresent(common_properties::molAtomMapNumber,
                              atomMapNumber)) {
    // XXX FIX ME->should we fail here? previously we would not assign
    // the atomMapNumber if it didn't exist which could result in garbage
    //  values.
    atomMapNumber = 0;
  }

  if (conf) {
    const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    x = pos.x;
    y = pos.y;
    z = pos.z;
    if (conf->is3D() && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        atom->getChiralTag() != Atom::CHI_OTHER && atom->getDegree() >= 3 &&
        atom->getTotalDegree() == 4) {
      parityFlag = getAtomParityFlag(atom, conf);
    }
  }
  if (hasNonDefaultValence(atom)) {
    if (atom->getTotalDegree() == 0) {
      // Specify zero valence for elements/metals without neighbors
      // or hydrogens (degree 0) instead of writing them as radicals.
      totValence = 15;
    } else {
      // write the total valence for other atoms
      totValence = atom->getTotalValence() % 15;
    }
  }
}

const std::string GetMolFileAtomLine(const Atom *atom, const Conformer *conf,
                                     boost::dynamic_bitset<> &queryListAtoms) {
  PRECONDITION(atom, "");
  std::string res;
  int totValence, atomMapNumber;
  unsigned int parityFlag;
  double x, y, z;
  GetMolFileAtomProperties(atom, conf, totValence, atomMapNumber, parityFlag, x,
                           y, z);

  int massDiff, chg, stereoCare, hCount, rxnComponentType, rxnComponentNumber,
      inversionFlag, exactChangeFlag;
  massDiff = 0;
  chg = 0;
  stereoCare = 0;
  hCount = 0;
  rxnComponentType = 0;
  rxnComponentNumber = 0;
  inversionFlag = 0;
  exactChangeFlag = 0;

  atom->getPropIfPresent(common_properties::molRxnRole, rxnComponentType);
  atom->getPropIfPresent(common_properties::molRxnComponent,
                         rxnComponentNumber);

  std::string symbol = AtomGetMolFileSymbol(atom, true, queryListAtoms);
#if 0
  const boost::format fmter(
      "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d");
  std::stringstream ss;
  ss << boost::format(fmter) % x % y % z % symbol.c_str() % massDiff % chg %
            parityFlag % hCount % stereoCare % totValence % rxnComponentType %
            rxnComponentNumber % atomMapNumber % inversionFlag %
            exactChangeFlag;
  res += ss.str();
#else
  // it feels ugly to use snprintf instead of boost::format, but at least of the
  // time of this writing (with boost 1.55), the snprintf version runs in 20% of
  // the time.
  char dest[128];
#ifndef _MSC_VER
  snprintf(dest, 128,
           "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d", x, y,
           z, symbol.c_str(), massDiff, chg, parityFlag, hCount, stereoCare,
           totValence, rxnComponentType, rxnComponentNumber, atomMapNumber,
           inversionFlag, exactChangeFlag);
#else
  // ok, technically we should be being more careful about this, but given that
  // the format string makes it impossible for this to overflow, I think we're
  // safe. I just used the snprintf above to prevent linters from complaining
  // about use of sprintf
  sprintf_s(dest, 128,
            "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d", x, y,
            z, symbol.c_str(), massDiff, chg, parityFlag, hCount, stereoCare,
            totValence, rxnComponentType, rxnComponentNumber, atomMapNumber,
            inversionFlag, exactChangeFlag);

#endif
  res += dest;
#endif
  return res;
};

int BondGetMolFileSymbol(const Bond *bond) {
  PRECONDITION(bond, "");
  // FIX: should eventually recognize queries
  int res = 0;
  if (bond->hasQuery()) {
    res = getQueryBondSymbol(bond);
  }
  if (!res) {
    switch (bond->getBondType()) {
      case Bond::SINGLE:
        if (bond->getIsAromatic()) {
          res = 4;
        } else {
          res = 1;
        }
        break;
      case Bond::DOUBLE:
        if (bond->getIsAromatic()) {
          res = 4;
        } else {
          res = 2;
        }
        break;
      case Bond::TRIPLE:
        res = 3;
        break;
      case Bond::AROMATIC:
        res = 4;
        break;
      case Bond::ZERO:
        res = 1;
        break;
      case Bond::DATIVE:
        // extension
        res = 9;
        break;
      default:
        break;
    }
  }
  return res;
  // return res.c_str();
}

const std::string GetMolFileBondLine(
    const Bond *bond,
    const std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> &wedgeBonds,
    const Conformer *conf, bool wasAromatic) {
  PRECONDITION(bond, "");

  int dirCode = 0;
  bool reverse = false;
  RDKit::Chirality::GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dirCode,
                                             reverse);
  // do not cross bonds which were aromatic before kekulization
  if (wasAromatic && dirCode == 3) {
    dirCode = 0;
  }
  int symbol = BondGetMolFileSymbol(bond);

  std::stringstream ss;
  if (reverse) {
    // switch the begin and end atoms on the bond line
    ss << std::setw(3) << bond->getEndAtomIdx() + 1;
    ss << std::setw(3) << bond->getBeginAtomIdx() + 1;
  } else {
    ss << std::setw(3) << bond->getBeginAtomIdx() + 1;
    ss << std::setw(3) << bond->getEndAtomIdx() + 1;
  }
  ss << std::setw(3) << symbol;
  ss << " " << std::setw(2) << dirCode;

  if (bond->hasQuery()) {
    int topol = getQueryBondTopology(bond);
    if (topol) {
      ss << " " << std::setw(2) << 0 << " " << std::setw(2) << topol;
    }
  }

  return ss.str();
}

const std::string GetV3000MolFileAtomLine(
    const Atom *atom, const Conformer *conf,
    boost::dynamic_bitset<> &queryListAtoms, unsigned int precision) {
  PRECONDITION(atom, "");
  int totValence, atomMapNumber;
  unsigned int parityFlag;
  double x, y, z;
  GetMolFileAtomProperties(atom, conf, totValence, atomMapNumber, parityFlag, x,
                           y, z);

  std::stringstream ss;
  ss << "M  V30 " << atom->getIdx() + 1;

  std::string symbol = AtomGetMolFileSymbol(atom, false, queryListAtoms);
  if (!isAtomListQuery(atom) || queryListAtoms[atom->getIdx()]) {
    ss << " " << symbol;
  } else {
    INT_VECT vals;
    getAtomListQueryVals(atom->getQuery(), vals);
    if (atom->getQuery()->getNegation()) {
      ss << " "
         << "\"NOT";
    }
    ss << " [";
    for (unsigned int i = 0; i < vals.size(); ++i) {
      if (i != 0) {
        ss << ",";
      }
      ss << PeriodicTable::getTable()->getElementSymbol(vals[i]);
    }
    ss << "]";
    if (atom->getQuery()->getNegation()) {
      ss << "\"";
    }
  }

  std::streamsize currentPrecision = ss.precision();
  ss << std::fixed;
  ss << std::setprecision(precision);
  ss << " " << x << " " << y << " " << z;
  ss << std::setprecision(currentPrecision);
  ss << std::defaultfloat;
  ss << " " << atomMapNumber;

  // Extra atom properties.
  int chg = atom->getFormalCharge();
  int isotope = atom->getIsotope();
  if (parityFlag != 0) {
    ss << " CFG=" << parityFlag;
  }
  if (chg != 0) {
    ss << " CHG=" << chg;
  }
  if (isotope != 0 && !isAtomRGroup(*atom)) {
    // the documentation for V3000 CTABs says that this should contain the
    // "absolute atomic weight" (whatever that means).
    // Online examples seem to have integer (isotope) values and Marvin won't
    // even read something that has a float.
    // We'll go with the int.
    int mass = static_cast<int>(std::round(atom->getMass()));
    // dummies may have an isotope set but they always have a mass of zero:
    if (!mass) {
      mass = isotope;
    }
    ss << " MASS=" << mass;
  }

  unsigned int nRadEs = atom->getNumRadicalElectrons();
  if (nRadEs != 0 && atom->getTotalDegree() != 0) {
    if (nRadEs % 2) {
      nRadEs = 2;
    } else {
      nRadEs = 3;  // we use triplets, not singlets:
    }
    ss << " RAD=" << nRadEs;
  }

  if (totValence != 0) {
    if (totValence == 15) {
      ss << " VAL=-1";
    } else {
      ss << " VAL=" << totValence;
    }
  }
  if (symbol == "R#") {
    unsigned int rLabel = 1;
    atom->getPropIfPresent(common_properties::_MolFileRLabel, rLabel);
    ss << " RGROUPS=(1 " << rLabel << ")";
  }

  {
    int iprop;
    if (atom->getPropIfPresent(common_properties::molAttachOrder, iprop) &&
        iprop) {
      ss << " ATTCHORD=" << iprop;
    }
    if (atom->getPropIfPresent(common_properties::molAttachPoint, iprop) &&
        iprop) {
      ss << " ATTCHPT=" << iprop;
    }
    if (atom->getPropIfPresent(common_properties::molAtomSeqId, iprop) &&
        iprop) {
      ss << " SEQID=" << iprop;
    }
    if (atom->getPropIfPresent(common_properties::molRxnExactChange, iprop) &&
        iprop) {
      ss << " EXACHG=" << iprop;
    }
    if (atom->getPropIfPresent(common_properties::molInversionFlag, iprop) &&
        iprop) {
      if (iprop == 1 || iprop == 2) {
        ss << " INVRET=" << iprop;
      }
    }
    if (atom->getPropIfPresent(common_properties::molStereoCare, iprop) &&
        iprop) {
      ss << " STBOX=" << iprop;
    }
    if (atom->getPropIfPresent(common_properties::molSubstCount, iprop) &&
        iprop) {
      ss << " SUBST=" << iprop;
    }
  }
  {
    std::string sprop;
    if (atom->getPropIfPresent(common_properties::molAtomClass, sprop)) {
      ss << " CLASS=" << sprop;
    }
  }
  // HCOUNT - *query* hydrogen count. Not written by this writer.

  return ss.str();
};

int GetV3000BondCode(const Bond *bond) {
  // JHJ: As far as I can tell, the V3000 bond codes are the same as for V2000.
  //      Except: The dative bond type is only supported in V3000.
  PRECONDITION(bond, "");
  int res = 0;
  // FIX: should eventually recognize queries
  if (bond->hasQuery()) {
    res = getQueryBondSymbol(bond);
  }
  if (!res) {
    switch (bond->getBondType()) {
      case Bond::SINGLE:
        if (bond->getIsAromatic()) {
          res = 4;
        } else {
          res = 1;
        }
        break;
      case Bond::DOUBLE:
        if (bond->getIsAromatic()) {
          res = 4;
        } else {
          res = 2;
        }
        break;
      case Bond::TRIPLE:
        res = 3;
        break;
      case Bond::AROMATIC:
        res = 4;
        break;
      case Bond::DATIVE:
        res = 9;
        break;
      case Bond::HYDROGEN:
        res = 10;
        break;
      case Bond::ZERO:
        res = 1;
        break;
      default:
        res = 0;
        break;
    }
  }
  return res;
}

int BondStereoCodeV2000ToV3000(int dirCode) {
  // The Any bond configuration (code 4 in v2000 ctabs) seems to be missing
  switch (dirCode) {
    case 0:
      return 0;
    case 1:
      return 1;  // V2000 Up       => Up.
    case 3:
      return 2;  // V2000 Unknown  => Either.
    case 4:
      return 2;  // V2000 Any      => Either.
    case 6:
      return 3;  // V2000 Down     => Down.
    default:
      return 0;
  }
}

namespace {
void createSMARTSQSubstanceGroups(ROMol &mol) {
  auto isRedundantQuery = [](const auto query) {
    if (query->getDescription() == "AtomAnd" &&
        (query->endChildren() - query->beginChildren() == 2) &&
        (*query->beginChildren())->getDescription() == "AtomAtomicNum" &&
        !(*query->beginChildren())->getNegation() &&
        !(*(query->beginChildren() + 1))->getNegation() &&
        ((*(query->beginChildren() + 1))->getDescription() == "AtomIsotope" ||
         (*(query->beginChildren() + 1))->getDescription() ==
             "AtomFormalCharge")) {
      return true;
    }
    return false;
  };
  for (const auto atom : mol.atoms()) {
    if (atom->hasQuery()) {
      std::string sma;

      if (!atom->getPropIfPresent(common_properties::MRV_SMA, sma) &&
          !isAtomListQuery(atom) &&
          atom->getQuery()->getDescription() != "AtomNull" &&
          // we may want to re-think this next one.
          // including AtomType queries will result in an entry
          // for every atom that comes from SMARTS, and I don't think
          // we want that.
          !boost::starts_with(atom->getQuery()->getDescription(), "AtomType") &&
          !boost::starts_with(atom->getQuery()->getDescription(),
                              "AtomAtomicNum") &&
          !isRedundantQuery(atom->getQuery())) {
        sma = SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(atom));
      }
      if (!sma.empty()) {
        SubstanceGroup sg(&mol, "DAT");
        sg.setProp("QUERYTYPE", "SMARTSQ");
        sg.setProp("QUERYOP", "=");
        std::vector<std::string> dataFields{sma};
        sg.setProp("DATAFIELDS", dataFields);
        sg.addAtomWithIdx(atom->getIdx());
        addSubstanceGroup(mol, sg);
      }
    }
  }
}

void createZBOSubstanceGroups(ROMol &mol) {
  SubstanceGroup bsg(&mol, "DAT");
  bsg.setProp("FIELDNAME", "ZBO");
  boost::dynamic_bitset<> atomsAffected(mol.getNumAtoms(), 0);
  for (const auto bond : mol.bonds()) {
    if (bond->getBondType() == Bond::ZERO) {
      bsg.addBondWithIdx(bond->getIdx());
      atomsAffected[bond->getBeginAtomIdx()] = 1;
      atomsAffected[bond->getEndAtomIdx()] = 1;
    }
  }
  if (atomsAffected.any()) {
    for (auto i = 0u; i < atomsAffected.size(); ++i) {
      if (atomsAffected[i]) {
        bsg.addAtomWithIdx(i);
      }
    }
    SubstanceGroup asg(&mol, "DAT");
    asg.setProp("FIELDNAME", "HYD");
    SubstanceGroup zsg(&mol, "DAT");
    zsg.setProp("FIELDNAME", "ZCH");
    std::string asgText;
    std::string zsgText;
    for (auto i = 0u; i < atomsAffected.size(); ++i) {
      if (atomsAffected[i]) {
        const Atom *atom = mol.getAtomWithIdx(i);
        asg.addAtomWithIdx(i);
        if (!asgText.empty()) {
          asgText += ";";
        }
        asgText += std::to_string(atom->getTotalNumHs());
        zsg.addAtomWithIdx(i);
        if (!zsgText.empty()) {
          zsgText += ";";
        }
        zsgText += std::to_string(atom->getFormalCharge());
      }
    }
    addSubstanceGroup(mol, bsg);

    std::vector<std::string> aDataFields{asgText};

    asg.setProp("DATAFIELDS", aDataFields);
    addSubstanceGroup(mol, asg);
    std::vector<std::string> zDataFields{zsgText};
    zsg.setProp("DATAFIELDS", zDataFields);
    addSubstanceGroup(mol, zsg);
  }
}
}  // namespace
namespace FileParserUtils {
void moveAdditionalPropertiesToSGroups(RWMol &mol) {
  GenericGroups::convertGenericQueriesToSubstanceGroups(mol);
  createSMARTSQSubstanceGroups(mol);
  createZBOSubstanceGroups(mol);
}
}  // namespace FileParserUtils
const std::string GetV3000MolFileBondLine(
    const Bond *bond,
    const std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> &wedgeBonds,
    const Conformer *conf, bool wasAromatic) {
  PRECONDITION(bond, "");

  int dirCode = 0;
  bool reverse = false;
  RDKit::Chirality::GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dirCode,
                                             reverse);
  // do not cross bonds which were aromatic before kekulization
  if (wasAromatic && dirCode == 3) {
    dirCode = 0;
  }

  std::stringstream ss;
  ss << "M  V30 " << bond->getIdx() + 1;
  ss << " " << GetV3000BondCode(bond);
  if (reverse) {
    // switch the begin and end atoms on the bond line
    ss << " " << bond->getEndAtomIdx() + 1;
    ss << " " << bond->getBeginAtomIdx() + 1;
  } else {
    ss << " " << bond->getBeginAtomIdx() + 1;
    ss << " " << bond->getEndAtomIdx() + 1;
  }
  if (dirCode != 0) {
    ss << " CFG=" << BondStereoCodeV2000ToV3000(dirCode);
  }
  if (bond->hasQuery()) {
    int topol = getQueryBondTopology(bond);
    if (topol) {
      ss << " TOPO=" << topol;
    }
  }

  {
    int iprop;
    if (bond->getPropIfPresent(common_properties::molReactStatus, iprop) &&
        iprop) {
      ss << " RXCTR=" << iprop;
    }
  }

  {
    std::string sprop;
    if (bond->getPropIfPresent(common_properties::molStereoCare, sprop) &&
        sprop != "0") {
      ss << " STBOX=" << sprop;
    }
    if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, sprop) &&
        sprop != "0") {
      ss << " ENDPTS=" << sprop;
    }
    if (bond->getPropIfPresent(common_properties::_MolFileBondAttach, sprop) &&
        sprop != "0") {
      ss << " ATTACH=" << sprop;
    }
  }

  return ss.str();
}

void appendEnhancedStereoGroups(
    std::string &res, const RWMol &tmol,
    std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> &wedgeBonds) {
  if (!tmol.getStereoGroups().empty()) {
    auto stereo_groups = tmol.getStereoGroups();
    assignStereoGroupIds(stereo_groups);
    res += "M  V30 BEGIN COLLECTION\n";
    std::string tmp;
    tmp.reserve(80);
    for (auto &&group : stereo_groups) {
      tmp += "M  V30 MDLV30/";
      switch (group.getGroupType()) {
        case RDKit::StereoGroupType::STEREO_ABSOLUTE:
          tmp += "STEABS";
          break;
        case RDKit::StereoGroupType::STEREO_OR:
          tmp += "STEREL";
          tmp += std::to_string(group.getWriteId());
          break;
        case RDKit::StereoGroupType::STEREO_AND:
          tmp += "STERAC";
          tmp += std::to_string(group.getWriteId());
          break;
      }
      tmp += " ATOMS=(";

      std::vector<unsigned int> atomIds;
      Atropisomers::getAllAtomIdsForStereoGroup(tmol, group, atomIds,
                                                wedgeBonds);

      tmp += std::to_string(atomIds.size());
      for (auto &&atom : atomIds) {
        tmp += ' ';
        // atoms are 1 indexed in molfiles
        auto idxStr = std::to_string(atom + 1);
        if (tmp.size() + idxStr.size() >= 78) {
          res += tmp + "-\n";
          tmp = "M  V30 ";
        }
        tmp += idxStr;
      }
      res += tmp + ")\n";
      tmp.clear();
    }
    res += tmp + "M  V30 END COLLECTION\n";
  }
}
namespace FileParserUtils {
std::string getV3000CTAB(const ROMol &tmol,
                         const boost::dynamic_bitset<> &wasAromatic, int confId,
                         unsigned int precision) {
  auto nAtoms = tmol.getNumAtoms();
  auto nBonds = tmol.getNumBonds();
  const auto &sgroups = getSubstanceGroups(tmol);
  auto nSGroups = sgroups.size();

  unsigned chiralFlag = 0;
  tmol.getPropIfPresent(common_properties::_MolFileChiralFlag, chiralFlag);

  const Conformer *conf = nullptr;
  if (confId >= 0 || tmol.getNumConformers()) {
    conf = &(tmol.getConformer(confId));
  }

  std::string res = "M  V30 BEGIN CTAB\n";
  std::stringstream ss;
  int num3DConstraints = 0;  //< not implemented
  ss << "M  V30 COUNTS " << nAtoms << " " << nBonds << " " << nSGroups << " "
     << num3DConstraints << " " << chiralFlag << "\n";

  res += ss.str();

  boost::dynamic_bitset<> queryListAtoms(tmol.getNumAtoms());
  res += "M  V30 BEGIN ATOM\n";
  for (ROMol::ConstAtomIterator atomIt = tmol.beginAtoms();
       atomIt != tmol.endAtoms(); ++atomIt) {
    res += GetV3000MolFileAtomLine(*atomIt, conf, queryListAtoms, precision);
    res += "\n";
  }
  res += "M  V30 END ATOM\n";

  auto wedgeBonds = Chirality::pickBondsToWedge(tmol, nullptr, conf);
  if (tmol.getNumBonds()) {
    res += "M  V30 BEGIN BOND\n";

    for (const auto bond : tmol.bonds()) {
      res += GetV3000MolFileBondLine(bond, wedgeBonds, conf,
                                     wasAromatic[bond->getIdx()]);
      res += "\n";
    }
    res += "M  V30 END BOND\n";
  }

  if (nSGroups > 0) {
    res += "M  V30 BEGIN SGROUP\n";
    unsigned int idx = 0;
    for (const auto &sgroup : sgroups) {
      res += GetV3000MolFileSGroupLines(++idx, sgroup);
    }
    res += "M  V30 END SGROUP\n";
  }

  if (tmol.hasProp(common_properties::molFileLinkNodes)) {
    auto pval = tmol.getProp<std::string>(common_properties::molFileLinkNodes);

    std::vector<std::string> linknodes;
    boost::split(linknodes, pval, boost::is_any_of("|"));
    for (const auto &linknode : linknodes) {
      res += "M  V30 LINKNODE " + linknode + "\n";
    }
  }
  appendEnhancedStereoGroups(res, tmol, wedgeBonds);

  res += "M  V30 END CTAB\n";
  return res;
}
}  // namespace FileParserUtils
enum class MolFileFormat { V2000, V3000, unspecified };

//------------------------------------------------
//
//  gets a mol block as a string
//
//------------------------------------------------
std::string outputMolToMolBlock(const RWMol &tmol, int confId,
                                MolFileFormat whichFormat,
                                unsigned int precision,
                                const boost::dynamic_bitset<> &aromaticBonds) {
  std::string res;
  unsigned int nAtoms, nBonds, nLists, chiralFlag, nsText, nRxnComponents;
  unsigned int nReactants, nProducts, nIntermediates;
  nAtoms = tmol.getNumAtoms();
  nBonds = tmol.getNumBonds();
  nLists = 0;

  const auto &sgroups = getSubstanceGroups(tmol);
  unsigned int nSGroups = sgroups.size();

  if (whichFormat == MolFileFormat::V2000 &&
      (nAtoms > 999 || nBonds > 999 || nSGroups > 999)) {
    throw ValueErrorException(
        "V2000 format does not support more than 999 atoms, bonds or SGroups.");
  }

  chiralFlag = 0;
  nsText = 0;
  nRxnComponents = 0;
  nReactants = 0;
  nProducts = 0;
  nIntermediates = 0;

  tmol.getPropIfPresent(common_properties::_MolFileChiralFlag, chiralFlag);

  const Conformer *conf;
  if (confId < 0 && tmol.getNumConformers() == 0) {
    conf = nullptr;
  } else {
    conf = &(tmol.getConformer(confId));
  }

  std::string text;
  if (tmol.getPropIfPresent(common_properties::_Name, text)) {
    res += text;
  }
  res += "\n";

  // info
  if (tmol.getPropIfPresent(common_properties::MolFileInfo, text)) {
    res += text;
  } else {
    std::stringstream ss;
    ss << "  " << std::setw(8) << "RDKit";
    ss << std::setw(10) << "";
    if (conf) {
      if (conf->is3D()) {
        ss << "3D";
      } else {
        ss << common_properties::TWOD;
      }
    }
    res += ss.str();
  }
  res += "\n";
  // comments
  if (tmol.getPropIfPresent(common_properties::MolFileComments, text)) {
    res += text;
  }
  res += "\n";

  bool hasDative = false;
  for (const auto bond : tmol.bonds()) {
    if (bond->getBondType() == Bond::DATIVE) {
      hasDative = true;
      break;
    }
  }

  bool isV3000 = false;
  if (whichFormat == MolFileFormat::V3000) {
    isV3000 = true;
  } else if (whichFormat == MolFileFormat::unspecified &&
             (hasDative || nAtoms > 999 || nBonds > 999 || nSGroups > 999 ||
              !tmol.getStereoGroups().empty())) {
    isV3000 = true;
  }

  // the counts line:
  std::stringstream ss;
  if (isV3000) {
    // All counts in the V3000 info line should be 0
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << "999 V3000\n";
  } else {
    ss << std::setw(3) << nAtoms;
    ss << std::setw(3) << nBonds;
    ss << std::setw(3) << nLists;
    ss << std::setw(3) << nSGroups;
    ss << std::setw(3) << chiralFlag;
    ss << std::setw(3) << nsText;
    ss << std::setw(3) << nRxnComponents;
    ss << std::setw(3) << nReactants;
    ss << std::setw(3) << nProducts;
    ss << std::setw(3) << nIntermediates;
    ss << "999 V2000\n";
  }
  res += ss.str();

  boost::dynamic_bitset<> queryListAtoms(tmol.getNumAtoms());
  if (!isV3000) {
    // V2000 output.
    for (ROMol::ConstAtomIterator atomIt = tmol.beginAtoms();
         atomIt != tmol.endAtoms(); ++atomIt) {
      res += GetMolFileAtomLine(*atomIt, conf, queryListAtoms);
      res += "\n";
    }

    auto wedgeBonds = Chirality::pickBondsToWedge(tmol, nullptr, conf);

    for (const auto bond : tmol.bonds()) {
      res += GetMolFileBondLine(bond, wedgeBonds, conf,
                                aromaticBonds[bond->getIdx()]);
      res += "\n";
    }

    res += GetMolFileChargeInfo(tmol);
    res += GetMolFileRGroupInfo(tmol);
    res += GetMolFileQueryInfo(tmol, queryListAtoms);
    res += GetMolFileAliasInfo(tmol);
    res += GetMolFileZBOInfo(tmol);

    res += GetMolFilePXAInfo(tmol);
    res += GetMolFileSGroupInfo(tmol);

    // FIX: R-group logic, SGroups and 3D features etc.
  } else {
    // V3000 output.
    res +=
        FileParserUtils::getV3000CTAB(tmol, aromaticBonds, confId, precision);
  }
  res += "M  END\n";
  return res;
}

void prepareMol(RWMol &trwmol, const MolWriterParams &params,
                boost::dynamic_bitset<> &aromaticBonds) {
  // NOTE: kekulize the molecule before writing it out
  // because of the way mol files handle aromaticity
  if (trwmol.needsUpdatePropertyCache()) {
    trwmol.updatePropertyCache(false);
  }
  if (params.kekulize && trwmol.getNumBonds()) {
    for (const auto bond : trwmol.bonds()) {
      if (bond->getIsAromatic()) {
        aromaticBonds.set(bond->getIdx());
      }
    }
    MolOps::Kekulize(trwmol);
  }

  if (params.includeStereo && !trwmol.getNumConformers()) {
    // generate coordinates so that the stereo we generate makes sense
    RDDepict::compute2DCoords(trwmol);
  }
  FileParserUtils::moveAdditionalPropertiesToSGroups(trwmol);
}

std::string MolToMolBlock(const ROMol &mol, const MolWriterParams &params,
                          int confId) {
  RDKit::Utils::LocaleSwitcher switcher;
  RWMol trwmol(mol);
  boost::dynamic_bitset<> aromaticBonds(trwmol.getNumBonds());
  prepareMol(trwmol, params, aromaticBonds);
  MolFileFormat whichFormat =
      params.forceV3000 ? MolFileFormat::V3000 : MolFileFormat::unspecified;
  return outputMolToMolBlock(trwmol, confId, whichFormat, params.precision,
                             aromaticBonds);
}

std::string MolToV2KMolBlock(const ROMol &mol, const MolWriterParams &params,
                             int confId) {
  RDKit::Utils::LocaleSwitcher switcher;
  RWMol trwmol(mol);
  boost::dynamic_bitset<> aromaticBonds(trwmol.getNumBonds());
  prepareMol(trwmol, params, aromaticBonds);
  return outputMolToMolBlock(trwmol, confId, MolFileFormat::V2000,
                             params.precision, aromaticBonds);
}

//------------------------------------------------
//
//  Dump a molecule to a file
//
//------------------------------------------------
void MolToMolFile(const ROMol &mol, const std::string &fName,
                  const MolWriterParams &params, int confId) {
  auto *outStream = new std::ofstream(fName.c_str());
  if (!(*outStream) || outStream->bad()) {
    delete outStream;
    std::ostringstream errout;
    errout << "Bad output file " << fName;
    throw BadFileException(errout.str());
  }
  std::string outString = MolToMolBlock(mol, params, confId);
  *outStream << outString;
  delete outStream;
}
}  // namespace RDKit
