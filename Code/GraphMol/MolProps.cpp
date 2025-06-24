//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <map>
#include <string>
#include <strstream>

namespace RDKit {
namespace MolOps {
double getAvgMolWt(const ROMol &mol, bool onlyHeavy) {
  double res = 0.0;
  const PeriodicTable *table = PeriodicTable::getTable();
  for (const auto &atom : mol.atoms()) {
    if (!onlyHeavy || atom->getAtomicNum() != 1) {
      res += atom->getMass();
    }

    // add our implicit Hs if we need to:
    if (!onlyHeavy) {
      res += atom->getTotalNumHs() * table->getAtomicWeight(1);
    }
  }
  return res;
}

double getExactMolWt(const ROMol &mol, bool onlyHeavy) {
  double res = 0.0;
  int nHsToCount = 0;
  const PeriodicTable *table = PeriodicTable::getTable();
  for (auto atom : mol.atoms()) {
    auto atNum = atom->getAtomicNum();
    if (atNum != 1 || !onlyHeavy) {
      if (!atom->getIsotope()) {
        res += table->getMostCommonIsotopeMass(atNum);
      } else {
        res += atom->getMass();
      }
      res -= constants::electronMass * atom->getFormalCharge();
    }

    // add our implicit Hs if we need to:
    if (!onlyHeavy) {
      nHsToCount += atom->getTotalNumHs(false);
    }
  }
  if (!onlyHeavy) {
    res += nHsToCount * table->getMostCommonIsotopeMass(1);
  }
  return res;
}

namespace {
bool HillCompare(const std::pair<unsigned int, std::string> &v1,
                 const std::pair<unsigned int, std::string> &v2) {
  bool nCompare = (v1.first < v2.first);

  // special cases: Cs, Hs, Ds, and Ts go at the beginning
  if (v1.second == "C") {
    if (v2.second != "C") {
      return true;
    }
    return nCompare;
  } else if (v2.second == "C") {
    return false;
  }

  if (v1.second == "H") {
    if (v2.second != "H") {
      return true;
    }
    return nCompare;
  } else if (v2.second == "H") {
    return false;
  }

  if (v1.second == "D") {
    return true;
  } else if (v2.second == "D") {
    return false;
  }

  if (v1.second == "T") {
    return true;
  } else if (v2.second == "T") {
    return false;
  }

  // finally, just compare the symbols and the isotopes:
  if (v1 != v2) {
    return v1 < v2;
  } else {
    return nCompare;
  }
}
}  // namespace

std::string getMolFormula(const ROMol &mol, bool separateIsotopes,
                          bool abbreviateHIsotopes) {
  std::ostringstream res;
  std::map<std::pair<unsigned int, std::string>, unsigned int> counts;
  int charge = 0;
  unsigned int nHs = 0;
  const PeriodicTable *table = PeriodicTable::getTable();
  for (const auto atom : mol.atoms()) {
    int atNum = atom->getAtomicNum();
    std::pair<unsigned int, std::string> key =
        std::make_pair(0, table->getElementSymbol(atNum));
    if (separateIsotopes) {
      unsigned int isotope = atom->getIsotope();
      if (abbreviateHIsotopes && atNum == 1 && (isotope == 2 || isotope == 3)) {
        if (isotope == 2) {
          key.second = "D";
        } else {
          key.second = "T";
        }
      } else {
        key.first = isotope;
      }
    }
    counts[key] += 1;
    nHs += atom->getTotalNumHs();
    charge += atom->getFormalCharge();
  }

  if (nHs) {
    counts[{0, "H"}] += nHs;
  }
  std::list<std::pair<unsigned int, std::string>> ks;
  for (const auto &pr : counts) {
    ks.push_back(pr.first);
  }
  ks.sort(HillCompare);

  for (const auto &key : ks) {
    if (key.first > 0) {
      res << "[" << key.first << key.second << "]";
    } else {
      res << key.second;
    }
    if (counts[key] > 1) {
      res << counts[key];
    }
  }
  if (charge > 0) {
    res << "+";
    if (charge > 1) {
      res << charge;
    }
  } else if (charge < 0) {
    res << "-";
    if (charge < -1) {
      res << -1 * charge;
    }
  }
  return res.str();
}

}  // namespace MolOps
}  // namespace RDKit