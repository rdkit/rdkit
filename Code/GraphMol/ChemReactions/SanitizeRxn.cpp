//
//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "SanitizeRxn.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>

namespace RDKit {
namespace RxnOps {

// molAtomMapNumber ==> int
// molFileRLabel ==> unsigned int
namespace {
template <class T>
T getMaxProp(ChemicalReaction &rxn, const std::string &prop) {
  T max_atom = (T)0;
  for (auto it = rxn.beginReactantTemplates(); it != rxn.endReactantTemplates();
       ++it) {
    for (auto atom : (*it)->atoms()) {
      T map;
      if (atom->getPropIfPresent<T>(prop, map)) {
        if (map > max_atom) {
          max_atom = map;
        }
      }
    }
  }

  for (auto it = rxn.beginAgentTemplates(); it != rxn.endAgentTemplates();
       ++it) {
    for (auto atom : (*it)->atoms()) {
      T map;
      if (atom->getPropIfPresent<T>(prop, map)) {
        if (map > max_atom) {
          max_atom = map;
        }
      }
    }
  }

  for (auto it = rxn.beginProductTemplates(); it != rxn.endProductTemplates();
       ++it) {
    for (auto atom : (*it)->atoms()) {
      T map;
      if (atom->getPropIfPresent<T>(prop, map)) {
        if (map > max_atom) {
          max_atom = map;
        }
      }
    }
  }

  return max_atom;
}

struct AtomInfo {
  Atom *atom;
  unsigned int templateIdx;
  unsigned int rlabel;
  int atomMap;
  int isotope;
  std::string dummyLabel;
  AtomInfo(Atom *at, unsigned int templateIdx)
      : atom(at),
        templateIdx(templateIdx),
        rlabel(0),
        atomMap(0),
        isotope(at->getIsotope()),
        dummyLabel() {
    atom->getPropIfPresent(common_properties::_MolFileRLabel, rlabel);
    atom->getPropIfPresent(common_properties::molAtomMapNumber, atomMap);
    atom->getPropIfPresent(common_properties::dummyLabel, dummyLabel);
    // std::cerr << atom->getIdx() << " : " << atom->getAtomicNum() << " " <<
    //    " rgroup: " << rlabel << " atomMap " << atomMap << " isotope " <<
    //    isotope <<
    //    " label " << dummyLabel <<
    //    std::endl;
  }

  bool NeedsRLabel() { return atom->getAtomicNum() == 0 && rlabel == 0; }

  unsigned int bestGuessRLabel() {
    if (rlabel) {
      return rlabel;
    }
    if (isotope) {
      return isotope;
    }
    if (atomMap) {
      return atomMap;
    }
    if (dummyLabel.size()) {
      try {
        return boost::lexical_cast<unsigned int>(
            dummyLabel.substr(1, dummyLabel.size() - 1));
      } catch (...) {
        return 0;
      }
    }
    return 0;
  }

  void setRLabel(unsigned int rlabel) {
    PRECONDITION(atom, "Internal error in SanitizeRxn - null atom");
    auto &mol = dynamic_cast<RWMol &>(atom->getOwningMol());

    QueryAtom qatom(*atom);
    qatom.setProp(common_properties::_MolFileRLabel, rlabel);
    std::string dLabel = "R" + std::to_string(rlabel);
    qatom.setProp(common_properties::dummyLabel, dLabel);
    if (rlabel > 0 && rlabel < 999) {
      qatom.setIsotope(rlabel);
    }
    qatom.setQuery(makeAtomNullQuery());
    unsigned int idx = atom->getIdx();
    mol.replaceAtom(idx, &qatom);
    atom = mol.getAtomWithIdx(idx);
  }

  void setAtomMap(int map) {
    atom->setProp(common_properties::molAtomMapNumber, map);
  }
};

std::string makeReactantErrorMessage(const std::string &error,
                                     const AtomInfo &at) {
  std::ostringstream str;
  str << error << " for reactant idx: " << at.templateIdx
      << " atom: " << at.atom->getIdx();
  return str.str();
}

std::string makeProductErrorMessage(const std::string &error,
                                    const AtomInfo &at) {
  std::ostringstream str;
  str << error << " for product idx: " << at.templateIdx
      << " atom: " << at.atom->getIdx();
  return str.str();
}
}  // namespace

// if we have query atoms without rlabels, make proper rlabels if possible
//  ensure that every rlabel in the reactant has one in the product
void fixRGroups(ChemicalReaction &rxn) {
  std::map<unsigned int, unsigned int> remapped;
  std::vector<AtomInfo> reactantAtomsToFix;
  std::vector<AtomInfo> productAtomsToFix;

  unsigned int templateIdx = 0;
  for (auto it = rxn.beginReactantTemplates(); it != rxn.endReactantTemplates();
       ++it, ++templateIdx) {
    for (auto atom : (*it)->atoms()) {
      AtomInfo at(atom, templateIdx);
      if (at.NeedsRLabel()) {
        reactantAtomsToFix.push_back(at);
      }
    }
  }

  templateIdx = 0;
  for (auto it = rxn.beginProductTemplates(); it != rxn.endProductTemplates();
       ++it, ++templateIdx) {
    for (auto atom : (*it)->atoms()) {
      AtomInfo at(atom, templateIdx);
      if (at.NeedsRLabel()) {
        productAtomsToFix.push_back(at);
      }
    }
  }

  if (!reactantAtomsToFix.size() && !productAtomsToFix.size()) {
    return;
  }

  if (reactantAtomsToFix.size() > productAtomsToFix.size()) {
    std::ostringstream str;
    str << "Mismatched potential rlabels: " << reactantAtomsToFix.size()
        << " unmapped reactant dummy atom rlabels," << productAtomsToFix.size()
        << " unmappped product dummy atom rlabels";
    BOOST_LOG(rdWarningLog) << str.str() << std::endl;
  }

  auto max_rlabel =
      getMaxProp<unsigned int>(rxn, common_properties::_MolFileRLabel);
  int max_atom_map = getMaxProp<int>(rxn, common_properties::molAtomMapNumber);

  for (auto &rat : reactantAtomsToFix) {
    bool found = false;
    unsigned int bestGuess = rat.bestGuessRLabel();
    if (!bestGuess) {
      continue;
    }

    for (auto &pat : productAtomsToFix) {
      if (!pat.atom) {
        continue;
      }

      if (rat.bestGuessRLabel() == pat.bestGuessRLabel()) {
        // if the atomMaps don't match, this is bad, no atomMap is ok(==0)
        if (rat.atomMap == pat.atomMap) {
          found = true;
          rat.setRLabel(max_rlabel + rat.bestGuessRLabel());
          pat.setRLabel(max_rlabel + pat.bestGuessRLabel());
          if (!rat.atomMap) {  // set atom mapping as well
            rat.setAtomMap(max_atom_map + rat.bestGuessRLabel());
            pat.setAtomMap(max_atom_map + rat.bestGuessRLabel());
          }
          pat.atom = nullptr;  // don't match again
          break;
        }
      }
    }

    if (!found) {
      BOOST_LOG(rdWarningLog)
          << "Could not find RLabel mapping for atom: " << rat.atom->getIdx()
          << " in template: " << rat.templateIdx << std::endl;
    }
  }
  return;
}

// if we have query atoms without rlabels, make proper rlabels if possible
//  ensure that every rlabel in the reactant has one in the product
void fixAtomMaps(ChemicalReaction &rxn) {
  int max_atom_map = getMaxProp<int>(rxn, common_properties::molAtomMapNumber);
  std::map<unsigned int, int> potential_mappings;

  unsigned int templateIdx = 0;

  for (auto it = rxn.beginReactantTemplates(); it != rxn.endReactantTemplates();
       ++it, ++templateIdx) {
    for (auto atom : (*it)->atoms()) {
      AtomInfo at(atom, templateIdx);
      if (at.rlabel && !at.atomMap) {
        if (potential_mappings.find(at.rlabel) != potential_mappings.end()) {
          throw RxnSanitizeException(std::string("Duplicated RLabels"));
        }
        int map = potential_mappings[at.rlabel] =
            rdcast<int>(at.rlabel) + max_atom_map;
        at.setAtomMap(map);
      }
    }
  }

  if (!potential_mappings.size()) {
    return;  // everything is ok!
  }

  templateIdx = 0;
  for (auto it = rxn.beginProductTemplates(); it != rxn.endProductTemplates();
       ++it, ++templateIdx) {
    for (auto atom : (*it)->atoms()) {
      AtomInfo at(atom, templateIdx);
      if (at.rlabel) {
        if (!at.atomMap) {
          at.setAtomMap(potential_mappings[at.rlabel]);
        } else {
          if (at.atomMap != potential_mappings[at.rlabel]) {
            throw RxnSanitizeException(makeProductErrorMessage(
                "RLabel is mapped in product but not in reactant", at));
          }
        }
      }
    }
  }
}

// might throw mol sanitization exception??? wrap in RxnSanitize?
void fixReactantTemplateAromaticity(ChemicalReaction &rxn) {
  unsigned int ops;
  for (auto it = rxn.beginReactantTemplates(); it != rxn.endReactantTemplates();
       ++it) {
    auto *rw = dynamic_cast<RWMol *>(it->get());
    if (rw) {
      sanitizeMol(*rw, ops, MolOps::SANITIZE_SETAROMATICITY);
    } else
      PRECONDITION(rw, "Oops, not really a RWMol?");
  }
}

void fixHs(ChemicalReaction &rxn) {
  {
    //  if mapped Hydrogens in reactants area mapped to heavy atoms
    //   keep mappings, in all other cases remove them.
    //   this allows us to merge query hydrogens atoms

    std::map<int, bool> mappedToNonHeavyProductAtom;

    for (auto it = rxn.beginProductTemplates(); it != rxn.endProductTemplates();
         ++it) {
      int atomMap = 0;
      for (auto atom : (*it)->atoms()) {
        if (atom->getAtomicNum() != 1) {  // hydrogen
          if (atom->getPropIfPresent(common_properties::molAtomMapNumber,
                                     atomMap)) {
            if (atomMap) {
              mappedToNonHeavyProductAtom[atomMap] = true;
            }
          }
        }
      }
    }

    for (auto it = rxn.beginReactantTemplates();
         it != rxn.endReactantTemplates(); ++it) {
      int atomMap = 0;
      for (auto atom : (*it)->atoms()) {
        if (atom->getAtomicNum() == 1) {  // hydrogen
          if (atom->getPropIfPresent(common_properties::molAtomMapNumber,
                                     atomMap)) {
            if (atomMap) {
              if (mappedToNonHeavyProductAtom.find(atomMap) ==
                  mappedToNonHeavyProductAtom.end()) {
                atom->clearProp(common_properties::molAtomMapNumber);
              } else {
                BOOST_LOG(rdWarningLog)
                    << "Reaction has explicit hydrogens, reactants will need "
                       "explicit hydrogens (addHs)"
                    << std::endl;
              }
            }
          }
        }
      }
    }
  }

  const bool mergeUnmappedOnly = true;
  for (auto it = rxn.beginReactantTemplates(); it != rxn.endReactantTemplates();
       ++it) {
    auto *rw = dynamic_cast<RWMol *>(it->get());
    if (rw) {
      MolOps::mergeQueryHs(*rw, mergeUnmappedOnly);
    } else
      PRECONDITION(rw, "Oops, not really an RWMol?");
  }
}

void adjustTemplates(const MOL_SPTR_VECT &templates,
                     const MolOps::AdjustQueryParameters &params) {
  for (auto &templ : templates) {
    auto *rw = dynamic_cast<RWMol *>(templ.get());
    if (rw) {
      adjustQueryProperties(*rw, &params);
    } else
      PRECONDITION(rw, "Oops, not really a RWMol?");
  }
}
void sanitizeRxn(ChemicalReaction &rxn, unsigned int &operationsThatFailed,
                 unsigned int ops,
                 const MolOps::AdjustQueryParameters &params) {
  operationsThatFailed = SANITIZE_NONE;

  if (ops & SANITIZE_RGROUP_NAMES) {
    operationsThatFailed = SANITIZE_RGROUP_NAMES;
    fixRGroups(rxn);
  }

  if (ops & SANITIZE_ATOM_MAPS) {
    operationsThatFailed = SANITIZE_ATOM_MAPS;
    fixAtomMaps(rxn);
  }
  if (ops & SANITIZE_ADJUST_REACTANTS) {
    operationsThatFailed = SANITIZE_ADJUST_REACTANTS;
    adjustTemplates(rxn.getReactants(), params);
  }
  if (ops & SANITIZE_ADJUST_PRODUCTS) {
    operationsThatFailed = SANITIZE_ADJUST_PRODUCTS;
    adjustTemplates(rxn.getProducts(), params);
  }
  if (ops & SANITIZE_MERGEHS) {
    operationsThatFailed = SANITIZE_MERGEHS;
    fixHs(rxn);
  }

  operationsThatFailed = SANITIZE_NONE;
}

void sanitizeRxn(ChemicalReaction &rxn,
                 const MolOps::AdjustQueryParameters &params) {
  unsigned int ops = 0;
  return sanitizeRxn(rxn, ops, SANITIZE_ALL, params);
}

void sanitizeRxnAsMols(ChemicalReaction &rxn, unsigned int sanitizeOps) {
  for (auto &mol : rxn.getReactants()) {
    unsigned int operationThatFailed;
    MolOps::sanitizeMol(*dynamic_cast<RWMol *>(mol.get()), operationThatFailed,
                        sanitizeOps);
  }
  for (auto &mol : rxn.getAgents()) {
    unsigned int operationThatFailed;
    MolOps::sanitizeMol(*dynamic_cast<RWMol *>(mol.get()), operationThatFailed,
                        sanitizeOps);
  }
  for (auto &mol : rxn.getProducts()) {
    unsigned int operationThatFailed;
    MolOps::sanitizeMol(*dynamic_cast<RWMol *>(mol.get()), operationThatFailed,
                        sanitizeOps);
  }
}

}  // namespace RxnOps
}  // namespace RDKit
