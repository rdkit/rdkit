//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <map>

#include "../MolOps.h"
#include "../Descriptors/MolDescriptors.h"
#include "StripSmallFragments.h"

namespace RDKit {
namespace StructureCheck {

static inline std::string getMolecularFormula(const ROMol &mol) {
  return RDKit::Descriptors::calcMolFormula(mol);
}

void AddMWMF(RWMol &mol,
             bool pre) {  // set formula & mass properties "MW_PRE" "MW_POST"
  double mass = 0.0;
  mass = RDKit::Descriptors::calcExactMW(mol);
  /*
          for (unsigned i = 0; i < mol.getNumAtoms(); i++) {
               const Atom& atom = *mol.getAtomWithIdx(i);
               mass += atom.getMass();
               mass += atom.getNumImplicitHs() * 1.0080; // and add implicit
     Hydrogens mass
           }
  */
  std::string formula = getMolecularFormula(mol);
  if (!formula.empty()) mol.setProp((pre ? "MF_PRE" : "MF_POST"), formula);
  char propertyValue[64];
  snprintf(propertyValue, sizeof(propertyValue), "%g", mass);
  mol.setProp((pre ? "MW_PRE" : "MW_POST"), mass);
}

bool StripSmallFragments(RWMol &mol) {
  bool removed = false;
  // there may be an argument about how much sense this makes, but it's
  // consistent with the avalon toolkit behavior
  // if (mol.hasProp(RDKit::common_properties::_MolFileChiralFlag)) {
  //   mol.clearProp(RDKit::common_properties::_MolFileChiralFlag);
  // }
  std::vector<int> frags;
  std::map<unsigned, unsigned> frag_count;

  unsigned int nFrags = RDKit::MolOps::getMolFrags(mol, frags);
  if (nFrags > 1) {
    unsigned maxFragSize = 0;
    unsigned maxSizeFragIdx = 0;
    for (unsigned i = 0; i < frags.size(); i++) {
      if (frag_count.find(frags[i]) != frag_count.end()) {
        frag_count.at(frags[i]) = frag_count.at(frags[i]) + 1;
      } else {
        frag_count.insert(std::pair<unsigned, unsigned>(frags[i], 1));
      }

      if (frag_count.at(frags[i]) > maxFragSize) {
        maxFragSize = frag_count.at(frags[i]);
        maxSizeFragIdx = frags[i];
      }
    }
    for (int i = frags.size() - 1; i >= 0; i--)
      if (frags[i] != maxSizeFragIdx) {
        mol.removeAtom(i);
      }
    removed = true;
  }
  return removed;
}

}  // namespace StructureCheck
}  // namespace RDKit
