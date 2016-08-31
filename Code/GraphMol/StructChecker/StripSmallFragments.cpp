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
#include "../SmilesParse/SmilesWrite.h"

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

bool StripSmallFragments(RWMol &mol, bool verbose) {
  const bool sanitize=false;
  std::vector<boost::shared_ptr<ROMol> > frags = MolOps::getMolFrags(mol, sanitize);
  size_t maxFragSize = 0;
  size_t maxFragIdx = 0;
  if (frags.size() == 1)
    return false;
  
  for(size_t i=0; i<frags.size(); ++i) {
    const unsigned int fragSize = frags[i].get()->getNumAtoms();
    if(fragSize >= maxFragSize) {
      maxFragSize = fragSize;
      maxFragIdx = i;
    }
  }

  if(verbose) {
    std::string name = "<no name>";
    mol.getPropIfPresent(common_properties::_Name, name);
    for(size_t i=0; i<frags.size(); ++i) {
      if (i != maxFragIdx) {
        BOOST_LOG(rdWarningLog) << name << " removed " << frags[i].get()->getNumAtoms()
                                << " atoms" << std::endl;
      }
    }

  }
  mol = *frags[maxFragIdx].get();
  BOOST_LOG(rdInfoLog) << MolToSmiles(mol) << "\n";
    
  return true;
}

}  // namespace StructureCheck
}  // namespace RDKit
