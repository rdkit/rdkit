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
#include "../FileParsers/MolFileStereochem.h"

// define snprintf for msvc
#if _MSC_VER
#if _MSC_VER < 1900
#define snprintf _snprintf
#endif
#endif

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
  const bool sanitize = false;
  std::vector<boost::shared_ptr<ROMol>> frags =
      MolOps::getMolFrags(mol, sanitize);
  if (frags.size() <= 1) return false;

  size_t maxFragSize = 0;
  size_t maxFragIdx = 0;

  for (size_t i = 0; i < frags.size(); ++i) {
    const unsigned int fragSize = frags[i].get()->getNumAtoms();
    if (fragSize >= maxFragSize) {
      maxFragSize = fragSize;
      maxFragIdx = i;
    }
  }

  if (verbose) {
    std::string name = "<no name>";
    mol.getPropIfPresent(common_properties::_Name, name);
    for (size_t i = 0; i < frags.size(); ++i) {
      if (i != maxFragIdx) {
        BOOST_LOG(rdWarningLog)
            << name << " removed fragment i=" << i << " with "
            << frags[i].get()->getNumAtoms() << " atoms" << std::endl;
      }
    }
  }

  // we need to save chirality for checking later
  bool checkChiral = false;
  if (mol.hasProp(RDKit::common_properties::_MolFileChiralFlag)) {
    unsigned int chiralflag =
        mol.getProp<unsigned int>(RDKit::common_properties::_MolFileChiralFlag);
    frags[maxFragIdx].get()->setProp<unsigned int>(
        RDKit::common_properties::_MolFileChiralFlag, chiralflag);
    checkChiral = chiralflag != 0;
  }

  mol = *frags[maxFragIdx].get();

  // We need to see if the mol file's chirality possibly came from this
  //  fragment.
  if (checkChiral) {
    bool ischiral = false;

    RWMol copy(mol);
    try {
      MolOps::sanitizeMol(copy);
      ClearSingleBondDirFlags(copy);
      MolOps::detectBondStereochemistry(copy);
      MolOps::assignStereochemistry(copy, true, true, true);
      for (ROMol::AtomIterator atIt = copy.beginAtoms();
           atIt != copy.endAtoms(); ++atIt) {
        if ((*atIt)->hasProp(common_properties::_ChiralityPossible)) {
          ischiral = true;
          checkChiral = false;
          break;
        }
      }
    } catch (...) {
    }

    // are chiral tags set
    if (checkChiral) {
      for (ROMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
           ++atIt) {
        if ((*atIt)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
            (*atIt)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
          ischiral = true;
          break;
        }
      }

      for (ROMol::BondIterator bondIt = mol.beginBonds();
           bondIt != mol.endBonds(); ++bondIt) {
        if ((*bondIt)->getBondDir() == Bond::BEGINDASH ||
            (*bondIt)->getBondDir() == Bond::BEGINWEDGE) {
          ischiral = true;
          break;
        }
      }
    }

    if (!ischiral) {
      mol.setProp<unsigned int>(RDKit::common_properties::_MolFileChiralFlag,
                                0);
    }
  }
  return true;
}

}  // namespace StructureCheck
}  // namespace RDKit
