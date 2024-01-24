//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Validate.h"
#include "Fragment.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentCatalogParams.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <iostream>
#include <vector>
#include <string>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace std;
using namespace RDKit;

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

std::vector<ValidationErrorInfo> CompositeValidation::validate(
  const ROMol &mol, bool reportAllFailures) const
{
  std::vector<ValidationErrorInfo> errors;
  for (const auto & method : validations) {
    auto partial = method->validate(mol, reportAllFailures);
    if (!partial.empty()) {
      std::copy(partial.begin(), partial.end(), std::back_inserter(errors));
      if (!reportAllFailures) {
        break;
      }
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> RDKitValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  ROMol molCopy = mol;
  std::vector<ValidationErrorInfo> errors;

  unsigned int na = mol.getNumAtoms();

  if (!na) {
    errors.emplace_back("ERROR: [NoAtomValidation] Molecule has no atoms");
  }

  // loop over atoms
  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    Atom *atom = molCopy.getAtomWithIdx(i);
    try {
      atom->calcExplicitValence();
    } catch (const MolSanitizeException &e) {
      errors.push_back("INFO: [ValenceValidation] " + std::string(e.what()));
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo>
NoAtomValidation::validate(const ROMol &mol, bool /*reportAllFailures*/) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();
  if (!na) {
    errors.emplace_back("ERROR: [NoAtomValidation] Molecule has no atoms");
  }
  return errors;
}

std::vector<ValidationErrorInfo>
FragmentValidation::validate(const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  // REVIEW: reportAllFailures is not being used here. is that correct?
  RDUNUSED_PARAM(reportAllFailures);
  std::shared_ptr<FragmentCatalogParams> fparams(new FragmentCatalogParams(""));
  FragmentCatalog fcat(fparams.get());

  const std::vector<std::shared_ptr<ROMol>> &fgrps = fparams->getFuncGroups();
  INT_VECT mapping;
  VECT_INT_VECT atom_mapping;
  std::vector<ROMOL_SPTR> frags =
      MolOps::getMolFrags(mol, true, &mapping, &atom_mapping);

  for (auto &fgrp : fgrps) {
    std::string fname;
    fgrp->getProp(common_properties::_Name, fname);
    std::vector<RDKit::MatchVectType> res;
    unsigned int matches = SubstructMatch(mol, *fgrp, res);
    //		std::cout << fname << " matches " << matches << std::endl;
    if (matches != 0 && frags.size() != 0) {
      VECT_INT_VECT substructmap;  // store idxs of frag from substructmatch
      for (const auto &match : res) {
        std::vector<int> vec;
        for (const auto &pair : match) {
          vec.push_back(pair.second);
        }
        substructmap.push_back(vec);
      }

      // to stop the same fragment being reported many times if present
      // multiple times in molecule
      bool fpresent = false;

      for (auto &molfragidx : atom_mapping) {
        std::sort(molfragidx.begin(), molfragidx.end());
        for (auto &substructidx : substructmap) {
          std::sort(substructidx.begin(), substructidx.end());
          //					// help to debug...
          //					std::cout << "molfragidx: "  <<
          // std::endl; 					for (const auto
          // &i : molfragidx)
          // {
          // std::cout << i; }; std::cout
          // << std::endl; std::cout <<
          //"substructidx: "  << std::endl;
          // for (const auto &i : substructidx) { std::cout << i; };
          // std::cout <<
          // std::endl;
          //					//
          if ((molfragidx == substructidx) && !fpresent) {
            std::string msg = fname + " is present";
            errors.push_back("INFO: [FragmentValidation] " + msg);
            fpresent = true;
          }
        }
      }
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo>
NeutralValidation::validate(const ROMol &mol, bool /*reportAllFailures*/) const {
  std::vector<ValidationErrorInfo> errors;
  int charge = RDKit::MolOps::getFormalCharge(mol);
  if (charge != 0) {
    std::string charge_str;
    if (charge > 0) {
      charge_str = "+" + std::to_string(charge);
    } else {
      charge_str = std::to_string(charge);
    }
    std::string msg = "Not an overall neutral system (" + charge_str + ')';
    errors.push_back("INFO: [NeutralValidation] " + msg);
  }
  return errors;
}

std::vector<ValidationErrorInfo>
IsotopeValidation::validate(const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();
  std::set<string> isotopes;

  // loop over atoms
  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    const Atom *atom = mol.getAtomWithIdx(i);
    unsigned int isotope = atom->getIsotope();
    if (isotope != 0) {
      std::string symbol = atom->getSymbol();
      isotopes.insert(std::to_string(isotope) + symbol);
    }
  }

  for (auto &isotope : isotopes) {
    errors.push_back("INFO: [IsotopeValidation] Molecule contains isotope " +
                     isotope);
  }
  return errors;
}

// constructor
MolVSValidation::MolVSValidation()
  : CompositeValidation({
      std::make_shared<NoAtomValidation>(),
      std::make_shared<FragmentValidation>(),
      std::make_shared<NeutralValidation>(),
      std::make_shared<IsotopeValidation>()
      })
{
}

// overloaded constructor
MolVSValidation::MolVSValidation(
    const std::vector<std::shared_ptr<ValidationMethod>> & validations)
  : CompositeValidation(validations)
{
}

std::vector<ValidationErrorInfo> AllowedAtomsValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();

  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    const Atom *qatom = mol.getAtomWithIdx(i);
    bool match = false;
    // checks to see qatom matches one of list of allowedAtoms
    for (const auto &allowedAtom : this->d_allowedList) {
      if (allowedAtom->Match(qatom)) {
        match = true;
        break;
      }
    }
    // if no match, append to list of errors.
    if (!match) {
      std::string symbol = qatom->getSymbol();
      errors.push_back("INFO: [AllowedAtomsValidation] Atom " + symbol +
                       " is not in allowedAtoms list");
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> DisallowedAtomsValidation::validate(
    const ROMol &mol, bool reportAllFailures) const {
  std::vector<ValidationErrorInfo> errors;
  unsigned int na = mol.getNumAtoms();

  for (size_t i = 0; i < na; ++i) {
    if (!reportAllFailures) {
      if (errors.size() >= 1) {
        break;
      }
    }
    const Atom *qatom = mol.getAtomWithIdx(i);
    bool match = false;
    // checks to see qatom matches one of list of allowedAtoms
    for (const auto &disallowedAtom : this->d_disallowedList) {
      if (disallowedAtom->Match(qatom)) {
        match = true;
      }
    }
    // if no match, append to list of errors.
    if (match) {
      std::string symbol = qatom->getSymbol();
      errors.push_back("INFO: [DisallowedAtomsValidation] Atom " + symbol +
                       " is in disallowedAtoms list");
    }
  }
  return errors;
}

std::vector<ValidationErrorInfo> validateSmiles(const std::string &smiles) {
  RWMOL_SPTR mol(SmilesToMol(smiles));
  if (!mol) {
    std::string message =
        "SMILES Parse Error: syntax error for input: " + smiles;
    throw ValueErrorException(message);
  }

  MolVSValidation vm;
  std::vector<ValidationErrorInfo> errors = vm.validate(*mol, true);

  return errors;
}

}  // namespace MolStandardize
}  // namespace RDKit
