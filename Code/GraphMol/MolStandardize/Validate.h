//
//  Copyright (C) 2018-2021 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Validate.h

        \brief Defines the ValidationErrorInfo class and four different
   validation methods: RDKitValidation, MolVSValidation, AllowedAtomsValidation,
   DisallowedAtomsValidation.

*/
#include <RDGeneral/export.h>
#ifndef RD_VALIDATE_H
#define RD_VALIDATE_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <iostream>
#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

//! The ValidationErrorInfo class is used to store the information returned by a
/// ValidationMethod validate.
using ValidationErrorInfo = std::string;

//! The ValidationMethod class is the abstract base class upon which all the
/// four different ValidationMethods inherit from.
class RDKIT_MOLSTANDARDIZE_EXPORT ValidationMethod {
 public:
  ValidationMethod() = default;
  virtual ~ValidationMethod() = default;

  virtual std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const = 0;
  virtual std::shared_ptr<ValidationMethod> copy() const = 0;
};

//! The CompositeValidation class provides a simple way to apply a collection of
// ValidationMethod instances in sequence
class RDKIT_MOLSTANDARDIZE_EXPORT CompositeValidation : public ValidationMethod {
 public:
  CompositeValidation(
    const std::vector<std::shared_ptr<ValidationMethod>> & validations)
    : validations(validations) {};

  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

  std::shared_ptr<ValidationMethod> copy() const override {
    return std::make_shared<CompositeValidation>(*this);
  }

 private:
   std::vector<std::shared_ptr<ValidationMethod>> validations;
};

//! The RDKitValidation class throws an error when there are no atoms in the
/// molecule or when there is incorrect atom valency.
/*!

  <b>Notes:</b>
    - RDKit automatically throws up atom valency issues but this class was made
  for completeness of the project.
*/
class RDKIT_MOLSTANDARDIZE_EXPORT RDKitValidation : public ValidationMethod {
 public:
  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

  std::shared_ptr<ValidationMethod> copy() const override {
    return std::make_shared<RDKitValidation>(*this);
  }
};

//////////////////////////////
/// MolVS Validations
//

//! The NoAtomValidation class throws an error if no atoms are present in the
/// molecule.
class RDKIT_MOLSTANDARDIZE_EXPORT NoAtomValidation : public ValidationMethod {
 public:
  std::vector<ValidationErrorInfo> validate(
    const ROMol &mol, bool reportAllFailures) const override;

  std::shared_ptr<ValidationMethod> copy() const override {
    return std::make_shared<NoAtomValidation>(*this);
  }
};

//! The FragmentValidation class logs if certain fragments are present.
class RDKIT_MOLSTANDARDIZE_EXPORT FragmentValidation : public ValidationMethod {
 public:
  std::vector<ValidationErrorInfo> validate(
    const ROMol &mol, bool reportAllFailures) const override;

  std::shared_ptr<ValidationMethod> copy() const override {
    return std::make_shared<FragmentValidation>(*this);
  }
};

//! The NeutralValidation class logs if not an overall neutral system.
class RDKIT_MOLSTANDARDIZE_EXPORT NeutralValidation : public ValidationMethod {
 public:
  std::vector<ValidationErrorInfo> validate(
    const ROMol &mol, bool reportAllFailures) const override;

  std::shared_ptr<ValidationMethod> copy() const override {
    return std::make_shared<NeutralValidation>(*this);
  }
};

//! The IsotopeValidation class logs if molecule contains isotopes.
class RDKIT_MOLSTANDARDIZE_EXPORT IsotopeValidation : public ValidationMethod {
 public:
  std::vector<ValidationErrorInfo> validate(
    const ROMol &mol, bool reportAllFailures) const override;

  std::shared_ptr<ValidationMethod> copy() const override {
    return std::make_shared<IsotopeValidation>(*this);
  }
};

////////////////////////////////
//! The MolVSValidation class includes most of the same validations as
/// molvs.validations, namely NoAtomValidation, FragmentValidation,
/// NeutralValidation, IsotopeValidation. MolVS also has IsNoneValidation and
/// DichloroethaneValidation but these were not included here (yet).
class RDKIT_MOLSTANDARDIZE_EXPORT MolVSValidation : public CompositeValidation {
 public:
  // constructor
  MolVSValidation();
  //! overloaded constructor to take in a user-defined list of ValidationMethod
  MolVSValidation(
      const std::vector<std::shared_ptr<ValidationMethod>> & validations);

  std::shared_ptr<ValidationMethod> copy() const override {
    return std::make_shared<MolVSValidation>(*this);
  }
};

//! The AllowedAtomsValidation class lets the user input a list of atoms,
//! anything not on
/// the list throws an error.
class RDKIT_MOLSTANDARDIZE_EXPORT AllowedAtomsValidation
    : public ValidationMethod {
 public:
  AllowedAtomsValidation(std::vector<std::shared_ptr<Atom>> atoms)
      : d_allowedList(std::move(atoms)) {}
  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

  std::shared_ptr<ValidationMethod> copy() const override {
    return std::make_shared<AllowedAtomsValidation>(*this);
  }

 private:
  std::vector<std::shared_ptr<Atom>> d_allowedList;
};

//! The DisallowedAtomsValidation class lets the user input a list of atoms and
//! as long
/// as there are no atoms from the list it is deemed acceptable.
class RDKIT_MOLSTANDARDIZE_EXPORT DisallowedAtomsValidation
    : public ValidationMethod {
 public:
  DisallowedAtomsValidation(std::vector<std::shared_ptr<Atom>> atoms)
      : d_disallowedList(std::move(atoms)) {}
  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

  std::shared_ptr<ValidationMethod> copy() const override {
    return std::make_shared<DisallowedAtomsValidation>(*this);
  }

 private:
  std::vector<std::shared_ptr<Atom>> d_disallowedList;
};

//! A convenience function for quickly validating a single SMILES string.
RDKIT_MOLSTANDARDIZE_EXPORT std::vector<ValidationErrorInfo> validateSmiles(
    const std::string &smiles);

}  // namespace MolStandardize
}  // namespace RDKit

#endif
