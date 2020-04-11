//
//  Copyright (C) 2018 Susan H. Leung
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
#ifndef __RD_VALIDATE_H__
#define __RD_VALIDATE_H__

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <iostream>
#include <exception>
#include <string>
#include <vector>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

//! The ValidationErrorInfo class is used to store the information returned by a
// ValidationMethod validate.
class RDKIT_MOLSTANDARDIZE_EXPORT ValidationErrorInfo : public std::exception {
 public:
  ValidationErrorInfo(const std::string &msg) : d_msg(msg) {
    BOOST_LOG(rdInfoLog) << d_msg << std::endl;
  };
  const char *what() const noexcept override { return d_msg.c_str(); };
  ~ValidationErrorInfo() noexcept {};

 private:
  std::string d_msg;
};  // class ValidationErrorInfo

//! The ValidationMethod class is the abstract base class upon which all the
// four different ValidationMethods inherit from.
class RDKIT_MOLSTANDARDIZE_EXPORT ValidationMethod {
 public:
  ValidationMethod() = default;
  virtual ~ValidationMethod() = default;

  virtual std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const = 0;
};

//! The RDKitValidation class throws an error when there are no atoms in the
// molecule or when there is incorrect atom valency.
/*!

  <b>Notes:</b>
    - RDKit automatically throws up atom valency issues but this class was made
  for completeness of the project.
*/
class RDKIT_MOLSTANDARDIZE_EXPORT RDKitValidation : public ValidationMethod {
 public:
  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;
};

//////////////////////////////
// MolVS Validations
//
//! The MolVSValidations class includes most of the same validations as
// molvs.validations, namely NoAtomValidation, FragmentValidation,
// NeutralValidation, IsotopeValidation. MolVS also has IsNoneValidation and
// DichloroethaneValidation but these were not included here (yet).
class RDKIT_MOLSTANDARDIZE_EXPORT MolVSValidations {
 public:
  virtual void run(const ROMol &mol, bool reportAllFailures,
                   std::vector<ValidationErrorInfo> &errors) const = 0;
  virtual boost::shared_ptr<MolVSValidations> copy() const = 0;
};

//! The NoAtomValidation class throws an error if no atoms are present in the
// molecule.
class RDKIT_MOLSTANDARDIZE_EXPORT NoAtomValidation final
    : public MolVSValidations {
 public:
  void run(const ROMol &mol, bool reportAllFailures,
           std::vector<ValidationErrorInfo> &errors) const override;
  //! makes a copy of NoAtomValidation object and returns a MolVSValidations
  //! pointer to it
  virtual boost::shared_ptr<MolVSValidations> copy() const override {
    return boost::make_shared<NoAtomValidation>(*this);
  };
};

//! The FragmentValidation class logs if certain fragments are present.
class RDKIT_MOLSTANDARDIZE_EXPORT FragmentValidation final
    : public MolVSValidations {
 public:
  void run(const ROMol &mol, bool reportAllFailures,
           std::vector<ValidationErrorInfo> &errors) const override;
  //! makes a copy of FragmentValidation object and returns a MolVSValidations
  //! pointer to it
  virtual boost::shared_ptr<MolVSValidations> copy() const override {
    return boost::make_shared<FragmentValidation>(*this);
  };
};

//! The NeutralValidation class logs if not an overall neutral system.
class RDKIT_MOLSTANDARDIZE_EXPORT NeutralValidation final
    : public MolVSValidations {
 public:
  void run(const ROMol &mol, bool reportAllFailures,
           std::vector<ValidationErrorInfo> &errors) const override;
  //! makes a copy of NeutralValidation object and returns a MolVSValidations
  //! pointer to it
  virtual boost::shared_ptr<MolVSValidations> copy() const override {
    return boost::make_shared<NeutralValidation>(*this);
  };
};

//! The IsotopeValidation class logs if molecule contains isotopes.
class RDKIT_MOLSTANDARDIZE_EXPORT IsotopeValidation final
    : public MolVSValidations {
 public:
  void run(const ROMol &mol, bool reportAllFailures,
           std::vector<ValidationErrorInfo> &errors) const override;
  //! makes a copy of IsotopeValidation object and returns a MolVSValidations
  //! pointer to it
  virtual boost::shared_ptr<MolVSValidations> copy() const override {
    return boost::make_shared<IsotopeValidation>(*this);
  };
};

////////////////////////////////

//! The MolVSValidation class can be used to perform all MolVSValidions.
class RDKIT_MOLSTANDARDIZE_EXPORT MolVSValidation : public ValidationMethod {
 public:
  // constructor
  MolVSValidation();
  //! overloaded constructor to take in a user-defined list of MolVSValidations
  MolVSValidation(
      const std::vector<boost::shared_ptr<MolVSValidations>> validations);
  MolVSValidation(const MolVSValidation &other);
  ~MolVSValidation();

  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

 private:
  std::vector<boost::shared_ptr<MolVSValidations>> d_validations;
};

//! The AllowedAtomsValidation class lets the user input a list of atoms,
//! anything not on
// the list throws an error.
class RDKIT_MOLSTANDARDIZE_EXPORT AllowedAtomsValidation
    : public ValidationMethod {
 public:
  AllowedAtomsValidation(const std::vector<std::shared_ptr<Atom>> &atoms)
      : d_allowedList(atoms){};
  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

 private:
  std::vector<std::shared_ptr<Atom>> d_allowedList;
};

//! The DisallowedAtomsValidation class lets the user input a list of atoms and
//! as long
// as there are no atoms from the list it is deemed acceptable.
class RDKIT_MOLSTANDARDIZE_EXPORT DisallowedAtomsValidation
    : public ValidationMethod {
 public:
  DisallowedAtomsValidation(const std::vector<std::shared_ptr<Atom>> &atoms)
      : d_disallowedList(atoms){};
  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

 private:
  std::vector<std::shared_ptr<Atom>> d_disallowedList;
};

//! A convenience function for quickly validating a single SMILES string.
RDKIT_MOLSTANDARDIZE_EXPORT std::vector<ValidationErrorInfo> validateSmiles(
    const std::string &smiles);

}  // namespace MolStandardize
}  // namespace RDKit

#endif
