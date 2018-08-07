//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
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

class ValidationErrorInfo : public std::exception {
 public:
  ValidationErrorInfo(const std::string &msg) : d_msg(msg) {
    BOOST_LOG(rdInfoLog) << d_msg << std::endl;
  };
  const char *message() const { return d_msg.c_str(); };
  ~ValidationErrorInfo() throw(){};

 private:
  std::string d_msg;
};  // class ValidationErrorInfo

class ValidationMethod {
 public:
  virtual std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const = 0;
};

class RDKitValidation : public ValidationMethod {
 public:
  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;
};

//////////////////////////////
// MolVS Validations
//
class MolVSValidations {
 public:
  virtual void run(const ROMol &mol, bool reportAllFailures,
                   std::vector<ValidationErrorInfo> &errors) const = 0;
  virtual MolVSValidations *copy() const = 0;
};

class NoAtomValidation : public MolVSValidations {
 public:
  void run(const ROMol &mol, bool reportAllFailures,
           std::vector<ValidationErrorInfo> &errors) const override;
  virtual MolVSValidations *copy() const override {
    return new NoAtomValidation(*this);
  };
};

class FragmentValidation : public MolVSValidations {
 public:
  void run(const ROMol &mol, bool reportAllFailures,
           std::vector<ValidationErrorInfo> &errors) const override;
  virtual MolVSValidations *copy() const override {
    return new FragmentValidation(*this);
  };
};

class NeutralValidation : public MolVSValidations {
 public:
  void run(const ROMol &mol, bool reportAllFailures,
           std::vector<ValidationErrorInfo> &errors) const override;
  virtual MolVSValidations *copy() const override {
    return new NeutralValidation(*this);
  };
};

class IsotopeValidation : public MolVSValidations {
 public:
  void run(const ROMol &mol, bool reportAllFailures,
           std::vector<ValidationErrorInfo> &errors) const override;
  virtual MolVSValidations *copy() const override {
    return new IsotopeValidation(*this);
  };
};

////////////////////////////////

class MolVSValidation : public ValidationMethod {
 public:
  // constructor
  MolVSValidation();
  // overloaded constructor
  MolVSValidation(const std::vector<MolVSValidations *> validations);
  MolVSValidation(const MolVSValidation &other);
  ~MolVSValidation();

  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

 private:
  std::vector<MolVSValidations *> d_validations;
};

class AllowedAtomsValidation : public ValidationMethod {
 public:
  AllowedAtomsValidation(const std::vector<std::shared_ptr<Atom>> &atoms)
      : d_allowedList(atoms){};
  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

 private:
  std::vector<std::shared_ptr<Atom>> d_allowedList;
  // void initializeDefaultAtoms; // TODO with filtersCatalog
  // stuff
};

class DisallowedAtomsValidation : public ValidationMethod {
 public:
  DisallowedAtomsValidation(const std::vector<std::shared_ptr<Atom>> &atoms)
      : d_disallowedList(atoms){};
  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

 private:
  std::vector<std::shared_ptr<Atom>> d_disallowedList;
  // void initializeDefaultAtoms; // TODO with filtersCatalog
  // stuff
};

}  // namespace MolStandardize
}  // namespace RDKit

#endif
