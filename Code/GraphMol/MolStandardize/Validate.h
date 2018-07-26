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

class MolVSValidation : public ValidationMethod {
 public:
  std::vector<ValidationErrorInfo> validate(
      const ROMol &mol, bool reportAllFailures) const override;

 private:
  void noAtomValidation(const ROMol &mol, bool reportAllFailures,
                        std::vector<ValidationErrorInfo> &errors) const;
  void fragmentValidation(const ROMol &mol, bool reportAllFailures,
                          std::vector<ValidationErrorInfo> &errors) const;
  void neutralValidation(const ROMol &mol, bool reportAllFailures,
                         std::vector<ValidationErrorInfo> &errors) const;
  void isotopeValidation(const ROMol &mol, bool reportAllFailures,
                         std::vector<ValidationErrorInfo> &errors) const;
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
