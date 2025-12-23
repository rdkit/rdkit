//
//  Copyright (C) 2024 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_SCSRMOL_H
#define RD_SCSRMOL_H

#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>

namespace RDKit {

class RDKIT_GRAPHMOL_EXPORT SCSRMol {
 private:
  std::unique_ptr<RWMol> p_mol;
  std::vector<std::unique_ptr<ROMol>> p_templates;

 public:
  SCSRMol() {};
  SCSRMol(const SCSRMol &other) = delete;
  SCSRMol(SCSRMol &&other) noexcept = delete;
  SCSRMol &operator=(SCSRMol &&other) noexcept = delete;

  SCSRMol &operator=(const SCSRMol &) = delete;  // disable assignment
  ~SCSRMol() {}

  void addTemplate(std::unique_ptr<ROMol> templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    p_templates.push_back(std::move(templateMol));
  }

  unsigned int getTemplateCount() const { return p_templates.size(); }

  ROMol *getTemplate(unsigned int index) { return p_templates[index].get(); };
  const ROMol *getTemplate(unsigned int index) const {
    return p_templates[index].get();
  };

  const ROMol *getMol() const { return p_mol.get(); }

  RWMol *getMol() { return p_mol.get(); }

  void setMol(std::unique_ptr<RWMol> mol) {
    PRECONDITION(mol, "bad molecule");
    p_mol = std::move(mol);
  }
};

typedef boost::shared_ptr<SCSRMol> SCSRMOL_SPTR;

}  // namespace RDKit

#endif
