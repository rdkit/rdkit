//
//  Copyright (C) 2003-2024 Tad Hurst and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file RWMol.h

  \brief Defines the editable molecule class \c RWMol

*/

#include <RDGeneral/export.h>

#ifndef RD_SCSRMOL_H
#define RD_SCSRMOL_H

// our stuff
#include "ROMol.h"
#include "RingInfo.h"

namespace RDKit {

//! SCSRMol contains an Self-contained Structure Representation (SCSR) of a
//! macro-molecule

class RDKIT_GRAPHMOL_EXPORT SCSRMol {
 private:
  std::unique_ptr<ROMol> p_mol;
  std::vector<std::unique_ptr<ROMol>> p_templates;

 public:
  SCSRMol();
  SCSRMol(const SCSRMol &other);
  // SCSRMol &operator=(const SCSRMol &);
  SCSRMol(SCSRMol &&other) noexcept;
  SCSRMol &operator=(SCSRMol &&other) noexcept;

  SCSRMol &operator=(const SCSRMol &) = delete;  // disable assignment
  virtual void destroyScsr();
  virtual ~SCSRMol() { destroyScsr(); }

  void addTemplate(std::unique_ptr<ROMol> templateMol) {
    PRECONDITION(templateMol, "bad template molecule");
    p_templates.push_back(std::move(templateMol));
  }

  unsigned int getTemplateCount() const { return p_templates.size(); }

  ROMol *getTemplate(unsigned int index) const {
    return p_templates[index].get();
  };

  const ROMol *getMol() const { return p_mol.get(); }

  ROMol *getMol() { return p_mol.get(); }
};

}  // namespace RDKit
#endif