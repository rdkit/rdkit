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
  std::unique_ptr<ROMol> mol;
  std::vector<std::unique_ptr<ROMol>> templates;

 public:
  SCSRMol();
  SCSRMol(const SCSRMol &other);
  SCSRMol &operator=(const SCSRMol &);
  SCSRMol(SCSRMol &&other) noexcept;
  SCSRMol &operator=(SCSRMol &&other) noexcept;

  virtual ~SCSRMol() {}

  void addTemplate(ROMol *templateMol) {
    templates.push_back(std::unique_ptr<ROMol>(templateMol));
  }

  unsigned int getTemplateCount() const { return templates.size(); }

  const ROMol *getTemplate(unsigned int index) const {
    return templates[index].get();
  };

  const ROMol *getMol() const { return mol.get(); }

  ROMol *getMol() { return mol.get(); }
};

typedef boost::shared_ptr<SCSRMol> SCSRMOL_SPTR;

}  // namespace RDKit
#endif