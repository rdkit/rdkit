//
//  Copyright (c) 2023, Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
#include <RDGeneral/export.h>
#ifndef XQMOL_H_MAY2023
#define XQMOL_H_MAY2023

#include <variant>
#include <memory>
#include <string>
#include <vector>

namespace RDKit {
class RWMol;
class MolBundle;
class TautomerQuery;

struct RDKIT_SUBSTRUCTMATCH_EXPORT ExtendedQueryMol {
  enum ExtendedQueryMolTypes : unsigned char {
    XQM_MOL = 1,
    XQM_MOLBUNDLE = 2,
    XQM_TAUTOMERQUERY = 3,
    XQM_TAUTOMERBUNDLE = 4
  };
  using ContainedType =
      std::variant<std::unique_ptr<RWMol>, std::unique_ptr<MolBundle>,
                   std::unique_ptr<TautomerQuery>,
                   std::vector<std::unique_ptr<TautomerQuery>>>;

  ContainedType xqmol;
  std::string toBinary() const;
};

}  // namespace RDKit
#endif
