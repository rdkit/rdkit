//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_MOLENUMERATOR_H
#define RDKIT_MOLENUMERATOR_H

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>

#include <vector>
#include <map>
#include <string>

namespace RDKit {
namespace MolEnumerator {

//! abstract base class
class RDKIT_MOLENUMERATOR_EXPORT MolEnumeratorOp {
 private:
  ROMol d_mol;

 public:
  MolEnumeratorOp(const ROMol &mol) : d_mol(mol){};
  MolEnumeratorOp(const MolEnumeratorOp &other) : d_mol(other.d_mol){};
  MolEnumeratorOp &operator=(const MolEnumeratorOp &other) {
    if (&other != this) {
      d_mol = ROMol(other.d_mol);
    }
    return *this;
  };
  ~MolEnumeratorOp() = default;
};

struct RDKIT_MOLENUMERATOR_EXPORT MolEnumeratorParams {
  bool sanitize = false;
  size_t maxToEnumerate = 0;
  bool doRandom = false;
  int randomSeed = -1;
  std::vector<MolEnumeratorOp> operations;
};

}  // namespace MolEnumerator
}  // namespace RDKit

#endif