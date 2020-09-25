//
//  Copyright (C) 2020 Manan Goel
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <ForceField/ForceField.h>
#include <ForceField/ANI/AtomicContrib.h>
#include "Builder.h"

namespace RDKit {
namespace ANI {
using namespace ForceFields::ANI;
namespace Tools {
void addANIContribs(const ROMol &mol, ForceFields::ForceField *field,
                    std::string model, int confId) {
  PRECONDITION(field, "bad ForceField");
  // get atomic number for each atom (speciesVec)
  auto numAtoms = mol.getNumAtoms();
  VectorXi speciesVec(numAtoms);
  for (unsigned int i = 0; i < numAtoms; i++) {
    auto atom = mol.getAtomWithIdx(i);
    speciesVec[i] = atom->getAtomicNum();
  }
  auto ac = new ForceFields::ANI::ANIAtomContrib(field, speciesVec,
                                                 model);
  field->contribs().push_back(ForceFields::ContribPtr(ac));
}
}  // namespace Tools

ForceFields::ForceField *constructForceField(ROMol &mol, std::string model,
                                             int confId) {
  auto *res = new ForceFields::ForceField();
  Conformer &conf = mol.getConformer(confId);
  for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    res->positions().push_back(&conf.getAtomPos(i));
  }
  Tools::addANIContribs(mol, res, model, confId);
  return res;
}

}  // namespace ANI
}  // namespace RDKit