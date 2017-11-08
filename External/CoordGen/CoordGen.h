//
//  Copyright (C) 2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include "coordgenlibs/sketcherMinimizer.h"

namespace RDKit {
template <typename T>
void addCoordsWithCoordGen(T& mol) {
  sketcherMinimizer minimizer;
  auto min_mol = new sketcherMinimizerMolecule();

  std::map<unsigned int, sketcherMinimizerAtom*> atomMap;
  for (auto atit = mol.beginAtoms(); atit != mol.endAtoms(); ++atit) {
    auto oatom = *atit;
    auto atom = min_mol->addNewAtom();
    atom->molecule = min_mol;  // seems like this should be in addNewAtom()
    atom->atomicNumber = oatom->getAtomicNum();
    atom->charge = oatom->getFormalCharge();
    atomMap[oatom->getIdx()] = atom;
  }
  for (auto bndit = mol.beginBonds(); bndit != mol.endBonds(); ++bndit) {
    auto obnd = *bndit;
    auto bnd = min_mol->addNewBond(atomMap[obnd->getBeginAtomIdx()],
                                   atomMap[obnd->getEndAtomIdx()]);
    // FIX: This is no doubt wrong
    switch (obnd->getBondType()) {
      case Bond::SINGLE:
        bnd->bondOrder = 1;
        break;
      case Bond::DOUBLE:
        bnd->bondOrder = 2;
        break;
      case Bond::TRIPLE:
        bnd->bondOrder = 3;
        break;
      case Bond::AROMATIC:
        bnd->bondOrder = 1;
        break;
      default:
        BOOST_LOG(rdWarningLog) << "unrecognized bond type";
    }
  }

  minimizer.initialize(min_mol);
  minimizer.runGenerateCoordinates();
  auto conf = new Conformer(mol.getNumAtoms());
  for (auto i = 0; i < mol.getNumAtoms(); ++i) {
    conf->setAtomPos(i, RDGeom::Point3D(atomMap[i]->coordinates.x(),
                                        atomMap[i]->coordinates.y(), 0.0));
    // std::cerr << atom->coordinates << std::endl;
  }
  conf->set3D(false);
  mol.clearConformers();
  mol.addConformer(conf, true);
}
}  // end of namespace RDKit
