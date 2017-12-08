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
namespace CoordGen {

struct CoordGenParams {
  RDGeom::INT_POINT2D_MAP coordMap;
};

template <typename T>
void addCoords(T& mol, const CoordGenParams* params = nullptr) {
  sketcherMinimizer minimizer;
  auto min_mol = new sketcherMinimizerMolecule();

  std::vector<sketcherMinimizerAtom*> ats(mol.getNumAtoms());
  for (auto atit = mol.beginAtoms(); atit != mol.endAtoms(); ++atit) {
    auto oatom = *atit;
    auto atom = min_mol->addNewAtom();
    atom->molecule = min_mol;  // seems like this should be in addNewAtom()
    atom->atomicNumber = oatom->getAtomicNum();
    atom->charge = oatom->getFormalCharge();
    if (params &&
        params->coordMap.find(oatom->getIdx()) != params->coordMap.end()) {
      atom->constrained = true;
      atom->fixed = true;
      const RDGeom::Point2D& coords =
          params->coordMap.find(oatom->getIdx())->second;
      atom->templateCoordinates = sketcherMinimizerPointF(coords.x, coords.y);
    }
    ats[oatom->getIdx()] = atom;
  }

  std::vector<sketcherMinimizerBond*> bnds(mol.getNumBonds());
  for (auto bndit = mol.beginBonds(); bndit != mol.endBonds(); ++bndit) {
    auto obnd = *bndit;
    auto bnd = min_mol->addNewBond(ats[obnd->getBeginAtomIdx()],
                                   ats[obnd->getEndAtomIdx()]);
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
    bnds[obnd->getIdx()] = bnd;
  }

  // setup double bond stereo
  min_mol->assignBondsAndNeighbors(ats, bnds);
  for (auto bndit = mol.beginBonds(); bndit != mol.endBonds(); ++bndit) {
    auto obnd = *bndit;
    if (obnd->getBondType() != Bond::DOUBLE ||
        obnd->getStereo() <= Bond::STEREOANY ||
        obnd->getStereo() > Bond::STEREOTRANS)
      continue;

    sketcherMinimizerBondStereoInfo sinfo;
    sinfo.atom1 = ats[obnd->getStereoAtoms()[0]];
    sinfo.atom2 = ats[obnd->getStereoAtoms()[1]];
    sinfo.stereo = (obnd->getStereo() == Bond::STEREOZ ||
                    obnd->getStereo() == Bond::STEREOCIS)
                       ? sketcherMinimizerBondStereoInfo::cis
                       : sketcherMinimizerBondStereoInfo::trans;
    auto bnd = bnds[obnd->getIdx()];
    bnd->setStereoChemistry(sinfo);
    bnd->setAbsoluteStereoFromStereoInfo();
  }

  minimizer.initialize(min_mol);
  minimizer.runGenerateCoordinates();
  auto conf = new Conformer(mol.getNumAtoms());
  for (size_t i = 0; i < mol.getNumAtoms(); ++i) {
    conf->setAtomPos(i, RDGeom::Point3D(ats[i]->coordinates.x(),
                                        ats[i]->coordinates.y(), 0.0));
    // std::cerr << atom->coordinates << std::endl;
  }
  conf->set3D(false);
  mol.clearConformers();
  mol.addConformer(conf, true);
}
}  // end of namespace CoordGen
}  // end of namespace RDKit
