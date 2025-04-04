//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>

namespace RDKit {
namespace MolOps {

ROMol *renumberAtoms(const ROMol &mol,
                     const std::vector<unsigned int> &newOrder) {
  unsigned int nAts = mol.getNumAtoms();
  PRECONDITION(newOrder.size() == nAts, "bad newOrder size");

  std::vector<unsigned int> revOrder(nAts);
  for (unsigned int nIdx = 0; nIdx < nAts; ++nIdx) {
    unsigned int oIdx = newOrder[nIdx];
    if (oIdx > nAts) {
      throw ValueErrorException("idx value exceeds numAtoms");
    }
    revOrder[oIdx] = nIdx;
  }

  // ------
  // newOrder[i] : which atom should be in position i of the new mol
  // revOrder[i] : where atom i of the original mol landed in the new mol
  auto *res = new RWMol();

  // copy over the atoms:
  for (unsigned int nIdx = 0; nIdx < nAts; ++nIdx) {
    unsigned int oIdx = newOrder[nIdx];
    const Atom *oAtom = mol.getAtomWithIdx(oIdx);
    Atom *nAtom = oAtom->copy();
    res->addAtom(nAtom, false, true);

    // take care of atom-numbering-dependent properties:
    INT_VECT nAtoms;
    if (nAtom->getPropIfPresent(common_properties::_ringStereoAtoms, nAtoms)) {
      // FIX: ought to be able to avoid this copy.
      for (auto &val : nAtoms) {
        if (val < 0) {
          val = -1 * (revOrder[(-val - 1)] + 1);
        } else {
          val = revOrder[val - 1] + 1;
        }
      }
      nAtom->setProp(common_properties::_ringStereoAtoms, nAtoms, true);
    }

    unsigned int otherAtom;
    if (nAtom->getPropIfPresent(common_properties::_ringStereoOtherAtom,
                                otherAtom)) {
      otherAtom = revOrder[otherAtom];
      nAtom->setProp(common_properties::_ringStereoOtherAtom, otherAtom, true);
    }
  }

  // now the bonds:
  for (ROMol::ConstBondIterator bi = mol.beginBonds(); bi != mol.endBonds();
       ++bi) {
    const Bond *oBond = (*bi);
    Bond *nBond = oBond->copy();
    nBond->setBeginAtomIdx(revOrder[oBond->getBeginAtomIdx()]);
    nBond->setEndAtomIdx(revOrder[oBond->getEndAtomIdx()]);
    res->addBond(nBond, true);
    // take care of atom-numbering-dependent properties:
    for (auto &idx : nBond->getStereoAtoms()) {
      idx = revOrder[idx];
    }
  }

  // Conformers:
  for (auto oConf = mol.beginConformers(); oConf != mol.endConformers();
       ++oConf) {
    auto *nConf = new Conformer(nAts);
    for (unsigned int i = 0; i < nAts; ++i) {
      nConf->setAtomPos(i, (*oConf)->getAtomPos(newOrder[i]));
    }
    nConf->setId((*oConf)->getId());
    nConf->set3D((*oConf)->is3D());
    res->addConformer(nConf);
  }

  // update the ring info:
  const RingInfo *oRings = mol.getRingInfo();
  if (oRings && oRings->isInitialized()) {
    RingInfo *nRings = res->getRingInfo();
    nRings->reset();
    nRings->initialize();
    for (unsigned int i = 0; i < oRings->numRings(); ++i) {
      const INT_VECT &oRing = oRings->atomRings()[i];
      INT_VECT nRing(oRing.size());
      for (unsigned int j = 0; j < oRing.size(); ++j) {
        nRing[j] = revOrder[oRing[j]];
      }
      nRings->addRing(nRing, oRings->bondRings()[i]);
    }
  }

  if (mol.getStereoGroups().size()) {
    std::vector<StereoGroup> nsgs;
    nsgs.reserve(mol.getStereoGroups().size());
    for (const auto &osg : mol.getStereoGroups()) {
      std::vector<Atom *> ats;
      std::vector<Bond *> bds;
      ats.reserve(osg.getAtoms().size());
      bds.reserve(osg.getBonds().size());
      for (const auto aptr : osg.getAtoms()) {
        ats.push_back(res->getAtomWithIdx(revOrder[aptr->getIdx()]));
      }
      for (const auto bptr : osg.getBonds()) {
        bds.push_back(
            res->getBondWithIdx(bptr->getIdx()));  // bonds do not change order
      }

      nsgs.emplace_back(osg.getGroupType(), ats, bds, osg.getReadId());
    }
    res->setStereoGroups(std::move(nsgs));
  }

  return dynamic_cast<ROMol *>(res);
}  // namespace

};  // end of namespace MolOps
};  // end of namespace RDKit
