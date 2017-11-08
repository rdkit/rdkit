//
//  Copyright (C) 2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>

#include <fstream>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <RDGeneral/RDLog.h>

#include "coordgenlibs/sketcherMinimizer.h"

using namespace RDKit;

void addCoordsWithCoordGen(ROMol& mol) {
  sketcherMinimizer minimizer;
  sketcherMinimizerMolecule* min_mol = new sketcherMinimizerMolecule();

  std::vector<sketcherMinimizerAtom*> ats;
  for (ROMol::AtomIterator atit = mol.beginAtoms(); atit != mol.endAtoms();
       ++atit) {
    sketcherMinimizerAtom* atom = new sketcherMinimizerAtom();
    atom->atomicNumber = (*atit)->getAtomicNum();
    atom->charge = (*atit)->getFormalCharge();
    ats.push_back(atom);
    min_mol->_atoms.push_back(atom);
  }
  std::vector<sketcherMinimizerBond*> bnds;
  for (ROMol::BondIterator bndit = mol.beginBonds(); bndit != mol.endBonds();
       ++bndit) {
    Bond* obnd = *bndit;
    sketcherMinimizerBond* bnd = new sketcherMinimizerBond();
    bnd->startAtom = ats[obnd->getBeginAtomIdx()];
    bnd->endAtom = ats[obnd->getEndAtomIdx()];
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

    bnds.push_back(bnd);
    min_mol->_bonds.push_back(bnd);
  }

  minimizer.initialize(min_mol);
  minimizer.runGenerateCoordinates();
  Conformer* conf = new Conformer(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    conf->setAtomPos(
        i, RDGeom::Point3D(minimizer._atoms[i]->coordinates.x(),
                           minimizer._atoms[i]->coordinates.y(), 0.0));
    // std::cerr << atom->coordinates << std::endl;
  }
  conf->set3D(false);
  mol.clearConformers();
  mol.addConformer(conf, true);
  std::string mb = MolToMolBlock(mol);
  std::cerr << mb << std::endl;
}

void test1() {
  BOOST_LOG(rdInfoLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << "test1: basics" << std::endl;
  {
    // ROMol* m = SmilesToMol("c1ccncc1");

    ROMol* m = SmilesToMol("CC(C)C");
    TEST_ASSERT(m);
    addCoordsWithCoordGen(*m);
    delete m;
  }

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}
int main(int argc, char* argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();
  test1();
}
