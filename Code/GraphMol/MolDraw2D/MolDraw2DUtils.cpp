//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

namespace RDKit {
namespace MolDraw2DUtils {
void prepareMolForDrawing(RWMol &mol, bool kekulize, bool addChiralHs,
                          bool wedgeBonds, bool forceCoords) {
  if (kekulize) {
    MolOps::Kekulize(mol, false);  // kekulize, but keep the aromatic flags!
  }
  if (addChiralHs) {
    std::vector<unsigned int> chiralAts;
    for (RWMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
         ++atIt) {
      if ((*atIt)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW ||
          (*atIt)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW) {
        chiralAts.push_back((*atIt)->getIdx());
      }
    }
    if (chiralAts.size()) {
      bool addCoords = false;
      if (!forceCoords && mol.getNumConformers()) addCoords = true;
      MolOps::addHs(mol, false, addCoords, &chiralAts);
    }
  }
  if (forceCoords || !mol.getNumConformers()) {
    // compute 2D coordinates in a standard orientation:
    const bool canonOrient = true;
    RDDepict::compute2DCoords(mol, NULL, canonOrient);
  }
  if (wedgeBonds) {
    WedgeMolBonds(mol, &mol.getConformer());
  }
}
}
}
