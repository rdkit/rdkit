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
#include <GraphMol/MolDraw2D/MolDraw2D.h>

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace MolDraw2DUtils {

namespace {
bool isAtomCandForChiralH(const RWMol &mol, const Atom *atom) {
  // conditions for needing a chiral H:
  //   - stereochem specified
  //   - in at least two rings
  if ((!mol.getRingInfo()->isInitialized() ||
       mol.getRingInfo()->numAtomRings(atom->getIdx()) > 1) &&
      (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW ||
       atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW)) {
    return true;
  }
  return false;
}
}  // end of anonymous namespace

void prepareMolForDrawing(RWMol &mol, bool kekulize, bool addChiralHs,
                          bool wedgeBonds, bool forceCoords) {
  if (kekulize) {
    MolOps::Kekulize(mol, false);  // kekulize, but keep the aromatic flags!
  }
  if (addChiralHs) {
    std::vector<unsigned int> chiralAts;
    for (RWMol::AtomIterator atIt = mol.beginAtoms(); atIt != mol.endAtoms();
         ++atIt) {
      if (isAtomCandForChiralH(mol, *atIt)) {
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

void updateDrawerParamsFromJSON(MolDraw2D &drawer, const char *json) {
  PRECONDITION(json, "no parameter string");
  updateDrawerParamsFromJSON(drawer, std::string(json));
};
#define PT_OPT_GET(opt) opts.opt = pt.get(#opt, opts.opt)

void updateDrawerParamsFromJSON(MolDraw2D &drawer, const std::string &json) {
  if (json == "") return;
  std::istringstream ss;
  ss.str(json);
  MolDrawOptions &opts = drawer.drawOptions();
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  PT_OPT_GET(atomLabelDeuteriumTritium);
  PT_OPT_GET(dummiesAreAttachments);
  //
  //
  // opts.atomLabelDeuteriumTritium =
  //     pt.get("atomLabelDeuteriumTritium", opts.atomLabelDeuteriumTritium);
}

}  // end of MolDraw2DUtils namespace
}  // end of RDKit namespace
