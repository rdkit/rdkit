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
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
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

void get_colour_option(boost::property_tree::ptree *pt, const char *pnm,
                       DrawColour &colour) {
  PRECONDITION(pnm && strlen(pnm), "bad property name");
  if (pt->find(pnm) == pt->not_found()) return;

  boost::property_tree::ptree::const_iterator itm = pt->get_child(pnm).begin();
  colour.get<0>() = itm->second.get_value<float>();
  ++itm;
  colour.get<1>() = itm->second.get_value<float>();
  ++itm;
  colour.get<2>() = itm->second.get_value<float>();
  ++itm;
}

void updateDrawerParamsFromJSON(MolDraw2D &drawer, const std::string &json) {
  if (json == "") return;
  std::istringstream ss;
  ss.str(json);
  MolDrawOptions &opts = drawer.drawOptions();
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  PT_OPT_GET(atomLabelDeuteriumTritium);
  PT_OPT_GET(dummiesAreAttachments);
  PT_OPT_GET(circleAtoms);
  PT_OPT_GET(continuousHighlight);
  PT_OPT_GET(flagCloseContactsDist);
  PT_OPT_GET(includeAtomTags);
  PT_OPT_GET(clearBackground);
  PT_OPT_GET(legendFontSize);
  PT_OPT_GET(multipleBondOffset);
  PT_OPT_GET(padding);
  PT_OPT_GET(additionalAtomLabelPadding);
  get_colour_option(&pt, "highlightColour", opts.highlightColour);
  get_colour_option(&pt, "backgroundColour", opts.backgroundColour);
  get_colour_option(&pt, "legendColour", opts.legendColour);
  if (pt.find("atomLabels") != pt.not_found()) {
    BOOST_FOREACH (boost::property_tree::ptree::value_type const &item,
                   pt.get_child("atomLabels")) {
      opts.atomLabels[boost::lexical_cast<int>(item.first)] =
          item.second.get_value<std::string>();
    }
  }
}

}  // end of MolDraw2DUtils namespace
}  // end of RDKit namespace
