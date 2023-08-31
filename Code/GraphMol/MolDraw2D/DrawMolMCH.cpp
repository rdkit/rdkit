//
//  Copyright (C) 2021-2023 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#include <GraphMol/RWMol.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/MolDraw2D/DrawMolMCH.h>

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************
DrawMolMCH::DrawMolMCH(
    const ROMol &mol, const std::string &legend, int width, int height,
    MolDrawOptions &drawOptions, DrawText &textDrawer,
    const std::map<int, std::vector<DrawColour>> &highlight_atom_map,
    const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
    const std::map<int, double> &highlight_radii,
    const std::map<int, int> &highlight_linewidth_multipliers, int confId)
    : DrawMol(mol, legend, width, height, drawOptions, textDrawer, nullptr,
              nullptr, nullptr, nullptr, nullptr, &highlight_radii, confId),
      mcHighlightAtomMap_(highlight_atom_map),
      mcHighlightBondMap_(highlight_bond_map),
      highlightLinewidthMultipliers_(highlight_linewidth_multipliers) {}

// ****************************************************************************
void DrawMolMCH::extractHighlights(double scale) {
  DrawMol::extractHighlights(scale);
  extractMCHighlights();
}

// ****************************************************************************
void DrawMolMCH::getAtomRadius(unsigned int atomIdx, double &xradius,
                               double &yradius) const {
  xradius = drawOptions_.highlightRadius;
  yradius = xradius;
  if (highlightRadii_.find(atomIdx) != highlightRadii_.end()) {
    xradius = highlightRadii_.find(atomIdx)->second;
    yradius = xradius;
  }
}
}  // namespace MolDraw2D_detail
}  // namespace RDKit