//
//  Copyright (C) 2021-2022 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//
// This class, derived from DrawMol, does multi-coloured highlights.
// Like DrawMol, it is not part of the public API and should only be
// used by MolDraw2D.

#ifndef RDKIT_DRAWMOLMCH_H
#define RDKIT_DRAWMOLMCH_H

#include <GraphMol/MolDraw2D/AtomSymbol.h>
#include <GraphMol/MolDraw2D/DrawMol.h>

namespace RDKit {
namespace MolDraw2D_detail {

class DrawMolMCH : public DrawMol {
 public:
  DrawMolMCH(const ROMol &mol, const std::string &legend, int width, int height,
             MolDrawOptions &drawOptions, DrawText &textDrawer,
             const std::map<int, std::vector<DrawColour>> &highlight_atom_map,
             const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
             const std::map<int, double> &highlight_radii,
             const std::map<int, int> &highlight_linewidth_multipliers,
             int confId = -1);
  DrawMolMCH(const DrawMol &) = delete;
  DrawMolMCH(DrawMol &&) = delete;
  DrawMolMCH &operator=(const DrawMol &) = delete;
  DrawMolMCH &operator=(DrawMol &&) = delete;

  void extractHighlights() override;
  void extractMCHighlights();
  void makeBondHighlights();
  void makeAtomHighlights();
  void adjustLineEndForHighlight(int at_idx, Point2D p1, Point2D &p2) const;
  void calcSymbolEllipse(unsigned int atomIdx, Point2D &centre, double &xradius,
                         double &yradius) const;

  const std::map<int, std::vector<DrawColour>> mcHighlightAtomMap_;
  const std::map<int, std::vector<DrawColour>> &mcHighlightBondMap_;
  const std::map<int, int> &highlightLinewidthMultipliers_;
};

}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif  // RDKIT_DRAWMOLMCH_H
