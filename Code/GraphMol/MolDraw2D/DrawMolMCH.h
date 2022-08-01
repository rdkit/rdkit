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
  /*!
   Make a DrawMol that does multi-coloured highlights.  Both atoms and
   bonds can have more than 1 highlight.  There's no maximum although more
   than 4 is very cluttered.
   \param mol             : the molecule to draw
   \param legend          : the legend (to be drawn under the molecule)
   \param width           : width (in pixels) of the rendering
   set this to -1 to have the canvas size set automatically
   \param height          : height (in pixels) of the rendering
   set this to -1 to have the canvas size set automatically
   \param drawOptions    : a MolDrawOptions object from the owning MolDraw2D
   \param textDrawer     : a DrawText object from the owning MolDraw2D
   \param highlight_atom_map : indexed on atom idx, the colours to be used to
   highlight atoms.  Not all atoms need to be mentioned.
   \param highlight_bond_map : indexed on bond idx, the colours to be used to
   highlight bonds.  Not all bonds need to be mentioned.
   \param highlightRadii : map from atomId -> radius (in molecule coordinates)
   for the radii of atomic highlights. If not provided for an atom, the default
   value from \c drawOptions() will be used.
   \param highlight_linewidth_multipliers : map from bondId -> int, to change
   the thickness of the lines used for the highlighting.
   \param confId         : (optional) conformer ID to be used for atomic
   coordinates
   */
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

  void extractHighlights(double scale) override;
  void extractMCHighlights();
  void makeBondHighlights(
      std::vector<std::unique_ptr<DrawShape>> &bondHighlights);
  void makeAtomHighlights(
      std::vector<std::unique_ptr<DrawShape>> &atomHighlights);
  // adjust p2 so that the line from p1 to p2 stops where it intersects
  // the label ellipse for at_idx.
  void adjustLineEndForHighlight(int at_idx, Point2D p1, Point2D &p2) const;
  void calcSymbolEllipse(unsigned int atomIdx, Point2D &centre, double &xradius,
                         double &yradius) const;
  void fixHighlightJoinProblems(
      std::vector<std::unique_ptr<DrawShape>> &atomHighlights,
      std::vector<std::unique_ptr<DrawShape>> &bondHighlights);

  const std::map<int, std::vector<DrawColour>> mcHighlightAtomMap_;
  const std::map<int, std::vector<DrawColour>> &mcHighlightBondMap_;
  const std::map<int, int> &highlightLinewidthMultipliers_;
};

// adjust p2 so that the line from p1 to p2 stops where it intersects
// the ellipse.
void adjustLineEndForEllipse(const Point2D &centre, double xradius,
                             double yradius, Point2D p1, Point2D &p2);
}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif  // RDKIT_DRAWMOLMCH_H
