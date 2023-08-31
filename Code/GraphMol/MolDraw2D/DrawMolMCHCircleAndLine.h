//
//  Copyright (C) 2023 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This is the original multi-coloured highlighting code, which is done
// by multi-coloured arcs and bars.

#ifndef DRAWMOLMCHCIRCLEANDLINE_H
#define DRAWMOLMCHCIRCLEANDLINE_H

#include <GraphMol/MolDraw2D/DrawMolMCH.h>

namespace RDKit {
namespace MolDraw2D_detail {

class DrawMolMCHCircleAndLine : public DrawMolMCH {
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
  DrawMolMCHCircleAndLine(
      const ROMol &mol, const std::string &legend, int width, int height,
      MolDrawOptions &drawOptions, DrawText &textDrawer,
      const std::map<int, std::vector<DrawColour>> &highlight_atom_map,
      const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
      const std::map<int, double> &highlight_radii,
      const std::map<int, int> &highlight_linewidth_multipliers,
      int confId = -1);
  DrawMolMCHCircleAndLine(const DrawMol &) = delete;
  DrawMolMCHCircleAndLine(DrawMol &&) = delete;
  DrawMolMCHCircleAndLine &operator=(const DrawMol &) = delete;
  DrawMolMCHCircleAndLine &operator=(DrawMol &&) = delete;

  void extractHighlights(double scale) override;
  void extractMCHighlights() override;

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
};
}  // namespace MolDraw2D_detail
}  // namespace RDKit
#endif  // DRAWMOLMCHCIRCLEANDLINE_H
