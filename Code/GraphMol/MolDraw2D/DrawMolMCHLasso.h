//
//  Copyright (C) 2023 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#ifndef DRAWMOLMCHLASSO_H
#define DRAWMOLMCHLASSO_H

#include <GraphMol/MolDraw2D/DrawMolMCH.h>

namespace RDKit {
namespace MolDraw2D_detail {

// Struct to hold a pair of lines that make up part of the
// lasso.  Each atom in the lasso will have one for each
// of its bonds that are in the lasso.  The angles are
// subtended by the line ends at the atom centres relative
// to the x axis.  The radius is also the distance from
// the nearest end of each line ot the atom.
struct LinePair {
  DrawShapeSimpleLine *line1;
  DrawShapeSimpleLine *line2;
  double angle1;
  double angle2;
  double radius;
  int atom;
};

class DrawMolMCHLasso : public DrawMolMCH {
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
 \param highlight_bond_map : If this is empty, then the lassos will be between
                             all bonded atoms in the highlight_atom_map, in
                             the appropriate colours.  If it is given, the lasso
                             will be in more than one part if there is a bond
                             between 2 highlighted atoms that isn't in the
                             highlight_bond_map.  The colours in
                             highlight_bond_map will be NOT used for the lines
                             between the atoms - the lasso will be one colour,
                             taken from the atoms.
 \param highlightRadii : ignored in this molecule representation.
 \param highlight_linewidth_multipliers : ignored in this molecule
 representation. \param confId         : (optional) conformer ID to be used for
 atomic coordinates
 */
  DrawMolMCHLasso(
      const ROMol &mol, const std::string &legend, int width, int height,
      MolDrawOptions &drawOptions, DrawText &textDrawer,
      const std::map<int, std::vector<DrawColour>> &highlight_atom_map,
      const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
      const std::map<int, double> &highlight_radii,
      const std::map<int, int> &highlight_linewidth_multipliers,
      int confId = -1);
  DrawMolMCHLasso(const DrawMol &) = delete;
  DrawMolMCHLasso(DrawMol &&) = delete;
  DrawMolMCHLasso &operator=(const DrawMol &) = delete;
  DrawMolMCHLasso &operator=(DrawMol &&) = delete;

  void extractHighlights(double scale) override;
  void extractMCHighlights() override;
  // turn mcHighlightAtomMap_ into a mapping of vectors of atoms of a particular
  // colour, indexed by the colour.
  void extractAtomColourLists(std::vector<DrawColour> &colours,
                              std::vector<std::vector<int>> &colourAtoms,
                              std::vector<std::vector<int>> &colourLists) const;
  void drawLasso(size_t lassoNum, const DrawColour &col,
                 const std::vector<int> &colAtoms);
  void extractBondLines(
      size_t lassoNum, const DrawColour &col, const std::vector<int> &colAtoms,
      std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines,
      std::vector<std::vector<LinePair>> &atomLines) const;
  void extractAtomArcs(std::vector<std::vector<LinePair>> &atomLines,
                       std::vector<std::unique_ptr<DrawShapeArc>> &arcs) const;
  void addNoLineArcs(
      const std::vector<int> &colAtoms, size_t lassoNum,
      const RDKit::DrawColour &col,
      const std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines,
      std::vector<std::unique_ptr<DrawShapeArc>> &arcs) const;
  void makeIntersectingArcs(
      const std::vector<std::vector<unsigned int>> &intersects, int lassoNum,
      const RDKit::DrawColour &col,
      std::vector<std::unique_ptr<DrawShapeArc>> &currArcs,
      std::vector<DrawShapeArc *> &arcs) const;
  // Put the start points of each pair of atoms in order so that going in
  // an anti-clockwise direction the bond that the lines straddle is in
  // between them such that the arc between the first and second is always
  // the bit that isn't drawn.  And then order the pairs so that the bonds
  // are moving round in an anti-clockwise direction.
  void orderAtomLines(std::vector<std::vector<LinePair>> &atomLines) const;
};
}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif  // DRAWMOLMCHLASSO_H
