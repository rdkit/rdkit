

#ifndef DRAWMOLMCHLASSO_H
#define DRAWMOLMCHLASSO_H

#include <GraphMol/MolDraw2D/DrawMolMCH.h>

namespace RDKit {
namespace MolDraw2D_detail {

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
 \param highlight_bond_map : ignored in this molecule representation.
 \param highlightRadii : ignored in this molecule representation[.
 \param highlight_linewidth_multipliers : map from bondId -> int, to change
 the thickness of the lines used for the highlighting.
 \param confId         : (optional) conformer ID to be used for atomic
 coordinates
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
                              std::vector<std::vector<int>> &colourAtoms) const;
  void drawLasso(size_t colNum, const DrawColour &col,
                 const std::vector<int> &colAtoms);
  void extractAtomArcs(size_t colNum, const DrawColour &col,
                       const std::vector<int> &colAtoms,
                       std::vector<std::unique_ptr<DrawShapeArc>> &arcs) const;
  void extractBondLines(
      size_t colNum, const DrawColour &col, const std::vector<int> &colAtoms,
      std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) const;
  void fixArcsAndLines(
      std::vector<std::unique_ptr<DrawShapeArc>> &arcs,
      std::vector<std::unique_ptr<DrawShapeSimpleLine>> &lines) const;
};
}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif  // DRAWMOLMCHLASSO_H
