//
//  Copyright (C) 2014-2021 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (AstraZeneca)
// 27th May 2014
//
// This class makes a 2D drawing of an RDKit molecule.
// It draws heavily on $RDBASE/GraphMol/MolDrawing/MolDrawing.h.
// One purpose of this is to make it easier to overlay annotations on top of
// the molecule drawing, which is difficult to do from the output of
// MolDrawing.h
// The class design philosophy echoes a standard one:
// a virtual base class defines the interface and does all
// the heavy lifting and concrete derived classes implement
// library-specific drawing code such as drawing lines, writing strings
// etc.

#include <RDGeneral/export.h>
#ifndef RDKITMOLDRAW2D_H
#define RDKITMOLDRAW2D_H

#include <vector>

#include <Geometry/point.h>
#include <Geometry/Transform2D.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

// ****************************************************************************
using RDGeom::Point2D;

namespace RDKit {

namespace MolDraw2D_detail {
class DrawMol;
class DrawText;
}  // namespace MolDraw2D_detail

//! MolDraw2D is the base class for doing 2D renderings of molecules
class RDKIT_MOLDRAW2D_EXPORT MolDraw2D {
 public:
  //! constructor for a particular size
  /*!
    \param width       : width (in pixels) of the rendering
    set this to -1 to have the canvas size set automatically
    \param height      : height (in pixels) of the rendering
    set this to -1 to have the canvas size set automatically
    \param panelWidth  : (optional) width (in pixels) of a single panel
    \param panelHeight : (optional) height (in pixels) of a single panel

    The \c panelWidth and \c panelHeight arguments are used to provide the
    sizes of the panels individual molecules are drawn in when
    \c drawMolecules() is called.
  */
  MolDraw2D(int width, int height, int panelWidth, int panelHeight);
  MolDraw2D(const MolDraw2D &rhs) = delete;
  MolDraw2D(MolDraw2D &&rhs) = delete;
  MolDraw2D &operator=(const MolDraw2D &rhs) = delete;
  MolDraw2D &operator=(MolDraw2D &&rhs) = delete;
  virtual ~MolDraw2D();

  //! draw a single molecule
  /*!
    \param mol             : the molecule to draw
    \param legend          : the legend (to be drawn under the molecule)
    \param highlight_atoms : (optional) vector of atom ids to highlight
    \param highlight_atoms : (optional) vector of bond ids to highlight
    \param highlight_atom_map   : (optional) map from atomId -> DrawColour
    providing the highlight colors. If not provided the default highlight colour
    from \c drawOptions() will be used.
    \param highlight_bond_map   : (optional) map from bondId -> DrawColour
    providing the highlight colors. If not provided the default highlight colour
    from \c drawOptions() will be used.
    \param highlight_radii : (optional) map from atomId -> radius (in molecule
    coordinates) for the radii of atomic highlights. If not provided the default
    value from \c drawOptions() will be used.
    \param confId          : (optional) conformer ID to be used for atomic
    coordinates

  */
  virtual void drawMolecule(
      const ROMol &mol, const std::string &legend,
      const std::vector<int> *highlight_atoms,
      const std::vector<int> *highlight_bonds,
      const std::map<int, DrawColour> *highlight_atom_map = nullptr,
      const std::map<int, DrawColour> *highlight_bond_map = nullptr,
      const std::map<int, double> *highlight_radii = nullptr, int confId = -1);

  //! \overload
  virtual void drawMolecule(
      const ROMol &mol, const std::vector<int> *highlight_atoms = nullptr,
      const std::map<int, DrawColour> *highlight_map = nullptr,
      const std::map<int, double> *highlight_radii = nullptr, int confId = -1);

  //! \overload
  virtual void drawMolecule(
      const ROMol &mol, const std::string &legend,
      const std::vector<int> *highlight_atoms = nullptr,
      const std::map<int, DrawColour> *highlight_map = nullptr,
      const std::map<int, double> *highlight_radii = nullptr, int confId = -1);

  //! \overload
  virtual void drawMolecule(
      const ROMol &mol, const std::vector<int> *highlight_atoms,
      const std::vector<int> *highlight_bonds,
      const std::map<int, DrawColour> *highlight_atom_map = nullptr,
      const std::map<int, DrawColour> *highlight_bond_map = nullptr,
      const std::map<int, double> *highlight_radii = nullptr, int confId = -1);

  //! draw molecule with multiple colours allowed per atom.
  /*!
    \param mol             : the molecule to draw
    \param legend          : the legend (to be drawn under the molecule)
    \param highlight_atom_map   : map from atomId -> DrawColours
    providing the highlight colours.
    \param highlight_bond_map   : map from bondId -> DrawColours
    providing the highlight colours.
    \param highlight_radii : map from atomId -> radius (in molecule
    coordinates) for the radii of atomic highlights. If not provided for an
    index, the default value from \c drawOptions() will be used.
    \param highlight_linewidth_multipliers : map from atomId -> int, used to
    vary the width the highlight lines.  Only active if
    drawOptions().fillHighlights is false.
    \param confId          : (optional) conformer ID to be used for atomic
    coordinates
  */
  virtual void drawMoleculeWithHighlights(
      const ROMol &mol, const std::string &legend,
      const std::map<int, std::vector<DrawColour>> &highlight_atom_map,
      const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
      const std::map<int, double> &highlight_radii,
      const std::map<int, int> &highlight_linewidth_multipliers,
      int confId = -1);

  //! draw multiple molecules in a grid
  /*!
    \param mols             : the molecules to draw
    \param legends          : (optional) the legends (to be drawn under the
    molecules)
    \param highlight_atoms  : (optional) vectors of atom ids to highlight
    \param highlight_atoms  : (optional) vectors of bond ids to highlight
    \param highlight_atom_map   : (optional) maps from atomId -> DrawColour
    providing the highlight colors. If not provided the default highlight colour
    from \c drawOptions() will be used.
    \param highlight_bond_map   : (optional) maps from bondId -> DrawColour
    providing the highlight colors. If not provided the default highlight colour
    from \c drawOptions() will be used.
    \param highlight_radii  : (optional) maps from atomId -> radius (in molecule
    coordinates) for the radii of atomic highlights. If not provided the default
    value from \c drawOptions() will be used.
    \param confId           : (optional) conformer IDs to be used for atomic
    coordinates

    The \c panelWidth and \c panelHeight values will be used to determine the
    number of rows and columns to be drawn. Theres not a lot of error checking
    here, so if you provide too many molecules for the number of panes things
    are likely to get screwed up.
    If the number of rows or columns ends up being <= 1, molecules will be
    being drawn in a single row/column.
  */
  virtual void drawMolecules(
      const std::vector<ROMol *> &mols,
      const std::vector<std::string> *legends = nullptr,
      const std::vector<std::vector<int>> *highlight_atoms = nullptr,
      const std::vector<std::vector<int>> *highlight_bonds = nullptr,
      const std::vector<std::map<int, DrawColour>> *highlight_atom_maps =
          nullptr,
      const std::vector<std::map<int, DrawColour>> *highlight_bond_maps =
          nullptr,
      const std::vector<std::map<int, double>> *highlight_radii = nullptr,
      const std::vector<int> *confIds = nullptr);

  //! draw a ChemicalReaction
  /*!
    \param rxn                 : the reaction to draw
    \param highlightByReactant : (optional) if this is set, atoms and bonds will
    be highlighted based on which reactant they come from. Atom map numbers
    will not be shown.
    \param highlightColorsReactants : (optional) provide a vector of colors for
    the
    reactant highlighting.
    \param confIds   : (optional) vector of confIds to use for rendering. These
    are numbered by reactants, then agents, then products.
  */
  virtual void drawReaction(
      const ChemicalReaction &rxn, bool highlightByReactant = false,
      const std::vector<DrawColour> *highlightColorsReactants = nullptr,
      const std::vector<int> *confIds = nullptr);

  //! clears the contents of the drawing
  virtual void clearDrawing() {
    if (needs_init_) {
      initDrawing();
      needs_init_ = false;
    }
  };
  //! draws a line from \c cds1 to \c cds2 using the current drawing style
  //! in atom coords.  If rawCoords is passed as true,
  //! the coordinates are used as is, if not they are assumed to be in
  //! the molecule coordinate frame and converted with getDrawCoords
  //! into canvas coords.
  virtual void drawLine(const Point2D &cds1, const Point2D &cds2,
                        bool rawCoords = false) = 0;
  //! draw a polygon.  Note that if fillPolys() returns false, it
  //! doesn't close the path.  If you want it to in that case, you
  //! do it explicitly yourself.  If rawCoords is passed as true,
  //! the coordinates are used as is, if not they are assumed to be in
  //! the molecule coordinate frame and converted with getDrawCoords
  //! into canvas coords.
  virtual void drawPolygon(const std::vector<Point2D> &cds,
                           bool rawCoords = false) = 0;
  //! @}

  //! A Whole bunch of drawing primitives.  They may be over-ridden
  //! by different renderers, or they may be implemented in terms of
  //! drawLine and drawPolygon above.  If rawCoords is passed as true,
  // the coordinates are used as is, if not they are assumed to be in
  // the molecule coordinate frame and converted with getDrawCoords
  // into canvas coords.
  //! draw a line where the ends are different colours
  virtual void drawLine(const Point2D &cds1, const Point2D &cds2,
                        const DrawColour &col1, const DrawColour &col2,
                        bool rawCoords = false);
  //! draw a triangle
  virtual void drawTriangle(const Point2D &cds1, const Point2D &cds2,
                            const Point2D &cds3, bool rawCoords = false);
  //! draw an ellipse
  virtual void drawEllipse(const Point2D &cds1, const Point2D &cds2,
                           bool rawCoords = false);
  // draw the arc of a circle between ang1 and ang2.  Note that 0 is
  // at 3 o-clock and 90 at 12 o'clock as you'd expect from your maths.
  // ang2 must be > ang1 - it won't draw backwards.  This is not enforced.
  // Angles in degrees.
  virtual void drawArc(const Point2D &centre, double radius, double ang1,
                       double ang2, bool rawCoords = false);
  // and a general ellipse form
  virtual void drawArc(const Point2D &centre, double xradius, double yradius,
                       double ang1, double ang2, bool rawCoords = false);
  //! draw a rectangle given two opposite corners
  virtual void drawRect(const Point2D &cds1, const Point2D &cds2,
                        bool rawCoords = false);
  //! draw a line indicating the presence of an attachment point (normally a
  //! squiggle line perpendicular to a bond)
  virtual void drawAttachmentLine(const Point2D &cds1, const Point2D &cds2,
                                  const DrawColour &col, double len = 1.0,
                                  unsigned int nSegments = 16,
                                  bool rawCoords = false);
  //! draw a wavy line like that used to indicate unknown stereochemistry
  virtual void drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                            const DrawColour &col1, const DrawColour &col2,
                            unsigned int nSegments = 16,
                            double vertOffset = 0.05, bool rawCoords = false);
  //! Draw an arrow with either lines or a filled head (when asPolygon is true)
  virtual void drawArrow(const Point2D &cds1, const Point2D &cds2,
                         bool asPolygon = false, double frac = 0.05,
                         double angle = M_PI / 6,
                         const DrawColour &col = DrawColour(0.0, 0.0, 0.0),
                         bool rawCoords = false);
  // draw a plus sign with lines at the given position.
  virtual void drawPlus(const Point2D &cds, int plusWidth,
                        const DrawColour &col, bool rawCoords = false);

  //! drawString centres the string on cds.
  virtual void drawString(const std::string &str, const Point2D &cds,
                          bool rawCoords = false);
  // unless the specific drawer over-rides this overload, it will just call
  // the first one.  SVG for one needs the alignment flag.
  virtual void drawString(const std::string &str, const Point2D &cds,
                          MolDraw2D_detail::TextAlignType align,
                          bool rawCoords = false);

  //! \name Transformations
  //! @{
  // transform a set of coords in the molecule's coordinate system
  // to drawing system coordinates and vice versa. Note that the coordinates
  // have
  // the origin in the top left corner, which is how Qt and Cairo have it, no
  // doubt a holdover from X Windows. This means that a higher y value will be
  // nearer the bottom of the screen. This doesn't really matter except when
  // doing text superscripts and subscripts.

  //! transform a point from the molecule coordinate system into the drawing
  //! coordinate system.
  //! Prefers globalDrawTrans_ if it exists, otherwise
  //! uses drawMols_[activeMolIdx_].
  virtual Point2D getDrawCoords(const Point2D &mol_cds) const;
  //! returns the drawing coordinates of a particular atom
  virtual Point2D getDrawCoords(int at_num) const;
  //! Prefers globalDrawTrans_ if it exists, otherwise
  //! uses drawMols_[activeMolIdx_].
  virtual Point2D getAtomCoords(const std::pair<int, int> &screen_cds) const;
  //! transform a point from drawing coordinates to the molecule coordinate
  //! system. Prefers globalDrawTrans_ if it exists, otherwise
  //! uses drawMols_[activeMolIdx_].
  virtual Point2D getAtomCoords(
      const std::pair<double, double> &screen_cds) const;
  //! returns the molecular coordinates of a particular atom.  at_num refers
  //! to the atom in activeMolIdx_.
  virtual Point2D getAtomCoords(int at_num) const;
  //! @}
  //! returns the coordinates of the atoms of the activeMolIdx_ molecule in
  //! molecular coordinates.
  const std::vector<Point2D> &atomCoords() const;
  //! returns the atomic symbols of the activeMolIdx_ molecule
  const std::vector<std::pair<std::string, MolDraw2D_detail::OrientType>>
      &atomSyms() const;

  //! return the width of the drawing area.
  virtual int width() const { return width_; }
  //! return the height of the drawing area.
  virtual int height() const { return height_; }
  //! return the width of the drawing panels.
  virtual int panelWidth() const { return panel_width_; }
  //! return the height of the drawing panels.
  virtual int panelHeight() const { return panel_height_; }
  virtual int drawHeight() const { return panel_height_ - legend_height_; }
  // returns the width to draw a line in draw coords.
  virtual double getDrawLineWidth() const;

  //! returns the drawing scale (conversion from molecular coords -> drawing
  /// coords)
  double scale() const;
  //! explicitly sets the scaling factors for the drawing
  void setScale(double newScale);
  void setScale(int width, int height, const Point2D &minv, const Point2D &maxv,
                const ROMol *mol = nullptr);
  //! sets the drawing offset (in drawing coords)
  void setOffset(int x, int y) {
    x_offset_ = x;
    y_offset_ = y;
  }
  //! returns the drawing offset (in drawing coords)
  Point2D offset() const { return Point2D(x_offset_, y_offset_); }

  //! returns the minimum point of the drawing (in molecular coords)
  Point2D minPt() const;
  //! returns the width and height of the grid (in molecular coords)
  Point2D range() const;

  //! font size in drawing coordinate units. That's probably pixels.
  virtual double fontSize() const;
  virtual void setFontSize(double new_size);

  //! sets the current draw color
  virtual void setColour(const DrawColour &col) { curr_colour_ = col; }
  //! returns the current draw color
  virtual DrawColour colour() const { return curr_colour_; }
  //! sets the current dash pattern
  virtual void setDash(const DashPattern &patt) { curr_dash_ = patt; }
  //! returns the current dash pattern
  virtual const DashPattern &dash() const { return curr_dash_; }

  //! sets the current line width
  virtual void setLineWidth(double width) { drawOptions().bondLineWidth = width; }
  //! returns the current line width
  virtual double lineWidth() const { return drawOptions().bondLineWidth; }

  //! using the current scale, work out the size of the label in molecule
  //! coordinates.
  /*!
     Bear in mind when implementing this, that, for example, NH2 will appear as
     NH<sub>2</sub> to convey that the 2 is a subscript, and this needs to
     accounted for in the width and height.
   */
  virtual void getStringSize(const std::string &label, double &label_width,
                             double &label_height) const;
  // get the overall size of the label, allowing for it being split
  // into pieces according to orientation.
  void getLabelSize(const std::string &label,
                    MolDraw2D_detail::OrientType orient, double &label_width,
                    double &label_height) const;
  // return extremes for string in molecule coords.
  void getStringExtremes(const std::string &label,
                         MolDraw2D_detail::OrientType orient,
                         const Point2D &cds, double &x_min, double &y_min,
                         double &x_max, double &y_max) const;

  //! adds additional information about the atoms to the output. Does not make
  //! sense for all renderers.
  virtual void tagAtoms(const ROMol &mol) { RDUNUSED_PARAM(mol); }
  //! set whether or not polygons are being filled
  virtual bool fillPolys() const { return fill_polys_; }
  //! returns either or not polygons should be filled
  virtual void setFillPolys(bool val) { fill_polys_ = val; }

  //! returns our current drawing options
  MolDrawOptions &drawOptions() { return options_; }
  //! \overload
  const MolDrawOptions &drawOptions() const { return options_; }

  virtual bool supportsAnnotations() const { return true; }

  void setActiveMolIdx(int newIdx);
  bool hasActiveAtmIdx() const { return activeAtmIdx1_ >= 0; }
  int getActiveAtmIdx1() const { return activeAtmIdx1_; }
  int getActiveAtmIdx2() const { return activeAtmIdx2_; }
  void setActiveAtmIdx(int at_idx1 = -1, int at_idx2 = -1);
  bool hasActiveBndIdx() const { return activeBndIdx_ >= 0; }
  int getActiveBndIdx() const { return activeBndIdx_; }
  void setActiveBndIdx(int bnd_idx = -1) {
    activeBndIdx_ = (bnd_idx < 0 ? -1 : bnd_idx);
  }
  void setActiveClass(std::string actClass = std::string("")) {
    d_activeClass = actClass;
  }
  std::string getActiveClass() const { return d_activeClass; }

 protected:
  std::unique_ptr<MolDraw2D_detail::DrawText> text_drawer_;
  std::string d_activeClass;
  bool needs_init_ = true;
  std::vector<std::pair<std::string, std::string>> d_metadata;
  unsigned int d_numMetadataEntries = 0;

 private:
  //! \name Methods that must be provided by child classes
  //! @{
  virtual void initDrawing() = 0;
  virtual void initTextDrawer(bool noFreetype) = 0;

  // if the width or height of the DrawMol was -1, the new dimensions need to be
  // transferred to MolDraw2D.
  void fixVariableDimensions(const MolDraw2D_detail::DrawMol &drawMol);

  // split the reaction up into the reagents, products and agents, each as
  // a separate entity with its own scale.
  void getReactionDrawMols(
      const ChemicalReaction &rxn, bool highlightByReactant,
      const std::vector<DrawColour> *highlightColorsReactants,
      const std::vector<int> *confIds,
      std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &reagents,
      std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &products,
      std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &agents,
      int &plusWidth);
  // take the given components from the reaction (bits will be either
  // reagents, products or agents) and create the corresponding DrawMols.
  void makeReactionComponents(
      std::vector<RDKit::ROMOL_SPTR> const &bits,
      const std::vector<int> *confIds, int heightToUse,
      std::map<int, DrawColour> &atomColours,
      std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &dms,
      double &minScale, double &minFontScale);
  // this puts a pointer to the DrawMol into _drawMols as well, hence the use
  // of shared_ptr for reagents, products and agents above.
  void makeReactionDrawMol(
      RWMol &mol, int confId, int molHeight,
      const std::vector<int> &highlightAtoms,
      const std::vector<int> &highlightBonds,
      const std::map<int, DrawColour> &highlightAtomMap,
      const std::map<int, DrawColour> &highlightBondMap,
      std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &mols);
  // Strictly speaking, this isn't actually a const function, although the
  // compiler can't spot it, because the scales of reagents etc may be changed,
  // and they are also in drawMols_.
  void calcReactionOffsets(
      std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &reagents,
      std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &products,
      std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &agents,
      int &plusWidth, std::vector<Point2D> &offsets, Point2D &arrowBeg,
      Point2D &arrowEnd);
  // returns the final offset. plusWidth of 0 means no pluses to be drawn.
  int drawReactionPart(
      std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> &reactBit,
      int plusWidth, int initOffset, const std::vector<Point2D> &offsets);
  // returns a map of colours indexed by the atomMapNum. Each reagent gives
  // a different colour.
  void findReactionHighlights(
      const ChemicalReaction &rxn, bool highlightByReactant,
      const std::vector<DrawColour> *highlightColorsReactants,
      std::map<int, DrawColour> &atomColours) const;

  int width_, height_, panel_width_, panel_height_, legend_height_;
  // if the user calls setScale() to explicitly force a scale on the
  // DrawMols, this is set to true.
  bool forceScale_ = false;
  double scale_, fontScale_;
  int x_offset_, y_offset_;  // translation in screen coordinates
  bool fill_polys_;
  int activeMolIdx_;
  int activeAtmIdx1_;
  int activeAtmIdx2_;
  int activeBndIdx_;
  // these are shared_ptr rather than unique_ptr because the reactions
  // keep their own copy.
  std::vector<std::shared_ptr<MolDraw2D_detail::DrawMol>> drawMols_;
  // this is for when we want to set up a MolDraw2D to a given scale and be
  // able to draw molecules and arbitrary lines, arcs etc. onto it all to the
  // same drawing transformation.  If present, it will always be applied to
  // any new drawMols_ before they are drawn.  A separate class might have
  // been better, but this is convenient.
  std::unique_ptr<MolDraw2D_detail::DrawMol> globalDrawTrans_;

  DrawColour curr_colour_;
  DashPattern curr_dash_;
  MolDrawOptions options_;

  // Do the drawing, the new way
  void startDrawing();
  void drawTheMolecule(MolDraw2D_detail::DrawMol &drawMol);
  void setupTextDrawer();

  virtual void updateMetadata(const ROMol &mol, int confId) {
    RDUNUSED_PARAM(mol);
    RDUNUSED_PARAM(confId);
  }
  virtual void updateMetadata(const ChemicalReaction &rxn) {
    RDUNUSED_PARAM(rxn);
  }
};

inline void setDarkMode(MolDrawOptions &opts) {
  assignDarkModePalette(opts.atomColourPalette);
  opts.backgroundColour = DrawColour{0, 0, 0, 1};
  opts.annotationColour = DrawColour{0.9, 0.9, 0.9, 1};
  opts.legendColour = DrawColour{0.9, 0.9, 0.9, 1};
  opts.symbolColour = DrawColour{0.9, 0.9, 0.9, 1};
  opts.variableAttachmentColour = DrawColour{0.3, 0.3, 0.3, 1};
}
inline void setDarkMode(MolDraw2D &d2d) { setDarkMode(d2d.drawOptions()); }
inline void setMonochromeMode(MolDrawOptions &opts, const DrawColour &fgColour,
                              const DrawColour &bgColour) {
  opts.atomColourPalette.clear();
  opts.atomColourPalette[-1] = fgColour;
  opts.backgroundColour = bgColour;
  opts.annotationColour = fgColour;
  opts.legendColour = fgColour;
  opts.symbolColour = fgColour;
  opts.variableAttachmentColour = fgColour;
}
inline void setMonochromeMode(MolDraw2D &drawer, const DrawColour &fgColour,
                              const DrawColour &bgColour) {
  setMonochromeMode(drawer.drawOptions(), fgColour, bgColour);
}

}  // namespace RDKit

#endif  // RDKITMOLDRAW2D_H
