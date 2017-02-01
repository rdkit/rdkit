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

#ifndef RDKITMOLDRAW2D_H
#define RDKITMOLDRAW2D_H

#include <vector>

#include <Geometry/point.h>
#include <GraphMol/RDKitBase.h>

#include <boost/tuple/tuple.hpp>

// ****************************************************************************
using RDGeom::Point2D;

namespace RDKit {

typedef boost::tuple<float, float, float> DrawColour;
typedef std::vector<unsigned int> DashPattern;

struct MolDrawOptions {
  bool atomLabelDeuteriumTritium;  // toggles replacing 2H with D and 3H with T
  bool dummiesAreAttachments;      // draws "breaks" at dummy atoms
  bool circleAtoms;                // draws circles under highlighted atoms
  DrawColour highlightColour;      // default highlight color
  bool continuousHighlight;  // highlight by drawing an outline *underneath* the
                             // molecule
  int flagCloseContactsDist;  // if positive, this will be used as a cutoff (in
                              // pixels) for highlighting close contacts
  bool includeAtomTags;  // toggles inclusion of atom tags in the output. does
                         // not make sense for all renderers.
  bool clearBackground;  // toggles clearing the background before drawing a
                         // molecule
  DrawColour
      backgroundColour;  // color to be used while clearing the background
  int legendFontSize;    // font size (in pixels) to be used for the legend (if
                         // present)
  DrawColour legendColour;    // color to be used for the legend (if present)
  double multipleBondOffset;  // offset (in Angstroms) for the extra lines in a
                              // multiple bond
  double padding;  // fraction of empty space to leave around the molecule
  double additionalAtomLabelPadding;  // additional padding to leave around atom
                                      // labels. Expressed as a fraction of the
                                      // font size.
  std::map<int, std::string> atomLabels;       // replacement labels for atoms
  std::vector<std::vector<int> > atomRegions;  // regions

  MolDrawOptions()
      : atomLabelDeuteriumTritium(false),
        dummiesAreAttachments(false),
        circleAtoms(true),
        highlightColour(1, .5, .5),
        continuousHighlight(true),
        flagCloseContactsDist(3),
        includeAtomTags(false),
        clearBackground(true),
        backgroundColour(1, 1, 1),
        legendFontSize(12),
        legendColour(0, 0, 0),
        multipleBondOffset(0.15),
        padding(0.05),
        additionalAtomLabelPadding(0.0){};
};

class MolDraw2D {
 public:
  typedef enum { C = 0, N, E, S, W } OrientType;
  typedef enum {
    TextDrawNormal = 0,
    TextDrawSuperscript,
    TextDrawSubscript
  } TextDrawType;

  MolDraw2D(int width, int height, int panelWidth = -1, int panelHeight = -1);

  virtual ~MolDraw2D() {}

  virtual void drawMolecule(
      const ROMol &mol, const std::vector<int> *highlight_atoms = NULL,
      const std::map<int, DrawColour> *highlight_map = NULL,
      const std::map<int, double> *highlight_radii = NULL, int confId = -1);

  virtual void drawMolecule(
      const ROMol &mol, const std::string &legend,
      const std::vector<int> *highlight_atoms = NULL,
      const std::map<int, DrawColour> *highlight_map = NULL,
      const std::map<int, double> *highlight_radii = NULL, int confId = -1);

  virtual void drawMolecule(
      const ROMol &mol, const std::vector<int> *highlight_atoms,
      const std::vector<int> *highlight_bonds,
      const std::map<int, DrawColour> *highlight_atom_map = NULL,
      const std::map<int, DrawColour> *highlight_bond_map = NULL,
      const std::map<int, double> *highlight_radii = NULL, int confId = -1);

  virtual void drawMolecule(
      const ROMol &mol, const std::string &legend,
      const std::vector<int> *highlight_atoms,
      const std::vector<int> *highlight_bonds,
      const std::map<int, DrawColour> *highlight_atom_map = NULL,
      const std::map<int, DrawColour> *highlight_bond_map = NULL,
      const std::map<int, double> *highlight_radii = NULL, int confId = -1);

  virtual void drawMolecules(
      const std::vector<ROMol *> &mols,
      const std::vector<std::string> *legends = NULL,
      const std::vector<std::vector<int> > *highlight_atoms = NULL,
      const std::vector<std::vector<int> > *highlight_bonds = NULL,
      const std::vector<std::map<int, DrawColour> > *highlight_atom_maps = NULL,
      const std::vector<std::map<int, DrawColour> > *highlight_bond_maps = NULL,
      const std::vector<std::map<int, double> > *highlight_radii = NULL,
      const std::vector<int> *confIds = NULL);

  // transform a set of coords in the molecule's coordinate system
  // to drawing system coordinates and vice versa. Note that the coordinates
  // have
  // the origin in the top left corner, which is how Qt and Cairo have it, no
  // doubt a holdover from X Windows. This means that a higher y value will be
  // nearer the bottom of the screen. This doesn't really matter except when
  // doing text superscripts and subscripts.
  virtual Point2D getDrawCoords(const Point2D &mol_cds) const;
  virtual Point2D getDrawCoords(int at_num) const;
  virtual Point2D getAtomCoords(const std::pair<int, int> &screen_cds) const;
  virtual Point2D getAtomCoords(
      const std::pair<double, double> &screen_cds) const;
  virtual Point2D getAtomCoords(int at_num) const;

  virtual int width() const { return width_; }
  virtual int height() const { return height_; }
  virtual int panelWidth() const { return panel_width_; }
  virtual int panelHeight() const { return panel_height_; }

  double scale() const { return scale_; }
  void calculateScale(int width, int height);
  void calculateScale() { calculateScale(panel_width_, panel_height_); };
  void setScale(int width, int height, const Point2D &minv,
                const Point2D &maxv);
  void setOffset(int x, int y) {
    x_offset_ = x;
    y_offset_ = y;
  }
  Point2D offset() { return Point2D(x_offset_, y_offset_); }

  virtual double fontSize() const { return font_size_; }
  // set font size in molecule coordinate units. That's probably Angstrom for
  // RDKit.
  virtual void setFontSize(double new_size);

  virtual void setColour(const DrawColour &col) { curr_colour_ = col; }
  virtual DrawColour colour() const { return curr_colour_; }

  virtual void setDash(const DashPattern &patt) { curr_dash_ = patt; }
  virtual const DashPattern &dash() const { return curr_dash_; }

  virtual void setLineWidth(int width) { curr_width_ = width; }
  virtual int lineWidth() const { return curr_width_; }

  // establishes whether to put string draw mode into super- or sub-script
  // mode based on contents of instring from i onwards. Increments i
  // appropriately
  // and returns true or false depending on whether it did something or not
  bool setStringDrawMode(const std::string &instring, TextDrawType &draw_mode,
                         int &i) const;

  virtual void clearDrawing() = 0;

  virtual void drawLine(const Point2D &cds1, const Point2D &cds2) = 0;

  // using the current scale, work out the size of the label in molecule
  // coordinates.
  // Bear in mind when implementing this, that, for example, NH2 will appear as
  // NH<sub>2</sub> to convey that the 2 is a subscript, and this needs to
  // accounted
  // for in the width and height.
  virtual void getStringSize(const std::string &label, double &label_width,
                             double &label_height) const = 0;
  // drawString centres the string on cds.
  virtual void drawString(const std::string &str, const Point2D &cds);

  // draw polygons:
  virtual void drawPolygon(const std::vector<Point2D> &cds) = 0;
  virtual void drawTriangle(const Point2D &cds1, const Point2D &cds2,
                            const Point2D &cds3);
  virtual void drawEllipse(const Point2D &cds1, const Point2D &cds2);
  virtual void drawRect(const Point2D &cds1, const Point2D &cds2);
  virtual void drawAttachmentLine(const Point2D &cds1, const Point2D &cds2,
                                  const DrawColour &col, double len = 1.0,
                                  unsigned int nSegments = 16);
  virtual void drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                            const DrawColour &col1, const DrawColour &col2,
                            unsigned int nSegments = 16,
                            double vertOffset = 0.05);
  virtual void tagAtoms(const ROMol &mol) { RDUNUSED_PARAM(mol); };

  virtual bool fillPolys() const { return fill_polys_; }
  virtual void setFillPolys(bool val) { fill_polys_ = val; }

  MolDrawOptions &drawOptions() { return options_; }
  const MolDrawOptions &drawOptions() const { return options_; }

  const std::vector<Point2D> &atomCoords() const {
    PRECONDITION(activeMolIdx_ >= 0, "no index");
    return at_cds_[activeMolIdx_];
  };
  const std::vector<std::pair<std::string, OrientType> > &atomSyms() const {
    PRECONDITION(activeMolIdx_ >= 0, "no index");
    return atom_syms_[activeMolIdx_];
  };

 private:
  bool needs_scale_;
  int width_, height_, panel_width_, panel_height_;
  double scale_;
  double x_min_, y_min_, x_range_, y_range_;
  double x_trans_, y_trans_;
  int x_offset_, y_offset_;  // translation in screen coordinates
  // font_size_ in molecule coordinate units. Default 0.5 (a bit bigger
  // than the default width of a double bond)
  double font_size_;
  int curr_width_;
  bool fill_polys_;
  int activeMolIdx_;

  DrawColour curr_colour_;
  DashPattern curr_dash_;
  MolDrawOptions options_;

  std::vector<std::vector<Point2D> > at_cds_;  // from mol
  std::vector<std::vector<int> > atomic_nums_;
  std::vector<std::vector<std::pair<std::string, OrientType> > > atom_syms_;
  Point2D bbox_[2];

  // draw the char, with the bottom left hand corner at cds
  virtual void drawChar(char c, const Point2D &cds) = 0;

  // return a DrawColour based on the contents of highlight_atoms or
  // highlight_map, falling back to atomic number by default
  DrawColour getColour(int atom_idx,
                       const std::vector<int> *highlight_atoms = NULL,
                       const std::map<int, DrawColour> *highlight_map = NULL);
  DrawColour getColourByAtomicNum(int atomic_num);

  void extractAtomCoords(const ROMol &mol, int confId, bool updateBBox);
  void extractAtomSymbols(const ROMol &mol);

  virtual void drawLine(const Point2D &cds1, const Point2D &cds2,
                        const DrawColour &col1, const DrawColour &col2);
  void drawBond(const ROMol &mol, const BOND_SPTR &bond, int at1_idx,
                int at2_idx, const std::vector<int> *highlight_atoms = NULL,
                const std::map<int, DrawColour> *highlight_atom_map = NULL,
                const std::vector<int> *highlight_bonds = NULL,
                const std::map<int, DrawColour> *highlight_bond_map = NULL);
  void drawWedgedBond(const Point2D &cds1, const Point2D &cds2,
                      bool draw_dashed, const DrawColour &col1,
                      const DrawColour &col2);
  void drawAtomLabel(int atom_num,
                     const std::vector<int> *highlight_atoms = NULL,
                     const std::map<int, DrawColour> *highlight_map = NULL);
  // cds1 and cds2 are 2 atoms in a ring.  Returns the perpendicular pointing
  // into
  // the ring.
  Point2D bondInsideRing(const ROMol &mol, const BOND_SPTR &bond,
                         const Point2D &cds1, const Point2D &cds2);
  // cds1 and cds2 are 2 atoms in a chain double bond.  Returns the
  // perpendicular
  // pointing into the inside of the bond
  Point2D bondInsideDoubleBond(const ROMol &mol, const BOND_SPTR &bond);
  // calculate normalised perpendicular to vector between two coords, such
  // that
  // it's inside the angle made between (1 and 2) and (2 and 3).
  Point2D calcInnerPerpendicular(const Point2D &cds1, const Point2D &cds2,
                                 const Point2D &cds3);

  // take the coords for atnum, with neighbour nbr_cds, and move cds out to
  // accommodate
  // the label associated with it.
  void adjustBondEndForLabel(int atnum, const Point2D &nbr_cds,
                             Point2D &cds) const;

  // adds LaTeX-like annotation for super- and sub-script.
  std::pair<std::string, OrientType> getAtomSymbolAndOrientation(
      const Atom &atom, const Point2D &nbr_sum);

 protected:
  virtual void doContinuousHighlighting(
      const ROMol &mol, const std::vector<int> *highlight_atoms,
      const std::vector<int> *highlight_bonds,
      const std::map<int, DrawColour> *highlight_atom_map,
      const std::map<int, DrawColour> *highlight_bond_map,
      const std::map<int, double> *highlight_radii);

  virtual void highlightCloseContacts();

  // calculate normalised perpendicular to vector between two coords
  Point2D calcPerpendicular(const Point2D &cds1, const Point2D &cds2);
};
}

#endif  // RDKITMOLDRAW2D_H
