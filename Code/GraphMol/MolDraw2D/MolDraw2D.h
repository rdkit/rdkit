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
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>

// ****************************************************************************
using RDGeom::Point2D;

namespace RDKit {

class DrawText;
enum class TextAlignType : unsigned char;
enum class OrientType : unsigned char;

struct DrawColour {
  double r = 0.0, g = 0.0, b = 0.0, a = 1.0;
  DrawColour() = default;
  DrawColour(double r, double g, double b, double a = 1.0)
      : r(r), g(g), b(b), a(a){};
  bool operator==(const DrawColour &other) const {
    return r == other.r && g == other.g && b == other.b && a == other.a;
  }
  bool feq(const DrawColour &other, double tol = 0.001,
           bool ignoreAlpha = true) const {
    return fabs(r - other.r) <= tol && fabs(g - other.g) <= tol &&
           fabs(b - other.b) <= tol &&
           (ignoreAlpha || fabs(a - other.a) <= tol);
  };
  DrawColour operator+(const DrawColour &other) const {
    return {r + other.r, g + other.g, b + other.b, a + other.a};
  }
  DrawColour operator-(const DrawColour &other) const {
    return {r - other.r, g - other.g, b - other.b, a - other.a};
  }
  DrawColour operator/(double v) const {
    PRECONDITION(v != 0.0, "divide by zero");
    return {r / v, g / v, b / v, a / v};
  }
  DrawColour operator*(double v) const { return {r * v, g * v, b * v, a * v}; }
};

// for holding dimensions of the rectangle round a string.
struct StringRect {
  Point2D trans_;     // Where to draw char relative to other chars in string
  Point2D offset_;    // offset for draw coords so char is centred correctly
  Point2D g_centre_;  // glyph centre relative to the origin of the char.
  double y_shift_;  // shift the whole thing in y by this. For multi-line text.
  double width_, height_;  // of the glyph itself, not the character cell
  double rect_corr_;  // because if we move a char one way, we need to move the
                      // rectangle the other.
  int clash_score_;   // rough measure of how badly it clashed with other things
                      // lower is better, 0 is no clash.
  StringRect()
      : trans_(0.0, 0.0),
        offset_(0.0, 0.0),
        g_centre_(offset_),
        y_shift_(0.0),
        width_(0.0),
        height_(0.0),
        rect_corr_(0.0),
        clash_score_(0) {}
  StringRect(const Point2D &offset, const Point2D &g_centre, double w, double h)
      : trans_(0.0, 0.0),
        offset_(offset),
        g_centre_(g_centre),
        y_shift_(0.0),
        width_(w),
        height_(h),
        rect_corr_(0.0),
        clash_score_(0) {}
  // tl is top, left; br is bottom, right of the glyph, relative to the
  // centre. Padding in draw coords.
  void calcCorners(Point2D &tl, Point2D &tr, Point2D &br, Point2D &bl,
                   double padding) const {
    double wb2 = padding + width_ / 2.0;
    double hb2 = padding + height_ / 2.0;
    Point2D c = trans_ + g_centre_ - offset_;
    c.y -= y_shift_;
    tl = Point2D(c.x - wb2, c.y - hb2);
    tr = Point2D(c.x + wb2, c.y - hb2);
    br = Point2D(c.x + wb2, c.y + hb2);
    bl = Point2D(c.x - wb2, c.y + hb2);
  }
  bool doesItIntersect(const StringRect &other) const {
    Point2D ttl, ttr, tbr, tbl;
    calcCorners(ttl, ttr, tbr, tbl, 0.0);
    // is +ve y up or down?
    if (ttl.y < tbl.y) {
      std::swap(ttl, tbl);
      std::swap(ttr, tbr);
    }
    Point2D otl, otr, obr, obl;
    other.calcCorners(otl, otr, obr, obl, 0.0);
    if (otl.y < obl.y) {
      std::swap(otl, obl);
      std::swap(otr, obr);
    }
    if ((otl.x >= ttl.x && otl.x <= ttr.x && otl.y >= tbl.y &&
         otl.y <= ttl.y) ||
        (otr.x >= ttl.x && otr.x <= ttr.x && otr.y >= tbl.y &&
         otr.y <= ttl.y) ||
        (obr.x >= ttl.x && obr.x <= ttr.x && obr.y >= tbl.y &&
         obr.y <= ttl.y) ||
        (obl.x >= ttl.x && obl.x <= ttr.x && obl.y >= tbl.y &&
         obl.y <= ttl.y)) {
      return true;
    }
    if ((ttl.x >= otl.x && ttl.x <= otr.x && ttl.y >= obl.y &&
         ttl.y <= otl.y) ||
        (ttr.x >= otl.x && ttr.x <= otr.x && ttr.y >= obl.y &&
         ttr.y <= otl.y) ||
        (tbr.x >= otl.x && tbr.x <= otr.x && tbr.y >= obl.y &&
         tbr.y <= otl.y) ||
        (tbl.x >= otl.x && tbl.x <= otr.x && tbl.y >= obl.y &&
         tbl.y <= otl.y)) {
      return true;
    }
    return false;
  }
};

typedef std::map<int, DrawColour> ColourPalette;
typedef std::vector<unsigned int> DashPattern;

inline void assignDefaultPalette(ColourPalette &palette) {
  palette.clear();
  palette[-1] = DrawColour(0, 0, 0);
  palette[0] = DrawColour(0.1, 0.1, 0.1);
  palette[1] = palette[6] = DrawColour(0.0, 0.0, 0.0);
  palette[7] = DrawColour(0.0, 0.0, 1.0);
  palette[8] = DrawColour(1.0, 0.0, 0.0);
  palette[9] = DrawColour(0.2, 0.8, 0.8);
  palette[15] = DrawColour(1.0, 0.5, 0.0);
  palette[16] = DrawColour(0.8, 0.8, 0.0);
  palette[17] = DrawColour(0.0, 0.802, 0.0);
  palette[35] = DrawColour(0.5, 0.3, 0.1);
  palette[53] = DrawColour(0.63, 0.12, 0.94);
};

inline void assignBWPalette(ColourPalette &palette) {
  palette.clear();
  palette[-1] = DrawColour(0, 0, 0);
};

struct RDKIT_MOLDRAW2D_EXPORT MolDrawOptions {
  bool atomLabelDeuteriumTritium =
      false;  // toggles replacing 2H with D and 3H with T
  bool dummiesAreAttachments = false;  // draws "breaks" at dummy atoms
  bool circleAtoms = true;             // draws circles under highlighted atoms
  DrawColour highlightColour{1, 0.5, 0.5};  // default highlight color
  bool continuousHighlight = true;          // highlight by drawing an outline
                                            // *underneath* the molecule
  bool fillHighlights = true;     // fill the areas used to highlight atoms and
                                  // atom regions
  double highlightRadius = 0.3;   // default if nothing given for a particular
                                  // atom. units are "Angstrom"
  int flagCloseContactsDist = 3;  // if positive, this will be used as a cutoff
                                  // (in pixels) for highlighting close contacts
  bool includeAtomTags =
      false;  // toggles inclusion of atom tags in the output. does
              // not make sense for all renderers.
  bool clearBackground = true;  // toggles clearing the background before
                                // drawing a molecule
  DrawColour backgroundColour{
      1, 1, 1};             // color to be used while clearing the background
  int legendFontSize = 16;  // font size (in pixels) to be used for the legend
                            // (if present)
  int maxFontSize = 40;  // maximum size in pixels for font in drawn molecule.
                         // -1 means no max.
  int minFontSize = 6;   // likewise for -1.
  double annotationFontScale = 0.5;  // scales font relative to atom labels for
                                     // atom and bond annotation.
  std::string fontFile = "";  // name of font for freetype rendering.  If given,
                              // over-rides default
  DrawColour legendColour{0, 0,
                          0};  // color to be used for the legend (if present)
  double multipleBondOffset = 0.15;  // offset (in Angstrom) for the extra lines
                                     // in a multiple bond
  double padding =
      0.05;  // fraction of empty space to leave around the molecule
  double additionalAtomLabelPadding = 0.0;  // additional padding to leave
                                            // around atom labels. Expressed as
                                            // a fraction of the font size.
  std::map<int, std::string> atomLabels;    // replacement labels for atoms
  std::vector<std::vector<int>> atomRegions;  // regions
  DrawColour symbolColour{
      0, 0, 0};  // color to be used for the symbols and arrows in reactions
  int bondLineWidth = 2;        // default line width when drawing bonds
  bool scaleBondWidth = false;  // whether to apply scale() to the bond width
  bool scaleHighlightBondWidth = true;   // likewise with bond highlights.
  int highlightBondWidthMultiplier = 8;  // what to multiply standard bond width
                                         // by for highlighting.
  bool prepareMolsBeforeDrawing = true;  // call prepareMolForDrawing() on each
                                         // molecule passed to drawMolecules()
  std::vector<DrawColour> highlightColourPalette;  // defining 10 default colors
  // for highlighting atoms and bonds
  // or reactants in a reactions
  ColourPalette atomColourPalette;  // the palette used to assign
                                    // colors to atoms based on
                                    // atomic number.
  double fixedScale =
      -1.0;  // fixes scale to this fraction of draw window width, so
             // an average bond is this fraction of the width.  If
             // scale comes out smaller than this, reduces scale, but
             // won't make it larger.  The default of -1.0 means no fix.
  double fixedBondLength =
      -1.0;             // fixes the bond length (and hence the scale) to
                        // always be this number of pixels.  Assuming a bond
                        // length in coordinates is 1, as is normal.  If
                        // scale comes out smaller than this, reduces scale,
                        // but won't make it larger.  The default -1.0 means no
                        // fix. If both fixedScale and fixedBondLength are >
                        // 0.0, fixedScale wins.
  double rotate = 0.0;  // angle in degrees to rotate coords by about centre
                        // before drawing.
  bool addAtomIndices = false;  // adds atom indices to drawings.
  bool addBondIndices = false;  // adds bond indices to drawings.

  bool addStereoAnnotation = false;       // adds E/Z and R/S to drawings.
  bool atomHighlightsAreCircles = false;  // forces atom highlights always to be
                                          // circles. Default (false) is to put
                                          // ellipses round longer labels.
  bool centreMoleculesBeforeDrawing = false;  // moves the centre of the drawn
                                              // molecule to (0,0)
  bool explicitMethyl = false;  // draw terminal methyl and related as CH3
  bool includeRadicals =
      true;  // include radicals in the drawing (it can be useful to turn this
             // off for reactions and queries)
  bool includeMetadata =
      true;  // when possible include metadata about molecules and reactions in
             // the output to allow them to be reconstructed

  MolDrawOptions() {
    highlightColourPalette.emplace_back(
        DrawColour(1., 1., .67));  // popcorn yellow
    highlightColourPalette.emplace_back(DrawColour(1., .8, .6));  // sand
    highlightColourPalette.emplace_back(
        DrawColour(1., .71, .76));  // light pink
    highlightColourPalette.emplace_back(
        DrawColour(.8, 1., .8));  // offwhitegreen
    highlightColourPalette.emplace_back(DrawColour(.87, .63, .87));  // plum
    highlightColourPalette.emplace_back(
        DrawColour(.76, .94, .96));  // pastel blue
    highlightColourPalette.emplace_back(
        DrawColour(.67, .67, 1.));  // periwinkle
    highlightColourPalette.emplace_back(DrawColour(.64, .76, .34));  // avocado
    highlightColourPalette.emplace_back(
        DrawColour(.56, .93, .56));  // light green
    highlightColourPalette.emplace_back(DrawColour(.20, .63, .79));  // peacock
    assignDefaultPalette(atomColourPalette);
  };
};

//! MolDraw2D is the base class for doing 2D renderings of molecules
class RDKIT_MOLDRAW2D_EXPORT MolDraw2D {
 public:
  //! constructor for a particular size
  /*!
    \param width       : width (in pixels) of the rendering
    \param height      : height (in pixels) of the rendering
    \param panelWidth  : (optional) width (in pixels) of a single panel
    \param panelHeight : (optional) height (in pixels) of a single panel

    The \c panelWidth and \c panelHeight arguments are used to provide the
    sizes of the panels individual molecules are drawn in when
    \c drawMolecules() is called.
  */
  MolDraw2D(int width, int height, int panelWidth, int panelHeight);
  virtual ~MolDraw2D();

  //! \name Methods that must be provided by child classes
  //@{
 private:
  virtual void initDrawing() = 0;
  virtual void initTextDrawer(bool noFreetype) = 0;

 public:
  //! clears the contents of the drawing
  virtual void clearDrawing() = 0;
  //! draws a line from \c cds1 to \c cds2 using the current drawing style
  // in atom coords.
  virtual void drawLine(const Point2D &cds1, const Point2D &cds2) = 0;
  //! draw a polygon.  Note that if fillPolys() returns false, it
  //! doesn't close the path.  If you want it to in that case, you
  //! do it explicitly yourself.
  virtual void drawPolygon(const std::vector<Point2D> &cds) = 0;
  //@}

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

  //! \name Transformations
  //@{
  // transform a set of coords in the molecule's coordinate system
  // to drawing system coordinates and vice versa. Note that the coordinates
  // have
  // the origin in the top left corner, which is how Qt and Cairo have it, no
  // doubt a holdover from X Windows. This means that a higher y value will be
  // nearer the bottom of the screen. This doesn't really matter except when
  // doing text superscripts and subscripts.

  //! transform a point from the molecule coordinate system into the drawing
  //! coordinate system
  virtual Point2D getDrawCoords(const Point2D &mol_cds) const;
  //! returns the drawing coordinates of a particular atom
  virtual Point2D getDrawCoords(int at_num) const;
  virtual Point2D getAtomCoords(const std::pair<int, int> &screen_cds) const;
  //! transform a point from drawing coordinates to the molecule coordinate
  //! system
  virtual Point2D getAtomCoords(
      const std::pair<double, double> &screen_cds) const;
  //! returns the molecular coordinates of a particular atom
  virtual Point2D getAtomCoords(int at_num) const;
  //@}
  //! return the width of the drawing area.
  virtual int width() const { return width_; }
  //! return the height of the drawing area.
  virtual int height() const { return height_; }
  //! return the width of the drawing panels.
  virtual int panelWidth() const { return panel_width_; }
  //! return the height of the drawing panels.
  virtual int panelHeight() const { return panel_height_; }
  virtual int drawHeight() const { return panel_height_ - legend_height_; }

  //! returns the drawing scale (conversion from molecular coords -> drawing
  // coords)
  double scale() const { return scale_; }
  //! calculates the drawing scale (conversion from molecular coords -> drawing
  // coords)
  void calculateScale(int width, int height, const ROMol &mol,
                      const std::vector<int> *highlight_atoms = nullptr,
                      const std::map<int, double> *highlight_radii = nullptr);
  //! overload
  // calculate a single scale that will suit all molecules.  For use by
  // drawMolecules primarily.
  void calculateScale(int width, int height, const std::vector<ROMol *> &mols,
                      const std::vector<std::vector<int>> *highlight_atoms,
                      const std::vector<std::map<int, double>> *highlight_radii,
                      const std::vector<int> *confIds,
                      std::vector<std::unique_ptr<RWMol>> &tmols);
  // set [xy]_trans_ to the middle of the draw area in molecule coords
  void centrePicture(int width, int height);

  //! explicitly sets the scaling factors for the drawing
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
  Point2D minPt() const { return Point2D(x_min_, y_min_); }
  //! returns the width and height of the grid (in molecular coords)
  Point2D range() const { return Point2D(x_range_, y_range_); }

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
  virtual void setLineWidth(int width) { drawOptions().bondLineWidth = width; }
  //! returns the current line width
  virtual int lineWidth() const { return drawOptions().bondLineWidth; }

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
  void getLabelSize(const std::string &label, OrientType orient,
                    double &label_width, double &label_height) const;
  // return extremes for string in molecule coords.
  void getStringExtremes(const std::string &label, OrientType orient,
                         const Point2D &cds, double &x_min, double &y_min,
                         double &x_max, double &y_max) const;

  //! drawString centres the string on cds.
  virtual void drawString(const std::string &str, const Point2D &cds);
  // unless the specific drawer over-rides this overload, it will just call
  // the first one.  SVG for one needs the alignment flag.
  virtual void drawString(const std::string &str, const Point2D &cds,
                          TextAlignType align);
  //! draw a triangle
  virtual void drawTriangle(const Point2D &cds1, const Point2D &cds2,
                            const Point2D &cds3);
  //! draw an ellipse
  virtual void drawEllipse(const Point2D &cds1, const Point2D &cds2);
  // draw the arc of a circle between ang1 and ang2.  Note that 0 is
  // at 3 o-clock and 90 at 12 o'clock as you'd expect from your maths.
  // ang2 must be > ang1 - it won't draw backwards.  This is not enforced.
  // Angles in degrees.
  virtual void drawArc(const Point2D &centre, double radius, double ang1,
                       double ang2);
  // and a general ellipse form
  virtual void drawArc(const Point2D &centre, double xradius, double yradius,
                       double ang1, double ang2);
  //! draw a rectangle
  virtual void drawRect(const Point2D &cds1, const Point2D &cds2);
  //! draw a line indicating the presence of an attachment point (normally a
  //! squiggle line perpendicular to a bond)
  virtual void drawAttachmentLine(const Point2D &cds1, const Point2D &cds2,
                                  const DrawColour &col, double len = 1.0,
                                  unsigned int nSegments = 16);
  //! draw a wavy line like that used to indicate unknown stereochemistry
  virtual void drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                            const DrawColour &col1, const DrawColour &col2,
                            unsigned int nSegments = 16,
                            double vertOffset = 0.05);
  //! adds additional information about the atoms to the output. Does not make
  //! sense for all renderers.
  virtual void tagAtoms(const ROMol &mol) { RDUNUSED_PARAM(mol); };
  //! set whether or not polygons are being filled
  virtual bool fillPolys() const { return fill_polys_; }
  //! returns either or not polygons should be filled
  virtual void setFillPolys(bool val) { fill_polys_ = val; }

  //! returns our current drawing options
  MolDrawOptions &drawOptions() { return options_; }
  //! \overload
  const MolDrawOptions &drawOptions() const { return options_; }

  //! returns the coordinates of the atoms of the current molecule in molecular
  //! coordinates
  const std::vector<Point2D> &atomCoords() const {
    PRECONDITION(activeMolIdx_ >= 0, "no index");
    return at_cds_[activeMolIdx_];
  };
  //! returns the atomic symbols of the current molecule
  const std::vector<std::pair<std::string, OrientType>> &atomSyms() const {
    PRECONDITION(activeMolIdx_ >= 0, "no index");
    return atom_syms_[activeMolIdx_];
  };
  //! Draw an arrow with either lines or a filled head (when asPolygon is true)
  virtual void drawArrow(const Point2D &cds1, const Point2D &cds2,
                         bool asPolygon = false, double frac = 0.05,
                         double angle = M_PI / 6);

  // reset to default values all the things the c'tor sets
  void tabulaRasa();

  virtual bool supportsAnnotations() { return true; }

 protected:
  std::unique_ptr<DrawText> text_drawer_;

 private:
  bool needs_scale_;
  int width_, height_, panel_width_, panel_height_, legend_height_;
  double scale_;
  double x_min_, y_min_, x_range_, y_range_;
  double x_trans_, y_trans_;
  int x_offset_, y_offset_;  // translation in screen coordinates
  bool fill_polys_;
  int activeMolIdx_;

  DrawColour curr_colour_;
  DashPattern curr_dash_;
  MolDrawOptions options_;

  std::vector<std::vector<Point2D>> at_cds_;  // from mol
  std::vector<std::vector<int>> atomic_nums_;
  std::vector<std::vector<std::pair<std::string, OrientType>>> atom_syms_;
  // by the time atom_notes_ and bonds_notes_ are drawn, we're only ever
  // using the trans_ member of the StringRect, but it is convenient to
  // keep the whole thing rather than just a StringPos for the position
  // for calculating the scale of the drawing.  Went a long way down
  // the rabbit hole before realising this, hence this note.
  std::vector<std::vector<std::shared_ptr<StringRect>>> atom_notes_;
  std::vector<std::vector<std::shared_ptr<StringRect>>> bond_notes_;
  std::vector<std::vector<std::pair<std::shared_ptr<StringRect>, OrientType>>>
      radicals_;

  Point2D bbox_[2];

  // return a DrawColour based on the contents of highlight_atoms or
  // highlight_map, falling back to atomic number by default
  DrawColour getColour(
      int atom_idx, const std::vector<int> *highlight_atoms = nullptr,
      const std::map<int, DrawColour> *highlight_map = nullptr);
  DrawColour getColourByAtomicNum(int atomic_num);

  // set the system up to draw the molecule including calculating the scale.
  std::unique_ptr<RWMol> setupDrawMolecule(
      const ROMol &mol, const std::vector<int> *highlight_atoms,
      const std::map<int, double> *highlight_radii, int confId, int width,
      int height);
  // copies of atom coords, atomic symbols etc. are stashed for convenience.
  // these put empty collections onto the stack and pop the off when done.
  void pushDrawDetails();
  void popDrawDetails();

  // do the initial setup bits for drawing a molecule.
  std::unique_ptr<RWMol> setupMoleculeDraw(
      const ROMol &mol, const std::vector<int> *highlight_atoms,
      const std::map<int, double> *highlight_radii, int confId = -1);
  void setupTextDrawer();

  // if bond_colours is given, it must have an entry for every bond, and it
  // trumps everything else.  First in pair is bonds begin atom, second is
  // end atom.
  void drawBonds(const ROMol &draw_mol,
                 const std::vector<int> *highlight_atoms = nullptr,
                 const std::map<int, DrawColour> *highlight_atom_map = nullptr,
                 const std::vector<int> *highlight_bonds = nullptr,
                 const std::map<int, DrawColour> *highlight_bond_map = nullptr,
                 const std::vector<std::pair<DrawColour, DrawColour>>
                     *bond_colours = nullptr);
  // do the finishing touches to the drawing
  void finishMoleculeDraw(const ROMol &draw_mol,
                          const std::vector<DrawColour> &atom_colours);
  void drawLegend(const std::string &legend);
  // draw a circle in the requested colour(s) around the atom.
  void drawHighlightedAtom(int atom_idx, const std::vector<DrawColour> &colours,
                           const std::map<int, double> *highlight_radii);
  // calculate the rectangle that goes round the string, taking its
  // orientation into account.  Centre of StringRect
  // won't be the same as label_coords, necessarily, as the string might
  // be offset according to orient.
  StringRect calcLabelRect(const std::string &label, OrientType orient,
                           const Point2D &label_coords) const;
  // calculate parameters for an ellipse that roughly goes round the label
  // of the given atom.
  void calcLabelEllipse(int atom_idx,
                        const std::map<int, double> *highlight_radii,
                        Point2D &centre, double &xradius,
                        double &yradius) const;
  // these both assume there is a note on the atom or bond.  That should
  // have been checked by the calling function. StringRect will have a
  // width of -1.0 if there's a problem.
  StringRect calcAnnotationPosition(const ROMol &mol, const Atom *atom);
  StringRect calcAnnotationPosition(const ROMol &mol, const Bond *bond);
  // find where to put the given annotation around an atom.  Starting
  // search at angle start_ang, in degrees.
  void calcAtomAnnotationPosition(const ROMol &mol, const Atom *atom,
                                  double start_ang, StringRect &rect);

  // draw 1 or more coloured line along bonds
  void drawHighlightedBonds(
      const ROMol &mol,
      const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
      const std::map<int, int> &highlight_linewidth_multipliers,
      const std::map<int, double> *highlight_radii);
  int getHighlightBondWidth(
      int bond_idx,
      const std::map<int, int> *highlight_linewidth_multipliers) const;
  // move p2 so that the line defined by p1 to p2 touches the ellipse for the
  // atom highlighted.
  void adjustLineEndForHighlight(int at_idx,
                                 const std::map<int, double> *highlight_radii,
                                 Point2D p1, Point2D &p2) const;

  void extractAtomCoords(const ROMol &mol, int confId, bool updateBBox);
  void extractAtomSymbols(const ROMol &mol);
  void extractAtomNotes(const ROMol &mol);
  void extractBondNotes(const ROMol &mol);
  void extractRadicals(const ROMol &mol);

  // coords in atom coords
  virtual void drawLine(const Point2D &cds1, const Point2D &cds2,
                        const DrawColour &col1, const DrawColour &col2);
  void drawWedgedBond(const Point2D &cds1, const Point2D &cds2,
                      bool draw_dashed, const DrawColour &col1,
                      const DrawColour &col2);
  // draw an arrow for a dative bond, with the arrowhead at cds2.
  void drawDativeBond(const Point2D &cds1, const Point2D &cds2,
                      const DrawColour &col1, const DrawColour &col2);
  void drawAtomLabel(int atom_num,
                     const std::vector<int> *highlight_atoms = nullptr,
                     const std::map<int, DrawColour> *highlight_map = nullptr);
  OrientType calcRadicalRect(const ROMol &mol, const Atom *atom,
                             StringRect &rad_rect);
  void drawRadicals(const ROMol &mol);
  // find a good starting point for scanning round the annotation
  // atom.  If we choose well, the first angle should be the one.
  // Returns angle in radians.
  double getNoteStartAngle(const ROMol &mol, const Atom *atom) const;
  // see if the note will clash with anything else drawn on the molecule.
  // note_vec should have unit length.  note_rad is the radius along
  // note_vec that the note will be drawn.
  bool doesAtomNoteClash(StringRect &note_rect,
                         const std::vector<std::shared_ptr<StringRect>> &rects,
                         const ROMol &mol, unsigned int atom_idx);
  bool doesBondNoteClash(StringRect &note_rect,
                         const std::vector<std::shared_ptr<StringRect>> &rects,
                         const ROMol &mol, const Bond *bond);
  // does the note_vec form an unacceptably acute angle with one of the
  // bonds from atom to its neighbours.
  bool doesNoteClashNbourBonds(
      const StringRect &note_rect,
      const std::vector<std::shared_ptr<StringRect>> &rects, const ROMol &mol,
      const Atom *atom) const;
  // does the note intersect with atsym, and if not, any other atom symbol.
  bool doesNoteClashAtomLabels(
      const StringRect &note_rect,
      const std::vector<std::shared_ptr<StringRect>> &rects, const ROMol &mol,
      unsigned int atom_idx) const;
  bool doesNoteClashOtherNotes(
      const StringRect &note_rect,
      const std::vector<std::shared_ptr<StringRect>> &rects) const;

  // cds1 and cds2 are 2 atoms in a ring.  Returns the perpendicular pointing
  // into the ring.
  Point2D bondInsideRing(const ROMol &mol, const Bond *bond,
                         const Point2D &cds1, const Point2D &cds2) const;
  // cds1 and cds2 are 2 atoms in a chain double bond.  Returns the
  // perpendicular pointing into the inside of the bond
  Point2D bondInsideDoubleBond(const ROMol &mol, const Bond *bond) const;
  // calculate normalised perpendicular to vector between two coords, such
  // that
  // it's inside the angle made between (1 and 2) and (2 and 3).
  Point2D calcInnerPerpendicular(const Point2D &cds1, const Point2D &cds2,
                                 const Point2D &cds3) const;

  // take the coords for atnum, with neighbour nbr_cds, and move cds out to
  // accommodate
  // the label associated with it.
  void adjustBondEndForLabel(int atnum, const Point2D &nbr_cds,
                             Point2D &cds) const;

  // adds LaTeX-like annotation for super- and sub-script.
  std::pair<std::string, OrientType> getAtomSymbolAndOrientation(
      const Atom &atom) const;
  std::string getAtomSymbol(const Atom &atom, OrientType orientation) const;
  OrientType getAtomOrientation(const Atom &atom) const;

  // things used by calculateScale.
  void adjustScaleForAtomLabels(const std::vector<int> *highlight_atoms,
                                const std::map<int, double> *highlight_radii);
  void adjustScaleForRadicals(const ROMol &mol);
  void adjustScaleForAnnotation(
      const std::vector<std::shared_ptr<StringRect>> &notes);

 private:
  virtual void updateMetadata(const ROMol &mol, int confId) {
    RDUNUSED_PARAM(mol);
    RDUNUSED_PARAM(confId);
  };
  virtual void updateMetadata(const ChemicalReaction &rxn) {
    RDUNUSED_PARAM(rxn);
  };

 protected:
  std::vector<std::pair<std::string, std::string>> d_metadata;
  unsigned int d_numMetadataEntries = 0;

  virtual void doContinuousHighlighting(
      const ROMol &mol, const std::vector<int> *highlight_atoms,
      const std::vector<int> *highlight_bonds,
      const std::map<int, DrawColour> *highlight_atom_map,
      const std::map<int, DrawColour> *highlight_bond_map,
      const std::map<int, double> *highlight_radii);

  virtual void highlightCloseContacts();
  // if bond_colours is given, it must have an entry for every bond, and it
  // trumps everything else.  First in pair is bonds begin atom, second is
  // end atom.
  virtual void drawBond(
      const ROMol &mol, const Bond *bond, int at1_idx, int at2_idx,
      const std::vector<int> *highlight_atoms = nullptr,
      const std::map<int, DrawColour> *highlight_atom_map = nullptr,
      const std::vector<int> *highlight_bonds = nullptr,
      const std::map<int, DrawColour> *highlight_bond_map = nullptr,
      const std::vector<std::pair<DrawColour, DrawColour>> *bond_colours =
          nullptr);
  virtual void drawAtomLabel(int atom_num, const DrawColour &draw_colour);
  virtual void drawAnnotation(const std::string &note,
                              const std::shared_ptr<StringRect> &note_rect);

  // calculate normalised perpendicular to vector between two coords
  Point2D calcPerpendicular(const Point2D &cds1, const Point2D &cds2) const;
  // assuming there's a double bond between atom1 and atom2, calculate
  // the ends of the 2 lines that should be used to draw it, distance
  // offset apart.  Includes bonds of type AROMATIC.
  void calcDoubleBondLines(const ROMol &mol, double offset, const Bond *bond,
                           const Point2D &at1_cds, const Point2D &at2_cds,
                           Point2D &l1s, Point2D &l1f, Point2D &l2s,
                           Point2D &l2f) const;
  // returns true if atom has degree 2 and both bonds are close to
  // linear.
  bool isLinearAtom(const Atom &atom) const;
  // and the same for triple bonds.  One line is from atom to atom,
  // so it doesn't need a separate return.
  void calcTripleBondLines(double offset, const Bond *bond,
                           const Point2D &at1_cds, const Point2D &at2_cds,
                           Point2D &l1s, Point2D &l1f, Point2D &l2s,
                           Point2D &l2f) const;

  // calculate the width to draw a line in draw coords.
  virtual double getDrawLineWidth() const;

  // sort out coords and scale for drawing reactions.
  void get2DCoordsForReaction(ChemicalReaction &rxn, Point2D &arrowBegin,
                              Point2D &arrowEnd, std::vector<double> &plusLocs,
                              double spacing, const std::vector<int> *confIds);
  // despite the name, this is only ever used for molecules in a reaction.
  void get2DCoordsMol(RWMol &mol, double &offset, double spacing, double &maxY,
                      double &minY, int confId, bool shiftAgents,
                      double coordScale);
};

// return true if the line l1s->l1f intersects line l2s->l2f.  If ip is not
// nullptr, the intersection point is stored in it.
RDKIT_MOLDRAW2D_EXPORT bool doLinesIntersect(const Point2D &l1s,
                                             const Point2D &l1f,
                                             const Point2D &l2s,
                                             const Point2D &l2f,
                                             Point2D *ip = nullptr);
// return true if line ls->lf intersects (or is fully inside) the
// rectangle of the string.
RDKIT_MOLDRAW2D_EXPORT bool doesLineIntersectLabel(const Point2D &ls,
                                                   const Point2D &lf,
                                                   const StringRect &lab_rect,
                                                   double padding = 0.0);

}  // namespace RDKit

#endif  // RDKITMOLDRAW2D_H
