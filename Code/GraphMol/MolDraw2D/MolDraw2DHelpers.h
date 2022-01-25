//
//  Copyright (C) 2014-2021 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//
// A load of helper classes used by MolDraw2D.

#ifndef RDKIT_MOLDRAW2DHELPERS_H
#define RDKIT_MOLDRAW2DHELPERS_H

#include <Geometry/point.h>

using RDGeom::Point2D;

namespace RDKit {

namespace MolDraw2D_detail {
// for aligning the drawing of text to the passed in coords.
enum class OrientType : unsigned char { C = 0, N, E, S, W };
enum class TextAlignType : unsigned char { MIDDLE = 0, START, END };
}

struct DrawColour {
  double r = 0.0, g = 0.0, b = 0.0, a = 1.0;
  DrawColour() = default;
  DrawColour(double r, double g, double b, double a = 1.0)
      : r(r), g(g), b(b), a(a) {}
  bool operator==(const DrawColour &other) const {
    return r == other.r && g == other.g && b == other.b && a == other.a;
  }
  bool feq(const DrawColour &other, double tol = 0.001,
           bool ignoreAlpha = true) const {
    return fabs(r - other.r) <= tol && fabs(g - other.g) <= tol &&
           fabs(b - other.b) <= tol &&
           (ignoreAlpha || fabs(a - other.a) <= tol);
  }
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

//! for annotating the type of the extra shapes
enum class MolDrawShapeType {
  Arrow,  // ordering of points is: start, end, p1, p2
  Polyline,
  Ellipse,
};

//! extra shape to add to canvas
struct MolDrawShape {
  MolDrawShapeType shapeType = MolDrawShapeType::Polyline;
  std::vector<Point2D> points;
  DrawColour lineColour{0, 0, 0};
  int lineWidth = 2;
  bool fill = false;
  bool scaleLineWidth = false;
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
    Point2D c;
    calcCentre(c);
    tl = Point2D(c.x - wb2, c.y - hb2);
    tr = Point2D(c.x + wb2, c.y - hb2);
    br = Point2D(c.x + wb2, c.y + hb2);
    bl = Point2D(c.x - wb2, c.y + hb2);
  }
  void calcCentre(Point2D &c) const {
    c = trans_ + g_centre_ - offset_;
    c.y -= y_shift_;
  }
  bool doesItIntersect(const StringRect &other, double padding = 0.0) const {
    Point2D ttl, ttr, tbr, tbl;
    calcCorners(ttl, ttr, tbr, tbl, padding);
    // is +ve y up or down?
    if (ttl.y < tbl.y) {
      std::swap(ttl, tbl);
      std::swap(ttr, tbr);
    }
    Point2D otl, otr, obr, obl;
    other.calcCorners(otl, otr, obr, obl, padding);
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
typedef std::vector<double> DashPattern;

//! use the RDKit's default palette r
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

//! use the color palette from the Avalon renderer
inline void assignAvalonPalette(ColourPalette &palette) {
  palette.clear();
  palette[-1] = DrawColour(0, 0, 0);
  palette[0] = DrawColour(0.1, 0.1, 0.1);
  palette[1] = palette[6] = DrawColour(0.0, 0.0, 0.0);
  palette[7] = DrawColour(0.0, 0.0, 1.0);
  palette[8] = DrawColour(1.0, 0.0, 0.0);
  palette[9] = DrawColour(0.0, 0.498, 0.0);
  palette[15] = DrawColour(0.498, 0.0, 0.498);
  palette[16] = DrawColour(0.498, 0.247, 0.0);
  palette[17] = DrawColour(0.0, 0.498, 0.0);
  palette[35] = DrawColour(0.0, 0.498, 0.0);
  palette[53] = DrawColour(0.247, 0.0, 0.498);
};

//! use (part of) the CDK color palette
/*!
  data source:
  https://github.com/cdk/cdk/blob/master/display/render/src/main/java/org/openscience/cdk/renderer/color/CDK2DAtomColors.java
*/
inline void assignCDKPalette(ColourPalette &palette) {
  palette.clear();
  palette[-1] = DrawColour(0, 0, 0);
  palette[0] = DrawColour(0.1, 0.1, 0.1);
  palette[1] = palette[6] = DrawColour(0.0, 0.0, 0.0);
  palette[7] = DrawColour(0.188, 0.314, 0.972);
  palette[8] = DrawColour(1.0, 0.051, 0.051);
  palette[9] = DrawColour(0.565, 0.878, 0.314);
  palette[15] = DrawColour(1.0, 0.5, 0.0);
  palette[16] = DrawColour(0.776, 0.776, 0.173);
  palette[17] = DrawColour(0.122, 0.498, 0.122);
  palette[35] = DrawColour(0.651, 0.161, 0.161);
  palette[53] = DrawColour(0.580, 0.0, 0.580);
  palette[5] = DrawColour(1.000, 0.710, 0.710);
};

inline void assignDarkModePalette(ColourPalette &palette) {
  palette.clear();
  palette[-1] = DrawColour(0.8, 0.8, 0.8);
  palette[0] = DrawColour(0.9, 0.9, 0.9);
  palette[1] = palette[6] = DrawColour(0.9, 0.9, 0.9);
  palette[7] = DrawColour(0.2, 0.2, 1.0);
  palette[8] = DrawColour(1.0, 0.2, 0.2);
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
  bool splitBonds = false;             // split bonds into per atom segments
                                                 // most useful for dynamic manipulation of drawing
                                                 // especially for svg
  DrawColour highlightColour{1, 0.5, 0.5, 1.0};  // default highlight color
  bool continuousHighlight = true;  // highlight by drawing an outline
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
      1, 1, 1, 1};          // color to be used while clearing the background
  int legendFontSize = 16;  // font size (in pixels) to be used for the legend
                                     // (if present)
  double legendFraction =
      0.1;  // fraction of the draw panel to be used for the legend if present
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
  bool noAtomLabels =
      false;  // disables inclusion of atom labels in the rendering
  std::vector<std::vector<int>> atomRegions;  // regions
  DrawColour symbolColour{
      0, 0, 0, 1};  // color to be used for the symbols and arrows in reactions
  DrawColour annotationColour{0, 0, 0, 1};  // color to be used for annotations
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
  bool addAtomIndices = false;     // adds atom indices to drawings.
  bool addBondIndices = false;     // adds bond indices to drawings.
  bool isotopeLabels = true;       // adds isotope to non-dummy atoms.
  bool dummyIsotopeLabels = true;  // adds isotope labels to dummy atoms.

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
  bool comicMode = false;  // simulate hand-drawn lines for bonds. When combined
                                         // with a font like Comic-Sans or Comic-Neue, this
                                         // gives xkcd-like drawings.
  int variableBondWidthMultiplier = 16;  // what to multiply standard bond width
                                         // by for variable attachment points.
  double variableAtomRadius = 0.4;  // radius value to use for atoms involved in
                                         // variable attachment points.
  DrawColour variableAttachmentColour = {
      0.8, 0.8, 0.8, 1.0};  // colour to use for variable attachment points
  bool includeChiralFlagLabel =
      false;  // add a molecule annotation with "ABS" if the chiral flag is set
  bool simplifiedStereoGroupLabel =
      false;  // if all specified stereocenters are in a single StereoGroup,
              // show a molecule-level annotation instead of the individual
              // labels
  bool singleColourWedgeBonds =
      false;  // if true wedged and dashed bonds are drawn
              // using symbolColour rather than inheriting
              // their colour from the atoms
  double scalingFactor = 20.0;  // scaling factor used for pixels->angstroms
                                // when auto scaling is being used
  double baseFontSize =
      -1.0;  // when > 0 this is used to set the baseFontSize used for text
             // drawing. As a reference point: the default value for
             // DrawText::baseFontSize  is 0.6

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
  }
};

} // namespace RDKit

#endif  // RDKIT_MOLDRAW2DHELPERS_H
