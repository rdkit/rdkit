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
}  // namespace MolDraw2D_detail

struct DrawColour {
  double r = 0.0, g = 0.0, b = 0.0, a = 1.0;
  DrawColour() = default;
  DrawColour(double r, double g, double b, double a = 1.0)
      : r(r), g(g), b(b), a(a) {}
  bool operator==(const DrawColour &other) const {
    return r == other.r && g == other.g && b == other.b && a == other.a;
  }
  bool operator!=(const DrawColour &other) const { return !(*this == other); }
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

typedef std::map<int, DrawColour> ColourPalette;
typedef std::vector<double> DashPattern;

// This is used to convert the line width into something that SVG and
// Cairo use.  It was picked by eye, and was formerly hidden in
// MolDraw2D::getDrawLineWidth().
static const double lineWidthScaleFactor = 0.02;

//! use the RDKit's default palette r
// 201 is for hydrogens when atom symbols are not being drawn.
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
  palette[201] = DrawColour(0.68, 0.85, 0.90);
};

//! use the color palette from the Avalon renderer
// 201 is for hydrogens when atom symbols are not being drawn.
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
  palette[201] = DrawColour(0.68, 0.85, 0.90);
};

//! use (part of) the CDK color palette
/*!
  data source:
  https://github.com/cdk/cdk/blob/master/display/render/src/main/java/org/openscience/cdk/renderer/color/CDK2DAtomColors.java
*/
// 201 is for hydrogens when atom symbols are not being drawn.
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
  palette[201] = DrawColour(0.68, 0.85, 0.90);
};

// 201 is for hydrogens when atom symbols are not being drawn.
inline void assignDarkModePalette(ColourPalette &palette) {
  palette.clear();
  palette[-1] = DrawColour(0.8, 0.8, 0.8);
  palette[0] = DrawColour(0.9, 0.9, 0.9);
  palette[1] = palette[6] = DrawColour(0.9, 0.9, 0.9);
  palette[7] = DrawColour(0.33, 0.41, 0.92);
  palette[8] = DrawColour(1.0, 0.2, 0.2);
  palette[9] = DrawColour(0.2, 0.8, 0.8);
  palette[15] = DrawColour(1.0, 0.5, 0.0);
  palette[16] = DrawColour(0.8, 0.8, 0.0);
  palette[17] = DrawColour(0.0, 0.802, 0.0);
  palette[35] = DrawColour(0.71, 0.4, 0.07);
  palette[53] = DrawColour(0.89, 0.004, 1);
  palette[201] = DrawColour(0.68, 0.85, 0.90);
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
  DrawColour highlightColour{1.0, 0.5, 0.5, 1.0};  // default highlight color
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
      1.0, 1.0, 1.0, 1.0};  // color to be used while clearing the background
  DrawColour queryColour{0.5, 0.5, 0.5,
                         1.0};  // color to be used for query bonds
  int legendFontSize = 16;  // font size (in pixels) to be used for the legend
                            // (if present)
  double legendFraction =
      0.1;  // fraction of the draw panel to be used for the legend if present
  int maxFontSize = 40;  // maximum size in pixels for font in drawn molecule.
                         // -1 means no max.
  int minFontSize = 6;   // likewise for -1.
  int fixedFontSize =
      -1;  // font size to use, in pixels.  Default -1 means not fixed.  If set,
           // always used irrespective of scale, minFontSize and maxFontSize.
  double annotationFontScale = 0.5;  // scales font relative to atom labels for
                                     // atom and bond annotation.
  std::string fontFile = "";  // name of font file for freetype rendering.  If
                              // given, over-rides default
                              // (BuiltinTelexRegular).  Can also be
                              // BuiltinRobotoRegular.
  DrawColour legendColour{0, 0,
                          0};  // color to be used for the legend (if present)
  double multipleBondOffset = 0.15;  // offset for the extra lines
                                     // in a multiple bond as a fraction of
                                     // mean bond length
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
      0.0, 0.0, 0.0,
      1.0};  // color to be used for the symbols and arrows in reactions
  DrawColour annotationColour{0.0, 0.0, 0.0,
                              1.0};  // color to be used for annotations
  double bondLineWidth = 2.0;        // default line width when drawing bonds
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
  bool unspecifiedStereoIsUnknown = false;  // if true, double bonds with
                                            // unspecified stereo are drawn
                                            // crossed, potential stereocenters
                                            // with unspecified stereo are drawn
                                            // with a wavy bond.
  bool singleColourWedgeBonds =
      false;                        // if true wedged and dashed bonds are drawn
                                    // using symbolColour rather than inheriting
                                    // their colour from the atoms
  bool useMolBlockWedging = false;  // If the molecule came from a MolBlock,
                                    // prefer the wedging information that
                                    // provides.  If false, use RDKit rules.
  double scalingFactor = 20.0;      // scaling factor used for pixels->angstrom
                                    // when auto scaling is being used
  double baseFontSize =
      -1.0;  // when > 0 this is used to set the baseFontSize used for text
             // drawing. As a reference point: the default value for
             // DrawText::baseFontSize  is 0.6
  bool drawMolsSameScale = true;  // when drawing multiple molecules with
                                  // DrawMolecules, forces them to use the same
                                  // scale.  Default is true.
  bool useComplexQueryAtomSymbols =
      true;  // replace any atom, any hetero, any halo queries
             // with complex query symbols A, Q, X, M, optionally followed
             // by H if hydrogen is included (except for AH, which stays *).
             // Default is true.

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

}  // namespace RDKit

#endif  // RDKIT_MOLDRAW2DHELPERS_H
