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
// This class is a helper class used by MolDraw2D to draw an ROMol.
// It is not part of the public API and is not intended to be used except
// by MolDraw2D.

#ifndef RDKIT_DRAWMOL_H
#define RDKIT_DRAWMOL_H

#include <map>
#include <string>
#include <vector>

#include <Geometry/point.h>
#include <GraphMol/MolDraw2D/AtomSymbol.h>
#include <GraphMol/MolDraw2D/DrawAnnotation.h>
#include <GraphMol/MolDraw2D/DrawShape.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

namespace RDKit {

class Atom;
class Bond;
class ROMol;
class RWMol;
class DrawText;

namespace MolDraw2D_detail {

class DrawMol {
 public:
  virtual ~DrawMol() = default;

  // Make the object, scaled to a given pixel size.
  /*!
    \param mol             : the molecule to draw
    \param legend          : the legend (to be drawn under the molecule)
    \param width           : width (in pixels) of the rendering
    set this to -1 to have the canvas size set automatically
    \param height          : height (in pixels) of the rendering
    set this to -1 to have the canvas size set automatically
    \param drawOptions    : a MolDrawOptions object from the owning MolDraw2D
    \param textDrawer     : a DrawText object from the owning MolDraw2D
    \param highlightAtoms : (optional) vector of atom ids to highlight
    \param highlightBonds : (optional) vector of bond ids to highlight
    \param highlightAtomMap   : (optional) map from atomId -> DrawColour
    providing the highlight colors. If not provided the default highlight colour
    from \c drawOptions() will be used.
    \param highlightBondMap   : (optional) map from bondId -> DrawColour
    providing the highlight colors. If not provided the default highlight colour
    from \c drawOptions() will be used.
    \param highlightRadii : (optional) map from atomId -> radius (in molecule
    coordinates) for the radii of atomic highlights. If not provided the default
    value from \c drawOptions() will be used.
    \param confId          : (optional) conformer ID to be used for atomic
    coordinates
  */
  DrawMol(const ROMol &mol, const std::string &legend, int width, int height,
          const MolDrawOptions &drawOptions, DrawText &textDrawer,
          const std::vector<int> *highlightAtoms = nullptr,
          const std::vector<int> *highlightBonds = nullptr,
          const std::map<int, DrawColour> *highlightAtomMap = nullptr,
          const std::map<int, DrawColour> *highlightBondMap = nullptr,
          const std::vector<std::pair<DrawColour, DrawColour>> *bondColours =
              nullptr,
          const std::map<int, double> *highlight_radii = nullptr,
          bool includeAnnotations = true, int confId = -1,
          bool isReactionMol = false);
  /*!
   Make a DrawMol when there's no molecule to draw, but we still want
   a DrawMol in the MolDraw2D for scale, conversion of atom coords to
   draw coords etc.  And so DrawMol starts sounding like a poor name for
   the class.
   \param width        : width (in pixels) of the rendering
   \param height       : height (in pixels) of the rendering
   \param drawOptions  : a MolDrawOptions object from the owning MolDraw2D
   \param textDrawer   : a DrawText object from the owning MolDraw2D
   \param xmin         : minimum value expected in X
   \param xmax         : miaximum value expected in X
   \param ymin         : minimum value expected in Y
   \param ymax         : maximum value expected in Y
   \param scale        : scale to use
   \param fontscale    : font scale to use
   */
  DrawMol(int width, int height, const MolDrawOptions &drawOptions,
          DrawText &textDrawer, double xmin, double xmax, double ymin,
          double ymax, double scale, double fontscale);
  DrawMol(const DrawMol &) = delete;
  DrawMol(DrawMol &&) = delete;
  DrawMol &operator=(const DrawMol &) = delete;
  DrawMol &operator=(DrawMol &&) = delete;

  // this must be called before a drawing can be done.
  void createDrawObjects();
  // common bits used by createDrawObjects and setScale.
  void finishCreateDrawObjects();
  void initDrawMolecule(const ROMol &mol);
  void extractAll(double scale);
  void extractAtomCoords();
  void extractAtomSymbols();
  void extractBonds();
  virtual void extractHighlights(double scale);
  void extractRegions();
  void extractAttachments();
  void extractMolNotes();
  void extractStereoGroups();
  void extractAtomNotes();
  void extractBondNotes();
  void extractRadicals();
  void extractSGroupData();
  void extractVariableBonds();
  void extractBrackets();
  void extractLinkNodes();
  // extractCloseContacts is to show where 2 atoms are drawn too close
  // together and so needs the final drawing coords.  It is therefore called
  // at the end of changeToDrawCoords() and any necessary DrawShapePolyLines
  // added to postShapes_ in drawCoords.
  void extractCloseContacts();
  void calculateScale();
  void findExtremes();
  void changeToDrawCoords();
  void draw(MolDraw2D &drawer) const;
  void drawRadicals(MolDraw2D &drawer) const;
  void resetEverything();
  // reduce width_ and height_ to just accomodate the Xrange_ and YRange_
  // at the current scale.  Recentres everything.  So the DrawMol takes up
  // no more screen real estate than it needs.
  void shrinkToFit(bool withPadding = true);

  // adds LaTeX-like annotation for super- and sub-script.
  std::pair<std::string, OrientType> getAtomSymbolAndOrientation(
      const Atom &atom) const;
  std::string getAtomSymbol(const Atom &atom, OrientType orientation) const;
  OrientType getAtomOrientation(const Atom &atom) const;

  // if there's a legend, partition the height_ to accommodate it
  void partitionForLegend();
  // extractLegend is different from all the other extract... functions
  // in that it needs to be called once the final scale has been found
  // by calculateScale.
  void extractLegend();

  void calcMeanBondLength();
  void makeStandardBond(Bond *bond, double doubleBondOffset);
  void makeQueryBond(Bond *bond, double doubleBondOffset);
  void makeDoubleBondLines(Bond *bond, double doubleBondOffset,
                           const std::pair<DrawColour, DrawColour> &cols);
  void makeTripleBondLines(Bond *bond, double doubleBondOffset,
                           const std::pair<DrawColour, DrawColour> &cols);
  void makeWedgedBond(Bond *bond,
                      const std::pair<DrawColour, DrawColour> &cols);
  void makeWavyBond(Bond *bond, double offset,
                    const std::pair<DrawColour, DrawColour> &cols);
  void makeDativeBond(Bond *bond, double offset,
                      const std::pair<DrawColour, DrawColour> &cols);
  void makeZeroBond(Bond *bond, const std::pair<DrawColour, DrawColour> &cols,
                    const DashPattern &dashPattern);
  void adjustBondEndsForLabels(int begAtIdx, int endAtIdx, Point2D &begCds,
                               Point2D &endCds) const;
  void newBondLine(const Point2D &pt1, const Point2D &pt2,
                   const DrawColour &col1, const DrawColour &col2, int atom1Idx,
                   int atom2Idx, int bondIdx, const DashPattern &dashPattern);
  std::pair<DrawColour, DrawColour> getBondColours(Bond *bond);
  void makeContinuousHighlights(double scale);
  void makeAtomCircleHighlights();
  void makeAtomEllipseHighlights(double lineWidth);
  // These are the lines for continuous highlights, that are
  // now filled trapezoids rather than fat simple lines.
  void makeBondHighlightLines(double lineWidth, double scale);
  void calcAnnotationPosition(const Atom *atom, DrawAnnotation &annot) const;
  void calcAnnotationPosition(const Bond *bond, DrawAnnotation &annot) const;
  double getNoteStartAngle(const Atom *atom) const;
  // see if the note will clash with anything else drawn on the molecule.
  // Returns 0 if no clash, 1-4 if there is a clash, denoting what clashed.
  int doesNoteClash(const DrawAnnotation &annot) const;
  int doesRectClash(const StringRect &rect, double padding) const;
  OrientType calcRadicalRect(const Atom *atom, StringRect &rad_rect) const;
  void getDrawTransformers(Point2D &trans, Point2D &scale,
                           Point2D &toCentre) const;
  // Given some coords in molecule space (angstrom, probably) return the
  // screen coords.
  Point2D getDrawCoords(const Point2D &atCds, const Point2D &trans,
                        const Point2D &scaleFactor,
                        const Point2D &toCentre) const;
  Point2D getDrawCoords(const Point2D &atCds) const;
  Point2D getDrawCoords(int atnum) const;
  // and the other way.
  Point2D getAtomCoords(const Point2D &screenCds) const;
  Point2D getAtomCoords(int atnum) const;
  double getScale() const { return scale_; }
  double getFontScale() const { return fontScale_; }
  // More often than not, newScale and newFontScale will be the same,
  // but not if minFontScale of maxFontScale have become involved.
  // The newFontScale will be used without checking the min and max.
  void setScale(double newScale, double newFontScale,
                bool ignoreFontLimits = true);
  // Set all the transformation details from the incoming DrawMol to this
  // one, so they can be overlaid properly.  Doesn't change the offsets.
  void setTransformation(const DrawMol &sourceMol);

  // For drawing into a grid, for example.  Must be set before
  // changeToDrawCoords is called for it to have effect.
  void setOffsets(double xOffset, double yOffset);
  // So we can add metadata later.  Most likely used after changeToDrawCoords
  // has been called.
  void tagAtomsWithCoords();
  // Apply the transformations to everything. trans and toCentre are added,
  // scale is multiplied.
  void transformAll(const Point2D *trans = nullptr, Point2D *scale = nullptr,
                    const Point2D *toCentre = nullptr);
  // Apply the transformations to the given point and return a new one.
  Point2D transformPoint(const Point2D &pt, const Point2D *trans = nullptr,
                         Point2D *scale = nullptr,
                         const Point2D *toCentre = nullptr) const;
  void calcDoubleBondLines(double offset, const Bond &bond, Point2D &l1s,
                           Point2D &l1f, Point2D &l2s, Point2D &l2f) const;
  void bondInsideRing(const Bond &bond, double offset, Point2D &l2s,
                      Point2D &l2f) const;
  void bondNonRing(const Bond &bond, double offset, Point2D &l2s,
                   Point2D &l2f) const;
  // deal with terminal double bond between at1 and at2, either to atoms of
  // degree 2 or 3.
  void doubleBondTerminal(Atom *at1, Atom *at2, double offset, Point2D &l1s,
                          Point2D &l1f, Point2D &l2s, Point2D &l2f) const;
  // assuming at[1-3] are atoms where at1 is bonded to at2 and at2 is bonded
  // to at3, find the position of the at2 end of a double bond between at2
  // and at3.  If trunc, it'll be along the vector that bisects the two bonds
  // on the inside, otherwise it's perpendicular to the bond from at1 to at2.
  Point2D doubleBondEnd(unsigned int at1, unsigned int at2, unsigned int at3,
                        double offset, bool trunc) const;
  void calcTripleBondLines(double offset, const Bond &bond, Point2D &l1s,
                           Point2D &l1f, Point2D &l2s, Point2D &l2f);
  // find the vectors of any atoms singly bonded to atom that aren't
  // otherAtom.
  void findOtherBondVecs(const Atom *atom, const Atom *otherAtom,
                         std::vector<Point2D> &otherBondVecs) const;
  void adjustBondsOnSolidWedgeEnds();
  void smoothBondJoins();
  // If doing a continuous highlight, add to points the 2 or 3 points that
  // will be for the end1 end of the highlight.  The final highlight will
  // be a 4-6 sided polygon formed by calling this twice, with the ends in
  // opposite order the second time.
  void makeHighlightEnd(const Atom *end1, const Atom *end2, double lineWidth,
                        const std::vector<Atom *> &end1HighNbrs,
                        std::vector<Point2D> &points);
  DrawColour getColour(int atom_idx) const;

  const MolDrawOptions &drawOptions_;
  DrawText &textDrawer_;
  // For drawing reactions, padding needs to be 0 irrespective
  // of what drawOptions_ does elsewhere, so it is copied from drawOptions_
  // on construction, and is then immune to changes in the outer program.
  double marginPadding_;
  std::vector<int> highlightAtoms_;
  std::vector<int> highlightBonds_;
  std::map<int, DrawColour> highlightAtomMap_;
  std::map<int, DrawColour> highlightBondMap_;
  std::vector<std::pair<DrawColour, DrawColour>> bondColours_;
  std::map<int, double> highlightRadii_;
  bool includeAnnotations_;
  bool isReactionMol_;
  std::string legend_;

  std::unique_ptr<RWMol> drawMol_;
  int confId_;
  // atCds_ are as extracted from the molecule, except that the y is
  // inverted and drawOptions_.rotate is applied.
  std::vector<Point2D> atCds_;
  std::vector<std::unique_ptr<DrawShape>> bonds_;
  std::vector<std::unique_ptr<DrawShape>> preShapes_;
  std::vector<std::unique_ptr<DrawShape>> postShapes_;
  std::vector<int> atomicNums_;
  std::vector<std::pair<std::string, OrientType>> atomSyms_;
  std::vector<std::unique_ptr<AtomSymbol>> atomLabels_;
  std::vector<std::unique_ptr<DrawShape>> highlights_;
  std::vector<std::unique_ptr<DrawAnnotation>> annotations_;
  std::vector<std::unique_ptr<DrawAnnotation>> legends_;
  std::vector<std::tuple<StringRect, OrientType, int>> radicals_;
  std::vector<int> singleBondLines_;

  // The total width and height of the canvas.  Either or both may be
  // -1 initially, in which case they will be calculated so the molecule
  // can be drawn within it at a fixed scale (the so-called flexiCanvas).
  int width_, height_;
  // The width and height of the drawing area within the canvas.  This is
  // width_ and height_ less the drawOptions._padding round the outside.
  // Will always be <= width_, height_.
  int drawWidth_, drawHeight_;
  // to allow for min and max font sizes, the font scale needs to be
  // independent of the main scale.
  double scale_, fontScale_;
  double xMin_, yMin_, xMax_, yMax_, xRange_, yRange_;
  // offsets are for drawing molecules in grids, for example.
  double xOffset_ = 0.0, yOffset_ = 0.0;
  double meanBondLength_ = 0.0;
  // if there's a legend, we reserve a bit for it.  molHeight_ is the
  // bit for drawing the molecule, legendHeight_ the bit under that
  // for the legend.  In pixels.
  int molHeight_, legendHeight_ = 0;
  bool drawingInitialised_ = false;
  // when drawing the atoms and bonds in an SVG, they are given a class
  // via MolDraw2D's activeAtmIdx[12]_ and activeBndIdx.  We don't always want
  // them to start from 0 for atom/bond 0.
  int activeAtmIdxOffset_ = 0, activeBndIdxOffset_ = 0;
  bool flexiCanvasX_ = false;
  bool flexiCanvasY_ = false;
};

void centerMolForDrawing(RWMol &mol, int confId = 1);
bool isLinearAtom(const Atom &atom, const std::vector<Point2D> &atCds);
std::string getAtomListText(const Atom &atom);
DrawColour getColourByAtomicNum(int atomicNum,
                                const MolDrawOptions &drawOptions);
DrawColour getHighlightBondColour(
    const Bond *bond, const MolDrawOptions &drawOptions,
    const std::vector<int> &highlightBonds,
    const std::map<int, DrawColour> &highlightBondMap,
    const std::vector<int> &highlightAtoms,
    const std::map<int, DrawColour> &highlightAtomMap);
double getHighlightBondWidth(
    const MolDrawOptions &drawOptions, int bond_idx,
    const std::map<int, int> *highlight_linewidth_multipliers);

Point2D calcPerpendicular(const Point2D &cds1, const Point2D &cds2);
Point2D calcInnerPerpendicular(const Point2D &cds1, const Point2D &cds2,
                               const Point2D &cds3);
// return a point that is moveEnd moved so as not to clash with any of the
// rects of a label.  moveEnd to end2 are the coords of 2 ends of a bond.
void adjustBondEndForString(
    const Point2D &end2, double padding,
    const std::vector<std::shared_ptr<StringRect>> &rects, Point2D &moveEnd);
void findRadicalExtremes(
    const std::vector<std::tuple<StringRect, OrientType, int>> &radicals,
    double &xmin, double &xmax, double &ymin, double &ymax);
void findRectExtremes(const StringRect &rect, const TextAlignType &align,
                      double &xmin, double &xmax, double &ymin, double &ymax);
void getBondHighlightsForAtoms(const ROMol &mol,
                               const std::vector<int> &highlight_atoms,
                               std::vector<int> &highlight_bonds);
// returns true if the vector at2->at1 points in roughly the opposite
// direction to at3->at4.  Basically, if the dot product is negative.
bool areBondsTrans(const Point2D &at1, const Point2D &at2, const Point2D &at3,
                   const Point2D &at4);
// returns true if the vector at2->at1 points is roughly linear with
// direction of at3->at4.  Basically, if the dot product is 1.0 within the
// given tolerance.
bool areBondsParallel(const Point2D &at1, const Point2D &at2,
                      const Point2D &at3, const Point2D &at4,
                      double tol = 1.0e-4);

// find the nborNum'th neighbour of firstAtom that isn't secondAtom
const Atom *otherNeighbor(const Atom *firstAtom, const Atom *secondAtom,
                          int nborNum, const ROMol &mol);
}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif  // RDKIT_DRAWMOL_H
