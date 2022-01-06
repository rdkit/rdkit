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
// This class is a helper class used by MolDraw2D to draw an ROMol.
// It is not part of the public API and is not intended to be used except
// by MolDraw2D.

#ifndef RDKIT_DRAWMOL_H
#define RDKIT_DRAWMOL_H

#include <map>
#include <string>
#include <vector>

#include <Geometry/point.h>
#include <GraphMol/MolDraw2D/DrawAnnotation.h>
#include <GraphMol/MolDraw2D/DrawShape.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

namespace RDKit {

class Atom;
class Bond;
class ROMol;
class RWMol;
class AtomSymbol;
class DrawText;

class DrawMol {
  friend class MolDraw2D;

  // everything's private because we don't want anyone using it.
 protected:
  // Make the object, scaled to a given pixel size.
  /*!
    \param mol             : the molecule to draw
    \param legend          : the legend (to be drawn under the molecule)
    \param width           : width (in pixels) of the rendering
    set this to -1 to have the canvas size set automatically
    \param height          : height (in pixels) of the rendering
    set this to -1 to have the canvas size set automatically
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
  DrawMol(const ROMol &mol, const std::string &legend, int width, int height,
          MolDrawOptions &drawOptions, DrawText &textDrawer,
          const std::vector<int> *highlightAtoms = nullptr,
          const std::vector<int> *highlightBonds = nullptr,
          const std::map<int, DrawColour> *highlightAtomMap = nullptr,
          const std::map<int, DrawColour> *highlightBondMap = nullptr,
          const std::vector<std::pair<DrawColour, DrawColour>> *bondColours =
              nullptr,
          const std::map<int, double> *highlight_radii = nullptr,
          bool includeAnnotations = true, int confId = -1);
  DrawMol(const DrawMol &) = delete;
  DrawMol(const DrawMol &&) = delete;
  DrawMol &operator=(const DrawMol &) = delete;

  // this must be called before a drawing can be done.
  void createDrawObjects();
  // common bits used by createDrawObjects and setScale.
  void finishCreateDrawObjects();
  void initDrawMolecule(const ROMol &mol);
  void extractAll();
  void extractAtomCoords();
  void extractAtomSymbols();
  void extractBonds();
  virtual void extractHighlights();
  void extractRegions();
  void extractAttachments();
  void extractMolNotes();
  void extractAtomNotes();
  void extractBondNotes();
  void extractRadicals();
  void extractSGroupData();
  void extractVariableBonds();
  void extractBrackets();
  void extractLinkNodes();
  // extractCloseContacts is to show where 2 atoms are drawn too close together
  // and so needs the final drawing coords.  It is therefore called at the end
  // of changeToDrawCoords() and any necessary DrawShapePolyLines added to
  // postShapes_ in drawCoords.
  void extractCloseContacts();
  void calculateScale();
  void findExtremes();
  void changeToDrawCoords();
  void draw(MolDraw2D &drawer) const;
  void drawRadicals(MolDraw2D &drawer) const;
  void resetEverything();

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

  void calcMeanBondLengthSquare();
  void makeStandardBond(Bond *bond, double doubleBondOffset);
  void makeQueryBond(Bond *bond, double doubleBondOffset);
  void makeDoubleBondLines(Bond *bond, double doubleBondOffset,
                           const std::pair<DrawColour, DrawColour> &cols);
  void makeTripleBondLines(Bond *bond, double doubleBondOffset,
                           const std::pair<DrawColour, DrawColour> &cols);
  void makeWedgedBond(Bond *bond,
                      const std::pair<DrawColour, DrawColour> &cols);
  void makeWavyBond(Bond *bond, const std::pair<DrawColour, DrawColour> &cols);
  void makeDativeBond(Bond *bond,
                      const std::pair<DrawColour, DrawColour> &cols);
  void makeZeroBond(Bond *bond, const std::pair<DrawColour, DrawColour> &cols,
                    const DashPattern &dashPattern);
  void adjustBondEndsForLabels(int begAtIdx, int endAtIdx, Point2D &begCds,
                               Point2D &endCds);
  void newBondLine(const Point2D &pt1, const Point2D &pt2,
                   const DrawColour &col1, const DrawColour &col2, int atom1Idx,
                   int atom2Idx, int bondIdx, const DashPattern &dashPattern);
  std::pair<DrawColour, DrawColour> getBondColours(Bond *bond);
  void makeContinuousHighlights();
  void makeAtomCircleHighlights();
  void makeAtomEllipseHighlights(int lineWidth);
  void makeBondHighlightLines(int lineWidth);
  void calcAnnotationPosition(const Atom *atom, DrawAnnotation &annot) const;
  void calcAnnotationPosition(const Bond *bond, DrawAnnotation &annot) const;
  double getNoteStartAngle(const Atom *atom) const;
  void calcMolNotePosition(const std::vector<Point2D> atCds,
                           DrawAnnotation &annot) const;
  // see if the note will clash with anything else drawn on the molecule.
  // Returns 0 if no clash, 1-4 if there is a clash, denoting what clashed.
  int doesNoteClash(const DrawAnnotation &annot) const;
  int doesRectClash(const StringRect &rect, double padding) const;
  OrientType calcRadicalRect(const Atom *atom, StringRect &rad_rect) const;
  void getDrawTransformers(Point2D &trans, Point2D &scale,
                           Point2D &toCentre) const;
  // given some coords in molecule space (angstrom, probably) return the
  // screen coords
  Point2D getDrawCoords(const Point2D &atCds, const Point2D &trans,
                        const Point2D &scaleFactor,
                        const Point2D &toCentre) const;
  Point2D getDrawCoords(const Point2D &atCds) const;
  // and the other way.
  Point2D getAtomCoords(const Point2D &screenCds) const;
  double getScale() const { return scale_; }
  void setScale(double newScale);
  // for drawing into a grid, for example.  Must be set before
  // changeToDrawCoords is called for it to have effect.
  void setOffsets(double xOffset, double yOffset);
  // so we can add metadata later.  Most likely used after changeToDrawCoords
  // has been called.
  void tagAtomsWithCoords();
  // apply the transformations to everything except the atomLabels_ and
  // annotations_. trans and toCentre are added, scale is multiplied.
  void transformAllButAtomLabels(const Point2D &trans, Point2D &scale,
                                 const Point2D &toCentre);

  MolDrawOptions &drawOptions_;
  DrawText &textDrawer_;
  const std::vector<int> *highlightAtoms_;
  const std::vector<int> *highlightBonds_;
  const std::map<int, DrawColour> *highlightAtomMap_;
  const std::map<int, DrawColour> *highlightBondMap_;
  const std::vector<std::pair<DrawColour, DrawColour>> *bondColours_;
  const std::map<int, double> *highlightRadii_;
  bool includeAnnotations_;
  std::string legend_;

  std::unique_ptr<RWMol> drawMol_;
  int confId_;
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

  int width_, height_;
  // to allow for min and max font sizes, the font scale needs to be
  // independent of the main scale.
  double scale_, fontScale_;
  double xMin_, yMin_, xMax_, yMax_, xRange_, yRange_;
  // offsets are for drawing molecules in grids, for example.
  double xOffset_ = 0.0, yOffset_ = 0.0;
  double meanBondLengthSquare_ = 0.0;
  // if there's a legend, we reserve a bit for it.
  int drawHeight_, legendHeight_;
  bool drawingInitialised_ = false;
};

void centerMolForDrawing(RWMol &mol, int confId = 1);
void prepareStereoGroups(RWMol &mol);
bool isLinearAtom(const Atom &atom, const std::vector<Point2D> &atCds);
std::string getAtomListText(const Atom &atom);
DrawColour getColour(int atom_idx, const MolDrawOptions &drawOptions,
                     const std::vector<int> &atomicNums,
                     const std::vector<int> *highlightAtoms,
                     const std::map<int, DrawColour> *highlightMap);
DrawColour getColourByAtomicNum(int atomicNum,
                                const MolDrawOptions &drawOptions);
int getHighlightBondWidth(
    MolDrawOptions &drawOptions, int bond_idx,
    const std::map<int, int> *highlight_linewidth_multipliers);

void calcDoubleBondLines(const ROMol &mol, double offset, const Bond &bond,
                         const std::vector<Point2D> &at_cds, Point2D &l1s,
                         Point2D &l1f, Point2D &l2s, Point2D &l2f);
void calcTripleBondLines(double offset, const Bond &bond,
                         const std::vector<Point2D> &at_cds, Point2D &l1s,
                         Point2D &l1f, Point2D &l2s, Point2D &l2f);
Point2D calcPerpendicular(const Point2D &cds1, const Point2D &cds2);
Point2D calcInnerPerpendicular(const Point2D &cds1, const Point2D &cds2,
                               const Point2D &cds3);
Point2D bondInsideRing(const ROMol &mol, const Bond &bond, const Point2D &cds1,
                       const Point2D &cds2, const std::vector<Point2D> &at_cds);
Point2D bondInsideDoubleBond(const ROMol &mol, const Bond &bond,
                             const std::vector<Point2D> &at_cds);
// return a point that is end1 moved so as not to clash with any of the
// rects of a label.  end1 to end2 and the coords of 2 ends of a bond.
void adjustBondEndForString(
    const Point2D &end1, const Point2D &end2, double padding,
    const std::vector<std::shared_ptr<StringRect>> &rects, Point2D &moveEnd);
void findRadicalExtremes(
    const std::vector<std::tuple<StringRect, OrientType, int>> &radicals,
    double &xmin, double &xmax, double &ymin, double &ymax);
void findRectExtremes(const StringRect &rect, const TextAlignType &align,
                      double &xmin, double &xmax, double &ymin, double &ymax);
}  // namespace RDKit

#endif  // RDKIT_DRAWMOL_H
