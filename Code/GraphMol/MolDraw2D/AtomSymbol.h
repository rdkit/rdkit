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
// It holds the information needed to draw an atom symbol, including
// all the extra bits like isotope labels.

#ifndef RDKIT_ATOMSYMBOL_H
#define RDKIT_ATOMSYMBOL_H

#include <string>

#include <GraphMol/MolDraw2D/DrawText.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

namespace RDKit {

class MolDraw2D;

namespace MolDraw2D_detail {

class DrawMol;
class DrawMolMCH;

class AtomSymbol {
  friend class DrawMol;
  friend class DrawMolMCH;

 public:
  ~AtomSymbol() = default;

  // everything's private because we don't want anyone using it.
 private:
  /*!
   *
   * @param symbol     : the full symbol
   * @param orient     : text orientation (up, down, left, right)
   * @param textDrawer : instance of DrawText to get the character sizes
   * etc.
   */
  AtomSymbol(const std::string &symbol, int atIdx, OrientType orient,
             const Point2D &cds, const DrawColour &colour,
             DrawText &textDrawer);

  AtomSymbol(const AtomSymbol &) = delete;
  AtomSymbol(AtomSymbol &&) = delete;
  AtomSymbol &operator=(const AtomSymbol &) = delete;
  AtomSymbol &operator=(AtomSymbol &&) = delete;

  std::string symbol_;
  int atIdx_;
  OrientType orient_;
  Point2D cds_;
  DrawColour colour_;
  DrawText &textDrawer_;

  std::vector<std::shared_ptr<StringRect>> rects_;
  std::vector<TextDrawType> drawModes_;
  std::vector<char> drawChars_;

  // expects xmin etc to be initialised to something sensible.
  virtual void findExtremes(double &xmin, double &xmax, double &ymin,
                            double &ymax) const;
  // scaleFactor moves the cds_, but the fontScaleFactor changes rects_, because
  // we might be scaling the font differently from the drawing as a whole.
  virtual void scale(const Point2D &scaleFactor);
  virtual void move(const Point2D &trans);
  void draw(MolDraw2D &molDrawer) const;
  bool doesRectClash(const StringRect &rect, double padding) const;

  // this is for debugging almost always.
  void drawRects(MolDraw2D &molDrawer) const;
};

}  // namespace MolDraw2D_detail
} // namespace RDKit

#endif  // RDKIT_ATOMSYMBOL_H
