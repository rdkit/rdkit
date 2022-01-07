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

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <GraphMol/MolDraw2D/AtomSymbol.h>
#include <GraphMol/MolDraw2D/DrawText.h>

namespace RDKit {

// ****************************************************************************
AtomSymbol::AtomSymbol(const std::string &symbol, int atIdx, OrientType orient,
                     const Point2D &cds, const DrawColour &colour,
                     DrawText &textDrawer)
    : symbol_(symbol),
      atIdx_(atIdx),
      orient_(orient),
      cds_(cds),
      colour_(colour),
      textDrawer_(textDrawer) {
  if (symbol_.empty()) {
    return;
  }
  textDrawer_.getStringRects(symbol_, orient_, rects_, draw_modes_,
                             draw_chars_, false, TextAlignType::MIDDLE);
}

// ****************************************************************************
void AtomSymbol::findExtremes(double &xmin, double &xmax, double &ymin,
                             double &ymax) const {
  Point2D bl, br, tl, tr, origTrans;
  for (auto rect : rects_) {
    origTrans = rect->trans_;
    rect->trans_ += cds_;
    rect->calcCorners(tl, tr, br, bl, 0.0);
    // sometimes the rect is in a coordinate frame where +ve y is down,
    // sometimes it's up.  For these purposes, we don't care so long as
    // the ymax is larger than the ymin.  We probably don't need to do
    // all the tests for xmin and xmax.
    xmin = std::min(bl.x, xmin);
    xmin = std::min(tr.x, xmin);
    ymin = std::min(bl.y, ymin);
    ymin = std::min(tr.y, ymin);
    xmax = std::max(bl.x, xmax);
    xmax = std::max(tr.x, xmax);
    ymax = std::max(bl.y, ymax);
    ymax = std::max(tr.y, ymax);
    rect->trans_ = origTrans;
  }
}

// ****************************************************************************
void AtomSymbol::scale(const Point2D &scaleFactor,
                       const Point2D &fontScaleFactor) {
  cds_.x *= scaleFactor.x;
  cds_.y *= scaleFactor.y;
  for (auto &rect : rects_) {
    rect->trans_.x *= fontScaleFactor.x;
    rect->trans_.y *= fontScaleFactor.y;
    rect->offset_.x *= fontScaleFactor.x;
    rect->offset_.y *= fontScaleFactor.y;
    rect->g_centre_.x *= fontScaleFactor.x;
    rect->g_centre_.y *= fontScaleFactor.y;
    rect->y_shift_ *= fontScaleFactor.y;
    rect->width_ *= fontScaleFactor.x;
    rect->height_ *= fontScaleFactor.y;
    rect->rect_corr_ *= fontScaleFactor.y;
  }
}

// ****************************************************************************
void AtomSymbol::move(const Point2D &trans) {
  cds_ += trans;
}

// ****************************************************************************
void AtomSymbol::draw(MolDraw2D &molDrawer) const {
  std::string o_class = molDrawer.getActiveClass();
  std::string actClass = o_class;
  if (!actClass.empty()) {
    actClass += " ";
  }
  actClass += boost::str(boost::format("atom-%d") % atIdx_);
  molDrawer.setActiveClass(actClass);
  textDrawer_.setColour(colour_);
  textDrawer_.drawString(symbol_, cds_, orient_);
  molDrawer.setActiveClass(o_class);
//  drawRects(molDrawer);
}

// ****************************************************************************
bool AtomSymbol::doesRectClash(const StringRect &rect, double padding) const {
  for (auto &alrect : rects_) {
    StringRect r{*alrect};
    r.trans_ += cds_;
    if (r.doesItIntersect(rect, padding)) {
      return true;
    }
  }
  return false;
}

// ****************************************************************************
void AtomSymbol::drawRects(MolDraw2D &molDrawer) const {
  Point2D tl, tr, br, bl, origTrans;
  for (auto &rect : rects_) {
    origTrans = rect->trans_;
    rect->trans_ += cds_;
    rect->calcCorners(tl, tr, br, bl, 0.0);
    molDrawer.setColour(DrawColour(1.0, 0.0, 0.0));
    molDrawer.drawLine(tl, tr, true);
    molDrawer.setColour(DrawColour(0.0, 1.0, 0.0));
    molDrawer.drawLine(tr, br, true);
    molDrawer.setColour(DrawColour(0.0, 0.0, 1.0));
    molDrawer.drawLine(br, bl, true);
    molDrawer.setColour(DrawColour(0.0, 0.95, 0.95));
    molDrawer.drawLine(bl, tl, true);
    rect->trans_ = origTrans;
  }
}

} // namespace RDKit

