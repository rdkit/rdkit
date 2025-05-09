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

#include <RDGeneral/BoostStartInclude.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/MolDraw2D/AtomSymbol.h>
#include <GraphMol/MolDraw2D/DrawText.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>

namespace RDKit {
namespace MolDraw2D_detail {

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
  textDrawer_.getStringRects(symbol_, orient_, rects_, drawModes_, drawChars_,
                             false, TextAlignType::MIDDLE);
  adjustColons();
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
    xmin = std::min({tr.x, bl.x, xmin});
    ymin = std::min({tr.y, bl.y, ymin});
    xmax = std::max({tr.x, bl.x, xmax});
    ymax = std::max({tr.y, bl.y, ymax});
    rect->trans_ = origTrans;
  }
}

// ****************************************************************************
void AtomSymbol::scale(const Point2D &scaleFactor) {
  cds_.x *= scaleFactor.x;
  cds_.y *= scaleFactor.y;
  recalculateRects();
}

// ****************************************************************************
void AtomSymbol::move(const Point2D &trans) { cds_ += trans; }

// ****************************************************************************
void AtomSymbol::recalculateRects() {
  // rebuild the rectangles, because the fontScale may be different,
  // and the widths etc might not scale by the same amount.
  rects_.clear();
  drawModes_.clear();
  drawChars_.clear();
  textDrawer_.getStringRects(symbol_, orient_, rects_, drawModes_, drawChars_,
                             false, TextAlignType::MIDDLE);
}

// ****************************************************************************
void AtomSymbol::draw(MolDraw2D &molDrawer) const {
  std::string o_class = molDrawer.getActiveClass();
  std::string actClass = o_class;
  if (!actClass.empty()) {
    actClass += " ";
  }
  actClass += (boost::format("atom-%d") % atIdx_).str();
  molDrawer.setActiveClass(actClass);
  textDrawer_.setColour(colour_);
  textDrawer_.drawString(symbol_, cds_, orient_);
  molDrawer.setActiveClass(o_class);
  //  drawRects(molDrawer);
}

// ****************************************************************************
bool AtomSymbol::doesRectClash(const StringRect &rect, double padding) const {
  for (auto &alrect : rects_) {
    auto oldTrans = alrect->trans_;
    alrect->trans_ += cds_;
    bool dii = alrect->doesItIntersect(rect, padding);
    alrect->trans_ = oldTrans;
    if (dii) {
      return true;
    }
  }
  return false;
}

// ****************************************************************************
void AtomSymbol::adjustColons() {
  if (symbol_.empty()) {
    return;  // but probably it's always got something in it.
  }
  size_t colonPos = symbol_.find(':');
  if (colonPos == std::string::npos) {
    return;
  }
  // we need to allow for markup in the symbol, such as <lit>[CH2;X2:4]</lit>
  // and the easiest way to do that is to use the fact that atomLabelToPieces
  // strips it out.
  std::string tmpSym = symbol_;
  while (true) {
    size_t ltPos = tmpSym.find('<');
    if (ltPos == std::string::npos) {
      break;
    }
    size_t gtPos = tmpSym.find('>');
    tmpSym = tmpSym.substr(0, ltPos) + tmpSym.substr(gtPos + 1);
  }
  colonPos = tmpSym.find(':');
  if (colonPos == std::string::npos) {
    return;
  }
  CHECK_INVARIANT(colonPos <= rects_.size(), "bad rects_ size");
  double leftHeight = colonPos ? rects_[colonPos - 1]->height_ : 0;
  double rightHeight =
      colonPos < symbol_.size() - 1 ? rects_[colonPos + 1]->height_ : 0;
  rects_[colonPos]->height_ = std::min(leftHeight, rightHeight);
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

}  // namespace MolDraw2D_detail
}  // namespace RDKit
