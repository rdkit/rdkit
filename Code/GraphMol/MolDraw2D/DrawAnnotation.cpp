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

#include <GraphMol/MolDraw2D/DrawAnnotation.h>
#include <GraphMol/MolDraw2D/DrawText.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************
DrawAnnotation::DrawAnnotation(const std::string &note,
                               const TextAlignType &align,
                               const std::string &cls, double relFontScale,
                               const Point2D &pos, const DrawColour &colour,
                               DrawText &textDrawer)
    : text_(note),
      align_(align),
      class_(cls),
      textDrawer_(textDrawer),
      pos_(pos),
      colour_(colour) {
  fontScale_ = relFontScale * textDrawer_.fontScale();
  extractRects();
}

// ****************************************************************************
void DrawAnnotation::findExtremes(double &xmin, double &xmax, double &ymin,
                                  double &ymax, double padding) const {
  if (text_.empty()) {
    return;
  }
  Point2D tl, tr, br, bl, otrans;
  for (auto r : rects_) {
    otrans = r->trans_;
    r->trans_ += pos_;
    r->calcCorners(tl, tr, br, bl, padding);
    // sometimes the rect is in a coordinate frame where +ve y is down,
    // sometimes it's up.  For these purposes, we don't care so long as
    // the ymax is larger than the ymin.  We probably don't need to do
    // all the tests for xmin and xmax;
    xmin = std::min({tr.x, bl.x, xmin});
    ymin = std::min({tr.y, bl.y, ymin});
    xmax = std::max({tr.x, bl.x, xmax});
    ymax = std::max({tr.y, bl.y, ymax});
    r->trans_ = otrans;
  }
}

// ****************************************************************************
void DrawAnnotation::getDimensions(double &width, double &height) const {
  double xMin, yMin, xMax, yMax;
  xMin = yMin = std::numeric_limits<double>::max();
  xMax = yMax = std::numeric_limits<double>::lowest();
  findExtremes(xMin, xMax, yMin, yMax);
  width = xMax - xMin;
  height = yMax - yMin;
}

// ****************************************************************************
void DrawAnnotation::extractRects() {
  // We don't need these for notes, which are always on 1 line and plain
  // text.
  std::vector<TextDrawType> drawModes;
  std::vector<char> drawChars;
  double ofs = textDrawer_.fontScale();
  textDrawer_.setFontScale(fontScale_, true);
  fontScale_ = textDrawer_.fontScale();
  textDrawer_.getStringRects(text_, OrientType::C, rects_, drawModes, drawChars,
                             true, align_);
  textDrawer_.setFontScale(ofs, true);
}

// ****************************************************************************
void DrawAnnotation::draw(MolDraw2D &molDrawer) const {
  std::string o_class = molDrawer.getActiveClass();
  std::string actClass = o_class;
  if (!actClass.empty()) {
    actClass += " ";
  }
  actClass += class_;
  molDrawer.setActiveClass(actClass);
  textDrawer_.setColour(colour_);
  double ofs = textDrawer_.fontScale();
  textDrawer_.setFontScale(fontScale_, true);
  textDrawer_.drawString(text_, pos_, align_);
  textDrawer_.setFontScale(ofs, true);
  molDrawer.setActiveClass(o_class);
  // drawRects(molDrawer);
}

// ****************************************************************************
void DrawAnnotation::drawRects(MolDraw2D &molDrawer) const {
  auto olw = molDrawer.lineWidth();
  molDrawer.setLineWidth(1);
  Point2D tl, tr, br, bl, origTrans;
  for (auto &rect : rects_) {
    origTrans = rect->trans_;
    rect->trans_ += pos_;
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
  molDrawer.setLineWidth(olw);
}

// ****************************************************************************
void DrawAnnotation::scale(const Point2D &scaleFactor) {
  pos_.x *= scaleFactor.x;
  pos_.y *= scaleFactor.y;
  // arbitrarily choose x scale for fonts.  It is highly unlikely that the
  // x and y are different, in any case.
  fontScale_ *= scaleFactor.x;
  // rebuild the rectangles, because the fontScale may be different,
  // and the widths etc might not scale by the same amount.
  rects_.clear();
  extractRects();
}

// ****************************************************************************
void DrawAnnotation::move(const Point2D &trans) {
  pos_ += trans;
  // the internals of the rects_ are all relative to pos_, so no need to
  // do anything to them.
}

// ****************************************************************************
bool DrawAnnotation::doesRectClash(const StringRect &rect,
                                   double padding) const {
  for (auto &alrect : rects_) {
    auto oldTrans = alrect->trans_;
    alrect->trans_ += pos_;
    bool dii = alrect->doesItIntersect(rect, padding);
    alrect->trans_ = oldTrans;
    if (dii) {
      return true;
    }
  }
  return false;
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit
