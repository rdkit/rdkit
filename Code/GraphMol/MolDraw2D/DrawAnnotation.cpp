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

#include <GraphMol/MolDraw2D/DrawAnnotation.h>
#include <GraphMol/MolDraw2D/DrawText.h>

namespace RDKit {

// ****************************************************************************
DrawAnnotation::DrawAnnotation(const std::string &note,
                               const TextAlignType &align,
                               const std::string &cls, double relFontScale,
                               const Point2D &pos, const DrawColour &colour,
                               DrawText &textDrawer)
    : text_(note),
      align_(align),
      class_(cls),
      relFontScale_(relFontScale),
      textDrawer_(textDrawer),
      pos_(pos),
      colour_(colour) {
  std::vector<TextDrawType> drawModes;
  std::vector<char> drawChars;
  double ofs = textDrawer_.fontScale();
  textDrawer_.setFontScale(relFontScale_ * textDrawer_.fontScale(), true);
  textDrawer_.getStringRects(text_, OrientType::C, rects_, drawModes, drawChars,
                             true, align_);
  textDrawer_.setFontScale(ofs, true);
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
    xmin = std::min(bl.x, xmin);
    xmin = std::min(tr.x, xmin);
    ymin = std::min(bl.y, ymin);
    ymin = std::min(tr.y, ymin);
    xmax = std::max(bl.x, xmax);
    xmax = std::max(tr.x, xmax);
    ymax = std::max(bl.y, ymax);
    ymax = std::max(tr.y, ymax);
    r->trans_ = otrans;
  }
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
    textDrawer_.setFontScale(relFontScale_ * ofs, true);
    textDrawer_.drawString(text_, pos_, align_);
    textDrawer_.setFontScale(ofs, true);
    molDrawer.setActiveClass(o_class);
//    drawRects(molDrawer);
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
void DrawAnnotation::scale(const Point2D &scaleFactor,
                           const Point2D &fontScaleFactor) {
  pos_.x *= scaleFactor.x;
  pos_.y *= scaleFactor.y;
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
void DrawAnnotation::move(const Point2D &trans) {
  pos_ += trans;
  // the internals of the rects_ are all relative to pos_, so no need to
  // do anything to them.
}

// ****************************************************************************
bool DrawAnnotation::doesRectClash(const StringRect &rect,
                                   double padding) const {
  for (auto &alrect : rects_) {
    StringRect r{*alrect};
    r.trans_ += pos_;
    if (r.doesItIntersect(rect, padding)) {
      return true;
    }
  }
  return false;
}

} // namespace RDKit
