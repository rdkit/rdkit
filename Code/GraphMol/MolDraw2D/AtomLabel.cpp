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

#include <GraphMol/MolDraw2D/AtomLabel.h>
#include <GraphMol/MolDraw2D/DrawText.h>

namespace RDKit {

// ****************************************************************************
AtomLabel::AtomLabel(const std::string &label, int atIdx, OrientType orient,
                     const Point2D &cds, const DrawColour &colour,
                     DrawText &textDrawer)
    : label_(label),
      atIdx_(atIdx),
      orient_(orient),
      cds_(cds),
      colour_(colour),
      textDrawer_(textDrawer) {
  if (label.empty()) {
    return;
  }
  textDrawer_.getStringRects(label_, orient_, rects_, draw_modes_,
                             draw_chars_, false, TextAlignType::MIDDLE);
}

// ****************************************************************************
void AtomLabel::draw(MolDraw2D &molDrawer) const {
  textDrawer_.setColour(colour_);
  molDrawer.setActiveAtmIdx(atIdx_);
  textDrawer_.drawString(label_, cds_, orient_);
  molDrawer.setActiveAtmIdx();
  drawRects(molDrawer);
}

// ****************************************************************************
void AtomLabel::drawRects(MolDraw2D &molDrawer) const {
  for (auto &rect : rects_) {
    Point2D origTrans = rect->trans_;
    rect->trans_ += cds_;
    Point2D tl, tr, br, bl;
    rect->calcCorners(tl, tr, br, bl, 0.0);
    molDrawer.setColour(DrawColour(1.0, 0.0, 0.0));
    molDrawer.drawLine(tl, tr);
    molDrawer.setColour(DrawColour(0.0, 1.0, 0.0));
    molDrawer.drawLine(tr, br);
    molDrawer.setColour(DrawColour(0.0, 0.0, 1.0));
    molDrawer.drawLine(br, bl);
    molDrawer.setColour(DrawColour(0.0, 0.95, 0.95));
    molDrawer.drawLine(bl, tl);
    rect->trans_ = origTrans;
  }
}

// ****************************************************************************
void AtomLabel::scale(const Point2D &scaleFactor) {
  cds_.x *= scaleFactor.x;
  cds_.y *= scaleFactor.y;
  for (auto &rect : rects_) {
    rect->trans_.x *= scaleFactor.x;
    rect->trans_.y *= scaleFactor.y;
    rect->offset_.x *= scaleFactor.x;
    rect->offset_.y *= scaleFactor.y;
    rect->g_centre_.x *= scaleFactor.x;
    rect->g_centre_.y *= scaleFactor.y;
    rect->y_shift_ *= scaleFactor.y;
    rect->width_ *= scaleFactor.x;
    rect->height_ *= scaleFactor.y;
    rect->rect_corr_ *= scaleFactor.y;
  }
}

// ****************************************************************************
void AtomLabel::move(const Point2D &trans) {
  cds_ += trans;
}

} // namespace RDKit

