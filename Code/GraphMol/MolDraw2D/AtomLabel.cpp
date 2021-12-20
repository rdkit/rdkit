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
  std::cout << label_ << " at " << cds_;
  for (auto &rect : rects_) {
    std::cout << " : " << rect->trans_ << " size " << rect->width_
              << " by " << rect->height_;    Point2D tl, tr, br, bl;
    rect->calcCorners(tl, tr, br, bl, 0.0);
    std::cout << "  rect width : " << tr.x - tl.x << " or " << br.x - bl.x << std::endl;
  }
  std::cout << std::endl;
}

// ****************************************************************************
void AtomLabel::draw(MolDraw2D &molDrawer) const {
  textDrawer_.setColour(colour_);
  molDrawer.setActiveAtmIdx(atIdx_);
  textDrawer_.drawString(label_, cds_, orient_);
  molDrawer.setActiveAtmIdx();
//  drawRects(molDrawer);
}

// ****************************************************************************
void AtomLabel::drawRects(MolDraw2D &molDrawer) const {
  for (auto &rect : rects_) {
    Point2D tl, tr, br, bl;
    rect->calcCorners(tl, tr, br, bl, 0.0);
    tl += cds_;
    tr += cds_;
    br += cds_;
    bl += cds_;
    molDrawer.setColour(DrawColour(1.0, 0.0, 0.0));
    molDrawer.drawLine(tl, tr);
    molDrawer.setColour(DrawColour(0.0, 1.0, 0.0));
    molDrawer.drawLine(tr, br);
    molDrawer.setColour(DrawColour(0.0, 0.0, 1.0));
    molDrawer.drawLine(br, bl);
    molDrawer.setColour(DrawColour(0.0, 0.95, 0.95));
    molDrawer.drawLine(bl, tl);
  }
}

// ****************************************************************************
void AtomLabel::scale(const Point2D &scale_factor) {
  cds_.x *= scale_factor.x;
  cds_.y *= scale_factor.y;

  for (auto &rect : rects_) {
    rect->trans_.x *= scale_factor.x;
    rect->trans_.y *= scale_factor.y;
    rect->offset_.x *= scale_factor.x;
    rect->offset_.y *= scale_factor.y;
    rect->g_centre_.x *= scale_factor.x;
    rect->g_centre_.y *= scale_factor.y;
    rect->y_shift_ *= scale_factor.y;
    // width and height of the character use the font scale not the global
    // drawing scale, which is what scale_factor is.  fontScale() needs to
    // be set before this scale() is called.
    rect->width_ *= textDrawer_.fontScale();
    rect->height_ *= textDrawer_.fontScale();
    rect->rect_corr_ *= scale_factor.y;
  }
}

// ****************************************************************************
void AtomLabel::move(const Point2D &trans) {
  cds_ += trans;
}

} // namespace RDKit

