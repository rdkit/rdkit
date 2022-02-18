//
//  Copyright (C) 2022 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// Original author: David Cosgrove (CozChemIx).
//

#include <GraphMol/MolDraw2D/DrawTextNotFT.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

namespace RDKit {
namespace MolDraw2D_detail {
// ****************************************************************************
DrawTextNotFT::DrawTextNotFT(double max_fnt_sz, double min_fnt_sz)
    : DrawText(max_fnt_sz, min_fnt_sz) {}

// ****************************************************************************
DrawTextNotFT::~DrawTextNotFT() {}

// ****************************************************************************
void DrawTextNotFT::alignString(
    TextAlignType talign, const std::vector<TextDrawType> &draw_modes,
    std::vector<std::shared_ptr<StringRect>> &rects) const {
  // std::string comes in with rects aligned with first char with its
  // left hand and bottom edges at 0 on y and x respectively.
  // Adjust relative to that so that the relative alignment point is at
  // (0,0).
  if (talign == TextAlignType::MIDDLE) {
    size_t num_norm = count(draw_modes.begin(), draw_modes.end(),
                            TextDrawType::TextDrawNormal);
    if (num_norm == 1) {
      talign = TextAlignType::START;
    }
  }

  Point2D align_trans, align_offset;
  if (talign == TextAlignType::START || talign == TextAlignType::END) {
    size_t align_char = 0;
    for (size_t i = 0; i < rects.size(); ++i) {
      if (draw_modes[i] == TextDrawType::TextDrawNormal) {
        align_char = i;
        if (talign == TextAlignType::START) {
          break;
        }
      }
    }
    align_trans = rects[align_char]->trans_;
    align_offset = rects[align_char]->offset_;
  } else {
    // centre on the middle of the Normal text.  The super- or subscripts
    // should be at the ends.
    double x_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::lowest();
    align_offset.x = align_offset.y = 0.0;
    int num_norm = 0;
    for (size_t i = 0; i < rects.size(); ++i) {
      if (draw_modes[i] == TextDrawType::TextDrawNormal) {
        Point2D tl, tr, br, bl;
        rects[i]->calcCorners(tl, tr, br, bl, 0.0);
        // sometimes the rect is in a coordinate frame where +ve y is down,
        // sometimes it's up.  For these purposes, we don't care so long as
        // the y_max is larger than the y_min.  We probably don't need to do
        // all the tests for x_min and x_max;
        x_min = std::min({bl.x, tr.x, x_min});
        x_max = std::max({bl.x, tr.x, x_max});
        align_offset += rects[i]->offset_;
        ++num_norm;
      }
    }
    align_trans.x = (x_max - x_min) / 2.0;
    align_trans.y = 0.0;
    align_offset /= num_norm;
  }

  for (auto r : rects) {
    r->trans_ -= align_trans;
    r->offset_ = align_offset;
  }
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit
