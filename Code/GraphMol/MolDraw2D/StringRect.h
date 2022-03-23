//
//  Copyright (C) 2021-2022 David Cosgrove and other RDKit contributors
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
// This is used by DrawText classes.  It is not intended for general use.

#ifndef RDKIT_STRINGRECT_H
#define RDKIT_STRINGRECT_H

#include <Geometry/point.h>

namespace RDKit {
namespace MolDraw2D_detail {

// for holding dimensions of the rectangle round a string.
struct StringRect {
  Point2D trans_;     // Where to draw char relative to other chars in string
  Point2D offset_;    // offset for draw coords so char is centred correctly
  Point2D g_centre_;  // glyph centre relative to the origin of the char.
  double y_shift_;  // shift the whole thing in y by this. For multi-line text.
  double width_, height_;  // of the glyph itself, not the character cell
  double rect_corr_;  // because if we move a char one way, we need to move the
                           // rectangle the other.
  int clash_score_;   // rough measure of how badly it clashed with other things
                           // lower is better, 0 is no clash.

  StringRect()
      : trans_(0.0, 0.0),
        offset_(0.0, 0.0),
        g_centre_(offset_),
        y_shift_(0.0),
        width_(0.0),
        height_(0.0),
        rect_corr_(0.0),
        clash_score_(0) {}
  StringRect(const Point2D &offset, const Point2D &g_centre, double w, double h)
      : trans_(0.0, 0.0),
        offset_(offset),
        g_centre_(g_centre),
        y_shift_(0.0),
        width_(w),
        height_(h),
        rect_corr_(0.0),
        clash_score_(0) {}
  // tl is top, left; br is bottom, right of the glyph, relative to the
  // centre. Padding in draw coords.
  void calcCorners(Point2D &tl, Point2D &tr, Point2D &br, Point2D &bl,
                   double padding) const {
    double wb2 = padding + width_ / 2.0;
    double hb2 = padding + height_ / 2.0;
    Point2D c;
    calcCentre(c);
    tl = Point2D(c.x - wb2, c.y - hb2);
    tr = Point2D(c.x + wb2, c.y - hb2);
    br = Point2D(c.x + wb2, c.y + hb2);
    bl = Point2D(c.x - wb2, c.y + hb2);
  }
  void calcCentre(Point2D &c) const {
    c = trans_ + g_centre_ - offset_;
    c.y -= y_shift_;
  }
  bool isPointInside(const Point2D &pt, double padding = 0.0) const {
    Point2D tl, tr, br, bl;
    calcCorners(tl, tr, br, bl, padding);
    // is +ve y up or down?
    if (tl.y < bl.y) {
      std::swap(tl, bl);
      std::swap(tr, br);
    }
    return pt.x >= tl.x && pt.x <= br.x && pt.y >= br.y && pt.y <= tl.y;
  }
  bool doesItIntersect(const StringRect &other, double padding = 0.0) const {
    Point2D ttl, ttr, tbr, tbl;
    calcCorners(ttl, ttr, tbr, tbl, padding);
    // is +ve y up or down?
    if (ttl.y < tbl.y) {
      std::swap(ttl, tbl);
      std::swap(ttr, tbr);
    }
    Point2D otl, otr, obr, obl;
    other.calcCorners(otl, otr, obr, obl, padding);
    if (otl.y < obl.y) {
      std::swap(otl, obl);
      std::swap(otr, obr);
    }
    // This could be done with isPointInside, but that would recalculate
    // the corners each time.
    if ((otl.x >= ttl.x && otl.x <= ttr.x && otl.y >= tbl.y &&
         otl.y <= ttl.y) ||
        (otr.x >= ttl.x && otr.x <= ttr.x && otr.y >= tbl.y &&
         otr.y <= ttl.y) ||
        (obr.x >= ttl.x && obr.x <= ttr.x && obr.y >= tbl.y &&
         obr.y <= ttl.y) ||
        (obl.x >= ttl.x && obl.x <= ttr.x && obl.y >= tbl.y &&
         obl.y <= ttl.y)) {
      return true;
    }
    if ((ttl.x >= otl.x && ttl.x <= otr.x && ttl.y >= obl.y &&
         ttl.y <= otl.y) ||
        (ttr.x >= otl.x && ttr.x <= otr.x && ttr.y >= obl.y &&
         ttr.y <= otl.y) ||
        (tbr.x >= otl.x && tbr.x <= otr.x && tbr.y >= obl.y &&
         tbr.y <= otl.y) ||
        (tbl.x >= otl.x && tbl.x <= otr.x && tbl.y >= obl.y &&
         tbl.y <= otl.y)) {
      return true;
    }
    return false;
  }
};

}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif  // RDKIT_STRINGRECT_H
