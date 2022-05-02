//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
#include "DrawTextFTQt.h"

#include <QPainter>
#include <QPainterPath>

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************
DrawTextFTQt::DrawTextFTQt(double max_fnt_sz, double min_fnt_sz,
                           const std::string &font_file, QPainter *qp)
    : DrawTextFT(max_fnt_sz, min_fnt_sz, font_file), d_qp(qp) {}

DrawTextFTQt::~DrawTextFTQt() = default;

// ****************************************************************************
double DrawTextFTQt::extractOutline() {
  dp_qpp.reset(new QPainterPath());
  double adv = DrawTextFT::extractOutline();
  const auto &col = colour();
  QColor this_col(int(255.0 * col.r), int(255.0 * col.g), int(255.0 * col.b),
                  int(255.0 * col.a));

  QPen pen(this_col);
  pen.setJoinStyle(Qt::RoundJoin);
  d_qp->setPen(pen);

  QBrush brush(this_col);
  brush.setStyle(Qt::SolidPattern);
  d_qp->setBrush(brush);

  d_qp->fillPath(*dp_qpp, brush);
  return adv;
}

// ****************************************************************************
int DrawTextFTQt::MoveToFunctionImpl(const FT_Vector *to) {
  PRECONDITION(dp_qpp, "no path");
  double dx, dy;
  fontPosToDrawPos(to->x, to->y, dx, dy);
  dp_qpp->moveTo(dx, dy);
  return 0;
}

// ****************************************************************************
int DrawTextFTQt::LineToFunctionImpl(const FT_Vector *to) {
  PRECONDITION(dp_qpp, "no path");
  double dx, dy;
  fontPosToDrawPos(to->x, to->y, dx, dy);
  dp_qpp->lineTo(dx, dy);

  return 0;
}

// ****************************************************************************
int DrawTextFTQt::ConicToFunctionImpl(const FT_Vector *control,
                                      const FT_Vector *to) {
  PRECONDITION(dp_qpp, "no path");
  // Qt doesn't have a quadratic bezier function, we need to promote
  // it to a cubic using formula from
  // https://lists.cairographics.org/archives/cairo/2010-April/019691.html
  // and correcting the typo on the last line.

  double x0, y0;
  auto currp = dp_qpp->currentPosition();
  x0 = currp.x();
  y0 = currp.y();

  double x1, y1;
  fontPosToDrawPos(control->x, control->y, x1, y1);

  double x2, y2;
  fontPosToDrawPos(to->x, to->y, x2, y2);

  dp_qpp->cubicTo((2.0 / 3.0) * x1 + (1.0 / 3.0) * x0,
                  (2.0 / 3.0) * y1 + (1.0 / 3.0) * y0,
                  (2.0 / 3.0) * x1 + (1.0 / 3.0) * x2,
                  (2.0 / 3.0) * y1 + (1.0 / 3.0) * y2, x2, y2);

  return 0;
}

// ****************************************************************************
int DrawTextFTQt::CubicToFunctionImpl(const FT_Vector *controlOne,
                                      const FT_Vector *controlTwo,
                                      const FT_Vector *to) {
  PRECONDITION(dp_qpp, "no path");
  double controlOneX, controlOneY;
  fontPosToDrawPos(controlOne->x, controlOne->y, controlOneX, controlOneY);
  double controlTwoX, controlTwoY;
  fontPosToDrawPos(controlTwo->x, controlTwo->y, controlTwoX, controlTwoY);

  double dx, dy;
  fontPosToDrawPos(to->x, to->y, dx, dy);

  dp_qpp->cubicTo(controlOneX, controlOneY, controlTwoX, controlTwoY, dx, dy);
  return 0;
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit