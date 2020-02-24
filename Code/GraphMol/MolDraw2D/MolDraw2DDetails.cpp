//
//  Copyright (C) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ****************************************************************************

namespace RDKit {
namespace MolDraw2D_detail {
// implementation from $RDBASE/rdkit/sping/pid.py
void arcPoints(const Point2D &cds1, const Point2D &cds2,
               std::vector<Point2D> &res, float startAng, float extent) {
  // Note: this implementation is simple and not particularly efficient.
  float xScale = (cds2.x - cds1.x) / 2.0;
  float yScale = (cds2.y - cds1.y) / 2.0;
  if (xScale < 0) {
    xScale *= -1;
  }
  if (yScale < 0) {
    yScale *= -1;
  }

  float x = std::min(cds1.x, cds2.x) + xScale;
  float y = std::min(cds1.y, cds2.y) + yScale;

  int steps = std::max(static_cast<int>(extent * 2), 5);
  float step = M_PI * extent / (180 * steps);
  float angle = M_PI * startAng / 180;
  for (int i = 0; i <= steps; ++i) {
    Point2D point(x + xScale * cos(angle), y - yScale * sin(angle));
    res.emplace_back(point);
    angle += step;
  }
}
}
}
