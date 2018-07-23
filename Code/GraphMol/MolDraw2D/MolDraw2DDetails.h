//
//  Copyright (C) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RDKITMOLDRAW2DDETAILS_H
#define RDKITMOLDRAW2DDETAILS_H

#include <vector>

#include <Geometry/point.h>
#include <GraphMol/RDKitBase.h>

#include <boost/tuple/tuple.hpp>

// ****************************************************************************
using RDGeom::Point2D;

namespace RDKit {
namespace MolDraw2D_detail {
// data taken from the helvetica font info in
// $RDBASE/rdkit/sping/PDF/pdfmetrics.py
const int char_widths[] = {
    0,   0,    0,   0,   0,   0,   0,    0,    0,    0,   0,    0,   0,    0,
    0,   0,    0,   0,   0,   0,   0,    0,    0,    0,   0,    0,   0,    0,
    0,   0,    0,   0,   278, 278, 355,  556,  556,  889, 667,  222, 333,  333,
    389, 584,  278, 333, 278, 278, 556,  556,  556,  556, 556,  556, 556,  556,
    556, 556,  278, 278, 584, 584, 584,  556,  1015, 667, 667,  722, 722,  667,
    611, 778,  722, 278, 500, 667, 556,  833,  722,  778, 667,  778, 722,  667,
    611, 722,  667, 944, 667, 667, 611,  278,  278,  278, 469,  556, 222,  556,
    556, 500,  556, 556, 278, 556, 556,  222,  222,  500, 222,  833, 556,  556,
    556, 556,  333, 500, 278, 556, 500,  722,  500,  500, 500,  334, 260,  334,
    584, 0,    0,   0,   0,   0,   0,    0,    0,    0,   0,    0,   0,    0,
    0,   0,    0,   0,   0,   0,   0,    0,    0,    0,   0,    0,   0,    0,
    0,   0,    0,   0,   0,   0,   0,    333,  556,  556, 167,  556, 556,  556,
    556, 191,  333, 556, 333, 333, 500,  500,  0,    556, 556,  556, 278,  0,
    537, 350,  222, 333, 333, 556, 1000, 1000, 0,    611, 0,    333, 333,  333,
    333, 333,  333, 333, 333, 0,   333,  333,  0,    333, 333,  333, 1000, 0,
    0,   0,    0,   0,   0,   0,   0,    0,    0,    0,   0,    0,   0,    0,
    0,   1000, 0,   370, 0,   0,   0,    0,    556,  778, 1000, 365, 0,    0,
    0,   0,    0,   889, 0,   0,   0,    278,  0,    0,   222,  611, 944,  611,
    0,   0,    834};

RDKIT_MOLDRAW2D_EXPORT void arcPoints(const Point2D &cds1, const Point2D &cds2,
               std::vector<Point2D> &res, float startAng = 0,
               float extent = 360);
}
}

#endif
