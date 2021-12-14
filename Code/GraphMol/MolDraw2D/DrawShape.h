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

// A set of shapes used in 2D drawing.d  Not part of the public API.

#ifndef RDKIT_DRAWSHAPE_H
#define RDKIT_DRAWSHAPE_H

#include <vector>

#include <Geometry/point.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

namespace RDKit {

class MolDraw2D;

class DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  virtual ~DrawShape() = default;

 protected:
  DrawShape(const std::vector<Point2D> &points, int lineWidth = 2,
            bool scaleLineWidth = false,
            DrawColour lineColour = DrawColour(0, 0, 0), bool fill = false);

  void draw(MolDraw2D &drawer);
  virtual void myDraw(MolDraw2D &drawer) const = 0;
  virtual void findExtremes(double &xmin, double &xmax,
                            double &ymin, double &ymax) const;
  virtual void scale(const Point2D &scale);
  virtual void move(const Point2D &trans);

  std::vector<Point2D> points_;
  DrawColour lineColour_;
  int lineWidth_;
  bool fill_;
  bool scaleLineWidth_;

};

class DrawShapeArrow: protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeArrow() = default;

 protected:
  DrawShapeArrow(const std::vector<Point2D> &points, int lineWidth = 2,
                 bool scaleLineWidth = false,
                 DrawColour lineColour = DrawColour(0, 0, 0),
                 bool fill = false);
  void myDraw(MolDraw2D &drawer) const;

};

class DrawShapeEllipse: protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeEllipse() = default;

 protected:
  DrawShapeEllipse(const std::vector<Point2D> &points, int lineWidth = 2,
                   bool scaleLineWidth = false,
                   DrawColour lineColour = DrawColour(0, 0, 0),
                   bool fill = false);
  void myDraw(MolDraw2D &drawer) const;

};

class DrawShapePolyline: protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapePolyline() = default;

 protected:
  DrawShapePolyline(const std::vector<Point2D> &points, int lineWidth = 2,
                    bool scaleLineWidth = false,
                    DrawColour lineColour = DrawColour(0, 0, 0),
                    bool fill = false);
  void myDraw(MolDraw2D &drawer) const;

};


} // namespace RDKit

#endif  // RDKIT_DRAWSHAPE_H
