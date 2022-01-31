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

// A set of shapes used in 2D drawing.d  Not part of the public API.

#ifndef RDKIT_DRAWSHAPE_H
#define RDKIT_DRAWSHAPE_H

#include <vector>

#include <Geometry/point.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

namespace RDKit {

class MolDraw2D;
const DashPattern noDash;
const DashPattern dots{2.0, 6.0};
const DashPattern dashes{6, 0, 6.0};
const DashPattern shortDashes{2.0, 2.0};

namespace MolDraw2D_detail {

struct StringRect;

class DrawShape {
  friend class RDKit::MolDraw2D;
  friend class DrawMol;
  friend class DrawMolMCH;

 public:
  virtual ~DrawShape() = default;

 protected:
  DrawShape(const std::vector<Point2D> &points, int lineWidth = 2,
            bool scaleLineWidth = false,
            DrawColour lineColour = DrawColour(0, 0, 0), bool fill = false,
            int atom1 = -1, int atom2 = -1, int bond = -1);
  DrawShape(const DrawShape &) = delete;
  DrawShape(DrawShape &&) = delete;
  DrawShape &operator=(const DrawShape &) = delete;
  DrawShape &operator=(DrawShape &&) = delete;

  void draw(MolDraw2D &drawer);
  virtual void myDraw(MolDraw2D &drawer) const = 0;
  virtual void findExtremes(double &xmin, double &xmax, double &ymin,
                            double &ymax) const;
  virtual void scale(const Point2D &scale_factor);
  virtual void move(const Point2D &trans);
  virtual bool doesRectClash(const StringRect &rect, double padding) const;

  std::vector<Point2D> points_;
  int lineWidth_;
  bool scaleLineWidth_;
  DrawColour lineColour_;
  bool fill_;
  int atom1_, atom2_, bond_;
};

class DrawShapeArrow : protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeArrow() = default;

 protected:
  DrawShapeArrow(const std::vector<Point2D> &points, int lineWidth = 2,
                 bool scaleLineWidth = false,
                 DrawColour lineColour = DrawColour(0, 0, 0), bool fill = false,
                 int atom1 = -1, int atom2 = -1, int bond = -1,
                 double frac = 0.2, double angle = M_PI / 6);
  DrawShapeArrow(const DrawShapeArrow &) = delete;
  DrawShapeArrow(DrawShapeArrow &&) = delete;
  DrawShapeArrow &operator=(const DrawShapeArrow &) = delete;
  DrawShapeArrow &operator=(DrawShapeArrow &&) = delete;
  void myDraw(MolDraw2D &drawer) const override;
  bool doesRectClash(const StringRect &rect, double padding) const override;

  double frac_;
  double angle_;
};

class DrawShapeEllipse : protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;
  friend class DrawMolMCH;

 public:
  ~DrawShapeEllipse() = default;

 protected:
  // points are the 2 foci of the ellipse
  DrawShapeEllipse(const std::vector<Point2D> &points, int lineWidth = 2,
                   bool scaleLineWidth = false,
                   DrawColour lineColour = DrawColour(0, 0, 0),
                   bool fill = false, int atom1 = -1);
  DrawShapeEllipse(const DrawShapeEllipse &) = delete;
  DrawShapeEllipse(DrawShapeEllipse &&) = delete;
  DrawShapeEllipse &operator=(const DrawShapeEllipse &) = delete;
  DrawShapeEllipse &operator=(DrawShapeEllipse &&) = delete;
  void myDraw(MolDraw2D &drawer) const override;
  void findExtremes(double &xmin, double &xmax, double &ymin,
                    double &ymax) const override;
  bool doesRectClash(const StringRect &rect, double padding) const override;
};

class DrawShapeSimpleLine : protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;
  friend class DrawMolMCH;

 public:
  ~DrawShapeSimpleLine() = default;

 protected:
  DrawShapeSimpleLine(const std::vector<Point2D> &points, int lineWidth = 2,
                      bool scaleLineWidth = false,
                      DrawColour lineColour = DrawColour(0, 0, 0),
                      int atom1 = -1, int atom2 = -1, int bond = -1,
                      DashPattern dashPattern = noDash);
  DrawShapeSimpleLine(const DrawShapeSimpleLine &) = delete;
  DrawShapeSimpleLine(DrawShapeSimpleLine &&) = delete;
  DrawShapeSimpleLine &operator=(const DrawShapeSimpleLine &) = delete;
  DrawShapeSimpleLine &operator=(DrawShapeSimpleLine &&) = delete;
  void myDraw(MolDraw2D &drawer) const override;
  bool doesRectClash(const StringRect &rect, double padding) const override;

  DashPattern dashPattern_;
};

class DrawShapePolyLine : protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;
  friend class DrawMolMCH;

 public:
  ~DrawShapePolyLine() = default;

 protected:
  DrawShapePolyLine(const std::vector<Point2D> &points, int lineWidth = 2,
                    bool scaleLineWidth = false,
                    DrawColour lineColour = DrawColour(0, 0, 0),
                    bool fill = false, int atom1 = -1, int atom2 = -1,
                    int bond = -1, DashPattern dashPattern = noDash);
  DrawShapePolyLine(const DrawShapePolyLine &) = delete;
  DrawShapePolyLine(DrawShapePolyLine &&) = delete;
  DrawShapePolyLine &operator=(const DrawShapePolyLine &) = delete;
  DrawShapePolyLine &operator=(DrawShapePolyLine &&) = delete;
  void myDraw(MolDraw2D &drawer) const override;
  bool doesRectClash(const StringRect &rect, double padding) const override;

  DashPattern dashPattern_;
};

class DrawShapeSolidWedge : protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeSolidWedge() = default;

 protected:
  DrawShapeSolidWedge(const std::vector<Point2D> points, const DrawColour &col1,
                      const DrawColour &col2, bool splitBonds, int atom1 = -1,
                      int atom2 = -1, int bond = -1);
  DrawShapeSolidWedge(const DrawShapeSolidWedge &) = delete;
  DrawShapeSolidWedge(DrawShapeSolidWedge &&) = delete;
  DrawShapeSolidWedge &operator=(const DrawShapeSolidWedge &) = delete;
  DrawShapeSolidWedge &operator=(DrawShapeSolidWedge &&) = delete;
  void buildTriangles();
  void myDraw(MolDraw2D &drawer) const override;
  bool doesRectClash(const StringRect &rect, double padding) const override;

  DrawColour col2_;
  bool splitBonds_;
};

class DrawShapeDashedWedge : protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeDashedWedge() = default;

 protected:
  DrawShapeDashedWedge(const std::vector<Point2D> points,
                       const DrawColour &col1, const DrawColour &col2,
                       int atom1 = -1, int atom2 = -1, int bond = -1);
  DrawShapeDashedWedge(const DrawShapeDashedWedge &) = delete;
  DrawShapeDashedWedge(DrawShapeDashedWedge &&) = delete;
  DrawShapeDashedWedge &operator=(const DrawShapeDashedWedge &) = delete;
  DrawShapeDashedWedge &operator=(DrawShapeDashedWedge &&) = delete;
  void buildLines();
  void myDraw(MolDraw2D &drawer) const override;
  void scale(const Point2D &scale_factor) override;
  void move(const Point2D &trans) override;
  bool doesRectClash(const StringRect &rect, double padding) const override;

  DrawColour col2_;
  std::vector<DrawColour> lineColours_;
  // for when we re-create the lines when it gets too wide, this is
  // the initial points[0] from the c'tor.
  Point2D at1Cds_;
};

class DrawShapeWavyLine : protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeWavyLine() = default;

 protected:
  DrawShapeWavyLine(const std::vector<Point2D> points, int lineWidth = 2,
                    bool scaleLineWidth = false,
                    const DrawColour &col1 = DrawColour(0, 0, 0),
                    const DrawColour &col2 = DrawColour(0, 0, 0),
                    double offset = 0.05, int atom1 = -1, int atom2 = -1,
                    int bond = -1);
  DrawShapeWavyLine(const DrawShapeWavyLine &) = delete;
  DrawShapeWavyLine(DrawShapeWavyLine &&) = delete;
  DrawShapeWavyLine &operator=(const DrawShapeWavyLine &) = delete;
  DrawShapeWavyLine &operator=(DrawShapeWavyLine &&) = delete;
  void myDraw(MolDraw2D &drawer) const override;
  void scale(const Point2D &scaleFactor) override;
  bool doesRectClash(const StringRect &rect, double padding) const override;

  DrawColour col2_;
  double offset_;
};

class DrawShapeArc : protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;
  friend class DrawMolMCH;

 public:
  ~DrawShapeArc() = default;

 protected:
  // draw the arc of an ellipse between ang1 and ang2.  Note that 0 is
  // at 3 o-clock and 90 at 12 o'clock as you'd expect from your maths.
  // ang2 must be > ang1 - it won't draw backwards. Angles in degrees,
  // between 0 and 360.0.
  // Points should be size 2 - the first entry is the centre, the second
  // gives the x and y radii of the ellipse.
  DrawShapeArc(const std::vector<Point2D> points, double ang1, double ang2,
               int lineWidth = 2, bool scaleLineWidth = false,
               const DrawColour &col1 = DrawColour(0, 0, 0), bool fill = false,
               int atom1 = -1);
  DrawShapeArc(const DrawShapeArc &) = delete;
  DrawShapeArc(DrawShapeArc &&) = delete;
  DrawShapeArc &operator=(const DrawShapeArc &) = delete;
  DrawShapeArc &operator=(DrawShapeArc &&) = delete;

  void myDraw(MolDraw2D &drawer) const override;
  void findExtremes(double &xmin, double &xmax, double &ymin,
                    double &ymax) const override;
  void move(const Point2D &trans) override;
  bool doesRectClash(const StringRect &rect, double padding) const override;

  double ang1_, ang2_;
};

std::vector<Point2D> calcScaledWedgePoints(const Point2D &point,
                                           const Point2D &end1,
                                           const Point2D &end2, double widthsq);
}  // namespace MolDraw2D_detail
}  // namespace RDKit

#endif  // RDKIT_DRAWSHAPE_H
