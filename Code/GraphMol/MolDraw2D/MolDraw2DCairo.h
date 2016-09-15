//
//  Copyright (C) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// derived from Dave Cosgrove's MolDraw2D
//
// This is a concrete class derived from MolDraw2D that uses RDKit to draw a
// molecule into a cairo drawing context

#ifndef MOLDRAW2DCAIRO_H
#define MOLDRAW2DCAIRO_H

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <cairo.h>

// ****************************************************************************

namespace RDKit {

class MolDraw2DCairo : public MolDraw2D {
 public:
  // does not take ownership of the drawing context
  MolDraw2DCairo(int width, int height, cairo_t *cr, int panelWidth = -1,
                 int panelHeight = -1)
      : MolDraw2D(width, height, panelWidth, panelHeight), dp_cr(cr) {
    cairo_reference(dp_cr);
    initDrawing();
  };
  MolDraw2DCairo(int width, int height, int panelWidth = -1,
                 int panelHeight = -1)
      : MolDraw2D(width, height, panelWidth, panelHeight) {
    cairo_surface_t *surf =
        cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
    dp_cr = cairo_create(surf);
    cairo_surface_destroy(surf);  // dp_cr has a reference to this now;
    initDrawing();
  };
  ~MolDraw2DCairo() {
    if (dp_cr) {
      if (cairo_get_reference_count(dp_cr) > 0) {
        cairo_destroy(dp_cr);
      }
      dp_cr = NULL;
    }
  }

  // set font size in molecule coordinate units. That's probably Angstrom for
  // RDKit. It will turned into drawing units using scale_, which might be
  // changed as a result, to make sure things still appear in the window.
  void setFontSize(double new_size);
  void setColour(const DrawColour &col);

  // not sure if this goes here or if we should do a dtor since initDrawing() is
  // called in the ctor,
  // but we'll start here
  void finishDrawing();

  void drawLine(const Point2D &cds1, const Point2D &cds2);
  void drawChar(char c, const Point2D &cds);
  // void drawString( const std::string &str, const Point2D &cds );
  void drawPolygon(const std::vector<Point2D> &cds);
  void clearDrawing();

  void drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                    const DrawColour &col1, const DrawColour &col2,
                    unsigned int nSegments = 16, double vertOffset = 0.05);

  // using the current scale, work out the size of the label in molecule
  // coordinates
  void getStringSize(const std::string &label, double &label_width,
                     double &label_height) const;

  // returns the PNG data in a string
  std::string getDrawingText() const;
  // writes the PNG data to a file
  void writeDrawingText(const std::string &fName) const;

 private:
  cairo_t *dp_cr;

  void initDrawing();
};
}
#endif  // MOLDRAW2DCAIRO_H
