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

#include <RDGeneral/export.h>
#ifndef MOLDRAW2DCAIRO_H
#define MOLDRAW2DCAIRO_H

#include <cairo.h>

#include <GraphMol/MolDraw2D/MolDraw2D.h>

// ****************************************************************************

namespace RDKit {

class RDKIT_MOLDRAW2D_EXPORT MolDraw2DCairo : public MolDraw2D {
 public:
  // does not take ownership of the drawing context
  MolDraw2DCairo(int width, int height, cairo_t *cr, int panelWidth = -1,
                 int panelHeight = -1, bool noFreetype = false)
      : MolDraw2D(width, height, panelWidth, panelHeight), dp_cr(cr) {
    cairo_reference(dp_cr);
    df_noFreetype = noFreetype;
    initTextDrawer(noFreetype);
  }
  MolDraw2DCairo(int width, int height, int panelWidth = -1,
                 int panelHeight = -1, bool noFreetype = false)
      : MolDraw2D(width, height, panelWidth, panelHeight), dp_cr(nullptr) {
    if (width > 0 && height > 0) {
      initDrawing();
      needs_init_ = false;
    }
    initTextDrawer(noFreetype);
    df_noFreetype = noFreetype;
  }
  MolDraw2DCairo(const MolDraw2DCairo &) = delete;
  MolDraw2DCairo(MolDraw2DCairo &&) = delete;
  MolDraw2DCairo &operator=(const MolDraw2DCairo &) = delete;
  MolDraw2DCairo &operator=(MolDraw2DCairo &&) = delete;
  ~MolDraw2DCairo() {
    if (dp_cr) {
      if (cairo_get_reference_count(dp_cr) > 0) {
        cairo_destroy(dp_cr);
      }
      dp_cr = nullptr;
    }
  }

  void setColour(const DrawColour &col) override;

  // not sure if this goes here or if we should do a dtor since initDrawing() is
  // called in the ctor,
  // but we'll start here
  void finishDrawing();

  void drawLine(const Point2D &cds1, const Point2D &cds2,
                bool rawCoords = false) override;
  // void drawString( const std::string &str, const Point2D &cds );
  void drawPolygon(const std::vector<Point2D> &cds,
                   bool rawCoords = false) override;
  void clearDrawing() override;

  void drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                    const DrawColour &col1, const DrawColour &col2,
                    unsigned int nSegments = 16, double vertOffset = 0.05,
                    bool rawCoords = false) override;

  // returns the PNG data in a string
  std::string getDrawingText() const;
  // writes the PNG data to a file
  void writeDrawingText(const std::string &fName) const;

#if defined(WIN32) && !defined(RDK_BUILD_FREETYPE_SUPPORT)
  bool supportsAnnotations() const override { return false; }
#endif

 private:
  cairo_t *dp_cr = nullptr;
  bool df_noFreetype = false;

  void updateMetadata(const ROMol &mol, int confId) override;
  void updateMetadata(const ChemicalReaction &rxn) override;

  void initDrawing() override;
  void initTextDrawer(bool noFreetype) override;
  std::string addMetadataToPNG(const std::string &png) const;
  void updateMetadata(const ROMol &mol) const;
  void updateMetadata(const ChemicalReaction &rxn) const;
};
}  // namespace RDKit
#endif  // MOLDRAW2DCAIRO_H
