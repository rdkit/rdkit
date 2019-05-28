//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Author: Igor Filippov based on the work of David Cosgrove (AstraZeneca)
//
// This is a concrete class derived from MolDraw2D that uses RDKit to draw a
// molecule into a wxDC

#include <RDGeneral/export.h>
#ifndef MOLDRAW2DWX_H
#define MOLDRAW2DWX_H

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <wx/dc.h>
#include <wx/font.h>
#include <wx/pen.h>
#include <wx/colour.h>
#include <wx/brush.h>

// ****************************************************************************

namespace RDKit {

class RDKIT_MOLDRAW2D_EXPORT MolDraw2Dwx : public MolDraw2D {
 public:
  MolDraw2Dwx(int width, int height, wxDC &dc, int panelWidth = -1,
              int panelHeight = -1)
      : MolDraw2D(width, height, panelWidth, panelHeight), m_dc(dc) {
    // m_dc.SetFont(wxFont(10, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL,
    // wxFONTWEIGHT_NORMAL));
  }

  // set font size in molecule coordinate units. That's probably Angstrom for
  // RDKit. It will turned into drawing units using scale_, which might be
  // changed as a result, to make sure things still appear in the window.

  void setFontSize(double new_size) {
    MolDraw2D::setFontSize(new_size);
    double font_size_in_points = fontSize() * scale();
    wxFont font = m_dc.GetFont();
    // font.SetPointSize(font_size_in_points);
    font.SetPixelSize(wxSize(0, font_size_in_points));
    m_dc.SetFont(font);
  }

  void setColour(const DrawColour &col) {
    MolDraw2D::setColour(col);
    double r = col.get<0>();
    double g = col.get<1>();
    double b = col.get<2>();
    wxColour colour(r * 255, g * 255, b * 255);
    m_dc.SetTextForeground(colour);
    m_dc.SetPen(wxPen(colour));
    m_dc.SetBrush(wxBrush(colour));
  }

  void drawLine(const Point2D &cds1, const Point2D &cds2) {
    Point2D c1 = getDrawCoords(cds1);
    Point2D c2 = getDrawCoords(cds2);
    m_dc.DrawLine(c1.x, c1.y, c2.x, c2.y);
  }

  void drawChar(char c, const Point2D &cds) {
    m_dc.DrawText(wxString(c), cds.x, cds.y);
  }

  void drawPolygon(const std::vector<Point2D> &cds) {
    PRECONDITION(cds.size() >= 3, "must have at least three points");
    wxPoint lines[cds.size()];
    for (unsigned int i = 0; i < cds.size(); ++i) {
      Point2D c1 = getDrawCoords(cds[i]);
      lines[i] = wxPoint(c1.x, c1.y);
    }
    // FIX: deal with toggling fills
    m_dc.DrawPolygon(cds.size(), lines);
  };

  void clearDrawing() {
    const wxBrush &brush = m_dc.GetBrush();
    const wxPen &pen = m_dc.GetPen();
    setColour(drawOptions.backgroundColour);
    m_dc.DrawRectangle(0, 0, width(), height());
    m_dc.SetBrush(brush);
    m_dc.SetPen(pen);
  }

  // using the current scale, work out the size of the label in molecule
  // coordinates
  void getStringSize(const std::string &label, double &label_width,
                     double &label_height) const {
    if (m_dc.CanGetTextExtent()) {
      wxCoord width, height;
      m_dc.GetTextExtent(wxString(label), &width, &height);
      label_width = double(width) / scale();
      label_height = double(height) / scale();
    }
  }

 private:
  wxDC &m_dc;
};
}  // namespace RDKit
#endif  // MOLDRAW2DWX_H
