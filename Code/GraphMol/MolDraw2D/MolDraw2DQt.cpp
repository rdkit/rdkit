//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (AstraZeneca)
// 19th June 2014
//

#include "MolDraw2DQt.h"
#ifdef RDK_BUILD_FREETYPE_SUPPORT
#include <GraphMol/MolDraw2D/DrawTextFTQt.h>
#else
#include <GraphMol/MolDraw2D/DrawTextQt.h>
#endif

#include <QPainter>
#include <QString>

using namespace boost;
using namespace std;

namespace RDKit {

// ****************************************************************************
MolDraw2DQt::MolDraw2DQt(int width, int height, QPainter &qp, int panelWidth,
                         int panelHeight, bool noFreetype)
    : MolDraw2D(width, height, panelWidth, panelHeight), d_qp(qp) {
  initDrawing();
  initTextDrawer(noFreetype);
}

void MolDraw2DQt::initDrawing() {
  d_qp.setRenderHint(QPainter::RenderHint::Antialiasing);
}
void MolDraw2DQt::initTextDrawer(bool noFreetype) {
  double max_fnt_sz = drawOptions().maxFontSize;
  double min_fnt_sz = drawOptions().minFontSize;

  if (noFreetype) {
    text_drawer_.reset(new DrawTextQt(max_fnt_sz, min_fnt_sz, d_qp));
  } else {
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    try {
      text_drawer_.reset(new DrawTextFTQt(max_fnt_sz, min_fnt_sz,
                                          drawOptions().fontFile, d_qp));
    } catch (std::runtime_error &e) {
      BOOST_LOG(rdWarningLog)
          << e.what() << std::endl
          << "Falling back to native Qt text handling." << std::endl;
      text_drawer_.reset(new DrawTextQt(max_fnt_sz, min_fnt_sz, d_qp));
    }
#else
    text_drawer_.reset(new DrawTextQt(max_fnt_sz, min_fnt_sz, d_qp));
#endif
  }
}

// ****************************************************************************
void MolDraw2DQt::setColour(const DrawColour &col) {
  MolDraw2D::setColour(col);
  QColor this_col(int(255.0 * col.r), int(255.0 * col.g), int(255.0 * col.b),
                  int(255.0 * col.a));

  QPen pen(this_col);
  pen.setJoinStyle(Qt::RoundJoin);
  d_qp.setPen(pen);

  QBrush brush(this_col);
  brush.setStyle(Qt::SolidPattern);
  d_qp.setBrush(brush);
}

// ****************************************************************************
void MolDraw2DQt::drawLine(const Point2D &cds1, const Point2D &cds2) {
  Point2D c1 = getDrawCoords(cds1);
  Point2D c2 = getDrawCoords(cds2);

  const DashPattern &dashes = dash();
  QPen pen = d_qp.pen();
  if (dashes.size()) {
    QVector<qreal> dd;
    for (unsigned int di = 0; di < dashes.size(); ++di) dd << dashes[di];
    pen.setDashPattern(dd);
  } else {
    pen.setStyle(Qt::SolidLine);
  }
  pen.setWidth(lineWidth());
  d_qp.setPen(pen);
  d_qp.drawLine(QPointF(c1.x, c1.y), QPointF(c2.x, c2.y));
}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void MolDraw2DQt::drawChar(char c, const Point2D &cds) {
  QRectF br = d_qp.boundingRect(0, 0, 100, 100, Qt::AlignLeft | Qt::AlignBottom,
                                QString(c));
  d_qp.drawText(QRectF(cds.x, cds.y - br.height(), br.width(), br.height()),
                Qt::AlignLeft | Qt::AlignBottom, QString(c), &br);
}

// ****************************************************************************
void MolDraw2DQt::drawPolygon(const vector<Point2D> &cds) {
  PRECONDITION(cds.size() >= 3, "must have at least three points");
#ifdef NOTYET
  QBrush brush("Black");
  brush.setStyle(Qt::SolidPattern);
  DrawColour cc = colour();
  brush.setColor(QColor(255.0 * cc.r, 255.0 * cc.g, 255.0 * cc.b));
#endif

  d_qp.save();
  QBrush brush = d_qp.brush();
  if (fillPolys())
    brush.setStyle(Qt::SolidPattern);
  else
    brush.setStyle(Qt::NoBrush);
  d_qp.setBrush(brush);

  QPointF points[cds.size()];
  for (unsigned int i = 0; i < cds.size(); ++i) {
    Point2D lc = getDrawCoords(cds[i]);
    points[i] = QPointF(lc.x, lc.y);
  }
  d_qp.drawConvexPolygon(points, cds.size());
  d_qp.restore();
}

// ****************************************************************************
void MolDraw2DQt::clearDrawing() {
  QColor this_col(int(255.0 * drawOptions().backgroundColour.r),
                  int(255.0 * drawOptions().backgroundColour.g),
                  int(255.0 * drawOptions().backgroundColour.b),
                  int(255.0 * drawOptions().backgroundColour.a));

  d_qp.setBackground(QBrush(this_col));
  d_qp.fillRect(offset().x, offset().y, width(), height(), this_col);
}

#if 0
// ****************************************************************************
void MolDraw2DQt::setFontSize(double new_size) {
  MolDraw2D::setFontSize(new_size);
  double font_size_in_points = fontSize() * scale();
#ifdef NOTYET
  cout << "initial font size in points : " << d_qp.font().pointSizeF() << endl;
  cout << "font_size_in_points : " << font_size_in_points << endl;
#endif
  QFont font(d_qp.font());
  font.setPointSizeF(font_size_in_points);
  d_qp.setFont(font);

  while (1) {
    double old_font_size_in_points = font_size_in_points;
    double font_size_in_points = fontSize() * scale();
    if (fabs(font_size_in_points - old_font_size_in_points) < 0.1) {
      break;
    }
    QFont font(d_qp.font());
    font.setPointSizeF(font_size_in_points);
    d_qp.setFont(font);
    calculateScale();
  }
}
#endif

#if 0
// ****************************************************************************
// using the current scale, work out the size of the label in molecule
// coordinates
void MolDraw2DQt::getStringSize(const string &label, double &label_width,
                                double &label_height) const {
  label_width = 0.0;
  label_height = 0.0;

  TextDrawType draw_mode =
      TextDrawNormal;  // 0 for normal, 1 for superscript, 2 for subscript
  QString next_char(" ");
  bool had_a_super = false;

  for (int i = 0, is = label.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == label[i] && setStringDrawMode(label, draw_mode, i)) {
      continue;
    }

    next_char[0] = label[i];
    QRectF br = d_qp.boundingRect(0, 0, 100, 100,
                                 Qt::AlignBottom | Qt::AlignLeft, next_char);
    label_height = br.height() / scale();
    double char_width = br.width() / scale();
    if (TextDrawSubscript == draw_mode) {
      char_width *= 0.5;
    } else if (TextDrawSuperscript == draw_mode) {
      char_width *= 0.5;
      had_a_super = true;
    }
    label_width += char_width;
  }

  // subscript keeps its bottom in line with the bottom of the bit chars,
  // superscript goes above the original char top by a quarter
  if (had_a_super) {
    label_height *= 1.25;
  }
}
#endif

}  // namespace RDKit
