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

#include <QPainter>
#include <QString>

using namespace boost;
using namespace std;

namespace RDKit {

// ****************************************************************************
MolDraw2DQt::MolDraw2DQt(int width, int height, QPainter &qp, int panelWidth,
                         int panelHeight)
    : MolDraw2D(width, height, panelWidth, panelHeight), qp_(qp) {}

// ****************************************************************************
void MolDraw2DQt::setColour(const DrawColour &col) {
  MolDraw2D::setColour(col);
  QColor this_col(int(255.0 * col.r), int(255.0 * col.g), int(255.0 * col.b),
                  int(255.0 * col.a));

  QPen pen(this_col);
  pen.setJoinStyle(Qt::RoundJoin);
  pen.setColor(this_col);
  qp_.setPen(pen);

  QBrush brush(this_col);
  brush.setStyle(Qt::SolidPattern);
  qp_.setBrush(brush);
}

// ****************************************************************************
void MolDraw2DQt::drawLine(const Point2D &cds1, const Point2D &cds2) {
  Point2D c1 = getDrawCoords(cds1);
  Point2D c2 = getDrawCoords(cds2);

  const DashPattern &dashes = dash();
  QPen pen = qp_.pen();
  if (dashes.size()) {
    QVector<qreal> dd;
    for (unsigned int di = 0; di < dashes.size(); ++di) dd << dashes[di];
    pen.setDashPattern(dd);
  } else {
    pen.setStyle(Qt::SolidLine);
  }
  pen.setWidth(lineWidth());
  qp_.setPen(pen);
  qp_.drawLine(QPointF(c1.x, c1.y), QPointF(c2.x, c2.y));
}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void MolDraw2DQt::drawChar(char c, const Point2D &cds) {
  QRectF br = qp_.boundingRect(0, 0, 100, 100, Qt::AlignLeft | Qt::AlignBottom,
                               QString(c));
  qp_.drawText(QRectF(cds.x, cds.y - br.height(), br.width(), br.height()),
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

  qp_.save();
  QBrush brush = qp_.brush();
  if (fillPolys())
    brush.setStyle(Qt::SolidPattern);
  else
    brush.setStyle(Qt::NoBrush);
  qp_.setBrush(brush);

  QPointF points[cds.size()];
  for (unsigned int i = 0; i < cds.size(); ++i) {
    Point2D lc = getDrawCoords(cds[i]);
    points[i] = QPointF(lc.x, lc.y);
  }
  qp_.drawConvexPolygon(points, cds.size());
  qp_.restore();
}

// ****************************************************************************
void MolDraw2DQt::clearDrawing() {
  QColor this_col(int(255.0 * drawOptions().backgroundColour.r),
                  int(255.0 * drawOptions().backgroundColour.g),
                  int(255.0 * drawOptions().backgroundColour.b),
                  int(255.0 * drawOptions().backgroundColour.a));

  qp_.setBackground(QBrush(this_col));
  qp_.fillRect(0, 0, width(), height(), this_col);
}

// ****************************************************************************
void MolDraw2DQt::setFontSize(double new_size) {
  MolDraw2D::setFontSize(new_size);
  double font_size_in_points = fontSize() * scale();
#ifdef NOTYET
  cout << "initial font size in points : " << qp_.font().pointSizeF() << endl;
  cout << "font_size_in_points : " << font_size_in_points << endl;
#endif
  QFont font(qp_.font());
  font.setPointSizeF(font_size_in_points);
  qp_.setFont(font);

  while (1) {
    double old_font_size_in_points = font_size_in_points;
    double font_size_in_points = fontSize() * scale();
    if (fabs(font_size_in_points - old_font_size_in_points) < 0.1) {
      break;
    }
    QFont font(qp_.font());
    font.setPointSizeF(font_size_in_points);
    qp_.setFont(font);
    calculateScale();
  }
}

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
    QRectF br = qp_.boundingRect(0, 0, 100, 100,
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

}  // namespace RDKit
