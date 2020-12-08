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
  d_qp.save();
  QBrush brush = d_qp.brush();
  if (fillPolys())
    brush.setStyle(Qt::SolidPattern);
  else
    brush.setStyle(Qt::NoBrush);
  d_qp.setBrush(brush);

  QPointF *points = new QPointF[cds.size()];
  for (unsigned int i = 0; i < cds.size(); ++i) {
    Point2D lc = getDrawCoords(cds[i]);
    points[i] = QPointF(lc.x, lc.y);
  }
  d_qp.drawConvexPolygon(points, cds.size());
  d_qp.restore();
  delete[] points;
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

}  // namespace RDKit
