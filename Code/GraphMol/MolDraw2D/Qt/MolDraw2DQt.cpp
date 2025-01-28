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
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <QPainter>
#include <QPainterPath>
#include <QString>

#ifdef RDK_BUILD_FREETYPE_SUPPORT
#include "DrawTextFTQt.h"
#else
#include "DrawTextQt.h"
#endif

namespace RDKit {

const char *rdkitQtVersion = RDK_QT_VERSION;

// ****************************************************************************
MolDraw2DQt::MolDraw2DQt(int width, int height, QPainter *qp, int panelWidth,
                         int panelHeight, bool noFreetype)
    : MolDraw2D(width, height, panelWidth, panelHeight), d_qp(qp) {
  PRECONDITION(width > 0, "bad width");
  PRECONDITION(height > 0, "bad height");
  initDrawing();
  initTextDrawer(noFreetype);
}

void MolDraw2DQt::initDrawing() {
  d_qp->setRenderHint(QPainter::RenderHint::Antialiasing);
}
void MolDraw2DQt::initTextDrawer(bool noFreetype) {
  double max_fnt_sz = drawOptions().maxFontSize;
  double min_fnt_sz = drawOptions().minFontSize;

  if (noFreetype) {
    text_drawer_.reset(
        new MolDraw2D_detail::DrawTextQt(max_fnt_sz, min_fnt_sz, d_qp));
  } else {
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    try {
      text_drawer_.reset(new MolDraw2D_detail::DrawTextFTQt(
          max_fnt_sz, min_fnt_sz, drawOptions().fontFile, d_qp));
    } catch (std::runtime_error &e) {
      BOOST_LOG(rdWarningLog)
          << e.what() << std::endl
          << "Falling back to native Qt text handling." << std::endl;
      text_drawer_.reset(
          new MolDraw2D_detail::DrawTextQt(max_fnt_sz, min_fnt_sz, d_qp));
    }
#else
    text_drawer_.reset(
        new MolDraw2D_detail::DrawTextQt(max_fnt_sz, min_fnt_sz, d_qp));
#endif
  }
  if (drawOptions().baseFontSize > 0.0) {
    text_drawer_->setBaseFontSize(drawOptions().baseFontSize);
  }
}

// ****************************************************************************
void MolDraw2DQt::setColour(const DrawColour &col) {
  MolDraw2D::setColour(col);
  QColor this_col(int(255.0 * col.r), int(255.0 * col.g), int(255.0 * col.b),
                  int(255.0 * col.a));

  QPen pen(this_col);
  pen.setJoinStyle(Qt::RoundJoin);
  d_qp->setPen(pen);

  QBrush brush(this_col);
  brush.setStyle(Qt::SolidPattern);
  d_qp->setBrush(brush);
}

namespace {
void setDashes(QPen &pen, const DashPattern &dashes) {
  if (dashes.size()) {
    QVector<qreal> dd;
    for (unsigned int di = 0; di < dashes.size(); ++di) dd << dashes[di];
    pen.setDashPattern(dd);
  } else {
    pen.setStyle(Qt::SolidLine);
    pen.setCapStyle(Qt::FlatCap);
  }
}
}  // namespace

// ****************************************************************************
void MolDraw2DQt::drawLine(const Point2D &cds1, const Point2D &cds2,
                           bool rawCoords) {
  Point2D c1 = rawCoords ? cds1 : getDrawCoords(cds1);
  Point2D c2 = rawCoords ? cds2 : getDrawCoords(cds2);

  QPen pen = d_qp->pen();
  setDashes(pen, dash());
  pen.setWidth(getDrawLineWidth());
  d_qp->setPen(pen);
  d_qp->drawLine(QPointF(c1.x, c1.y), QPointF(c2.x, c2.y));
}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void MolDraw2DQt::drawChar(char c, const Point2D &cds) {
  QRectF br = d_qp->boundingRect(0, 0, 100, 100,
                                 Qt::AlignLeft | Qt::AlignBottom, QString(c));
  d_qp->drawText(QRectF(cds.x, cds.y - br.height(), br.width(), br.height()),
                 Qt::AlignLeft | Qt::AlignBottom, QString(c), &br);
}

// ****************************************************************************
void MolDraw2DQt::drawPolygon(const std::vector<Point2D> &cds, bool rawCoords) {
  PRECONDITION(cds.size() >= 3, "must have at least three points");
  QPointF *points = new QPointF[cds.size()];
  for (unsigned int i = 0; i < cds.size(); ++i) {
    Point2D lc = rawCoords ? cds[i] : getDrawCoords(cds[i]);
    points[i] = QPointF(lc.x, lc.y);
  }
  QPen pen = d_qp->pen();
  setDashes(pen, dash());
  pen.setStyle(Qt::SolidLine);
  pen.setCapStyle(Qt::FlatCap);
  pen.setWidth(getDrawLineWidth());
  if (fillPolys()) {
    d_qp->save();
    d_qp->setPen(pen);
    QBrush brush = d_qp->brush();
    brush.setStyle(Qt::SolidPattern);
    d_qp->setBrush(brush);

    d_qp->drawConvexPolygon(points, cds.size());
    d_qp->restore();
  } else {
    QPainterPath path;
    path.moveTo(points[0]);
    for (unsigned int i = 0; i < cds.size(); ++i) {
      path.lineTo(points[i]);
    }
    d_qp->strokePath(path, pen);
  }
  delete[] points;
}

// ****************************************************************************
void MolDraw2DQt::drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                               const DrawColour &col1, const DrawColour &,
                               unsigned int nSegments, double vertOffset,
                               bool rawCoords) {
  PRECONDITION(nSegments > 1, "too few segments");

  auto segments =
      MolDraw2D_detail::getWavyLineSegments(cds1, cds2, nSegments, vertOffset);

  QPen pen = d_qp->pen();
  pen.setStyle(Qt::SolidLine);
  pen.setCapStyle(Qt::FlatCap);
  pen.setWidth(getDrawLineWidth());

  setColour(col1);

  QPainterPath path;

  auto c1 = std::get<0>(segments[0]);
  c1 = rawCoords ? c1 : getDrawCoords(c1);

  path.moveTo(QPointF(c1.x, c1.y));
  for (unsigned int i = 0; i < nSegments; ++i) {
    auto cpt1 = std::get<1>(segments[i]);
    cpt1 = rawCoords ? cpt1 : getDrawCoords(cpt1);
    auto cpt2 = std::get<2>(segments[i]);
    cpt2 = rawCoords ? cpt2 : getDrawCoords(cpt2);
    auto segpt = std::get<3>(segments[i]);
    segpt = rawCoords ? segpt : getDrawCoords(segpt);
    path.cubicTo(QPointF(cpt1.x, cpt1.y), QPointF(cpt2.x, cpt2.y),
                 QPointF(segpt.x, segpt.y));
  }
  d_qp->strokePath(path, pen);
}

// ****************************************************************************
void MolDraw2DQt::clearDrawing() {
  MolDraw2D::clearDrawing();

  QColor this_col(int(255.0 * drawOptions().backgroundColour.r),
                  int(255.0 * drawOptions().backgroundColour.g),
                  int(255.0 * drawOptions().backgroundColour.b),
                  int(255.0 * drawOptions().backgroundColour.a));

  d_qp->setBackground(QBrush(this_col));
  d_qp->fillRect(offset().x, offset().y, width(), height(), this_col);
}
}  // namespace RDKit
