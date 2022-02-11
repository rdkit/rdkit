//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
#include "DrawTextQt.h"

#include <QPainter>
#include <QCoreApplication>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>

using namespace std;

namespace RDKit {
namespace MolDraw2D_detail {

// ****************************************************************************
DrawTextQt::DrawTextQt(double max_fnt_sz, double min_fnt_sz, QPainter *qp)
    : DrawTextNotFT(max_fnt_sz, min_fnt_sz), d_qp(qp) {
  PRECONDITION(
      QCoreApplication::instance(),
      "need a global QGuiApplication instance to use the Qt font system");
  QFont font("Sans Serif", fontSize());
  d_qp->setFont(font);
}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void DrawTextQt::drawChar(char c, const Point2D &cds) {
  auto font = d_qp->font();
  font.setPixelSize(fontSize());
  d_qp->setFont(font);
  QPointF loc(cds.x, cds.y);
  QChar qc(c);
  QString text(qc);
  const auto col = colour();
  QColor qcol(int(255.0 * col.r), int(255.0 * col.g), int(255.0 * col.b),
              int(255.0 * col.a));
  d_qp->setPen(qcol);
  d_qp->drawText(loc, text);
}

// ****************************************************************************
void DrawTextQt::getStringRects(const string &text,
                                vector<shared_ptr<StringRect>> &rects,
                                vector<TextDrawType> &draw_modes,
                                vector<char> &draw_chars) const {
  double running_x = 0.0;
  double act_font_size = fontSize();
  double char_height;
  double max_width = 0.0;
  TextDrawType draw_mode = TextDrawType::TextDrawNormal;
  QFontMetrics qfm(d_qp->font());
  for (size_t i = 0; i < text.length(); ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == text[i] && setStringDrawMode(text, draw_mode, i)) {
      continue;
    }
    draw_modes.push_back(draw_mode);
    draw_chars.push_back(text[i]);

    max_width = std::max(
        max_width,
        static_cast<double>(qfm.boundingRect(QChar(draw_chars[i])).width()));
  }
  for (size_t i = 0; i < draw_chars.size(); ++i) {
    auto br = qfm.boundingRect(QChar(draw_chars[i]));
    double char_width = 0.6 * act_font_size * br.width() / max_width;
    // Absent real formatted string handling this is something of an empirical
    // bodge.
    if (draw_chars[i] == '+') {
      char_height = 0.6 * act_font_size;
    } else if (draw_chars[i] == '-') {
      char_height = 0.4 * act_font_size;
    } else {
      char_height = 0.8 * act_font_size;
    }
    double cscale = selectScaleFactor(draw_chars[i], draw_modes[i]);
    char_height *= cscale;
    char_width *= cscale;
    Point2D offset(char_width / 2, char_height / 2);
    if (draw_chars[i] == '+' || draw_chars[i] == '-') {
      offset.y /= 2.0;
    }
    Point2D g_centre(char_width / 2, char_height / 2);
    rects.push_back(std::shared_ptr<StringRect>(
        new StringRect(offset, g_centre, char_width, char_height)));
    rects.back()->trans_.x += running_x;
    // empirical spacing.
    if (draw_modes[i] != TextDrawType::TextDrawNormal) {
      running_x += char_width * 1.05;
    } else {
      running_x += char_width * 1.15;
    }
  }
  for (auto r : rects) {
    r->g_centre_.y = act_font_size - r->g_centre_.y;
    r->offset_.y = act_font_size / 2.0;
  }

  adjustStringRectsForSuperSubScript(draw_modes, rects);
}

}  // namespace MolDraw2D_detail
}  // namespace RDKit