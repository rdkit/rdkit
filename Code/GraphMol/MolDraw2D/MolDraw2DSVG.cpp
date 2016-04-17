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

#include "MolDraw2DSVG.h"
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <sstream>

namespace RDKit {
namespace {
std::string DrawColourToSVG(const DrawColour &col) {
  const char *convert = "0123456789ABCDEF";
  std::string res(7, ' ');
  res[0] = '#';
  unsigned int v;
  unsigned int i = 1;
  v = int(255 * col.get<0>());
  res[i++] = convert[v / 16];
  res[i++] = convert[v % 16];
  v = int(255 * col.get<1>());
  res[i++] = convert[v / 16];
  res[i++] = convert[v % 16];
  v = int(255 * col.get<2>());
  res[i++] = convert[v / 16];
  res[i++] = convert[v % 16];
  return res;
}
}

void MolDraw2DSVG::initDrawing() {
  d_os << "<?xml version='1.0' encoding='iso-8859-1'?>\n";
  d_os << "<svg:svg version='1.1' baseProfile='full'\n      \
        xmlns:svg='http://www.w3.org/2000/svg'\n              \
        xmlns:rdkit='http://www.rdkit.org/xml'\n              \
        xmlns:xlink='http://www.w3.org/1999/xlink'\n          \
        xml:space='preserve'\n";
  d_os << "width='" << width() << "px' height='" << height() << "px' >\n";
  // d_os<<"<svg:g transform='translate("<<width()*.05<<","<<height()*.05<<")
  // scale(.85,.85)'>";
}

// ****************************************************************************
void MolDraw2DSVG::finishDrawing() {
  // d_os << "</svg:g>";
  d_os << "</svg:svg>\n";
}

// ****************************************************************************
void MolDraw2DSVG::setColour(const DrawColour &col) {
  MolDraw2D::setColour(col);
}

// ****************************************************************************
void MolDraw2DSVG::drawLine(const Point2D &cds1, const Point2D &cds2) {
  Point2D c1 = getDrawCoords(cds1);
  Point2D c2 = getDrawCoords(cds2);
  std::string col = DrawColourToSVG(colour());
  unsigned int width = lineWidth();
  std::string dashString = "";
  const DashPattern &dashes = dash();
  if (dashes.size()) {
    std::stringstream dss;
    dss << ";stroke-dasharray:";
    std::copy(dashes.begin(), dashes.end() - 1,
              std::ostream_iterator<unsigned int>(dss, ","));
    dss << dashes.back();
    dashString = dss.str();
  }
  d_os << "<svg:path ";
  d_os << "d='M " << c1.x << "," << c1.y << " " << c2.x << "," << c2.y << "' ";
  d_os << "style='fill:none;fill-rule:evenodd;stroke:" << col
       << ";stroke-width:" << width
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       << dashString << "'";
  d_os << " />\n";
}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void MolDraw2DSVG::drawChar(char c, const Point2D &cds) {
  unsigned int fontSz = scale() * fontSize();
  std::string col = DrawColourToSVG(colour());

  d_os << "<svg:text";
  d_os << " x='" << cds.x;
  d_os << "' y='" << cds.y << "'";
  d_os << " style='font-size:" << fontSz
       << "px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;"
          "font-family:sans-serif;text-anchor:start;"
       << "fill:" << col << "'";
  d_os << " >";
  d_os << c;
  d_os << "</svg:text>";
}

// ****************************************************************************
void MolDraw2DSVG::drawPolygon(const std::vector<Point2D> &cds) {
  PRECONDITION(cds.size() >= 3, "must have at least three points");

  std::string col = DrawColourToSVG(colour());
  unsigned int width = lineWidth();
  std::string dashString = "";
  d_os << "<svg:path ";
  d_os << "d='M";
  Point2D c0 = getDrawCoords(cds[0]);
  d_os << " " << c0.x << "," << c0.y;
  for (unsigned int i = 1; i < cds.size(); ++i) {
    Point2D ci = getDrawCoords(cds[i]);
    d_os << " " << ci.x << "," << ci.y;
  }
  d_os << " " << c0.x << "," << c0.y;
  d_os << "' style='";
  if (fillPolys())
    d_os << "fill:" << col << ";fill-rule:evenodd";
  else
    d_os << "fill:none;";

  d_os << "stroke:" << col << ";stroke-width:" << width
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       << dashString << "'";
  d_os << " />\n";
}

void MolDraw2DSVG::drawEllipse(const Point2D &cds1, const Point2D &cds2) {
  Point2D c1 = getDrawCoords(cds1);
  Point2D c2 = getDrawCoords(cds2);
  double w = c2.x - c1.x;
  double h = c2.y - c1.y;
  double cx = c1.x + w / 2;
  double cy = c1.y + h / 2;
  w = w > 0 ? w : -1 * w;
  h = h > 0 ? h : -1 * h;

  std::string col = DrawColourToSVG(colour());
  unsigned int width = lineWidth();
  std::string dashString = "";
  d_os << "<svg:ellipse"
       << " cx='" << cx << "'"
       << " cy='" << cy << "'"
       << " rx='" << w / 2 << "'"
       << " ry='" << h / 2 << "'";

  d_os << " style='";
  if (fillPolys())
    d_os << "fill:" << col << ";fill-rule:evenodd";
  else
    d_os << "fill:none;";

  d_os << "stroke:" << col << ";stroke-width:" << width
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       << dashString << "'";
  d_os << " />\n";
}

// ****************************************************************************
void MolDraw2DSVG::clearDrawing() {
  std::string col = DrawColourToSVG(drawOptions().backgroundColour);
  d_os << "<svg:rect";
  d_os << " style='opacity:1.0;fill:" << col << ";stroke:none'";
  d_os << " width='" << width() << "' height='" << height() << "'";
  d_os << " x='0' y='0'";
  d_os << "> </svg:rect>\n";
}

// ****************************************************************************
void MolDraw2DSVG::setFontSize(double new_size) {
  MolDraw2D::setFontSize(new_size);
  // double font_size_in_points = fontSize() * scale();
}

// ****************************************************************************
// using the current scale, work out the size of the label in molecule
// coordinates
void MolDraw2DSVG::getStringSize(const std::string &label, double &label_width,
                                 double &label_height) const {
  label_width = 0.0;
  label_height = 0.0;

  TextDrawType draw_mode = TextDrawNormal;

  bool had_a_super = false;
  bool had_a_sub = false;

  for (int i = 0, is = label.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == label[i] && setStringDrawMode(label, draw_mode, i)) {
      continue;
    }

    label_height = fontSize();
    double char_width =
        fontSize() *
        static_cast<double>(MolDraw2D_detail::char_widths[(int)label[i]]) /
        MolDraw2D_detail::char_widths[(int)'M'];
    if (TextDrawSubscript == draw_mode) {
      char_width *= 0.5;
      had_a_sub =true;
    } else if (TextDrawSuperscript == draw_mode) {
      char_width *= 0.5;
      had_a_super = true;
    }
    label_width += char_width;
  }

  // subscript keeps its bottom in line with the bottom of the bit chars,
  // superscript goes above the original char top by a bit (empirical)
  if (had_a_super) {
    label_height *= 1.1;
  }
  if (had_a_sub) {
    label_height *= 1.1;
  }
}

// ****************************************************************************
// draws the string centred on cds
void MolDraw2DSVG::drawString(const std::string &str, const Point2D &cds) {
  unsigned int fontSz = scale() * fontSize();

  double string_width, string_height;
  getStringSize(str, string_width, string_height);

  double draw_x = cds.x - string_width / 2.0;
  double draw_y = cds.y - string_height / 2.0;

#if 0
  // for debugging text output
  DrawColour tcolour =colour();
  setColour(DrawColour(.8,.8,.8));
  std::vector<Point2D> poly;
  poly.push_back(Point2D(draw_x,draw_y));
  poly.push_back(Point2D(draw_x+string_width,draw_y));
  poly.push_back(Point2D(draw_x+string_width,draw_y+string_height));
  poly.push_back(Point2D(draw_x,draw_y+string_height));
  drawPolygon(poly);
  setColour(tcolour);
#endif
  std::string col = DrawColourToSVG(colour());

  Point2D draw_coords = getDrawCoords(Point2D(draw_x, draw_y));

  d_os << "<svg:text";
  d_os << " x='" << draw_coords.x;

  d_os << "' y='" << draw_coords.y << "'";

  d_os << " style='font-size:" << fontSz
       << "px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;"
          "font-family:sans-serif;text-anchor:start;"
       << "fill:" << col << "'";
  d_os << " >";

  TextDrawType draw_mode =
      TextDrawNormal;  // 0 for normal, 1 for superscript, 2 for subscript
  std::string span;
  bool first_span = true;
  for (int i = 0, is = str.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == str[i] && setStringDrawMode(str, draw_mode, i)) {
      if (!first_span) {
        d_os << span << "</svg:tspan>";
        span = "";
      }
      first_span = false;
      d_os << "<svg:tspan";
      switch (draw_mode) {
        case TextDrawSuperscript:
          d_os << " style='baseline-shift:super;font-size:" << fontSz * 0.75
               << "px;"
               << "'";
          break;
        case TextDrawSubscript:
          d_os << " style='baseline-shift:sub;font-size:" << fontSz * 0.75
               << "px;"
               << "'";
          break;
        default:
          break;
      }
      d_os << ">";
      continue;
    }
    if (first_span) {
      first_span = false;
      d_os << "<svg:tspan>";
      span = "";
    }
    span += str[i];
  }
  d_os << span << "</svg:tspan>";
  d_os << "</svg:text>\n";
}

void MolDraw2DSVG::tagAtoms(const ROMol &mol) {
  PRECONDITION(d_os, "no output stream");
  ROMol::VERTEX_ITER this_at, end_at;
  boost::tie(this_at, end_at) = mol.getVertices();
  while (this_at != end_at) {
    int this_idx = mol[*this_at]->getIdx();
    ++this_at;
    Point2D pos = getDrawCoords(atomCoords()[this_idx]);
    std::string lbl = atomSyms()[this_idx].first;

    d_os << "<rdkit:atom"
         << " idx=\"" << this_idx + 1 << "\"";
    if (lbl != "") {
      d_os << " label=\"" << lbl << "\"";
    }
    d_os << " x=\"" << pos.x << "\""
         << " y=\"" << pos.y << "\""
         << " />" << std::endl;
  }
}
}  // EO namespace RDKit
