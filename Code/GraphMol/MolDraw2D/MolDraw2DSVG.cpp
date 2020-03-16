//
//  Copyright (C) 2015-2019 Greg Landrum
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
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <Geometry/point.h>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <sstream>

namespace RDKit {
namespace {
std::string DrawColourToSVG(const DrawColour &col) {
  const char *convert = "0123456789ABCDEF";
  std::string res(7, ' ');
  res[0] = '#';
  unsigned int v;
  unsigned int i = 1;
  v = int(255 * col.r);
  if (v > 255) {
    throw ValueErrorException(
        "elements of the color should be between 0 and 1");
  }
  res[i++] = convert[v / 16];
  res[i++] = convert[v % 16];
  v = int(255 * col.g);
  if (v > 255) {
    throw ValueErrorException(
        "elements of the color should be between 0 and 1");
  }
  res[i++] = convert[v / 16];
  res[i++] = convert[v % 16];
  v = int(255 * col.b);
  if (v > 255) {
    throw ValueErrorException(
        "elements of the color should be between 0 and 1");
  }
  res[i++] = convert[v / 16];
  res[i++] = convert[v % 16];
  return res;
}
}  // namespace

void MolDraw2DSVG::initDrawing() {
  d_os << "<?xml version='1.0' encoding='iso-8859-1'?>\n";
  d_os << "<svg version='1.1' baseProfile='full'\n      \
        xmlns='http://www.w3.org/2000/svg'\n              \
        xmlns:rdkit='http://www.rdkit.org/xml'\n              \
        xmlns:xlink='http://www.w3.org/1999/xlink'\n          \
        xml:space='preserve'\n";
  d_os << boost::format{"width='%1%px' height='%2%px' viewBox='0 0 %1% %2%'>\n"}
      % width() % height();
  d_os << "<!-- END OF HEADER -->\n";

  // d_os<<"<g transform='translate("<<width()*.05<<","<<height()*.05<<")
  // scale(.85,.85)'>";
}

// ****************************************************************************
void MolDraw2DSVG::finishDrawing() {
  // d_os << "</g>";
  d_os << "</svg>\n";
}

// ****************************************************************************
void MolDraw2DSVG::setColour(const DrawColour &col) {
  MolDraw2D::setColour(col);
}

void MolDraw2DSVG::drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                                const DrawColour &col1, const DrawColour &col2,
                                unsigned int nSegments, double vertOffset) {
  PRECONDITION(nSegments > 1, "too few segments");
  RDUNUSED_PARAM(col2);

  if (nSegments % 2) {
    ++nSegments;  // we're going to assume an even number of segments
  }
  setColour(col1);

  Point2D perp = calcPerpendicular(cds1, cds2);
  Point2D delta = (cds2 - cds1);
  perp *= vertOffset;
  delta /= nSegments;

  Point2D c1 = getDrawCoords(cds1);

  std::string col = DrawColourToSVG(colour());
  unsigned int width = getDrawLineWidth();
  d_os << "<path ";
  if (d_activeClass != "") {
    d_os << "class='" << d_activeClass << "' ";
  }
  d_os << "d='M" << c1.x << "," << c1.y;
  for (unsigned int i = 0; i < nSegments; ++i) {
    Point2D startpt = cds1 + delta * i;
    Point2D segpt = getDrawCoords(startpt + delta);
    Point2D cpt1 =
        getDrawCoords(startpt + delta / 3. + perp * (i % 2 ? -1 : 1));
    Point2D cpt2 =
        getDrawCoords(startpt + delta * 2. / 3. + perp * (i % 2 ? -1 : 1));
    d_os << " C" << cpt1.x << "," << cpt1.y << " " << cpt2.x << "," << cpt2.y
         << " " << segpt.x << "," << segpt.y;
  }
  d_os << "' ";

  d_os << "style='fill:none;stroke:" << col << ";stroke-width:" << width
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       << "'";
  d_os << " />\n";
}

void MolDraw2DSVG::drawBond(
    const ROMol &mol, const Bond *bond, int at1_idx, int at2_idx,
    const std::vector<int> *highlight_atoms,
    const std::map<int, DrawColour> *highlight_atom_map,
    const std::vector<int> *highlight_bonds,
    const std::map<int, DrawColour> *highlight_bond_map,
    const std::vector<std::pair<DrawColour, DrawColour> > *bond_colours) {
  PRECONDITION(bond, "bad bond");
  std::string o_class = d_activeClass;
  if (!d_activeClass.empty()) {
    d_activeClass += " ";
  }
  d_activeClass += boost::str(boost::format("bond-%d") % bond->getIdx());
  MolDraw2D::drawBond(mol, bond, at1_idx, at2_idx, highlight_atoms,
                      highlight_atom_map, highlight_bonds, highlight_bond_map,
                      bond_colours);
  d_activeClass = o_class;
};

// ****************************************************************************
void MolDraw2DSVG::drawLine(const Point2D &cds1, const Point2D &cds2) {
  Point2D c1 = getDrawCoords(cds1);
  Point2D c2 = getDrawCoords(cds2);
  std::string col = DrawColourToSVG(colour());
  unsigned int width = getDrawLineWidth();
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
  d_os << "<path ";
  if (!d_activeClass.empty()) {
    d_os << "class='" << d_activeClass << "' ";
  }
  d_os << "d='M " << c1.x << "," << c1.y << " L " << c2.x << "," << c2.y
       << "' ";
  d_os << "style='fill:none;fill-rule:evenodd;stroke:" << col
       << ";stroke-width:" << width
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       << dashString << "'";
  d_os << " />\n";
}

// ****************************************************************************
// draw the char, with the bottom left hand corner at cds
void MolDraw2DSVG::drawChar(char c, const Point2D &cds) {
  unsigned int fontSz = drawFontSize();
  std::string col = DrawColourToSVG(colour());

  d_os << "<text";
  d_os << " x='" << cds.x;
  d_os << "' y='" << cds.y << "'";
  d_os << " style='font-size:" << fontSz
       << "px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;"
          "font-family:sans-serif;text-anchor:start;"
       << "fill:" << col << "'";
  d_os << " >";
  d_os << c;
  d_os << "</text>";
}

// ****************************************************************************
void MolDraw2DSVG::drawPolygon(const std::vector<Point2D> &cds) {
  PRECONDITION(cds.size() >= 3, "must have at least three points");

  std::string col = DrawColourToSVG(colour());
  unsigned int width = getDrawLineWidth();
  std::string dashString = "";
  d_os << "<path ";
  if (d_activeClass != "") {
    d_os << "class='" << d_activeClass << "' ";
  }
  d_os << "d='M";
  Point2D c0 = getDrawCoords(cds[0]);
  d_os << " " << c0.x << "," << c0.y;
  for (unsigned int i = 1; i < cds.size(); ++i) {
    Point2D ci = getDrawCoords(cds[i]);
    d_os << " L " << ci.x << "," << ci.y;
  }
  if (fillPolys()) {
    // the Z closes the path which we don't want for unfilled polygons
    d_os << " Z' style='fill:" << col << ";fill-rule:evenodd;fill-opacity=" << colour().a
         << ";";
  } else {
    d_os << "' style='fill:none;";
  }

  d_os << "stroke:" << col << ";stroke-width:" << width
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:"
       << colour().a << ";" << dashString << "'";
  d_os << " />\n";
}

// ****************************************************************************
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
  unsigned int width = getDrawLineWidth();
  std::string dashString = "";
  d_os << "<ellipse"
       << " cx='" << cx << "'"
       << " cy='" << cy << "'"
       << " rx='" << w / 2 << "'"
       << " ry='" << h / 2 << "'";

  if (d_activeClass != "") {
    d_os << " class='" << d_activeClass << "'";
  }
  d_os << " style='";
  if (fillPolys()) {
    d_os << "fill:" << col << ";fill-rule:evenodd;";
  } else {
    d_os << "fill:none;";
  }

  d_os << "stroke:" << col << ";stroke-width:" << width
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       << dashString << "'";
  d_os << " />\n";
}

// ****************************************************************************
void MolDraw2DSVG::clearDrawing() {
  std::string col = DrawColourToSVG(drawOptions().backgroundColour);
  d_os << "<rect";
  d_os << " style='opacity:1.0;fill:" << col << ";stroke:none'";
  d_os << " width='" << width() << "' height='" << height() << "'";
  d_os << " x='0' y='0'";
  d_os << "> </rect>\n";
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

  double act_font_size = drawFontSize() / scale();
  for (int i = 0, is = label.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == label[i] && setStringDrawMode(label, draw_mode, i)) {
      continue;
    }

    label_height = act_font_size;
    double char_width =
        act_font_size *
        static_cast<double>(MolDraw2D_detail::char_widths[(int)label[i]]) /
        MolDraw2D_detail::char_widths[(int)'M'];
    if (TextDrawSubscript == draw_mode) {
      char_width *= 0.5;
      had_a_sub = true;
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

namespace {
void escape_xhtml(std::string &data) {
  boost::algorithm::replace_all(data, "&", "&amp;");
  boost::algorithm::replace_all(data, "\"", "&quot;");
  boost::algorithm::replace_all(data, "\'", "&apos;");
  boost::algorithm::replace_all(data, "<", "&lt;");
  boost::algorithm::replace_all(data, ">", "&gt;");
}
}  // namespace

// ****************************************************************************
void MolDraw2DSVG::drawString(const std::string &str, const Point2D &cds) {

  drawString(str, cds, MIDDLE);

}

// ****************************************************************************
// draws the string aligned as requested.
void MolDraw2DSVG::drawString(const std::string &str, const Point2D &cds,
                              AlignType align) {
  unsigned int fontSz = drawFontSize();
  double string_width, string_height;
  getStringSize(str, string_width, string_height);

  double draw_x = cds.x;
  double draw_y = cds.y;

  // for debugging text output
#if 0
  draw_x = cds.x - string_width / 2.0;
  draw_y = cds.y - string_height / 2.0;
  DrawColour tcolour =colour();
  setColour(DrawColour(.8,.8,.8));
  std::vector<Point2D> poly;
  poly.push_back(Point2D(draw_x,draw_y));
  poly.push_back(Point2D(draw_x+string_width,draw_y));
  poly.push_back(Point2D(draw_x+string_width,draw_y+string_height));
  poly.push_back(Point2D(draw_x,draw_y+string_height));
  drawPolygon(poly);
  setColour(tcolour);
  draw_x = cds.x;
  draw_y = cds.y;
#endif
  std::string col = DrawColourToSVG(colour());

  std::string text_anchor = "middle";
  double tmult = 0.0;
  if(align == END) {
    text_anchor = "end";
    tmult = -1.0;
  } else if(align == START) {
    text_anchor = "start";
    tmult = 1.0;
  }
  Point2D draw_coords = getDrawCoords(Point2D(draw_x, draw_y));
  // fonts are laid out with room for wider letters like W and hanging bits like g.
  // Very few atomic symbols need to care about this, and common ones look a bit
  // out of line.  For example O sits to the left of a double bond.  This is an
  // empirical tweak to push it back a bit.  Use the string_height in x and y
  // as it's just a correction factor based on the scaled font size.
  draw_coords.x += string_height * tmult * 0.1 * scale();
  draw_coords.y += string_height * 0.15  *scale();

  d_os << "<text dominant-baseline=\"central\" text-anchor=\""
       << text_anchor << "\"";
  d_os << " x='" << draw_coords.x;
  d_os << "' y='" << draw_coords.y << "'";

  if (!d_activeClass.empty()) {
    d_os << " class='" << d_activeClass << "'";
  }
  d_os << " style='font-size:" << fontSz
       << "px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;"
          "font-family:sans-serif;"
       << "fill:" << col << "'";
  d_os << " >";

  TextDrawType draw_mode =
      TextDrawNormal;  // 0 for normal, 1 for superscript, 2 for subscript
  std::string span;
  bool first_span = true;
  auto write_span = [&]() {
    if (!first_span) {
      escape_xhtml(span);
      d_os << span << "</tspan>";
      span = "";
    }
    first_span = false;
  };

  for (int i = 0, is = str.length(); i < is; ++i) {
    // setStringDrawMode moves i along to the end of any <sub> or <sup>
    // markup
    if ('<' == str[i] && setStringDrawMode(str, draw_mode, i)) {
      write_span();
      d_os << "<tspan";
      switch (draw_mode) {
        // To save people time later - on macOS Catalina, at least, Firefox
        // renders the superscript as a subscript.  It's fine on Safari.
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
      d_os << "<tspan>";
      span = "";
    }
    span += str[i];
  }
  escape_xhtml(span);
  d_os << span << "</tspan>";
  d_os << "</text>\n";
}

// ****************************************************************************
void MolDraw2DSVG::alignString(const std::string &str, const std::string &align_char,
                               int align, const Point2D &in_cds,
                               Point2D &out_cds) const {

  RDUNUSED_PARAM(str);

  if(align != 0 && align != 1) {
    out_cds = in_cds;
    return;
  }

  double ac_width, ac_height;
  getStringSize(align_char, ac_width, ac_height);
  // align == 0 is left align - first char to go at in_cds.
  double dir = align == 0 ? 1.0 : -1.0;
  // this works with SVG, so long as we use the correct text anchor -
  // W => end, E => start, N, S => middle
  out_cds.x = in_cds.x - dir * 0.5 * ac_width;
  // and this assumes dominant-baseline="central"
  out_cds.y = in_cds.y;

}

// ****************************************************************************
static const char *RDKIT_SVG_VERSION = "0.9";
void MolDraw2DSVG::addMoleculeMetadata(const ROMol &mol, int confId) const {
  PRECONDITION(d_os, "no output stream");
  d_os << "<metadata>" << std::endl;
  d_os << "<rdkit:mol"
       << " xmlns:rdkit = \"http://www.rdkit.org/xml\""
       << " version=\"" << RDKIT_SVG_VERSION << "\""
       << ">" << std::endl;
  for (const auto atom : mol.atoms()) {
    d_os << "<rdkit:atom idx=\"" << atom->getIdx() + 1 << "\"";
    bool doKekule = false, allHsExplicit = true, isomericSmiles = true;
    d_os << " atom-smiles=\""
         << SmilesWrite::GetAtomSmiles(atom, doKekule, nullptr, allHsExplicit,
                                       isomericSmiles)
         << "\"";
    auto tag = boost::str(boost::format("_atomdrawpos_%d") % confId);

    const Conformer &conf = mol.getConformer(confId);
    RDGeom::Point3D pos = conf.getAtomPos(atom->getIdx());

    Point2D dpos(pos.x, pos.y);
    if (atom->hasProp(tag)) {
      dpos = atom->getProp<Point2D>(tag);
    } else {
      dpos = getDrawCoords(dpos);
    }
    d_os << " drawing-x=\"" << dpos.x << "\""
         << " drawing-y=\"" << dpos.y << "\"";
    d_os << " x=\"" << pos.x << "\""
         << " y=\"" << pos.y << "\""
         << " z=\"" << pos.z << "\"";

    d_os << " />" << std::endl;
  }
  for (const auto bond : mol.bonds()) {
    d_os << "<rdkit:bond idx=\"" << bond->getIdx() + 1 << "\"";
    d_os << " begin-atom-idx=\"" << bond->getBeginAtomIdx() + 1 << "\"";
    d_os << " end-atom-idx=\"" << bond->getEndAtomIdx() + 1 << "\"";
    bool doKekule = false, allBondsExplicit = true;
    d_os << " bond-smiles=\""
         << SmilesWrite::GetBondSmiles(bond, -1, doKekule, allBondsExplicit)
         << "\"";
    d_os << " />" << std::endl;
  }
  d_os << "</rdkit:mol></metadata>" << std::endl;
}

void MolDraw2DSVG::addMoleculeMetadata(const std::vector<ROMol *> &mols,
                                       const std::vector<int> confIds) const {
  for (unsigned int i = 0; i < mols.size(); ++i) {
    int confId = -1;
    if (confIds.size() == mols.size()) {
      confId = confIds[i];
    }
    addMoleculeMetadata(*(mols[i]), confId);
  }
};

void MolDraw2DSVG::tagAtoms(const ROMol &mol, double radius,
                            const std::map<std::string, std::string> &events) {
  PRECONDITION(d_os, "no output stream");
  // first bonds so that they are under the atoms
  for (const auto &bond : mol.bonds()) {
    auto this_idx = bond->getIdx();
    auto a1pos = getDrawCoords(atomCoords()[bond->getBeginAtomIdx()]);
    auto a2pos = getDrawCoords(atomCoords()[bond->getEndAtomIdx()]);
    auto width = 2 + lineWidth();
    d_os << "<path "
         << " d='M " << a1pos.x << "," << a1pos.y << " L " << a2pos.x << ","
         << a2pos.y << "'";
    d_os << " class='bond-selector bond-" << this_idx;
    if (d_activeClass != "") {
      d_os << " " << d_activeClass;
    }
    d_os << "'";
    d_os << " style='fill:#fff;stroke:#fff;stroke-width:" << width
         << "px;fill-opacity:0;"
            "stroke-opacity:0' ";
    d_os << "/>\n";
  }
  for (const auto &at : mol.atoms()) {
    auto this_idx = at->getIdx();
    auto pos = getDrawCoords(atomCoords()[this_idx]);
    d_os << "<circle "
         << " cx='" << pos.x << "'"
         << " cy='" << pos.y << "'"
         << " r='" << (scale() * radius) << "'";
    d_os << " class='atom-selector atom-" << this_idx;
    if (d_activeClass != "") {
      d_os << " " << d_activeClass;
    }
    d_os << "'";
    d_os << " style='fill:#fff;stroke:#fff;stroke-width:1px;fill-opacity:0;"
            "stroke-opacity:0' ";
    for (const auto &event : events) {
      d_os << " " << event.first << "='" << event.second << "(" << this_idx
           << ");"
           << "'";
    }
    d_os << "/>\n";
  }
}
}  // namespace RDKit
