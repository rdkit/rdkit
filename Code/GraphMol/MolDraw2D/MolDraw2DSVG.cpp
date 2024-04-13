//
//  Copyright (C) 2015-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// derived from Dave Cosgrove's MolDraw2D
//

#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/DrawText.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <Geometry/point.h>
#include <GraphMol/MolDraw2D/DrawTextSVG.h>
#ifdef RDK_BUILD_FREETYPE_SUPPORT
#include <GraphMol/MolDraw2D/DrawTextFTSVG.h>
#endif

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <sstream>

namespace RDKit {
namespace {
template <class t_obj>
void outputTagClasses(const t_obj *obj, std::ostream &d_os,
                      const std::string &d_activeClass) {
  if (!d_activeClass.empty()) {
    d_os << " " << d_activeClass;
  }
  if (obj->hasProp("_tagClass")) {
    std::string value;
    obj->getProp("_tagClass", value);
    std::replace(value.begin(), value.end(), '\"', '_');
    std::replace(value.begin(), value.end(), '\'', '_');
    std::replace(value.begin(), value.end(), '.', '_');
    d_os << " " << value;
  }
}

template <class t_obj>
void outputMetaData(const t_obj *obj, std::ostream &d_os) {
  std::string value;
  for (const auto &prop : obj->getPropList()) {
    if (prop.length() < 11 || prop.rfind("_metaData-", 0) != 0) {
      continue;
    }
    obj->getProp(prop, value);
    boost::replace_all(value, "\"", "&quot;");
    d_os << " " << prop.substr(10) << "=\"" << value << "\"";
  }
}
}  // namespace

std::string DrawColourToSVG(const DrawColour &col) {
  const char *convert = "0123456789ABCDEF";
  bool hasAlpha = 1.0 - col.a > 1e-3;
  std::string res(hasAlpha ? 9 : 7, ' ');
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
  if (hasAlpha) {
    v = int(255 * col.a);
    if (v > 255) {
      throw ValueErrorException(
          "elements of the color should be between 0 and 1");
    }
    res[i++] = convert[v / 16];
    res[i++] = convert[v % 16];
  }
  return res;
}

// ****************************************************************************
MolDraw2DSVG::MolDraw2DSVG(int width, int height, std::ostream &os,
                           int panelWidth, int panelHeight, bool noFreetype)
    : MolDraw2D(width, height, panelWidth, panelHeight), d_os(os) {
  if (width > 0 && height > 0) {
    initDrawing();
    needs_init_ = false;
  }
  initTextDrawer(noFreetype);
}

// ****************************************************************************
MolDraw2DSVG::MolDraw2DSVG(int width, int height, int panelWidth,
                           int panelHeight, bool noFreetype)
    : MolDraw2D(width, height, panelWidth, panelHeight), d_os(d_ss) {
  if (width > 0 && height > 0) {
    initDrawing();
    needs_init_ = false;
  }
  initTextDrawer(noFreetype);
}

// ****************************************************************************
void MolDraw2DSVG::initDrawing() {
  d_os << "<?xml version='1.0' encoding='iso-8859-1'?>\n";
  d_os << "<svg version='1.1' baseProfile='full'\n      \
        xmlns='http://www.w3.org/2000/svg'\n              \
        xmlns:rdkit='http://www.rdkit.org/xml'\n              \
        xmlns:xlink='http://www.w3.org/1999/xlink'\n          \
        xml:space='preserve'\n";
  d_os
      << boost::format{"width='%1%px' height='%2%px' viewBox='0 0 %1% %2%'>\n"} %
             width() % height();
  d_os << "<!-- END OF HEADER -->\n";
}

// ****************************************************************************
void MolDraw2DSVG::initTextDrawer(bool noFreetype) {
  double max_fnt_sz = drawOptions().maxFontSize;
  double min_fnt_sz = drawOptions().minFontSize;

  if (noFreetype) {
    text_drawer_.reset(new MolDraw2D_detail::DrawTextSVG(max_fnt_sz, min_fnt_sz,
                                                         d_os, d_activeClass));
  } else {
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    try {
      text_drawer_.reset(new MolDraw2D_detail::DrawTextFTSVG(
          max_fnt_sz, min_fnt_sz, drawOptions().fontFile, d_os, d_activeClass));
    } catch (std::runtime_error &e) {
      BOOST_LOG(rdWarningLog)
          << e.what() << std::endl
          << "Falling back to native SVG text handling." << std::endl;
      text_drawer_.reset(new MolDraw2D_detail::DrawTextSVG(
          max_fnt_sz, min_fnt_sz, d_os, d_activeClass));
    }
#else
    text_drawer_.reset(new MolDraw2D_detail::DrawTextSVG(max_fnt_sz, min_fnt_sz,
                                                         d_os, d_activeClass));
#endif
  }
  if (drawOptions().baseFontSize > 0.0) {
    text_drawer_->setBaseFontSize(drawOptions().baseFontSize);
  }
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

// ****************************************************************************
void MolDraw2DSVG::drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                                const DrawColour &col1, const DrawColour &,
                                unsigned int nSegments, double vertOffset,
                                bool rawCoords) {
  PRECONDITION(nSegments > 1, "too few segments");

  setColour(col1);

  auto segments =
      MolDraw2D_detail::getWavyLineSegments(cds1, cds2, nSegments, vertOffset);

  std::string col = DrawColourToSVG(colour());
  double width = getDrawLineWidth();
  d_os << "<path ";
  outputClasses();

  auto c1 = std::get<0>(segments[0]);
  c1 = rawCoords ? c1 : getDrawCoords(c1);
  d_os << "d='M" << c1.x << "," << c1.y;
  for (unsigned int i = 0; i < nSegments; ++i) {
    auto cpt1 = std::get<1>(segments[i]);
    cpt1 = rawCoords ? cpt1 : getDrawCoords(cpt1);
    auto cpt2 = std::get<2>(segments[i]);
    cpt2 = rawCoords ? cpt2 : getDrawCoords(cpt2);
    auto segpt = std::get<3>(segments[i]);
    segpt = rawCoords ? segpt : getDrawCoords(segpt);
    d_os << " C" << cpt1.x << "," << cpt1.y << " " << cpt2.x << "," << cpt2.y
         << " " << segpt.x << "," << segpt.y;
  }
  d_os << "' ";

  d_os << "style='fill:none;stroke:" << col
       << ";stroke-width:" << boost::format("%.1f") % width
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       << "'";
  d_os << " />\n";
}

namespace {
std::string getDashString(const DashPattern &dashes) {
  std::string res;
  if (dashes.size()) {
    std::stringstream dss;
    dss << ";stroke-dasharray:";
    std::copy(dashes.begin(), dashes.end() - 1,
              std::ostream_iterator<double>(dss, ","));
    dss << dashes.back();
    res = dss.str();
  }
  return res;
}
}  // namespace

// ****************************************************************************
void MolDraw2DSVG::drawLine(const Point2D &cds1, const Point2D &cds2,
                            bool rawCoords) {
  Point2D c1 = rawCoords ? cds1 : getDrawCoords(cds1);
  Point2D c2 = rawCoords ? cds2 : getDrawCoords(cds2);
  std::string col = DrawColourToSVG(colour());
  double width = getDrawLineWidth();
  std::string dashString = getDashString(dash());
  d_os << "<path ";
  outputClasses();
  d_os << "d='M " << MolDraw2D_detail::formatDouble(c1.x) << ","
       << MolDraw2D_detail::formatDouble(c1.y) << " L "
       << MolDraw2D_detail::formatDouble(c2.x) << ","
       << MolDraw2D_detail::formatDouble(c2.y) << "' ";
  d_os << "style='fill:none;fill-rule:evenodd;stroke:" << col
       << ";stroke-width:" << MolDraw2D_detail::formatDouble(width)
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       << dashString << "'";
  d_os << " />\n";
}

// ****************************************************************************
void MolDraw2DSVG::drawPolygon(const std::vector<Point2D> &cds,
                               bool rawCoords) {
  PRECONDITION(cds.size() >= 3, "must have at least three points");

  std::string col = DrawColourToSVG(colour());
  double width = getDrawLineWidth();
  std::string dashString = getDashString(dash());

  d_os << "<path ";
  outputClasses();
  d_os << "d='M";
  Point2D c0 = rawCoords ? cds[0] : getDrawCoords(cds[0]);
  d_os << " " << MolDraw2D_detail::formatDouble(c0.x) << ","
       << MolDraw2D_detail::formatDouble(c0.y);
  for (unsigned int i = 1; i < cds.size(); ++i) {
    Point2D ci = rawCoords ? cds[i] : getDrawCoords(cds[i]);
    d_os << " L " << MolDraw2D_detail::formatDouble(ci.x) << ","
         << MolDraw2D_detail::formatDouble(ci.y);
  }
  if (fillPolys()) {
    // the Z closes the path which we don't want for unfilled polygons
    d_os << " Z' style='fill:" << col
         << ";fill-rule:evenodd;fill-opacity:" << colour().a << ";";
  } else {
    d_os << "' style='fill:none;";
  }

  d_os << "stroke:" << col
       << ";stroke-width:" << MolDraw2D_detail::formatDouble(width)
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:"
       << colour().a << ";" << dashString << "'";
  d_os << " />\n";
}

// ****************************************************************************
void MolDraw2DSVG::drawEllipse(const Point2D &cds1, const Point2D &cds2,
                               bool rawCoords) {
  Point2D c1 = rawCoords ? cds1 : getDrawCoords(cds1);
  Point2D c2 = rawCoords ? cds2 : getDrawCoords(cds2);
  double w = c2.x - c1.x;
  double h = c2.y - c1.y;
  double cx = c1.x + w / 2;
  double cy = c1.y + h / 2;
  w = w > 0 ? w : -1 * w;
  h = h > 0 ? h : -1 * h;

  std::string col = DrawColourToSVG(colour());
  double width = getDrawLineWidth();
  std::string dashString = getDashString(dash());
  d_os << "<ellipse"
       << " cx='" << MolDraw2D_detail::formatDouble(cx) << "'"
       << " cy='" << MolDraw2D_detail::formatDouble(cy) << "'"
       << " rx='" << MolDraw2D_detail::formatDouble(w / 2) << "'"
       << " ry='" << MolDraw2D_detail::formatDouble(h / 2) << "' ";
  outputClasses();
  d_os << " style='";
  if (fillPolys()) {
    d_os << "fill:" << col << ";fill-rule:evenodd;";
  } else {
    d_os << "fill:none;";
  }

  d_os << "stroke:" << col << ";stroke-width:" << boost::format("%.1f") % width
       << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       << dashString << "'";
  d_os << " />\n";
}

// ****************************************************************************
void MolDraw2DSVG::clearDrawing() {
  MolDraw2D::clearDrawing();

  std::string col = DrawColourToSVG(drawOptions().backgroundColour);
  d_os << "<rect";
  d_os << " style='opacity:1.0;fill:" << col << ";stroke:none'";
  d_os << " width='" << MolDraw2D_detail::formatDouble(width()) << "' height='"
       << MolDraw2D_detail::formatDouble(height()) << "'";
  d_os << " x='" << MolDraw2D_detail::formatDouble(offset().x) << "' y='"
       << MolDraw2D_detail::formatDouble(offset().y) << "'";
  d_os << "> </rect>\n";
}

// ****************************************************************************
static const char *RDKIT_SVG_VERSION = "0.9";
void MolDraw2DSVG::addMoleculeMetadata(const ROMol &mol, int confId) const {
  PRECONDITION(d_os, "no output stream");
  d_os << "<metadata>"
       << "\n";
  d_os << "<rdkit:mol"
       << " xmlns:rdkit = \"http://www.rdkit.org/xml\""
       << " version=\"" << RDKIT_SVG_VERSION << "\""
       << ">"
       << "\n";
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
    d_os << " drawing-x=\"" << MolDraw2D_detail::formatDouble(dpos.x) << "\""
         << " drawing-y=\"" << MolDraw2D_detail::formatDouble(dpos.y) << "\"";
    d_os << " x=\"" << pos.x << "\""
         << " y=\"" << pos.y << "\""
         << " z=\"" << pos.z << "\"";

    outputMetaData(atom, d_os);

    d_os << " />"
         << "\n";
  }
  for (const auto bond : mol.bonds()) {
    d_os << "<rdkit:bond idx=\"" << bond->getIdx() + 1 << "\"";
    d_os << " begin-atom-idx=\"" << bond->getBeginAtomIdx() + 1 << "\"";
    d_os << " end-atom-idx=\"" << bond->getEndAtomIdx() + 1 << "\"";
    bool doKekule = false, allBondsExplicit = true;
    d_os << " bond-smiles=\""
         << SmilesWrite::GetBondSmiles(bond, -1, doKekule, allBondsExplicit)
         << "\"";

    outputMetaData(bond, d_os);

    d_os << " />"
         << "\n";
  }
  d_os << "</rdkit:mol></metadata>"
       << "\n";
}

// ****************************************************************************
void MolDraw2DSVG::addMoleculeMetadata(const std::vector<ROMol *> &mols,
                                       const std::vector<int> confIds) {
  for (unsigned int i = 0; i < mols.size(); ++i) {
    int confId = -1;
    if (confIds.size() == mols.size()) {
      confId = confIds[i];
    }
    setActiveMolIdx(i);
    addMoleculeMetadata(*(mols[i]), confId);
  }
};

// ****************************************************************************
void MolDraw2DSVG::tagAtoms(const ROMol &mol, double radius,
                            const std::map<std::string, std::string> &events) {
  PRECONDITION(d_os, "no output stream");
  // first bonds so that they are under the atoms
  for (const auto &bond : mol.bonds()) {
    const auto this_idx = bond->getIdx();
    const auto a_idx1 = bond->getBeginAtomIdx();
    const auto a_idx2 = bond->getEndAtomIdx();
    const auto a1pos = getDrawCoords(bond->getBeginAtomIdx());
    const auto a2pos = getDrawCoords(bond->getEndAtomIdx());
    const auto width = 2 + lineWidth();
    if (drawOptions().splitBonds) {
      const auto midp = (a1pos + a2pos) / 2;
      // from begin to mid
      d_os << "<path "
           << " d='M " << MolDraw2D_detail::formatDouble(a1pos.x) << ","
           << MolDraw2D_detail::formatDouble(a1pos.y) << " L "
           << MolDraw2D_detail::formatDouble(midp.x) << ","
           << MolDraw2D_detail::formatDouble(midp.y) << "'";
      d_os << " class='bond-selector bond-" << this_idx << " atom-" << a_idx1;
      outputTagClasses(bond, d_os, d_activeClass);
      d_os << "'";
      d_os << " style='fill:#fff;stroke:#fff;stroke-width:"
           << MolDraw2D_detail::formatDouble(width)
           << "px;fill-opacity:0;"
              "stroke-opacity:0' ";
      d_os << "/>\n";
      // mid to end
      d_os << "<path "
           << " d='M " << MolDraw2D_detail::formatDouble(midp.x) << ","
           << MolDraw2D_detail::formatDouble(midp.y) << " L "
           << MolDraw2D_detail::formatDouble(a2pos.x) << ","
           << MolDraw2D_detail::formatDouble(a2pos.y) << "'";
      d_os << " class='bond-selector bond-" << this_idx << " atom-" << a_idx2;
      outputTagClasses(bond, d_os, d_activeClass);
      d_os << "'";
      d_os << " style='fill:#fff;stroke:#fff;stroke-width:"
           << MolDraw2D_detail::formatDouble(width)
           << "px;fill-opacity:0;"
              "stroke-opacity:0' ";
      d_os << "/>\n";
    } else {
      d_os << "<path "
           << " d='M " << MolDraw2D_detail::formatDouble(a1pos.x) << ","
           << MolDraw2D_detail::formatDouble(a1pos.y) << " L "
           << MolDraw2D_detail::formatDouble(a2pos.x) << ","
           << MolDraw2D_detail::formatDouble(a2pos.y) << "'";
      d_os << " class='bond-selector bond-" << this_idx << " atom-" << a_idx1
           << " atom-" << a_idx2;
      outputTagClasses(bond, d_os, d_activeClass);
      d_os << "'";
      d_os << " style='fill:#fff;stroke:#fff;stroke-width:"
           << MolDraw2D_detail::formatDouble(width)
           << "px;fill-opacity:0;"
              "stroke-opacity:0' ";
      d_os << "/>\n";
    }
  }
  for (const auto &at : mol.atoms()) {
    auto this_idx = at->getIdx();
    auto pos = getDrawCoords(this_idx);
    d_os << "<circle "
         << " cx='" << MolDraw2D_detail::formatDouble(pos.x) << "'"
         << " cy='" << MolDraw2D_detail::formatDouble(pos.y) << "'"
         << " r='" << MolDraw2D_detail::formatDouble(scale() * radius) << "'";
    d_os << " class='atom-selector atom-" << this_idx;
    outputTagClasses(at, d_os, d_activeClass);
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

// ****************************************************************************
void MolDraw2DSVG::outputClasses() {
  if (d_activeClass.empty() && !hasActiveAtmIdx()) {
    return;
  }

  d_os << "class='";
  if (!d_activeClass.empty()) {
    d_os << d_activeClass;
  }
  if (hasActiveBndIdx()) {
    d_os << "bond-" << getActiveBndIdx();
  }

  if (!hasActiveAtmIdx()) {
    d_os << "' ";
    return;
  }
  if (hasActiveBndIdx()) {
    d_os << " ";
  }
  d_os << (!d_activeClass.empty() ? " " : "");
  const auto aidx1 = getActiveAtmIdx1();
  if (aidx1 >= 0) {
    d_os << "atom-" << aidx1;
  }
  const auto aidx2 = getActiveAtmIdx2();
  if (aidx2 >= 0 && aidx2 != aidx1) {
    d_os << " atom-" << aidx2;
  }
  d_os << "' ";
}

}  // namespace RDKit
