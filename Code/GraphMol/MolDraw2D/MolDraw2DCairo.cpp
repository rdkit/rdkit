//
//  Copyright (C) 2015-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// derived from Dave Cosgrove's MolDraw2D
//

#include <cairo.h>
#include <fstream>
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#ifdef RDK_BUILD_FREETYPE_SUPPORT
#include <GraphMol/MolDraw2D/DrawTextFTCairo.h>
#else
#include <GraphMol/MolDraw2D/DrawTextCairo.h>
#endif
#include <GraphMol/FileParsers/PNGParser.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>

#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <boost/format.hpp>

namespace RDKit {
void MolDraw2DCairo::initDrawing() {
  PRECONDITION(dp_cr, "no draw context");
  cairo_set_line_cap(dp_cr, CAIRO_LINE_CAP_BUTT);
  //  drawOptions().backgroundColour = DrawColour(0.9, 0.9, 0.0);
}

void MolDraw2DCairo::initTextDrawer(bool noFreetype) {
  double max_fnt_sz = drawOptions().maxFontSize;
  double min_fnt_sz = drawOptions().minFontSize;

  if (noFreetype) {
    text_drawer_.reset(new DrawTextCairo(max_fnt_sz, min_fnt_sz, dp_cr));
  } else {
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    try {
      text_drawer_.reset(new DrawTextFTCairo(max_fnt_sz, min_fnt_sz,
                                             drawOptions().fontFile, dp_cr));
    } catch (std::runtime_error &e) {
      BOOST_LOG(rdWarningLog)
          << e.what() << std::endl
          << "Falling back to native Cairo text handling." << std::endl;
      text_drawer_.reset(new DrawTextCairo(max_fnt_sz, min_fnt_sz, dp_cr));
    }
#else
    text_drawer_.reset(new DrawTextCairo(max_fnt_sz, min_fnt_sz, dp_cr));
#endif
  }
}

// ****************************************************************************
void MolDraw2DCairo::finishDrawing() {}

// ****************************************************************************
void MolDraw2DCairo::setColour(const DrawColour &col) {
  PRECONDITION(dp_cr, "no draw context");
  MolDraw2D::setColour(col);
  cairo_set_source_rgb(dp_cr, col.r, col.g, col.b);
}

// ****************************************************************************
void MolDraw2DCairo::drawLine(const Point2D &cds1, const Point2D &cds2) {
  PRECONDITION(dp_cr, "no draw context");
  Point2D c1 = getDrawCoords(cds1);
  Point2D c2 = getDrawCoords(cds2);

  double width = getDrawLineWidth();
  std::string dashString = "";

  cairo_set_line_width(dp_cr, width);

  const DashPattern &dashes = dash();
  if (dashes.size()) {
    auto *dd = new double[dashes.size()];
    std::copy(dashes.begin(), dashes.end(), dd);
    cairo_set_dash(dp_cr, dd, dashes.size(), 0);
    delete[] dd;
  } else {
    cairo_set_dash(dp_cr, nullptr, 0, 0);
  }

  cairo_move_to(dp_cr, c1.x, c1.y);
  cairo_line_to(dp_cr, c2.x, c2.y);
  cairo_stroke(dp_cr);
}

void MolDraw2DCairo::drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                                  const DrawColour &col1,
                                  const DrawColour &col2,
                                  unsigned int nSegments, double vertOffset) {
  PRECONDITION(dp_cr, "no draw context");
  PRECONDITION(nSegments > 1, "too few segments");
  RDUNUSED_PARAM(col2);

  if (nSegments % 2) {
    ++nSegments;  // we're going to assume an even number of segments
  }

  Point2D perp = calcPerpendicular(cds1, cds2);
  Point2D delta = (cds2 - cds1);
  perp *= vertOffset;
  delta /= nSegments;

  Point2D c1 = getDrawCoords(cds1);

  double width = getDrawLineWidth();
  cairo_set_line_width(dp_cr, width);
  cairo_set_dash(dp_cr, nullptr, 0, 0);
  setColour(col1);
  cairo_move_to(dp_cr, c1.x, c1.y);
  for (unsigned int i = 0; i < nSegments; ++i) {
    Point2D startpt = cds1 + delta * i;
    Point2D segpt = getDrawCoords(startpt + delta);
    Point2D cpt1 =
        getDrawCoords(startpt + delta / 3. + perp * (i % 2 ? -1 : 1));
    Point2D cpt2 =
        getDrawCoords(startpt + delta * 2. / 3. + perp * (i % 2 ? -1 : 1));
    // if (i == nSegments / 2 && col2 != col1) setColour(col2);
    cairo_curve_to(dp_cr, cpt1.x, cpt1.y, cpt2.x, cpt2.y, segpt.x, segpt.y);
  }
  cairo_stroke(dp_cr);
}

// ****************************************************************************
void MolDraw2DCairo::drawPolygon(const std::vector<Point2D> &cds) {
  PRECONDITION(dp_cr, "no draw context");
  PRECONDITION(cds.size() >= 3, "must have at least three points");

  unsigned int width = getDrawLineWidth();

  cairo_line_cap_t olinecap = cairo_get_line_cap(dp_cr);
  cairo_line_join_t olinejoin = cairo_get_line_join(dp_cr);

  cairo_set_line_cap(dp_cr, CAIRO_LINE_CAP_BUTT);
  cairo_set_line_join(dp_cr, CAIRO_LINE_JOIN_BEVEL);
  cairo_set_dash(dp_cr, nullptr, 0, 0);
  cairo_set_line_width(dp_cr, width);

  if (!cds.empty()) {
    Point2D lc = getDrawCoords(cds.front());
    cairo_move_to(dp_cr, lc.x, lc.y);
    for (unsigned int i = 1; i < cds.size(); ++i) {
      Point2D lci = getDrawCoords(cds[i]);
      cairo_line_to(dp_cr, lci.x, lci.y);
    }
  }

  if (fillPolys()) {
    cairo_close_path(dp_cr);
    cairo_fill_preserve(dp_cr);
  }
  cairo_stroke(dp_cr);
  cairo_set_line_cap(dp_cr, olinecap);
  cairo_set_line_join(dp_cr, olinejoin);
}

// ****************************************************************************
void MolDraw2DCairo::clearDrawing() {
  PRECONDITION(dp_cr, "no draw context");
  setColour(drawOptions().backgroundColour);
  cairo_rectangle(dp_cr, offset().x, offset().y, width(), height());
  cairo_fill(dp_cr);
}

// ****************************************************************************
namespace {
cairo_status_t grab_str(void *closure, const unsigned char *data,
                        unsigned int len) {
  auto *str_ptr = (std::string *)closure;
  (*str_ptr) += std::string((const char *)data, len);
  return CAIRO_STATUS_SUCCESS;
}
}  // namespace
std::string MolDraw2DCairo::getDrawingText() const {
  PRECONDITION(dp_cr, "no draw context");
  std::string res = "";
  cairo_surface_t *surf = cairo_get_target(dp_cr);
  cairo_surface_write_to_png_stream(surf, &grab_str, (void *)&res);
  res = addMetadataToPNG(res);
  return res;
};

void MolDraw2DCairo::writeDrawingText(const std::string &fName) const {
  PRECONDITION(dp_cr, "no draw context");
  auto png = getDrawingText();
  std::ofstream outs(fName.c_str(), std::ios_base::binary | std::ios_base::out);
  if (!outs || outs.bad()) {
    BOOST_LOG(rdErrorLog) << "Failed to write PNG file " << fName << std::endl;
    return;
  }
  outs.write(png.c_str(), png.size());
};

std::string MolDraw2DCairo::addMetadataToPNG(const std::string &png) const {
  if (d_metadata.empty()) {
    return png;
  }

  return RDKit::addMetadataToPNGString(png, d_metadata);
}

// ****************************************************************************
namespace {
void addMoleculeMetadata(
    const ROMol &mol, int confId,
    std::vector<std::pair<std::string, std::string>> &metadata, unsigned idx) {
  // it's legal to repeat a tag, but a lot of software doesn't show all values
  // so we append a suffix to disambiguate when necessary
  std::string suffix = "";
  if (idx) {
    suffix = (boost::format("%d") % idx).str();
  }
  std::string pkl;
  MolPickler::pickleMol(mol, pkl);
  metadata.push_back(
      std::make_pair(augmentTagName(PNGData::pklTag + suffix), pkl));

  bool includeStereo = true;
  if (mol.getNumConformers()) {
    bool kekulize = false;
    auto molb = MolToMolBlock(mol, includeStereo, confId, kekulize);
    metadata.push_back(
        std::make_pair(augmentTagName(PNGData::molTag + suffix), molb));
  }
  // MolToCXSmiles() is missing the feature that lets us specify confIds
  auto smiles = MolToCXSmiles(mol);
  metadata.push_back(
      std::make_pair(augmentTagName(PNGData::smilesTag + suffix), smiles));
}
void addReactionMetadata(
    const ChemicalReaction &rxn,
    std::vector<std::pair<std::string, std::string>> &metadata) {
  std::string pkl;
  ReactionPickler::pickleReaction(rxn, pkl);
  metadata.push_back(std::make_pair(augmentTagName(PNGData::rxnPklTag), pkl));
  metadata.push_back(std::make_pair(augmentTagName(PNGData::rxnSmartsTag),
                                    ChemicalReactionToRxnSmarts(rxn)));
}
}  // namespace

void MolDraw2DCairo::updateMetadata(const ROMol &mol, int confId) {
  addMoleculeMetadata(mol, confId, d_metadata, d_numMetadataEntries);
  ++d_numMetadataEntries;
}
void MolDraw2DCairo::updateMetadata(const ChemicalReaction &rxn) {
  addReactionMetadata(rxn, d_metadata);
}

}  // namespace RDKit
