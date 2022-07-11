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
#include <GraphMol/MolDraw2D/DrawTextCairo.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#ifdef RDK_BUILD_FREETYPE_SUPPORT
#include <GraphMol/MolDraw2D/DrawTextFTCairo.h>
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
  if (dp_cr) {
    if (cairo_get_reference_count(dp_cr) > 0) {
      cairo_destroy(dp_cr);
    }
    dp_cr = nullptr;
  }
  cairo_surface_t *surf =
      cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width(), height());
  dp_cr = cairo_create(surf);
  cairo_surface_destroy(surf);  // dp_cr has a reference to this now;
  cairo_set_line_cap(dp_cr, CAIRO_LINE_CAP_BUTT);
  if (!text_drawer_) {
    initTextDrawer(df_noFreetype);
  } else {
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    if (df_noFreetype) {
      dynamic_cast<MolDraw2D_detail::DrawTextCairo *>(text_drawer_.get())
          ->setCairoContext(dp_cr);
    } else {
      MolDraw2D_detail::DrawTextFTCairo *dt =
          dynamic_cast<MolDraw2D_detail::DrawTextFTCairo *>(text_drawer_.get());
      dt->setCairoContext(dp_cr);
    }
#else
    dynamic_cast<MolDraw2D_detail::DrawTextCairo *>(text_drawer_.get())
        ->setCairoContext(dp_cr);
#endif
  }
}

void MolDraw2DCairo::initTextDrawer(bool noFreetype) {
  double max_fnt_sz = drawOptions().maxFontSize;
  double min_fnt_sz = drawOptions().minFontSize;

  if (noFreetype) {
    text_drawer_.reset(
        new MolDraw2D_detail::DrawTextCairo(max_fnt_sz, min_fnt_sz, dp_cr));
  } else {
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    try {
      text_drawer_.reset(new MolDraw2D_detail::DrawTextFTCairo(
          max_fnt_sz, min_fnt_sz, drawOptions().fontFile, dp_cr));
    } catch (std::runtime_error &e) {
      BOOST_LOG(rdWarningLog)
          << e.what() << std::endl
          << "Falling back to native Cairo text handling." << std::endl;
      text_drawer_.reset(
          new MolDraw2D_detail::DrawTextCairo(max_fnt_sz, min_fnt_sz, dp_cr));
    }
#else
    text_drawer_.reset(
        new MolDraw2D_detail::DrawTextCairo(max_fnt_sz, min_fnt_sz, dp_cr));
#endif
  }
  if (drawOptions().baseFontSize > 0.0) {
    text_drawer_->setBaseFontSize(drawOptions().baseFontSize);
  }
}

// ****************************************************************************
void MolDraw2DCairo::finishDrawing() {}

// ****************************************************************************
void MolDraw2DCairo::setColour(const DrawColour &col) {
  PRECONDITION(dp_cr, "no draw context");
  MolDraw2D::setColour(col);
  cairo_set_source_rgba(dp_cr, col.r, col.g, col.b, col.a);
}

// ****************************************************************************
void MolDraw2DCairo::drawLine(const Point2D &cds1, const Point2D &cds2,
                              bool rawCoords) {
  PRECONDITION(dp_cr, "no draw context");
  Point2D c1 = rawCoords ? cds1 : getDrawCoords(cds1);
  Point2D c2 = rawCoords ? cds2 : getDrawCoords(cds2);

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

// ****************************************************************************
void MolDraw2DCairo::drawWavyLine(const Point2D &cds1, const Point2D &cds2,
                                  const DrawColour &col1, const DrawColour &,
                                  unsigned int nSegments, double vertOffset,
                                  bool rawCoords) {
  PRECONDITION(dp_cr, "no draw context");
  PRECONDITION(nSegments > 1, "too few segments");

  auto segments =
      MolDraw2D_detail::getWavyLineSegments(cds1, cds2, nSegments, vertOffset);

  double width = getDrawLineWidth();
  cairo_set_line_width(dp_cr, width);
  cairo_set_dash(dp_cr, nullptr, 0, 0);
  setColour(col1);

  auto c1 = std::get<0>(segments[0]);
  c1 = rawCoords ? c1 : getDrawCoords(c1);

  cairo_move_to(dp_cr, c1.x, c1.y);
  for (unsigned int i = 0; i < nSegments; ++i) {
    auto cpt1 = std::get<1>(segments[i]);
    cpt1 = rawCoords ? cpt1 : getDrawCoords(cpt1);
    auto cpt2 = std::get<2>(segments[i]);
    cpt2 = rawCoords ? cpt2 : getDrawCoords(cpt2);
    auto segpt = std::get<3>(segments[i]);
    segpt = rawCoords ? segpt : getDrawCoords(segpt);
    cairo_curve_to(dp_cr, cpt1.x, cpt1.y, cpt2.x, cpt2.y, segpt.x, segpt.y);
  }
  cairo_stroke(dp_cr);
}

// ****************************************************************************
void MolDraw2DCairo::drawPolygon(const std::vector<Point2D> &cds,
                                 bool rawCoords) {
  PRECONDITION(dp_cr, "no draw context");
  PRECONDITION(cds.size() >= 3, "must have at least three points");

  double width = getDrawLineWidth();

  cairo_line_cap_t olinecap = cairo_get_line_cap(dp_cr);
  cairo_line_join_t olinejoin = cairo_get_line_join(dp_cr);

  cairo_set_line_cap(dp_cr, CAIRO_LINE_CAP_BUTT);
  cairo_set_line_join(dp_cr, CAIRO_LINE_JOIN_BEVEL);
  cairo_set_dash(dp_cr, nullptr, 0, 0);
  cairo_set_line_width(dp_cr, width);

  if (!cds.empty()) {
    Point2D lc = rawCoords ? cds.front() : getDrawCoords(cds.front());
    cairo_move_to(dp_cr, lc.x, lc.y);
    for (unsigned int i = 1; i < cds.size(); ++i) {
      Point2D lci = rawCoords ? cds[i] : getDrawCoords(cds[i]);
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
  MolDraw2D::clearDrawing();
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
    auto molb = MolToV3KMolBlock(mol, includeStereo, confId, kekulize);
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
