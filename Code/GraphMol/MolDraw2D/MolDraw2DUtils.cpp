//
//  Copyright (C) 2016-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Chirality.h>
#ifdef RDK_BUILD_CAIRO_SUPPORT
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#endif
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <limits>
#include <cmath>
#include <sys/stat.h>
#include <Numerics/Conrec.h>

namespace RDKit {
namespace MolDraw2DUtils {

namespace {
bool isAtomCandForChiralH(const RWMol &mol, const Atom *atom) {
  // conditions for needing a chiral H:
  //   - stereochem specified
  //   - in at least two rings
  return mol.getRingInfo()->isInitialized() &&
         mol.getRingInfo()->numAtomRings(atom->getIdx()) > 1u &&
         (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW ||
          atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW);
}
}  // end of anonymous namespace

void prepareMolForDrawing(RWMol &mol, bool kekulize, bool addChiralHs,
                          bool wedgeBonds, bool forceCoords, bool wavyBonds) {
  if (kekulize) {
    RDLog::LogStateSetter blocker;
    MolOps::KekulizeIfPossible(
        mol, false);  // kekulize, but keep the aromatic flags!
  }
  if (addChiralHs) {
    std::vector<unsigned int> chiralAts;
    for (auto atom : mol.atoms()) {
      if (isAtomCandForChiralH(mol, atom)) {
        chiralAts.push_back(atom->getIdx());
      }
    }
    if (chiralAts.size()) {
      bool addCoords = false;
      if (!forceCoords && mol.getNumConformers()) {
        addCoords = true;
      }
      MolOps::addHs(mol, false, addCoords, &chiralAts);
    }
  }
  if (forceCoords || !mol.getNumConformers()) {
    // compute 2D coordinates in a standard orientation:
    const bool canonOrient = true;
    RDDepict::compute2DCoords(mol, nullptr, canonOrient);
  }
  if (wedgeBonds) {
    WedgeMolBonds(mol, &mol.getConformer());
  }
  if (wavyBonds) {
    addWavyBondsForStereoAny(mol);
  }
}

void prepareAndDrawMolecule(MolDraw2D &drawer, const ROMol &mol,
                            const std::string &legend,
                            const std::vector<int> *highlight_atoms,
                            const std::vector<int> *highlight_bonds,
                            const std::map<int, DrawColour> *highlight_atom_map,
                            const std::map<int, DrawColour> *highlight_bond_map,
                            const std::map<int, double> *highlight_radii,
                            int confId, bool kekulize) {
  RWMol cpy(mol);
  prepareMolForDrawing(cpy, kekulize);
  // having done the prepare, we don't want to do it again in drawMolecule.
  bool old_prep_mol = drawer.drawOptions().prepareMolsBeforeDrawing;
  drawer.drawOptions().prepareMolsBeforeDrawing = false;
  drawer.drawMolecule(cpy, legend, highlight_atoms, highlight_bonds,
                      highlight_atom_map, highlight_bond_map, highlight_radii,
                      confId);
  drawer.drawOptions().prepareMolsBeforeDrawing = old_prep_mol;
}

void updateMolDrawOptionsFromJSON(MolDrawOptions &opts, const char *json) {
  PRECONDITION(json, "no parameter string");
  updateMolDrawOptionsFromJSON(opts, std::string(json));
};

RDKIT_MOLDRAW2D_EXPORT void updateDrawerParamsFromJSON(MolDraw2D &drawer,
                                                       const char *json) {
  updateMolDrawOptionsFromJSON(drawer.drawOptions(), json);
}

#define PT_OPT_GET(opt) opts.opt = pt.get(#opt, opts.opt)

void get_rgba(const boost::property_tree::ptree &node, DrawColour &colour) {
  boost::property_tree::ptree::const_iterator itm = node.begin();
  colour.r = itm->second.get_value<float>();
  ++itm;
  colour.g = itm->second.get_value<float>();
  ++itm;
  colour.b = itm->second.get_value<float>();
  ++itm;
  if (itm != node.end()) {
    colour.a = itm->second.get_value<float>();
    ++itm;
  }
}

void get_colour_option(boost::property_tree::ptree *pt, const char *pnm,
                       DrawColour &colour) {
  PRECONDITION(pnm && strlen(pnm), "bad property name");
  if (pt->find(pnm) == pt->not_found()) {
    return;
  }

  const auto &node = pt->get_child(pnm);
  get_rgba(node, colour);
}

void get_colour_palette_option(boost::property_tree::ptree *pt, const char *pnm,
                               ColourPalette &palette) {
  PRECONDITION(pnm && strlen(pnm), "bad property name");
  if (pt->find(pnm) == pt->not_found()) {
    return;
  }

  for (const auto &atomicNumNodeIt : pt->get_child(pnm)) {
    int atomicNum = boost::lexical_cast<int>(atomicNumNodeIt.first);
    DrawColour colour;
    get_rgba(atomicNumNodeIt.second, colour);
    palette[atomicNum] = colour;
  }
}

void updateMolDrawOptionsFromJSON(MolDrawOptions &opts,
                                  const std::string &json) {
  if (json == "") {
    return;
  }
  std::istringstream ss;
  ss.str(json);
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  PT_OPT_GET(atomLabelDeuteriumTritium);
  PT_OPT_GET(dummiesAreAttachments);
  PT_OPT_GET(circleAtoms);
  PT_OPT_GET(splitBonds);
  PT_OPT_GET(continuousHighlight);
  PT_OPT_GET(fillHighlights);
  PT_OPT_GET(highlightRadius);
  PT_OPT_GET(flagCloseContactsDist);
  PT_OPT_GET(includeAtomTags);
  PT_OPT_GET(clearBackground);
  PT_OPT_GET(legendFontSize);
  PT_OPT_GET(legendFraction);
  PT_OPT_GET(maxFontSize);
  PT_OPT_GET(minFontSize);
  PT_OPT_GET(fixedFontSize);
  PT_OPT_GET(baseFontSize);
  PT_OPT_GET(annotationFontScale);
  PT_OPT_GET(fontFile);
  PT_OPT_GET(multipleBondOffset);
  PT_OPT_GET(padding);
  PT_OPT_GET(additionalAtomLabelPadding);
  PT_OPT_GET(noAtomLabels);
  PT_OPT_GET(bondLineWidth);
  PT_OPT_GET(scaleBondWidth);
  PT_OPT_GET(scaleHighlightBondWidth);
  PT_OPT_GET(highlightBondWidthMultiplier);
  PT_OPT_GET(prepareMolsBeforeDrawing);
  PT_OPT_GET(fixedScale);
  PT_OPT_GET(fixedBondLength);
  PT_OPT_GET(rotate);
  PT_OPT_GET(addAtomIndices);
  PT_OPT_GET(addBondIndices);
  PT_OPT_GET(isotopeLabels);
  PT_OPT_GET(dummyIsotopeLabels);
  PT_OPT_GET(addStereoAnnotation);
  PT_OPT_GET(atomHighlightsAreCircles);
  PT_OPT_GET(centreMoleculesBeforeDrawing);
  PT_OPT_GET(explicitMethyl);
  PT_OPT_GET(includeMetadata);
  PT_OPT_GET(includeRadicals);
  PT_OPT_GET(comicMode);
  PT_OPT_GET(variableBondWidthMultiplier);
  PT_OPT_GET(variableAtomRadius);
  PT_OPT_GET(includeChiralFlagLabel);
  PT_OPT_GET(simplifiedStereoGroupLabel);
  PT_OPT_GET(unspecifiedStereoIsUnknown);
  PT_OPT_GET(singleColourWedgeBonds);
  PT_OPT_GET(useMolBlockWedging);
  PT_OPT_GET(scalingFactor);
  PT_OPT_GET(drawMolsSameScale);

  get_colour_option(&pt, "highlightColour", opts.highlightColour);
  get_colour_option(&pt, "backgroundColour", opts.backgroundColour);
  get_colour_option(&pt, "legendColour", opts.legendColour);
  get_colour_option(&pt, "symbolColour", opts.symbolColour);
  get_colour_option(&pt, "annotationColour", opts.annotationColour);
  get_colour_option(&pt, "variableAttachmentColour",
                    opts.variableAttachmentColour);
  get_colour_palette_option(&pt, "atomColourPalette", opts.atomColourPalette);
  if (pt.find("atomLabels") != pt.not_found()) {
    for (const auto &item : pt.get_child("atomLabels")) {
      opts.atomLabels[boost::lexical_cast<int>(item.first)] =
          item.second.get_value<std::string>();
    }
  }
}

RDKIT_MOLDRAW2D_EXPORT void updateDrawerParamsFromJSON(
    MolDraw2D &drawer, const std::string &json) {
  updateMolDrawOptionsFromJSON(drawer.drawOptions(), json);
}

void contourAndDrawGrid(MolDraw2D &drawer, const double *grid,
                        const std::vector<double> &xcoords,
                        const std::vector<double> &ycoords, size_t nContours,
                        std::vector<double> &levels,
                        const ContourParams &params, const ROMol *mol) {
  PRECONDITION(grid, "no data");
  PRECONDITION(params.colourMap.size() > 1,
               "colourMap must have at least two entries");

  if (params.setScale) {
    Point2D minP = {xcoords[0], ycoords[0]};
    Point2D maxP = {xcoords.back(), ycoords.back()};
    drawer.setScale(drawer.width(), drawer.height(), minP, maxP, mol);
  }

  size_t nX = xcoords.size();
  size_t nY = ycoords.size();
  double minV = std::numeric_limits<double>::max();
  double maxV = std::numeric_limits<double>::lowest();
  if (!levels.size() || params.fillGrid) {
    for (size_t i = 0; i < nX; ++i) {
      for (size_t j = 0; j < nY; ++j) {
        minV = std::min(minV, grid[i * nY + j]);
        maxV = std::max(maxV, grid[i * nY + j]);
      }
    }
    if (!levels.size()) {
      levels.resize(nContours);
      for (size_t i = 0; i < nContours; ++i) {
        levels[i] = minV + i * (maxV - minV) / (nContours - 1);
      }
    }
  }
  if (maxV <= minV) {
    return;
  }

  const auto olw = drawer.lineWidth();
  const auto odash = drawer.dash();
  const auto ocolor = drawer.colour();
  const auto ofill = drawer.fillPolys();
  const auto owidth = drawer.lineWidth();
  if (params.fillGrid) {
    drawer.setFillPolys(true);
    drawer.setLineWidth(1);
    auto delta = (maxV - minV);
    if (params.colourMap.size() > 2) {
      // need to find how fractionally far we are from zero, not the min
      if (-minV > maxV) {
        delta = -minV;
      } else {
        delta = maxV;
      }
    }
    for (size_t i = 0; i < nX - 1; ++i) {
      for (size_t j = 0; j < nY - 1; ++j) {
        auto gridV = grid[i * nY + j];
        auto fracV = (gridV - minV) / delta;
        if (params.colourMap.size() > 2) {
          // need to find how fractionally far we are from zero, not the min
          fracV = gridV / delta;
          if (fracV < 0) {
            fracV *= -1;
          }
        }
        auto c1 = (gridV < 0 || params.colourMap.size() == 2)
                      ? params.colourMap[1]
                      : params.colourMap[1];
        auto c2 = (gridV < 0 || params.colourMap.size() == 2)
                      ? params.colourMap[0]
                      : params.colourMap[2];
        auto c = c1 + (c2 - c1) * fracV;
        // don't bother drawing boxes that are the same as the background color:
        double tol = 0.01;
        if (c.feq(drawer.drawOptions().backgroundColour, tol)) {
          continue;
        }
        drawer.setColour(c);
        Point2D p1 = {xcoords[i], ycoords[j]};
        Point2D p2 = {xcoords[i + 1], ycoords[j + 1]};
        drawer.drawRect(p1, p2);
      }
    }
  }

  if (nContours) {
    if (nContours > levels.size()) {
      throw ValueErrorException(
          "nContours larger than the size of the level list");
    }
    std::vector<conrec::ConrecSegment> segs;
    conrec::Contour(grid, 0, nX - 1, 0, nY - 1, xcoords.data(), ycoords.data(),
                    nContours, levels.data(), segs);
    static DashPattern negDash = {2, 6};
    static DashPattern posDash;
    drawer.setColour(params.contourColour);
    drawer.setLineWidth(params.contourWidth);
    for (const auto &seg : segs) {
      if (params.dashNegative && seg.isoVal < 0) {
        drawer.setDash(negDash);
      } else {
        drawer.setDash(posDash);
      }
      drawer.drawLine(seg.p1, seg.p2);
    }
  }

  drawer.setDash(odash);
  drawer.setLineWidth(olw);
  drawer.setColour(ocolor);
  drawer.setFillPolys(ofill);
  drawer.setLineWidth(owidth);
};

void contourAndDrawGaussians(MolDraw2D &drawer,
                             const std::vector<Point2D> &locs,
                             const std::vector<double> &weights,
                             const std::vector<double> &widths,
                             size_t nContours, std::vector<double> &levels,
                             const ContourParams &params, const ROMol *mol) {
  PRECONDITION(locs.size() == weights.size(), "size mismatch");
  PRECONDITION(locs.size() == widths.size(), "size mismatch");

  // start by setting up the grid
  if (params.setScale) {
    Point2D minP, maxP;
    minP.x = minP.y = std::numeric_limits<double>::max();
    maxP.x = maxP.y = std::numeric_limits<double>::lowest();
    for (const auto &loc : locs) {
      minP.x = std::min(loc.x, minP.x);
      minP.y = std::min(loc.y, minP.y);
      maxP.x = std::max(loc.x, maxP.x);
      maxP.y = std::max(loc.y, maxP.y);
    }
    Point2D dims = maxP - minP;
    minP.x -= drawer.drawOptions().padding * dims.x;
    minP.y -= drawer.drawOptions().padding * dims.y;
    maxP.x += drawer.drawOptions().padding * dims.x;
    maxP.y += drawer.drawOptions().padding * dims.y;

    if (params.extraGridPadding > 0) {
      minP.x -= params.extraGridPadding;
      minP.y -= params.extraGridPadding;
      maxP.x += params.extraGridPadding;
      maxP.y += params.extraGridPadding;
    }
    drawer.setScale(drawer.width(), drawer.height(), minP, maxP, mol);
  }

  size_t nx = (size_t)ceil(drawer.range().x / params.gridResolution) + 1;
  size_t ny = (size_t)ceil(drawer.range().y / params.gridResolution) + 1;
  std::vector<double> xcoords(nx);
  for (size_t i = 0; i < nx; ++i) {
    xcoords[i] = drawer.minPt().x + i * params.gridResolution;
  }
  std::vector<double> ycoords(ny);
  for (size_t i = 0; i < ny; ++i) {
    ycoords[i] = drawer.minPt().y + i * params.gridResolution;
  }
  std::unique_ptr<double[]> grid(new double[nx * ny]);

  // populate the grid from the gaussians:
  for (size_t ix = 0; ix < nx; ++ix) {
    auto px = drawer.minPt().x + ix * params.gridResolution;
    for (size_t iy = 0; iy < ny; ++iy) {
      auto py = drawer.minPt().y + iy * params.gridResolution;
      Point2D pt(px, py);
      double accum = 0.0;
      for (size_t ig = 0; ig < locs.size(); ++ig) {
        auto d2 = (pt - locs[ig]).lengthSq();
        auto contrib = weights[ig] / widths[ig] *
                       exp(-0.5 * d2 / (widths[ig] * widths[ig]));
        accum += contrib;
      }
      grid[ix * ny + iy] = accum / (2 * M_PI);
    }
  }

  // and render it:
  ContourParams paramsCopy = params;
  paramsCopy.setScale = false;  // if scaling was needed, we did it already
  contourAndDrawGrid(drawer, grid.get(), xcoords, ycoords, nContours, levels,
                     paramsCopy);
};

// ****************************************************************************
void drawMolACS1996(MolDraw2D &drawer, const ROMol &mol,
                    const std::string &legend,
                    const std::vector<int> *highlight_atoms,
                    const std::vector<int> *highlight_bonds,
                    const std::map<int, DrawColour> *highlight_atom_map,
                    const std::map<int, DrawColour> *highlight_bond_map,
                    const std::map<int, double> *highlight_radii, int confId) {
  if (drawer.width() != -1 || drawer.height() != -1) {
    BOOST_LOG(rdWarningLog)
        << "ACS drawing mode works best with a flexiCanvas i.e. a drawer"
        << " created with width and height of -1.  The scale will be fixed,"
        << " and that may not look great with a pre-determined size."
        << std::endl;
  }
  double meanBondLen = MolDraw2DUtils::meanBondLength(mol, confId);
  setACS1996Options(drawer.drawOptions(), meanBondLen);
  drawer.drawMolecule(mol, legend, highlight_atoms, highlight_bonds,
                      highlight_atom_map, highlight_bond_map, highlight_radii,
                      confId);
}

// ****************************************************************************
void setACS1996Options(MolDrawOptions &opts, double meanBondLen) {
  opts.bondLineWidth = 0.6;
  opts.scaleBondWidth = false;
  // the guideline is for a bond length of 14.4px, and we set things up
  // in pixels per Angstrom.
  opts.scalingFactor = 14.4 / meanBondLen;
  // setting the fixedBondLength means the drawing won't be scaled
  // up in a drawer of defined size, so the bond length won't exceed
  // 14.4 pixels.
  opts.fixedBondLength = 14.4 / meanBondLen;
  // offset for multiple bonds is 18% of the bond length.
  opts.multipleBondOffset = 0.18;
  opts.highlightBondWidthMultiplier = 32;
  setMonochromeMode(opts, DrawColour(0.0, 0.0, 0.0), DrawColour(1.0, 1.0, 1.0));

  opts.fixedFontSize = 10;
  opts.additionalAtomLabelPadding = 0.066;
  // The guidelines say Arial font, which is not a free font.  A close
  // approximation is FreeSans, but that is under GPL v3.0, so can't be
  // embedded.  Use it if it's there, but fall back on the Roboto font
  // which uses an Apache 2.0 license and is also fairly close to Arial.
  // It is up to the user to put the FreeSans.ttf in the right place.
  // If the user has already specified a fontFile, assume they know
  // what they're doing and use it.
  if (opts.fontFile.empty()) {
    const char *rdbase = getenv("RDBASE");
    bool have_free_sans = false;
    if (rdbase) {
      opts.fontFile = std::string(rdbase) + "/Data/Fonts/FreeSans.ttf";
      struct stat buffer;
      have_free_sans = (stat(opts.fontFile.c_str(), &buffer) == 0);
    }
    if (!rdbase || !have_free_sans) {
      opts.fontFile = "BuiltinRobotoRegular";
    }
  }
}

// ****************************************************************************
double meanBondLength(const ROMol &mol, int confId) {
  double bondLen = 0.0;
  if (mol.getNumBonds()) {
    auto conf = mol.getConformer(confId);
    for (auto bond : mol.bonds()) {
      bondLen += MolTransforms::getBondLength(conf, bond->getBeginAtomIdx(),
                                              bond->getEndAtomIdx());
    }
    bondLen /= mol.getNumBonds();
  }
  return bondLen;
}

}  // namespace MolDraw2DUtils
}  // namespace RDKit
