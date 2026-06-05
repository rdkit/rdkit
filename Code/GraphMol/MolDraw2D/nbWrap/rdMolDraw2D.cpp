//
//  Copyright (C) 2015-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef USE_BETTER_ENUMS
#define USE_BETTER_ENUMS
#endif
#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>

#include <numpy/arrayobject.h>
#include <RDBoost/Wrap_nb.h>
#include <RDBoost/import_array.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <Geometry/point.h>
#ifdef RDK_BUILD_CAIRO_SUPPORT
#include <cairo.h>
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#endif

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

struct IntStringMap {
  std::map<int, std::string> *dp_map;
};

void tagAtomHelper(MolDraw2DSVG &self, const ROMol &mol, double radius,
                   nb::object pyo) {
  std::map<std::string, std::string> events;
  if (!pyo.is_none()) {
    auto tDict = nb::cast<nb::dict>(pyo);
    for (auto item : tDict) {
      events[nb::cast<std::string>(item.first)] =
          nb::cast<std::string>(item.second);
    }
  }
  self.tagAtoms(mol, radius, events);
}

void pyDictToColourMap(nb::object pyo, ColourPalette &res) {
  auto tDict = nb::cast<nb::dict>(pyo);
  for (auto item : tDict) {
    auto tpl = nb::cast<nb::tuple>(item.second);
    float r = nb::cast<float>(tpl[0]);
    float g = nb::cast<float>(tpl[1]);
    float b = nb::cast<float>(tpl[2]);
    float a = 1.0;
    if (nb::len(tpl) > 3) {
      a = nb::cast<float>(tpl[3]);
    }
    DrawColour clr(r, g, b, a);
    res[nb::cast<int>(item.first)] = clr;
  }
}

ColourPalette *pyDictToColourMap(nb::object pyo) {
  ColourPalette *res = nullptr;
  if (!pyo.is_none()) {
    res = new ColourPalette;
    pyDictToColourMap(pyo, *res);
  }
  return res;
}

void pyDictToDoubleMap(nb::object pyo, std::map<int, double> &res) {
  auto tDict = nb::cast<nb::dict>(pyo);
  for (auto item : tDict) {
    res[nb::cast<int>(item.first)] = nb::cast<double>(item.second);
  }
}

std::map<int, double> *pyDictToDoubleMap(nb::object pyo) {
  std::map<int, double> *res = nullptr;
  if (!pyo.is_none()) {
    res = new std::map<int, double>;
    pyDictToDoubleMap(pyo, *res);
  }
  return res;
}

void pyDictToIntMap(nb::object pyo, std::map<int, int> &res) {
  auto tDict = nb::cast<nb::dict>(pyo);
  for (auto item : tDict) {
    res[nb::cast<int>(item.first)] = nb::cast<int>(item.second);
  }
}

std::map<int, int> *pyDictToIntMap(nb::object pyo) {
  std::map<int, int> *res = nullptr;
  if (!pyo.is_none()) {
    res = new std::map<int, int>;
    pyDictToIntMap(pyo, *res);
  }
  return res;
}

DrawColour pyTupleToDrawColour(const nb::tuple tpl) {
  float r = nb::cast<float>(tpl[0]);
  if (r > 1 || r < 0) {
    throw ValueErrorException("RGBA color value needs to be between 0 and 1.");
  }
  float g = nb::cast<float>(tpl[1]);
  if (g > 1 || g < 0) {
    throw ValueErrorException("RGBA color value needs to be between 0 and 1.");
  }
  float b = nb::cast<float>(tpl[2]);
  if (b > 1 || b < 0) {
    throw ValueErrorException("RGBA color value needs to be between 0 and 1.");
  }
  float a = 1;
  if (nb::len(tpl) > 3) {
    a = nb::cast<float>(tpl[3]);
    if (a > 1 || a < 0) {
      throw ValueErrorException(
          "RGBA color value needs to be between 0 and 1.");
    }
  }
  return DrawColour(r, g, b, a);
}

void pyListToColourVec(nb::object pyo, std::vector<DrawColour> &res) {
  res.clear();
  auto tList = nb::cast<nb::list>(pyo);
  for (size_t i = 0; i < nb::len(tList); ++i) {
    auto tpl = nb::cast<nb::tuple>(tList[i]);
    res.push_back(pyTupleToDrawColour(tpl));
  }
}

void pyDictToMapColourVec(nb::object pyo,
                          std::map<int, std::vector<DrawColour>> &res) {
  auto tDict = nb::cast<nb::dict>(pyo);
  for (auto item : tDict) {
    auto pl = nb::cast<nb::list>(item.second);
    std::vector<DrawColour> v;
    pyListToColourVec(pl, v);
    res[nb::cast<int>(item.first)] = v;
  }
}

std::map<int, std::vector<DrawColour>> *pyDictToMapColourVec(nb::object pyo) {
  std::map<int, std::vector<DrawColour>> *res = nullptr;
  if (!pyo.is_none()) {
    res = new std::map<int, std::vector<DrawColour>>;
    pyDictToMapColourVec(pyo, *res);
  }
  return res;
}

void drawMoleculeHelper1(MolDraw2D &self, const ROMol &mol,
                         nb::object highlight_atoms,
                         nb::object highlight_atom_map,
                         nb::object highlight_atom_radii, int confId,
                         std::string legend) {
  std::unique_ptr<std::vector<int>> highlightAtoms =
      pythonObjectToVect(highlight_atoms, static_cast<int>(mol.getNumAtoms()));
  ColourPalette *ham = pyDictToColourMap(highlight_atom_map);
  std::map<int, double> *har = pyDictToDoubleMap(highlight_atom_radii);

  self.drawMolecule(mol, legend, highlightAtoms.get(), ham, har, confId);

  delete ham;
  delete har;
}

void drawMoleculeHelper2(MolDraw2D &self, const ROMol &mol,
                         nb::object highlight_atoms,
                         nb::object highlight_bonds,
                         nb::object highlight_atom_map,
                         nb::object highlight_bond_map,
                         nb::object highlight_atom_radii, int confId,
                         std::string legend) {
  std::unique_ptr<std::vector<int>> highlightAtoms =
      pythonObjectToVect(highlight_atoms, static_cast<int>(mol.getNumAtoms()));
  std::unique_ptr<std::vector<int>> highlightBonds =
      pythonObjectToVect(highlight_bonds, static_cast<int>(mol.getNumBonds()));
  ColourPalette *ham = pyDictToColourMap(highlight_atom_map);
  ColourPalette *hbm = pyDictToColourMap(highlight_bond_map);
  std::map<int, double> *har = pyDictToDoubleMap(highlight_atom_radii);

  self.drawMolecule(mol, legend, highlightAtoms.get(), highlightBonds.get(),
                    ham, hbm, har, confId);

  delete ham;
  delete hbm;
  delete har;
}

nb::tuple getMolSizeHelper(MolDraw2D &self, const ROMol &mol,
                           nb::object highlight_atoms,
                           nb::object highlight_bonds,
                           nb::object highlight_atom_map,
                           nb::object highlight_bond_map,
                           nb::object highlight_atom_radii, int confId,
                           std::string legend) {
  std::unique_ptr<std::vector<int>> highlightAtoms =
      pythonObjectToVect(highlight_atoms, static_cast<int>(mol.getNumAtoms()));
  std::unique_ptr<std::vector<int>> highlightBonds =
      pythonObjectToVect(highlight_bonds, static_cast<int>(mol.getNumBonds()));
  ColourPalette *ham = pyDictToColourMap(highlight_atom_map);
  ColourPalette *hbm = pyDictToColourMap(highlight_bond_map);
  std::map<int, double> *har = pyDictToDoubleMap(highlight_atom_radii);

  auto sz = self.getMolSize(mol, legend, highlightAtoms.get(),
                            highlightBonds.get(), ham, hbm, har, confId);

  delete ham;
  delete hbm;
  delete har;
  return nb::make_tuple(sz.first, sz.second);
}

void drawMoleculeWithHighlightsHelper(
    MolDraw2D &self, const ROMol &mol, std::string legend,
    nb::object highlight_atom_map, nb::object highlight_bond_map,
    nb::object highlight_atom_radii,
    nb::object highlight_linewidth_multipliers, int confId) {
  std::map<int, std::vector<DrawColour>> *ham =
      pyDictToMapColourVec(highlight_atom_map);
  if (!ham) {
    ham = new std::map<int, std::vector<DrawColour>>();
  }
  std::map<int, std::vector<DrawColour>> *hbm =
      pyDictToMapColourVec(highlight_bond_map);
  if (!hbm) {
    hbm = new std::map<int, std::vector<DrawColour>>();
  }
  std::map<int, double> *har = pyDictToDoubleMap(highlight_atom_radii);
  if (!har) {
    har = new std::map<int, double>();
  }
  std::map<int, int> *hlm = pyDictToIntMap(highlight_linewidth_multipliers);
  if (!hlm) {
    hlm = new std::map<int, int>();
  }
  self.drawMoleculeWithHighlights(mol, legend, *ham, *hbm, *har, *hlm, confId);

  delete ham;
  delete hbm;
  delete har;
  delete hlm;
}

void prepareAndDrawMoleculeHelper(
    MolDraw2D &drawer, const ROMol &mol, std::string legend,
    nb::object highlight_atoms, nb::object highlight_bonds,
    nb::object highlight_atom_map, nb::object highlight_bond_map,
    nb::object highlight_atom_radii, int confId, bool kekulize) {
  std::unique_ptr<std::vector<int>> highlightAtoms =
      pythonObjectToVect(highlight_atoms, static_cast<int>(mol.getNumAtoms()));
  std::unique_ptr<std::vector<int>> highlightBonds =
      pythonObjectToVect(highlight_bonds, static_cast<int>(mol.getNumBonds()));
  ColourPalette *ham = pyDictToColourMap(highlight_atom_map);
  ColourPalette *hbm = pyDictToColourMap(highlight_bond_map);
  std::map<int, double> *har = pyDictToDoubleMap(highlight_atom_radii);
  MolDraw2DUtils::prepareAndDrawMolecule(
      drawer, mol, legend, highlightAtoms.get(), highlightBonds.get(), ham, hbm,
      har, confId, kekulize);

  delete ham;
  delete hbm;
  delete har;
}

void drawMoleculeACS1996Helper(
    MolDraw2D &drawer, const ROMol &mol, std::string legend,
    nb::object highlight_atoms, nb::object highlight_bonds,
    nb::object highlight_atom_map, nb::object highlight_bond_map,
    nb::object highlight_atom_radii, int confId) {
  std::unique_ptr<std::vector<int>> highlightAtoms =
      pythonObjectToVect(highlight_atoms, static_cast<int>(mol.getNumAtoms()));
  std::unique_ptr<std::vector<int>> highlightBonds =
      pythonObjectToVect(highlight_bonds, static_cast<int>(mol.getNumBonds()));
  std::unique_ptr<ColourPalette> ham{pyDictToColourMap(highlight_atom_map)};
  std::unique_ptr<ColourPalette> hbm{pyDictToColourMap(highlight_bond_map)};
  std::unique_ptr<std::map<int, double>> har{
      pyDictToDoubleMap(highlight_atom_radii)};
  MolDraw2DUtils::drawMolACS1996(drawer, mol, legend, highlightAtoms.get(),
                                 highlightBonds.get(), ham.get(), hbm.get(),
                                 har.get(), confId);
}

void drawMoleculesHelper2(MolDraw2D &self, nb::object pmols,
                          nb::object highlight_atoms,
                          nb::object highlight_bonds,
                          nb::object highlight_atom_map,
                          nb::object highlight_bond_map,
                          nb::object highlight_atom_radii,
                          nb::object pconfIds, nb::object plegends) {
  std::unique_ptr<std::vector<ROMol *>> mols =
      pythonObjectToVect<ROMol *>(pmols);
  if (mols == nullptr || !mols->size()) {
    return;
  }
  unsigned int nThere = mols->size();
  std::unique_ptr<std::vector<std::vector<int>>> highlightAtoms;
  if (!highlight_atoms.is_none()) {
    if (nb::len(highlight_atoms) != nThere) {
      throw ValueErrorException(
          "If highlightAtoms is provided it must be the same length as the "
          "molecule list.");
    }
    highlightAtoms.reset(new std::vector<std::vector<int>>(nThere));
    for (unsigned int i = 0; i < nThere; ++i) {
      pythonObjectToVect(highlight_atoms[i], (*highlightAtoms)[i]);
    }
  }
  std::unique_ptr<std::vector<std::vector<int>>> highlightBonds;
  if (!highlight_bonds.is_none()) {
    if (nb::len(highlight_bonds) != nThere) {
      throw ValueErrorException(
          "If highlightBonds is provided it must be the same length as the "
          "molecule list.");
    }
    highlightBonds.reset(new std::vector<std::vector<int>>(nThere));
    for (unsigned int i = 0; i < nThere; ++i) {
      pythonObjectToVect(highlight_bonds[i], (*highlightBonds)[i]);
    }
  }

  std::unique_ptr<std::vector<ColourPalette>> highlightAtomMap;
  if (!highlight_atom_map.is_none()) {
    if (nb::len(highlight_atom_map) != nThere) {
      throw ValueErrorException(
          "If highlightAtomMap is provided it must be the same length as the "
          "molecule list.");
    }
    highlightAtomMap.reset(new std::vector<ColourPalette>(nThere));
    for (unsigned int i = 0; i < nThere; ++i) {
      pyDictToColourMap(highlight_atom_map[i], (*highlightAtomMap)[i]);
    }
  }
  std::unique_ptr<std::vector<ColourPalette>> highlightBondMap;
  if (!highlight_bond_map.is_none()) {
    if (nb::len(highlight_bond_map) != nThere) {
      throw ValueErrorException(
          "If highlightBondMap is provided it must be the same length as the "
          "molecule list.");
    }
    highlightBondMap.reset(new std::vector<ColourPalette>(nThere));
    for (unsigned int i = 0; i < nThere; ++i) {
      pyDictToColourMap(highlight_bond_map[i], (*highlightBondMap)[i]);
    }
  }
  std::unique_ptr<std::vector<std::map<int, double>>> highlightRadii;
  if (!highlight_atom_radii.is_none()) {
    if (nb::len(highlight_atom_radii) != nThere) {
      throw ValueErrorException(
          "If highlightAtomRadii is provided it must be the same length as the "
          "molecule list.");
    }
    highlightRadii.reset(new std::vector<std::map<int, double>>(nThere));
    for (unsigned int i = 0; i < nThere; ++i) {
      pyDictToDoubleMap(highlight_atom_radii[i], (*highlightRadii)[i]);
    }
  }
  std::unique_ptr<std::vector<int>> confIds = pythonObjectToVect<int>(pconfIds);
  std::unique_ptr<std::vector<std::string>> legends =
      pythonObjectToVect<std::string>(plegends);

  self.drawMolecules(*mols, legends.get(), highlightAtoms.get(),
                     highlightBonds.get(), highlightAtomMap.get(),
                     highlightBondMap.get(), highlightRadii.get(),
                     confIds.get());
}

void drawReactionHelper(MolDraw2D &self, const ChemicalReaction &rxn,
                        bool highlightByReactant,
                        nb::object phighlightColorsReactants,
                        nb::object pconfIds) {
  std::unique_ptr<std::vector<DrawColour>> highlightColorsReactants;
  if (!phighlightColorsReactants.is_none()) {
    highlightColorsReactants.reset(new std::vector<DrawColour>);
    pyListToColourVec(phighlightColorsReactants, *highlightColorsReactants);
  }

  std::unique_ptr<std::vector<int>> confIds = pythonObjectToVect<int>(pconfIds);

  self.drawReaction(rxn, highlightByReactant, highlightColorsReactants.get(),
                    confIds.get());
}

#ifdef RDK_BUILD_CAIRO_SUPPORT
nb::bytes getCairoDrawingText(const RDKit::MolDraw2DCairo &self) {
  std::string res = self.getDrawingText();
  return nb::bytes(res.c_str(), res.size());
}
#endif

ROMol *prepMolForDrawing(nb::object mol, bool kekulize, bool addChiralHs,
                         bool wedgeBonds, bool forceCoords, bool wavyBonds) {
  if (mol.is_none()) {
    throw std::runtime_error("molecule must not be None");
  }
  const ROMol *m = nb::cast<const ROMol *>(mol);
  auto *res = new RWMol(*m);
  MolDraw2DUtils::prepareMolForDrawing(*res, kekulize, addChiralHs, wedgeBonds,
                                       forceCoords, wavyBonds);
  return static_cast<ROMol *>(res);
}

nb::tuple colourToPyTuple(const DrawColour &clr) {
  return nb::make_tuple(clr.r, clr.g, clr.b, clr.a);
}

nb::tuple getBgColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.backgroundColour);
}
nb::tuple getQyColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.queryColour);
}
nb::tuple getHighlightColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.highlightColour);
}
void setBgColour(RDKit::MolDrawOptions &self, nb::tuple tpl) {
  self.backgroundColour = pyTupleToDrawColour(tpl);
}
void setQyColour(RDKit::MolDrawOptions &self, nb::tuple tpl) {
  self.queryColour = pyTupleToDrawColour(tpl);
}
void setHighlightColour(RDKit::MolDrawOptions &self, nb::tuple tpl) {
  self.highlightColour = pyTupleToDrawColour(tpl);
}
nb::tuple getSymbolColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.symbolColour);
}
void setSymbolColour(RDKit::MolDrawOptions &self, nb::tuple tpl) {
  self.symbolColour = pyTupleToDrawColour(tpl);
}
nb::tuple getLegendColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.legendColour);
}
void setLegendColour(RDKit::MolDrawOptions &self, nb::tuple tpl) {
  self.legendColour = pyTupleToDrawColour(tpl);
}
nb::tuple getAnnotationColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.annotationColour);
}
void setAnnotationColour(RDKit::MolDrawOptions &self, nb::tuple tpl) {
  self.annotationColour = pyTupleToDrawColour(tpl);
}
void setAtomNoteColour(RDKit::MolDrawOptions &self, nb::tuple tpl) {
  self.atomNoteColour = pyTupleToDrawColour(tpl);
}
nb::tuple getAtomNoteColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.atomNoteColour);
}
void setBondNoteColour(RDKit::MolDrawOptions &self, nb::tuple tpl) {
  self.bondNoteColour = pyTupleToDrawColour(tpl);
}
nb::tuple getBondNoteColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.bondNoteColour);
}
nb::tuple getVariableAttachmentColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.variableAttachmentColour);
}
void setVariableAttachmentColour(RDKit::MolDrawOptions &self, nb::tuple tpl) {
  self.variableAttachmentColour = pyTupleToDrawColour(tpl);
}
void useDefaultAtomPalette(RDKit::MolDrawOptions &self) {
  assignDefaultPalette(self.atomColourPalette);
}
void useBWAtomPalette(RDKit::MolDrawOptions &self) {
  assignBWPalette(self.atomColourPalette);
}
void useAvalonAtomPalette(RDKit::MolDrawOptions &self) {
  assignAvalonPalette(self.atomColourPalette);
}
void useCDKAtomPalette(RDKit::MolDrawOptions &self) {
  assignCDKPalette(self.atomColourPalette);
}
void updateAtomPalette(RDKit::MolDrawOptions &self, nb::object cmap) {
  pyDictToColourMap(cmap, self.atomColourPalette);
}
void setAtomPalette(RDKit::MolDrawOptions &self, nb::object cmap) {
  self.atomColourPalette.clear();
  updateAtomPalette(self, cmap);
}
nb::dict getAtomPalette(const RDKit::MolDrawOptions &self) {
  nb::dict res;
  for (const auto &pair : self.atomColourPalette) {
    res[nb::cast(pair.first)] = colourToPyTuple(pair.second);
  }
  return res;
}

void setMonochromeMode_helper1(RDKit::MolDrawOptions &options, nb::tuple fg,
                               nb::tuple bg) {
  auto fgc = pyTupleToDrawColour(fg);
  auto bgc = pyTupleToDrawColour(bg);
  RDKit::setMonochromeMode(options, fgc, bgc);
}

void setMonochromeMode_helper2(RDKit::MolDraw2D &d2d, nb::tuple fg,
                               nb::tuple bg) {
  auto fgc = pyTupleToDrawColour(fg);
  auto bgc = pyTupleToDrawColour(bg);
  RDKit::setMonochromeMode(d2d, fgc, bgc);
}

void contourAndDrawGaussiansHelper(
    RDKit::MolDraw2D &drawer, nb::object pylocs, nb::object pyheights,
    nb::object pywidths, unsigned int nContours, nb::object pylevels,
    const MolDraw2DUtils::ContourParams &params, nb::object mol) {
  std::unique_ptr<std::vector<RDGeom::Point2D>> locs =
      pythonObjectToVect<RDGeom::Point2D>(pylocs);
  if (!locs) {
    throw nb::value_error("locs argument must be non-empty");
  }
  std::unique_ptr<std::vector<double>> heights =
      pythonObjectToVect<double>(pyheights);
  if (!heights) {
    throw nb::value_error("heights argument must be non-empty");
  }
  std::unique_ptr<std::vector<double>> widths =
      pythonObjectToVect<double>(pywidths);
  if (!widths) {
    throw nb::value_error("widths argument must be non-empty");
  }
  std::unique_ptr<std::vector<double>> levels;
  if (!pylevels.is_none()) {
    levels = pythonObjectToVect<double>(pylevels);
  } else {
    levels = std::unique_ptr<std::vector<double>>(new std::vector<double>);
  }
  ROMol *mol_p = nullptr;
  if (!mol.is_none()) {
    mol_p = nb::cast<ROMol *>(mol);
  }
  MolDraw2DUtils::contourAndDrawGaussians(drawer, *locs, *heights, *widths,
                                          nContours, *levels, params, mol_p);
}

void contourAndDrawGridHelper(
    RDKit::MolDraw2D &drawer,
    nb::ndarray<nb::numpy, double, nb::ndim<2>, nb::c_contig> data,
    nb::object pyxcoords, nb::object pyycoords, unsigned int nContours,
    nb::object pylevels, const MolDraw2DUtils::ContourParams &params,
    nb::object mol) {
  std::unique_ptr<std::vector<double>> xcoords =
      pythonObjectToVect<double>(pyxcoords);
  if (!xcoords) {
    throw nb::value_error("xcoords argument must be non-empty");
  }
  std::unique_ptr<std::vector<double>> ycoords =
      pythonObjectToVect<double>(pyycoords);
  if (!ycoords) {
    throw nb::value_error("ycoords argument must be non-empty");
  }
  std::unique_ptr<std::vector<double>> levels;
  if (!pylevels.is_none()) {
    levels = pythonObjectToVect<double>(pylevels);
  } else {
    levels = std::unique_ptr<std::vector<double>>(new std::vector<double>);
  }

  if (data.shape(0) != xcoords->size()) {
    throw nb::value_error(
        "data array and xcoords sizes do not match.\n"
        "Did you forget to call np.transpose() on the array?");
  }

  if (data.shape(1) != ycoords->size()) {
    throw nb::value_error("data array and ycoords sizes do not match");
  }

  ROMol *mol_p = nullptr;
  if (!mol.is_none()) {
    mol_p = nb::cast<RDKit::ROMol *>(mol);
  }
  MolDraw2DUtils::contourAndDrawGrid(drawer, data.data(), *xcoords, *ycoords,
                                     nContours, *levels, params, mol_p);
}

void setColoursHelper(RDKit::MolDraw2DUtils::ContourParams &params,
                      nb::object pycolors) {
  std::vector<RDKit::DrawColour> cs;
  for (size_t i = 0; i < nb::len(pycolors); ++i) {
    cs.push_back(pyTupleToDrawColour(nb::cast<nb::tuple>(pycolors[i])));
  }
  params.colourMap = cs;
}

nb::tuple getColoursHelper(
    const RDKit::MolDraw2DUtils::ContourParams &params) {
  nb::list res;
  for (const auto &clr : params.colourMap) {
    res.append(colourToPyTuple(clr));
  }
  return nb::tuple(res);
}

void setContourColour(RDKit::MolDraw2DUtils::ContourParams &params,
                      nb::tuple tpl) {
  params.contourColour = pyTupleToDrawColour(tpl);
}

nb::tuple getContourColour(
    const RDKit::MolDraw2DUtils::ContourParams &params) {
  return colourToPyTuple(params.contourColour);
}

void drawPolygonHelper(RDKit::MolDraw2D &self, nb::object py_cds,
                       bool rawCoords) {
  std::unique_ptr<std::vector<RDGeom::Point2D>> cds =
      pythonObjectToVect<RDGeom::Point2D>(py_cds);
  if (!cds) {
    throw nb::value_error("cds argument must be non-empty");
  }

  self.drawPolygon(*cds, rawCoords);
}

void drawAttachmentLineHelper(RDKit::MolDraw2D &self, const Point2D &cds1,
                              const Point2D &cds2, nb::tuple pycol,
                              double len, unsigned int nSegments,
                              bool rawCoords) {
  auto col = pyTupleToDrawColour(pycol);
  self.drawAttachmentLine(cds1, cds2, col, len, nSegments, rawCoords);
}

void drawWavyLineHelper(RDKit::MolDraw2D &self, const Point2D &cds1,
                        const Point2D &cds2, nb::tuple pycol1,
                        nb::tuple pycol2, unsigned int nSegments,
                        double vertOffset, bool rawCoords) {
  auto col1 = pyTupleToDrawColour(pycol1);
  auto col2 = pyTupleToDrawColour(pycol2);
  self.drawWavyLine(cds1, cds2, col1, col2, nSegments, vertOffset, rawCoords);
}

void drawArrowHelper(RDKit::MolDraw2D &self, const Point2D &cds1,
                     const Point2D &cds2, bool asPolygon, double frac,
                     double angle, nb::object pycol, bool rawCoords) {
  DrawColour col{0.0, 0.0, 0.0};
  if (!pycol.is_none()) {
    col = pyTupleToDrawColour(nb::cast<nb::tuple>(pycol));
  }
  self.drawArrow(cds1, cds2, asPolygon, frac, angle, col, rawCoords);
}

void setDrawOptions(RDKit::MolDraw2D &self, const MolDrawOptions &opts) {
  self.drawOptions() = opts;
}

void setDrawerColour(RDKit::MolDraw2D &self, nb::tuple tpl) {
  self.setColour(pyTupleToDrawColour(tpl));
}

void updateMolDrawOptionsHelper(RDKit::MolDrawOptions &obj, std::string json) {
  MolDraw2DUtils::updateMolDrawOptionsFromJSON(obj, json);
}

void updateDrawerParamsHelper(RDKit::MolDraw2D &obj, std::string json) {
  MolDraw2DUtils::updateDrawerParamsFromJSON(obj, json);
}

std::string molToSVG(const ROMol &mol, unsigned int width, unsigned int height,
                     nb::object pyHighlightAtoms, bool kekulize,
                     unsigned int lineWidthMult, bool includeAtomCircles,
                     int confId) {
  RDUNUSED_PARAM(kekulize);
  std::unique_ptr<std::vector<int>> highlightAtoms =
      pythonObjectToVect(pyHighlightAtoms, static_cast<int>(mol.getNumAtoms()));
  std::stringstream outs;
  MolDraw2DSVG drawer(width, height, outs);
  drawer.setLineWidth(drawer.lineWidth() * lineWidthMult);
  drawer.drawOptions().circleAtoms = includeAtomCircles;
  drawer.drawOptions().prepareMolsBeforeDrawing = false;
  drawer.drawMolecule(mol, highlightAtoms.get(), nullptr, nullptr, confId);
  drawer.finishDrawing();
  return outs.str();
}

std::string molToACS1996SVG(const ROMol &mol, std::string legend,
                            nb::object highlight_atoms,
                            nb::object highlight_bonds,
                            nb::object highlight_atom_map,
                            nb::object highlight_bond_map,
                            nb::object highlight_atom_radii, int confId) {
  std::stringstream outs;
  MolDraw2DSVG drawer(-1, -1, outs);
  drawMoleculeACS1996Helper(drawer, mol, legend, highlight_atoms,
                            highlight_bonds, highlight_atom_map,
                            highlight_bond_map, highlight_atom_radii, confId);
  drawer.finishDrawing();
  return outs.str();
}

void drawStringHelper(MolDraw2D &self, std::string text, const Point2D &loc,
                      int align, bool rawCoords) {
  MolDraw2D_detail::TextAlignType talign =
      MolDraw2D_detail::TextAlignType::MIDDLE;
  switch (align) {
    case 0:
      talign = MolDraw2D_detail::TextAlignType::MIDDLE;
      break;
    case 1:
      talign = MolDraw2D_detail::TextAlignType::START;
      break;
    case 2:
      talign = MolDraw2D_detail::TextAlignType::END;
      break;
    default:
      throw nb::value_error("align must be 0, 1, or 2");
  }
  self.drawString(text, loc, talign, rawCoords);
}

void setScaleHelper(MolDraw2D &self, int width, int height, const Point2D &minv,
                    const Point2D &maxv, nb::object mol) {
  ROMol *mol_p = nullptr;
  if (!mol.is_none()) {
    mol_p = nb::cast<RDKit::ROMol *>(mol);
  }
  self.setScale(width, height, minv, maxv, mol_p);
}

}  // namespace

NB_MODULE(rdMolDraw2D, m) {
  m.doc() = "Module containing a C++ implementation of 2D molecule drawing";

  rdkit_import_array();

  nb::enum_<RDKit::MultiColourHighlightStyle>(m, "MultiColourHighlightStyle")
      .value("CircleAndLine", RDKit::MultiColourHighlightStyle::CIRCLEANDLINE)
      .value("Lasso", RDKit::MultiColourHighlightStyle::LASSO)
      .export_values();

  nb::enum_<RDKit::MolDrawOptions::LegendPosition>(m, "LegendPosition")
      .value("Bottom", RDKit::MolDrawOptions::LegendPosition::Bottom)
      .value("Top", RDKit::MolDrawOptions::LegendPosition::Top)
      .value("Left", RDKit::MolDrawOptions::LegendPosition::Left)
      .value("Right", RDKit::MolDrawOptions::LegendPosition::Right)
      .export_values();

  {
    auto drawElementEnum =
        nb::enum_<RDKit::DrawElement::_enumerated>(m, "DrawElement", nb::is_arithmetic());
    for (const auto *key : RDKit::DrawElement::_names()) {
      drawElementEnum.value(key, RDKit::DrawElement::_from_string(key));
    }
    drawElementEnum.export_values();
  }

  nb::class_<IntStringMap>(m, "IntStringMap")
      .def("__getitem__",
           [](const IntStringMap &self, int key) {
             auto it = self.dp_map->find(key);
             if (it == self.dp_map->end()) throw nb::key_error("key not found");
             return it->second;
           })
      .def("__setitem__",
           [](IntStringMap &self, int key, const std::string &val) {
             (*self.dp_map)[key] = val;
           })
      .def("__delitem__",
           [](IntStringMap &self, int key) {
             if (!self.dp_map->erase(key)) throw nb::key_error("key not found");
           })
      .def("__len__", [](const IntStringMap &self) { return self.dp_map->size(); })
      .def("__contains__",
           [](const IntStringMap &self, int key) {
             return self.dp_map->count(key) > 0;
           })
      .def("__repr__",
           [](const IntStringMap &self) {
             std::string s = "{";
             bool first = true;
             for (const auto &p : *self.dp_map) {
               if (!first) s += ", ";
               s += std::to_string(p.first) + ": '" + p.second + "'";
               first = false;
             }
             return s + "}";
           });

  nb::class_<RDKit::MolDrawOptions>(m, "MolDrawOptions", "Drawing options")
      .def(nb::init<>())
      .def_rw("dummiesAreAttachments",
               &RDKit::MolDrawOptions::dummiesAreAttachments)
      .def_rw("circleAtoms", &RDKit::MolDrawOptions::circleAtoms)
      .def_rw("splitBonds", &RDKit::MolDrawOptions::splitBonds)
      .def_prop_rw(
          "backgroundColour", &getBgColour, &setBgColour,
          "the background colour as an (R,G,B,A) tuple, values should be between 0 and 1")
      .def_prop_rw(
          "queryColour", &getQyColour, &setQyColour,
          "the query colour as an (R,G,B,A) tuple, values should be between 0 and 1")
      .def_prop_rw(
          "highlightColour", &getHighlightColour, &setHighlightColour,
          "the highlight colour as an (R,G,B,A) tuple, values should be between 0 and 1")
      .def_prop_rw(
          "symbolColour", &getSymbolColour, &setSymbolColour,
          "the symbol colour as an (R,G,B,A) tuple, values should be between 0 and 1")
      .def_prop_rw(
          "annotationColour", &getAnnotationColour, &setAnnotationColour,
          "the annotation colour as an (R,G,B,A) tuple, values should be between 0 and 1")
      .def_prop_rw(
          "atomNoteColour", &getAtomNoteColour, &setAtomNoteColour,
          "the atom note colour as an (R,G,B,A) tuple, values should be between 0 and 1")
      .def_prop_rw(
          "bondNoteColour", &getBondNoteColour, &setBondNoteColour,
          "the bond note colour as an (R,G,B,A) tuple, values should be between 0 and 1")
      .def_prop_rw(
          "legendColour", &getLegendColour, &setLegendColour,
          "the legend colour as an (R,G,B,A) tuple, values should be between 0 and 1")
      .def_prop_rw(
          "variableAttachmentColour", &getVariableAttachmentColour,
          &setVariableAttachmentColour,
          "the variable attachment colour as an (R,G,B,A) tuple, values should be between 0 and 1")
      .def("getBackgroundColour", &getBgColour, "method returning the background colour")
      .def("getQueryColour", &getQyColour, "method returning the query colour")
      .def("getHighlightColour", &getHighlightColour, "method returning the highlight colour")
      .def("setBackgroundColour", &setBgColour, "tpl"_a,
           "method for setting the background colour")
      .def("setQueryColour", &setQyColour, "tpl"_a,
           "method for setting the query colour")
      .def("setHighlightColour", &setHighlightColour, "tpl"_a,
           "method for setting the highlight colour")
      .def("getSymbolColour", &getSymbolColour, "method returning the symbol colour")
      .def("setSymbolColour", &setSymbolColour, "tpl"_a,
           "method for setting the symbol colour")
      .def("getAnnotationColour", &getAnnotationColour, "method returning the annotation colour")
      .def("setAnnotationColour", &setAnnotationColour, "tpl"_a,
           "method for setting the annotation colour")
      .def("setAtomNoteColour", &setAtomNoteColour, "tpl"_a,
           "method for setting the atom note colour")
      .def("getAtomNoteColour", &getAtomNoteColour, "method returning the atom note colour")
      .def("setBondNoteColour", &setBondNoteColour, "tpl"_a,
           "method for setting the bond note colour")
      .def("getBondNoteColour", &getBondNoteColour, "method returning the bond note colour")
      .def("getLegendColour", &getLegendColour, "method returning the legend colour")
      .def("setLegendColour", &setLegendColour, "tpl"_a,
           "method for setting the legend colour")
      .def("useDefaultAtomPalette", &useDefaultAtomPalette, "use the default colour palette for atoms and bonds")
      .def("useBWAtomPalette", &useBWAtomPalette, "use a black and white palette for atoms and bonds")
      .def("useAvalonAtomPalette", &useAvalonAtomPalette, "use the Avalon renderer palette for atoms and bonds")
      .def("useCDKAtomPalette", &useCDKAtomPalette, "use the CDK palette for atoms and bonds")
      .def("updateAtomPalette", &updateAtomPalette, "cmap"_a,
           "updates the palette for atoms and bonds from a dictionary mapping "
           "ints to 3-tuples")
      .def("setAtomPalette", &setAtomPalette, "cmap"_a,
           "sets the palette for atoms and bonds from a dictionary mapping ints "
           "to 3-tuples")
      .def("getAtomPalette", &getAtomPalette, "returns the current atom palette as a dictionary mapping ints "
           "to 4-tuples")
      .def_prop_rw(
          "atomLabels",
          [](RDKit::MolDrawOptions &self) {
            return IntStringMap{&self.atomLabels};
          },
          [](RDKit::MolDrawOptions &self, nb::dict d) {
            self.atomLabels.clear();
            for (auto item : d) {
              self.atomLabels[nb::cast<int>(item.first)] =
                  nb::cast<std::string>(item.second);
            }
          },
          nb::rv_policy::reference_internal,
          "maps indices to atom labels")
      .def_rw("atomLabelDeuteriumTritium",
               &RDKit::MolDrawOptions::atomLabelDeuteriumTritium,
               "labels deuterium as D and tritium as T")
      .def_rw("continuousHighlight",
               &RDKit::MolDrawOptions::continuousHighlight)
      .def_rw("fillHighlights", &RDKit::MolDrawOptions::fillHighlights)
      .def_rw("highlightRadius", &RDKit::MolDrawOptions::highlightRadius,
               "Default radius for highlight circles.")
      .def_rw("flagCloseContactsDist",
               &RDKit::MolDrawOptions::flagCloseContactsDist)
      .def_rw("atomRegions", &RDKit::MolDrawOptions::atomRegions,
               "regions to outline")
      .def_rw("includeAtomTags", &RDKit::MolDrawOptions::includeAtomTags,
               "include atom tags in output")
      .def_rw("clearBackground", &RDKit::MolDrawOptions::clearBackground,
               "clear the background before drawing a molecule")
      .def_rw("legendFontSize", &RDKit::MolDrawOptions::legendFontSize,
               "font size in pixels of the legend (if drawn)")
      .def_rw("legendFraction", &RDKit::MolDrawOptions::legendFraction,
               "fraction of the draw panel to be used for the legend if present")
      .def_rw("legendPosition", &RDKit::MolDrawOptions::legendPosition,
               R"DOC(legend position enum. Default=Bottom.
Values: LegendPosition.Bottom, LegendPosition.Top, LegendPosition.Left, LegendPosition.Right.)DOC")
      .def_rw("legendVerticalText", &RDKit::MolDrawOptions::legendVerticalText,
               "when legend is Left or Right, draw text vertically (one char per line)")
      .def_rw("maxFontSize", &RDKit::MolDrawOptions::maxFontSize,
               "maximum font size in pixels. default=40, -1 means no maximum.")
      .def_rw("minFontSize", &RDKit::MolDrawOptions::minFontSize,
               "minimum font size in pixels. default=6, -1 means no minimum.")
      .def_rw("fixedFontSize", &RDKit::MolDrawOptions::fixedFontSize,
               R"DOC(font size in pixels. default=-1 means not fixed.  If set,
always used irrespective of scale, minFontSize and maxFontSize.)DOC")
      .def_rw("baseFontSize", &RDKit::MolDrawOptions::baseFontSize,
               "relative size of font.  Defaults to 0.6.  -1 means use default.")
      .def_rw("annotationFontScale", &RDKit::MolDrawOptions::annotationFontScale,
               R"DOC(Scale of font for atom and bond annotation relative to atom
label font.  Default=0.75.)DOC")
      .def_rw("fontFile", &RDKit::MolDrawOptions::fontFile,
               R"DOC(Font file for use with FreeType text drawer.  Can also be
BuiltinTelexRegular (the default) or BuiltinRobotoRegular.)DOC")
      .def_rw("multipleBondOffset", &RDKit::MolDrawOptions::multipleBondOffset,
               R"DOC(offset for the extra lines in a multiple bond as a fraction of mean bond length)DOC")
      .def_rw("padding", &RDKit::MolDrawOptions::padding,
               "Fraction of empty space to leave around molecule.  Default=0.05.")
      .def_rw("reagentPadding", &RDKit::MolDrawOptions::componentPadding,
               R"DOC(Fraction of empty space to leave around each component
of a reaction drawing.  Default=0.0.)DOC")
      .def_rw("bondLineWidth", &RDKit::MolDrawOptions::bondLineWidth,
               "if positive, this overrides the default line width for bonds")
      .def_rw("scaleBondWidth", &RDKit::MolDrawOptions::scaleBondWidth,
               "Scales the width of drawn bonds using image scaling.")
      .def_rw("scaleHighlightBondWidth",
               &RDKit::MolDrawOptions::scaleHighlightBondWidth,
               "Scales the width of drawn highlighted bonds using image scaling.")
      .def_rw("highlightBondWidthMultiplier",
               &RDKit::MolDrawOptions::highlightBondWidthMultiplier,
               "What to multiply default bond width by for highlighting bonds. Default-8.")
      .def_rw("prepareMolsBeforeDrawing",
               &RDKit::MolDrawOptions::prepareMolsBeforeDrawing,
               "call prepareMolForDrawing() on each molecule passed to DrawMolecules()")
      .def_rw("fixedScale", &RDKit::MolDrawOptions::fixedScale,
               R"DOC(If > 0.0, fixes scale to that fraction of width of
draw window unless that would make it too big.  Default -1.0 means adjust scale to fit.)DOC")
      .def_rw("fixedBondLength", &RDKit::MolDrawOptions::fixedBondLength,
               R"DOC(If > 0.0, fixes bond length to this number of pixels
unless that would make it too big.  Default -1.0 means
no fix.  If both set, fixedScale takes precedence.)DOC")
      .def_rw("rotate", &RDKit::MolDrawOptions::rotate,
               "Rotates molecule about centre by this number of degrees,")
      .def_rw("addStereoAnnotation",
               &RDKit::MolDrawOptions::addStereoAnnotation,
               "adds R/S and E/Z to drawings. Default False.")
      .def_rw("showAllCIPCodes", &RDKit::MolDrawOptions::showAllCIPCodes,
               "show all defined CIP codes (no hiding!). Default False.")
      .def_rw("addAtomIndices", &RDKit::MolDrawOptions::addAtomIndices,
               "adds atom indices to drawings. Default False.")
      .def_rw("addBondIndices", &RDKit::MolDrawOptions::addBondIndices,
               "adds bond indices to drawings. Default False.")
      .def_rw("isotopeLabels", &RDKit::MolDrawOptions::isotopeLabels,
               "adds isotope labels on non-dummy atoms. Default True.")
      .def_rw("dummyIsotopeLabels",
               &RDKit::MolDrawOptions::dummyIsotopeLabels,
               "adds isotope labels on dummy atoms. Default True.")
      .def_rw("atomHighlightsAreCircles",
               &RDKit::MolDrawOptions::atomHighlightsAreCircles,
               R"DOC(forces atom highlights always to be circles.
Default (false) is to put ellipses round longer labels.)DOC")
      .def_rw("multiColourHighlightStyle",
               &RDKit::MolDrawOptions::multiColourHighlightStyle,
               R"DOC(Either 'CircleAndLine' or 'Lasso', to control style of
multi-coloured highlighting in DrawMoleculeWithHighlights.
Default is CircleAndLine.)DOC")
      .def_rw("centreMoleculesBeforeDrawing",
               &RDKit::MolDrawOptions::centreMoleculesBeforeDrawing,
               "Moves the centre of the drawn molecule to (0,0). Default False.")
      .def_rw("additionalAtomLabelPadding",
               &RDKit::MolDrawOptions::additionalAtomLabelPadding,
               R"DOC(additional padding to leave around atom labels.
Expressed as a fraction of the font size.)DOC")
      .def_rw("noAtomLabels", &RDKit::MolDrawOptions::noAtomLabels,
               "disables inclusion of atom labels in the rendering")
      .def_rw("explicitMethyl", &RDKit::MolDrawOptions::explicitMethyl,
               "Draw terminal methyls explictly.  Default is false.")
      .def_rw("includeMetadata", &RDKit::MolDrawOptions::includeMetadata,
               R"DOC(When possible, include metadata about molecules and reactions to
allow them to be reconstructed. Default is true.)DOC")
      .def_rw("includeRadicals", &RDKit::MolDrawOptions::includeRadicals,
               R"DOC(include radicals in the drawing (it can be useful to turn this off
for reactions and queries). Default is true.)DOC")
      .def_rw("comicMode", &RDKit::MolDrawOptions::comicMode,
               R"DOC(simulate hand-drawn lines for bonds. When combined with
a font like Comic-Sans or Comic-Neue, this gives
xkcd-like drawings. Default is false.)DOC")
      .def_rw("variableBondWidthMultiplier",
               &RDKit::MolDrawOptions::variableBondWidthMultiplier,
               "what to multiply standard bond width by for variable attachment points.")
      .def_rw("variableAtomRadius",
               &RDKit::MolDrawOptions::variableAtomRadius,
               "radius value to use for atoms involved in variable attachment points.")
      .def_rw("includeChiralFlagLabel",
               &RDKit::MolDrawOptions::includeChiralFlagLabel,
               R"DOC(add a molecule annotation with "ABS" if the chiral
flag is set. Default is false.)DOC")
      .def_rw("simplifiedStereoGroupLabel",
               &RDKit::MolDrawOptions::simplifiedStereoGroupLabel,
               R"DOC(if all specified stereocenters are in a single
StereoGroup, show a molecule-level annotation instead of
the individual labels. Default is false.)DOC")
      .def_rw("unspecifiedStereoIsUnknown",
               &RDKit::MolDrawOptions::unspecifiedStereoIsUnknown,
               R"DOC(if true, double bonds with unspecified stereo are drawn
crossed, potential stereocenters with unspecified stereo
are drawn with a wavy bond. Default is false.)DOC")
      .def_rw("singleColourWedgeBonds",
               &RDKit::MolDrawOptions::singleColourWedgeBonds,
               R"DOC(if true wedged and dashed bonds are drawn using symbolColour
rather than inheriting their colour from the atoms.
Default is false.)DOC")
      .def_rw("singleColourBonds",
               &RDKit::MolDrawOptions::singleColourBonds,
               "if true all bonds are drawn using symbolColour "
               "rather than inheriting their colour from the atoms. "
               "Default is false.")
      .def_rw("useMolBlockWedging",
               &RDKit::MolDrawOptions::useMolBlockWedging,
               R"DOC(If the molecule came from a MolBlock, prefer the wedging
information that provides.  If false, use RDKit rules.
Default false)DOC")
      .def_rw("scalingFactor", &RDKit::MolDrawOptions::scalingFactor,
               R"DOC(scaling factor for pixels->angstrom when auto scaling
being used.  Default is 20.)DOC")
      .def_rw("drawMolsSameScale",
               &RDKit::MolDrawOptions::drawMolsSameScale,
               R"DOC(when drawing multiple molecules with DrawMolecules,
forces them to use the same scale.  Default is true.)DOC")
      .def_rw("useComplexQueryAtomSymbols",
               &RDKit::MolDrawOptions::useComplexQueryAtomSymbols,
               R"DOC(replace any atom, any hetero, any halo queries
with complex query symbols A, Q, X, M, optionally followed
by H if hydrogen is included (except for AH, which stays *).
Default is true)DOC")
      .def_rw("bracketsAroundAtomLists",
               &RDKit::MolDrawOptions::bracketsAroundAtomLists,
               R"DOC(Whether to put brackets round atom lists in query atoms.
Default is true.)DOC")
      .def_rw("standardColoursForHighlightedAtoms",
               &RDKit::MolDrawOptions::standardColoursForHighlightedAtoms,
               R"DOC(If true, highlighted hetero atoms are drawn in standard colours
rather than black.  Default=False)DOC")
      .def_rw("drawingExtentsInclude",
               &RDKit::MolDrawOptions::drawingExtentsInclude,
               R"DOC(Drawing extents are computed taking into account only selected
DrawElement items.  Default=DrawElement.ALL)DOC")
      .def("getVariableAttachmentColour", &getVariableAttachmentColour,
           "method for getting the colour of variable attachment points")
      .def("setVariableAttachmentColour", &setVariableAttachmentColour,
           "tpl"_a,
           "method for setting the colour of variable attachment points")
      .def("__setattr__", &safeSetattr);

  nb::class_<RDKit::MolDraw2D>(m, "MolDraw2D", "Drawer abstract base class")
      .def("SetFontSize", &RDKit::MolDraw2D::setFontSize, "new_size"_a,
           "change the default font size. The units are, roughly, pixels.")
      .def("FontSize", &RDKit::MolDraw2D::fontSize, "get the default font size. The units are, roughly, pixels.")
      .def("DrawMolecule", &drawMoleculeHelper1,
           "mol"_a,
           "highlightAtoms"_a = nb::none(),
           "highlightAtomColors"_a = nb::none(),
           "highlightAtomRadii"_a = nb::none(),
           "confId"_a = -1, "legend"_a = std::string(""),
           "renders a molecule\n")
      .def("DrawMolecule", &drawMoleculeHelper2,
           "mol"_a,
           "highlightAtoms"_a,
           "highlightBonds"_a,
           "highlightAtomColors"_a = nb::none(),
           "highlightBondColors"_a = nb::none(),
           "highlightAtomRadii"_a = nb::none(),
           "confId"_a = -1, "legend"_a = std::string(""),
           "renders a molecule\n")
      .def("GetMolSize", &getMolSizeHelper,
           "mol"_a,
           "highlightAtoms"_a = nb::none(),
           "highlightBonds"_a = nb::none(),
           "highlightAtomColors"_a = nb::none(),
           "highlightBondColors"_a = nb::none(),
           "highlightAtomRadii"_a = nb::none(),
           "confId"_a = -1, "legend"_a = std::string(""),
           "returns the width and height required to draw a molecule at the current size")
      .def("DrawMoleculeWithHighlights", &drawMoleculeWithHighlightsHelper,
           "mol"_a, "legend"_a,
           "highlight_atom_map"_a,
           "highlight_bond_map"_a, "highlight_radii"_a,
           "highlight_linewidth_multipliers"_a,
           "confId"_a = -1,
           "renders a molecule with multiple highlight colours\n")
      .def("DrawMolecules", &drawMoleculesHelper2,
           "mols"_a,
           "highlightAtoms"_a = nb::none(),
           "highlightBonds"_a = nb::none(),
           "highlightAtomColors"_a = nb::none(),
           "highlightBondColors"_a = nb::none(),
           "highlightAtomRadii"_a = nb::none(),
           "confIds"_a = nb::none(),
           "legends"_a = nb::none(),
           "renders multiple molecules\n")
      .def("DrawReaction", &drawReactionHelper,
           "rxn"_a,
           "highlightByReactant"_a = false,
           "highlightColorsReactants"_a = nb::none(),
           "confIds"_a = nb::none(),
           "renders a reaction\n")
      .def("Width", &RDKit::MolDraw2D::width, "get the width of the drawing canvas")
      .def("Height", &RDKit::MolDraw2D::height, "get the height of the drawing canvas")
      .def("SetOffset", &RDKit::MolDraw2D::setOffset, "x"_a, "y"_a,
           "set the offset (in drawing coordinates) for the drawing")
      .def("Offset", &RDKit::MolDraw2D::offset, "returns the offset (in drawing coordinates) for the drawing")
      .def("SetScale", &setScaleHelper,
           "width"_a, "height"_a,
           "minv"_a, "maxv"_a,
           "mol"_a = nb::none(),
           "uses the values provided to set the drawing scaling")
      .def("FlexiMode", &RDKit::MolDraw2D::flexiMode, "returns whether or not FlexiMode is being used")
      .def("SetFlexiMode", &RDKit::MolDraw2D::setFlexiMode, "mode"_a,
           R"DOC(when FlexiMode is set, molecules will always been drawn with the default values for bond length, font size, etc.)DOC")
      .def("SetLineWidth", &RDKit::MolDraw2D::setLineWidth, "width"_a,
           "set the line width being used")
      .def("SetColour", &setDrawerColour, "tpl"_a,
           "set the color being used fr drawing and filling")
      .def("LineWidth", &RDKit::MolDraw2D::lineWidth, "returns the line width being used")
      .def("SetFillPolys", &RDKit::MolDraw2D::setFillPolys, "val"_a,
           "sets whether or not polygons are filled")
      .def("FillPolys", &RDKit::MolDraw2D::fillPolys, "returns whether or not polygons are being filled")
      .def("DrawLine",
           (void(RDKit::MolDraw2D::*)(const Point2D &, const Point2D &, bool)) &
               RDKit::MolDraw2D::drawLine,
           "cds1"_a, "cds2"_a, "rawCoords"_a = false,
           "draws a line with the current drawing style. The coordinates "
           "are in the molecule frame unless rawCoords is true, "
           "in which case the coordinates are in pixels.")
      .def("DrawArrow", &drawArrowHelper,
           "cds1"_a, "cds2"_a,
           "asPolygon"_a = false, "frac"_a = 0.05,
           "angle"_a = M_PI / 6,
           "color"_a = nb::none(),
           "rawCoords"_a = false,
           R"DOC(draws an arrow with the current drawing style. The coordinates
are in the molecule frame unless rawCoords is true,
in which case the coordinates are in pixels.
If asPolygon is true the head of the
arrow will be drawn as a triangle, otherwise two lines are used.
The fraction of the arrow length to use for the head is given by
frac. The angle of the arrowhead
(the angle between the main line and each arrowhead line) is given by angle.
The color is a tuple of 3 floats (0-1) in red, green, blue (RGB) order.)DOC")
      .def("DrawTriangle", &RDKit::MolDraw2D::drawTriangle,
           "cds1"_a, "cds2"_a, "cds3"_a, "rawCoords"_a = false,
           "draws a triangle with the current drawing style. The coordinates "
           "are in the molecule frame unless rawCoords is true, "
           "in which case the coordinates are in pixels.")
      .def("DrawPolygon", &drawPolygonHelper,
           "cds"_a, "rawCoords"_a = false,
           "draws a polygon with the current drawing style. The coordinates "
           "are in the molecule frame unless rawCoords is true, "
           "in which case the coordinates are in pixels.")
      .def("DrawEllipse", &RDKit::MolDraw2D::drawEllipse,
           "cds1"_a, "cds2"_a, "rawCoords"_a = false,
           "draws a triangle with the current drawing style in the rectangle "
           "defined by the two points. The coordinates "
           "are in the molecule frame unless rawCoords is true, "
           "in which case the coordinates are in pixels.")
      .def("DrawRect", &RDKit::MolDraw2D::drawRect,
           "cds1"_a, "cds2"_a, "rawCoords"_a = false,
           "draws a rectangle with the current drawing style in the rectangle "
           "defined by the two points. The coordinates "
           "are in the molecule frame unless rawCoords is true, "
           "in which case the coordinates are in pixels.")
      .def("DrawArc",
           (void(RDKit::MolDraw2D::*)(const Point2D &, double, double, double,
                                      bool)) &
               RDKit::MolDraw2D::drawArc,
           "center"_a, "radius"_a, "angle1"_a, "angle2"_a,
           "rawCoords"_a = false,
           R"DOC(draws an arc with the current drawing style. The coordinates
are in the molecule frame unless rawCoords is true,
in which case the coordinates are in pixels.
The angles are in degrees; angle2 should be > angle1.)DOC")
      .def("DrawAttachmentLine", &drawAttachmentLineHelper,
           "cds1"_a, "cds2"_a,
           "color"_a, "len"_a = 1.0,
           "nSegments"_a = 16, "rawCoords"_a = false,
           R"DOC(draw a line indicating the presence of an attachment point
(normally a squiggle line perpendicular to a bond).
The coordinates
are in the molecule frame unless rawCoords is true,
in which case the coordinates are in pixels.)DOC")
      .def("DrawWavyLine", &drawWavyLineHelper,
           "cds1"_a, "cds2"_a,
           "color1"_a, "color2"_a,
           "nSegments"_a = 16, "vertOffset"_a = 0.05,
           "rawCoords"_a = false,
           R"DOC(draw a line indicating the presence of an attachment point
(normally a squiggle line perpendicular to a bond).
The coordinates
are in the molecule frame unless rawCoords is true,
in which case the coordinates are in pixels.)DOC")
      .def("DrawString",
           (void(RDKit::MolDraw2D::*)(const std::string &,
                                      const RDGeom::Point2D &, bool)) &
               RDKit::MolDraw2D::drawString,
           "string"_a, "pos"_a, "rawCoords"_a = false,
           "add text to the canvas. The coordinates "
           "are in the molecule frame unless rawCoords is true, "
           "in which case the coordinates are in pixels.")
      .def("DrawString", &drawStringHelper,
           "string"_a, "pos"_a,
           "align"_a, "rawCoords"_a = false,
           R"DOC(add aligned text to the canvas. The align argument can be 0
(=MIDDLE), 1 (=START), or 2 (=END).
The coordinates
are in the molecule frame unless rawCoords is true,
in which case the coordinates are in pixels.)DOC")
      .def("GetDrawCoords",
           (RDGeom::Point2D(RDKit::MolDraw2D::*)(const RDGeom::Point2D &)
                const) &
               RDKit::MolDraw2D::getDrawCoords,
           "point"_a,
           "get the coordinates in drawing space for a particular point in "
           "molecule space")
      .def("GetDrawCoords",
           (RDGeom::Point2D(RDKit::MolDraw2D::*)(int) const) &
               RDKit::MolDraw2D::getDrawCoords,
           "atomIndex"_a,
           "get the coordinates in drawing space for a particular atom")
      .def("ClearDrawing", &RDKit::MolDraw2D::clearDrawing, "clears the drawing by filling it with the background color")
      .def("drawOptions",
           (RDKit::MolDrawOptions & (RDKit::MolDraw2D::*)()) &
               RDKit::MolDraw2D::drawOptions,
           nb::rv_policy::reference_internal,
           "Returns a modifiable version of the current drawing options")
      .def("SetDrawOptions", &setDrawOptions, "opts"_a,
           "Copies the drawing options passed in over our drawing options");

  nb::class_<RDKit::MolDraw2DSVG, RDKit::MolDraw2D>(m, "MolDraw2DSVG",
                                                     "SVG molecule drawer")
      .def(nb::init<int, int, int, int, bool>(),
           "width"_a, "height"_a,
           "panelWidth"_a = -1, "panelHeight"_a = -1,
           "noFreetype"_a = false)
      .def("FinishDrawing", &RDKit::MolDraw2DSVG::finishDrawing, "add the last bits of SVG to finish the drawing")
      .def("AddMoleculeMetadata",
           (void(RDKit::MolDraw2DSVG::*)(const RDKit::ROMol &, int) const) &
               RDKit::MolDraw2DSVG::addMoleculeMetadata,
           "mol"_a, "confId"_a = -1,
           "add RDKit-specific information to the bottom of the drawing")
      .def("TagAtoms", &tagAtomHelper,
           "mol"_a,
           "radius"_a = 0.2,
           "events"_a = nb::none(),
           "allow atom selection in the SVG")
      .def("GetDrawingText", &RDKit::MolDraw2DSVG::getDrawingText, "return the SVG");

#ifdef RDK_BUILD_CAIRO_SUPPORT
  nb::class_<RDKit::MolDraw2DCairo, RDKit::MolDraw2D>(m, "MolDraw2DCairo",
                                                       "Cairo molecule drawer")
      .def(nb::init<int, int, int, int, bool>(),
           "width"_a, "height"_a,
           "panelWidth"_a = -1, "panelHeight"_a = -1,
           "noFreetype"_a = false)
      .def("FinishDrawing", &RDKit::MolDraw2DCairo::finishDrawing, "add the last bits to finish the drawing")
      .def("GetDrawingText", &getCairoDrawingText, "return the PNG data as a string")
      .def("WriteDrawingText", &RDKit::MolDraw2DCairo::writeDrawingText,
           "fName"_a,
           "write the PNG data to the named file");
#endif

  m.def("PrepareMolForDrawing", &prepMolForDrawing,
        "mol"_a.none(), "kekulize"_a = true,
        "addChiralHs"_a = true, "wedgeBonds"_a = true,
        "forceCoords"_a = false, "wavyBonds"_a = false,
        R"DOC(Does some cleanup operations on the molecule to prepare it to draw nicely.
The operations include: kekulization, addition of chiral Hs (so that we can draw
wedges to them), wedging of bonds at chiral centers, and generation of a 2D
conformation if the molecule does not already have a conformation

Returns a modified copy of the molecule.)DOC",
        nb::rv_policy::take_ownership);

  m.def("PrepareAndDrawMolecule", &prepareAndDrawMoleculeHelper,
        "drawer"_a, "mol"_a, "legend"_a = "",
        "highlightAtoms"_a = nb::none(),
        "highlightBonds"_a = nb::none(),
        "highlightAtomColors"_a = nb::none(),
        "highlightBondColors"_a = nb::none(),
        "highlightAtomRadii"_a = nb::none(),
        "confId"_a = -1, "kekulize"_a = true,
        "Preps a molecule for drawing and actually draws it\n");

  m.def("DrawMoleculeACS1996", &drawMoleculeACS1996Helper,
        "drawer"_a, "mol"_a, "legend"_a = "",
        "highlightAtoms"_a = nb::none(),
        "highlightBonds"_a = nb::none(),
        "highlightAtomColors"_a = nb::none(),
        "highlightBondColors"_a = nb::none(),
        "highlightAtomRadii"_a = nb::none(),
        "confId"_a = -1,
        "Draws molecule in ACS 1996 mode.");

  nb::class_<RDKit::MolDraw2DUtils::ContourParams>(m, "ContourParams",
                                                    "Parameters for drawing contours")
      .def(nb::init<>())
      .def_rw("setScale",
               &RDKit::MolDraw2DUtils::ContourParams::setScale,
               "set the scale of the drawing object (useful if you draw "
               "the grid/contours first)")
      .def_rw("dashNegative",
               &RDKit::MolDraw2DUtils::ContourParams::dashNegative,
               "use a dashed line for negative contours")
      .def_rw("fillGrid",
               &RDKit::MolDraw2DUtils::ContourParams::fillGrid,
               "colors the grid in addition to drawing contours")
      .def_rw("gridResolution",
               &RDKit::MolDraw2DUtils::ContourParams::gridResolution,
               "set the resolution of the grid")
      .def_rw("contourWidth",
               &RDKit::MolDraw2DUtils::ContourParams::contourWidth,
               "line width of the contours")
      .def_rw("extraGridPadding",
               &RDKit::MolDraw2DUtils::ContourParams::extraGridPadding,
               "extra space (in molecule coords) around the grid")
      .def_rw("drawAsLines", &RDKit::MolDraw2DUtils::ContourParams::drawAsLines,
               "draw the contours as continuous lines isntead of line segments")
      .def_rw("coordScaleForQuantization",
               &RDKit::MolDraw2DUtils::ContourParams::coordScaleForQuantization,
               "scaling factor used to convert coordinates to ints when forming the continuous lines")
      .def_rw("isovalScaleForQuantization",
               &RDKit::MolDraw2DUtils::ContourParams::isovalScaleForQuantization,
               "scaling factor used to convert isovalues to ints when forming the continuous lines")
      .def_rw("useFillThreshold",
               &RDKit::MolDraw2DUtils::ContourParams::useFillThreshold,
               "use a magnitude threshold to determine if a grid point is filled")
      .def_rw("fillThreshold",
               &RDKit::MolDraw2DUtils::ContourParams::fillThreshold,
               "magnitude threshold to determine if a grid point is filled")
      .def_rw("fillThresholdIsFraction",
               &RDKit::MolDraw2DUtils::ContourParams::fillThresholdIsFraction,
               "if true, fillThreshold is a fraction of the range of the data")
      .def_prop_rw("colourMap", &getColoursHelper, &setColoursHelper,
                   "the color map to use when filling the grid")
      .def_prop_rw("contourColour", &getContourColour, &setContourColour,
                   "the color to use for drawing the contours")
      .def("setContourColour", &setContourColour, "colour"_a)
      .def("setColourMap", &setColoursHelper, "colours"_a)
      .def("__setattr__", &safeSetattr);

  m.def("ContourAndDrawGaussians", &contourAndDrawGaussiansHelper,
        "drawer"_a, "locs"_a, "heights"_a,
        "widths"_a, "nContours"_a = 10,
        "levels"_a = nb::none(),
        "params"_a = RDKit::MolDraw2DUtils::ContourParams(),
        "mol"_a = nb::none(),
        R"DOC(Generates and draws contours for a set of gaussians

- drawer: the MolDraw2D object to use
- locs: locations of the gaussians
- heights: the heights (or weights) of the gaussians
- widths: the standard deviations of the gaussians
- nContours: the number of contours to draw
- levels: the contours to use
- ps: additional parameters controlling the contouring.
- mol: molecule used to help set scale.

The values are calculated on a grid with spacing params.gridResolution.
If params.setScale  is set, the grid size will be calculated based on the
locations of the gaussians and params.extraGridPadding. Otherwise the current
size of the viewport will be used.

If the levels argument is empty, the contour levels will be determined
automatically from the max and min values on the grid and levels will
be updated to include the contour levels.

If params.fillGrid is set, the data on the grid will also be drawn using
the color scheme in params.colourMap

If mol is not 0, uses the molecule to help set the scale, assuming that
it will be drawn over the plot, so needs to fit on it.)DOC");

  m.def("ContourAndDrawGrid", &contourAndDrawGridHelper,
        "drawer"_a, "data"_a, "xcoords"_a,
        "ycoords"_a, "nContours"_a = 10,
        "levels"_a = nb::none(),
        "params"_a = RDKit::MolDraw2DUtils::ContourParams(),
        "mol"_a = nb::none(),
        R"DOC(Generates and draws contours for data on a grid

- drawer: the MolDraw2D object to use
- data: numpy array with the data to be contoured
- xcoords: the x coordinates of the grid
- ycoords: the y coordinates of the grid
- nContours: the number of contours to draw
- levels: the contours to use
- ps: additional parameters controlling the contouring
- mol: molecule used to help set scale.

The values are calculated on a grid with spacing params.gridResolution.
If params.setScale  is set, the grid size will be calculated based on the
locations of the gaussians and params.extraGridPadding. Otherwise the current
size of the viewport will be used.

If the levels argument is empty, the contour levels will be determined
automatically from the max and min values on the grid and levels will
be updated to include the contour levels.

If params.fillGrid is set, the data on the grid will also be drawn using
the color scheme in params.colourMap

If mol is not 0, uses the molecule to help set the scale, assuming that
it will be drawn over the plot, so needs to fit on it.)DOC");

  m.def("UpdateMolDrawOptionsFromJSON", &updateMolDrawOptionsHelper,
        "opts"_a, "json"_a);
  m.def("UpdateDrawerParamsFromJSON", &updateDrawerParamsHelper,
        "drawer"_a, "json"_a);

  m.def("MolToSVG", &molToSVG,
        "mol"_a, "width"_a = 300,
        "height"_a = 300,
        "highlightAtoms"_a = nb::none(),
        "kekulize"_a = true, "lineWidthMult"_a = 1,
        "includeAtomCircles"_a = true, "confId"_a = -1,
        "Returns svg for a molecule");

  m.def("MolToACS1996SVG", &molToACS1996SVG,
        "mol"_a, "legend"_a = "",
        "highlightAtoms"_a = nb::none(),
        "highlightBonds"_a = nb::none(),
        "highlightAtomColors"_a = nb::none(),
        "highlightBondColors"_a = nb::none(),
        "highlightAtomRadii"_a = nb::none(),
        "confId"_a = -1,
        "Returns ACS 1996 mode svg for a molecule");

  m.def("SetACS1996Mode", &RDKit::MolDraw2DUtils::setACS1996Options,
        "drawOptions"_a, "meanBondLength"_a,
        R"DOC(Set the draw options to produce something as close as possible to
the ACS 1996 guidelines as described at
https://en.wikipedia.org/wiki/Wikipedia:Manual_of_Style/Chemistry/Structure_drawing

- MolDrawOptions opt - the options what will be changed
- float meanBondLength - mean bond length of the molecule

Works best if the MolDraw2D object is created with width and height -1 (a
flexiCanvas).
The mean bond length may be calculated with MeanBondLength.
It is used to calculate the offset for the lines in multiple bonds.

Options changed are:
  bondLineWidth = 0.6
  scaleBondWidth = false
  scalingFactor = 14.4 / meanBondLen
  multipleBondOffset = 0.18
  highlightBondWidthMultiplier = 32
  setMonochromeMode - black and white
  fixedFontSize = 10
  additionalAtomLabelPadding = 0.066
  fontFile - if it isn't set already, then if RDBASE is set and the file
             exists, uses $RDBASE/Data/Fonts/FreeSans.ttf.  Otherwise uses
             BuiltinRobotoRegular.)DOC");

  m.def("MeanBondLength", &RDKit::MolDraw2DUtils::meanBondLength,
        "mol"_a, "confId"_a = -1,
        "Calculate the mean bond length for the molecule.");

  m.def("SetDarkMode",
        (void (*)(RDKit::MolDrawOptions &)) & RDKit::setDarkMode,
        "d2d"_a, "set dark mode for a MolDrawOptions object");
  m.def("SetDarkMode",
        (void (*)(RDKit::MolDraw2D &)) & RDKit::setDarkMode,
        "d2d"_a, "set dark mode for a MolDraw2D object");
  m.def("SetMonochromeMode", &setMonochromeMode_helper1,
        "options"_a, "fgColour"_a, "bgColour"_a,
        "set monochrome mode for a MolDrawOptions object");
  m.def("SetMonochromeMode", &setMonochromeMode_helper2,
        "drawer"_a, "fgColour"_a, "bgColour"_a,
        "set monochrome mode for a MolDraw2D object");
}
