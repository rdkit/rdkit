//
//  Copyright (C) 2015-2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdmoldraw2d_array_API
#include <RDBoost/python.h>
#include <RDBoost/boost_numpy.h>
#include <numpy/arrayobject.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>

#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <Geometry/point.h>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#ifdef RDK_BUILD_CAIRO_SUPPORT
#include <cairo.h>
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#endif

namespace python = boost::python;

namespace RDKit {
namespace {
void tagAtomHelper(MolDraw2DSVG &self, const ROMol &mol, double radius,
                   python::object pyo) {
  std::map<std::string, std::string> events;
  if (pyo) {
    python::dict tDict = python::extract<python::dict>(pyo);
    python::list keys = tDict.keys();
    python::list vals = tDict.values();
    for (unsigned int i = 0;
         i < python::extract<unsigned int>(keys.attr("__len__")()); ++i) {
      events[python::extract<std::string>(keys[i])] =
          python::extract<std::string>(vals[i]);
    }
  }
  self.tagAtoms(mol, radius, events);
}
void pyDictToColourMap(python::object pyo, ColourPalette &res) {
  python::dict tDict = python::extract<python::dict>(pyo);
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(tDict.keys().attr("__len__")()); ++i) {
    python::tuple tpl = python::extract<python::tuple>(tDict.values()[i]);
    float r = python::extract<float>(tpl[0]);
    float g = python::extract<float>(tpl[1]);
    float b = python::extract<float>(tpl[2]);
    DrawColour clr(r, g, b);
    res[python::extract<int>(tDict.keys()[i])] = clr;
  }
}
ColourPalette *pyDictToColourMap(python::object pyo) {
  ColourPalette *res = nullptr;
  if (pyo) {
    res = new ColourPalette;
    pyDictToColourMap(pyo, *res);
  }
  return res;
}
void pyDictToDoubleMap(python::object pyo, std::map<int, double> &res) {
  python::dict tDict = python::extract<python::dict>(pyo);
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(tDict.keys().attr("__len__")()); ++i) {
    double r = python::extract<double>(tDict.values()[i]);
    res[python::extract<int>(tDict.keys()[i])] = r;
  }
}
std::map<int, double> *pyDictToDoubleMap(python::object pyo) {
  std::map<int, double> *res = nullptr;
  if (pyo) {
    res = new std::map<int, double>;
    pyDictToDoubleMap(pyo, *res);
  }
  return res;
}
void pyDictToIntMap(python::object pyo, std::map<int, int> &res) {
  python::dict tDict = python::extract<python::dict>(pyo);
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(tDict.keys().attr("__len__")()); ++i) {
    int r = python::extract<int>(tDict.values()[i]);
    res[python::extract<int>(tDict.keys()[i])] = r;
  }
}
std::map<int, int> *pyDictToIntMap(python::object pyo) {
  std::map<int, int> *res = nullptr;
  if (pyo) {
    res = new std::map<int, int>;
    pyDictToIntMap(pyo, *res);
  }
  return res;
}

DrawColour pyTupleToDrawColour(const python::tuple tpl) {
  float r = python::extract<float>(tpl[0]);
  if (r > 1 || r < 0) {
    throw ValueErrorException("RGBA color value needs to be between 0 and 1.");
  }
  float g = python::extract<float>(tpl[1]);
  if (g > 1 || g < 0) {
    throw ValueErrorException("RGBA color value needs to be between 0 and 1.");
  }
  float b = python::extract<float>(tpl[2]);
  if (b > 1 || b < 0) {
    throw ValueErrorException("RGBA color value needs to be between 0 and 1.");
  }
  float a = 1;
  if (python::extract<unsigned int>(tpl.attr("__len__")()) > 3) {
    a = python::extract<float>(tpl[3]);
    if (a > 1 || a < 0) {
      throw ValueErrorException(
          "RGBA color value needs to be between 0 and 1.");
    }
  }
  DrawColour clr(r, g, b, a);
  return clr;
}
void pyListToColourVec(python::object pyo, std::vector<DrawColour> &res) {
  res.clear();
  python::list tList = python::extract<python::list>(pyo);
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(tList.attr("__len__")()); ++i) {
    python::tuple tpl = python::extract<python::tuple>(tList[i]);
    res.push_back(pyTupleToDrawColour(tpl));
  }
}
}  // namespace

void pyDictToMapColourVec(python::object pyo,
                          std::map<int, std::vector<DrawColour>> &res) {
  python::dict tDict = python::extract<python::dict>(pyo);
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(tDict.keys().attr("__len__")()); ++i) {
    python::list pl = python::extract<python::list>(tDict.values()[i]);
    std::vector<DrawColour> v;
    pyListToColourVec(pl, v);
    res[python::extract<int>(tDict.keys()[i])] = v;
  }
}
std::map<int, std::vector<DrawColour>> *pyDictToMapColourVec(
    python::object pyo) {
  std::map<int, std::vector<DrawColour>> *res = nullptr;
  if (pyo) {
    res = new std::map<int, std::vector<DrawColour>>;
    pyDictToMapColourVec(pyo, *res);
  }
  return res;
}

void drawMoleculeHelper1(MolDraw2D &self, const ROMol &mol,
                         python::object highlight_atoms,
                         python::object highlight_atom_map,
                         python::object highlight_atom_radii, int confId,
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
                         python::object highlight_atoms,
                         python::object highlight_bonds,
                         python::object highlight_atom_map,
                         python::object highlight_bond_map,
                         python::object highlight_atom_radii, int confId,
                         std::string legend) {
  std::unique_ptr<std::vector<int>> highlightAtoms =
      pythonObjectToVect(highlight_atoms, static_cast<int>(mol.getNumAtoms()));
  std::unique_ptr<std::vector<int>> highlightBonds =
      pythonObjectToVect(highlight_bonds, static_cast<int>(mol.getNumBonds()));
  // FIX: support these
  ColourPalette *ham = pyDictToColourMap(highlight_atom_map);
  ColourPalette *hbm = pyDictToColourMap(highlight_bond_map);
  std::map<int, double> *har = pyDictToDoubleMap(highlight_atom_radii);

  self.drawMolecule(mol, legend, highlightAtoms.get(), highlightBonds.get(),
                    ham, hbm, har, confId);

  delete ham;
  delete hbm;
  delete har;
}

void drawMoleculeWithHighlightsHelper(
    MolDraw2D &self, const ROMol &mol, std::string legend,
    python::object highlight_atom_map, python::object highlight_bond_map,
    python::object highlight_atom_radii,
    python::object highlight_linewidth_multipliers, int confId) {
  // highlight_atom_map and highlight_bond_map come in as a dict of
  // lists of tuples of floats (the R, G, B values for the colours),
  // and need to be changed to a map of vectors of DrawColour.
  // All of the dict to map converters return nullptr if the dict
  // was empty.  We need real objects in all cases for
  // drawMoleculeWithHighlights.
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
    python::object highlight_atoms, python::object highlight_bonds,
    python::object highlight_atom_map, python::object highlight_bond_map,
    python::object highlight_atom_radii, int confId) {
  std::unique_ptr<std::vector<int>> highlightAtoms =
      pythonObjectToVect(highlight_atoms, static_cast<int>(mol.getNumAtoms()));
  std::unique_ptr<std::vector<int>> highlightBonds =
      pythonObjectToVect(highlight_bonds, static_cast<int>(mol.getNumBonds()));
  // FIX: support these
  ColourPalette *ham = pyDictToColourMap(highlight_atom_map);
  ColourPalette *hbm = pyDictToColourMap(highlight_bond_map);
  std::map<int, double> *har = pyDictToDoubleMap(highlight_atom_radii);
  MolDraw2DUtils::prepareAndDrawMolecule(
      drawer, mol, legend, highlightAtoms.get(), highlightBonds.get(), ham, hbm,
      har, confId);

  delete ham;
  delete hbm;
  delete har;
}

void drawMoleculesHelper2(MolDraw2D &self, python::object pmols,
                          python::object highlight_atoms,
                          python::object highlight_bonds,
                          python::object highlight_atom_map,
                          python::object highlight_bond_map,
                          python::object highlight_atom_radii,
                          python::object pconfIds, python::object plegends) {
  std::unique_ptr<std::vector<ROMol *>> mols =
      pythonObjectToVect<ROMol *>(pmols);
  if (mols == nullptr || !mols->size()) {
    return;
  }
  unsigned int nThere = mols->size();
  std::unique_ptr<std::vector<std::vector<int>>> highlightAtoms;
  if (highlight_atoms) {
    if (python::extract<unsigned int>(highlight_atoms.attr("__len__")()) !=
        nThere) {
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
  if (highlight_bonds) {
    if (python::extract<unsigned int>(highlight_bonds.attr("__len__")()) !=
        nThere) {
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
  if (highlight_atom_map) {
    if (python::extract<unsigned int>(highlight_atom_map.attr("__len__")()) !=
        nThere) {
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
  if (highlight_bond_map) {
    if (python::extract<unsigned int>(highlight_bond_map.attr("__len__")()) !=
        nThere) {
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
  if (highlight_atom_radii) {
    if (python::extract<unsigned int>(highlight_atom_radii.attr("__len__")()) !=
        nThere) {
      throw ValueErrorException(
          "If highlightAtomRadii is provided it must be the same length as the "
          "molecule list.");
    }
    highlightRadii.reset(new std::vector<std::map<int, double>>(nThere));
    for (unsigned int i = 0; i < nThere; ++i) {
      pyDictToDoubleMap(highlight_atom_radii[i], (*highlightRadii)[i]);
    }
  }
  // std::unique_ptr<std::vector<int> > highlightAtoms =
  //     pythonObjectToVect(highlight_atoms,
  //     static_cast<int>(mol.getNumAtoms()));
  // std::unique_ptr<std::vector<int> > highlightBonds =
  //     pythonObjectToVect(highlight_bonds,
  //     static_cast<int>(mol.getNumBonds()));
  // FIX: support these
  // std::map<int, DrawColour> *ham = pyDictToColourMap(highlight_atom_map);
  // std::map<int, DrawColour> *hbm = pyDictToColourMap(highlight_bond_map);
  // std::map<int, double> *har = pyDictToDoubleMap(highlight_atom_radii);
  //
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
                        python::object phighlightColorsReactants,
                        python::object pconfIds) {
  std::unique_ptr<std::vector<DrawColour>> highlightColorsReactants;
  if (phighlightColorsReactants) {
    highlightColorsReactants.reset(new std::vector<DrawColour>);
    pyListToColourVec(phighlightColorsReactants, *highlightColorsReactants);
  }

  std::unique_ptr<std::vector<int>> confIds = pythonObjectToVect<int>(pconfIds);

  self.drawReaction(rxn, highlightByReactant, highlightColorsReactants.get(),
                    confIds.get());
}

#ifdef RDK_BUILD_CAIRO_SUPPORT
python::object getCairoDrawingText(const RDKit::MolDraw2DCairo &self) {
  std::string res = self.getDrawingText();
  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}
#endif
ROMol *prepMolForDrawing(const ROMol *m, bool kekulize = true,
                         bool addChiralHs = true, bool wedgeBonds = true,
                         bool forceCoords = false) {
  auto *res = new RWMol(*m);
  MolDraw2DUtils::prepareMolForDrawing(*res, kekulize, addChiralHs, wedgeBonds,
                                       forceCoords);
  return static_cast<ROMol *>(res);
}

python::tuple colourToPyTuple(const DrawColour &clr) {
  python::list res;
  res.append(clr.r);
  res.append(clr.g);
  res.append(clr.b);
  res.append(clr.a);

  return python::tuple(res);
}
python::object getBgColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.backgroundColour);
}
python::object getHighlightColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.highlightColour);
}
void setBgColour(RDKit::MolDrawOptions &self, python::tuple tpl) {
  self.backgroundColour = pyTupleToDrawColour(tpl);
}
void setHighlightColour(RDKit::MolDrawOptions &self, python::tuple tpl) {
  self.highlightColour = pyTupleToDrawColour(tpl);
}
python::object getSymbolColour(const RDKit::MolDrawOptions &self) {
  return colourToPyTuple(self.symbolColour);
}
void setSymbolColour(RDKit::MolDrawOptions &self, python::tuple tpl) {
  self.symbolColour = pyTupleToDrawColour(tpl);
}
void useDefaultAtomPalette(RDKit::MolDrawOptions &self) {
  assignDefaultPalette(self.atomColourPalette);
}
void useBWAtomPalette(RDKit::MolDrawOptions &self) {
  assignBWPalette(self.atomColourPalette);
}
void updateAtomPalette(RDKit::MolDrawOptions &self, python::object cmap) {
  pyDictToColourMap(cmap, self.atomColourPalette);
}
void setAtomPalette(RDKit::MolDrawOptions &self, python::object cmap) {
  self.atomColourPalette.clear();
  updateAtomPalette(self, cmap);
}

void contourAndDrawGaussiansHelper(
    RDKit::MolDraw2D &drawer, python::object pylocs, python::object pyheights,
    python::object pywidths, unsigned int nContours, python::object pylevels,
    const MolDraw2DUtils::ContourParams &params, python::object mol) {
  std::unique_ptr<std::vector<RDGeom::Point2D>> locs =
      pythonObjectToVect<RDGeom::Point2D>(pylocs);
  std::unique_ptr<std::vector<double>> heights =
      pythonObjectToVect<double>(pyheights);
  std::unique_ptr<std::vector<double>> widths =
      pythonObjectToVect<double>(pywidths);
  std::unique_ptr<std::vector<double>> levels;
  if (pylevels) {
    levels = pythonObjectToVect<double>(pylevels);
  } else {
    levels = std::unique_ptr<std::vector<double>>(new std::vector<double>);
  }
  ROMol *mol_p = nullptr;
  if (mol) {
    mol_p = python::extract<ROMol *>(mol);
  }
  MolDraw2DUtils::contourAndDrawGaussians(drawer, *locs, *heights, *widths,
                                          nContours, *levels, params, mol_p);
}

void contourAndDrawGridHelper(RDKit::MolDraw2D &drawer, python::object &data,
                              python::object &pyxcoords,
                              python::object &pyycoords, unsigned int nContours,
                              python::object &pylevels,
                              const MolDraw2DUtils::ContourParams &params,
                              python::object mol) {
  if (!PyArray_Check(data.ptr())) {
    throw_value_error("data argument must be a numpy array");
  }
  auto *dataArr = reinterpret_cast<PyArrayObject *>(
      PyArray_ContiguousFromObject(data.ptr(), NPY_DOUBLE, 2, 2));
  if (!dataArr) {
    throw_value_error("could not convert data argument");
  }

  std::unique_ptr<std::vector<double>> xcoords =
      pythonObjectToVect<double>(pyxcoords);
  std::unique_ptr<std::vector<double>> ycoords =
      pythonObjectToVect<double>(pyycoords);
  std::unique_ptr<std::vector<double>> levels;
  if (pylevels) {
    levels = pythonObjectToVect<double>(pylevels);
  } else {
    levels = std::unique_ptr<std::vector<double>>(new std::vector<double>);
  }

  if (PyArray_DIM(dataArr, 0) != static_cast<int>(xcoords->size())) {
    throw_value_error(
        "data array and xcoords sizes do not match.\n"
        "Did you forget to call np.transpose() on the array?");
  }

  if (PyArray_DIM(dataArr, 1) != static_cast<int>(ycoords->size())) {
    throw_value_error("data array and ycoords sizes do not match");
  }

  ROMol *mol_p = nullptr;
  if (mol) {
    mol_p = python::extract<RDKit::ROMol *>(mol);
  }
  MolDraw2DUtils::contourAndDrawGrid(drawer, (double *)PyArray_DATA(dataArr),
                                     *xcoords, *ycoords, nContours, *levels,
                                     params, mol_p);

  Py_DECREF(dataArr);
}

void setColoursHelper(RDKit::MolDraw2DUtils::ContourParams &params,
                      python::object pycolors) {
  std::vector<RDKit::DrawColour> cs;
  for (size_t i = 0; i < python::extract<size_t>(pycolors.attr("__len__")());
       ++i) {
    cs.push_back(
        pyTupleToDrawColour(python::extract<python::tuple>(pycolors[i])));
  }
  params.colourMap = cs;
}

void setContourColour(RDKit::MolDraw2DUtils::ContourParams &params,
                      python::tuple tpl) {
  params.contourColour = pyTupleToDrawColour(tpl);
}
void drawPolygonHelper(RDKit::MolDraw2D &self, python::object py_cds) {
  std::unique_ptr<std::vector<RDGeom::Point2D>> cds =
      pythonObjectToVect<RDGeom::Point2D>(py_cds);

  self.drawPolygon(*cds);
}

void drawAttachmentLineHelper(RDKit::MolDraw2D &self, const Point2D &cds1,
                              const Point2D &cds2, python::tuple &pycol,
                              double len, unsigned int nSegments) {
  auto col = pyTupleToDrawColour(pycol);
  self.drawAttachmentLine(cds1, cds2, col, len, nSegments);
}

void drawWavyLineHelper(RDKit::MolDraw2D &self, const Point2D &cds1,
                        const Point2D &cds2, python::tuple &pycol1,
                        python::tuple &pycol2, unsigned int nSegments,
                        double vertOffset) {
  auto col1 = pyTupleToDrawColour(pycol1);
  auto col2 = pyTupleToDrawColour(pycol2);
  self.drawWavyLine(cds1, cds2, col1, col2, nSegments, vertOffset);
}

void setDrawOptions(RDKit::MolDraw2D &self, const MolDrawOptions &opts) {
  self.drawOptions() = opts;
}

}  // namespace RDKit

BOOST_PYTHON_MODULE(rdMolDraw2D) {
  python::scope().attr("__doc__") =
      "Module containing a C++ implementation of 2D molecule drawing";

  rdkit_import_array();

  python::class_<std::map<int, std::string>>("IntStringMap")
      .def(python::map_indexing_suite<std::map<int, std::string>, true>());

  std::string docString = "Drawing options";
  python::class_<RDKit::MolDrawOptions, boost::noncopyable>("MolDrawOptions",
                                                            docString.c_str())
      .def_readwrite("dummiesAreAttachments",
                     &RDKit::MolDrawOptions::dummiesAreAttachments)
      .def_readwrite("circleAtoms", &RDKit::MolDrawOptions::circleAtoms)
      .def("getBackgroundColour", &RDKit::getBgColour,
           "method returning the background colour")
      .def("getHighlightColour", &RDKit::getHighlightColour,
           "method returning the highlight colour")
      .def("setBackgroundColour", &RDKit::setBgColour,
           "method for setting the background colour")
      .def("setHighlightColour", &RDKit::setHighlightColour,
           "method for setting the highlight colour")
      .def("getSymbolColour", &RDKit::getSymbolColour,
           "method returning the symbol colour")
      .def("setSymbolColour", &RDKit::setSymbolColour,
           "method for setting the symbol colour")
      .def("useDefaultAtomPalette", &RDKit::useDefaultAtomPalette,
           "use the default colour palette for atoms and bonds")
      .def("useBWAtomPalette", &RDKit::useBWAtomPalette,
           "use the black & white palette for atoms and bonds")
      .def("updateAtomPalette", &RDKit::updateAtomPalette,
           "updates the palette for atoms and bonds from a dictionary mapping "
           "ints to 3-tuples")
      .def(
          "setAtomPalette", &RDKit::setAtomPalette,
          "sets the palette for atoms and bonds from a dictionary mapping ints "
          "to 3-tuples")
      .def_readwrite("atomLabels", &RDKit::MolDrawOptions::atomLabels,
                     "maps indices to atom labels")
      .def_readwrite("atomLabelDeuteriumTritium",
                     &RDKit::MolDrawOptions::atomLabelDeuteriumTritium,
                     "labels deuterium as D and tritium as T")
      .def_readwrite("continuousHighlight",
                     &RDKit::MolDrawOptions::continuousHighlight)
      .def_readwrite("fillHighlights", &RDKit::MolDrawOptions::fillHighlights)
      .def_readwrite("highlightRadius", &RDKit::MolDrawOptions::highlightRadius,
                     "Default radius for highlight circles.")
      .def_readwrite("flagCloseContactsDist",
                     &RDKit::MolDrawOptions::flagCloseContactsDist)
      .def_readwrite("atomRegions", &RDKit::MolDrawOptions::atomRegions,
                     "regions to outline")
      .def_readwrite("includeAtomTags", &RDKit::MolDrawOptions::includeAtomTags,
                     "include atom tags in output")
      .def_readwrite("clearBackground", &RDKit::MolDrawOptions::clearBackground,
                     "clear the background before drawing a molecule")
      .def_readwrite("legendFontSize", &RDKit::MolDrawOptions::legendFontSize,
                     "font size in pixels of the legend (if drawn)")
      .def_readwrite("maxFontSize", &RDKit::MolDrawOptions::maxFontSize,
                     "maximum font size in pixels. default=40, -1 means no"
                     " maximum.")
      .def_readwrite("minFontSize", &RDKit::MolDrawOptions::minFontSize,
                     "minimum font size in pixels. default=12, -1 means no"
                     " minimum.")
      .def_readwrite(
          "annotationFontScale", &RDKit::MolDrawOptions::annotationFontScale,
          "Scale of font for atom and bond annotation relative to atom"
          "label font.  Default=0.75.")
      .def_readwrite("fontFile", &RDKit::MolDrawOptions::fontFile,
                     "Font file for use with FreeType text drawer.")
      .def_readwrite(
          "multipleBondOffset", &RDKit::MolDrawOptions::multipleBondOffset,
          "offset (in Angstroms) for the extra lines in a multiple bond")
      .def_readwrite("padding", &RDKit::MolDrawOptions::padding,
                     "fraction of empty space to leave around molecule")
      .def_readwrite(
          "bondLineWidth", &RDKit::MolDrawOptions::bondLineWidth,
          "if positive, this overrides the default line width for bonds")
      .def_readwrite("scaleBondWidth", &RDKit::MolDrawOptions::scaleBondWidth,
                     "Scales the width of drawn bonds using image scaling.")
      .def_readwrite("scaleHighlightBondWidth",
                     &RDKit::MolDrawOptions::scaleHighlightBondWidth,
                     "Scales the width of drawn highlighted bonds using"
                     " image scaling.")
      .def_readwrite("highlightBondWidthMultiplier",
                     &RDKit::MolDrawOptions::highlightBondWidthMultiplier,
                     "What to multiply default bond width by for highlighting "
                     "bonds. Default-8.")
      .def_readwrite("prepareMolsBeforeDrawing",
                     &RDKit::MolDrawOptions::prepareMolsBeforeDrawing,
                     "call prepareMolForDrawing() on each molecule passed to "
                     "DrawMolecules()")
      .def_readwrite("fixedScale", &RDKit::MolDrawOptions::fixedScale,
                     "If > 0.0, fixes scale to that fraction of width of"
                     "draw window unless that would make it too big."
                     "  Default -1.0 means adjust scale to fit.")
      .def_readwrite("fixedBondLength", &RDKit::MolDrawOptions::fixedBondLength,
                     "If > 0.0, fixes bond length to this number of pixels"
                     "unless that would make it too big.  Default -1.0 means"
                     "no fix.  If both set, fixedScale takes precedence.")
      .def_readwrite("fixedScale", &RDKit::MolDrawOptions::fixedScale,
                     "If > 0.0, fixes scale to that fraction of width of"
                     "draw window.  Default -1.0 means adjust scale to fit.")
      .def_readwrite("rotate", &RDKit::MolDrawOptions::rotate,
                     "Rotates molecule about centre by this number of degrees,")
      .def_readwrite("addStereoAnnotation",
                     &RDKit::MolDrawOptions::addStereoAnnotation,
                     "adds R/S and E/Z to drawings. Default False.")
      .def_readwrite("addAtomIndices", &RDKit::MolDrawOptions::addAtomIndices,
                     "adds atom indices drawings. Default False.")
      .def_readwrite("addBondIndices", &RDKit::MolDrawOptions::addBondIndices,
                     "adds bond indices drawings. Default False.")
      .def_readwrite("atomHighlightsAreCircles",
                     &RDKit::MolDrawOptions::atomHighlightsAreCircles,
                     "forces atom highlights always to be circles."
                     "Default (false) is to put ellipses round"
                     "longer labels.")
      .def_readwrite("centreMoleculesBeforeDrawing",
                     &RDKit::MolDrawOptions::centreMoleculesBeforeDrawing,
                     "Moves the centre of the drawn molecule to (0,0)."
                     "Default False.")
      .def_readwrite("additionalAtomLabelPadding",
                     &RDKit::MolDrawOptions::additionalAtomLabelPadding,
                     "additional padding to leave around atom labels. "
                     "Expressed as a fraction of the font size.")
      .def_readwrite("explicitMethyl", &RDKit::MolDrawOptions::explicitMethyl,
                     "Draw terminal methyls explictly.  Default is false.")
      .def_readwrite(
          "includeMetadata", &RDKit::MolDrawOptions::includeMetadata,
          "When possible, include metadata about molecules and reactions to "
          "allow them to be reconstructed. Default is true.");
  docString = "Drawer abstract base class";
  python::class_<RDKit::MolDraw2D, boost::noncopyable>(
      "MolDraw2D", docString.c_str(), python::no_init)
      .def("SetFontSize", &RDKit::MolDraw2D::setFontSize,
           "change the default font size. The units are, roughly, pixels.")
      .def("FontSize", &RDKit::MolDraw2D::fontSize,
           "get the default font size. The units are, roughly, pixels.")
      .def(
          "DrawMolecule", RDKit::drawMoleculeHelper1,
          (python::arg("self"), python::arg("mol"),
           python::arg("highlightAtoms") = python::object(),
           python::arg("highlightAtomColors") = python::object(),
           python::arg("highlightAtomRadii") = python::object(),
           python::arg("confId") = -1, python::arg("legend") = std::string("")),
          "renders a molecule\n")
      .def(
          "DrawMolecule", RDKit::drawMoleculeHelper2,
          (python::arg("self"), python::arg("mol"),
           python::arg("highlightAtoms"), python::arg("highlightBonds"),
           python::arg("highlightAtomColors") = python::object(),
           python::arg("highlightBondColors") = python::object(),
           python::arg("highlightAtomRadii") = python::object(),
           python::arg("confId") = -1, python::arg("legend") = std::string("")),
          "renders a molecule\n")
      .def("DrawMoleculeWithHighlights",
           RDKit::drawMoleculeWithHighlightsHelper,
           (python::arg("self"), python::arg("mol"), python::arg("legend"),
            python::arg("highlight_atom_map"),
            python::arg("highlight_bond_map"), python::arg("highlight_radii"),
            python::arg("highlight_linewidth_multipliers"),
            python::arg("confId") = -1),
           "renders a molecule with multiple highlight colours\n")
      .def("DrawMolecules", RDKit::drawMoleculesHelper2,
           (python::arg("self"), python::arg("mols"),
            python::arg("highlightAtoms") = python::object(),
            python::arg("highlightBonds") = python::object(),
            python::arg("highlightAtomColors") = python::object(),
            python::arg("highlightBondColors") = python::object(),
            python::arg("highlightAtomRadii") = python::object(),
            python::arg("confIds") = python::object(),
            python::arg("legends") = python::object()),
           "renders multiple molecules\n")
      .def("DrawReaction", RDKit::drawReactionHelper,
           (python::arg("self"), python::arg("rxn"),
            python::arg("highlightByReactant") = false,
            python::arg("highlightColorsReactants") = python::object(),
            python::arg("confIds") = python::object()),
           "renders a reaction\n")
      .def("Width", &RDKit::MolDraw2D::width,
           "get the width of the drawing canvas")
      .def("Height", &RDKit::MolDraw2D::height,
           "get the height of the drawing canvas")
      .def("SetOffset", &RDKit::MolDraw2D::setOffset,
           "set the offset (in drawing coordinates) for the drawing")
      .def("Offset", &RDKit::MolDraw2D::offset,
           "returns the offset (in drawing coordinates) for the drawing")
      .def("SetScale", &RDKit::MolDraw2D::setScale,
           "uses the values provided to set the drawing scaling")
      .def("SetLineWidth", &RDKit::MolDraw2D::setLineWidth,
           "set the line width being used")
      .def("LineWidth", &RDKit::MolDraw2D::lineWidth,
           "returns the line width being used")
      .def("SetFillPolys", &RDKit::MolDraw2D::setFillPolys,
           "sets whether or not polygons are filled")
      .def("FillPolys", &RDKit::MolDraw2D::fillPolys,
           "returns whether or not polygons are being filled")
      .def("DrawLine",
           (void (RDKit::MolDraw2D::*)(const Point2D &, const Point2D &)) &
               RDKit::MolDraw2D::drawLine,
           (python::arg("self"), python::arg("cds1"), python::arg("cds2")),
           "draws a line with the current drawing style. The coordinates "
           "are in the molecule frame")
      .def("DrawArrow", &RDKit::MolDraw2D::drawArrow,
           (python::arg("self"), python::arg("cds1"), python::arg("cds2"),
            python::arg("asPolygon") = false, python::arg("frac") = 0.05,
            python::arg("angle") = M_PI / 6),
           "draws an arrow with the current drawing style. The coordinates "
           "are in the molecule frame. If asPolygon is true the head of the "
           "arrow will be drawn as a triangle, otherwise two lines are used.")
      .def("DrawTriangle", &RDKit::MolDraw2D::drawTriangle,
           (python::arg("self"), python::arg("cds1"), python::arg("cds2"),
            python::arg("cds3")),
           "draws a triangle with the current drawing style. The coordinates "
           "are in the molecule frame")
      .def("DrawPolygon", RDKit::drawPolygonHelper,
           (python::arg("self"), python::arg("cds")),
           "draws a polygon with the current drawing style. The coordinates "
           "are in the molecule frame")
      .def("DrawEllipse", &RDKit::MolDraw2D::drawEllipse,
           (python::arg("self"), python::arg("cds1"), python::arg("cds2")),
           "draws a triangle with the current drawing style in the rectangle "
           "defined by the two points. The coordinates "
           "are in the molecule frame")
      .def("DrawRect", &RDKit::MolDraw2D::drawRect,
           (python::arg("self"), python::arg("cds1"), python::arg("cds2")),
           "draws a rectangle with the current drawing style in the rectangle "
           "defined by the two points. The coordinates "
           "are in the molecule frame")
      .def("DrawAttachmentLine", &RDKit::drawAttachmentLineHelper,
           (python::arg("self"), python::arg("cds1"), python::arg("cds2"),
            python::arg("color"), python::arg("len") = 1.0,
            python::arg("nSegments") = 16),
           "draw a line indicating the presence of an attachment point "
           "(normally a squiggle line perpendicular to a bond)")
      .def("DrawWavyLine", &RDKit::drawWavyLineHelper,
           (python::arg("self"), python::arg("cds1"), python::arg("cds2"),
            python::arg("color1"), python::arg("color2"),
            python::arg("nSegments") = 16, python::arg("vertOffset") = 0.05),
           "draw a line indicating the presence of an attachment point "
           "(normally a squiggle line perpendicular to a bond)")
      .def("DrawString",
           (void (RDKit::MolDraw2D::*)(const std::string &,
                                       const RDGeom::Point2D &)) &
               RDKit::MolDraw2D::drawString,
           (python::arg("self"), python::arg("string"), python::arg("pos")),
           "add text to the canvas")
      .def("DrawString",
           (void (RDKit::MolDraw2D::*)(const std::string &,
                                       const RDGeom::Point2D &,
                                       RDKit::TextAlignType)) &
               RDKit::MolDraw2D::drawString,
           (python::arg("self"), python::arg("string"), python::arg("pos"),
            python::arg("align")),
           "add aligned text to the canvas")
      .def("GetDrawCoords",
           (RDGeom::Point2D(RDKit::MolDraw2D::*)(const RDGeom::Point2D &)
                const) &
               RDKit::MolDraw2D::getDrawCoords,
           (python::arg("self"), python::arg("point")),
           "get the coordinates in drawing space for a particular point in "
           "molecule space")
      .def("GetDrawCoords",
           (RDGeom::Point2D(RDKit::MolDraw2D::*)(int) const) &
               RDKit::MolDraw2D::getDrawCoords,
           (python::arg("self"), python::arg("atomIndex")),
           "get the coordinates in drawing space for a particular atom")
      .def("ClearDrawing", &RDKit::MolDraw2D::clearDrawing,
           (python::arg("self")),
           "clears the drawing by filling it with the background color")
      .def("drawOptions",
           (RDKit::MolDrawOptions & (RDKit::MolDraw2D::*)()) &
               RDKit::MolDraw2D::drawOptions,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>(),
           "Returns a modifiable version of the current drawing options")
      .def("SetDrawOptions", RDKit::setDrawOptions,
           "Copies the drawing options passed in over our drawing options");

  docString = "SVG molecule drawer";
  python::class_<RDKit::MolDraw2DSVG, python::bases<RDKit::MolDraw2D>,
                 boost::noncopyable>("MolDraw2DSVG", docString.c_str(),
                                     python::init<int, int>())
      .def(python::init<int, int, int, int>())
      .def(python::init<int, int, int, int, bool>())
      .def("FinishDrawing", &RDKit::MolDraw2DSVG::finishDrawing,
           "add the last bits of SVG to finish the drawing")
      .def("AddMoleculeMetadata",
           (void (RDKit::MolDraw2DSVG::*)(const RDKit::ROMol &, int) const) &
               RDKit::MolDraw2DSVG::addMoleculeMetadata,
           (python::arg("mol"), python::arg("confId") = -1),
           "add RDKit-specific information to the bottom of the drawing")
      .def("TagAtoms", RDKit::tagAtomHelper,
           (python::arg("mol"), python::arg("radius") = 0.2,
            python::arg("events") = python::object()),
           "allow atom selection in the SVG")
      .def("GetDrawingText", &RDKit::MolDraw2DSVG::getDrawingText,
           "return the SVG");

#ifdef RDK_BUILD_CAIRO_SUPPORT
  docString = "Cairo molecule drawer";
  python::class_<RDKit::MolDraw2DCairo, python::bases<RDKit::MolDraw2D>,
                 boost::noncopyable>("MolDraw2DCairo", docString.c_str(),
                                     python::init<int, int>())
      .def(python::init<int, int, int, int>())
      .def("FinishDrawing", &RDKit::MolDraw2DCairo::finishDrawing,
           "add the last bits to finish the drawing")
      .def("GetDrawingText", &RDKit::getCairoDrawingText,
           "return the PNG data as a string")
      .def("WriteDrawingText", &RDKit::MolDraw2DCairo::writeDrawingText,
           "write the PNG data to the named file");
#endif
  docString =
      "Does some cleanup operations on the molecule to prepare it to draw "
      "nicely.\n"
      "The operations include: kekulization, addition of chiral Hs (so that we "
      "can draw\n"
      "wedges to them), wedging of bonds at chiral centers, and generation of "
      "a 2D\n"
      "conformation if the molecule does not already have a conformation\n"
      "\nReturns a modified copy of the molecule.\n";
  python::def(
      "PrepareMolForDrawing", &RDKit::prepMolForDrawing,
      (python::arg("mol"), python::arg("kekulize") = true,
       python::arg("addChiralHs") = true, python::arg("wedgeBonds") = true,
       python::arg("forceCoords") = false),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());
  python::def(
      "PrepareAndDrawMolecule", &RDKit::prepareAndDrawMoleculeHelper,
      (python::arg("drawer"), python::arg("mol"), python::arg("legend") = "",
       python::arg("highlightAtoms") = python::object(),
       python::arg("highlightBonds") = python::object(),
       python::arg("highlightAtomColors") = python::object(),
       python::arg("highlightBondColors") = python::object(),
       python::arg("highlightAtomRadii") = python::object(),
       python::arg("confId") = -1),
      "Preps a molecule for drawing and actually draws it\n");
  docString = "Parameters for drawing contours";
  python::class_<RDKit::MolDraw2DUtils::ContourParams>(
      "ContourParams", docString.c_str(), python::init<>())
      .def_readwrite("setScale",
                     &RDKit::MolDraw2DUtils::ContourParams::setScale,
                     "set the scale of the drawing object (useful if you draw "
                     "the grid/contours first)")
      .def_readwrite("dashNegative",
                     &RDKit::MolDraw2DUtils::ContourParams::dashNegative,
                     "use a dashed line for negative contours")
      .def_readwrite("fillGrid",
                     &RDKit::MolDraw2DUtils::ContourParams::fillGrid,
                     "colors the grid in addition to drawing contours")
      .def_readwrite("gridResolution",
                     &RDKit::MolDraw2DUtils::ContourParams::gridResolution,
                     "set the resolution of the grid")
      .def_readwrite("contourWidth",
                     &RDKit::MolDraw2DUtils::ContourParams::contourWidth,
                     "line width of the contours")
      .def_readwrite("extraGridPadding",
                     &RDKit::MolDraw2DUtils::ContourParams::extraGridPadding,
                     "extra space (in molecule coords) around the grid")
      .def("setContourColour", &RDKit::setContourColour,
           (python::arg("self"), python::arg("colour")))
      .def("setColourMap", &RDKit::setColoursHelper,
           (python::arg("self"), python::arg("colours")));
  docString = R"DOC(Generates and draws contours for a set of gaussians

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
  it will be drawn over the plot, so needs to fit on it.
*/)DOC";
  python::def(
      "ContourAndDrawGaussians", &RDKit::contourAndDrawGaussiansHelper,
      (python::arg("drawer"), python::arg("locs"), python::arg("heights"),
       python::arg("widths"), python::arg("nContours") = 10,
       python::arg("levels") = python::object(),
       python::arg("params") = RDKit::MolDraw2DUtils::ContourParams(),
       python::arg("mol") = python::object()),
      docString.c_str());
  docString = R"DOC(Generates and draws contours for data on a grid

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
  it will be drawn over the plot, so needs to fit on it.
*/)DOC";
  python::def(
      "ContourAndDrawGrid", &RDKit::contourAndDrawGridHelper,
      (python::arg("drawer"), python::arg("data"), python::arg("xcoords"),
       python::arg("ycoords"), python::arg("nContours") = 10,
       python::arg("levels") = python::object(),
       python::arg("params") = RDKit::MolDraw2DUtils::ContourParams(),
       python::arg("mol") = python::object()),
      docString.c_str());
}
