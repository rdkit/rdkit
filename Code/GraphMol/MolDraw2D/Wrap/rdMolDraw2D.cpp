//
//  Copyright (C) 2015-2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdmoldraw2d_array_API
#include <RDBoost/python.h>
#include <GraphMol/ROMol.h>
#include <RDBoost/Wrap.h>

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

DrawColour pyTupleToDrawColour(const python::tuple tpl) {
  float r = python::extract<float>(tpl[0]);
  if (r > 1 || r < 0) {
    throw ValueErrorException("RGB color value needs to be between 0 and 1.");
  }
  float g = python::extract<float>(tpl[1]);
  if (g > 1 || g < 0) {
    throw ValueErrorException("RGB color value needs to be between 0 and 1.");
  }
  float b = python::extract<float>(tpl[2]);
  if (b > 1 || b < 0) {
    throw ValueErrorException("RGB color value needs to be between 0 and 1.");
  }
  DrawColour clr(r, g, b);
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
  res.append(clr.get<0>());
  res.append(clr.get<1>());
  res.append(clr.get<2>());

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

void addMoleculeMetadata(const RDKit::MolDraw2DSVG &self, const RDKit::ROMol &m,
                         int confId) {
  self.addMoleculeMetadata(m, confId);
}

}  // namespace RDKit

BOOST_PYTHON_MODULE(rdMolDraw2D) {
  python::scope().attr("__doc__") =
      "Module containing a C++ implementation of 2D molecule drawing";
  std::string docString;

  python::class_<std::map<int, std::string>>("IntStringMap")
      .def(python::map_indexing_suite<std::map<int, std::string>, true>());

  docString = "Drawing options";
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
      .def_readwrite(
          "multipleBondOffset", &RDKit::MolDrawOptions::multipleBondOffset,
          "offset (in Angstroms) for the extra lines in a multiple bond")
      .def_readwrite("padding", &RDKit::MolDrawOptions::padding,
                     "fraction of empty space to leave around molecule")
      .def_readwrite(
          "bondLineWidth", &RDKit::MolDrawOptions::bondLineWidth,
          "if positive, this overrides the default line width for bonds")
      .def_readwrite("additionalAtomLabelPadding",
                     &RDKit::MolDrawOptions::additionalAtomLabelPadding,
                     "additional padding to leave around atom labels. "
                     "Expressed as a fraction of the font size.");
  docString = "Drawer abstract base class";
  python::class_<RDKit::MolDraw2D, boost::noncopyable>(
      "MolDraw2D", docString.c_str(), python::no_init)
      .def("SetFontSize", &RDKit::MolDraw2D::setFontSize,
           "change the default font size")
      .def("FontSize", &RDKit::MolDraw2D::fontSize, "get the default font size")
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
      .def("DrawString", &RDKit::MolDraw2D::drawString,
           (python::arg("self"), python::arg("string"), python::arg("pos")),
           "add text to the canvas")
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
      .def("drawOptions",
           (RDKit::MolDrawOptions & (RDKit::MolDraw2D::*)()) &
               RDKit::MolDraw2D::drawOptions,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1>>(),
           "Returns a modifiable version of the current drawing options");
  docString = "SVG molecule drawer";
  python::class_<RDKit::MolDraw2DSVG, python::bases<RDKit::MolDraw2D>,
                 boost::noncopyable>("MolDraw2DSVG", docString.c_str(),
                                     python::init<int, int>())
      .def(python::init<int, int, int, int>())
      .def("FinishDrawing", &RDKit::MolDraw2DSVG::finishDrawing,
           "add the last bits of SVG to finish the drawing")
      .def("AddMoleculeMetadata", RDKit::addMoleculeMetadata,
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
}
