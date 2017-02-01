//
//  Copyright (C) 2015 Greg Landrum
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
#ifdef RDK_CAIRO_BUILD
#include <cairo.h>
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#endif

namespace python = boost::python;

namespace RDKit {
namespace {
std::map<int, DrawColour> *pyDictToColourMap(python::object pyo) {
  std::map<int, DrawColour> *res = NULL;
  if (pyo) {
    res = new std::map<int, DrawColour>;
    python::dict tDict = python::extract<python::dict>(pyo);
    for (unsigned int i = 0;
         i < python::extract<unsigned int>(tDict.keys().attr("__len__")());
         ++i) {
      python::tuple tpl = python::extract<python::tuple>(tDict.values()[i]);
      float r = python::extract<float>(tpl[0]);
      float g = python::extract<float>(tpl[1]);
      float b = python::extract<float>(tpl[2]);
      DrawColour clr(r, g, b);
      (*res)[python::extract<int>(tDict.keys()[i])] = clr;
    }
  }
  return res;
}
std::map<int, double> *pyDictToDoubleMap(python::object pyo) {
  std::map<int, double> *res = NULL;
  if (pyo) {
    res = new std::map<int, double>;
    python::dict tDict = python::extract<python::dict>(pyo);
    for (unsigned int i = 0;
         i < python::extract<unsigned int>(tDict.keys().attr("__len__")());
         ++i) {
      double r = python::extract<double>(tDict.values()[i]);
      (*res)[python::extract<int>(tDict.keys()[i])] = r;
    }
  }
  return res;
}
}
void drawMoleculeHelper1(MolDraw2D &self, const ROMol &mol,
                         python::object highlight_atoms,
                         python::object highlight_atom_map,
                         python::object highlight_atom_radii, int confId,
                         std::string legend) {
  rdk_auto_ptr<std::vector<int> > highlightAtoms =
      pythonObjectToVect(highlight_atoms, static_cast<int>(mol.getNumAtoms()));
  std::map<int, DrawColour> *ham = pyDictToColourMap(highlight_atom_map);
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
  rdk_auto_ptr<std::vector<int> > highlightAtoms =
      pythonObjectToVect(highlight_atoms, static_cast<int>(mol.getNumAtoms()));
  rdk_auto_ptr<std::vector<int> > highlightBonds =
      pythonObjectToVect(highlight_bonds, static_cast<int>(mol.getNumBonds()));
  // FIX: support these
  std::map<int, DrawColour> *ham = pyDictToColourMap(highlight_atom_map);
  std::map<int, DrawColour> *hbm = pyDictToColourMap(highlight_bond_map);
  std::map<int, double> *har = pyDictToDoubleMap(highlight_atom_radii);

  self.drawMolecule(mol, legend, highlightAtoms.get(), highlightBonds.get(),
                    ham, hbm, har, confId);

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
  rdk_auto_ptr<std::vector<ROMol *> > mols = pythonObjectToVect<ROMol *>(pmols);
  // rdk_auto_ptr<std::vector<int> > highlightAtoms =
  //     pythonObjectToVect(highlight_atoms,
  //     static_cast<int>(mol.getNumAtoms()));
  // rdk_auto_ptr<std::vector<int> > highlightBonds =
  //     pythonObjectToVect(highlight_bonds,
  //     static_cast<int>(mol.getNumBonds()));
  // FIX: support these
  // std::map<int, DrawColour> *ham = pyDictToColourMap(highlight_atom_map);
  // std::map<int, DrawColour> *hbm = pyDictToColourMap(highlight_bond_map);
  // std::map<int, double> *har = pyDictToDoubleMap(highlight_atom_radii);
  //
  rdk_auto_ptr<std::vector<int> > confIds = pythonObjectToVect<int>(pconfIds);
  rdk_auto_ptr<std::vector<std::string> > legends =
      pythonObjectToVect<std::string>(plegends);

  self.drawMolecules(*mols, legends.get(), NULL, NULL, NULL, NULL, NULL,
                     confIds.get());
}

#ifdef RDK_CAIRO_BUILD
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
  RWMol *res = new RWMol(*m);
  MolDraw2DUtils::prepareMolForDrawing(*res, kekulize, addChiralHs, wedgeBonds,
                                       forceCoords);
  return static_cast<ROMol *>(res);
}
}

BOOST_PYTHON_MODULE(rdMolDraw2D) {
  python::scope().attr("__doc__") =
      "Module containing a C++ implementation of 2D molecule drawing";

  std::string docString;

  python::class_<std::map<int, std::string> >("IntStringMap")
      .def(python::map_indexing_suite<std::map<int, std::string>, true>());

  docString = "Drawing options";
  python::class_<RDKit::MolDrawOptions, boost::noncopyable>("MolDrawOptions",
                                                            docString.c_str())
      .def_readwrite("dummiesAreAttachments",
                     &RDKit::MolDrawOptions::dummiesAreAttachments)
      .def_readwrite("circleAtoms", &RDKit::MolDrawOptions::circleAtoms)
      //.def_readwrite("highlightColour",&RDKit::MolDrawOptions::highlightColour)
      .def_readwrite("atomLabels", &RDKit::MolDrawOptions::atomLabels,
                     "maps indices to atom labels")
      .def_readwrite("atomLabelDeuteriumTritium",
                     &RDKit::MolDrawOptions::atomLabelDeuteriumTritium,
                     "labels deuterium as D and tritium as T")
      .def_readwrite("continuousHighlight",
                     &RDKit::MolDrawOptions::continuousHighlight)
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
      .def("DrawString", &RDKit::MolDraw2D::drawString,
           (python::arg("self"), python::arg("string"), python::arg("pos")),
           "add text to the canvas")
      .def("drawOptions", (RDKit::MolDrawOptions & (RDKit::MolDraw2D::*)()) &
                              RDKit::MolDraw2D::drawOptions,
           python::return_internal_reference<
               1, python::with_custodian_and_ward_postcall<0, 1> >(),
           "Returns a modifiable version of the current drawing options");
  docString = "SVG molecule drawer";
  python::class_<RDKit::MolDraw2DSVG, python::bases<RDKit::MolDraw2D>,
                 boost::noncopyable>("MolDraw2DSVG", docString.c_str(),
                                     python::init<int, int>())
      .def(python::init<int, int, int, int>())
      .def("FinishDrawing", &RDKit::MolDraw2DSVG::finishDrawing,
           "add the last bits of SVG to finish the drawing")
      .def("GetDrawingText", &RDKit::MolDraw2DSVG::getDrawingText,
           "return the SVG");

#ifdef RDK_CAIRO_BUILD
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
}
