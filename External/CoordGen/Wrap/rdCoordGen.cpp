//
//  Copyright (C) 2017 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <boost/python/list.hpp>
#include <RDGeneral/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <CoordGen/CoordGen.h>

namespace python = boost::python;

namespace RDKit {

namespace {
void SetCoordMap(CoordGen::CoordGenParams *self, python::dict &coordMap) {
  self->coordMap.clear();
  python::list ks = coordMap.keys();
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(ks.attr("__len__")()); i++) {
    unsigned int id = python::extract<unsigned int>(ks[i]);
    self->coordMap[id] = python::extract<RDGeom::Point2D>(coordMap[id]);
  }
}
void addCoordsHelper(ROMol &mol, python::object &params) {
  CoordGen::CoordGenParams *ps = nullptr;
  if (params != python::object()) {
    ps = python::extract<CoordGen::CoordGenParams *>(params);
  }
  CoordGen::addCoords(mol, ps);
}
void SetTemplateMol(CoordGen::CoordGenParams *self, const ROMol *templ) {
  self->templateMol = templ;
}
void SetDefaultTemplateFileDir(const std::string &dir) {
  CoordGen::defaultParams.templateFileDir = dir;
}
}  // end of anonymous namespace

struct coordgen_wrapper {
  static void wrap() {
    std::string docString = "";

    python::class_<CoordGen::CoordGenParams>(
        "CoordGenParams", "Parameters controlling coordinate generation")
        .def(
            "SetCoordMap", SetCoordMap,
            "expects a dictionary of Point2D objects with template coordinates")
        .def("SetTemplateMol", SetTemplateMol,
             python::with_custodian_and_ward<1, 2>(),
             "sets a molecule to be used as the template")
        .def_readwrite("coordgenScaling",
                       &CoordGen::CoordGenParams::coordgenScaling,
                       "scaling factor for a single bond")
        .def_readwrite("dbg_useConstrained",
                       &CoordGen::CoordGenParams::dbg_useConstrained,
                       "for debugging use")
        .def_readwrite("dbg_useFixed", &CoordGen::CoordGenParams::dbg_useFixed,
                       "for debugging use")
        .def_readwrite("templateFileDir",
                       &CoordGen::CoordGenParams::templateFileDir,
                       "directory containing the templates.mae file")
        .def_readonly("sketcherBestPrecision",
                      &CoordGen::CoordGenParams::sketcherBestPrecision,
                      "highest quality (and slowest) precision setting")
        .def_readonly("sketcherStandardPrecision",
                      &CoordGen::CoordGenParams::sketcherStandardPrecision,
                      "standard quality precision setting, the default for the "
                      "coordgen project")
        .def_readonly("sketcherQuickPrecision",
                      &CoordGen::CoordGenParams::sketcherQuickPrecision,
                      "faster precision setting")
        .def_readonly(
            "sketcherCoarsePrecision",
            &CoordGen::CoordGenParams::sketcherCoarsePrecision,
            "\"coarse\" (fastest) precision setting, produces good-quality "
            "coordinates"
            " most of the time, this is the default setting for the RDKit")
        .def_readwrite("minimizerPrecision",
                       &CoordGen::CoordGenParams::minimizerPrecision,
                       "controls sketcher precision")
        .def_readwrite(
            "treatNonterminalBondsToMetalAsZOBs",
            &CoordGen::CoordGenParams::treatNonterminalBondsToMetalAsZeroOrder);
    python::def("SetDefaultTemplateFileDir", SetDefaultTemplateFileDir);
    docString =
        "Add 2D coordinates.\n"
        "ARGUMENTS:\n"
        "   - mol: molecule to modify\n"
        "   - params: (optional) parameters controlling the coordinate "
        "generation\n"
        "\n";
    python::def("AddCoords", addCoordsHelper,
                (python::arg("mol"), python::arg("params") = python::object()),
                docString.c_str());
  }
};

}  // end of namespace RDKit
BOOST_PYTHON_MODULE(rdCoordGen) {
  python::scope().attr("__doc__") =
      "Module containing interface to the CoordGen library.";
  RDKit::coordgen_wrapper::wrap();
}
