//
//  Copyright (C) 2017-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <GraphMol/RDKitBase.h>
#include <CoordGen/CoordGen.h>
#include <RDBoost/Wrap_nb.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {

namespace {

void SetCoordMap(CoordGen::CoordGenParams *self, nb::dict coordMap) {
  self->coordMap.clear();
  for (auto item : coordMap) {
    unsigned int id = nb::cast<unsigned int>(item.first);
    self->coordMap[id] = nb::cast<RDGeom::Point2D>(item.second);
  }
}

void addCoordsHelper(ROMol &mol, nb::object params) {
  CoordGen::CoordGenParams *ps = nullptr;
  if (!params.is_none()) {
    ps = nb::cast<CoordGen::CoordGenParams *>(params);
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

}  // end of namespace RDKit

NB_MODULE(rdCoordGen, m) {
  m.doc() = "Module containing interface to the CoordGen library.";

  nb::class_<RDKit::CoordGen::CoordGenParams>(m, "CoordGenParams",
                                               "Parameters controlling coordinate generation")
      .def(nb::init<>())
      .def("SetCoordMap", &RDKit::SetCoordMap, "coordMap"_a,
           "expects a dictionary of Point2D objects with template coordinates")
      .def("SetTemplateMol", &RDKit::SetTemplateMol, nb::keep_alive<1, 2>(),
           "templ"_a, "sets a molecule to be used as the template")
      .def_rw("coordgenScaling",
              &RDKit::CoordGen::CoordGenParams::coordgenScaling,
              "scaling factor for a single bond")
      .def_rw("dbg_useConstrained",
              &RDKit::CoordGen::CoordGenParams::dbg_useConstrained,
              "for debugging use")
      .def_rw("dbg_useFixed", &RDKit::CoordGen::CoordGenParams::dbg_useFixed,
              "for debugging use")
      .def_rw("templateFileDir",
              &RDKit::CoordGen::CoordGenParams::templateFileDir,
              "directory containing the templates.mae file")
      .def_ro("sketcherBestPrecision",
              &RDKit::CoordGen::CoordGenParams::sketcherBestPrecision,
              "highest quality (and slowest) precision setting")
      .def_ro("sketcherStandardPrecision",
              &RDKit::CoordGen::CoordGenParams::sketcherStandardPrecision,
              "standard quality precision setting, the default for the coordgen project")
      .def_ro("sketcherQuickPrecision",
              &RDKit::CoordGen::CoordGenParams::sketcherQuickPrecision,
              "faster precision setting")
      .def_ro("sketcherCoarsePrecision",
              &RDKit::CoordGen::CoordGenParams::sketcherCoarsePrecision,
              R"DOC("coarse" (fastest) precision setting, produces good-quality
coordinates most of the time, this is the default setting for the RDKit)DOC")
      .def_rw("minimizerPrecision",
              &RDKit::CoordGen::CoordGenParams::minimizerPrecision,
              "controls sketcher precision")
      .def_rw("treatNonterminalBondsToMetalAsZOBs",
              &RDKit::CoordGen::CoordGenParams::treatNonterminalBondsToMetalAsZeroOrder)
      .def("__setattr__", &safeSetattr);

  m.def("SetDefaultTemplateFileDir", &RDKit::SetDefaultTemplateFileDir,
        "dir"_a);

  m.def(
      "AddCoords", &RDKit::addCoordsHelper,
      "mol"_a, "params"_a = nb::none(),
      R"DOC(Add 2D coordinates.
ARGUMENTS:
   - mol: molecule to modify
   - params: (optional) parameters controlling the coordinate generation
)DOC");
}
