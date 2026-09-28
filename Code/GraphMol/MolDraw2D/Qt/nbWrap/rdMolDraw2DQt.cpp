//
//  Copyright (C) 2015-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>

#include <GraphMol/MolDraw2D/Qt/MolDraw2DQt.h>
#include <QPainter>

namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
static MolDraw2DQt *moldrawFromQPainter(int width, int height, size_t ptr,
                                        int panelWidth, int panelHeight) {
  if (!ptr) {
    throw nb::value_error("QPainter pointer is null");
  }
  QPainter *qptr = reinterpret_cast<QPainter *>(ptr);
  return new MolDraw2DQt(width, height, qptr, panelWidth, panelHeight);
}

}  // namespace RDKit

NB_MODULE(rdMolDraw2DQt, m) {
  m.doc() =
      "Module containing a C++ implementation of 2D molecule drawing using Qt";

  m.attr("rdkitQtVersion") = RDKit::rdkitQtVersion;

  nb::class_<RDKit::MolDraw2DQt, RDKit::MolDraw2D>(m, "MolDraw2DQt",
                                                    "Qt molecule drawer",
                                                    nb::dynamic_attr());

  m.def(
      "MolDraw2DFromQPainter_", &RDKit::moldrawFromQPainter,
      "width"_a, "height"_a, "pointer_to_QPainter"_a,
      "panelWidth"_a = -1, "panelHeight"_a = -1,
      R"DOC(Returns a MolDraw2DQt instance set to use a QPainter.
Use sip.unwrapinstance(qptr) to get the required pointer information.
Please note that this is somewhat fragile.)DOC",
      nb::rv_policy::take_ownership);
}
