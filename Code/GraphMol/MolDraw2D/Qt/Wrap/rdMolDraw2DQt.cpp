//
//  Copyright (C) 2015-2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>

#include <GraphMol/MolDraw2D/Qt/MolDraw2DQt.h>
#include <QPainter>

namespace python = boost::python;

namespace RDKit {
MolDraw2DQt *moldrawFromQPainter(int width, int height, size_t ptr,
                                 int panelWidth, int panelHeight) {
  if (!ptr) {
    throw_value_error("QPainter pointer is null");
  }
  QPainter *qptr = reinterpret_cast<QPainter *>(ptr);
  return new MolDraw2DQt(width, height, qptr, panelWidth, panelHeight);
}

}  // namespace RDKit

BOOST_PYTHON_MODULE(rdMolDraw2DQt) {
  python::scope().attr("__doc__") =
      "Module containing a C++ implementation of 2D molecule drawing using Qt";

  python::scope().attr("rdkitQtVersion") = RDKit::rdkitQtVersion;

  std::string docString = "Qt molecule drawer";
  python::class_<RDKit::MolDraw2DQt, python::bases<RDKit::MolDraw2D>,
                 boost::noncopyable>("MolDraw2DQt", docString.c_str(),
                                     python::no_init);
  python::def("MolDraw2DFromQPainter_", RDKit::moldrawFromQPainter,
              (python::arg("width"), python::arg("height"),
               python::arg("pointer_to_QPainter"),
               python::arg("panelWidth") = -1, python::arg("panelHeight") = -1),
              "Returns a MolDraw2DQt instance set to use a QPainter.\nUse "
              "sip.unwrapinstance(qptr) to get the required pointer "
              "information. Please note that this is somewhat fragile.",
              python::return_value_policy<python::manage_new_object>());
}
