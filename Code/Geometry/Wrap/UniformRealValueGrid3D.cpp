//
//  Copyright (c) 2014-2024, Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <boost/python.hpp>
#include <RDBoost/Wrap.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <DataStructs/RealValueVect.h>
#include <Geometry/point.h>
#include <Geometry/UniformRealValueGrid3D.h>
#include <Geometry/GridUtils.h>
namespace python = boost::python;

using namespace RDKit;

namespace RDGeom {
struct urvg3d_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const UniformRealValueGrid3D &self) {
    auto res = self.toString();
    python::object retval(
        python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
    return python::make_tuple(retval);
  }
};

UniformRealValueGrid3D *makeUniformRealValueGrid3D(
    double dimX, double dimY, double dimZ, double spacing = 0.5,
    const Point3D *offSet = nullptr) {
  UniformRealValueGrid3D *grd =
      new UniformRealValueGrid3D(dimX, dimY, dimZ, spacing, offSet);
  return grd;
}

double getValPoint(const UniformRealValueGrid3D &grid, const Point3D &pt) {
  return grid.getVal(pt);
}

double getValIndex(const UniformRealValueGrid3D &grid, unsigned int id) {
  return grid.getVal(id);
}

void setValIndex(UniformRealValueGrid3D &grid, unsigned int id, double val) {
  grid.setVal(id, val);
}

void setValPoint(UniformRealValueGrid3D &grid, const Point3D &pt, double val) {
  grid.setVal(pt, val);
}

python::tuple getGridIndicesWrap(const UniformRealValueGrid3D &grid,
                                 unsigned int idx) {
  unsigned int xi, yi, zi;
  grid.getGridIndices(idx, xi, yi, zi);
  return python::make_tuple(xi, yi, zi);
}

std::string urvGridClassDoc =
    "Class to represent a uniform three-dimensional\n\
    cubic grid. Each grid point can store a floating point value. \n";

struct urvGrid3D_wrapper {
  static void wrap() {
    python::class_<UniformRealValueGrid3D>(
        "UniformRealValueGrid3D", urvGridClassDoc.c_str(),
        python::init<std::string>("Pickle constructor"))
        .def(python::init<>("Default constructor"))
        .def(python::init<RDGeom::UniformRealValueGrid3D>("Copy constructor"))
        .def("__init__",
             python::make_constructor(
                 makeUniformRealValueGrid3D, python::default_call_policies(),
                 (python::arg("dimX"), python::arg("dimY"), python::arg("dimZ"),
                  python::arg("spacing") = 0.5,
                  python::arg("offSet") = python::object())),
             "Constructor")
        .def("GetGridPointIndex", &UniformRealValueGrid3D::getGridPointIndex,
             "Get the index to the grid point closest to the specified point")
        .def(
            "GetGridIndex", &UniformRealValueGrid3D::getGridIndex,
            "Get the index to the grid point with the three integer indices provided")
        .def("GetGridIndices", &getGridIndicesWrap,
             "Returns the integer indices of the grid index provided.")
        .def("GetValPoint", getValPoint,
             "Get the value at the closest grid point")
        .def("GetVal", getValIndex, "Get the value at the specified grid point")
        .def("SetVal", setValIndex, "Set the value at the specified grid point")
        .def("SetValPoint", setValPoint,
             "Set the value at grid point closest to the specified point")
        .def("GetGridPointLoc", &UniformRealValueGrid3D::getGridPointLoc,
             "Get the location of the specified grid point")
        .def("GetSize", &UniformRealValueGrid3D::getSize,
             "Get the size of the grid (number of grid points)")
        .def("GetNumX", &UniformRealValueGrid3D::getNumX,
             "Get the number of grid points along x-axis")
        .def("GetNumY", &UniformRealValueGrid3D::getNumY,
             "Get the number of grid points along y-axis")
        .def("GetNumZ", &UniformRealValueGrid3D::getNumZ,
             "Get the number of grid points along z-axis")
        .def("GetOffset", &UniformRealValueGrid3D::getOffset,
             python::return_value_policy<python::copy_const_reference>(),
             "Get the location of the center of the grid")
        .def("GetSpacing", &UniformRealValueGrid3D::getSpacing,
             "Get the grid spacing")
        .def("GetOccupancyVect", &UniformRealValueGrid3D::getOccupancyVect,
             python::return_value_policy<python::reference_existing_object>(),
             "Get the occupancy vector for the grid")
        .def("CompareVectors", &UniformRealValueGrid3D::compareVectors,
             "Compare the vector values between two grid objects.")
        .def("CompareParams", &UniformRealValueGrid3D::compareParams,
             "Compare the parameters between two grid object.")
        .def("CompareGrids", &UniformRealValueGrid3D::compareGrids,
             "Compare the parameters and values between two grid objects.")
        .def(python::self &= python::self)
        .def(python::self |= python::self)
        .def(python::self += python::self)
        .def(python::self -= python::self)
        .def_pickle(RDGeom::urvg3d_pickle_suite());
  }
};
}  // namespace RDGeom

void wrap_uniformrealvalueGrid() { RDGeom::urvGrid3D_wrapper::wrap(); }
