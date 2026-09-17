//
//  Copyright (c) 2026 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>

#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <DataStructs/RealValueVect.h>
#include <Geometry/point.h>
#include <Geometry/UniformRealValueGrid3D.h>
#include <Geometry/GridUtils.h>
namespace nb = nanobind;
using namespace nb::literals;

using namespace RDKit;

namespace RDGeom {
#if 0
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
#endif

std::string urvGridClassDoc =
    "Class to represent a uniform three-dimensional\n\
    cubic grid. Each grid point can store a floating point value. \n";

struct urvGrid3D_wrapper {
  static void wrap(nb::module_ &m) {
    nb::class_<UniformRealValueGrid3D>(m, "UniformRealValueGrid3D")
        .def(nb::init<>(), "Default constructor")
        .def(nb::init<RDGeom::UniformRealValueGrid3D>(), "Copy constructor")
        .def(
            nb::init<double, double, double, double, const RDGeom::Point3D *>(),
            "dimX"_a, "dimY"_a, "dimZ"_a, "spacing"_a = 0.5,
            "offSet"_a = (const RDGeom::Point3D *)nullptr, "Constructor")
        .def("__init__",
             [](UniformRealValueGrid3D *t, nb::bytes b) {
               new (t) UniformRealValueGrid3D(
                   std::string(static_cast<const char *>(b.data()),
                               static_cast<size_t>(b.size())));
             })
        .def("GetGridPointIndex", &UniformRealValueGrid3D::getGridPointIndex,
             "Get the index to the grid point closest to the specified point")
        .def(
            "GetGridIndex", &UniformRealValueGrid3D::getGridIndex,
            "Get the index to the grid point with the three integer indices provided")
        .def(
            "GetGridIndices",
            [](const UniformRealValueGrid3D &grid, unsigned int idx) {
              unsigned int xi, yi, zi;
              grid.getGridIndices(idx, xi, yi, zi);
              return std::make_tuple(xi, yi, zi);
            },
            "idx"_a, "Returns the integer indices of the grid index provided.")
        .def("GetValPoint",
             nb::overload_cast<const Point3D &>(&UniformRealValueGrid3D::getVal,
                                                nb::const_),
             "pt"_a, "Get the value at the closest grid point")
        .def("GetVal",
             nb::overload_cast<unsigned int>(&UniformRealValueGrid3D::getVal,
                                             nb::const_),
             "id"_a, "Get the value at the specified grid index")
        .def("SetVal",
             nb::overload_cast<unsigned int, double>(
                 &UniformRealValueGrid3D::setVal),
             "id"_a, "val"_a, "Set the value at the specified grid index")
        .def("SetValPoint",
             nb::overload_cast<const Point3D &, double>(
                 &UniformRealValueGrid3D::setVal),
             "pt"_a, "val"_a,
             "Set the value at grid point closest to the specified point")
        .def("GetGridPointLoc", &UniformRealValueGrid3D::getGridPointLoc,
             "pointId"_a, "Get the location of the specified grid point")
        .def("GetSize", &UniformRealValueGrid3D::getSize,
             "Get the size of the grid (number of grid points)")
        .def("GetNumX", &UniformRealValueGrid3D::getNumX,
             "Get the number of grid points along x-axis")
        .def("GetNumY", &UniformRealValueGrid3D::getNumY,
             "Get the number of grid points along y-axis")
        .def("GetNumZ", &UniformRealValueGrid3D::getNumZ,
             "Get the number of grid points along z-axis")
        .def("GetOffset", &UniformRealValueGrid3D::getOffset,
             nb::rv_policy::copy, "Get the location of the center of the grid")
        .def("GetSpacing", &UniformRealValueGrid3D::getSpacing,
             "Get the grid spacing")
        .def("GetOccupancyVect", &UniformRealValueGrid3D::getOccupancyVect,
             nb::rv_policy::reference_internal,
             "Get the occupancy vector for the grid")
        .def("CompareVectors", &UniformRealValueGrid3D::compareVectors,
             "Compare the vector values between two grid objects.")
        .def("CompareParams", &UniformRealValueGrid3D::compareParams,
             "Compare the parameters between two grid object.")
        .def("CompareGrids", &UniformRealValueGrid3D::compareGrids,
             "Compare the parameters and values between two grid objects.")
        .def(nb::self &= nb::self)
        .def(nb::self |= nb::self)
        .def(nb::self += nb::self)
        .def(nb::self -= nb::self)
        .def("__getstate__",
             [](const UniformRealValueGrid3D &grid) {
               const auto pkl = grid.toString();
               return std::make_tuple(nb::bytes(pkl.data(), pkl.size()));
             })
        .def("__setstate__",
             [](UniformRealValueGrid3D &grid,
                const std::tuple<nb::bytes> &state) {
               std::string pkl = std::string(
                   static_cast<const char *>(std::get<0>(state).data()),
                   static_cast<size_t>(std::get<0>(state).size()));
               new (&grid) UniformRealValueGrid3D(pkl);
             })
        .doc() = urvGridClassDoc.c_str();
  }
};
}  // namespace RDGeom

void wrap_uniformrealvalueGrid(nb::module_ &m) {
  RDGeom::urvGrid3D_wrapper::wrap(m);
}
