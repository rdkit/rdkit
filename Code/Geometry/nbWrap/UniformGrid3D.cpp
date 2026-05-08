//
//  Copyright (C) 2026 Greg Landrum
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
#include <DataStructs/DiscreteValueVect.h>
#include <Geometry/point.h>
#include <Geometry/UniformGrid3D.h>
#include <Geometry/GridUtils.h>
namespace nb = nanobind;
using namespace nb::literals;

using namespace RDKit;

namespace RDGeom {

#if 0
UniformGrid3D *makeUnformGrid3D(double dimX, double dimY, double dimZ,
                                double spacing = 0.5,
                                DiscreteValueVect::DiscreteValueType valType =
                                    DiscreteValueVect::TWOBITVALUE,
                                const Point3D *offSet = nullptr) {
  auto *grd = new UniformGrid3D(dimX, dimY, dimZ, spacing, valType, offSet);
  return grd;
}

int getValPoint(const UniformGrid3D &grid, const Point3D &pt) {
  return grid.getVal(pt);
}

int getValIndex(const UniformGrid3D &grid, unsigned int id) {
  return grid.getVal(id);
}

void setValIndex(UniformGrid3D &grid, unsigned int id, unsigned int val) {
  grid.setVal(id, val);
}

void setValPoint(UniformGrid3D &grid, const Point3D &pt, unsigned int val) {
  grid.setVal(pt, val);
}

python::tuple computeGridCentroidWrap(const UniformGrid3D &grid,
                                      const Point3D &pt, double windowRadius) {
  double weightSum;
  Point3D centroid = computeGridCentroid(grid, pt, windowRadius, weightSum);
  return python::make_tuple(weightSum, centroid);
}
python::tuple getGridIndicesWrap(const UniformGrid3D &grid, unsigned int idx) {
     unsigned int xi, yi, zi;
     grid.getGridIndices(idx, xi, yi, zi);
     python::list pyRes;
     pyRes.append(xi);
     pyRes.append(yi);
     pyRes.append(zi);
     return python::tuple(pyRes);
}
python::tuple findGridTerminalPointsWrap(const UniformGrid3D &grid,
     double windowRadius,
     double inclusionFraction) {
          std::vector<Point3D> res =
          findGridTerminalPoints(grid, windowRadius, inclusionFraction);
          python::list pyRes;
          for (auto &re : res) {
               pyRes.append(re);
          }
          return python::tuple(pyRes);
     }
#endif

std::string uGridClassDoc =
    "Class to represent a uniform three-dimensional\n\
    cubic grid. Each grid point can store a positive integer value. For the sake\n\
    of efficiency these value can either be binary or fit in 2, 4, 8 or 16 bits\n";

struct uGrid3D_wrapper {
  static void wrap(nb::module_ &m) {
    nb::class_<UniformGrid3D>(m, "UniformGrid3D")
        .def(nb::init<double, double, double, double,
                      DiscreteValueVect::DiscreteValueType, const Point3D *>(),
             "dimX"_a, "dimY"_a, "dimZ"_a, "spacing"_a = 0.5,
             "valType"_a = DiscreteValueVect::TWOBITVALUE,
             "offSet"_a = (const Point3D *)nullptr,
             "Constructor for a UniformGrid3D object")
        .def(nb::init<std::string>(), "pkl"_a, "pickle constructor")
        .def("__init__",
             [](UniformGrid3D *t, nb::bytes b) {
               new (t) UniformGrid3D(
                   std::string(static_cast<const char *>(b.data()),
                               static_cast<size_t>(b.size())));
             })
        .def("GetGridPointIndex", &UniformGrid3D::getGridPointIndex, "point"_a,
             "Get the index to the grid point closest to the specified point")
        .def("GetGridIndex", &UniformGrid3D::getGridIndex, "xi"_a, "yi"_a,
             "zi"_a,
             "Get the index to the grid point with the three integer indices "
             "provided")
        .def(
            "GetGridIndices",
            [](const UniformGrid3D &grid, unsigned int idx) {
              unsigned int xi, yi, zi;
              grid.getGridIndices(idx, xi, yi, zi);
              return std::make_tuple(xi, yi, zi);
            },
            "idx"_a, "Returns the integer indices of the grid index provided.")
        .def("GetValPoint",
             nb::overload_cast<const Point3D &>(&UniformGrid3D::getVal,
                                                nb::const_),
             "pt"_a, "Get the value at the closest grid point")
        .def(
            "GetVal",
            nb::overload_cast<unsigned int>(&UniformGrid3D::getVal, nb::const_),
            "id"_a, "Get the value at the specified grid index")
        .def("SetVal",
             nb::overload_cast<unsigned int, unsigned int>(
                 &UniformGrid3D::setVal),
             "id"_a, "val"_a, "Set the value at the specified grid index")
        .def("SetValPoint",
             nb::overload_cast<const Point3D &, unsigned int>(
                 &UniformGrid3D::setVal),
             "pt"_a, "val"_a,
             "Set the value at grid point closest to the specified point")
        .def("GetGridPointLoc", &UniformGrid3D::getGridPointLoc, "pointId"_a,
             "Get the location of the specified grid point")
        .def("GetSize", &UniformGrid3D::getSize,
             "Get the size of the grid (number of grid points)")
        .def("GetNumX", &UniformGrid3D::getNumX,
             "Get the number of grid points along x-axis")
        .def("GetNumY", &UniformGrid3D::getNumY,
             "Get the number of grid points along y-axis")
        .def("GetNumZ", &UniformGrid3D::getNumZ,
             "Get the number of grid points along z-axis")
        .def("GetOffset", &UniformGrid3D::getOffset, nb::rv_policy::copy,
             "Get the location of the center of the grid")
        .def("GetSpacing", &UniformGrid3D::getSpacing, "Get the grid spacing")
        .def("GetOccupancyVect", &UniformGrid3D::getOccupancyVect,
             nb::rv_policy::reference_internal,
             "Get the occupancy vector for the grid")
        .def("CompareParams", &UniformGrid3D::compareParams, "other"_a,
             "Compare the parameters between two grid object")
        .def("SetSphereOccupancy", &UniformGrid3D::setSphereOccupancy,
             "center"_a, "radius"_a, "stepSize"_a, "maxLayers"_a = -1,
             "ignoreOutOfBound"_a = true,
             "Set the occupancy on the grid for a sphere or specified radius\n"
             " and multiple layers around this sphere, with decreasing values "
             "of \n"
             "occupancy\n")

        .def(nb::self += nb::self)
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
        .def(nb::self &=
             nb::self)  // clang warns incorrectly on these constructs
        .def(nb::self |= nb::self)
        .def(nb::self -= nb::self)
#ifdef __clang__
#pragma GCC diagnostic pop
#endif
        .def("__getstate__",
             [](const UniformGrid3D &grd) {
               const auto pkl = grd.toString();
               return std::make_tuple(nb::bytes(pkl.c_str(), pkl.size()));
             })
        .def("__setstate__",
             [](UniformGrid3D &grd, const std::tuple<nb::bytes> &state) {
               std::string pkl = std::string(
                   static_cast<const char *>(std::get<0>(state).data()),
                   static_cast<size_t>(std::get<0>(state).size()));
               new (&grd) UniformGrid3D(pkl);
             })
        .doc() = uGridClassDoc.c_str();

    m.def("WriteGridToFile", writeGridToFile, "grid"_a, "filename"_a,
          "Write the grid to a grid file");

    m.def("TverskyIndex", tverskyIndex<UniformGrid3D>, "grid1"_a, "grid2"_a,
          "alpha"_a, "beta"_a,
          "Compute the tversky index between two grid objects");
    m.def("TanimotoDistance", tanimotoDistance<UniformGrid3D>, "grid1"_a,
          "grid2"_a, "Compute the tanimoto distance between two grid objects");
    m.def("ProtrudeDistance", protrudeDistance<UniformGrid3D>, "grid1"_a,
          "grid2"_a, "Compute the protrude distance between two grid objects");
    m.def(
        "ComputeGridCentroid",
        [](const UniformGrid3D &grid, const Point3D &pt, double windowRadius) {
          double weightSum;
          Point3D centroid =
              computeGridCentroid(grid, pt, windowRadius, weightSum);
          return std::make_tuple(weightSum, centroid);
        },
        "Compute the grid point at the center of sphere around a Point3D");
    m.def("FindGridTerminalPoints", &findGridTerminalPoints, "grid"_a,
          "windowRadius"_a, "inclusionFraction"_a,
          "Find a grid's terminal points (defined in the subshape algorithm).");
  }
};
}  // namespace RDGeom

void wrap_uniformGrid(nb::module_ &m) { RDGeom::uGrid3D_wrapper::wrap(m); }
