// $Id$
//
//  Copyright (C) 2005 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>

#include <RDBoost/Wrap.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <DataStructs/DiscreteValueVect.h>
#include <Geometry/point.h>
#include <Geometry/UniformGrid3D.h>
#include <Geometry/GridUtils.h>
namespace python = boost::python;

using namespace RDKit;

namespace RDGeom {
struct ug3d_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const UniformGrid3D &self) {
    std::string res = self.toString();
    python::object retval = python::object(
        python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
    return python::make_tuple(retval);
  };
};

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

std::string uGridClassDoc =
    "Class to represent a uniform three-dimensional\n\
    cubic grid. Each grid point can store a poisitive integer value. For the sake\n\
    of efficiency these value can either be binary, fit in 2, 4, 8 or 16 bits\n";

struct uGrid3D_wrapper {
  static void wrap() {
    python::class_<UniformGrid3D>(
        "UniformGrid3D_", uGridClassDoc.c_str(),
        python::init<std::string>(python::args("self", "pkl"),
                                  "pickle constructor"))
        .def("GetGridPointIndex", &UniformGrid3D::getGridPointIndex,
             python::args("self", "point"),
             "Get the index to the grid point closest to the specified point")
        .def("GetGridIndex", &UniformGrid3D::getGridIndex,
             python::args("self", "xi", "yi", "zi"),
             "Get the index to the grid point with the three integer indices "
             "provided")
        .def("GetGridIndices", &getGridIndicesWrap, python::args("self", "idx"),
             "Returns the integer indices of the grid index provided.")
        .def("GetValPoint", getValPoint, python::args("self", "pt"),
             "Get the value at the closest grid point")
        .def("GetVal", getValIndex, python::args("self", "id"),
             "Get the value at the specified grid point")
        .def("SetVal", setValIndex, python::args("self", "id", "val"),
             "Set the value at the specified grid point")
        .def("SetValPoint", setValPoint, python::args("self", "pt", "val"),
             "Set the value at grid point closest to the specified point")
        .def("GetGridPointLoc", &UniformGrid3D::getGridPointLoc,
             python::args("self", "pointId"),
             "Get the location of the specified grid point")
        .def("GetSize", &UniformGrid3D::getSize, python::args("self"),
             "Get the size of the grid (number of grid points)")
        .def("GetNumX", &UniformGrid3D::getNumX, python::args("self"),
             "Get the number of grid points along x-axis")
        .def("GetNumY", &UniformGrid3D::getNumY, python::args("self"),
             "Get the number of grid points along y-axis")
        .def("GetNumZ", &UniformGrid3D::getNumZ, python::args("self"),
             "Get the number of grid points along z-axis")
        .def("GetOffset", &UniformGrid3D::getOffset,
             python::return_value_policy<python::copy_const_reference>(),
             python::args("self"), "Get the location of the center of the grid")
        .def("GetSpacing", &UniformGrid3D::getSpacing, python::args("self"),
             "Get the grid spacing")
        .def("GetOccupancyVect", &UniformGrid3D::getOccupancyVect,
             python::return_value_policy<python::reference_existing_object>(),
             python::args("self"), "Get the occupancy vector for the grid")
        .def("CompareParams", &UniformGrid3D::compareParams,
             python::args("self", "other"),
             "Compare the parameters between two grid object")
        .def("SetSphereOccupancy", &UniformGrid3D::setSphereOccupancy,
             (python::arg("self"), python::arg("center"), python::arg("radius"),
              python::arg("stepSize"), python::arg("maxLayers") = -1,
              python::arg("ignoreOutOfBound") = true),
             "Set the occupancy on the grid for a sphere or specified radius\n"
             " and multiple layers around this sphere, with decreasing values "
             "of \n"
             "occupancy\n")

        .def(python::self += python::self)
#ifdef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
        .def(python::self &=
             python::self)  // clang warns incorrectly on these constructs
        .def(python::self |= python::self)
        .def(python::self -= python::self)
#ifdef __clang__
#pragma GCC diagnostic pop
#endif

        .def_pickle(RDGeom::ug3d_pickle_suite())

        ;

    python::def("UniformGrid3D", makeUnformGrid3D,
                (python::arg("dimX"), python::arg("dimY"), python::arg("dimZ"),
                 python::arg("spacing") = 0.5,
                 python::arg("valType") = DiscreteValueVect::TWOBITVALUE,
                 python::arg("offSet") = (const Point3D *)nullptr),
                "Faking the constructor",
                python::return_value_policy<python::manage_new_object>());

    python::def("WriteGridToFile", writeGridToFile,
                python::args("grid", "filename"),
                "Write the grid to a grid file");

    python::def("TverskyIndex", tverskyIndex<UniformGrid3D>,
                python::args("grid1", "grid2", "alpha", "beta"),
                "Compute the tversky index between two grid objects");
    python::def("TanimotoDistance", tanimotoDistance<UniformGrid3D>,
                python::args("grid1", "grid2"),
                "Compute the tanimoto distance between two grid objects");
    python::def("ProtrudeDistance", protrudeDistance<UniformGrid3D>,
                python::args("grid1", "grid2"),
                "Compute the protrude distance between two grid objects");
    python::def(
        "ComputeGridCentroid", computeGridCentroidWrap,
        python::args("grid", "pt", "windowRadius"),
        "Compute the grid point at the center of sphere around a Point3D");
    python::def(
        "FindGridTerminalPoints", findGridTerminalPointsWrap,
        python::args("grid", "windowRadius", "inclusionFraction"),
        "Find a grid's terminal points (defined in the subshape algorithm).");
  }
};
}  // namespace RDGeom

void wrap_uniformGrid() { RDGeom::uGrid3D_wrapper::wrap(); }
