//
//   Copyright (C) 2005-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/tuple.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>

#include <Geometry/Transform3D.h>
#include <Geometry/UniformGrid3D.h>

#include <GraphMol/ShapeHelpers/ShapeEncoder.h>
#include <GraphMol/ShapeHelpers/ShapeUtils.h>
#include <DataStructs/DiscreteValueVect.h>
#include <Geometry/point.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

static void copyTransform(
    const nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig>
        &transMat,
    RDGeom::Transform3D &trans) {
  unsigned int nrows = transMat.shape(0);
  unsigned int ncols = transMat.shape(1);
  if ((nrows != 4) || (ncols != 4)) {
    throw nb::value_error("The transform has to be square matrix, of size 4x4");
  }

  unsigned int dSize = nrows * nrows;
  const double *inData = transMat.data();
  double *tData = trans.getData();
  memcpy(static_cast<void *>(tData), static_cast<const void *>(inData),
         dSize * sizeof(double));
}

}  // namespace

NB_MODULE(rdShapeHelpers, m) {
  m.doc() = R"DOC(Module containing functions to encode and compare the shapes of molecules)DOC";

  m.def(
      "EncodeShape",
      [](const ROMol &mol, RDGeom::UniformGrid3D &grid, int confId,
         nb::object trans, double vdwScale, double stepSize, int maxLayers,
         bool ignoreHs) {
        if (!trans.is_none()) {
          auto transMat = nb::cast<
              nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig>>(
              trans);
          RDGeom::Transform3D ctrans;
          copyTransform(transMat, ctrans);
          MolShapes::EncodeShape(mol, grid, confId, &ctrans, vdwScale, stepSize,
                                 maxLayers, ignoreHs);
        } else {
          MolShapes::EncodeShape(mol, grid, confId, nullptr, vdwScale, stepSize,
                                 maxLayers, ignoreHs);
        }
      },
      "mol"_a, "grid"_a, "confId"_a = -1, "trans"_a = nb::none(),
      "vdwScale"_a = 0.8, "stepSize"_a = 0.25, "maxLayers"_a = -1,
      "ignoreHs"_a = true,
      R"DOC(Encode the shape of a molecule (one of its conformer) onto a grid

ARGUMENTS:

   - mol : the molecule of interest
   - grid : grid onto which the encoding is written
   - confId : id of the conformation of interest on mol (defaults to the first one)
   - trans : any transformation that needs to be used to encode onto the grid (note the molecule remains unchanged)
   - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
                used in the encoding - grid points inside this sphere carry the maximum occupancy
   - setpSize : thickness of the layers outside the base radius, the occupancy value is decreased
                from layer to layer from the maximum value
   - maxLayers : the maximum number of layers - defaults to the number of bits
                 used per grid point - e.g. two bits per grid point will allow 3 layers
   - ignoreHs : when set, the contribution of Hs to the shape will be ignored
)DOC");

  m.def(
      "ShapeTverskyIndex",
      nb::overload_cast<const ROMol &, const ROMol &, double, double, int, int,
                        double, DiscreteValueVect::DiscreteValueType, double,
                        double, int, bool>(&MolShapes::tverskyIndex),
      "mol1"_a, "mol2"_a, "alpha"_a, "beta"_a, "confId1"_a = -1,
      "confId2"_a = -1, "gridSpacing"_a = 0.5,
      "bitsPerPoint"_a = DiscreteValueVect::TWOBITVALUE,
      "vdwScale"_a = 0.8, "stepSize"_a = 0.25, "maxLayers"_a = -1,
      "ignoreHs"_a = true,
      R"DOC(Compute the shape tversky index between two molecule based on a predefined alignment

ARGUMENTS:
   - mol1 : The first molecule of interest
   - mol2 : The second molecule of interest
   - alpha : first parameter of the Tversky index
   - beta : second parameter of the Tversky index
   - confId1 : Conformer in the first molecule (defaults to first conformer)
   - confId2 : Conformer in the second molecule (defaults to first conformer)
   - gridSpacing : resolution of the grid used to encode the molecular shapes
   - bitsPerPoint : number of bits used to encode the occupancy at each grid point
                         defaults to two bits per grid point
   - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
               used in the encoding - grid points inside this sphere carry the maximum occupancy
   - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased
                from layer to layer from the maximum value
   - maxLayers : the maximum number of layers - defaults to the number of bits
                 used per grid point - e.g. two bits per grid point will allow 3 layers
   - ignoreHs : when set, the contribution of Hs to the shape will be ignored
)DOC");

  m.def(
      "ShapeTanimotoDist",
      nb::overload_cast<const ROMol &, const ROMol &, int, int, double,
                        DiscreteValueVect::DiscreteValueType, double, double,
                        int, bool>(&MolShapes::tanimotoDistance),
      "mol1"_a, "mol2"_a, "confId1"_a = -1, "confId2"_a = -1,
      "gridSpacing"_a = 0.5,
      "bitsPerPoint"_a = DiscreteValueVect::TWOBITVALUE,
      "vdwScale"_a = 0.8, "stepSize"_a = 0.25, "maxLayers"_a = -1,
      "ignoreHs"_a = true,
      R"DOC(Compute the shape tanimoto distance between two molecule based on a predefined alignment

ARGUMENTS:
   - mol1 : The first molecule of interest
   - mol2 : The second molecule of interest
   - confId1 : Conformer in the first molecule (defaults to first conformer)
   - confId2 : Conformer in the second molecule (defaults to first conformer)
   - gridSpacing : resolution of the grid used to encode the molecular shapes
   - bitsPerPoint : number of bits used to encode the occupancy at each grid point
                         defaults to two bits per grid point
   - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
               used in the encoding - grid points inside this sphere carry the maximum occupancy
   - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased
                from layer to layer from the maximum value
   - maxLayers : the maximum number of layers - defaults to the number of bits
                 used per grid point - e.g. two bits per grid point will allow 3 layers
   - ignoreHs : when set, the contribution of Hs to the shape will be ignored
)DOC");

  m.def(
      "ShapeProtrudeDist",
      nb::overload_cast<const ROMol &, const ROMol &, int, int, double,
                        DiscreteValueVect::DiscreteValueType, double, double,
                        int, bool, bool>(&MolShapes::protrudeDistance),
      "mol1"_a, "mol2"_a, "confId1"_a = -1, "confId2"_a = -1,
      "gridSpacing"_a = 0.5,
      "bitsPerPoint"_a = DiscreteValueVect::TWOBITVALUE,
      "vdwScale"_a = 0.8, "stepSize"_a = 0.25, "maxLayers"_a = -1,
      "ignoreHs"_a = true, "allowReordering"_a = true,
      R"DOC(Compute the shape protrude distance between two molecule based on a predefined alignment

ARGUMENTS:
   - mol1 : The first molecule of interest
   - mol2 : The second molecule of interest
   - confId1 : Conformer in the first molecule (defaults to first conformer)
   - confId2 : Conformer in the second molecule (defaults to first conformer)
   - gridSpacing : resolution of the grid used to encode the molecular shapes
   - bitsPerPoint : number of bit used to encode the occupancy at each grid point
                         defaults to two bits per grid point
   - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
               used in the encoding - grid points inside this sphere carry the maximum occupancy
   - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased
                from layer to layer from the maximum value
   - maxLayers : the maximum number of layers - defaults to the number of bits
                 used per grid point - e.g. two bits per grid point will allow 3 layers
   - ignoreHs : when set, the contribution of Hs to the shape will be ignored
   - allowReordering : when set, the order will be automatically updated so that the value calculated
                       is the protrusion of the smaller shape from the larger one.
)DOC");

  m.def(
      "ComputeConfDimsAndOffset",
      [](const Conformer &conf, nb::object trans, double padding) {
        RDGeom::Point3D dims, offSet;
        if (!trans.is_none()) {
          auto transMat = nb::cast<
              nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig>>(
              trans);
          RDGeom::Transform3D ctrans;
          copyTransform(transMat, ctrans);
          MolShapes::computeConfDimsAndOffset(conf, dims, offSet, &ctrans,
                                              padding);
        } else {
          MolShapes::computeConfDimsAndOffset(conf, dims, offSet, nullptr,
                                              padding);
        }
        return nb::make_tuple(dims, offSet);
      },
      "conf"_a, "trans"_a = nb::none(), "padding"_a = 2.0,
      R"DOC(Compute the size of the box that can fit the conformations, and offset
of the box from the origin)DOC");

  m.def(
      "ComputeConfBox",
      [](const Conformer &conf, nb::object trans, double padding) {
        RDGeom::Point3D lowerCorner, upperCorner;
        if (!trans.is_none()) {
          auto transMat = nb::cast<
              nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig>>(
              trans);
          RDGeom::Transform3D ctrans;
          copyTransform(transMat, ctrans);
          MolShapes::computeConfBox(conf, lowerCorner, upperCorner, &ctrans,
                                    padding);
        } else {
          MolShapes::computeConfBox(conf, lowerCorner, upperCorner, nullptr,
                                    padding);
        }
        return nb::make_tuple(lowerCorner, upperCorner);
      },
      "conf"_a, "trans"_a = nb::none(), "padding"_a = 2.0,
      "Compute the lower and upper corners of a cuboid that will fit the conformer");

  m.def(
      "ComputeUnionBox",
      [](nb::tuple box1, nb::tuple box2) {
        if (nb::len(box1) != 2 || nb::len(box2) != 2) {
          throw nb::value_error(
              "In correct format for one of the box: expecting a tuple of two "
              "Point3D");
        }
        RDGeom::Point3D lC1 = nb::cast<RDGeom::Point3D>(box1[0]);
        RDGeom::Point3D uC1 = nb::cast<RDGeom::Point3D>(box1[1]);
        RDGeom::Point3D lC2 = nb::cast<RDGeom::Point3D>(box2[0]);
        RDGeom::Point3D uC2 = nb::cast<RDGeom::Point3D>(box2[1]);
        RDGeom::Point3D lowerCorner, upperCorner;
        MolShapes::computeUnionBox(lC1, uC1, lC2, uC2, lowerCorner,
                                   upperCorner);
        return nb::make_tuple(lowerCorner, upperCorner);
      },
      "box1"_a, "box2"_a,
      R"DOC(Compute the union of two boxes, so that all the points in both boxes are
contained in the new box)DOC");
}
