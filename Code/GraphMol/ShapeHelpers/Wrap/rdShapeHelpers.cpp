// $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdshapehelpers_array_API
#include <boost/python.hpp>
#include "numpy/oldnumeric.h"
#include <GraphMol/ROMol.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/import_array.h>

#include <Geometry/Transform3D.h>
#include <Geometry/UniformGrid3D.h>

#include <GraphMol/ShapeHelpers/ShapeEncoder.h>
#include <GraphMol/ShapeHelpers/ShapeUtils.h>
#include <DataStructs/DiscreteValueVect.h>
#include <Geometry/point.h>

namespace python = boost::python;

namespace RDKit {
  void _copyTransform(const PyArrayObject *transMat, RDGeom::Transform3D &trans) {
    unsigned int nrows = transMat->dimensions[0];
    unsigned int ncols = transMat->dimensions[1];
    if ((nrows != 4) || (ncols != 4)) {
      throw_value_error("The transform has to be square matrix, of size 4x4");
    }
    if (transMat->descr->type_num != PyArray_DOUBLE)
      throw_value_error("Only double arrays allowed for transform object ");
    
    unsigned int dSize = nrows*nrows;
    const double *inData = reinterpret_cast<const double *>(transMat->data);
    double *tData  = trans.getData();
    memcpy(static_cast<void *>(tData), static_cast<const void *>(inData), dSize*sizeof(double));
  }

  python::tuple getConformerDimsAndOffset(const Conformer &conf, python::object trans=python::object(),
                                          double padding=2.5) {
    RDGeom::Point3D dims, offSet;
    PyObject *transObj = trans.ptr();
    if (PyArray_Check(transObj)) {
      PyArrayObject *transMat = reinterpret_cast<PyArrayObject *>(transObj);
      RDGeom::Transform3D ctrans;
      _copyTransform(transMat, ctrans);
      MolShapes::computeConfDimsAndOffset(conf, dims, offSet, &ctrans, padding);
    } else {
      MolShapes::computeConfDimsAndOffset(conf, dims, offSet, 0, padding);
    }
 
    python::tuple res = python::make_tuple(dims, offSet);
    return res;
  }
    
  python::tuple getConfBox(const Conformer &conf, python::object trans=python::object(),
                           double padding=2.5) {
    RDGeom::Point3D lowerCorner, upperCorner;
    PyObject *transObj = trans.ptr();
    if (PyArray_Check(transObj)) {
      PyArrayObject *transMat = reinterpret_cast<PyArrayObject *>(transObj);
      RDGeom::Transform3D ctrans;
      _copyTransform(transMat, ctrans);
      MolShapes::computeConfBox(conf, lowerCorner, upperCorner, &ctrans, padding);
    } else {
      MolShapes::computeConfBox(conf, lowerCorner, upperCorner, 0, padding);
    }
    python::tuple res = python::make_tuple(lowerCorner, upperCorner);
    return res;
  }

  python::tuple getUnionOfTwoBox(python::tuple box1, python::tuple box2) {
    unsigned int len1 = python::extract<unsigned int>(box1.attr("__len__")());
    unsigned int len2 = python::extract<unsigned int>(box2.attr("__len__")());
    if ((len1 != 2) || (len2 != 2)) {
      throw_value_error("In correct format for one of the box: expecting a tuple of two Point3D");
    }
    RDGeom::Point3D lC1 = python::extract<RDGeom::Point3D>(box1.attr("__getitem__")(0));
    RDGeom::Point3D uC1 = python::extract<RDGeom::Point3D>(box1.attr("__getitem__")(1));

    RDGeom::Point3D lC2 = python::extract<RDGeom::Point3D>(box2.attr("__getitem__")(0));
    RDGeom::Point3D uC2 = python::extract<RDGeom::Point3D>(box2.attr("__getitem__")(1));

    RDGeom::Point3D lowerCorner, upperCorner;
    MolShapes::computeUnionBox(lC1, uC1, lC2, uC2, lowerCorner, upperCorner);
    python::tuple res = python::make_tuple(lowerCorner, upperCorner);
    return res;
  }

  void EncodeMolShape(const ROMol &mol, RDGeom::UniformGrid3D &grid, 
                      int confId=-1, python::object trans=python::object(), //PyObject *trans=0,
                      double vdwScale=0.8, double stepSize=0.25, 
                      int maxLayers=-1, bool ignoreHs=true) {
    PyObject *transObj = trans.ptr();
    
    if (PyArray_Check(transObj)) {
      PyArrayObject *transMat = reinterpret_cast<PyArrayObject *>(transObj);
      RDGeom::Transform3D ctrans;
      _copyTransform(transMat, ctrans);
      MolShapes::EncodeShape(mol, grid, confId, &ctrans, vdwScale, stepSize, maxLayers, ignoreHs);
    } else {
      MolShapes::EncodeShape(mol, grid, confId, 0, vdwScale, stepSize, maxLayers, ignoreHs);
    }
  }
  
  double tanimotoMolShapes(const ROMol &mol1, const ROMol &mol2, int confId1=-1, int confId2=-1,
                           double gridSpacing=0.5, 
                           DiscreteValueVect::DiscreteValueType bitsPerPoint=DiscreteValueVect::TWOBITVALUE, 
                           double vdwScale=0.8, double stepSize=0.25, int maxLayers=-1, bool ignoreHs=true) {
    return MolShapes::tanimotoDistance(mol1, mol2, confId1, confId2, gridSpacing, bitsPerPoint,
                                       vdwScale, stepSize, maxLayers, ignoreHs);
  }
  double protrudeMolShapes(const ROMol &mol1, const ROMol &mol2, int confId1=-1, int confId2=-1,
                           double gridSpacing=0.5, 
                           DiscreteValueVect::DiscreteValueType bitsPerPoint=DiscreteValueVect::TWOBITVALUE, 
                           double vdwScale=0.8, double stepSize=0.25, int maxLayers=-1, bool ignoreHs=true,
                           bool allowReordering=true) {
    return MolShapes::protrudeDistance(mol1, mol2, confId1, confId2, gridSpacing, bitsPerPoint,
                                       vdwScale, stepSize, maxLayers, ignoreHs, allowReordering);
  }
}

BOOST_PYTHON_MODULE(rdShapeHelpers) {
  python::scope().attr("__doc__") =
    "Module containing functions to encode and compare the shapes of molecules"
    ;
  
  rdkit_import_array();
  
  //RegisterListConverter<RDKit::Atom*>();
  
  std::string docString = "Encode the shape of a molecule (one of its conformer) onto a grid\n\n\
 \n\
 ARGUMENTS:\n\n\
    - mol : the molecule of interest\n\
    - grid : grid onto which the encoding is written \n\
    - confId : id of the conformation of interest on mol (defaults to the first one) \n\
    - trans : any transformation that needs to be used to encode onto the grid (note the molecule remains unchanged) \n\
    - vdwScale : Scaling factor for the radius of the atoms to determine the base radius \n\
                 used in the encoding - grid points inside this sphere carry the maximum occupancy \n\
    - setpSize : thickness of the layers outside the base radius, the occupancy value is decreased \n\
                 from layer to layer from the maximum value \n\
    - maxLayers : the maximum number of layers - defaults to the number of bits \n\
                  used per grid point - e.g. two bits per grid point will allow 3 layers\n\
    - ignoreHs : when set, the contribution of Hs to the shape will be ignored\n";
  python::def("EncodeShape", RDKit::EncodeMolShape,
              (python::arg("mol"), python::arg("grid"),
               python::arg("confId")=-1, python::arg("trans")=python::object(),
               python::arg("vdwScale")=0.8, python::arg("stepSize")=0.25,
               python::arg("maxLayers")=-1, python::arg("ignoreHs")=true),
              docString.c_str());

  docString = "Compute the shape tanimoto distance between two molecule based on a predefined alignment\n\
  \n\
  ARGUMENTS:\n\
    - mol1 : The first molecule of interest \n\
    - mol2 : The second molecule of interest \n\
    - confId1 : Conformer in the first molecule (defaults to first conformer) \n\
    - confId2 : Conformer in the second molecule (defaults to first conformer) \n\
    - gridSpacing : resolution of the grid used to encode the molecular shapes \n\
    - bitsPerPoint : number of bits used to encode the occupancy at each grid point \n\
                          defaults to two bits per grid point \n\
    - vdwScale : Scaling factor for the radius of the atoms to determine the base radius \n\
                used in the encoding - grid points inside this sphere carry the maximum occupancy \n\
    - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased \n\
                 from layer to layer from the maximum value \n\
    - maxLayers : the maximum number of layers - defaults to the number of bits \n\
                  used per grid point - e.g. two bits per grid point will allow 3 layers \n\
    - ignoreHs : when set, the contribution of Hs to the shape will be ignored\n";
                  
  python::def("ShapeTanimotoDist", RDKit::tanimotoMolShapes,
              (python::arg("mol1"), python::arg("mol2"), 
               python::arg("confId1")=-1, python::arg("confId2")=-1,
               python::arg("gridSpacing")=0.5, 
               python::arg("bitsPerPoint")=RDKit::DiscreteValueVect::TWOBITVALUE,
               python::arg("vdwScale")=0.8, python::arg("stepSize")=0.25,
               python::arg("maxLayers")=-1, python::arg("ignoreHs")=true),
              docString.c_str());

  docString = "Compute the shape protrude distance between two molecule based on a predefined alignment\n\
  \n\
  ARGUMENTS:\n\
    - mol1 : The first molecule of interest \n\
    - mol2 : The second molecule of interest \n\
    - confId1 : Conformer in the first molecule (defaults to first conformer) \n\
    - confId2 : Conformer in the second molecule (defaults to first conformer) \n\
    - gridSpacing : resolution of the grid used to encode the molecular shapes \n\
    - bitsPerPoint : number of bit used to encode the occupancy at each grid point \n\
                          defaults to two bits per grid point \n\
    - vdwScale : Scaling factor for the radius of the atoms to determine the base radius \n\
                used in the encoding - grid points inside this sphere carry the maximum occupancy \n\
    - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased \n\
                 from layer to layer from the maximum value \n\
    - maxLayers : the maximum number of layers - defaults to the number of bits \n\
                  used per grid point - e.g. two bits per grid point will allow 3 layers \n\
    - ignoreHs : when set, the contribution of Hs to the shape will be ignored\n\
    - allowReordering : when set, the order will be automatically updated so that the value calculated\n\
                        is the protrusion of the smaller shape from the larger one.\n";
  python::def("ShapeProtrudeDist", RDKit::protrudeMolShapes,
              (python::arg("mol1"), python::arg("mol2"), 
               python::arg("confId1")=-1, python::arg("confId2")=-1,
               python::arg("gridSpacing")=0.5, 
               python::arg("bitsPerPoint")=RDKit::DiscreteValueVect::TWOBITVALUE,
               python::arg("vdwScale")=0.8, python::arg("stepSize")=0.25,
               python::arg("maxLayers")=-1, python::arg("ignoreHs")=true,
               python::arg("allowReordering")=true),
              docString.c_str());

  
  docString = "Compute the size of the box that can fit the conformations, and offset \n\
   of the box from the origin\n";
  python::def("ComputeConfDimsAndOffset", RDKit::getConformerDimsAndOffset,
              (python::arg("conf"), python::arg("trans")=python::object(),
               python::arg("padding")=2.0), docString.c_str()); 

  docString = "Compute the lower and upper corners of a cuboid that will fit the conformer";
  python::def("ComputeConfBox", RDKit::getConfBox, 
              (python::arg("conf"), python::arg("trans")=python::object(),
               python::arg("padding")=2.0), docString.c_str());
  
  docString = "Compute the union of two boxes, so that all the points in both boxes are \n\
    contained in the new box";
  python::def("ComputeUnionBox", RDKit::getUnionOfTwoBox, docString.c_str());
              
}
              
