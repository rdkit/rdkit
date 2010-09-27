// $Id$
//
//  Copyright (C) 2005-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdmoltransforms_array_API
#include <boost/python.hpp>
#include "numpy/arrayobject.h"
#include <GraphMol/ROMol.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>

namespace python = boost::python;

namespace RDKit {
  PyObject *computeCanonTrans(const Conformer &conf, const RDGeom::Point3D *center=0,
                              bool normalizeCovar=false, bool ignoreHs=true) {
    RDGeom::Transform3D *trans;
    trans = MolTransforms::computeCanonicalTransform(conf, center, 
                                                     normalizeCovar, ignoreHs);
    npy_intp dims[2];
    dims[0] = 4;
    dims[1] = 4;
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    double *resData=reinterpret_cast<double *>(res->data);
    const double *tdata = trans->getData();
    memcpy(static_cast<void *>(resData), static_cast<const void *>(tdata), 4*4*sizeof(double));
    delete trans;
    return PyArray_Return(res);
  }

  void transConformer(Conformer &conf, python::object trans) {
    PyObject *transObj = trans.ptr();
    if (!PyArray_Check(transObj)) {
      throw_value_error("Expecting a numeric array for transformation");
    }
    PyArrayObject *transMat = reinterpret_cast<PyArrayObject *>(transObj);
    unsigned int nrows = transMat->dimensions[0];
    unsigned int dSize = nrows*nrows;
    double *inData = reinterpret_cast<double *>(transMat->data);
    RDGeom::Transform3D transform;
    double *tData  = transform.getData();
    memcpy(static_cast<void *>(tData), static_cast<void *>(inData), dSize*sizeof(double));
    MolTransforms::transformConformer(conf, transform);
  }
}

BOOST_PYTHON_MODULE(rdMolTransforms) {
  
  python::scope().attr("__doc__") =
    "Module containing functions to perform 3D operations like rotate and translate conformations";
   
  import_array();

  std::string docString = "Compute the centroid of the conformation - hydrogens are ignored and no attention\n\
                           if paid to the difference in sizes of the heavy atoms\n";
  python::def("ComputeCentroid", MolTransforms::computeCentroid,
              (python::arg("conf"), python::arg("ignoreHs")=true),
              docString.c_str());

  docString = "Compute the transformation required aligna conformer so that\n\
               the the principal axes align up with the x,y, z axes\n\
               The conformer itself is left unchanged\n\
  ARGUMENTS:\n\
    - conf : the conformer of interest\n\
    - center : optional center point to compute the principal axes around (defaults to the centroid)\n\
    - normalizeCovar : optionally normalize the covariance matrix by the number of atoms\n";
  python::def("ComputeCanonicalTransform", RDKit::computeCanonTrans,
              (python::arg("conf"), python::arg("center")=(RDGeom::Point3D *)(0),
               python::arg("normalizeCovar")=false, python::arg("ignoreHs")=true),
              docString.c_str());

  python::def("TransformConformer", RDKit::transConformer,
              "Transform the coordinates of a conformer");


  docString = "Canonicalize the orientation of a conformer so that its principal axes\n\
               around the specified center point coincide with the x, y, z axes\n\
  \n\
  ARGUMENTS:\n\
    - conf : conformer of interest \n\
    - center : optionally center point about which the principal axes are computed \n\
                          if not specified the centroid of the conformer will be used\n\
    - normalizeCovar : Optionally normalize the covariance matrix by the number of atoms\n";
  python::def("CanonicalizeConformer", MolTransforms::canonicalizeConformer,
              (python::arg("conf"), python::arg("center")=(RDGeom::Point3D *)(0), 
               python::arg("normalizeCovar")=false, python::arg("ignoreHs")=true),
              docString.c_str());
  
  python::def("CanonicalizeMol", MolTransforms::canonicalizeMol,
              (python::arg("mol"), python::arg("normalizeCovar")=false, python::arg("ignoreHs")=true),
              "Loop over the conformers in a molecule and canonicalize their orientation");
}

              
