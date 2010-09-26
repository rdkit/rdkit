// $Id$
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdmolalign_array_API
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include "numpy/arrayobject.h"
#include <GraphMol/MolAlign/AlignMolecules.h>

#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/ROMol.h>

namespace python = boost::python;

namespace RDKit {
  MatchVectType *_translateAtomMap(python::object atomMap) {
    PySequenceHolder<python::object> aMapSeq(atomMap);
    MatchVectType *aMap;
    aMap = 0;
    unsigned int i, nAtms = aMapSeq.size();
    if (nAtms > 0) {
      aMap = new MatchVectType;
      for (i = 0; i < nAtms; ++i) {
        PySequenceHolder<int> item(aMapSeq[i]);
        if (item.size() != 2) {
          delete aMap;
          aMap = 0;
          throw_value_error("Incorrect format for atomMap");
        }
        aMap->push_back(std::pair<int, int>(item[0], item[1]));
      }
    }
    return aMap;
  }

  RDNumeric::DoubleVector *_translateWeights(python::object weights) {
    PySequenceHolder<double> wts(weights);
    unsigned int nwts = wts.size();
    RDNumeric::DoubleVector *wtsVec;
    wtsVec = 0;
    unsigned int i;
    if (nwts > 0) {
      wtsVec = new RDNumeric::DoubleVector(nwts);
      for ( i = 0; i < nwts; i++) {
        wtsVec->setVal(i, wts[i]);
      }
    }
    return wtsVec;
  }

  std::vector<unsigned int>* _translateIds(python::object ids) {
    PySequenceHolder<unsigned int> idsSeq(ids);
    std::vector<unsigned int>* ivec = 0;
    if (idsSeq.size() > 0) {
      ivec = new std::vector<unsigned int>;
      for(unsigned int i = 0; i < idsSeq.size(); ++i) {
        ivec->push_back(idsSeq[i]);
      }
    }
    return ivec;
  }

  void alignMolConfs(ROMol &mol, python::object atomIds=python::list(),
                     python::object confIds=python::list(), 
                     python::object weights=python::list(),
                     bool reflect=false, unsigned int maxIters=50) {
    
    RDNumeric::DoubleVector *wtsVec = _translateWeights(weights);
    std::vector<unsigned int> *aIds = _translateIds(atomIds);
    std::vector<unsigned int> *cIds = _translateIds(confIds);
    MolAlign::alignMolConformers(mol, aIds, cIds, wtsVec, reflect, maxIters);
    if (wtsVec) {
      delete wtsVec;
    }
    if (aIds) {
      delete aIds;
    }
    if (cIds) {
      delete cIds;
    }
  }
    
  PyObject* getMolAlignTransform(const ROMol &prbMol, const ROMol &refMol,
                                      int prbCid=-1, int refCid=-1, 
                                      python::object atomMap=python::list(),
                                      python::object weights=python::list(),
                                      bool reflect=false, unsigned int maxIters=50) {
    MatchVectType *aMap = _translateAtomMap(atomMap);
    unsigned int i, nAtms;
    if (aMap) {
      nAtms = aMap->size();
    } else {
      nAtms = prbMol.getNumAtoms();
    }
    RDNumeric::DoubleVector *wtsVec = _translateWeights(weights);
    if (wtsVec) {
      if (wtsVec->size() != nAtms) {
        throw_value_error("Incorrect number of weights specified");
      }
    }
    RDGeom::Transform3D trans;
    double rmsd = MolAlign::getAlignmentTransform(prbMol, refMol, trans, prbCid, refCid, aMap, 
                                                  wtsVec, reflect, maxIters);
    npy_intp dims[2];
    dims[0] = 4;
    dims[1] = 4;
    PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    double *resData=reinterpret_cast<double *>(res->data);
    unsigned int j, itab;
    const double *tdata = trans.getData();
    for(i=0; i < trans.numRows(); ++i){
      itab = i*4;
      for (j = 0; j < trans.numRows(); ++j) {
        resData[itab + j] = tdata[itab+j];
      }
    }
    if (aMap) {
      delete aMap;
    } 
    if (wtsVec) {
      delete wtsVec;
    }

    PyObject *resTup = PyTuple_New(2);
    PyObject *rmsdItem = PyFloat_FromDouble(rmsd);
    PyTuple_SetItem(resTup,0,rmsdItem);
    PyTuple_SetItem(resTup,1,PyArray_Return(res));
    return resTup;
  }

  double AlignMolecule(ROMol &prbMol, const ROMol &refMol,
                       int prbCid=-1, int refCid=-1, 
                       python::object atomMap=python::list(),
                       python::object weights=python::list(),
                       bool reflect=false, unsigned int maxIters=50) {
    MatchVectType *aMap = _translateAtomMap(atomMap);
    unsigned int nAtms;
    if (aMap) {
      nAtms = aMap->size();
    } else {
      nAtms = prbMol.getNumAtoms();
    }
    RDNumeric::DoubleVector *wtsVec = _translateWeights(weights);
    if (wtsVec) {
      if (wtsVec->size() != nAtms) {
        throw_value_error("Incorrect number of weights specified");
      }
    }
    double rmsd = MolAlign::alignMol(prbMol, refMol, prbCid, refCid, aMap, 
                                     wtsVec, reflect, maxIters);
    if (aMap) {
      delete aMap;
    } 
    if (wtsVec) {
      delete wtsVec;
    }
    return rmsd;
  }
}

BOOST_PYTHON_MODULE(rdMolAlign) {
  import_array();
  python::scope().attr("__doc__") =
    "Module containing functions to align a molecule to a second molecule";
    
  std::string docString = "Compute the transformation required to align a molecule\n\
     \n\
      The 3D transformation required to align the specied conformation in the probe molecule\n\
      to a specified conformation in the reference molecule is computed so that the root mean\n\
      squared distance between a specified set of atoms is minimized\n\
     \n\
     ARGUMENTS\n\
      - prbMol    molecule that is to be aligned\n\
      - refMol    molecule used as the reference for the alignment\n\
      - prbCid    ID of the conformation in the probe to be used \n\
                       for the alignment (defaults to first conformation)\n\
      - refCid    ID of the conformation in the ref molecule to which \n\
                       the alignment is computed (defaults to first conformation)\n\
      - atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)\n\
                       used to compute the alignments. If this mapping is \n\
                       not specified an attempt is made to generate on by\n\
                       substructure matching\n\
      - weights   Optionally specify weights for each of the atom pairs\n\
      - reflect   if true reflect the conformation of the probe molecule\n\
      - maxIters  maximum number of iteration used in mimizing the RMSD\n\
       \n\
      RETURNS\n\
      a tuple of (RMSD value, transform matrix) \n\
    \n";
  python::def("GetAlignmentTransform", RDKit::getMolAlignTransform,
              (python::arg("prbMol"), python::arg("refMol"),
               python::arg("prbCid")=-1,python::arg("refCid")=-1,
               python::arg("atomMap")=python::list(), python::arg("weights")=python::list(),
               python::arg("reflect")=false, python::arg("maxIters")=50),
              docString.c_str());

  docString = "Optimally (minimum RMSD) align a molecule to another molecule\n\
     \n\
      The 3D transformation required to align the specied conformation in the probe molecule\n\
      to a specified conformation in the reference molecule is computed so that the root mean\n\
      squared distance between a specified set of atoms is minimized. \n\
      This transform is then applied to the specified conformation in the probe molecule\n\
     \n\
     ARGUMENTS\n\
      - prbMol    molecule that is to be aligned\n\
      - refMol    molecule used as the reference for the alignment\n\
      - prbCid    ID of the conformation in the probe to be used \n\
                       for the alignment (defaults to first conformation)\n\
      - refCid    ID of the conformation in the ref molecule to which \n\
                       the alignment is computed (defaults to first conformation)\n\
      - atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)\n\
                       used to compute the alignments. If this mapping is \n\
                       not specified an attempt is made to generate on by\n\
                       substructure matching\n\
      - weights   Optionally specify weights for each of the atom pairs\n\
      - reflect   if true reflect the conformation of the probe molecule\n\
      - maxIters  maximum number of iteration used in mimizing the RMSD\n\
       \n\
      RETURNS\n\
      RMSD value\n\
    \n";
  python::def("AlignMol", RDKit::AlignMolecule,
              (python::arg("prbMol"), python::arg("refMol"),
               python::arg("prbCid")=-1,python::arg("refCid")=-1,
               python::arg("atomMap")=python::list(), python::arg("weights")=python::list(),
               python::arg("reflect")=false, python::arg("maxIters")=50),
              docString.c_str());

  docString = "Alignment conformations in a molecule to each other\n\
     \n\
      The first conformation in the molecule is used as the reference\n\
     \n\
     ARGUMENTS\n\
      - mol       molecule of interest\n\
      - atomIds   List of atom ids to use a points for alingment - defaults to all atoms\n\
      - confIds   Ids of conformations to align - defaults to all conformers \n\
      - weights   Optionally specify weights for each of the atom pairs\n\
      - reflect   if true reflect the conformation of the probe molecule\n\
      - maxIters  maximum number of iteration used in mimizing the RMSD\n\
       \n\
      RETURNS\n\
      RMSD value\n\
    \n";
  python::def("AlignMolConformers", RDKit::alignMolConfs,
              (python::arg("mol"), python::arg("atomIds")=python::list(), 
               python::arg("confIds")=python::list(),
               python::arg("weights")=python::list(),
               python::arg("reflect")=false, python::arg("maxIters")=50),
              docString.c_str());
}
    
  
