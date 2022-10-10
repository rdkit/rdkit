// $Id$
//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//  Copyright (C) 2013 Paolo Tosco
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define PY_ARRAY_UNIQUE_SYMBOL rdmolalign_array_API
#include <RDBoost/python.h>
#include <RDBoost/import_array.h>
#include <RDBoost/boost_numpy.h>
#include <utility>
#include "numpy/arrayobject.h"
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolAlign/O3AAlignMolecules.h>
#include <ForceField/Wrap/PyForceField.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <RDBoost/PySequenceHolder.h>
#include <RDBoost/Wrap.h>
#include <GraphMol/ROMol.h>

namespace python = boost::python;

namespace RDKit {
void alignMolConfs(ROMol &mol, python::object atomIds, python::object confIds,
                   python::object weights, bool reflect, unsigned int maxIters,
                   python::object RMSlist) {
  std::unique_ptr<RDNumeric::DoubleVector> wtsVec(translateDoubleSeq(weights));
  std::unique_ptr<std::vector<unsigned int>> aIds(translateIntSeq(atomIds));
  std::unique_ptr<std::vector<unsigned int>> cIds(translateIntSeq(confIds));
  std::unique_ptr<std::vector<double>> RMSvector;
  if (RMSlist != python::object()) {
    RMSvector.reset(new std::vector<double>());
  }
  {
    NOGIL gil;
    MolAlign::alignMolConformers(mol, aIds.get(), cIds.get(), wtsVec.get(),
                                 reflect, maxIters, RMSvector.get());
  }
  if (RMSvector) {
    auto &pyl = static_cast<python::list &>(RMSlist);
    for (double i : *RMSvector) {
      pyl.append(i);
    }
  }
}

PyObject *generateRmsdTransMatchPyTuple(double rmsd,
                                        const RDGeom::Transform3D &trans,
                                        const MatchVectType *match = nullptr) {
  npy_intp dims[2];
  dims[0] = 4;
  dims[1] = 4;
  auto *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  auto *resData = reinterpret_cast<double *>(PyArray_DATA(res));
  unsigned int i, j, itab;
  const double *tdata = trans.getData();
  for (i = 0; i < trans.numRows(); ++i) {
    itab = i * 4;
    for (j = 0; j < trans.numRows(); ++j) {
      resData[itab + j] = tdata[itab + j];
    }
  }
  PyObject *resTup = PyTuple_New(2 + (match ? 1 : 0));
  PyObject *rmsdItem = PyFloat_FromDouble(rmsd);
  PyTuple_SetItem(resTup, 0, rmsdItem);
  PyTuple_SetItem(resTup, 1, PyArray_Return(res));
  if (match) {
    python::list pairList;
    for (const auto &pair : *match) {
      pairList.append(python::make_tuple(pair.first, pair.second));
    }
    auto pairTup = new python::tuple(pairList);
    PyTuple_SetItem(resTup, 2, pairTup->ptr());
  }
  return resTup;
}

PyObject *getMolAlignTransform(const ROMol &prbMol, const ROMol &refMol,
                               int prbCid = -1, int refCid = -1,
                               python::object atomMap = python::list(),
                               python::object weights = python::list(),
                               bool reflect = false,
                               unsigned int maxIters = 50) {
  std::unique_ptr<MatchVectType> aMap(translateAtomMap(atomMap));
  unsigned int nAtms;
  if (aMap) {
    nAtms = aMap->size();
  } else {
    nAtms = prbMol.getNumAtoms();
  }
  std::unique_ptr<RDNumeric::DoubleVector> wtsVec(translateDoubleSeq(weights));
  if (wtsVec) {
    if (wtsVec->size() != nAtms) {
      throw_value_error("Incorrect number of weights specified");
    }
  }
  RDGeom::Transform3D trans;
  double rmsd;
  {
    NOGIL gil;
    rmsd = MolAlign::getAlignmentTransform(prbMol, refMol, trans, prbCid,
                                           refCid, aMap.get(), wtsVec.get(),
                                           reflect, maxIters);
  }

  return generateRmsdTransMatchPyTuple(rmsd, trans);
}

PyObject *getBestMolAlignTransform(const ROMol &prbMol, const ROMol &refMol,
                                   int prbCid = -1, int refCid = -1,
                                   python::object map = python::list(),
                                   int maxMatches = 1000000,
                                   bool symmetrizeTerminalGroups = true,
                                   python::object weights = python::list(),
                                   bool reflect = false,
                                   unsigned int maxIters = 50) {
  std::vector<MatchVectType> aMapVec;
  unsigned int nAtms = 0;
  if (map != python::object()) {
    aMapVec = translateAtomMapSeq(map);
    if (!aMapVec.empty()) {
      nAtms = aMapVec.front().size();
    }
  }
  std::unique_ptr<RDNumeric::DoubleVector> wtsVec(translateDoubleSeq(weights));
  if (wtsVec) {
    if (wtsVec->size() != nAtms) {
      throw_value_error("Incorrect number of weights specified");
    }
  }
  RDGeom::Transform3D bestTrans;
  MatchVectType bestMatch;
  double rmsd;
  {
    NOGIL gil;
    rmsd = MolAlign::getBestAlignmentTransform(
        prbMol, refMol, bestTrans, bestMatch, prbCid, refCid, aMapVec,
        maxMatches, symmetrizeTerminalGroups, wtsVec.get(), reflect, maxIters);
  }

  return generateRmsdTransMatchPyTuple(rmsd, bestTrans, &bestMatch);
}

double AlignMolecule(ROMol &prbMol, const ROMol &refMol, int prbCid = -1,
                     int refCid = -1, python::object atomMap = python::list(),
                     python::object weights = python::list(),
                     bool reflect = false, unsigned int maxIters = 50) {
  std::unique_ptr<MatchVectType> aMap(translateAtomMap(atomMap));
  unsigned int nAtms;
  if (aMap) {
    nAtms = aMap->size();
  } else {
    nAtms = prbMol.getNumAtoms();
  }
  std::unique_ptr<RDNumeric::DoubleVector> wtsVec(translateDoubleSeq(weights));
  if (wtsVec) {
    if (wtsVec->size() != nAtms) {
      throw_value_error("Incorrect number of weights specified");
    }
  }

  double rmsd;
  {
    NOGIL gil;
    rmsd = MolAlign::alignMol(prbMol, refMol, prbCid, refCid, aMap.get(),
                              wtsVec.get(), reflect, maxIters);
  }
  return rmsd;
}

double GetBestRMS(ROMol &prbMol, ROMol &refMol, int prbId, int refId,
                  python::object map, int maxMatches,
                  bool symmetrizeTerminalGroups,
                  python::object weights = python::list()) {
  std::vector<MatchVectType> aMapVec;
  if (map != python::object()) {
    aMapVec = translateAtomMapSeq(map);
  }
  std::unique_ptr<RDNumeric::DoubleVector> wtsVec(translateDoubleSeq(weights));
  double rmsd;
  {
    NOGIL gil;
    rmsd =
        MolAlign::getBestRMS(prbMol, refMol, prbId, refId, aMapVec, maxMatches,
                             symmetrizeTerminalGroups, wtsVec.get());
  }
  return rmsd;
}

double CalcRMS(ROMol &prbMol, ROMol &refMol, int prbCid, int refCid,
               python::object map, int maxMatches,
               bool symmetrizeTerminalGroups,
               python::object weights = python::list()) {
  std::vector<MatchVectType> aMapVec;
  if (map != python::object()) {
    aMapVec = translateAtomMapSeq(map);
  }
  RDNumeric::DoubleVector *wtsVec = translateDoubleSeq(weights);
  double rmsd;
  {
    NOGIL gil;
    rmsd = MolAlign::CalcRMS(prbMol, refMol, prbCid, refCid, aMapVec,
                             maxMatches, symmetrizeTerminalGroups, wtsVec);
  }
  return rmsd;
}

namespace MolAlign {
class PyO3A {
 public:
  PyO3A(O3A *o) : o3a(o){};
  PyO3A(boost::shared_ptr<O3A> o) : o3a(std::move(o)){};
  ~PyO3A() = default;
  double align() { return o3a.get()->align(); };
  PyObject *trans() {
    RDGeom::Transform3D trans;
    double rmsd = o3a.get()->trans(trans);
    return RDKit::generateRmsdTransMatchPyTuple(rmsd, trans);
  };
  double score() { return o3a.get()->score(); };
  boost::python::list matches() {
    boost::python::list matchList;
    const RDKit::MatchVectType *o3aMatchVect = o3a->matches();

    for (const auto &i : *o3aMatchVect) {
      boost::python::list match;
      match.append(i.first);
      match.append(i.second);
      matchList.append(match);
    }

    return matchList;
  };
  boost::python::list weights() {
    boost::python::list weightList;
    const RDNumeric::DoubleVector *o3aWeights = o3a->weights();

    for (unsigned int i = 0; i < o3aWeights->size(); ++i) {
      weightList.append((*o3aWeights)[i]);
    }

    return weightList;
  };
  boost::shared_ptr<O3A> o3a;
};
PyO3A *getMMFFO3A(ROMol &prbMol, ROMol &refMol, python::object prbProps,
                  python::object refProps, int prbCid = -1, int refCid = -1,
                  bool reflect = false, unsigned int maxIters = 50,
                  unsigned int options = 0,
                  python::list constraintMap = python::list(),
                  python::list constraintWeights = python::list()) {
  std::unique_ptr<MatchVectType> cMap;
  if (python::len(constraintMap)) {
    cMap.reset(translateAtomMap(constraintMap));
  }
  std::unique_ptr<RDNumeric::DoubleVector> cWts;
  if (cMap) {
    cWts.reset(translateDoubleSeq(constraintWeights));
    if (cWts) {
      if (cMap->size() != cWts->size()) {
        throw_value_error(
            "The number of weights should match the number of constraints");
      }
    }
    for (const auto &i : *cMap) {
      if ((i.first < 0) || (i.first >= rdcast<int>(prbMol.getNumAtoms())) ||
          (i.second < 0) || (i.second >= rdcast<int>(refMol.getNumAtoms()))) {
        throw_value_error("Constrained atom idx out of range");
      }
      if ((prbMol[i.first]->getAtomicNum() == 1) ||
          (refMol[i.second]->getAtomicNum() == 1)) {
        throw_value_error("Constrained atoms must be heavy atoms");
      }
    }
  }
  std::unique_ptr<MMFF::MMFFMolProperties> prbMolProps;
  MMFF::MMFFMolProperties *prbMolPropsPtr = nullptr;
  std::unique_ptr<MMFF::MMFFMolProperties> refMolProps;
  MMFF::MMFFMolProperties *refMolPropsPtr = nullptr;

  if (prbProps != python::object()) {
    ForceFields::PyMMFFMolProperties *prbPyMMFFMolProperties =
        python::extract<ForceFields::PyMMFFMolProperties *>(prbProps);
    prbMolPropsPtr = prbPyMMFFMolProperties->mmffMolProperties.get();
  } else {
    prbMolProps.reset(new MMFF::MMFFMolProperties(prbMol));
    if (!prbMolProps->isValid()) {
      throw_value_error("missing MMFF94 parameters for probe molecule");
    }
    prbMolPropsPtr = prbMolProps.get();
  }
  if (refProps != python::object()) {
    ForceFields::PyMMFFMolProperties *refPyMMFFMolProperties =
        python::extract<ForceFields::PyMMFFMolProperties *>(refProps);
    refMolPropsPtr = refPyMMFFMolProperties->mmffMolProperties.get();
  } else {
    refMolProps.reset(new MMFF::MMFFMolProperties(refMol));
    if (!refMolProps->isValid()) {
      throw_value_error("missing MMFF94 parameters for reference molecule");
    }
    refMolPropsPtr = refMolProps.get();
  }
  O3A *o3a;
  {
    NOGIL gil;
    o3a = new MolAlign::O3A(prbMol, refMol, prbMolPropsPtr, refMolPropsPtr,
                            MolAlign::O3A::MMFF94, prbCid, refCid, reflect,
                            maxIters, options, cMap.get(), cWts.get());
  }

  return new PyO3A(o3a);
}

python::tuple getMMFFO3AForConfs(
    ROMol &prbMol, ROMol &refMol, int numThreads, python::object prbProps,
    python::object refProps, int refCid = -1, bool reflect = false,
    unsigned int maxIters = 50, unsigned int options = 0,
    python::list constraintMap = python::list(),
    python::list constraintWeights = python::list()) {
  std::unique_ptr<MatchVectType> cMap;
  if (python::len(constraintMap)) {
    cMap.reset(translateAtomMap(constraintMap));
  }
  std::unique_ptr<RDNumeric::DoubleVector> cWts;
  if (cMap) {
    cWts.reset(translateDoubleSeq(constraintWeights));
    if (cWts) {
      if (cMap->size() != cWts->size()) {
        throw_value_error(
            "The number of weights should match the number of constraints");
      }
    }
    for (const auto &i : *cMap) {
      if ((i.first < 0) || (i.first >= rdcast<int>(prbMol.getNumAtoms())) ||
          (i.second < 0) || (i.second >= rdcast<int>(refMol.getNumAtoms()))) {
        throw_value_error("Constrained atom idx out of range");
      }
      if ((prbMol[i.first]->getAtomicNum() == 1) ||
          (refMol[i.second]->getAtomicNum() == 1)) {
        throw_value_error("Constrained atoms must be heavy atoms");
      }
    }
  }
  std::unique_ptr<MMFF::MMFFMolProperties> prbMolProps;
  MMFF::MMFFMolProperties *prbMolPropsPtr = nullptr;
  std::unique_ptr<MMFF::MMFFMolProperties> refMolProps;
  MMFF::MMFFMolProperties *refMolPropsPtr = nullptr;

  if (prbProps != python::object()) {
    ForceFields::PyMMFFMolProperties *prbPyMMFFMolProperties =
        python::extract<ForceFields::PyMMFFMolProperties *>(prbProps);
    prbMolPropsPtr = prbPyMMFFMolProperties->mmffMolProperties.get();
  } else {
    prbMolProps.reset(new MMFF::MMFFMolProperties(prbMol));
    if (!prbMolProps->isValid()) {
      throw_value_error("missing MMFF94 parameters for probe molecule");
    }
    prbMolPropsPtr = prbMolProps.get();
  }
  if (refProps != python::object()) {
    ForceFields::PyMMFFMolProperties *refPyMMFFMolProperties =
        python::extract<ForceFields::PyMMFFMolProperties *>(refProps);
    refMolPropsPtr = refPyMMFFMolProperties->mmffMolProperties.get();
  } else {
    refMolProps.reset(new MMFF::MMFFMolProperties(refMol));
    if (!refMolProps->isValid()) {
      throw_value_error("missing MMFF94 parameters for reference molecule");
    }
    refMolPropsPtr = refMolProps.get();
  }
  std::vector<boost::shared_ptr<O3A>> res;
  {
    NOGIL gil;
    getO3AForProbeConfs(prbMol, refMol, prbMolPropsPtr, refMolPropsPtr, res,
                        numThreads, MolAlign::O3A::MMFF94, refCid, reflect,
                        maxIters, options, cMap.get(), cWts.get());
  }

  python::list pyres;
  for (auto &i : res) {
    pyres.append(new PyO3A(i));
  }

  return python::tuple(pyres);
}

PyO3A *getCrippenO3A(ROMol &prbMol, ROMol &refMol,
                     python::list prbCrippenContribs,
                     python::list refCrippenContribs, int prbCid = -1,
                     int refCid = -1, bool reflect = false,
                     unsigned int maxIters = 50, unsigned int options = 0,
                     python::list constraintMap = python::list(),
                     python::list constraintWeights = python::list()) {
  std::unique_ptr<MatchVectType> cMap;
  if (python::len(constraintMap)) {
    cMap.reset(translateAtomMap(constraintMap));
  }
  std::unique_ptr<RDNumeric::DoubleVector> cWts;
  if (cMap) {
    cWts.reset(translateDoubleSeq(constraintWeights));
    if (cWts) {
      if (cMap->size() != cWts->size()) {
        throw_value_error(
            "The number of weights should match the number of constraints");
      }
    }
    for (const auto &i : *cMap) {
      if ((i.first < 0) || (i.first >= rdcast<int>(prbMol.getNumAtoms())) ||
          (i.second < 0) || (i.second >= rdcast<int>(refMol.getNumAtoms()))) {
        throw_value_error("Constrained atom idx out of range");
      }
      if ((prbMol[i.first]->getAtomicNum() == 1) ||
          (refMol[i.second]->getAtomicNum() == 1)) {
        throw_value_error("Constrained atoms must be heavy atoms");
      }
    }
  }
  unsigned int prbNAtoms = prbMol.getNumAtoms();
  std::vector<double> prbLogpContribs(prbNAtoms);
  unsigned int refNAtoms = refMol.getNumAtoms();
  std::vector<double> refLogpContribs(refNAtoms);

  if ((prbCrippenContribs != python::list()) &&
      (python::len(prbCrippenContribs) == prbNAtoms)) {
    for (unsigned int i = 0; i < prbNAtoms; ++i) {
      python::tuple logpMRTuple =
          python::extract<python::tuple>(prbCrippenContribs[i]);
      prbLogpContribs[i] = python::extract<double>(logpMRTuple[0]);
    }
  } else {
    std::vector<double> prbMRContribs(prbNAtoms);
    std::vector<unsigned int> prbAtomTypes(prbNAtoms);
    std::vector<std::string> prbAtomTypeLabels(prbNAtoms);
    Descriptors::getCrippenAtomContribs(prbMol, prbLogpContribs, prbMRContribs,
                                        true, &prbAtomTypes,
                                        &prbAtomTypeLabels);
  }
  if ((refCrippenContribs != python::list()) &&
      (python::len(refCrippenContribs) == refNAtoms)) {
    for (unsigned int i = 0; i < refNAtoms; ++i) {
      python::tuple logpMRTuple =
          python::extract<python::tuple>(refCrippenContribs[i]);
      refLogpContribs[i] = python::extract<double>(logpMRTuple[0]);
    }
  } else {
    std::vector<double> refMRContribs(refNAtoms);
    std::vector<unsigned int> refAtomTypes(refNAtoms);
    std::vector<std::string> refAtomTypeLabels(refNAtoms);
    Descriptors::getCrippenAtomContribs(refMol, refLogpContribs, refMRContribs,
                                        true, &refAtomTypes,
                                        &refAtomTypeLabels);
  }
  O3A *o3a;
  {
    NOGIL gil;
    o3a = new MolAlign::O3A(prbMol, refMol, &prbLogpContribs, &refLogpContribs,
                            MolAlign::O3A::CRIPPEN, prbCid, refCid, reflect,
                            maxIters, options, cMap.get(), cWts.get());
  }

  return new PyO3A(o3a);
}

python::tuple getCrippenO3AForConfs(
    ROMol &prbMol, ROMol &refMol, int numThreads,
    python::list prbCrippenContribs, python::list refCrippenContribs,
    int refCid = -1, bool reflect = false, unsigned int maxIters = 50,
    unsigned int options = 0, python::list constraintMap = python::list(),
    python::list constraintWeights = python::list()) {
  std::unique_ptr<MatchVectType> cMap;
  if (python::len(constraintMap)) {
    cMap.reset(translateAtomMap(constraintMap));
  }
  std::unique_ptr<RDNumeric::DoubleVector> cWts;
  if (cMap) {
    cWts.reset(translateDoubleSeq(constraintWeights));
    if (cWts) {
      if (cMap->size() != cWts->size()) {
        throw_value_error(
            "The number of weights should match the number of constraints");
      }
    }
    for (const auto &i : *cMap) {
      if ((i.first < 0) || (i.first >= rdcast<int>(prbMol.getNumAtoms())) ||
          (i.second < 0) || (i.second >= rdcast<int>(refMol.getNumAtoms()))) {
        throw_value_error("Constrained atom idx out of range");
      }
      if ((prbMol[i.first]->getAtomicNum() == 1) ||
          (refMol[i.second]->getAtomicNum() == 1)) {
        throw_value_error("Constrained atoms must be heavy atoms");
      }
    }
  }
  unsigned int prbNAtoms = prbMol.getNumAtoms();
  std::vector<double> prbLogpContribs(prbNAtoms);
  unsigned int refNAtoms = refMol.getNumAtoms();
  std::vector<double> refLogpContribs(refNAtoms);

  if ((prbCrippenContribs != python::list()) &&
      (python::len(prbCrippenContribs) == prbNAtoms)) {
    for (unsigned int i = 0; i < prbNAtoms; ++i) {
      python::tuple logpMRTuple =
          python::extract<python::tuple>(prbCrippenContribs[i]);
      prbLogpContribs[i] = python::extract<double>(logpMRTuple[0]);
    }
  } else {
    std::vector<double> prbMRContribs(prbNAtoms);
    std::vector<unsigned int> prbAtomTypes(prbNAtoms);
    std::vector<std::string> prbAtomTypeLabels(prbNAtoms);
    Descriptors::getCrippenAtomContribs(prbMol, prbLogpContribs, prbMRContribs,
                                        true, &prbAtomTypes,
                                        &prbAtomTypeLabels);
  }
  if ((refCrippenContribs != python::list()) &&
      (python::len(refCrippenContribs) == refNAtoms)) {
    for (unsigned int i = 0; i < refNAtoms; ++i) {
      python::tuple logpMRTuple =
          python::extract<python::tuple>(refCrippenContribs[i]);
      refLogpContribs[i] = python::extract<double>(logpMRTuple[0]);
    }
  } else {
    std::vector<double> refMRContribs(refNAtoms);
    std::vector<unsigned int> refAtomTypes(refNAtoms);
    std::vector<std::string> refAtomTypeLabels(refNAtoms);
    Descriptors::getCrippenAtomContribs(refMol, refLogpContribs, refMRContribs,
                                        true, &refAtomTypes,
                                        &refAtomTypeLabels);
  }
  std::vector<boost::shared_ptr<O3A>> res;
  {
    NOGIL gil;
    getO3AForProbeConfs(prbMol, refMol, &prbLogpContribs, &refLogpContribs, res,
                        numThreads, MolAlign::O3A::CRIPPEN, refCid, reflect,
                        maxIters, options, cMap.get(), cWts.get());
  }
  python::list pyres;
  for (auto &re : res) {
    pyres.append(new PyO3A(re));
  }

  return python::tuple(pyres);
}
}  // end of namespace MolAlign
}  // end of namespace RDKit

BOOST_PYTHON_MODULE(rdMolAlign) {
  rdkit_import_array();
  python::scope().attr("__doc__") =
      "Module containing functions to align a molecule to a second molecule";

  std::string docString =
      "Compute the transformation required to align a molecule\n\
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
                       not specified an attempt is made to generate one by\n\
                       substructure matching\n\
      - weights   Optionally specify weights for each of the atom pairs\n\
      - reflect   if true reflect the conformation of the probe molecule\n\
      - maxIters  maximum number of iterations used in minimizing the RMSD\n\
       \n\
      RETURNS\n\
      a tuple of (RMSD value, transform matrix) \n\
    \n";
  python::def(
      "GetAlignmentTransform", RDKit::getMolAlignTransform,
      (python::arg("prbMol"), python::arg("refMol"), python::arg("prbCid") = -1,
       python::arg("refCid") = -1, python::arg("atomMap") = python::list(),
       python::arg("weights") = python::list(), python::arg("reflect") = false,
       python::arg("maxIters") = 50),
      docString.c_str());

  docString =
      "Compute the optimal RMS, transformation and atom map for aligning\n\
      two molecules, taking symmetry into account. Molecule coordinates\n\
      are left unaltered.\n\
    \n\
      This function will attempt to align all permutations of matching atom\n\
      orders in both molecules, for some molecules it will lead to 'combinatorial\n\
      explosion' especially if hydrogens are present.\n\
      Use 'GetAlignmentTransform' to align molecules without changing the atom order.\n\
    \n\
     ARGUMENTS\n\
      - prbMol      molecule that is to be aligned\n\
      - refMol      molecule used as the reference for the alignment\n\
      - prbCid      ID of the conformation in the probe to be used \n\
                    for the alignment (defaults to first conformation)\n\
      - refCid      ID of the conformation in the ref molecule to which \n\
                    the alignment is computed (defaults to first conformation)\n\
      - map:        (optional) a list of lists of (probeAtomId, refAtomId)\n\
                    tuples with the atom-atom mappings of the two\n\
                    molecules. If not provided, these will be generated\n\
                    using a substructure search.\n\
      - maxMatches  (optional) if atomMap is empty, this will be the max number of\n\
                    matches found in a SubstructMatch().\n\
      - symmetrizeConjugatedTerminalGroups (optional) if set, conjugated\n\
                    terminal functional groups (like nitro or carboxylate)\n\
                    will be considered symmetrically.\n\
      - weights     Optionally specify weights for each of the atom pairs\n\
      - reflect     if true reflect the conformation of the probe molecule\n\
      - maxIters    maximum number of iterations used in minimizing the RMSD\n\
       \n\
      RETURNS\n\
      a tuple of (RMSD value, best transform matrix, best atom map)\n\
    \n";

  python::def(
      "GetBestAlignmentTransform", RDKit::getBestMolAlignTransform,
      (python::arg("prbMol"), python::arg("refMol"), python::arg("prbCid") = -1,
       python::arg("refCid") = -1, python::arg("map") = python::list(),
       python::arg("maxMatches") = 1000000,
       python::arg("symmetrizeConjugatedTerminalGroups") = true,
       python::arg("weights") = python::list(), python::arg("reflect") = false,
       python::arg("maxIters") = 50),
      docString.c_str());

  docString =
      "Optimally (minimum RMSD) align a molecule to another molecule\n\
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
                       not specified an attempt is made to generate one by\n\
                       substructure matching\n\
      - weights   Optionally specify weights for each of the atom pairs\n\
      - reflect   if true reflect the conformation of the probe molecule\n\
      - maxIters  maximum number of iterations used in minimizing the RMSD\n\
       \n\
      RETURNS\n\
      RMSD value\n\
    \n";
  python::def(
      "AlignMol", RDKit::AlignMolecule,
      (python::arg("prbMol"), python::arg("refMol"), python::arg("prbCid") = -1,
       python::arg("refCid") = -1, python::arg("atomMap") = python::list(),
       python::arg("weights") = python::list(), python::arg("reflect") = false,
       python::arg("maxIters") = 50),
      docString.c_str());

  docString =
      "Returns the optimal RMS for aligning two molecules, taking\n\
       symmetry into account. As a side-effect, the probe molecule is\n\
       left in the aligned state.\n\
      \n\
       Note:\n\
       This function will attempt to align all permutations of matching atom\n\
       orders in both molecules, for some molecules it will lead to\n\
       'combinatorial explosion' especially if hydrogens are present.\n\
       Use 'rdkit.Chem.AllChem.AlignMol' to align molecules without changing\n\
       the atom order.\n\
      \n\
       ARGUMENTS\n\
        - prbMol:      the molecule to be aligned to the reference\n\
        - refMol:      the reference molecule\n\
        - prbId:       (optional) probe conformation to use\n\
        - refId:       (optional) reference conformation to use\n\
        - map:         (optional) a list of lists of (probeAtomId,refAtomId)\n\
                       tuples with the atom-atom mappings of the two\n\
                       molecules. If not provided, these will be generated\n\
                       using a substructure search.\n\
        - maxMatches:  (optional) if map isn't specified, this will be\n\
                       the max number of matches found in a SubstructMatch()\n\
        - symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugated\n\
                       terminal functional groups (like nitro or carboxylate)\n\
                       will be considered symmetrically\n\
        - weights:     (optional) weights for mapping\n\
       \n\
      RETURNS\n\
      The best RMSD found\n\
    \n";
  python::def(
      "GetBestRMS", RDKit::GetBestRMS,
      (python::arg("prbMol"), python::arg("refMol"), python::arg("prbId") = -1,
       python::arg("refId") = -1, python::arg("map") = python::object(),
       python::arg("maxMatches") = 1000000,
       python::arg("symmetrizeConjugatedTerminalGroups") = true,
       python::arg("weights") = python::list()),
      docString.c_str());

  docString =
      "Returns the RMS between two molecules, taking symmetry into account.\n\
       In contrast to getBestRMS, the RMS is computed 'in place', i.e.\n\
       probe molecules are not aligned to the reference ahead of the\n\
       RMS calculation. This is useful, for example, to compute\n\
       the RMSD between docking poses and the co-crystallized ligand.\n\
      \n\
       Note:\n\
       This function will attempt to match all permutations of matching atom\n\
       orders in both molecules, for some molecules it will lead to\n\
       'combinatorial explosion' especially if hydrogens are present.\n\
      \n\
       ARGUMENTS\n\
        - prbMol:      the molecule to be aligned to the reference\n\
        - refMol:      the reference molecule\n\
        - prbCId:      (optional) probe conformation to use\n\
        - refCId:      (optional) reference conformation to use\n\
        - map:         (optional) a list of lists of (probeAtomId, refAtomId)\n\
                       tuples with the atom-atom mappings of the two\n\
                       molecules. If not provided, these will be generated\n\
                       using a substructure search.\n\
        - maxMatches:  (optional) if map isn't specified, this will be\n\
                       the max number of matches found in a SubstructMatch()\n\
        - symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugated\n\
                       terminal functional groups (like nitro or carboxylate)\n\
                       will be considered symmetrically\n\
        - weights:     (optional) weights for mapping\n\
       \n\
      RETURNS\n\
      The best RMSD found\n\
    \n";
  python::def(
      "CalcRMS", RDKit::CalcRMS,
      (python::arg("prbMol"), python::arg("refMol"), python::arg("prbId") = -1,
       python::arg("refId") = -1, python::arg("map") = python::object(),
       python::arg("maxMatches") = 1000000,
       python::arg("symmetrizeConjugatedTerminalGroups") = true,
       python::arg("weights") = python::list()),
      docString.c_str());

  docString =
      "Align conformations in a molecule to each other\n\
     \n\
      The first conformation in the molecule is used as the reference\n\
     \n\
     ARGUMENTS\n\
      - mol          molecule of interest\n\
      - atomIds      List of atom ids to use a points for alingment - defaults to all atoms\n\
      - confIds      Ids of conformations to align - defaults to all conformers \n\
      - weights      Optionally specify weights for each of the atom pairs\n\
      - reflect      if true reflect the conformation of the probe molecule\n\
      - maxIters     maximum number of iterations used in minimizing the RMSD\n\
      - RMSlist      if provided, fills in the RMS values between the reference\n\
		     conformation and the other aligned conformations\n\
       \n\
    \n";
  python::def(
      "AlignMolConformers", RDKit::alignMolConfs,
      (python::arg("mol"), python::arg("atomIds") = python::list(),
       python::arg("confIds") = python::list(),
       python::arg("weights") = python::list(), python::arg("reflect") = false,
       python::arg("maxIters") = 50, python::arg("RMSlist") = python::object()),
      docString.c_str());

  docString =
      "Perform a random transformation on a molecule\n\
     \n\
     ARGUMENTS\n\
      - mol    molecule that is to be transformed\n\
      - cid    ID of the conformation in the mol to be transformed\n\
               (defaults to first conformation)\n\
      - seed   seed used to initialize the random generator\n\
               (defaults to -1, that is no seeding)\n\
       \n\
    \n";
  python::def(
      "RandomTransform", RDKit::MolAlign::randomTransform,
      (python::arg("mol"), python::arg("cid") = -1, python::arg("seed") = -1),
      docString.c_str());

  python::class_<RDKit::MolAlign::PyO3A,
                 boost::shared_ptr<RDKit::MolAlign::PyO3A>>(
      "O3A", "Open3DALIGN object", python::no_init)
      .def("Align", &RDKit::MolAlign::PyO3A::align, (python::arg("self")),
           "aligns probe molecule onto reference molecule")
      .def("Trans", &RDKit::MolAlign::PyO3A::trans, (python::arg("self")),
           "returns the transformation which aligns probe molecule onto "
           "reference molecule")
      .def("Score", &RDKit::MolAlign::PyO3A::score, (python::arg("self")),
           "returns the O3AScore of the alignment")
      .def("Matches", &RDKit::MolAlign::PyO3A::matches, (python::arg("self")),
           "returns the AtomMap as found by Open3DALIGN")
      .def("Weights", &RDKit::MolAlign::PyO3A::weights, (python::arg("self")),
           "returns the weight vector as found by Open3DALIGN");

  docString =
      "Get an O3A object with atomMap and weights vectors to overlay\n\
      the probe molecule onto the reference molecule based on\n\
      MMFF atom types and charges\n\
     \n\
     ARGUMENTS\n\
      - prbMol                   molecule that is to be aligned\n\
      - refMol                   molecule used as the reference for the alignment\n\
      - prbPyMMFFMolProperties   PyMMFFMolProperties object for the probe molecule as returned\n\
                                 by SetupMMFFForceField()\n\
      - refPyMMFFMolProperties   PyMMFFMolProperties object for the reference molecule as returned\n\
                                 by SetupMMFFForceField()\n\
      - prbCid                   ID of the conformation in the probe to be used \n\
                                 for the alignment (defaults to first conformation)\n\
      - refCid                   ID of the conformation in the ref molecule to which \n\
                                 the alignment is computed (defaults to first conformation)\n\
      - reflect                  if true reflect the conformation of the probe molecule\n\
                                 (defaults to false)\n\
      - maxIters                 maximum number of iterations used in minimizing the RMSD\n\
                                 (defaults to 50)\n\
      - options                  least 2 significant bits encode accuracy\n\
                                 (0: maximum, 3: minimum; defaults to 0)\n\
                                 bit 3 triggers local optimization of the alignment\n\
                                 (no computation of the cost matrix; defaults: off)\n\
      - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)\n\
                                 which shall be used for the alignment (defaults to [])\n\
      - constraintWeights        optionally specify weights for each of the constraints\n\
                                 (weights default to 100.0)\n\
       \n\
      RETURNS\n\
      The O3A object\n\
    \n";
  python::def("GetO3A", RDKit::MolAlign::getMMFFO3A,
              (python::arg("prbMol"), python::arg("refMol"),
               python::arg("prbPyMMFFMolProperties") = python::object(),
               python::arg("refPyMMFFMolProperties") = python::object(),
               python::arg("prbCid") = -1, python::arg("refCid") = -1,
               python::arg("reflect") = false, python::arg("maxIters") = 50,
               python::arg("options") = 0,
               python::arg("constraintMap") = python::list(),
               python::arg("constraintWeights") = python::list()),
              python::return_value_policy<python::manage_new_object>(),
              docString.c_str());
  docString =
      "Get an O3A object with atomMap and weights vectors to overlay\n\
      the probe molecule onto the reference molecule based on\n\
      Crippen logP atom contributions\n\
     \n\
     ARGUMENTS\n\
      - prbMol                   molecule that is to be aligned\n\
      - refMol                   molecule used as the reference for the alignment\n\
      - prbCrippenContribs       Crippen atom contributions for the probe molecule\n\
                                 as a list of (logp, mr) tuples, as returned\n\
                                 by _CalcCrippenContribs()\n\
      - refCrippenContribs       Crippen atom contributions for the reference molecule\n\
                                 as a list of (logp, mr) tuples, as returned\n\
                                 by _CalcCrippenContribs()\n\
      - prbCid                   ID of the conformation in the probe to be used \n\
                                 for the alignment (defaults to first conformation)\n\
      - refCid                   ID of the conformation in the ref molecule to which \n\
                                 the alignment is computed (defaults to first conformation)\n\
      - reflect                  if true reflect the conformation of the probe molecule\n\
                                 (defaults to false)\n\
      - maxIters                 maximum number of iterations used in minimizing the RMSD\n\
                                 (defaults to 50)\n\
      - options                  least 2 significant bits encode accuracy\n\
                                 (0: maximum, 3: minimum; defaults to 0)\n\
                                 bit 3 triggers local optimization of the alignment\n\
                                 (no computation of the cost matrix; defaults: off)\n\
      - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)\n\
                                 which shall be used for the alignment (defaults to [])\n\
      - constraintWeights        optionally specify weights for each of the constraints\n\
                                 (weights default to 100.0)\n\
       \n\
      RETURNS\n\
      The O3A object\n\
    \n";
  python::def("GetCrippenO3A", RDKit::MolAlign::getCrippenO3A,
              (python::arg("prbMol"), python::arg("refMol"),
               python::arg("prbCrippenContribs") = python::list(),
               python::arg("refCrippenContribs") = python::list(),
               python::arg("prbCid") = -1, python::arg("refCid") = -1,
               python::arg("reflect") = false, python::arg("maxIters") = 50,
               python::arg("options") = 0,
               python::arg("constraintMap") = python::list(),
               python::arg("constraintWeights") = python::list()),
              python::return_value_policy<python::manage_new_object>(),
              docString.c_str());

  docString =
      "Get a vector of O3A objects for the overlay of all \n\
      the probe molecule's conformations onto the reference molecule based on\n\
      MMFF atom types and charges\n\
     \n\
     ARGUMENTS\n\
      - prbMol                   molecule that is to be aligned\n\
      - refMol                   molecule used as the reference for the alignment\n\
      - numThreads :             the number of threads to use, only has an effect if\n\
                                 the RDKit was built with thread support (defaults to 1)\n\
                                 If set to zero, the max supported by the system will be used.\n\
      - prbPyMMFFMolProperties   PyMMFFMolProperties object for the probe molecule as returned\n\
                                 by SetupMMFFForceField()\n\
      - refPyMMFFMolProperties   PyMMFFMolProperties object for the reference molecule as returned\n\
                                 by SetupMMFFForceField()\n\
      - refCid                   ID of the conformation in the ref molecule to which \n\
                                 the alignment is computed (defaults to first conformation)\n\
      - reflect                  if true reflect the conformation of the probe molecule\n\
                                 (defaults to false)\n\
      - maxIters                 maximum number of iterations used in minimizing the RMSD\n\
                                 (defaults to 50)\n\
      - options                  least 2 significant bits encode accuracy\n\
                                 (0: maximum, 3: minimum; defaults to 0)\n\
                                 bit 3 triggers local optimization of the alignment\n\
                                 (no computation of the cost matrix; defaults: off)\n\
      - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)\n\
                                 which shall be used for the alignment (defaults to [])\n\
      - constraintWeights        optionally specify weights for each of the constraints\n\
                                 (weights default to 100.0)\n\
       \n\
      RETURNS\n\
      A vector of O3A objects\n\
    \n";
  python::def("GetO3AForProbeConfs", RDKit::MolAlign::getMMFFO3AForConfs,
              (python::arg("prbMol"), python::arg("refMol"),
               python::arg("numThreads") = 1,
               python::arg("prbPyMMFFMolProperties") = python::object(),
               python::arg("refPyMMFFMolProperties") = python::object(),
               python::arg("refCid") = -1, python::arg("reflect") = false,
               python::arg("maxIters") = 50, python::arg("options") = 0,
               python::arg("constraintMap") = python::list(),
               python::arg("constraintWeights") = python::list()),
              docString.c_str());

  docString =
      "Get a vector of O3A objects for the overlay of all \n\
      the probe molecule's conformations onto the reference molecule based on\n\
      MMFF atom types and charges\n\
     \n\
     ARGUMENTS\n\
      - prbMol                   molecule that is to be aligned\n\
      - refMol                   molecule used as the reference for the alignment\n\
      - numThreads :             the number of threads to use, only has an effect if\n\
                                 the RDKit was built with thread support (defaults to 1)\n\
      - prbCrippenContribs       Crippen atom contributions for the probe molecule\n\
                                 as a list of (logp, mr) tuples, as returned\n\
                                 by _CalcCrippenContribs()\n\
      - refCrippenContribs       Crippen atom contributions for the reference molecule\n\
                                 as a list of (logp, mr) tuples, as returned\n\
                                 by _CalcCrippenContribs()\n\
      - refCid                   ID of the conformation in the ref molecule to which \n\
                                 the alignment is computed (defaults to first conformation)\n\
      - reflect                  if true reflect the conformation of the probe molecule\n\
                                 (defaults to false)\n\
      - maxIters                 maximum number of iterations used in minimizing the RMSD\n\
                                 (defaults to 50)\n\
      - options                  least 2 significant bits encode accuracy\n\
                                 (0: maximum, 3: minimum; defaults to 0)\n\
                                 bit 3 triggers local optimization of the alignment\n\
                                 (no computation of the cost matrix; defaults: off)\n\
      - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)\n\
                                 which shall be used for the alignment (defaults to [])\n\
      - constraintWeights        optionally specify weights for each of the constraints\n\
                                 (weights default to 100.0)\n\
       \n\
      RETURNS\n\
      A vector of O3A objects\n\
    \n";
  python::def("GetCrippenO3AForProbeConfs",
              RDKit::MolAlign::getCrippenO3AForConfs,
              (python::arg("prbMol"), python::arg("refMol"),
               python::arg("numThreads") = 1,
               python::arg("prbCrippenContribs") = python::list(),
               python::arg("refCrippenContribs") = python::list(),
               python::arg("refCid") = -1, python::arg("reflect") = false,
               python::arg("maxIters") = 50, python::arg("options") = 0,
               python::arg("constraintMap") = python::list(),
               python::arg("constraintWeights") = python::list()),
              docString.c_str());
}
