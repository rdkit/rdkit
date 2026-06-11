//
//  Copyright (C) 2004-2026 Greg Landrum, Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolAlign/O3AAlignMolecules.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace {

// -- Sequence translation helpers --

MatchVectType nbTranslateAtomMap(nb::object obj) {
  MatchVectType result;
  if (obj.is_none()) return result;
  nb::sequence seq = nb::cast<nb::sequence>(obj);
  size_t n = nb::len(seq);
  for (size_t i = 0; i < n; ++i) {
    nb::sequence pair = nb::cast<nb::sequence>(seq[i]);
    if (nb::len(pair) != 2) {
      throw nb::value_error("Incorrect format for an atomMap");
    }
    result.push_back({nb::cast<int>(pair[0]), nb::cast<int>(pair[1])});
  }
  return result;
}

std::vector<MatchVectType> nbTranslateAtomMapSeq(nb::object obj) {
  std::vector<MatchVectType> result;
  if (obj.is_none()) return result;
  nb::sequence seq = nb::cast<nb::sequence>(obj);
  size_t n = nb::len(seq);
  for (size_t i = 0; i < n; ++i) {
    result.push_back(nbTranslateAtomMap(nb::cast<nb::object>(seq[i])));
  }
  return result;
}

std::vector<double> nbTranslateDoubleVec(nb::object obj) {
  std::vector<double> result;
  if (obj.is_none()) return result;
  nb::sequence seq = nb::cast<nb::sequence>(obj);
  size_t n = nb::len(seq);
  for (size_t i = 0; i < n; ++i) {
    result.push_back(nb::cast<double>(seq[i]));
  }
  return result;
}

std::vector<unsigned int> nbTranslateUIntVec(nb::object obj) {
  std::vector<unsigned int> result;
  if (obj.is_none()) return result;
  nb::sequence seq = nb::cast<nb::sequence>(obj);
  size_t n = nb::len(seq);
  for (size_t i = 0; i < n; ++i) {
    result.push_back(nb::cast<unsigned int>(seq[i]));
  }
  return result;
}

std::unique_ptr<RDNumeric::DoubleVector> makeDoubleVector(nb::object obj) {
  if (obj.is_none()) return nullptr;
  nb::sequence seq = nb::cast<nb::sequence>(obj);
  size_t n = nb::len(seq);
  if (n == 0) return nullptr;
  auto dv = std::make_unique<RDNumeric::DoubleVector>(n);
  for (size_t i = 0; i < n; ++i) {
    dv->setVal(i, nb::cast<double>(seq[i]));
  }
  return dv;
}

// -- Return value helpers --

auto makeTransform4x4(const RDGeom::Transform3D &trans) {
  double *resData = new double[16];
  memcpy(resData, trans.getData(), 16 * sizeof(double));
  nb::capsule owner(resData, [](void *f) noexcept {
    delete[] reinterpret_cast<double *>(f);
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<2>>(resData, {4, 4}, owner);
}

nb::object makeRmsdTransResult(double rmsd, const RDGeom::Transform3D &trans,
                               const MatchVectType *match = nullptr) {
  auto transArr = makeTransform4x4(trans);
  if (!match) {
    return nb::make_tuple(rmsd, transArr);
  }
  nb::list matchList;
  for (const auto &p : *match) {
    matchList.append(nb::make_tuple(p.first, p.second));
  }
  return nb::make_tuple(rmsd, transArr, matchList);
}

// -- Constraint parsing and validation --

std::pair<std::unique_ptr<MatchVectType>, std::unique_ptr<RDNumeric::DoubleVector>>
parseConstraints(nb::object constraintMap, nb::object constraintWeights,
                 const ROMol &prbMol, const ROMol &refMol) {
  std::unique_ptr<MatchVectType> cMap;
  std::unique_ptr<RDNumeric::DoubleVector> cWts;
  if (!constraintMap.is_none() && nb::len(constraintMap) > 0) {
    cMap = std::make_unique<MatchVectType>(nbTranslateAtomMap(constraintMap));
    cWts = makeDoubleVector(constraintWeights);
    if (cWts && cMap->size() != cWts->size()) {
      throw nb::value_error(
          "The number of weights should match the number of constraints");
    }
    for (const auto &p : *cMap) {
      if (p.first < 0 || p.first >= static_cast<int>(prbMol.getNumAtoms()) ||
          p.second < 0 || p.second >= static_cast<int>(refMol.getNumAtoms())) {
        throw nb::value_error("Constrained atom idx out of range");
      }
      if (prbMol[p.first]->getAtomicNum() == 1 ||
          refMol[p.second]->getAtomicNum() == 1) {
        throw nb::value_error("Constrained atoms must be heavy atoms");
      }
    }
  }
  return {std::move(cMap), std::move(cWts)};
}

// -- NbMMFFMolProperties --

struct NbMMFFMolProperties {
  std::shared_ptr<MMFF::MMFFMolProperties> dp_props;
  bool isValid() const { return dp_props != nullptr; }
};

// -- NbBestAlignmentParams --

struct NbBestAlignmentParams {
  int maxMatches = 1000000;
  bool symmetrizeConjugatedTerminalGroups = true;
  bool ignoreHs = true;
  int numThreads = 1;
  std::vector<MatchVectType> map;
  std::vector<double> weights;

  std::pair<MolAlign::BestAlignmentParams,
            std::unique_ptr<RDNumeric::DoubleVector>>
  toNative() const {
    MolAlign::BestAlignmentParams params;
    params.maxMatches = maxMatches;
    params.symmetrizeConjugatedTerminalGroups = symmetrizeConjugatedTerminalGroups;
    params.ignoreHs = ignoreHs;
    params.numThreads = numThreads;
    params.map = map;
    std::unique_ptr<RDNumeric::DoubleVector> weightsVec;
    if (!weights.empty()) {
      weightsVec = std::make_unique<RDNumeric::DoubleVector>(weights.size());
      for (unsigned int i = 0; i < weights.size(); ++i) {
        weightsVec->setVal(i, weights[i]);
      }
      params.weights = weightsVec.get();
    }
    return {params, std::move(weightsVec)};
  }
};

// -- NbO3A --

struct NbO3A {
  std::shared_ptr<MolAlign::O3A> o3a;

  double align() { return o3a->align(); }

  nb::object trans() {
    RDGeom::Transform3D t;
    double rmsd = o3a->trans(t);
    return makeRmsdTransResult(rmsd, t);
  }

  double score() { return o3a->score(); }

  nb::list matches() {
    nb::list result;
    const MatchVectType *m = o3a->matches();
    for (const auto &p : *m) {
      nb::list pair;
      pair.append(p.first);
      pair.append(p.second);
      result.append(pair);
    }
    return result;
  }

  nb::list weights() {
    nb::list result;
    const RDNumeric::DoubleVector *w = o3a->weights();
    for (unsigned int i = 0; i < w->size(); ++i) {
      result.append((*w)[i]);
    }
    return result;
  }
};

// -- Module functions --

nb::object getMolAlignTransform(const ROMol &prbMol, const ROMol &refMol,
                                int prbCid, int refCid, nb::object atomMap,
                                nb::object weights, bool reflect,
                                unsigned int maxIters) {
  auto aMap = nbTranslateAtomMap(atomMap);
  MatchVectType *aMapPtr = aMap.empty() ? nullptr : &aMap;
  auto wtsVec = makeDoubleVector(weights);
  if (wtsVec) {
    unsigned int nAtms = aMapPtr ? aMapPtr->size() : prbMol.getNumAtoms();
    if (wtsVec->size() != nAtms) {
      throw nb::value_error("Incorrect number of weights specified");
    }
  }
  RDGeom::Transform3D trans;
  double rmsd;
  {
    nb::gil_scoped_release release;
    rmsd = MolAlign::getAlignmentTransform(prbMol, refMol, trans, prbCid,
                                           refCid, aMapPtr, wtsVec.get(),
                                           reflect, maxIters);
  }
  return makeRmsdTransResult(rmsd, trans);
}

nb::object getBestMolAlignTransform(const ROMol &prbMol, const ROMol &refMol,
                                    int prbCid, int refCid, nb::object map,
                                    int maxMatches, bool symmetrize,
                                    nb::object weights, bool reflect,
                                    unsigned int maxIters, int numThreads) {
  NbBestAlignmentParams nbParams;
  nbParams.maxMatches = maxMatches;
  nbParams.symmetrizeConjugatedTerminalGroups = symmetrize;
  nbParams.ignoreHs = false;
  nbParams.numThreads = numThreads;
  if (!map.is_none()) {
    nbParams.map = nbTranslateAtomMapSeq(map);
  }
  if (!weights.is_none()) {
    nbParams.weights = nbTranslateDoubleVec(weights);
  }
  auto [params, weightsOwner] = nbParams.toNative();
  RDGeom::Transform3D bestTrans;
  MatchVectType bestMatch;
  double rmsd;
  {
    nb::gil_scoped_release release;
    rmsd = MolAlign::getBestAlignmentTransform(prbMol, refMol, bestTrans,
                                               bestMatch, params, prbCid,
                                               refCid, reflect, maxIters);
  }
  return makeRmsdTransResult(rmsd, bestTrans, &bestMatch);
}

nb::object getBestMolAlignTransformParams(const ROMol &prbMol,
                                          const ROMol &refMol,
                                          const NbBestAlignmentParams &nbParams,
                                          int prbCid, int refCid, bool reflect,
                                          unsigned int maxIters) {
  auto [params, weightsOwner] = nbParams.toNative();
  RDGeom::Transform3D bestTrans;
  MatchVectType bestMatch;
  double rmsd;
  {
    nb::gil_scoped_release release;
    rmsd = MolAlign::getBestAlignmentTransform(prbMol, refMol, bestTrans,
                                               bestMatch, params, prbCid,
                                               refCid, reflect, maxIters);
  }
  return makeRmsdTransResult(rmsd, bestTrans, &bestMatch);
}

double alignMolecule(ROMol &prbMol, const ROMol &refMol, int prbCid,
                     int refCid, nb::object atomMap, nb::object weights,
                     bool reflect, unsigned int maxIters) {
  auto aMap = nbTranslateAtomMap(atomMap);
  MatchVectType *aMapPtr = aMap.empty() ? nullptr : &aMap;
  auto wtsVec = makeDoubleVector(weights);
  if (wtsVec) {
    unsigned int nAtms = aMapPtr ? aMapPtr->size() : prbMol.getNumAtoms();
    if (wtsVec->size() != nAtms) {
      throw nb::value_error("Incorrect number of weights specified");
    }
  }
  double rmsd;
  {
    nb::gil_scoped_release release;
    rmsd = MolAlign::alignMol(prbMol, refMol, prbCid, refCid, aMapPtr,
                              wtsVec.get(), reflect, maxIters);
  }
  return rmsd;
}

double getBestRMS(ROMol &prbMol, ROMol &refMol, int prbId, int refId,
                  nb::object map, int maxMatches, bool symmetrize,
                  nb::object weights, int numThreads) {
  NbBestAlignmentParams nbParams;
  nbParams.maxMatches = maxMatches;
  nbParams.symmetrizeConjugatedTerminalGroups = symmetrize;
  nbParams.ignoreHs = false;
  nbParams.numThreads = numThreads;
  if (!map.is_none()) {
    nbParams.map = nbTranslateAtomMapSeq(map);
  }
  if (!weights.is_none()) {
    nbParams.weights = nbTranslateDoubleVec(weights);
  }
  auto [params, weightsOwner] = nbParams.toNative();
  double rmsd;
  {
    nb::gil_scoped_release release;
    rmsd = MolAlign::getBestRMS(prbMol, refMol, params, prbId, refId);
  }
  return rmsd;
}

double getBestRMSParams(ROMol &prbMol, ROMol &refMol,
                        const NbBestAlignmentParams &nbParams, int prbId,
                        int refId) {
  auto [params, weightsOwner] = nbParams.toNative();
  double rmsd;
  {
    nb::gil_scoped_release release;
    rmsd = MolAlign::getBestRMS(prbMol, refMol, params, prbId, refId);
  }
  return rmsd;
}

nb::tuple getAllConformerBestRMS(ROMol &mol, int numThreads, nb::object map,
                                int maxMatches, bool symmetrize,
                                nb::object weights) {
  NbBestAlignmentParams nbParams;
  nbParams.maxMatches = maxMatches;
  nbParams.symmetrizeConjugatedTerminalGroups = symmetrize;
  nbParams.ignoreHs = true;
  nbParams.numThreads = numThreads;
  if (!map.is_none()) {
    nbParams.map = nbTranslateAtomMapSeq(map);
  }
  if (!weights.is_none()) {
    nbParams.weights = nbTranslateDoubleVec(weights);
  }
  auto [params, weightsOwner] = nbParams.toNative();
  std::vector<double> rmsds;
  {
    nb::gil_scoped_release release;
    rmsds = MolAlign::getAllConformerBestRMS(mol, params);
  }
  nb::list res;
  for (double v : rmsds) {
    res.append(v);
  }
  return nb::tuple(res);
}

nb::tuple getAllConformerBestRMSParams(ROMol &mol,
                                      const NbBestAlignmentParams &nbParams) {
  auto [params, weightsOwner] = nbParams.toNative();
  std::vector<double> rmsds;
  {
    nb::gil_scoped_release release;
    rmsds = MolAlign::getAllConformerBestRMS(mol, params);
  }
  nb::list res;
  for (double v : rmsds) {
    res.append(v);
  }
  return nb::tuple(res);
}

double calcRMS(ROMol &prbMol, ROMol &refMol, int prbCid, int refCid,
               nb::object map, int maxMatches, bool symmetrize,
               nb::object weights) {
  std::vector<MatchVectType> aMapVec;
  if (!map.is_none()) {
    aMapVec = nbTranslateAtomMapSeq(map);
  }
  auto wtsVec = makeDoubleVector(weights);
  double rmsd;
  {
    nb::gil_scoped_release release;
    rmsd = MolAlign::CalcRMS(prbMol, refMol, prbCid, refCid, aMapVec,
                             maxMatches, symmetrize, wtsVec.get());
  }
  return rmsd;
}

void alignMolConfs(ROMol &mol, nb::object atomIds, nb::object confIds,
                   nb::object weights, bool reflect, unsigned int maxIters,
                   nb::object RMSlist) {
  auto aIds = nbTranslateUIntVec(atomIds);
  auto cIds = nbTranslateUIntVec(confIds);
  auto wtsVec = makeDoubleVector(weights);
  const std::vector<unsigned int> *aIdsPtr = aIds.empty() ? nullptr : &aIds;
  const std::vector<unsigned int> *cIdsPtr = cIds.empty() ? nullptr : &cIds;
  std::unique_ptr<std::vector<double>> RMSvector;
  if (!RMSlist.is_none()) {
    RMSvector = std::make_unique<std::vector<double>>();
  }
  {
    nb::gil_scoped_release release;
    MolAlign::alignMolConformers(mol, aIdsPtr, cIdsPtr, wtsVec.get(), reflect,
                                 maxIters, RMSvector.get());
  }
  if (RMSvector) {
    for (double v : *RMSvector) {
      RMSlist.attr("append")(v);
    }
  }
}

NbO3A getMMFFO3A(ROMol &prbMol, ROMol &refMol, nb::object prbProps,
                 nb::object refProps, int prbCid, int refCid, bool reflect,
                 unsigned int maxIters, unsigned int options,
                 nb::object constraintMap, nb::object constraintWeights) {
  auto [cMap, cWts] =
      parseConstraints(constraintMap, constraintWeights, prbMol, refMol);
  MMFF::MMFFMolProperties *prbMolPropsPtr = nullptr;
  MMFF::MMFFMolProperties *refMolPropsPtr = nullptr;
  std::unique_ptr<MMFF::MMFFMolProperties> prbMolProps;
  std::unique_ptr<MMFF::MMFFMolProperties> refMolProps;

  if (!prbProps.is_none()) {
    prbMolPropsPtr = nb::cast<NbMMFFMolProperties &>(prbProps).dp_props.get();
  } else {
    prbMolProps = std::make_unique<MMFF::MMFFMolProperties>(prbMol);
    if (!prbMolProps->isValid()) {
      throw nb::value_error("missing MMFF94 parameters for probe molecule");
    }
    prbMolPropsPtr = prbMolProps.get();
  }
  if (!refProps.is_none()) {
    refMolPropsPtr = nb::cast<NbMMFFMolProperties &>(refProps).dp_props.get();
  } else {
    refMolProps = std::make_unique<MMFF::MMFFMolProperties>(refMol);
    if (!refMolProps->isValid()) {
      throw nb::value_error("missing MMFF94 parameters for reference molecule");
    }
    refMolPropsPtr = refMolProps.get();
  }

  std::shared_ptr<MolAlign::O3A> o3a;
  {
    nb::gil_scoped_release release;
    o3a = std::shared_ptr<MolAlign::O3A>(
        new MolAlign::O3A(prbMol, refMol, prbMolPropsPtr, refMolPropsPtr,
                          MolAlign::O3A::MMFF94, prbCid, refCid, reflect,
                          maxIters, options, cMap.get(), cWts.get()));
  }
  return NbO3A{std::move(o3a)};
}

nb::tuple getMMFFO3AForConfs(ROMol &prbMol, ROMol &refMol, int numThreads,
                             nb::object prbProps, nb::object refProps,
                             int refCid, bool reflect, unsigned int maxIters,
                             unsigned int options, nb::object constraintMap,
                             nb::object constraintWeights) {
  auto [cMap, cWts] =
      parseConstraints(constraintMap, constraintWeights, prbMol, refMol);
  MMFF::MMFFMolProperties *prbMolPropsPtr = nullptr;
  MMFF::MMFFMolProperties *refMolPropsPtr = nullptr;
  std::unique_ptr<MMFF::MMFFMolProperties> prbMolProps;
  std::unique_ptr<MMFF::MMFFMolProperties> refMolProps;

  if (!prbProps.is_none()) {
    prbMolPropsPtr = nb::cast<NbMMFFMolProperties &>(prbProps).dp_props.get();
  } else {
    prbMolProps = std::make_unique<MMFF::MMFFMolProperties>(prbMol);
    if (!prbMolProps->isValid()) {
      throw nb::value_error("missing MMFF94 parameters for probe molecule");
    }
    prbMolPropsPtr = prbMolProps.get();
  }
  if (!refProps.is_none()) {
    refMolPropsPtr = nb::cast<NbMMFFMolProperties &>(refProps).dp_props.get();
  } else {
    refMolProps = std::make_unique<MMFF::MMFFMolProperties>(refMol);
    if (!refMolProps->isValid()) {
      throw nb::value_error("missing MMFF94 parameters for reference molecule");
    }
    refMolPropsPtr = refMolProps.get();
  }

  std::vector<boost::shared_ptr<MolAlign::O3A>> res;
  {
    nb::gil_scoped_release release;
    MolAlign::getO3AForProbeConfs(prbMol, refMol, prbMolPropsPtr,
                                  refMolPropsPtr, res, numThreads,
                                  MolAlign::O3A::MMFF94, refCid, reflect,
                                  maxIters, options, cMap.get(), cWts.get());
  }
  nb::list pyres;
  for (auto &i : res) {
    pyres.append(NbO3A{std::shared_ptr<MolAlign::O3A>(i.get(), [b = i](MolAlign::O3A *) {})});
  }
  return nb::tuple(pyres);
}

NbO3A getCrippenO3A(ROMol &prbMol, ROMol &refMol,
                    nb::object prbCrippenContribs,
                    nb::object refCrippenContribs, int prbCid, int refCid,
                    bool reflect, unsigned int maxIters, unsigned int options,
                    nb::object constraintMap, nb::object constraintWeights) {
  auto [cMap, cWts] =
      parseConstraints(constraintMap, constraintWeights, prbMol, refMol);
  unsigned int prbNAtoms = prbMol.getNumAtoms();
  unsigned int refNAtoms = refMol.getNumAtoms();
  std::vector<double> prbLogpContribs(prbNAtoms);
  std::vector<double> refLogpContribs(refNAtoms);

  if (!prbCrippenContribs.is_none() &&
      nb::len(prbCrippenContribs) == prbNAtoms) {
    nb::sequence prbSeq = nb::cast<nb::sequence>(prbCrippenContribs);
    for (unsigned int i = 0; i < prbNAtoms; ++i) {
      nb::sequence tup = nb::cast<nb::sequence>(prbSeq[i]);
      prbLogpContribs[i] = nb::cast<double>(tup[0]);
    }
  } else {
    std::vector<double> prbMRContribs(prbNAtoms);
    std::vector<unsigned int> prbAtomTypes(prbNAtoms);
    std::vector<std::string> prbAtomTypeLabels(prbNAtoms);
    Descriptors::getCrippenAtomContribs(prbMol, prbLogpContribs, prbMRContribs,
                                        true, &prbAtomTypes,
                                        &prbAtomTypeLabels);
  }
  if (!refCrippenContribs.is_none() &&
      nb::len(refCrippenContribs) == refNAtoms) {
    nb::sequence refSeq = nb::cast<nb::sequence>(refCrippenContribs);
    for (unsigned int i = 0; i < refNAtoms; ++i) {
      nb::sequence tup = nb::cast<nb::sequence>(refSeq[i]);
      refLogpContribs[i] = nb::cast<double>(tup[0]);
    }
  } else {
    std::vector<double> refMRContribs(refNAtoms);
    std::vector<unsigned int> refAtomTypes(refNAtoms);
    std::vector<std::string> refAtomTypeLabels(refNAtoms);
    Descriptors::getCrippenAtomContribs(refMol, refLogpContribs, refMRContribs,
                                        true, &refAtomTypes,
                                        &refAtomTypeLabels);
  }

  std::shared_ptr<MolAlign::O3A> o3a;
  {
    nb::gil_scoped_release release;
    o3a = std::shared_ptr<MolAlign::O3A>(
        new MolAlign::O3A(prbMol, refMol, &prbLogpContribs, &refLogpContribs,
                          MolAlign::O3A::CRIPPEN, prbCid, refCid, reflect,
                          maxIters, options, cMap.get(), cWts.get()));
  }
  return NbO3A{std::move(o3a)};
}

nb::tuple getCrippenO3AForConfs(ROMol &prbMol, ROMol &refMol, int numThreads,
                                nb::object prbCrippenContribs,
                                nb::object refCrippenContribs, int refCid,
                                bool reflect, unsigned int maxIters,
                                unsigned int options, nb::object constraintMap,
                                nb::object constraintWeights) {
  auto [cMap, cWts] =
      parseConstraints(constraintMap, constraintWeights, prbMol, refMol);
  unsigned int prbNAtoms = prbMol.getNumAtoms();
  unsigned int refNAtoms = refMol.getNumAtoms();
  std::vector<double> prbLogpContribs(prbNAtoms);
  std::vector<double> refLogpContribs(refNAtoms);

  if (!prbCrippenContribs.is_none() &&
      nb::len(prbCrippenContribs) == prbNAtoms) {
    nb::sequence prbSeq = nb::cast<nb::sequence>(prbCrippenContribs);
    for (unsigned int i = 0; i < prbNAtoms; ++i) {
      nb::sequence tup = nb::cast<nb::sequence>(prbSeq[i]);
      prbLogpContribs[i] = nb::cast<double>(tup[0]);
    }
  } else {
    std::vector<double> prbMRContribs(prbNAtoms);
    std::vector<unsigned int> prbAtomTypes(prbNAtoms);
    std::vector<std::string> prbAtomTypeLabels(prbNAtoms);
    Descriptors::getCrippenAtomContribs(prbMol, prbLogpContribs, prbMRContribs,
                                        true, &prbAtomTypes,
                                        &prbAtomTypeLabels);
  }
  if (!refCrippenContribs.is_none() &&
      nb::len(refCrippenContribs) == refNAtoms) {
    nb::sequence refSeq = nb::cast<nb::sequence>(refCrippenContribs);
    for (unsigned int i = 0; i < refNAtoms; ++i) {
      nb::sequence tup = nb::cast<nb::sequence>(refSeq[i]);
      refLogpContribs[i] = nb::cast<double>(tup[0]);
    }
  } else {
    std::vector<double> refMRContribs(refNAtoms);
    std::vector<unsigned int> refAtomTypes(refNAtoms);
    std::vector<std::string> refAtomTypeLabels(refNAtoms);
    Descriptors::getCrippenAtomContribs(refMol, refLogpContribs, refMRContribs,
                                        true, &refAtomTypes,
                                        &refAtomTypeLabels);
  }

  std::vector<boost::shared_ptr<MolAlign::O3A>> res;
  {
    nb::gil_scoped_release release;
    MolAlign::getO3AForProbeConfs(prbMol, refMol, &prbLogpContribs,
                                  &refLogpContribs, res, numThreads,
                                  MolAlign::O3A::CRIPPEN, refCid, reflect,
                                  maxIters, options, cMap.get(), cWts.get());
  }
  nb::list pyres;
  for (auto &i : res) {
    pyres.append(NbO3A{std::shared_ptr<MolAlign::O3A>(i.get(), [b = i](MolAlign::O3A *) {})});
  }
  return nb::tuple(pyres);
}

NbMMFFMolProperties mmffGetMoleculeProperties(ROMol &mol,
                                              const std::string &mmffVariant,
                                              unsigned int mmffVerbosity) {
  auto *p = new MMFF::MMFFMolProperties(mol, mmffVariant, mmffVerbosity);
  NbMMFFMolProperties result;
  if (p->isValid()) {
    result.dp_props.reset(p);
  } else {
    delete p;
  }
  return result;
}

}  // namespace

NB_MODULE(rdMolAlign, m) {
  m.doc() = "Module containing functions to align a molecule to a second molecule";

  nb::class_<NbMMFFMolProperties>(m, "MMFFMolProperties",
                                  "MMFF molecular properties for O3A alignment")
      .def("IsValid", &NbMMFFMolProperties::isValid,
           "Returns True if MMFF parameters are available for this molecule");

  nb::class_<NbBestAlignmentParams>(m, "BestAlignmentParams",
                                    "Parameters controlling RMSD alignment")
      .def(nb::init<>())
      .def_rw("maxMatches", &NbBestAlignmentParams::maxMatches,
              "maximum number of substructure matches to consider")
      .def_rw("symmetrizeConjugatedTerminalGroups",
              &NbBestAlignmentParams::symmetrizeConjugatedTerminalGroups,
              R"DOC(if true, conjugated terminal functional groups (like nitro or carboxylate)
will be considered symmetrically.)DOC")
      .def_rw("ignoreHs", &NbBestAlignmentParams::ignoreHs,
              "if true, hydrogens will be ignored in the alignment")
      .def_rw("numThreads", &NbBestAlignmentParams::numThreads,
              "number of threads to use")
      .def_prop_rw(
          "map",
          [](const NbBestAlignmentParams &p) {
            nb::list result;
            for (const auto &matchVect : p.map) {
              nb::list matchList;
              for (const auto &pair : matchVect) {
                matchList.append(nb::make_tuple(pair.first, pair.second));
              }
              result.append(matchList);
            }
            return nb::tuple(result);
          },
          [](NbBestAlignmentParams &p, nb::object obj) {
            p.map = nbTranslateAtomMapSeq(obj);
          },
          "the atom-atom mapping(s) used in the alignment")
      .def_prop_rw(
          "weights",
          [](const NbBestAlignmentParams &p) {
            nb::list result;
            for (double w : p.weights) {
              result.append(w);
            }
            return nb::tuple(result);
          },
          [](NbBestAlignmentParams &p, nb::object obj) {
            p.weights = nbTranslateDoubleVec(obj);
          },
          "the weights used in the alignment");

  nb::class_<NbO3A>(m, "O3A", "Open3DALIGN object")
      .def("Align", &NbO3A::align,
           "aligns probe molecule onto reference molecule")
      .def("Trans", &NbO3A::trans,
           "returns the transformation which aligns probe molecule onto "
           "reference molecule")
      .def("Score", &NbO3A::score, "returns the O3AScore of the alignment")
      .def("Matches", &NbO3A::matches,
           "returns the AtomMap as found by Open3DALIGN")
      .def("Weights", &NbO3A::weights,
           "returns the weight vector as found by Open3DALIGN");

  m.def(
      "MMFFGetMoleculeProperties", &mmffGetMoleculeProperties, "mol"_a,
      "mmffVariant"_a = "MMFF94", "mmffVerbosity"_a = 0u,
      R"DOC(Get MMFF molecule properties for use with GetO3A.

ARGUMENTS:
 - mol          : the molecule of interest
 - mmffVariant  : "MMFF94" or "MMFF94s" (defaults to "MMFF94")
 - mmffVerbosity: verbosity level (defaults to 0)

RETURNS:
An MMFFMolProperties object, or one with IsValid()==False if parameters are unavailable.)DOC");

  m.def(
      "GetAlignmentTransform", getMolAlignTransform,
      "prbMol"_a, "refMol"_a, "prbCid"_a = -1, "refCid"_a = -1,
      "atomMap"_a = nb::none(), "weights"_a = nb::none(),
      "reflect"_a = false, "maxIters"_a = 50u,
      R"DOC(Compute the transformation required to align a molecule

The 3D transformation required to align the specied conformation in the probe molecule
to a specified conformation in the reference molecule is computed so that the root mean
squared distance between a specified set of atoms is minimized

ARGUMENTS
 - prbMol    molecule that is to be aligned
 - refMol    molecule used as the reference for the alignment
 - prbCid    ID of the conformation in the probe to be used
                  for the alignment (defaults to first conformation)
 - refCid    ID of the conformation in the ref molecule to which
                  the alignment is computed (defaults to first conformation)
 - atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                  used to compute the alignments. If this mapping is
                  not specified an attempt is made to generate one by
                  substructure matching
 - weights   Optionally specify weights for each of the atom pairs
 - reflect   if true reflect the conformation of the probe molecule
 - maxIters  maximum number of iterations used in minimizing the RMSD

RETURNS
a tuple of (RMSD value, transform matrix))DOC");

  m.def(
      "GetBestAlignmentTransform", getBestMolAlignTransform,
      "prbMol"_a, "refMol"_a, "prbCid"_a = -1, "refCid"_a = -1,
      "map"_a = nb::none(), "maxMatches"_a = 1000000,
      "symmetrizeConjugatedTerminalGroups"_a = true,
      "weights"_a = nb::none(), "reflect"_a = false,
      "maxIters"_a = 50u, "numThreads"_a = 1,
      R"DOC(Compute the optimal RMS, transformation and atom map for aligning
two molecules, taking symmetry into account. Molecule coordinates
are left unaltered.

This function will attempt to align all permutations of matching atom
orders in both molecules, for some molecules it will lead to 'combinatorial
explosion' especially if hydrogens are present.
Use 'GetAlignmentTransform' to align molecules without changing the atom order.

ARGUMENTS
 - prbMol      molecule that is to be aligned
 - refMol      molecule used as the reference for the alignment
 - prbCid      ID of the conformation in the probe to be used
               for the alignment (defaults to first conformation)
 - refCid      ID of the conformation in the ref molecule to which
               the alignment is computed (defaults to first conformation)
 - map:        (optional) a list of lists of (probeAtomId, refAtomId)
               tuples with the atom-atom mappings of the two
               molecules. If not provided, these will be generated
               using a substructure search.
 - maxMatches  (optional) if atomMap is empty, this will be the max number of
               matches found in a SubstructMatch().
 - symmetrizeConjugatedTerminalGroups (optional) if set, conjugated
               terminal functional groups (like nitro or carboxylate)
               will be considered symmetrically.
 - weights     Optionally specify weights for each of the atom pairs
 - reflect     if true reflect the conformation of the probe molecule
 - maxIters    maximum number of iterations used in minimizing the RMSD
 - numThreads  (optional) number of threads to use

RETURNS
a tuple of (RMSD value, best transform matrix, best atom map))DOC");

  m.def(
      "GetBestAlignmentTransform", getBestMolAlignTransformParams,
      "prbMol"_a, "refMol"_a, "params"_a, "prbCid"_a = -1, "refCid"_a = -1,
      "reflect"_a = false, "maxIters"_a = 50u,
      R"DOC(Compute the optimal RMS, transformation and atom map for aligning
two molecules using a BestAlignmentParams object.

RETURNS
a tuple of (RMSD value, best transform matrix, best atom map))DOC");

  m.def(
      "AlignMol", alignMolecule,
      "prbMol"_a, "refMol"_a, "prbCid"_a = -1, "refCid"_a = -1,
      "atomMap"_a = nb::none(), "weights"_a = nb::none(),
      "reflect"_a = false, "maxIters"_a = 50u,
      R"DOC(Optimally (minimum RMSD) align a molecule to another molecule

The 3D transformation required to align the specied conformation in the probe molecule
to a specified conformation in the reference molecule is computed so that the root mean
squared distance between a specified set of atoms is minimized.
This transform is then applied to the specified conformation in the probe molecule

ARGUMENTS
 - prbMol    molecule that is to be aligned
 - refMol    molecule used as the reference for the alignment
 - prbCid    ID of the conformation in the probe to be used
                  for the alignment (defaults to first conformation)
 - refCid    ID of the conformation in the ref molecule to which
                  the alignment is computed (defaults to first conformation)
 - atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                  used to compute the alignments. If this mapping is
                  not specified an attempt is made to generate one by
                  substructure matching
 - weights   Optionally specify weights for each of the atom pairs
 - reflect   if true reflect the conformation of the probe molecule
 - maxIters  maximum number of iterations used in minimizing the RMSD

RETURNS
RMSD value)DOC");

  m.def(
      "GetBestRMS", getBestRMS,
      "prbMol"_a, "refMol"_a, "prbId"_a = -1, "refId"_a = -1,
      "map"_a = nb::none(), "maxMatches"_a = 1000000,
      "symmetrizeConjugatedTerminalGroups"_a = true,
      "weights"_a = nb::none(), "numThreads"_a = 1,
      R"DOC(Returns the optimal RMS for aligning two molecules, taking
symmetry into account. As a side-effect, the probe molecule is
left in the aligned state.

Note:
This function will attempt to align all permutations of matching atom
orders in both molecules, for some molecules it will lead to
'combinatorial explosion' especially if hydrogens are present.
Use 'rdkit.Chem.AllChem.AlignMol' to align molecules without changing
the atom order.

ARGUMENTS
 - prbMol:      the molecule to be aligned to the reference
 - refMol:      the reference molecule
 - prbId:       (optional) probe conformation to use
 - refId:       (optional) reference conformation to use
 - map:         (optional) a list of lists of (probeAtomId,refAtomId)
               tuples with the atom-atom mappings of the two
               molecules. If not provided, these will be generated
               using a substructure search.
 - maxMatches:  (optional) if map isn't specified, this will be
               the max number of matches found in a SubstructMatch()
 - symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugated
               terminal functional groups (like nitro or carboxylate)
               will be considered symmetrically
 - weights:     (optional) weights for mapping
 - numThreads:  (optional) number of threads to use

RETURNS
The best RMSD found)DOC");

  m.def(
      "GetBestRMS", getBestRMSParams,
      "prbMol"_a, "refMol"_a, "params"_a, "prbId"_a = -1, "refId"_a = -1,
      R"DOC(Returns the optimal RMS for aligning two molecules using a BestAlignmentParams object.

RETURNS
The best RMSD found)DOC");

  m.def(
      "GetAllConformerBestRMS", getAllConformerBestRMS,
      "mol"_a, "numThreads"_a = 1, "map"_a = nb::none(),
      "maxMatches"_a = 1000000, "symmetrizeConjugatedTerminalGroups"_a = true,
      "weights"_a = nb::none(),
      R"DOC(Returns the symmetric distance matrix between the conformers of a molecule.
getBestRMS() is used to calculate the inter-conformer distances

ARGUMENTS
 - mol:       the molecule to be considered
 - numThreads:  (optional) number of threads to use
 - map:         (optional) a list of lists of (probeAtomId,refAtomId)
               tuples with the atom-atom mappings of the two
               molecules. If not provided, these will be generated
               using a substructure search.
 - maxMatches:  (optional) if map isn't specified, this will be
               the max number of matches found in a SubstructMatch()
 - symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugated
               terminal functional groups (like nitro or carboxylate)
               will be considered symmetrically
 - weights:     (optional) weights for mapping

RETURNS
A tuple with the best RMSDS. The ordering is [(1,0),(2,0),(2,1),(3,0),... etc])DOC");

  m.def(
      "GetAllConformerBestRMS", getAllConformerBestRMSParams,
      "mol"_a, "params"_a,
      R"DOC(Returns the symmetric distance matrix between the conformers of a molecule
using a BestAlignmentParams object.

RETURNS
A tuple with the best RMSDS. The ordering is [(1,0),(2,0),(2,1),(3,0),... etc])DOC");

  m.def(
      "CalcRMS", calcRMS,
      "prbMol"_a, "refMol"_a, "prbId"_a = -1, "refId"_a = -1,
      "map"_a = nb::none(), "maxMatches"_a = 1000000,
      "symmetrizeConjugatedTerminalGroups"_a = true,
      "weights"_a = nb::none(),
      R"DOC(Returns the RMS between two molecules, taking symmetry into account.
In contrast to getBestRMS, the RMS is computed 'in place', i.e.
probe molecules are not aligned to the reference ahead of the
RMS calculation. This is useful, for example, to compute
the RMSD between docking poses and the co-crystallized ligand.

Note:
This function will attempt to match all permutations of matching atom
orders in both molecules, for some molecules it will lead to
'combinatorial explosion' especially if hydrogens are present.

ARGUMENTS
 - prbMol:      the molecule to be aligned to the reference
 - refMol:      the reference molecule
 - prbCId:      (optional) probe conformation to use
 - refCId:      (optional) reference conformation to use
 - map:         (optional) a list of lists of (probeAtomId, refAtomId)
               tuples with the atom-atom mappings of the two
               molecules. If not provided, these will be generated
               using a substructure search.
 - maxMatches:  (optional) if map isn't specified, this will be
               the max number of matches found in a SubstructMatch()
 - symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugated
               terminal functional groups (like nitro or carboxylate)
               will be considered symmetrically
 - weights:     (optional) weights for mapping

RETURNS
The best RMSD found)DOC");

  m.def(
      "AlignMolConformers", alignMolConfs,
      "mol"_a, "atomIds"_a = nb::none(), "confIds"_a = nb::none(),
      "weights"_a = nb::none(), "reflect"_a = false, "maxIters"_a = 50u,
      "RMSlist"_a = nb::none(),
      R"DOC(Align conformations in a molecule to each other

The first conformation in the molecule is used as the reference

ARGUMENTS
 - mol          molecule of interest
 - atomIds      List of atom ids to use a points for alignment - defaults to all atoms
 - confIds      Ids of conformations to align - defaults to all conformers
 - weights      Optionally specify weights for each of the atom pairs
 - reflect      if true reflect the conformation of the probe molecule
 - maxIters     maximum number of iterations used in minimizing the RMSD
 - RMSlist      if provided, fills in the RMS values between the reference
                conformation and the other aligned conformations)DOC");

  m.def(
      "RandomTransform", MolAlign::randomTransform,
      "mol"_a, "cid"_a = -1, "seed"_a = -1,
      R"DOC(Perform a random transformation on a molecule

ARGUMENTS
 - mol    molecule that is to be transformed
 - cid    ID of the conformation in the mol to be transformed
          (defaults to first conformation)
 - seed   seed used to initialize the random generator
          (defaults to -1, that is no seeding))DOC");

  m.def(
      "GetO3A", getMMFFO3A,
      "prbMol"_a, "refMol"_a,
      "prbPyMMFFMolProperties"_a = nb::none(),
      "refPyMMFFMolProperties"_a = nb::none(),
      "prbCid"_a = -1, "refCid"_a = -1, "reflect"_a = false,
      "maxIters"_a = 50u, "options"_a = 0u,
      "constraintMap"_a = nb::none(), "constraintWeights"_a = nb::none(),
      R"DOC(Get an O3A object with atomMap and weights vectors to overlay
the probe molecule onto the reference molecule based on
MMFF atom types and charges

ARGUMENTS
 - prbMol                   molecule that is to be aligned
 - refMol                   molecule used as the reference for the alignment
 - prbPyMMFFMolProperties   MMFFMolProperties object for the probe molecule as returned
                            by MMFFGetMoleculeProperties()
 - refPyMMFFMolProperties   MMFFMolProperties object for the reference molecule as returned
                            by MMFFGetMoleculeProperties()
 - prbCid                   ID of the conformation in the probe to be used
                            for the alignment (defaults to first conformation)
 - refCid                   ID of the conformation in the ref molecule to which
                            the alignment is computed (defaults to first conformation)
 - reflect                  if true reflect the conformation of the probe molecule
                            (defaults to false)
 - maxIters                 maximum number of iterations used in minimizing the RMSD
                            (defaults to 50)
 - options                  least 2 significant bits encode accuracy
                            (0: maximum, 3: minimum; defaults to 0)
                            bit 3 triggers local optimization of the alignment
                            (no computation of the cost matrix; defaults: off)
 - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                            which shall be used for the alignment (defaults to [])
 - constraintWeights        optionally specify weights for each of the constraints
                            (weights default to 100.0)

RETURNS
The O3A object)DOC");

  m.def(
      "GetCrippenO3A", getCrippenO3A,
      "prbMol"_a, "refMol"_a,
      "prbCrippenContribs"_a = nb::none(),
      "refCrippenContribs"_a = nb::none(),
      "prbCid"_a = -1, "refCid"_a = -1, "reflect"_a = false,
      "maxIters"_a = 50u, "options"_a = 0u,
      "constraintMap"_a = nb::none(), "constraintWeights"_a = nb::none(),
      R"DOC(Get an O3A object with atomMap and weights vectors to overlay
the probe molecule onto the reference molecule based on
Crippen logP atom contributions

ARGUMENTS
 - prbMol                   molecule that is to be aligned
 - refMol                   molecule used as the reference for the alignment
 - prbCrippenContribs       Crippen atom contributions for the probe molecule
                            as a list of (logp, mr) tuples, as returned
                            by _CalcCrippenContribs()
 - refCrippenContribs       Crippen atom contributions for the reference molecule
                            as a list of (logp, mr) tuples, as returned
                            by _CalcCrippenContribs()
 - prbCid                   ID of the conformation in the probe to be used
                            for the alignment (defaults to first conformation)
 - refCid                   ID of the conformation in the ref molecule to which
                            the alignment is computed (defaults to first conformation)
 - reflect                  if true reflect the conformation of the probe molecule
                            (defaults to false)
 - maxIters                 maximum number of iterations used in minimizing the RMSD
                            (defaults to 50)
 - options                  least 2 significant bits encode accuracy
                            (0: maximum, 3: minimum; defaults to 0)
                            bit 3 triggers local optimization of the alignment
                            (no computation of the cost matrix; defaults: off)
 - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                            which shall be used for the alignment (defaults to [])
 - constraintWeights        optionally specify weights for each of the constraints
                            (weights default to 100.0)

RETURNS
The O3A object)DOC");

  m.def(
      "GetO3AForProbeConfs", getMMFFO3AForConfs,
      "prbMol"_a, "refMol"_a, "numThreads"_a = 1,
      "prbPyMMFFMolProperties"_a = nb::none(),
      "refPyMMFFMolProperties"_a = nb::none(),
      "refCid"_a = -1, "reflect"_a = false, "maxIters"_a = 50u,
      "options"_a = 0u, "constraintMap"_a = nb::none(),
      "constraintWeights"_a = nb::none(),
      R"DOC(Get a vector of O3A objects for the overlay of all
the probe molecule's conformations onto the reference molecule based on
MMFF atom types and charges

ARGUMENTS
 - prbMol                   molecule that is to be aligned
 - refMol                   molecule used as the reference for the alignment
 - numThreads :             the number of threads to use, only has an effect if
                            the RDKit was built with thread support (defaults to 1)
                            If set to zero, the max supported by the system will be used.
 - prbPyMMFFMolProperties   MMFFMolProperties object for the probe molecule as returned
                            by MMFFGetMoleculeProperties()
 - refPyMMFFMolProperties   MMFFMolProperties object for the reference molecule as returned
                            by MMFFGetMoleculeProperties()
 - refCid                   ID of the conformation in the ref molecule to which
                            the alignment is computed (defaults to first conformation)
 - reflect                  if true reflect the conformation of the probe molecule
                            (defaults to false)
 - maxIters                 maximum number of iterations used in minimizing the RMSD
                            (defaults to 50)
 - options                  least 2 significant bits encode accuracy
                            (0: maximum, 3: minimum; defaults to 0)
                            bit 3 triggers local optimization of the alignment
                            (no computation of the cost matrix; defaults: off)
 - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                            which shall be used for the alignment (defaults to [])
 - constraintWeights        optionally specify weights for each of the constraints
                            (weights default to 100.0)

RETURNS
A vector of O3A objects)DOC");

  m.def(
      "GetCrippenO3AForProbeConfs", getCrippenO3AForConfs,
      "prbMol"_a, "refMol"_a, "numThreads"_a = 1,
      "prbCrippenContribs"_a = nb::none(),
      "refCrippenContribs"_a = nb::none(),
      "refCid"_a = -1, "reflect"_a = false, "maxIters"_a = 50u,
      "options"_a = 0u, "constraintMap"_a = nb::none(),
      "constraintWeights"_a = nb::none(),
      R"DOC(Get a vector of O3A objects for the overlay of all
the probe molecule's conformations onto the reference molecule based on
Crippen logP atom contributions

ARGUMENTS
 - prbMol                   molecule that is to be aligned
 - refMol                   molecule used as the reference for the alignment
 - numThreads :             the number of threads to use, only has an effect if
                            the RDKit was built with thread support (defaults to 1)
 - prbCrippenContribs       Crippen atom contributions for the probe molecule
                            as a list of (logp, mr) tuples, as returned
                            by _CalcCrippenContribs()
 - refCrippenContribs       Crippen atom contributions for the reference molecule
                            as a list of (logp, mr) tuples, as returned
                            by _CalcCrippenContribs()
 - refCid                   ID of the conformation in the ref molecule to which
                            the alignment is computed (defaults to first conformation)
 - reflect                  if true reflect the conformation of the probe molecule
                            (defaults to false)
 - maxIters                 maximum number of iterations used in minimizing the RMSD
                            (defaults to 50)
 - options                  least 2 significant bits encode accuracy
                            (0: maximum, 3: minimum; defaults to 0)
                            bit 3 triggers local optimization of the alignment
                            (no computation of the cost matrix; defaults: off)
 - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                            which shall be used for the alignment (defaults to [])
 - constraintWeights        optionally specify weights for each of the constraints
                            (weights default to 100.0)

RETURNS
A vector of O3A objects)DOC");
}
