//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <string>
#include <cmath>

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Substruct/SubstructUtils.h>
#include "substructmethods.h"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/GenericGroups/GenericGroups.h>
#include <GraphMol/Subset.h>
#include <RDBoost/Wrap_nb.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/bind_map.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/SmilesParse/CanonicalizeStereoGroups.h>

#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;
NB_MAKE_OPAQUE(std::vector<bool>);
NB_MAKE_OPAQUE(std::map<unsigned int, unsigned int>);

namespace RDKit {

nb::tuple computeAtomCIPRanksHelper(ROMol &mol) {
  UINT_VECT atomRanks;
  Chirality::assignAtomCIPRanks(mol, atomRanks);
  nb::list res;
  for (auto rank : atomRanks) {
    res.append(rank);
  }
  return nb::tuple(res);
}

nb::tuple fragmentOnSomeBondsHelper(const ROMol &mol, nb::object pyBondIndices,
                                    unsigned int nToBreak, bool addDummies,
                                    nb::object pyDummyLabels,
                                    nb::object pyBondTypes,
                                    bool returnCutsPerAtom) {
  auto bondIndices = pythonObjectToVect(pyBondIndices, mol.getNumBonds());
  if (!bondIndices.get()) {
    throw ValueErrorException("empty bond indices");
  }

  std::vector<std::pair<unsigned int, unsigned int>> *dummyLabels = nullptr;
  if (pyDummyLabels) {
    unsigned int nVs = nb::len(pyDummyLabels);
    dummyLabels = new std::vector<std::pair<unsigned int, unsigned int>>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      unsigned int v1 = nb::cast<unsigned int>(pyDummyLabels[i][0]);
      unsigned int v2 = nb::cast<unsigned int>(pyDummyLabels[i][1]);
      (*dummyLabels)[i] = std::make_pair(v1, v2);
    }
  }
  std::vector<Bond::BondType> *bondTypes = nullptr;
  if (pyBondTypes) {
    unsigned int nVs = nb::len(pyBondTypes);
    if (nVs != bondIndices->size()) {
      throw ValueErrorException("bondTypes shorter than bondIndices");
    }
    bondTypes = new std::vector<Bond::BondType>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      (*bondTypes)[i] = nb::cast<Bond::BondType>(pyBondTypes[i]);
    }
  }
  std::vector<std::vector<unsigned int>> *cutsPerAtom = nullptr;
  if (returnCutsPerAtom) {
    cutsPerAtom = new std::vector<std::vector<unsigned int>>;
  }

  std::vector<ROMOL_SPTR> frags;
  MolFragmenter::fragmentOnSomeBonds(mol, *bondIndices, frags, nToBreak,
                                     addDummies, dummyLabels, bondTypes,
                                     cutsPerAtom);
  nb::list res;
  for (auto &frag : frags) {
    res.append(frag);
  }
  delete dummyLabels;
  delete bondTypes;
  if (cutsPerAtom) {
    nb::list pyCutsPerAtom;
    for (auto &cut : *cutsPerAtom) {
      nb::list localL;
      for (unsigned int j = 0; j < mol.getNumAtoms(); ++j) {
        localL.append(cut[j]);
      }
      pyCutsPerAtom.append(nb::tuple(localL));
    }
    delete cutsPerAtom;
    nb::list tres;
    tres.append(nb::tuple(res));
    tres.append(nb::tuple(pyCutsPerAtom));
    return nb::tuple(tres);
  } else {
    return nb::tuple(res);
  }
}

nb::tuple getShortestPathHelper(const ROMol &mol, int aid1, int aid2) {
  if (aid1 < 0 || aid1 >= rdcast<int>(mol.getNumAtoms()) || aid2 < 0 ||
      aid2 >= rdcast<int>(mol.getNumAtoms())) {
    throw ValueErrorException("bad atom index");
  }
  nb::list res;
  for (const auto atomIdx : MolOps::getShortestPath(mol, aid1, aid2)) {
    res.append(atomIdx);
  }
  return nb::steal<nb::tuple>(PySequence_Tuple(res.ptr()));
}

ROMol *fragmentOnBondsHelper(const ROMol &mol, nb::object pyBondIndices,
                             bool addDummies, nb::object pyDummyLabels,
                             nb::object pyBondTypes, nb::list pyCutsPerAtom) {
  auto bondIndices = pythonObjectToVect(pyBondIndices, mol.getNumBonds());
  if (!bondIndices.get()) {
    throw ValueErrorException("empty bond indices");
  }
  std::vector<std::pair<unsigned int, unsigned int>> *dummyLabels = nullptr;
  if (pyDummyLabels) {
    unsigned int nVs = nb::len(pyDummyLabels);
    dummyLabels = new std::vector<std::pair<unsigned int, unsigned int>>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      unsigned int v1 = nb::cast<unsigned int>(pyDummyLabels[i][0]);
      unsigned int v2 = nb::cast<unsigned int>(pyDummyLabels[i][1]);
      (*dummyLabels)[i] = std::make_pair(v1, v2);
    }
  }
  std::vector<Bond::BondType> *bondTypes = nullptr;
  if (pyBondTypes) {
    unsigned int nVs = nb::len(pyBondTypes);
    if (nVs != bondIndices->size()) {
      throw ValueErrorException("bondTypes shorter than bondIndices");
    }
    bondTypes = new std::vector<Bond::BondType>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      (*bondTypes)[i] = nb::cast<Bond::BondType>(pyBondTypes[i]);
    }
  }
  std::vector<unsigned int> *cutsPerAtom = nullptr;
  if (pyCutsPerAtom) {
    cutsPerAtom = new std::vector<unsigned int>;
    unsigned int nAts = nb::len(pyCutsPerAtom);
    if (nAts < mol.getNumAtoms()) {
      throw ValueErrorException("cutsPerAtom shorter than the number of atoms");
    }
    cutsPerAtom->resize(nAts);
  }

  ROMol *res = MolFragmenter::fragmentOnBonds(
      mol, *bondIndices, addDummies, dummyLabels, bondTypes, cutsPerAtom);
  if (cutsPerAtom) {
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      pyCutsPerAtom[i] = (*cutsPerAtom)[i];
    }
    delete cutsPerAtom;
  }

  delete dummyLabels;
  delete bondTypes;
  return res;
}

ROMol *renumberAtomsHelper(const ROMol &mol, nb::object &pyNewOrder) {
  if (nb::len(pyNewOrder) < mol.getNumAtoms()) {
    throw ValueErrorException("atomCounts shorter than the number of atoms");
  }
  auto newOrder = pythonObjectToVect(pyNewOrder, mol.getNumAtoms());
  if (!newOrder) {
    throw ValueErrorException("newOrder argument must be non-empty");
  }
  ROMol *res = MolOps::renumberAtoms(mol, *newOrder);
  return res;
}

namespace {
std::string getResidue(const ROMol &, const Atom *at) {
  auto monomerInfo = at->getMonomerInfo();
  if (!monomerInfo ||
      monomerInfo->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    return "";
  }
  return static_cast<const AtomPDBResidueInfo *>(monomerInfo)->getResidueName();
}
std::string getChainId(const ROMol &, const Atom *at) {
  auto monomerInfo = at->getMonomerInfo();
  if (!monomerInfo ||
      monomerInfo->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
    return "";
  }
  return static_cast<const AtomPDBResidueInfo *>(monomerInfo)->getChainId();
}
}  // namespace
nb::dict splitMolByPDBResidues(const ROMol &mol, nb::object pyWhiteList,
                               bool negateList) {
  std::vector<std::string> *whiteList = nullptr;
  if (pyWhiteList) {
    unsigned int nVs = nb::len(pyWhiteList);
    whiteList = new std::vector<std::string>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      (*whiteList)[i] = nb::cast<std::string>(pyWhiteList[i]);
    }
  }
  std::map<std::string, boost::shared_ptr<ROMol>> res =
      MolOps::getMolFragsWithQuery(mol, getResidue, false, whiteList,
                                   negateList);
  delete whiteList;

  nb::dict pyres;
  for (std::map<std::string, boost::shared_ptr<ROMol>>::const_iterator iter =
           res.begin();
       iter != res.end(); ++iter) {
    pyres[iter->first.c_str()] = iter->second;
  }
  return pyres;
}
nb::dict splitMolByPDBChainId(const ROMol &mol, nb::object pyWhiteList,
                              bool negateList) {
  std::vector<std::string> *whiteList = nullptr;
  if (pyWhiteList) {
    unsigned int nVs = nb::len(pyWhiteList);
    whiteList = new std::vector<std::string>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      (*whiteList)[i] = nb::cast<std::string>(pyWhiteList[i]);
    }
  }
  std::map<std::string, boost::shared_ptr<ROMol>> res =
      MolOps::getMolFragsWithQuery(mol, getChainId, false, whiteList,
                                   negateList);
  delete whiteList;

  nb::dict pyres;
  for (std::map<std::string, boost::shared_ptr<ROMol>>::const_iterator iter =
           res.begin();
       iter != res.end(); ++iter) {
    pyres[iter->first.c_str()] = iter->second;
  }
  return pyres;
}

nb::dict parseQueryDefFileHelper(nb::object &input, bool standardize,
                                 std::string delimiter, std::string comment,
                                 unsigned int nameColumn,
                                 unsigned int smartsColumn) {
  std::string input_text;
  std::map<std::string, ROMOL_SPTR> queryDefs;

  if (nb::isinstance<nb::str>(input)) {
    parseQueryDefFile(nb::cast<std::string>(input), queryDefs, standardize,
                      delimiter, comment, nameColumn, smartsColumn);
  } else {
    input_text = nb::cast<std::string>(input.attr("read")());
    std::istringstream istr(input_text);
    parseQueryDefFile(&istr, queryDefs, standardize, delimiter, comment,
                      nameColumn, smartsColumn);
  }

  nb::dict res;
  for (std::map<std::string, ROMOL_SPTR>::const_iterator iter =
           queryDefs.begin();
       iter != queryDefs.end(); ++iter) {
    res[iter->first.c_str()] = iter->second;
  }

  return res;
}

void addRecursiveQueriesHelper(ROMol &mol, nb::dict replDict,
                               std::string propName) {
  std::map<std::string, ROMOL_SPTR> replacements;
  const auto items = replDict.items();
  for (unsigned int i = 0; i < nb::len(items); ++i) {
    const auto item = items[i];
    ROMol *m = nb::cast<ROMol *>(item[1]);
    ROMOL_SPTR nm(new ROMol(*m));
    std::string k = nb::cast<std::string>(item[0]);
    replacements[k] = nm;
  }
  addRecursiveQueries(mol, replacements, propName);
}

ROMol *addHs2(const ROMol &orig, MolOps::AddHsParameters params,
              nb::object onlyOnAtoms) {
  std::unique_ptr<std::vector<unsigned int>> onlyOn;
  if (!onlyOnAtoms.is_none()) {
    onlyOn = pythonObjectToVect(onlyOnAtoms, orig.getNumAtoms());
  }
  auto res = std::make_unique<RWMol>(orig);
  MolOps::addHs(*res, params, onlyOn.get());
  return static_cast<ROMol *>(res.release());
}

ROMol *addHs(const ROMol &orig, bool explicitOnly, bool addCoords,
             nb::object onlyOnAtoms, bool addResidueInfo) {
  MolOps::AddHsParameters params{explicitOnly, addCoords, addResidueInfo};
  return addHs2(orig, params, onlyOnAtoms);
}

VECT_INT_VECT getSSSR(ROMol &mol, bool includeDativeBonds,
                      bool includeHydrogenBonds) {
  VECT_INT_VECT rings;
  MolOps::findSSSR(mol, rings, includeDativeBonds, includeHydrogenBonds);
  return rings;
}

nb::tuple replaceSubstructures(const ROMol &orig, const ROMol &query,
                               const ROMol &replacement,
                               bool replaceAll = false,
                               unsigned int replacementConnectionPoint = 0,
                               bool useChirality = false) {
  std::vector<ROMOL_SPTR> v =
      replaceSubstructs(orig, query, replacement, replaceAll,
                        replacementConnectionPoint, useChirality);
  nb::list res;
  for (const auto &mv : v) {
    res.append(mv);
  }
  return nb::steal<nb::tuple>(PySequence_Tuple(res.ptr()));
}

std::vector<MatchVectType> seqOfSeqsToMatchVectTypeVect(
    const nb::object &matches) {
  if (nb::len(matches) == 0) {
    throw ValueErrorException("matches must not be empty");
  }
  std::vector<MatchVectType> matchVectVect;
  for (size_t matchNum = 0; matchNum < nb::len(matches); ++matchNum) {
    auto matchObj = nb::cast<nb::object>(matches[matchNum]);
    auto match = pythonObjectToVect<unsigned int>(matchObj);
    if (!match) {
      throw ValueErrorException("tuples in matches must not be empty");
    }
    MatchVectType matchVect(match->size());
    for (unsigned int i = 0; i < match->size(); ++i) {
      matchVect[i] = std::make_pair(static_cast<int>(i), match->at(i));
    }
    matchVectVect.push_back(std::move(matchVect));
  }
  return matchVectVect;
}

std::vector<int> getMostSubstitutedCoreMatchHelper(const ROMol &mol,
                                                   const ROMol &core,
                                                   const nb::object &matches) {
  auto matchVectVect = seqOfSeqsToMatchVectTypeVect(matches);
  const auto &best = getMostSubstitutedCoreMatch(mol, core, matchVectVect);
  std::vector<int> res(best.size());
  for (const auto &pair : best) {
    res[pair.first] = pair.second;
  }
  return res;
}

std::vector<std::vector<int>> sortMatchesByDegreeOfCoreSubstitutionHelper(
    const ROMol &mol, const ROMol &core, const nb::object &matches) {
  auto matchVectVect = seqOfSeqsToMatchVectTypeVect(matches);
  auto sortedMatches =
      sortMatchesByDegreeOfCoreSubstitution(mol, core, matchVectVect);
  std::vector<std::vector<int>> res;
  res.reserve(sortedMatches.size());
  for (const auto &mv : sortedMatches) {
    std::vector<int> local(mv.size());
    for (const auto &pair : mv) {
      local[pair.first] = pair.second;
    }
    res.push_back(std::move(local));
  }
  return res;
}

void addRecursiveQuery(ROMol &mol, const ROMol &query, unsigned int atomIdx,
                       bool preserveExistingQuery) {
  if (atomIdx >= mol.getNumAtoms()) {
    throw ValueErrorException("atom index exceeds mol.GetNumAtoms()");
  }
  auto *q = new RecursiveStructureQuery(new ROMol(query));

  Atom *oAt = mol.getAtomWithIdx(atomIdx);
  if (!oAt->hasQuery()) {
    QueryAtom qAt(*oAt);
    static_cast<RWMol &>(mol).replaceAtom(atomIdx, &qAt);
    oAt = mol.getAtomWithIdx(atomIdx);
  }

  if (!preserveExistingQuery) {
    oAt->setQuery(q);
  } else {
    oAt->expandQuery(q, Queries::COMPOSITE_AND);
  }
}

void reapplyWedging(ROMol &mol, bool allBondTypes, bool verify) {
  auto &wmol = static_cast<RWMol &>(mol);
  RDKit::Chirality::reapplyMolBlockWedging(wmol, allBondTypes, verify);
}

MolOps::SanitizeFlags sanitizeMol(ROMol &mol, boost::uint64_t sanitizeOps,
                                  bool catchErrors) {
  auto &wmol = static_cast<RWMol &>(mol);
  unsigned int operationThatFailed;
  if (catchErrors) {
    try {
      MolOps::sanitizeMol(wmol, operationThatFailed, sanitizeOps);
    } catch (const MolSanitizeException &) {
      // this really should not be necessary, but at some point it
      // started to be required with VC++17. Doesn't seem like it does
      // any harm.
    } catch (...) {
    }
  } else {
    MolOps::sanitizeMol(wmol, operationThatFailed, sanitizeOps);
  }
  return static_cast<MolOps::SanitizeFlags>(operationThatFailed);
}

RWMol *getEditable(const ROMol &mol) {
  auto *res = new RWMol(mol, false);
  return res;
}

ROMol *getNormal(const RWMol &mol) {
  auto *res = static_cast<ROMol *>(new RWMol(mol));
  return res;
}

void kekulizeMol(ROMol &mol, bool clearAromaticFlags = false,
                 bool canonical = true) {
  auto &wmol = static_cast<RWMol &>(mol);
  MolOps::Kekulize(wmol, clearAromaticFlags, canonical);
}
void kekulizeMolIfPossible(ROMol &mol, bool clearAromaticFlags = false,
                           bool canonical = true) {
  auto &wmol = static_cast<RWMol &>(mol);
  MolOps::KekulizeIfPossible(wmol, clearAromaticFlags, canonical);
}

void cleanupMol(ROMol &mol) {
  auto &rwmol = static_cast<RWMol &>(mol);
  MolOps::cleanUp(rwmol);
}

void cleanUpOrganometallicsMol(ROMol &mol) {
  auto &rwmol = static_cast<RWMol &>(mol);
  MolOps::cleanUpOrganometallics(rwmol);
}

void setAromaticityMol(ROMol &mol, MolOps::AromaticityModel model) {
  auto &wmol = static_cast<RWMol &>(mol);
  MolOps::setAromaticity(wmol, model);
}

void setConjugationMol(ROMol &mol) {
  auto &wmol = static_cast<RWMol &>(mol);
  MolOps::setConjugation(wmol);
}

void assignRadicalsMol(ROMol &mol) {
  auto &wmol = static_cast<RWMol &>(mol);
  MolOps::assignRadicals(wmol);
}
namespace {
ROMol *hapticBondsToDativeHelper(const ROMol &mol) {
  ROMol *res = MolOps::hapticBondsToDative(mol);
  return res;
}
ROMol *dativeBondsToHapticHelper(const ROMol &mol) {
  ROMol *res = MolOps::dativeBondsToHaptic(mol);
  return res;
}
}  // namespace

void setHybridizationMol(ROMol &mol) {
  auto &wmol = static_cast<RWMol &>(mol);
  MolOps::setHybridization(wmol);
}

void cleanupChiralityMol(ROMol &mol) {
  auto &rwmol = static_cast<RWMol &>(mol);
  MolOps::cleanupChirality(rwmol);
}

void cleanupAtropisomersMol(ROMol &mol) {
  auto &rwmol = static_cast<RWMol &>(mol);
  MolOps::cleanupAtropisomers(rwmol);
}

VECT_INT_VECT getSymmSSSR(ROMol &mol, bool includeDativeBonds,
                          bool includeHydrogenBonds) {
  VECT_INT_VECT rings;
  MolOps::symmetrizeSSSR(mol, rings, includeDativeBonds, includeHydrogenBonds);
  return rings;
}
auto getDistanceMatrix(ROMol &mol, bool useBO = false, bool useAtomWts = false,
                       bool force = false, const char *prefix = nullptr) {
  const auto nats = mol.getNumAtoms();
  double *distMat =
      MolOps::getDistanceMat(mol, useBO, useAtomWts, force, prefix);
  double *dmCopy = new double[nats * nats];
  memcpy(static_cast<void *>(dmCopy), static_cast<void *>(distMat),
         nats * nats * sizeof(double));
  // pattern from
  // https://nanobind.readthedocs.io/en/latest/ndarray.html#data-ownership
  nb::capsule owner(dmCopy, [](void *f) noexcept {
    double *data = reinterpret_cast<double *>(f);
    delete[] data;
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<2>>(
      /* data = */ dmCopy,
      /* shape = */ {nats, nats},
      /* owner = */ owner);
}
auto get3DDistanceMatrix(ROMol &mol, int confId = -1, bool useAtomWts = false,
                         bool force = false, const char *prefix = nullptr) {
  const auto nats = mol.getNumAtoms();
  double *distMat =
      MolOps::get3DDistanceMat(mol, confId, useAtomWts, force, prefix);

  double *dmCopy = new double[nats * nats];
  memcpy(static_cast<void *>(dmCopy), static_cast<void *>(distMat),
         nats * nats * sizeof(double));
  if (prefix == nullptr || std::string(prefix) == "") {
    delete[] distMat;
  }
  // pattern from
  // https://nanobind.readthedocs.io/en/latest/ndarray.html#data-ownership
  nb::capsule owner(dmCopy, [](void *f) noexcept {
    double *data = reinterpret_cast<double *>(f);
    delete[] data;
  });
  return nb::ndarray<nb::numpy, double, nb::ndim<2>>(
      /* data = */ dmCopy,
      /* shape = */ {nats, nats},
      /* owner = */ owner);
}

nb::object getAdjacencyMatrix(ROMol &mol, bool useBO = false, int emptyVal = 0,
                              bool force = false,
                              const char *prefix = nullptr) {
  const auto nats = mol.getNumAtoms();

  double *tmpMat =
      MolOps::getAdjacencyMatrix(mol, useBO, emptyVal, force, prefix);

  if (useBO) {
    double *resMat = new double[nats * nats];
    memcpy(static_cast<void *>(resMat), static_cast<void *>(tmpMat),
           nats * nats * sizeof(double));
    // pattern from
    // https://nanobind.readthedocs.io/en/latest/ndarray.html#data-ownership
    nb::capsule owner(resMat, [](void *f) noexcept {
      double *data = reinterpret_cast<double *>(f);
      delete[] data;
    });
    auto res = nb::ndarray<nb::numpy, double, nb::ndim<2>>(
        /* data = */ resMat,
        /* shape = */ {nats, nats},
        /* owner = */ owner);
    return res.cast();
  } else {
    int *resMat = new int[nats * nats];
    for (unsigned int i = 0; i < nats; i++) {
      for (unsigned int j = 0; j < nats; j++) {
        resMat[i * nats + j] = (int)std::round(tmpMat[i * nats + j]);
      }
    }
    // pattern from
    // https://nanobind.readthedocs.io/en/latest/ndarray.html#data-ownership
    nb::capsule owner(resMat, [](void *f) noexcept {
      int *data = reinterpret_cast<int *>(f);
      delete[] data;
    });
    auto res = nb::ndarray<nb::numpy, int, nb::ndim<2>>(
        /* data = */ resMat,
        /* shape = */ {nats, nats},
        /* owner = */ owner);
    return res.cast();
  }
}

nb::tuple GetMolFragsWithMapping(const ROMol &mol, bool asMols,
                                 bool sanitizeFrags,
                                 nb::object frags = nb::none(),
                                 nb::object fragsMolAtomMapping = nb::none()) {
  nb::list res;

  if (!asMols) {
    VECT_INT_VECT fragsVec;
    MolOps::getMolFrags(mol, fragsVec);

    for (auto &i : fragsVec) {
      nb::list tpl;
      for (unsigned int j = 0; j < i.size(); ++j) {
        tpl.append(i[j]);
      }
      res.append(nb::tuple(tpl));
    }
  } else {
    std::vector<std::vector<int>> fragsMolAtomMappingVec;
    std::vector<int> fragsVec;
    std::vector<boost::shared_ptr<ROMol>> molFrags;
    auto &fragsList = reinterpret_cast<nb::list &>(frags);
    auto &fragsMolAtomMappingList =
        reinterpret_cast<nb::list &>(fragsMolAtomMapping);
    bool hasFrags = !fragsList.is_none();
    bool hasFragsMolAtomMapping = !fragsMolAtomMappingList.is_none();
    molFrags =
        hasFrags || hasFragsMolAtomMapping
            ? MolOps::getMolFrags(
                  mol, sanitizeFrags, hasFrags ? &fragsVec : nullptr,
                  hasFragsMolAtomMapping ? &fragsMolAtomMappingVec : nullptr)
            : MolOps::getMolFrags(mol, sanitizeFrags);
    if (hasFrags) {
      for (int i : fragsVec) {
        fragsList.append(i);
      }
    }
    if (hasFragsMolAtomMapping) {
      for (auto &i : fragsMolAtomMappingVec) {
        nb::list perFragMolAtomMappingTpl;
        for (unsigned int j = 0; j < i.size(); ++j) {
          perFragMolAtomMappingTpl.append(i[j]);
        }
        fragsMolAtomMappingList.append(nb::tuple(perFragMolAtomMappingTpl));
      }
    }
    for (const auto &molFrag : molFrags) {
      res.append(molFrag);
    }
  }
  return nb::tuple(res);
}

nb::tuple GetMolFrags(const ROMol &mol, bool asMols, bool sanitizeFrags) {
  return GetMolFragsWithMapping(mol, asMols, sanitizeFrags);
}

ExplicitBitVect *wrapLayeredFingerprint(
    const ROMol &mol, unsigned int layerFlags, unsigned int minPath,
    unsigned int maxPath, unsigned int fpSize, nb::list atomCounts,
    ExplicitBitVect *includeOnlyBits, bool branchedPaths,
    nb::object fromAtoms) {
  std::unique_ptr<std::vector<unsigned int>> lFromAtoms =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::unique_ptr<std::vector<unsigned int>> atomCountsV;
  if (atomCounts) {
    atomCountsV.reset(new std::vector<unsigned int>);
    unsigned int nAts = nb::len(atomCounts);
    if (nAts < mol.getNumAtoms()) {
      throw ValueErrorException("atomCounts shorter than the number of atoms");
    }
    atomCountsV->resize(nAts);
    for (unsigned int i = 0; i < nAts; ++i) {
      (*atomCountsV)[i] = nb::cast<unsigned int>(atomCounts[i]);
    }
  }

  ExplicitBitVect *res;
  res = RDKit::LayeredFingerprintMol(mol, layerFlags, minPath, maxPath, fpSize,
                                     atomCountsV.get(), includeOnlyBits,
                                     branchedPaths, lFromAtoms.get());

  if (atomCountsV) {
    for (unsigned int i = 0; i < atomCountsV->size(); ++i) {
      atomCounts[i] = (*atomCountsV)[i];
    }
  }

  return res;
}
ExplicitBitVect *wrapPatternFingerprint(const ROMol &mol, unsigned int fpSize,
                                        nb::list atomCounts,
                                        ExplicitBitVect *includeOnlyBits,
                                        bool tautomerFingerprints) {
  std::vector<unsigned int> *atomCountsV = nullptr;
  if (atomCounts) {
    atomCountsV = new std::vector<unsigned int>;
    unsigned int nAts = nb::len(atomCounts);
    if (nAts < mol.getNumAtoms()) {
      throw ValueErrorException("atomCounts shorter than the number of atoms");
    }
    atomCountsV->resize(nAts);
    for (unsigned int i = 0; i < nAts; ++i) {
      (*atomCountsV)[i] = nb::cast<unsigned int>(atomCounts[i]);
    }
  }

  ExplicitBitVect *res;
  res = RDKit::PatternFingerprintMol(mol, fpSize, atomCountsV, includeOnlyBits,
                                     tautomerFingerprints);

  if (atomCountsV) {
    for (unsigned int i = 0; i < atomCountsV->size(); ++i) {
      atomCounts[i] = (*atomCountsV)[i];
    }
    delete atomCountsV;
  }

  return res;
}
ExplicitBitVect *wrapPatternFingerprintBundle(const MolBundle &bundle,
                                              unsigned int fpSize,
                                              ExplicitBitVect *includeOnlyBits,
                                              bool tautomerFingerprints) {
  ExplicitBitVect *res;
  res = RDKit::PatternFingerprintMol(bundle, fpSize, includeOnlyBits,
                                     tautomerFingerprints);
  return res;
}

ExplicitBitVect *wrapRDKFingerprintMol(
    const ROMol &mol, unsigned int minPath, unsigned int maxPath,
    unsigned int fpSize, unsigned int nBitsPerHash, bool useHs,
    double tgtDensity, unsigned int minSize, bool branchedPaths,
    bool useBondOrder, nb::object atomInvariants, nb::object fromAtoms,
    nb::object atomBits, nb::object bitInfo) {
  std::unique_ptr<std::vector<unsigned int>> lAtomInvariants =
      pythonObjectToVect<unsigned int>(atomInvariants);
  std::unique_ptr<std::vector<unsigned int>> lFromAtoms =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::vector<std::vector<std::uint32_t>> *lAtomBits = nullptr;
  std::map<std::uint32_t, std::vector<std::vector<int>>> *lBitInfo = nullptr;
  // if(!(atomBits.is_none())){
  if (!atomBits.is_none()) {
    lAtomBits = new std::vector<std::vector<std::uint32_t>>(mol.getNumAtoms());
  }
  if (!bitInfo.is_none()) {
    lBitInfo = new std::map<std::uint32_t, std::vector<std::vector<int>>>;
  }
  ExplicitBitVect *res;
  res = RDKit::RDKFingerprintMol(mol, minPath, maxPath, fpSize, nBitsPerHash,
                                 useHs, tgtDensity, minSize, branchedPaths,
                                 useBondOrder, lAtomInvariants.get(),
                                 lFromAtoms.get(), lAtomBits, lBitInfo);

  if (lAtomBits) {
    auto &pyl = static_cast<nb::list &>(atomBits);
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      nb::list tmp;
      for (auto v : (*lAtomBits)[i]) {
        tmp.append(v);
      }
      pyl.append(tmp);
    }
    delete lAtomBits;
  }
  if (lBitInfo) {
    auto &pyd = static_cast<nb::dict &>(bitInfo);
    for (auto &it : (*lBitInfo)) {
      nb::list temp;
      std::vector<std::vector<int>>::iterator itset;
      for (itset = it.second.begin(); itset != it.second.end(); ++itset) {
        nb::list temp2;
        for (int &i : *itset) {
          temp2.append(i);
        }
        temp.append(temp2);
      }
      pyd[it.first] = temp;
    }
    delete lBitInfo;
  }

  return res;
}

SparseIntVect<boost::uint64_t> *wrapUnfoldedRDKFingerprintMol(
    const ROMol &mol, unsigned int minPath, unsigned int maxPath, bool useHs,
    bool branchedPaths, bool useBondOrder, nb::object atomInvariants,
    nb::object fromAtoms, nb::object atomBits, nb::object bitInfo) {
  std::unique_ptr<std::vector<unsigned int>> lAtomInvariants =
      pythonObjectToVect<unsigned int>(atomInvariants);
  std::unique_ptr<std::vector<unsigned int>> lFromAtoms =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::vector<std::vector<boost::uint64_t>> *lAtomBits = nullptr;
  std::map<boost::uint64_t, std::vector<std::vector<int>>> *lBitInfo = nullptr;

  // if(!(atomBits.is_none())){
  if (!atomBits.is_none()) {
    lAtomBits =
        new std::vector<std::vector<boost::uint64_t>>(mol.getNumAtoms());
  }
  if (!bitInfo.is_none()) {
    lBitInfo = new std::map<boost::uint64_t, std::vector<std::vector<int>>>;
  }

  SparseIntVect<boost::uint64_t> *res;
  res = getUnfoldedRDKFingerprintMol(
      mol, minPath, maxPath, useHs, branchedPaths, useBondOrder,
      lAtomInvariants.get(), lFromAtoms.get(), lAtomBits, lBitInfo);

  if (lAtomBits) {
    auto &pyl = static_cast<nb::list &>(atomBits);
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      nb::list tmp;
      for (auto v : (*lAtomBits)[i]) {
        tmp.append(v);
      }
      pyl.append(tmp);
    }
    delete lAtomBits;
  }
  if (lBitInfo) {
    auto &pyd = static_cast<nb::dict &>(bitInfo);
    for (auto &it : (*lBitInfo)) {
      nb::list temp;
      std::vector<std::vector<int>>::iterator itset;
      for (itset = it.second.begin(); itset != it.second.end(); ++itset) {
        nb::list temp2;
        for (int &i : *itset) {
          temp2.append(i);
        }
        temp.append(temp2);
      }
      pyd[it.first] = temp;
    }
    delete lBitInfo;
  }

  return res;
}

nb::object findAllSubgraphsOfLengthsMtoNHelper(const ROMol &mol,
                                               unsigned int lowerLen,
                                               unsigned int upperLen,
                                               bool useHs = false,
                                               int rootedAtAtom = -1) {
  if (lowerLen > upperLen) {
    throw ValueErrorException("lowerLen > upperLen");
  }

  INT_PATH_LIST_MAP oMap = findAllSubgraphsOfLengthsMtoN(
      mol, lowerLen, upperLen, useHs, rootedAtAtom);
  nb::list res;
  for (unsigned int i = lowerLen; i <= upperLen; ++i) {
    nb::list tmp;
    const PATH_LIST &pth = oMap[i];
    for (const auto &pthit : pth) {
      tmp.append(nb::cast(pthit));
    }
    res.append(tmp);
  }
  return nb::steal<nb::tuple>(PySequence_Tuple(res.ptr()));
};

PATH_TYPE findAtomEnvironmentOfRadiusNHelper(const ROMol &mol,
                                             unsigned int radius,
                                             unsigned int rootedAtAtom,
                                             bool useHs, bool enforceSize,
                                             nb::object atomMap) {
  PATH_TYPE path;
  if (atomMap.is_none()) {
    path = findAtomEnvironmentOfRadiusN(mol, radius, rootedAtAtom, useHs,
                                        enforceSize);
  } else {
    std::unordered_map<unsigned int, unsigned int> cAtomMap;
    path = findAtomEnvironmentOfRadiusN(mol, radius, rootedAtAtom, useHs,
                                        enforceSize, &cAtomMap);
    // make sure the optional argument (atomMap) is actually a dictionary
    nb::dict typecheck = nb::cast<nb::dict>(atomMap);
    atomMap.attr("clear")();
    for (auto pair : cAtomMap) {
      atomMap[pair.first] = pair.second;
    }
  }
  return path;
}

ROMol *pathToSubmolHelper(const ROMol &mol, nb::object &path, bool useQuery,
                          nb::object atomMap) {
  ROMol *result;
  PATH_TYPE pth;
  for (unsigned int i = 0; i < nb::len(path); ++i) {
    pth.push_back(nb::cast<unsigned int>(path[i]));
  }
  std::map<int, int> mapping;
  result = Subgraphs::pathToSubmol(mol, pth, useQuery, mapping);
  if (!atomMap.is_none()) {
    // make sure the optional argument actually was a dictionary
    nb::dict typecheck = nb::cast<nb::dict>(atomMap);
    atomMap.attr("clear")();
    for (std::map<int, int>::const_iterator mIt = mapping.begin();
         mIt != mapping.end(); ++mIt) {
      atomMap[mIt->first] = mIt->second;
    }
  }
  return result;
}

ROMol *adjustQueryPropertiesHelper(const ROMol &mol, nb::object pyparams) {
  MolOps::AdjustQueryParameters params;
  if (!pyparams.is_none()) {
    params = nb::cast<MolOps::AdjustQueryParameters>(pyparams);
  }
  return MolOps::adjustQueryProperties(mol, &params);
}

ROMol *adjustQueryPropertiesWithGenericGroupsHelper(const ROMol &mol,
                                                    nb::object pyparams) {
  MolOps::AdjustQueryParameters params;
  if (!pyparams.is_none()) {
    params = nb::cast<MolOps::AdjustQueryParameters>(pyparams);
  }
  return GenericGroups::adjustQueryPropertiesWithGenericGroups(mol, &params);
}

nb::tuple detectChemistryProblemsHelper(const ROMol &mol,
                                        unsigned int sanitizeOps) {
  auto probs = MolOps::detectChemistryProblems(mol, sanitizeOps);
  nb::list res;
  for (const auto &exc_ptr : probs) {
    res.append(boost::shared_ptr<MolSanitizeException>(exc_ptr->copy()));
  }
  return nb::tuple(res);
}

ROMol *canonicalizeStereoGroupsHelper(
    ROMol &mol, RDKit::StereoGroupAbsOptions stereoGroupAbsOptions,
    unsigned int maxStereoGroups) {
  auto mol_uptr = std::unique_ptr<ROMol>(new ROMol(mol));

  RDKit::canonicalizeStereoGroups(mol_uptr, stereoGroupAbsOptions,
                                  maxStereoGroups);
  return mol_uptr.release();
  ;
}

ROMol *replaceCoreHelper(const ROMol &mol, const ROMol &core, nb::object match,
                         bool replaceDummies, bool labelByIndex,
                         bool requireDummyMatch = false) {
  // convert input to MatchVect
  MatchVectType matchVect;

  unsigned int length = nb::len(match);
  for (unsigned int i = 0; i < length; ++i) {
    // This is what boost::nb::len() does internally
    auto pyObj = static_cast<nb::object>(match[i]).ptr();
    unsigned int sz = PyObject_Length(pyObj);
    if (PyErr_Occurred()) {
      PyErr_Clear();
      sz = 1;
    }

    int v1, v2;
    switch (sz) {
      case 1:
        if (length != core.getNumAtoms()) {
          std::string entries = core.getNumAtoms() == 1 ? " entry" : " entries";

          std::stringstream ss;
          ss << std::string(
                    "When using input vector of type (molecule_atom_idx,...) "
                    "supplied core requires ")
             << core.getNumAtoms() << entries;
          throw ValueErrorException(ss.str());
        }
        v1 = (int)i;
        v2 = nb::cast<int>(match[i]);
        break;
      case 2:
        v1 = nb::cast<int>(match[i][0]);
        v2 = nb::cast<int>(match[i][1]);
        break;
      default:
        throw ValueErrorException(
            "Input not a vector of (core_atom_idx,molecule_atom_idx) or "
            "(molecule_atom_idx,...) entries");
    }
    matchVect.push_back(std::make_pair(v1, v2));
  }

  return replaceCore(mol, core, matchVect, replaceDummies, labelByIndex,
                     requireDummyMatch);
}

void setDoubleBondNeighborDirectionsHelper(ROMol &mol, nb::object confObj) {
  Conformer *conf = nullptr;
  if (confObj) {
    conf = nb::cast<Conformer *>(confObj);
  }
  MolOps::setDoubleBondNeighborDirections(mol, conf);
}

void setAtomSymbols(MolzipParams &p, nb::object symbols) {
  p.atomSymbols.clear();
  if (symbols) {
    unsigned int nVs = nb::len(symbols);
    for (unsigned int i = 0; i < nVs; ++i) {
      p.atomSymbols.push_back(nb::cast<std::string>(symbols[i]));
    }
  }
}

ROMol *molzip_new(const ROMol &a, const ROMol &b, const MolzipParams &p) {
  return molzip(a, b, p).release();
}

ROMol *molzip_new(const ROMol &a, const MolzipParams &p) {
  return molzip(a, p).release();
}

ROMol *molzipHelper(nb::object &pmols, const MolzipParams &p) {
  auto mols = pythonObjectToVect<ROMOL_SPTR>(pmols);
  if (mols == nullptr || mols->empty()) {
    return nullptr;
  }
  return molzip(*mols, p).release();
}

ROMol *rgroupRowZipHelper(nb::dict row, const MolzipParams &p) {
  std::map<std::string, ROMOL_SPTR> rgroup_row;
  nb::list items = row.items();
  for (size_t i = 0; i < (size_t)nb::len(items); ++i) {
    nb::object key = items[i][0];
    nb::object value = items[i][1];
    try {
      auto rgroup_key = nb::cast<std::string>(key);
      auto mol = nb::cast<ROMOL_SPTR>(value);
      rgroup_row[rgroup_key] = mol;
    } catch (...) {
      throw ValueErrorException(
          "Unable to retrieve rgroup key and molecule from dictionary");
    }
  }

  return molzip(rgroup_row, p).release();
}

nb::tuple hasQueryHsHelper(const ROMol &m) {
  nb::list res;
  auto hashs = MolOps::hasQueryHs(m);
  res.append(hashs.first);
  res.append(hashs.second);
  return nb::tuple(res);
}

// we can really only set some of these types from C++ which means
//  we need a helper function for testing that we can read them
//  correctly.
void _testSetProps(RDProps &props, const std::string &prefix) {
  props.setProp<bool>(prefix + "bool", true);
  props.setProp<unsigned int>(prefix + "uint", -1);
  props.setProp<double>(prefix + "double", 3.14159);

  std::vector<int> svint;
  svint.push_back(0);
  svint.push_back(1);
  svint.push_back(2);
  svint.push_back(-2);

  props.setProp<std::vector<int>>(prefix + "svint", svint);

  std::vector<unsigned int> svuint;
  svuint.push_back(0);
  svuint.push_back(1);
  svuint.push_back(2);
  svuint.push_back(-2);

  props.setProp<std::vector<unsigned int>>(prefix + "svuint", svuint);

  std::vector<double> svdouble;
  svdouble.push_back(0.);
  svdouble.push_back(1.);
  svdouble.push_back(2.);
  props.setProp<std::vector<double>>(prefix + "svdouble", svdouble);

  std::vector<std::string> svstring;
  svstring.push_back("The");
  svstring.push_back("RDKit");

  props.setProp<std::vector<std::string>>(prefix + "svstring", svstring);
}

void testSetProps(ROMol &mol) {
  _testSetProps(mol, "mol_");
  for (auto &atom : mol.atoms()) {
    _testSetProps(*atom, std::string("atom_") + std::to_string(atom->getIdx()));
  }
  for (auto &bond : mol.bonds()) {
    _testSetProps(*bond, std::string("bond_") + std::to_string(bond->getIdx()));
  }
  for (unsigned conf_idx = 0; conf_idx < mol.getNumConformers(); ++conf_idx) {
    _testSetProps(mol.getConformer(conf_idx),
                  "conf_" + std::to_string(conf_idx));
  }
}

void expandAttachmentPointsHelper(ROMol &mol, bool addAsQueries,
                                  bool addCoords) {
  MolOps::expandAttachmentPoints(static_cast<RWMol &>(mol), addAsQueries,
                                 addCoords);
}

void collapseAttachmentPointsHelper(ROMol &mol, bool markedOnly) {
  MolOps::collapseAttachmentPoints(static_cast<RWMol &>(mol), markedOnly);
}

nb::object findMesoHelper(const ROMol &mol, bool includeIsotopes,
                          bool includeAtomMaps) {
  auto meso = Chirality::findMesoCenters(mol, includeIsotopes, includeAtomMaps);
  nb::list res;
  for (const auto &pr : meso) {
    nb::list tpl;
    tpl.append(pr.first);
    tpl.append(pr.second);
    res.append(nb::tuple(tpl));
  }
  return nb::tuple(res);
}

ROMol *copyMolSubsetHelper1(const ROMol &mol, nb::object pyAtomIndices,
                            nb::object pyBondIndices,
                            const SubsetOptions &options = SubsetOptions()) {
  auto atomIndices = pythonObjectToVect<unsigned int>(pyAtomIndices);
  auto bondIndices = pythonObjectToVect<unsigned int>(pyBondIndices);
  if (!atomIndices.get()) {
    atomIndices = std::make_unique<std::vector<unsigned int>>();
  }
  if (!bondIndices.get()) {
    bondIndices = std::make_unique<std::vector<unsigned int>>();
  }

  return copyMolSubset(mol, *atomIndices, *bondIndices, options).release();
}

ROMol *copyMolSubsetHelper2(const ROMol &mol, nb::object pyAtomIndices,
                            nb::object pyBondIndices, SubsetInfo &info,
                            const SubsetOptions &options = SubsetOptions()) {
  auto atomIndices = pythonObjectToVect<unsigned int>(pyAtomIndices);
  auto bondIndices = pythonObjectToVect<unsigned int>(pyBondIndices);
  if (!atomIndices.get()) {
    atomIndices = std::make_unique<std::vector<unsigned int>>();
  }
  if (!bondIndices.get()) {
    bondIndices = std::make_unique<std::vector<unsigned int>>();
  }

  return copyMolSubset(mol, *atomIndices, *bondIndices, info, options)
      .release();
}

ROMol *copyMolSubsetHelper3(const ROMol &mol, nb::object path,
                            const SubsetOptions &options = SubsetOptions()) {
  auto pathvect = pythonObjectToVect<unsigned int>(path);
  if (!pathvect.get()) {
    pathvect = std::make_unique<std::vector<unsigned int>>();
  }
  return copyMolSubset(mol, *pathvect, options).release();
}

ROMol *copyMolSubsetHelper4(const ROMol &mol, nb::object path,
                            SubsetInfo &selectionInfo,
                            const SubsetOptions &options = SubsetOptions()) {
  auto pathvect = pythonObjectToVect<unsigned int>(path);
  if (!pathvect.get()) {
    pathvect = std::make_unique<std::vector<unsigned int>>();
  }
  return copyMolSubset(mol, *pathvect, selectionInfo, options).release();
}

struct molops_wrapper {
  static void wrap(nb::module_ &m) {
    std::string docString;
    nb::enum_<MolOps::SanitizeFlags>(m, "SanitizeFlags")
        .value("SANITIZE_NONE", MolOps::SANITIZE_NONE)
        .value("SANITIZE_CLEANUP", MolOps::SANITIZE_CLEANUP)
        .value("SANITIZE_PROPERTIES", MolOps::SANITIZE_PROPERTIES)
        .value("SANITIZE_SYMMRINGS", MolOps::SANITIZE_SYMMRINGS)
        .value("SANITIZE_KEKULIZE", MolOps::SANITIZE_KEKULIZE)
        .value("SANITIZE_FINDRADICALS", MolOps::SANITIZE_FINDRADICALS)
        .value("SANITIZE_SETAROMATICITY", MolOps::SANITIZE_SETAROMATICITY)
        .value("SANITIZE_SETCONJUGATION", MolOps::SANITIZE_SETCONJUGATION)
        .value("SANITIZE_SETHYBRIDIZATION", MolOps::SANITIZE_SETHYBRIDIZATION)
        .value("SANITIZE_CLEANUPCHIRALITY", MolOps::SANITIZE_CLEANUPCHIRALITY)
        .value("SANITIZE_CLEANUPATROPISOMERS",
               MolOps::SANITIZE_CLEANUPATROPISOMERS)
        .value("SANITIZE_ADJUSTHS", MolOps::SANITIZE_ADJUSTHS)
        .value("SANITIZE_CLEANUP_ORGANOMETALLICS",
               MolOps::SANITIZE_CLEANUP_ORGANOMETALLICS)
        .value("SANITIZE_ALL", MolOps::SANITIZE_ALL)
        .export_values();
    ;

    // ------------------------------------------------------------------------
    docString =
        "Assign stereochemistry to bonds based on coordinates and a conformer.\n\
        DEPRECATED\n\
        \n\
  ARGUMENTS:\n\
  \n\
    - mol: the molecule to be modified\n\
    - conformer: Conformer providing the coordinates\n\
\n";
    m.def("DetectBondStereoChemistry", DetectBondStereoChemistry, "mol"_a,
          "conformer"_a, docString.c_str());
    docString =
        "DEPRECATED\n\
    - mol: the molecule to be modified\n\
    - confId: Conformer to use for the coordinates\n\
\n";
    m.def("DetectBondStereochemistry", MolOps::detectBondStereochemistry,
          "mol"_a, "confId"_a = -1, docString.c_str());

    docString =
        "Uses the stereo info on double bonds to set the directions of neighboring single bonds\n\
        \n\
  ARGUMENTS:\n\
  \n\
    - mol: the molecule to be modified\n\
\n";
    m.def("SetDoubleBondNeighborDirections",
          setDoubleBondNeighborDirectionsHelper, "mol"_a, "conf"_a = nb::none(),
          docString.c_str());

    docString =
        "Uses the directions of neighboring bonds to set cis/trans stereo on double bonds.\n\
        \n\
  ARGUMENTS:\n\
  \n\
    - mol: the molecule to be modified\n\
\n";
    m.def("SetBondStereoFromDirections", MolOps::setBondStereoFromDirections,
          "mol"_a, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Kekulize, check valencies, set aromaticity, conjugation and hybridization\n\
\n\
    - The molecule is modified in place.\n\
\n\
    - If sanitization fails, an exception will be thrown unless catchErrors is set\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
    - sanitizeOps: (optional) sanitization operations to be carried out\n\
      these should be constructed by or'ing together the\n\
      operations in rdkit.Chem.SanitizeFlags\n\
    - catchErrors: (optional) if provided, instead of raising an exception\n\
      when sanitization fails (the default behavior), the \n\
      first operation that failed (as defined in rdkit.Chem.SanitizeFlags)\n\
      is returned. Zero is returned on success.\n\
\n";
    m.def("SanitizeMol", sanitizeMol, "mol"_a,
          "sanitizeOps"_a = static_cast<boost::uint64_t>(MolOps::SANITIZE_ALL),
          "catchErrors"_a = false, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Get the smallest set of simple rings for a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use.\n\
    - includeDativeBonds: whether or not dative bonds should be included in the ring finding.\n\
    - includeHydrogenBonds: whether or not hydrogen bonds should be included in the ring finding.\n\
\n\
  RETURNS: a sequence of sequences containing the rings found as atom ids\n\
         The length of this will be equal to NumBonds-NumAtoms+1 for single-fragment molecules.\n\
\n";
    m.def("GetSSSR", getSSSR, "mol"_a, "includeDativeBonds"_a = false,
          "includeHydrogenBonds"_a = false, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Get a symmetrized SSSR for a molecule.\n\
\n\
  The symmetrized SSSR is at least as large as the SSSR for a molecule.\n\
  In certain highly-symmetric cases (e.g. cubane), the symmetrized SSSR can be\n\
  a bit larger (i.e. the number of symmetrized rings is >= NumBonds-NumAtoms+1).\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use.\n\
    - includeDativeBonds: whether or not dative bonds should be included in the ring finding.\n\
    - includeHydrogenBonds: whether or not hydrogen bonds should be included in the ring finding.\n\
\n\
  RETURNS: a sequence of sequences containing the rings found as atom ids\n\
\n";
    m.def("GetSymmSSSR", getSymmSSSR, "mol"_a, "includeDativeBonds"_a = false,
          "includeHydrogenBonds"_a = false, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Sets Cartesian coordinates for a terminal atom.\n\
\n\
  Useful for growing an atom off a molecule with sensible \n\
  coordinates based on the geometry of the neighbor.\n\
\n\
  NOTE: this sets the appropriate coordinates in all of the molecule's conformers \n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule the atoms belong to.\n\
    - idx: index of the terminal atom whose coordinates are set.\n\
    - mol: index of the bonded neighbor atom.\n\
\n\
  RETURNS: Nothing\n\
\n";
    m.def("SetTerminalAtomCoords", MolOps::setTerminalAtomCoords, "mol"_a,
          "idx"_a, "otherIdx"_a, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Does a non-SSSR ring finding for a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use.\n\
\n\
  RETURNS: Nothing\n\
\n";
    m.def("FastFindRings", MolOps::fastFindRings, docString.c_str(), "mol"_a);

    docString =
        "Generate Unique Ring Families.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use.\n\
    - includeDativeBonds: whether or not dative bonds should be included in the ring families finding.\n\
    - includeHydrogenBonds: whether or not hydrogen bonds should be included in the ring families .\n\
\n\
  RETURNS: Nothing\n\
\n";
    m.def("FindRingFamilies", MolOps::findRingFamilies, "mol"_a,
          "includeDativeBonds"_a = false, "includeHydrogenBonds"_a = false,
          docString.c_str());

    // ------------------------------------------------------------------------
    docString = R"DOC(Parameters controlling H addition.)DOC";
    nb::class_<MolOps::AddHsParameters>(m, "AddHsParameters", docString.c_str())
        .def(nb::init<>())
        .def_rw("explicitOnly", &MolOps::AddHsParameters::explicitOnly,
                "only add explict Hs")
        .def_rw("addCoords", &MolOps::AddHsParameters::addCoords,
                "add coordinates for the Hs")
        .def_rw("addResidueInfo", &MolOps::AddHsParameters::addResidueInfo,
                "add residue info to the Hs")
        .def_rw("skipQueries", &MolOps::AddHsParameters::skipQueries,
                "do not add Hs to query atoms or atoms with query bonds");

    // ------------------------------------------------------------------------
    docString =
        R"DOC(Adds hydrogens to the graph of a molecule.

  ARGUMENTS:

    - mol: the molecule to be modified

    - params: AddHsParameters object controlling the addition.

    - onlyOnAtoms: (optional) if this sequence is provided, only these atoms will be
      considered to have Hs added to them

  RETURNS: a new molecule with added Hs

  NOTES:

    - The original molecule is *not* modified.

    - Much of the code assumes that Hs are not included in the molecular
      topology, so be *very* careful with the molecule that comes back from
      this function.\n)DOC";
    m.def("AddHs", addHs2, "mol"_a, "params"_a, "onlyOnAtoms"_a = nb::none(),
          docString.c_str(), nb::rv_policy::take_ownership);

    // ------------------------------------------------------------------------
    docString =
        "Adds hydrogens to the graph of a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - explicitOnly: (optional) if this toggle is set, only explicit Hs will\n\
      be added to the molecule.  Default value is 0 (add implicit and explicit Hs).\n\
\n\
    - addCoords: (optional) if this toggle is set, The Hs will have 3D coordinates\n\
      set.  Default value is 0 (no 3D coords).\n\
\n\
    - onlyOnAtoms: (optional) if this sequence is provided, only these atoms will be\n\
      considered to have Hs added to them\n\
\n\
    - addResidueInfo: (optional) if this is true, add residue info to\n\
      hydrogen atoms (useful for PDB files).\n\
\n\
  RETURNS: a new molecule with added Hs\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
    - Much of the code assumes that Hs are not included in the molecular\n\
      topology, so be *very* careful with the molecule that comes back from\n\
      this function.\n\
\n";
    m.def("AddHs", addHs, "mol"_a, "explicitOnly"_a = false,
          "addCoords"_a = false, "onlyOnAtoms"_a = nb::none(),
          "addResidueInfo"_a = false, docString.c_str(),
          nb::rv_policy::take_ownership);

    // ------------------------------------------------------------------------
    docString =
        "Removes any hydrogens from the graph of a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - implicitOnly: (optional) if this toggle is set, only implicit Hs will\n\
      be removed from the graph.  Default value is 0 (remove implicit and explicit Hs).\n\
\n\
    - updateExplicitCount: (optional) if this toggle is set, the explicit H count on atoms with \n\
      Hs will be updated. Default value is 0 (do not update explicit H count).\n\
\n\
    - sanitize: (optional) if this toggle is set, the molecule will be sanitized after the Hs\n\
      are removed. Default value is 1 (do sanitize).\n\
\n\
  RETURNS: a new molecule with the Hs removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
    - Hydrogens which aren't connected to a heavy atom will not be\n\
      removed.  This prevents molecules like [H][H] from having\n\
      all atoms removed.\n\
    - Labelled hydrogen (e.g. atoms with atomic number=1, but isotope > 1),\n\
      will not be removed.\n\
    - two coordinate Hs, like the central H in C[H-]C, will not be removed\n\
    - Hs connected to dummy atoms will not be removed\n\
    - Hs that are part of the definition of double bond Stereochemistry\n\
      will not be removed\n\
    - Hs that are not connected to anything else will not be removed\n\
\n ";
#if defined(__GNUC__) or defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    m.def("RemoveHs",
          (ROMol * (*)(const ROMol &, bool, bool, bool)) MolOps::removeHs,
          "mol"_a, "implicitOnly"_a = false, "updateExplicitCount"_a = false,
          "sanitize"_a = true, docString.c_str(),
          nb::rv_policy::take_ownership);
#if defined(__GNUC__) or defined(__clang__)
#pragma GCC diagnostic pop
#endif

    // ------------------------------------------------------------------------
    docString = R"DOC(Parameters controlling which Hs are removed.)DOC";
    nb::class_<MolOps::RemoveHsParameters>(m, "RemoveHsParameters",
                                           docString.c_str())
        .def_rw("removeDegreeZero",
                &MolOps::RemoveHsParameters::removeDegreeZero,
                "hydrogens that have no bonds")
        .def_rw("removeHigherDegrees",
                &MolOps::RemoveHsParameters::removeHigherDegrees,
                "hydrogens with two (or more) bonds")
        .def_rw("removeOnlyHNeighbors",
                &MolOps::RemoveHsParameters::removeOnlyHNeighbors,
                "hydrogens with bonds only to other hydrogens")
        .def_rw("removeIsotopes", &MolOps::RemoveHsParameters::removeIsotopes,
                "hydrogens with non-default isotopes")
        .def_rw("removeAndTrackIsotopes",
                &MolOps::RemoveHsParameters::removeAndTrackIsotopes,
                "hydrogens with non-default isotopes and store "
                "them in the _isotopicHs atom property such "
                "that AddHs() can add the same isotope at "
                "a later stage")
        .def_rw("removeDummyNeighbors",
                &MolOps::RemoveHsParameters::removeDummyNeighbors,
                "hydrogens with at least one dummy-atom neighbor")
        .def_rw("removeDefiningBondStereo",
                &MolOps::RemoveHsParameters::removeDefiningBondStereo,
                "hydrogens defining bond stereochemistry")
        .def_rw("removeWithWedgedBond",
                &MolOps::RemoveHsParameters::removeWithWedgedBond,
                "hydrogens with wedged bonds to them")
        .def_rw("removeWithQuery", &MolOps::RemoveHsParameters::removeWithQuery,
                "hydrogens with queries defined")
        .def_rw("removeMapped", &MolOps::RemoveHsParameters::removeMapped,
                "mapped hydrogens")
        .def_rw("removeInSGroups", &MolOps::RemoveHsParameters::removeInSGroups,
                "hydrogens involved in SubstanceGroups")
        .def_rw("removeNonimplicit",
                &MolOps::RemoveHsParameters::removeNonimplicit, "DEPRECATED")
        .def_rw("removeHydrides", &MolOps::RemoveHsParameters::removeHydrides,
                "hydrogens with formal charge -1")
        .def_rw("removeNontetrahedralNeighbors",
                &MolOps::RemoveHsParameters::removeNontetrahedralNeighbors,
                "hydrogens with neighbors that have non-tetrahedral "
                "stereochemistry")
        .def_rw("showWarnings", &MolOps::RemoveHsParameters::showWarnings,
                "display warning messages for some classes of removed Hs")
        .def_rw("updateExplicitCount",
                &MolOps::RemoveHsParameters::updateExplicitCount, "DEPRECATED");
    m.def(
        "RemoveHs",
        (ROMol * (*)(const ROMol &, const MolOps::RemoveHsParameters &, bool)) &
            MolOps::removeHs,
        "mol"_a, "params"_a, "sanitize"_a = true,
        "Returns a copy of the molecule with Hs removed. Which Hs are "
        "removed is controlled by the params argument",
        nb::rv_policy::take_ownership);
    m.def("RemoveAllHs",
          (ROMol * (*)(const ROMol &, bool)) & MolOps::removeAllHs, "mol"_a,
          "sanitize"_a = true,
          "Returns a copy of the molecule with all Hs removed.",
          nb::rv_policy::take_ownership);
    m.def("MergeQueryHs",
          (ROMol * (*)(const ROMol &, bool, bool)) & MolOps::mergeQueryHs,
          "mol"_a, "mergeUnmappedOnly"_a = false, "mergeIsotopes"_a = false,
          "merges hydrogens into their neighboring atoms as queries",
          nb::rv_policy::take_ownership);

    docString =
        "Check to see if the molecule has query Hs, this is normally used on query molecules\n\
such as those returned from MolFromSmarts\n\
Example: \n\
      (hasQueryHs, hasUnmergeableQueryHs) = HasQueryHs(mol)\n\
\n\
if hasUnmergeableQueryHs, these query hs cannot be removed by calling\n\
MergeQueryHs";
    m.def("HasQueryHs", hasQueryHsHelper, "mol"_a, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Removes atoms matching a substructure query from a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as a substructure query\n\
\n\
    - onlyFrags: (optional) if this toggle is set, atoms will only be removed if\n\
      the entire fragment in which they are found is matched by the query.\n\
      See below for examples.\n\
      Default value is 0 (remove the atoms whether or not the entire fragment matches)\n\
\n\
    - useChirality: (optional) match the substructure query using chirality\n\
\n\
  RETURNS: a new molecule with the substructure removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - DeleteSubstructs('CCOC','OC') -> 'CC'\n\
\n\
    - DeleteSubstructs('CCOC','OC',1) -> 'CCOC'\n\
\n\
    - DeleteSubstructs('CCOCCl.Cl','Cl',1) -> 'CCOCCl'\n\
\n\
    - DeleteSubstructs('CCOCCl.Cl','Cl') -> 'CCOC'\n\
\n";
    m.def("DeleteSubstructs", deleteSubstructs, "mol"_a, "query"_a,
          "onlyFrags"_a = false, "useChirality"_a = false, docString.c_str(),
          nb::rv_policy::take_ownership);
    docString = R"DOC(Do a Murcko decomposition and return the scaffold)DOC";
    m.def("MurckoDecompose", MurckoDecompose, "mol"_a, docString.c_str(),
          nb::rv_policy::take_ownership);
    docString =
        R"DOC(Combine the atoms from two molecules to produce a third)DOC";
    m.def("CombineMols", combineMols, "mol1"_a, "mol2"_a,
          "offset"_a = RDGeom::Point3D(0, 0, 0), docString.c_str(),
          nb::rv_policy::take_ownership);

    // ------------------------------------------------------------------------
    docString =
        "Replaces atoms matching a substructure query in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as a substructure query\n\
\n\
    - replacement: the molecule to be used as the replacement\n\
\n\
    - replaceAll: (optional) if this toggle is set, all substructures matching\n\
      the query will be replaced in a single result, otherwise each result will\n\
      contain a separate replacement.\n\
      Default value is False (return multiple replacements)\n\
    - replacementConnectionPoint: (optional) index of the atom in the replacement that\n\
      the bond should be made to.\n\
    - useChirality: (optional) match the substructure query using chirality\n\
\n\
  RETURNS: a tuple of new molecules with the substructures replaced removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
    - A bond is only formed to the remaining atoms, if any, that were bonded \n\
      to the first atom in the substructure query. (For finer control over\n\
      substructure replacement, consider using ChemicalReaction.)\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceSubstructs('CCOC','O[CH3]','NC') -> ('CCNC',)\n\
\n\
    - ReplaceSubstructs('COCCOC','O[CH3]','NC') -> ('COCCNC','CNCCOC')\n\
\n\
    - ReplaceSubstructs('COCCOC','O[CH3]','NC',True) -> ('CNCCNC',)\n\
\n\
    - ReplaceSubstructs('COCCOC','O[CH3]','CN',True,1) -> ('CNCCNC',)\n\
\n\
    - ReplaceSubstructs('CCOC','[CH3]O','NC') -> ('CC.CN',)\n\
\n";
    m.def("ReplaceSubstructs", replaceSubstructures, "mol"_a, "query"_a,
          "replacement"_a, "replaceAll"_a = false,
          "replacementConnectionPoint"_a = 0, "useChirality"_a = false,
          docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Postprocesses the results of a mol.GetSubstructMatches(core) call \n\
where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups). \n\
It returns the match with the largest number of non-hydrogen matches to \n\
the terminal dummy atoms.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule GetSubstructMatches was run on\n\
\n\
    - core: the molecule used as a substructure query\n\
\n\
    - matches: the result returned by GetSubstructMatches\n\
\n\
  RETURNS: the tuple where terminal dummy atoms in the core match the largest \n\
           number of non-hydrogen atoms in mol\n";
    m.def("GetMostSubstitutedCoreMatch", getMostSubstitutedCoreMatchHelper,
          "mol"_a, "core"_a, "matches"_a, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Postprocesses the results of a mol.GetSubstructMatches(core) call \n\
where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups). \n\
It returns a copy of matches sorted by decreasing number of non-hydrogen matches \n\
to the terminal dummy atoms.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule GetSubstructMatches was run on\n\
\n\
    - core: the molecule used as a substructure query\n\
\n\
    - matches: the result returned by GetSubstructMatches\n\
\n\
  RETURNS: a copy of matches sorted by decreasing number of non-hydrogen matches \n\
           to the terminal dummy atoms\n";
    m.def("SortMatchesByDegreeOfCoreSubstitution",
          sortMatchesByDegreeOfCoreSubstitutionHelper, "mol"_a, "core"_a,
          "matches"_a, docString.c_str());

    // ------------------------------------------------------------------------
    docString = R"DOC(Adds named recursive queries to atoms
  )DOC";
    m.def("MolAddRecursiveQueries", addRecursiveQueriesHelper, "mol"_a,
          "queries"_a, "propName"_a, docString.c_str());

    docString = R"DOC(reads query definitions from a simply formatted file
  )DOC";
    m.def("ParseMolQueryDefFile", parseQueryDefFileHelper, "fileobj"_a,
          "standardize"_a = true, "delimiter"_a = "\t", "comment"_a = "//",
          "nameColumn"_a = 0, "smartsColumn"_a = 1, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Returns the molecule's topological distance matrix.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - useBO: (optional) toggles use of bond orders in calculating the distance matrix.\n\
      Default value is 0.\n\
\n\
    - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the\n\
      matrix (to return a \"Balaban\" distance matrix).\n\
      Default value is 0.\n\
\n\
    - force: (optional) forces the calculation to proceed, even if there is a cached value.\n\
      Default value is 0.\n\
\n\
    - prefix: (optional, internal use) sets the prefix used in the property cache\n\
      Default value is "
        ".\n\
\n\
  RETURNS: a Numeric array of floats with the distance matrix\n\
\n";
    m.def("GetDistanceMatrix", getDistanceMatrix, "mol"_a, "useBO"_a = false,
          "useAtomWts"_a = false, "force"_a = false, "prefix"_a = "",
          docString.c_str());
    // ------------------------------------------------------------------------
    docString =
        "Returns the molecule's 3D distance matrix.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - confId: (optional) chooses the conformer Id to use\n\
      Default value is -1.\n\
\n\
    - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the\n\
      matrix (to return a \"Balaban\" distance matrix).\n\
      Default value is 0.\n\
\n\
    - force: (optional) forces the calculation to proceed, even if there is a cached value.\n\
      Default value is 0.\n\
\n\
    - prefix: (optional, internal use) sets the prefix used in the property cache\n\
      Default value is "
        ".\n\
\n\
  RETURNS: a Numeric array of floats with the distance matrix\n\
\n";
    m.def("Get3DDistanceMatrix", get3DDistanceMatrix, "mol"_a, "confId"_a = -1,
          "useAtomWts"_a = false, "force"_a = false, "prefix"_a = "",
          docString.c_str());
    // ------------------------------------------------------------------------
    docString =
        "Returns the molecule's adjacency matrix.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - useBO: (optional) toggles use of bond orders in calculating the matrix.\n\
      Default value is 0.\n\
\n\
    - emptyVal: (optional) sets the elements of the matrix between non-adjacent atoms\n\
      Default value is 0.\n\
\n\
    - force: (optional) forces the calculation to proceed, even if there is a cached value.\n\
      Default value is 0.\n\
\n\
    - prefix: (optional, internal use) sets the prefix used in the property cache\n\
      Default value is "
        ".\n\
\n\
  RETURNS: a Numeric array of floats containing the adjacency matrix\n\
\n";
    m.def("GetAdjacencyMatrix", getAdjacencyMatrix, "mol"_a, "useBO"_a = false,
          "emptyVal"_a = 0, "force"_a = false, "prefix"_a = "",
          docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        R"DOC(Kekulizes the molecule

  ARGUMENTS:

    - mol: the molecule to use

    - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the
      molecule will be marked non-aromatic following the kekulization.
      Default value is False.

    - canonical: (optional) if true, uses canonical atom ranking so
      that the kekulization result is independent of the atom ordering in the
      molecule.  Set to false to skip the ranking step for better performance
      when deterministic output is not required (e.g. during sanitization).
      Note, this "canonical" order only really makes sense when the molecule's
      chemistry is sane, like after sanitization. If stereochemistry hasn't been
      perceived, the chemistry of the molecule is inconsistent, and
      "canonical" atom ranks are only a technical artifact.

  NOTES:

    - The molecule is modified in place.

    - this does not modify query bonds which have bond type queries (like those
      which come from SMARTS) or rings containing them.

    - even if clearAromaticFlags is False the BondType for all modified
      aromatic bonds will be changed from AROMATIC to SINGLE or DOUBLE
      Kekulization.

)DOC";
    m.def("Kekulize", kekulizeMol, "mol"_a, "clearAromaticFlags"_a = false,
          "canonical"_a = true, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        R"DOC(Kekulizes the molecule if possible. Otherwise the molecule is not modified

  ARGUMENTS:

    - mol: the molecule to use

    - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the 
      molecule will be marked non-aromatic if the kekulization succeds.
      Default value is False.

    - canonical: (optional) if true  uses canonical atom ranking so
      that the kekulization result is independent of the atom ordering in the
      molecule.  Set to false to skip the ranking step for better performance
      when deterministic output is not required (e.g. during sanitization).
      Note, this "canonical" order only really makes sense when the molecule's
      chemistry is sane, like after sanitization. If stereochemistry hasn't been
      perceived, the chemistry of the molecule is inconsistent, and
      "canonical" atom ranks are only a technical artifact.

\n\
    - canonical: (optional) if True, uses canonical atom ranking so that the\n\
      kekulization result is independent of the atom ordering in the molecule.\n\
      Default value is False.\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
    )DOC";
    m.def("KekulizeIfPossible", kekulizeMolIfPossible, "mol"_a,
          "clearAromaticFlags"_a = false, "canonical"_a = true,
          docString.c_str());
    // ------------------------------------------------------------------------
    docString =
        R"DOC(cleans up certain common bad functionalities in the molecule

  ARGUMENTS:

    - mol: the molecule to use

  NOTES:

    - The molecule is modified in place.
)DOC";
    m.def("Cleanup", cleanupMol, "mol"_a, docString.c_str());
    // ------------------------------------------------------------------------
    docString =
        R"DOC(removes bogus chirality markers (e.g. tetrahedral flags on non-sp3 centers)

  ARGUMENTS:

    - mol: the molecule to use

  NOTES:

    - The molecule is modified in place.
)DOC";
    m.def("CleanupChirality", cleanupChiralityMol, "mol"_a, docString.c_str());
    // ------------------------------------------------------------------------
    docString =
        R"DOC(removes bogus atropisomeric markers (e.g. those without sp2 begin and end atoms)

  ARGUMENTS:

    - mol: the molecule to use

  NOTES:

    - The molecule is modified in place.
)DOC";
    m.def("CleanupAtropisomers", cleanupAtropisomersMol, "mol"_a,
          docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "cleans up certain common bad functionalities in the organometallic molecule\n\
\n\
  Note that this function is experimental and may either change in behavior\n\
  or be replaced with something else in future releases.\n\
\n\
        ARGUMENTS :\n\
\n - mol : the molecule to use\n\
\n NOTES :\n\
\n - The molecule is modified in place.\n\
\n ";
    m.def("CleanupOrganometallics", cleanUpOrganometallicsMol, "mol"_a,
          docString.c_str());

    nb::enum_<MolOps::AromaticityModel>(m, "AromaticityModel")
        .value("AROMATICITY_DEFAULT", MolOps::AROMATICITY_DEFAULT)
        .value("AROMATICITY_RDKIT", MolOps::AROMATICITY_RDKIT)
        .value("AROMATICITY_SIMPLE", MolOps::AROMATICITY_SIMPLE)
        .value("AROMATICITY_MDL", MolOps::AROMATICITY_MDL)
        .value("AROMATICITY_MMFF94", MolOps::AROMATICITY_MMFF94)
        .value("AROMATICITY_CUSTOM", MolOps::AROMATICITY_CUSTOM)
        .export_values();

    // ------------------------------------------------------------------------
    docString =
        "does aromaticity perception\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - model: the model to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
    m.def("SetAromaticity", setAromaticityMol, "mol"_a,
          "model"_a = MolOps::AROMATICITY_DEFAULT, docString.c_str());

    docString =
        "finds conjugated bonds\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
    m.def("SetConjugation", setConjugationMol, "mol"_a, docString.c_str());
    docString =
        "Assigns hybridization states to atoms\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
    m.def("SetHybridization", setHybridizationMol, "mol"_a, docString.c_str());
    docString =
        "Assigns radical counts to atoms\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
    m.def("AssignRadicals", assignRadicalsMol, "mol"_a, docString.c_str());

    docString =
        R"DOC(One way of showing haptic bonds (such as cyclopentadiene to
iron in ferrocene) is to use a dummy atom with a dative bond to the
iron atom with the bond labelled with the atoms involved in the
organic end of the bond.  Another way is to have explicit dative
bonds from the atoms of the haptic group to the metal atom.  This
function converts the former representation to the latter.

ARGUMENTS:

  - mol: the molecule to use

RETURNS:
  a modified copy of the molecule)DOC";
    m.def("HapticBondsToDative", hapticBondsToDativeHelper, "mol"_a,
          docString.c_str(), nb::rv_policy::take_ownership);

    docString =
        R"DOC(Does the reverse of hapticBondsToDative.  If there are multiple
contiguous atoms attached by dative bonds to an atom (probably a metal
atom), the dative bonds will be replaced by a dummy atom in their
centre attached to the (metal) atom by a dative bond, which is
labelled with ENDPTS of the atoms that had the original dative bonds.

ARGUMENTS:

  - mol: the molecule to use

RETURNS:
  a modified copy of the molecule)DOC";
    m.def("DativeBondsToHaptic", dativeBondsToHapticHelper, "mol"_a,
          docString.c_str(), nb::rv_policy::take_ownership);

    // ------------------------------------------------------------------------
    docString =
        "Finds all subgraphs of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target number of bonds for the subgraphs.\n\
\n\
    - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph\n\
      should be included in the results.\n\
      Defaults to 0.\n\
\n\
    - rootedAtAtom: (optional) if nonzero, only subgraphs from the specified\n\
      atom will be returned.\n\
\n\
  RETURNS: a tuple of 2-tuples with bond IDs\n\
\n\
  NOTES: \n\
\n\
   - Difference between _subgraphs_ and _paths_ :: \n\
\n\
       Subgraphs are potentially branched, whereas paths (in our \n\
       terminology at least) cannot be.  So, the following graph: \n\
\n\
            C--0--C--1--C--3--C\n\
                  |\n\
                  2\n\
                  |\n\
                  C\n\
  has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)\n\
  but only 2 _paths_ of length 3: (0,1,3),(2,1,3)\n\
\n";
    m.def("FindAllSubgraphsOfLengthN", &findAllSubgraphsOfLengthN, "mol"_a,
          "length"_a, "useHs"_a = false, "rootedAtAtom"_a = -1,
          docString.c_str());
    // ------------------------------------------------------------------------
    docString =
        "Finds all subgraphs of a particular length in a molecule\n\
  See documentation for FindAllSubgraphsOfLengthN for definitions\n\
\n";
    m.def("FindAllSubgraphsOfLengthMToN", &findAllSubgraphsOfLengthsMtoNHelper,
          "mol"_a, "min"_a, "max"_a, "useHs"_a = false, "rootedAtAtom"_a = -1,
          docString.c_str());
    // ------------------------------------------------------------------------
    docString =
        "Finds unique subgraphs of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target number of bonds for the subgraphs.\n\
\n\
    - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph\n\
      should be included in the results.\n\
      Defaults to 0.\n\
\n\
    - useBO: (optional) Toggles use of bond orders in distinguishing one subgraph from\n\
      another.\n\
      Defaults to 1.\n\
\n\
    - rootedAtAtom: (optional) if nonzero, only subgraphs from the specified\n\
      atom will be returned.\n\
\n\
  RETURNS: a tuple of tuples with bond IDs\n\
\n\
\n";
    m.def("FindUniqueSubgraphsOfLengthN", &findUniqueSubgraphsOfLengthN,
          "mol"_a, "length"_a, "useHs"_a = false, "useBO"_a = true,
          "rootedAtAtom"_a = -1, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Finds all paths of a particular length in a molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - length: an integer with the target length for the paths.\n\
\n\
    - useBonds: (optional) toggles the use of bond indices in the paths.\n\
      Otherwise atom indices are used.  *Note* this behavior is different\n\
      from that for subgraphs.\n\
      Defaults to 1.\n\
\n\
    - rootedAtAtom: (optional) if nonzero, only paths from the specified\n\
      atom will be returned.\n\
\n\
    - onlyShortestPaths: (optional) if set then only paths which are <= the shortest\n\
      path between the begin and end atoms will be included in the results\n\
\n\
  RETURNS: a tuple of tuples with IDs for the bonds.\n\
\n\
  NOTES: \n\
\n\
   - Difference between _subgraphs_ and _paths_ :: \n\
\n\
       Subgraphs are potentially branched, whereas paths (in our \n\
       terminology at least) cannot be.  So, the following graph: \n\
\n\
            C--0--C--1--C--3--C\n\
                  |\n\
                  2\n\
                  |\n\
                  C\n\
\n\
       has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)\n\
       but only 2 _paths_ of length 3: (0,1,3),(2,1,3)\n\
\n";
    m.def("FindAllPathsOfLengthN", &findAllPathsOfLengthN, "mol"_a, "length"_a,
          "useBonds"_a = true, "useHs"_a = false, "rootedAtAtom"_a = -1,
          "onlyShortestPaths"_a = false, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Find bonds of a particular radius around an atom. \n\
         Return empty result if there is no bond at the requested radius.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - radius: an integer with the target radius for the environment.\n\
\n\
    - rootedAtAtom: the atom to consider\n\
\n\
    - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph\n\
      should be included in the results.\n\
      Defaults to 0.\n\
\n\
    - enforceSize (optional) If set to False, all bonds within the requested radius is \n\
      collected. Defaults to 1. \n\
\n\
    - atomMap: (optional) If provided, it will measure the minimum distance of the atom \n\
      from the rooted atom (start with 0 from the rooted atom). The result is a pair of \n\
      the atom ID and the distance. \n\
\n\
  RETURNS: a vector of bond IDs\n\
\n";
    m.def("FindAtomEnvironmentOfRadiusN", findAtomEnvironmentOfRadiusNHelper,
          "mol"_a, "radius"_a, "rootedAtAtom"_a, "useHs"_a = false,
          "enforceSize"_a = true, "atomMap"_a = nb::none(), docString.c_str());

    m.def("PathToSubmol", pathToSubmolHelper, "mol"_a, "path"_a,
          "useQuery"_a = false, "atomMap"_a = nb::none(), "",
          nb::rv_policy::take_ownership);

    // ------------------------------------------------------------------------
    docString =
        "Finds the disconnected fragments from a molecule.\n\
\n\
  For example, for the molecule 'CC(=O)[O-].[NH3+]C' GetMolFrags() returns\n\
  ((0, 1, 2, 3), (4, 5))\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - asMols: (optional) if this is provided and true, the fragments\n\
      will be returned as molecules instead of atom ids.\n\
    - sanitizeFrags: (optional) if this is provided and true, the fragments\n\
      molecules will be sanitized before returning them.\n\
    - frags: (optional, defaults to None) if asMols is true and this is provided\n\
       as an empty list, the result will be mol.GetNumAtoms() long on return and\n\
       will contain the fragment assignment for each Atom\n\
    - fragsMolAtomMapping: (optional, defaults to None) if asMols is true and this\n\
      is provided as an empty list, the result will be numFrags long on \n\
      return, and each entry will contain the indices of the Atoms in that fragment:\n\
      [(0, 1, 2, 3), (4, 5)]\n\
\n\
  RETURNS: a tuple of tuples with IDs for the atoms in each fragment\n\
           or a tuple of molecules.\n\
\n";
    m.def("GetMolFrags", &GetMolFragsWithMapping, "mol"_a, "asMols"_a = false,
          "sanitizeFrags"_a = true, "frags"_a = nb::none(),
          "fragsMolAtomMapping"_a = nb::none(), docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Splits a molecule into pieces based on PDB residue information.\n\
          \n\
          ARGUMENTS:\n\
          \n\
          - mol: the molecule to use\n\
          - whiteList: only residues in this list will be returned\n\
          - negateList: if set, negates the white list inclusion logic\n\
          \n\
          RETURNS: a dictionary keyed by residue name with molecules as the values\n\
          \n";
    m.def("SplitMolByPDBResidues", &splitMolByPDBResidues, "mol"_a,
          "whiteList"_a = nb::none(), "negateList"_a = false,
          docString.c_str());
    // ------------------------------------------------------------------------
    docString =
        "Splits a molecule into pieces based on PDB chain information.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - whiteList: only residues in this list will be returned\n\
    - negateList: if set, negates the white list inclusion logic\n\
\n\
  RETURNS: a dictionary keyed by chain id with molecules as the values\n\
\n";
    m.def("SplitMolByPDBChainId", &splitMolByPDBChainId, "mol"_a,
          "whiteList"_a = nb::none(), "negateList"_a = false,
          docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Returns the formal charge for the molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n";
    m.def("GetFormalCharge", &MolOps::getFormalCharge, docString.c_str(),
          "mol"_a);

    // ------------------------------------------------------------------------
    docString =
        "Find the shortest path between two atoms using the Bellman-Ford algorithm.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - idx1: index of the first atom\n\
    - idx2: index of the second atom\n\
\n";
    m.def("GetShortestPath", getShortestPathHelper, docString.c_str(), "mol"_a,
          "aid1"_a, "aid2"_a);

    // ------------------------------------------------------------------------
    docString =
        R"DOC(Assign stereochemistry tags to atoms and bonds.
  If useLegacyStereoPerception is true, it also does the CIP stereochemistry
  assignment for the molecule's atoms (R/S) and double bonds (Z/E).
  This assignment is based on legacy code which is fast, but is
  known to incorrectly assign CIP labels in some cases.
  instead, to assign CIP labels based on an accurate, though slower,
  implementation of the CIP rules, call CIPLabeler::assignCIPLabels().
  Chiral atoms will have a property '_CIPCode' indicating their chiral code.

  ARGUMENTS:

    - mol: the molecule to use
    - cleanIt: (optional) if provided, any existing values of the property `_CIPCode`
        will be cleared, atoms with a chiral specifier that aren't
      actually chiral (e.g. atoms with duplicate substituents or only 2 substituents,
      etc.) will have their chiral code set to CHI_UNSPECIFIED. Bonds with 
      STEREOCIS/STEREOTRANS specified that have duplicate substituents based upon the CIP 
      atom ranks will be marked STEREONONE. 
    - force: (optional) causes the calculation to be repeated, even if it has already
      been done
    - flagPossibleStereoCenters (optional)   set the _ChiralityPossible property on
      atoms that are possible stereocenters
)DOC";
    m.def("AssignStereochemistry", MolOps::assignStereochemistry, "mol"_a,
          "cleanIt"_a = false, "force"_a = false,
          "flagPossibleStereoCenters"_a = false, docString.c_str());

    m.def("ComputeAtomCIPRanks", computeAtomCIPRanksHelper, "mol"_a,
          R"DOC(Computes the CIP ranks for the atoms in a molecule.
  The ranks are stored as an atom property '_CIPRank' and returned as a tuple.
)DOC");

    // ------------------------------------------------------------------------
    docString =
        "Uses bond directions to assign ChiralTypes to a molecule's atoms.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - confId: (optional) the conformation to use \n\
    - replaceExistingTags: (optional) replace any existing information about stereochemistry\n\
\n";
    m.def("AssignChiralTypesFromBondDirs",
          MolOps::assignChiralTypesFromBondDirs, "mol"_a, "confId"_a = -1,
          "replaceExistingTags"_a = true, docString.c_str());
    // ------------------------------------------------------------------------
    docString =
        "Uses a conformer (should be 3D) to assign ChiralTypes to a molecule's atoms\n\
        and stereo flags to its bonds\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - confId: (optional) the conformation to use \n\
    - replaceExistingTags: (optional) replace any existing information about stereochemistry\n\
\n";
    m.def("AssignStereochemistryFrom3D", MolOps::assignStereochemistryFrom3D,
          "mol"_a, "confId"_a = -1, "replaceExistingTags"_a = true,
          docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Find bonds than can be cis/trans in a molecule and mark them as 'any'.\n\
         This function finds any double bonds that can potentially be part\n\
         of a cis/trans system. No attempt is made here to mark them cis or trans\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - cleanIt: (optional) if this option is set to true, any previous marking of _CIPCode\n\
               on the bond is cleared - otherwise it is left untouched\n\
\n";
    m.def("FindPotentialStereoBonds", MolOps::findPotentialStereoBonds, "mol"_a,
          "cleanIt"_a = false, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Removes all stereochemistry info from the molecule.\n\
\n";
    m.def("RemoveStereochemistry", MolOps::removeStereochemistry, "mol"_a,
          docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Sets the chiral tags on a molecule's atoms based on\n\
  a 3D conformation.\n\
  NOTE that this does not check to see if atoms are chiral centers (i.e. all\n\
  substituents are different), it merely sets the chiral type flags based on the\n\
  coordinates and atom ordering. Use AssignStereochemistryFrom3D() if you\n\
  want chiral flags only on actual stereocenters.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - confId: the conformer id to use, -1 for the default \n\
    - replaceExistingTags: if True, existing stereochemistry information will be cleared\n\
    before running the calculation.\n\
\n";
    m.def("AssignAtomChiralTagsFromStructure", MolOps::assignChiralTypesFrom3D,
          "mol"_a, "confId"_a = -1, "replaceExistingTags"_a = true,
          docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Sets the chiral tags on a molecule's atoms based on\n\
  the molParity atom property.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - replaceExistingTags: if True, existing stereochemistry information will be cleared\n\
    before running the calculation.\n\
\n";
    m.def("AssignAtomChiralTagsFromMolParity",
          MolOps::assignChiralTypesFromMolParity, "mol"_a,
          "replaceExistingTags"_a = true, docString.c_str());

    docString = R"DOC(returns the meso centers in a molecule (if any).
    
  ARGUMENTS:
    
    - mol: the molecule to use
    - includeIsotopes: (optional) toggles whether or not isotopes should be included in the
      calculation.
    - includeAtomMaps: (optional) toggles whether or not atom maps should be included in the
      calculation.
    )DOC";
    m.def("FindMesoCenters", findMesoHelper, "mol"_a,
          "includeIsotopes"_a = true, "includeAtomMaps"_a = false,
          docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Returns an RDKit topological fingerprint for a molecule\n\
\n\
  Explanation of the algorithm below.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - minPath: (optional) minimum number of bonds to include in the subgraphs\n\
      Defaults to 1.\n\
\n\
    - maxPath: (optional) maximum number of bonds to include in the subgraphs\n\
      Defaults to 7.\n\
\n\
    - fpSize: (optional) number of bits in the fingerprint\n\
      Defaults to 2048.\n\
\n\
    - nBitsPerHash: (optional) number of bits to set per path\n\
      Defaults to 2.\n\
\n\
    - useHs: (optional) include paths involving Hs in the fingerprint if the molecule\n\
      has explicit Hs.\n\
      Defaults to True.\n\
\n\
    - tgtDensity: (optional) fold the fingerprint until this minimum density has\n\
      been reached\n\
      Defaults to 0.\n\
\n\
    - minSize: (optional) the minimum size the fingerprint will be folded to when\n\
      trying to reach tgtDensity\n\
      Defaults to 128.\n\
\n\
    - branchedPaths: (optional) if set both branched and unbranched paths will be\n\
      used in the fingerprint.\n\
      Defaults to True.\n\
\n\
    - useBondOrder: (optional) if set both bond orders will be used in the path hashes\n\
      Defaults to True.\n\
\n\
    - atomInvariants: (optional) a sequence of atom invariants to use in the path hashes\n\
      Defaults to empty.\n\
\n\
    - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs \n\
      starting from these atoms will be used.\n\
      Defaults to empty.\n\
\n\
    - atomBits: (optional) an empty list. If provided, the result will contain a list \n\
      containing the bits each atom sets.\n\
      Defaults to empty.\n\
\n\
    - bitInfo: (optional) an empty dict. If provided, the result will contain a dict \n\
      with bits as keys and corresponding bond paths as values.\n\
      Defaults to empty.\n\
\n\
  RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits\n\
\n\
  ALGORITHM:\n\
\n\
   This algorithm functions by find all subgraphs between minPath and maxPath in\n\
   length.  For each subgraph:\n\
\n\
     1) A hash is calculated.\n\
\n\
     2) The hash is used to seed a random-number generator\n\
\n\
     3) _nBitsPerHash_ random numbers are generated and used to set the corresponding\n\
        bits in the fingerprint\n\
\n\
\n";
    m.def("RDKFingerprint", wrapRDKFingerprintMol, "mol"_a, "minPath"_a = 1,
          "maxPath"_a = 7, "fpSize"_a = 2048, "nBitsPerHash"_a = 2,
          "useHs"_a = true, "tgtDensity"_a = 0.0, "minSize"_a = 128,
          "branchedPaths"_a = true, "useBondOrder"_a = true,
          "atomInvariants"_a = 0, "fromAtoms"_a = 0, "atomBits"_a = nb::none(),
          "bitInfo"_a = nb::none(), docString.c_str(),
          nb::rv_policy::take_ownership);
    m.attr("_RDKFingerprint_version") = RDKit::RDKFingerprintMolVersion;

    docString =
        "Returns an unfolded count-based version of the RDKit fingerprint for a molecule\n\
\n\
ARGUMENTS:\n\
    \n\
        - mol: the molecule to use\n\
    \n\
        - minPath: (optional) minimum number of bonds to include in the subgraphs\n\
          Defaults to 1.\n\
    \n\
        - maxPath: (optional) maximum number of bonds to include in the subgraphs\n\
          Defaults to 7.\n\
    \n\
        - useHs: (optional) include paths involving Hs in the fingerprint if the molecule\n\
          has explicit Hs.\n\
          Defaults to True.\n\
    \n\
        - branchedPaths: (optional) if set both branched and unbranched paths will be\n\
          used in the fingerprint.\n\
          Defaults to True.\n\
    \n\
        - useBondOrder: (optional) if set both bond orders will be used in the path hashes\n\
          Defaults to True.\n\
    \n\
        - atomInvariants: (optional) a sequence of atom invariants to use in the path hashes\n\
          Defaults to empty.\n\
    \n\
        - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs \n\
          starting from these atoms will be used.\n\
          Defaults to empty.\n\
    \n\
        - atomBits: (optional) an empty list. If provided, the result will contain a list \n\
          containing the bits each atom sets.\n\
          Defaults to empty.\n\
    \n\
        - bitInfo: (optional) an empty dict. If provided, the result will contain a dict \n\
          with bits as keys and corresponding bond paths as values.\n\
          Defaults to empty.\n\
     \n\
     \n";

    m.def("UnfoldedRDKFingerprintCountBased", wrapUnfoldedRDKFingerprintMol,
          "mol"_a, "minPath"_a = 1, "maxPath"_a = 7, "useHs"_a = true,
          "branchedPaths"_a = true, "useBondOrder"_a = true,
          "atomInvariants"_a = 0, "fromAtoms"_a = 0, "atomBits"_a = nb::none(),
          "bitInfo"_a = nb::none(), docString.c_str(),
          nb::rv_policy::take_ownership);

    // ------------------------------------------------------------------------
    docString =
        "Returns a layered fingerprint for a molecule\n\
\n\
  NOTE: This function is experimental. The API or results may change from\n\
    release to release.\n\
\n\
  Explanation of the algorithm below.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - layerFlags: (optional) which layers to include in the fingerprint\n\
      See below for definitions. Defaults to all.\n\
\n\
    - minPath: (optional) minimum number of bonds to include in the subgraphs\n\
      Defaults to 1.\n\
\n\
    - maxPath: (optional) maximum number of bonds to include in the subgraphs\n\
      Defaults to 7.\n\
\n\
    - fpSize: (optional) number of bits in the fingerprint\n\
      Defaults to 2048.\n\
\n\
    - atomCounts: (optional) \n\
      if provided, this should be a list at least as long as the number of atoms\n\
      in the molecule. It will be used to provide the count of the number \n\
      of paths that set bits each atom is involved in.\n\
      NOTE: the list is not zeroed out here.\n\
\n\
    - setOnlyBits: (optional) \n\
      if provided, only bits that are set in this bit vector will be set\n\
      in the result. This is essentially the same as doing:\n\
           res &= setOnlyBits\n\
      but also has an impact on the atomCounts (if being used)\n\
\n\
    - branchedPaths: (optional) if set both branched and unbranched paths will be\n\
      used in the fingerprint.\n\
      Defaults to True.\n\
\n\
    - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs \n\
      starting from these atoms will be used.\n\
      Defaults to empty.\n\
\n\
  RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits\n\
\n\
  Layer definitions:\n\
     - 0x01: pure topology\n\
     - 0x02: bond order\n\
     - 0x04: atom types\n\
     - 0x08: presence of rings\n\
     - 0x10: ring sizes\n\
     - 0x20: aromaticity\n\
\n\
\n";
    m.def("LayeredFingerprint", wrapLayeredFingerprint, "mol"_a,
          "layerFlags"_a = 0xFFFFFFFF, "minPath"_a = 1, "maxPath"_a = 7,
          "fpSize"_a = 2048, "atomCounts"_a = nb::list(),
          "setOnlyBits"_a = (ExplicitBitVect *)nullptr,
          "branchedPaths"_a = true, "fromAtoms"_a = 0, docString.c_str(),
          nb::rv_policy::take_ownership);
    m.attr("_LayeredFingerprint_version") = RDKit::LayeredFingerprintMolVersion;
    m.attr("LayeredFingerprint_substructLayers") = RDKit::substructLayers;

    // ------------------------------------------------------------------------
    docString =
        "A fingerprint using SMARTS patterns \n\
\n\
  NOTE: This function is experimental. The API or results may change from\n\
    release to release.\n";
    m.def("PatternFingerprint", wrapPatternFingerprint, "mol"_a,
          "fpSize"_a = 2048, "atomCounts"_a = nb::list(),
          "setOnlyBits"_a = (ExplicitBitVect *)nullptr,
          "tautomerFingerprints"_a = false, docString.c_str(),
          nb::rv_policy::take_ownership);
    m.attr("_PatternFingerprint_version") = RDKit::PatternFingerprintMolVersion;
    docString =
        "A fingerprint using SMARTS patterns \n\
\n\
  NOTE: This function is experimental. The API or results may change from\n\
    release to release.\n";
    m.def("PatternFingerprint", wrapPatternFingerprintBundle, "mol"_a,
          "fpSize"_a = 2048, "setOnlyBits"_a = (ExplicitBitVect *)nullptr,
          "tautomerFingerprints"_a = false, docString.c_str(),
          nb::rv_policy::take_ownership);

    nb::class_<Chirality::BondWedgingParameters>(
        m, "BondWedgingParameters",
        "Parameters controlling how bond wedging is done.")
        .def_rw("wedgeTwoBondsIfPossible",
                &Chirality::BondWedgingParameters::wedgeTwoBondsIfPossible,
                R"DOC(If this is enabled then two bonds will be wedged at chiral
  centers subject to the following constraints:
    1. ring bonds will not be wedged
    2. bonds to chiral centers will not be wedged
    3. bonds separated by more than 120 degrees will not be
        wedged)DOC");
    docString =
        "Set the wedging on single bonds in a molecule.\n\
   The wedging scheme used is that from Mol files.\n\
\n\
  ARGUMENTS:\n\
\n\
    - molecule: the molecule to update\n\
    - conformer: the conformer to use to determine wedge direction\n\
\n\
\n";
    m.def("WedgeMolBonds", Chirality::wedgeMolBonds, "mol"_a, "conformer"_a,
          "params"_a = nb::none(), docString.c_str());

    docString =
        "Set the wedging to that which was read from the original\n\
     MolBlock, over-riding anything that was originally there.\n\
\n\
          ARGUMENTS:\n\
        \n\
            - molecule: the molecule to update\n\
            - allBondTypes: reapply the wedging also on bonds other\n\
              than single and aromatic ones\n\
            - verify: if true, the function will check that the wedges are only\n\
              applied in sensible places (i.e.single bonds connected to chiral\n\
              centers or atropisomeric bonds)\n\
        \n\
        \n";
    m.def("ReapplyMolBlockWedging", reapplyWedging, "mol"_a,
          "allBondTypes"_a = true, "verify"_a = false, docString.c_str());

    docString =
        "Remove chiral markings that were derived from a 3D mol but were not \n\
        explicity marked in the mol block. (wedge bond or CFG indication\n\
        \n\
          ARGUMENTS:\n\
        \n\
            - molecule: the molecule to update\n\
        \n\
        \n";
    m.def("RemoveNonExplicit3DChirality",
          Chirality::removeNonExplicit3DChirality, "mol"_a, docString.c_str());

    nb::enum_<RDKit::StereoGroupAbsOptions>(m, "StereoGroupAbsOptions")
        .value("OnlyIncludeWhenOtherGroupsExist",
               RDKit::StereoGroupAbsOptions::OnlyIncludeWhenOtherGroupsExist)
        .value("NeverInclude", RDKit::StereoGroupAbsOptions::NeverInclude)
        .value("AlwaysInclude", RDKit::StereoGroupAbsOptions::AlwaysInclude);

    docString =
        "Rationalize Enhanced Stereo indications to a canonical form \n\
        \n\
          ARGUMENTS:\n\
        \n\
            - molecule: the molecule to update\n\
            -StereoGroupAbsOptions outputAbsoluteGroups: controls output of abs groups: \n\
              one of: OnlyIncludeWhenOtherGroupsExist, NeverInclude, AlwaysInclude \n\
             maxStereoGroups: maximm number of OR or AND stereo groups to process (default is 12): \n\
        \n\
        \n ";
    m.def("CanonicalizeStereoGroups", canonicalizeStereoGroupsHelper, "mol"_a,
          "outputAbsoluteGroups"_a =
              RDKit::StereoGroupAbsOptions::OnlyIncludeWhenOtherGroupsExist,
          "maxStereoGroups"_a = 12, docString.c_str(),
          nb::rv_policy::take_ownership);

    docString =
        R"DOC(Constants used to set the thresholds for which single bonds can be made wavy.)DOC";
    nb::class_<StereoBondThresholds>(m, "StereoBondThresholds",
                                     docString.c_str())
        .def_ro_static("DBL_BOND_NO_STEREO",
                       &StereoBondThresholds::DBL_BOND_NO_STEREO,
                       "neighboring double bond without stereo info")
        .def_ro_static("DBL_BOND_SPECIFIED_STEREO",
                       &StereoBondThresholds::DBL_BOND_SPECIFIED_STEREO,
                       "neighboring double bond with stereo specified")
        .def_ro_static("CHIRAL_ATOM", &StereoBondThresholds::CHIRAL_ATOM,
                       "atom with specified chirality")
        .def_ro_static("DIRECTION_SET", &StereoBondThresholds::DIRECTION_SET,
                       "single bond with the direction already set");

    docString = R"DOC(set wavy bonds around double bonds with STEREOANY stereo
  ARGUMENTS :
    - molecule : the molecule to update\n -
    - conformer : the conformer to use to determine wedge direction
)DOC";
    m.def("AddWavyBondsForStereoAny", addWavyBondsForStereoAny, "mol"_a,
          "clearDoubleBondFlags"_a = true,
          "addWhenImpossible"_a = StereoBondThresholds::DBL_BOND_NO_STEREO,
          docString.c_str());

    docString =
        R"DOC(Set the wedging on an individual bond from a molecule.
   The wedging scheme used is that from Mol files.
  ARGUMENTS:
    - bond: the bond to update
    - atom ID: the atom from which to do the wedging
    - conformer: the conformer to use to determine wedge direction
)DOC";
    m.def("WedgeBond", Chirality::wedgeBond, docString.c_str(), "bond"_a,
          "fromAtomIdx"_a, "conf"_a);

    // ------------------------------------------------------------------------
    docString =
        "Replaces sidechains in a molecule with dummy atoms for their attachment points.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - coreQuery: the molecule to be used as a substructure query for recognizing the core\n\
\n\
    - useChirality: (optional) match the substructure query using chirality\n\
\n\
  RETURNS: a new molecule with the sidechains removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceSidechains('CCC1CCC1','C1CCC1') -> '[Xa]C1CCC1'\n\
\n\
    - ReplaceSidechains('CCC1CC1','C1CCC1') -> ''\n\
\n\
    - ReplaceSidechains('C1CC2C1CCC2','C1CCC1') -> '[Xa]C1CCC1[Xb]'\n\
\n";
    m.def("ReplaceSidechains", replaceSidechains, "mol"_a, "coreQuery"_a,
          "useChirality"_a = false, docString.c_str(),
          nb::rv_policy::take_ownership);

    // ------------------------------------------------------------------------
    docString =
        "Removes the core of a molecule and labels the sidechains with dummy atoms based on\n\
The matches indices given in the matching vector matches.\n\
Calling:\n\
  ReplaceCore(mol,core,mol.GetSubstructMatch(core))\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - coreQuery: the molecule to be used as a substructure query for recognizing the core\n\
\n\
    - matches: a matching vector of the type returned by mol.GetSubstructMatch(...)\n\
\n\
    - replaceDummies: toggles replacement of atoms that match dummies in the query\n\
\n\
    - labelByIndex: toggles labeling the attachment point dummy atoms with \n\
      the index of the core atom they're attached to.\n\
\n\
    - requireDummyMatch: if the molecule has side chains that attach at points not\n\
      flagged with a dummy, it will be rejected (None is returned)\n\
\n\
  RETURNS: a new molecule with the core removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
EXAMPLES:\n\n\
    >>> from rdkit.Chem import MolToSmiles, MolFromSmiles, ReplaceCore\n\
    >>> mol = MolFromSmiles('C1ONNCC1')\n\
    >>> core = MolFromSmiles('NN')\n\
\n\
    >>> MolToSmiles(ReplaceCore(mol, core, mol.GetSubstructMatch(core)))\n\
    '[1*]OCCC[2*]'\n\
\n\
    Since NN is symmetric, we should actually get two matches here if we don't\n\
    uniquify the matches.\n\n\
    >>> [MolToSmiles(ReplaceCore(mol, core, match))\n\
    ...     for match in mol.GetSubstructMatches(core, uniquify=False)]\n\
    ['[1*]OCCC[2*]', '[1*]CCCO[2*]']\n\
\n\
";
    m.def("ReplaceCore", replaceCoreHelper, "mol"_a, "core"_a, "matches"_a,
          "replaceDummies"_a = true, "labelByIndex"_a = false,
          "requireDummyMatch"_a = false, docString.c_str(),
          nb::rv_policy::take_ownership);
    // ------------------------------------------------------------------------
    docString =
        "Removes the core of a molecule and labels the sidechains with dummy atoms.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - coreQuery: the molecule to be used as a substructure query for recognizing the core\n\
\n\
    - replaceDummies: toggles replacement of atoms that match dummies in the query\n\
\n\
    - labelByIndex: toggles labeling the attachment point dummy atoms with \n\
      the index of the core atom they're attached to.\n\
\n\
    - requireDummyMatch: if the molecule has side chains that attach at points not\n\
      flagged with a dummy, it will be rejected (None is returned)\n\
\n\
    - useChirality: use chirality matching in the coreQuery\n\
\n\
  RETURNS: a new molecule with the core removed\n\
\n\
  NOTES:\n\
\n\
    - The original molecule is *not* modified.\n\
\n\
  EXAMPLES:\n\
\n\
   >>> from rdkit.Chem import MolToSmiles, MolFromSmiles, MolFromSmarts, ReplaceCore\n\
\n\
   Basic usage: remove a core as specified by SMILES (or another molecule).\n\
   To get the atom labels which are stored as an isotope of the matched atom, \n\
   the output must be written as isomeric smiles.  \n\
   A small confusion is that atom isotopes of 0 aren't shown in smiles strings.\n\
\n\
   Here we remove a ring and leave the decoration (r-group) behind.\n\
\n\
   >>> MolToSmiles(ReplaceCore(MolFromSmiles('CCCC1CCC1'),MolFromSmiles('C1CCC1')),\n\
   ...             isomericSmiles=True)\n\
   '[1*]CCC'\n\
\n\
   The isotope label by default is matched by the first connection found. In order to\n\
   indicate which atom the decoration is attached in the core query, use labelByIndex=True.\n\
   Here the attachment is from the third atom in the smiles string, which is indexed by 3\n\
   in the core, like all good computer scientists expect, atoms indices start at 0.\n\n\
   >>> MolToSmiles(ReplaceCore(MolFromSmiles('CCN1CCC1'),MolFromSmiles('C1CCN1'),\n\
   ...                         labelByIndex=True),\n\
   ...   isomericSmiles=True)\n\
   '[3*]CC'\n\
\n\
   Non-core matches just return None\n\n\
   >>> ReplaceCore(MolFromSmiles('CCC1CC1'),MolFromSmiles('C1CCC1'))\n\
\n\
   The bond between atoms are considered part of the core and are removed as well\n\n\
   >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CC2C1CCC2'),MolFromSmiles('C1CCC1')),\n\
   ...             isomericSmiles=True)\n\
   '[1*]CCC[2*]'\n\
   >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmiles('N')),\n\
   ...             isomericSmiles=True)\n\
   '[1*]CCCC[2*]'\n\
\n\
   When using dummy atoms, cores should be read in as SMARTS.  When read as SMILES\n\
   dummy atoms only match other dummy atoms.\n\
   The replaceDummies flag indicates whether matches to the dummy atoms should be considered as part\n\
   of the core or as part of the decoration (r-group)\n\n\
   >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmarts('[*]N[*]'),\n\
   ...                         replaceDummies=True),\n\
   ...             isomericSmiles=True)\n\
   '[1*]CC[2*]'\n\
   >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmarts('[*]N[*]'),\n\
   ...                         replaceDummies=False),\n\
   ...             isomericSmiles=True)\n\
   '[1*]CCCC[2*]'\n\
\n\
\n\
   >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CCC1CN'),MolFromSmarts('C1CCC1[*]'),\n\
   ...                         replaceDummies=False),\n\
   ...             isomericSmiles=True)\n\
   '[1*]CN'\n\
\n\
\n";
    m.def("ReplaceCore",
          (ROMol * (*)(const ROMol &, const ROMol &, bool, bool, bool, bool))
              replaceCore,
          "mol"_a, "coreQuery"_a, "replaceDummies"_a = true,
          "labelByIndex"_a = false, "requireDummyMatch"_a = false,
          "useChirality"_a = false, docString.c_str(),
          nb::rv_policy::take_ownership);

    docString = R"DOC(Return a new molecule with all BRICS bonds broken)DOC";
    m.def("FragmentOnBRICSBonds", MolFragmenter::fragmentOnBRICSBonds, "mol"_a,
          docString.c_str(), nb::rv_policy::take_ownership);

    docString =
        "Return a new molecule with all specified bonds broken\n\
\n\
  ARGUMENTS:\n\
\n\
      - mol            - the molecule to be modified\n\
      - bondIndices    - indices of the bonds to be broken\n\
      - addDummies  - toggles addition of dummy atoms to indicate where \n\
        bonds were broken\n\
      - dummyLabels - used to provide the labels to be used for the dummies.\n\
        the first element in each pair is the label for the dummy\n\
        that replaces the bond's beginAtom, the second is for the \n\
        dummy that replaces the bond's endAtom. If not provided, the\n\
        dummies are labeled with atom indices.\n\
      - bondTypes - used to provide the bond type to use between the\n\
        fragments and the dummy atoms. If not provided, defaults to single. \n\
      - cutsPerAtom - used to return the number of cuts made at each atom. \n\
\n\
  RETURNS:\n\
      a new Mol with the modifications\n\
";
    m.def("FragmentOnBonds", fragmentOnBondsHelper, "mol"_a, "bondIndices"_a,
          "addDummies"_a = true, "dummyLabels"_a = nb::none(),
          "bondTypes"_a = nb::none(), "cutsPerAtom"_a = nb::list(),
          docString.c_str(), nb::rv_policy::take_ownership);
    docString = R"DOC(fragment on some bonds)DOC";
    m.def("FragmentOnSomeBonds", fragmentOnSomeBondsHelper, "mol"_a,
          "bondIndices"_a, "numToBreak"_a = 1, "addDummies"_a = true,
          "dummyLabels"_a = nb::none(), "bondTypes"_a = nb::none(),
          "returnCutsPerAtom"_a = false, docString.c_str());

    nb::enum_<MolzipLabel>(m, "MolzipLabel")
        .value("AtomMapNumber", MolzipLabel::AtomMapNumber)
        .value("Isotope", MolzipLabel::Isotope)
        .value("FragmentOnBonds", MolzipLabel::FragmentOnBonds)
        .value("AtomType", MolzipLabel::AtomType);

    docString =
        "Parameters controlling how to zip molecules together\n\
\n\
  OPTIONS:\n\
      label : set the MolzipLabel option [default MolzipLabel.AtomMapNumber]\n\
\n\
  MolzipLabel.AtomMapNumber: atom maps are on dummy atoms, zip together the corresponding\n\
     attached atoms, i.e.  zip 'C[*:1]' 'N[*:1]' results in 'CN'\n\
\n\
  MolzipLabel.Isotope: isotope labels are on dummy atoms, zip together the corresponding\n\
     attached atoms, i.e.  zip 'C[1*]' 'N[1*]' results in 'CN'\n\
\n\
  MolzipLabel.FragmentOnBonds: zip together molecules generated by fragment on bonds.\n\
    Note the atom indices cannot change or be reordered from the output of fragmentOnBonds\n\
\n\
  MolzipLabel.AtomTypes: choose the atom types to act as matching dummy atoms.\n\
    i.e.  'C[V]' and 'N[Xe]' with atoms pairs [('V', 'Xe')] results in 'CN'\n\
";

    nb::class_<MolzipParams>(m, "MolzipParams", docString.c_str())
        .def(nb::init<>())
        .def_rw("label", &MolzipParams::label,
                "Set the atom labeling system to zip together")
        .def_rw("enforceValenceRules", &MolzipParams::enforceValenceRules,
                "If true (default) enforce valences after zipping\n\
Setting this to false allows assembling chemically incorrect fragments.")
        .def_rw(
            "generateCoordinates", &MolzipParams::generateCoordinates,
            "If true will add depiction coordinates to input molecules and\n\
zipped molecule (for molzipFragments only)")
        .def_rw(
            "alignCoordinates", &MolzipParams::alignCoordinates,
            "if true and the input fragments have coordinates, the fragments\n\
will be aligned along connection vectors in the output molecule")
        .def("setAtomSymbols", &RDKit::setAtomSymbols, "symbols"_a,
             "Set the atom symbols used to zip mols together when using "
             "AtomType labeling")
        .def("__setattr__", &safeSetattr);

    docString =
        "molzip: zip molecules together preserving bond and atom stereochemistry.\n\
\n\
This is useful when dealing with results from fragmentOnBonds, RGroupDecomposition and MMPs.\n\
\n\
Example:\n\
    >>> from rdkit.Chem import MolFromSmiles,  MolToSmiles, molzip\n\
    >>> a = MolFromSmiles('C=C[*:1]')\n\
    >>> b = MolFromSmiles('O/C=N/[*:1]')\n\
    >>> c = molzip(a,b)\n\
    >>> MolToSmiles(c)\n\
    'C=C/N=C/O'\n\
\n\
    >>> from rdkit.Chem import CombineMols, MolFromSmiles,  MolToSmiles, molzip\n\
    >>> a = MolFromSmiles('C=C[*:1]')\n\
    >>> b = MolFromSmiles('O/C=N/[*:1]')\n\
    >>> combined = CombineMols(a, b)\n\
    >>> c = molzip(combined)\n\
    >>> MolToSmiles(c)\n\
    'C=C/N=C/O'\n\
\n\
The atoms to zip can be specified with the MolzipParams class.\n\
    >>> from rdkit.Chem import MolzipParams, MolzipLabel\n\
    >>> a = MolFromSmiles('C=C[1*]')\n\
    >>> b = MolFromSmiles('O/C=N/[1*]')\n\
    >>> p = MolzipParams()\n\
    >>> p.label = MolzipLabel.Isotope\n\
    >>> c = molzip(a,b, p)\n\
    >>> MolToSmiles(c)\n\
    'C=C/N=C/O'\n\
";
    m.def("molzip",
          (ROMol * (*)(const ROMol &, const ROMol &, const MolzipParams &)) &
              molzip_new,
          "a"_a, "b"_a, "params"_a = MolzipParams(),
          "zip together two molecules using the given matching parameters",
          nb::rv_policy::take_ownership);

    m.def(
        "molzip",
        (ROMol * (*)(const ROMol &, const MolzipParams &)) & molzip_new, "a"_a,
        "params"_a = MolzipParams(),
        "zip together multiple molecules within a combined molecule using the given matching parameters",
        nb::rv_policy::take_ownership);

    m.def("molzipFragments",
          (ROMol * (*)(nb::object &, const MolzipParams &)) & molzipHelper,
          "mols"_a, "params"_a = MolzipParams(),
          "zip together multiple molecules from an R group decomposition \n\
using the given matching parameters.  The first molecule in the list\n\
must be the core",
          nb::rv_policy::take_ownership);

    docString =
        "zip an RGroupRow together to recreate the original molecule.  This correctly handles\n"
        "broken cycles that can occur in decompositions.\n"
        " example:\n\n"
        "  >>> from rdkit import Chem\n"
        "  >>> from rdkit.Chem import rdRGroupDecomposition as rgd\n"
        "  >>> core = Chem.MolFromSmiles('CO')\n"
        "  >>> mols = [Chem.MolFromSmiles('C1NNO1')]\n"
        "  >>> rgroups, unmatched = rgd.RGroupDecompose(core, mols)\n"
        "  >>> for rgroup in rgroups:\n"
        "  ...     mol = rgd.molzip(rgroup)\n"
        "\n";
    m.def("molzip",
          (ROMol * (*)(nb::dict, const MolzipParams &)) & rgroupRowZipHelper,
          "row"_a, "params"_a = MolzipParams(), docString.c_str(),
          nb::rv_policy::take_ownership);

    // ------------------------------------------------------------------------
    docString =
        "Adds a recursive query to an atom\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - query: the molecule to be used as the recursive query (this will be copied)\n\
\n\
    - atomIdx: the atom to modify\n\
\n\
    - preserveExistingQuery: (optional) if this is set, existing query information on the atom will be preserved\n\
\n\
  RETURNS: None\n\
\n";
    m.def("AddRecursiveQuery", addRecursiveQuery, "mol"_a, "query"_a,
          "atomIdx"_a, "preserveExistingQuery"_a = true, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Returns a copy of a molecule with renumbered atoms\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to be modified\n\
\n\
    - newOrder: the new ordering the atoms (should be numAtoms long)\n\
      for example: if newOrder is [3,2,0,1], then atom 3 in the original \n\
      molecule will be atom 0 in the new one\n\
\n\
\n";
    m.def("RenumberAtoms", renumberAtomsHelper, "mol"_a, "newOrder"_a,
          docString.c_str(), nb::rv_policy::take_ownership);

    docString =
        R"DOC(Possible values:
  - ADJUST_IGNORENONE: nothing will be ignored
  - ADJUST_IGNORECHAINS: non-ring atoms/bonds will be ignored
  - ADJUST_IGNORERINGS: ring atoms/bonds will be ignored
  - ADJUST_IGNOREDUMMIES: dummy atoms will be ignored
  - ADJUST_IGNORENONDUMMIES: non-dummy atoms will be ignored
  - ADJUST_IGNOREMAPPED: mapped atoms will be ignored
  - ADJUST_IGNOREALL: everything will be ignored
)DOC";
    nb::enum_<MolOps::AdjustQueryWhichFlags>(m, "AdjustQueryWhichFlags")
        .value("ADJUST_IGNORENONE", MolOps::ADJUST_IGNORENONE)
        .value("ADJUST_IGNORECHAINS", MolOps::ADJUST_IGNORECHAINS)
        .value("ADJUST_IGNORERINGS", MolOps::ADJUST_IGNORERINGS)
        .value("ADJUST_IGNOREDUMMIES", MolOps::ADJUST_IGNOREDUMMIES)
        .value("ADJUST_IGNORENONDUMMIES", MolOps::ADJUST_IGNORENONDUMMIES)
        .value("ADJUST_IGNOREMAPPED", MolOps::ADJUST_IGNOREMAPPED)
        .value("ADJUST_IGNOREALL", MolOps::ADJUST_IGNOREALL)
        .export_values();

    docString =
        R"DOC(Parameters controlling which components of the query atoms/bonds are adjusted.

Note that some of the options here are either directly contradictory or make
  no sense when combined with each other. We generally assume that client code
  is doing something sensible and don't attempt to detect possible conflicts or
  problems.

A note on the flags controlling which atoms/bonds are modified: 
   These generally limit the set of atoms/bonds to be modified.
   For example:
       - ADJUST_IGNORERINGS atoms/bonds in rings will not be modified.
       - ADJUST_IGNORENONE causes all atoms/bonds to be modified
       - ADJUST_IGNOREALL no atoms/bonds will be modified
   Some of the options obviously make no sense for bonds
)DOC";
    nb::class_<MolOps::AdjustQueryParameters>(m, "AdjustQueryParameters",
                                              docString.c_str())
        .def_rw("adjustDegree", &MolOps::AdjustQueryParameters::adjustDegree,
                "add degree queries")
        .def_rw("adjustDegreeFlags",
                &MolOps::AdjustQueryParameters::adjustDegreeFlags,
                "controls which atoms have their degree queries changed")
        .def_rw("adjustHeavyDegree",
                &MolOps::AdjustQueryParameters::adjustHeavyDegree,
                "adjust the heavy-atom degree")
        .def_rw("adjustHeavyDegreeFlags",
                &MolOps::AdjustQueryParameters::adjustHeavyDegreeFlags,
                "controls which atoms have their heavy-atom degree "
                "queries changed")
        .def_rw("adjustRingCount",
                &MolOps::AdjustQueryParameters::adjustRingCount,
                "add ring-count queries")
        .def_rw("adjustRingCountFlags",
                &MolOps::AdjustQueryParameters::adjustRingCountFlags,
                "controls which atoms have ring-count queries added")
        .def_rw(
            "makeDummiesQueries",
            &MolOps::AdjustQueryParameters::makeDummiesQueries,
            "convert dummy atoms without isotope labels to any-atom queries")
        .def_rw("aromatizeIfPossible",
                &MolOps::AdjustQueryParameters::aromatizeIfPossible,
                "perceive and set aromaticity")
        .def_rw("makeBondsGeneric",
                &MolOps::AdjustQueryParameters::makeBondsGeneric,
                "converts bonds to generic queries (any bonds)")
        .def_rw("makeBondsGenericFlags",
                &MolOps::AdjustQueryParameters::makeBondsGenericFlags,
                "controls which bonds are converted to generic queries")
        .def_rw("makeAtomsGeneric",
                &MolOps::AdjustQueryParameters::makeAtomsGeneric,
                "convert atoms to generic queries (any atoms)")
        .def_rw("makeAtomsGenericFlags",
                &MolOps::AdjustQueryParameters::makeAtomsGenericFlags,
                "controls which atoms are converted to generic queries")
        .def_rw("adjustRingChain",
                &MolOps::AdjustQueryParameters::adjustRingChain,
                "add ring-chain queries to atoms")
        .def_rw("adjustRingChainFlags",
                &MolOps::AdjustQueryParameters::adjustRingChainFlags,
                "controls which atoms have ring-chain queries added")
        .def_rw(
            "useStereoCareForBonds",
            &MolOps::AdjustQueryParameters::useStereoCareForBonds,
            "if this is set sterochemistry information will be removed from "
            "double bonds that do not have the stereoCare property set")
        .def_rw("adjustConjugatedFiveRings",
                &MolOps::AdjustQueryParameters::adjustConjugatedFiveRings,
                "set bond queries in conjugated five-rings to "
                "SINGLE|DOUBLE|AROMATIC")
        .def_rw("setMDLFiveRingAromaticity",
                &MolOps::AdjustQueryParameters::setMDLFiveRingAromaticity,
                "uses the 5-ring aromaticity behavior of the (former) MDL "
                "software "
                "as documented in the Chemical Representation Guide")
        .def_rw("adjustSingleBondsToDegreeOneNeighbors",
                &MolOps::AdjustQueryParameters::
                    adjustSingleBondsToDegreeOneNeighbors,
                "set single bonds bewteen aromatic or conjugated atoms "
                "and degree-one neighbors to SINGLE|AROMATIC")
        .def_rw("adjustSingleBondsBetweenAromaticAtoms",
                &MolOps::AdjustQueryParameters::
                    adjustSingleBondsBetweenAromaticAtoms,
                "sets non-ring single bonds between two aromatic or "
                "conjugated atoms to SINGLE|AROMATIC")
        .def_static(
            "NoAdjustments", &MolOps::AdjustQueryParameters::noAdjustments,
            "Returns an AdjustQueryParameters object with all parameters set "
            "to false");

    docString =
        "Returns a new molecule where the query properties of atoms have "
        "been modified.";
    m.def("AdjustQueryProperties", adjustQueryPropertiesHelper, "mol"_a,
          "params"_a = nb::none(), docString.c_str(),
          nb::rv_policy::take_ownership);
    docString =
        "Returns a new molecule where the query properties of atoms have "
        "been modified and generic group queries have been prepared.";
    m.def("AdjustQueryPropertiesWithGenericGroups",
          adjustQueryPropertiesWithGenericGroupsHelper, "mol"_a,
          "params"_a = nb::none(), docString.c_str(),
          nb::rv_policy::take_ownership);
    docString = R"DOC(checks for chemistry problems)DOC";
    m.def("DetectChemistryProblems", detectChemistryProblemsHelper, "mol"_a,
          "sanitizeOps"_a = MolOps::SANITIZE_ALL, docString.c_str());
    m.def("SetGenericQueriesFromProperties",
          GenericGroups::setGenericQueriesFromProperties, "mol"_a,
          "useAtomLabels"_a = true, "useSGroups"_a = true, "documentation");
    m.def("ConvertGenericQueriesToSubstanceGroups",
          GenericGroups::convertGenericQueriesToSubstanceGroups, "mol"_a,
          "documentation");
    m.def(
        "SetAllowNontetrahedralChirality",
        Chirality::setAllowNontetrahedralChirality, "val"_a,
        "toggles recognition of non-tetrahedral chirality from 3D structures");
    m.def("GetAllowNontetrahedralChirality",
          Chirality::getAllowNontetrahedralChirality,
          "returns whether or not recognition of non-tetrahedral "
          "chirality from 3D structures is enabled");
    m.def("SetUseLegacyStereoPerception",
          Chirality::setUseLegacyStereoPerception, "val"_a,
          "sets usage of the legacy stereo perception code");
    m.def("GetUseLegacyStereoPerception",
          Chirality::getUseLegacyStereoPerception,
          "returns whether or not the legacy stereo perception code is "
          "being used");
    m.def(
        "TranslateChiralFlagToStereoGroups", translateChiralFlagToStereoGroups,
        "mol"_a, "zeroFlagGroupType"_a = StereoGroupType::STEREO_AND,
        R"DOC(Generate enhanced stereo groups based on the status of the chiral flag property.
            
            Arguments:
            - mol: molecule to be modified
            - zeroFlagGroupType: how to handle non-grouped stereo centers when the
            chiral flag is set to zero
            
            If the chiral flag is set to a value of 1 then all specified tetrahedral
            chiral centers which are not already in StereoGroups will be added to an
            ABS StereoGroup.
            
            If the chiral flag is set to a value of 0 then all specified tetrahedral
            chiral centers will be added to a StereoGroup of the type zeroFlagGroupType
            
            If there is no chiral flag set (i.e. the property is not present), the
            molecule will not be modified.)DOC");
#if 1

    m.def(
        "ExpandAttachmentPoints", expandAttachmentPointsHelper, "mol"_a,
        "addAsQueries"_a = true, "addCoords"_a = true,
        R"DOC(attachment points encoded as attachPt properties are added to the graph as dummy atoms

  Arguments:
   - mol: molecule to be modified
   - addAsQueries: if true, the dummy atoms will be added as null queries
        (i.e. they will match any atom in a substructure search)
   - addCoords: if true and the molecule has one or more conformers, 
        positions for the attachment points will be added to the conformer(s)
)DOC");

    m.def(
        "CollapseAttachmentPoints", collapseAttachmentPointsHelper, "mol"_a,
        "markedOnly"_a = true,
        R"DOC(dummy atoms in the graph are removed and replaced with attachment point annotations on the attached atoms

  Arguments:
   - mol: molecule to be modified
   - markedOnly: if true, only dummy atoms with the _fromAttachPoint
     property will be collapsed

  In order for a dummy atom to be considered for collapsing it must have:
   - degree 1 with a single or unspecified bond
   - the bond to it can not be wedged
   - either no query or be an AtomNullQuery
)DOC");
    m.def("AddStereoAnnotations", Chirality::addStereoAnnotations, "mol"_a,
          "absLabel"_a = "abs ({cip})", "orLabel"_a = "or{id}",
          "andLabel"_a = "and{id}", "cipLabel"_a = "({cip})",
          "bondLabel"_a = "({cip})",
          R"DOC(add R/S, relative stereo, and E/Z annotations to atoms and bonds

  Arguments:
   - mol: molecule to modify
   - absLabel: label for atoms in an ABS stereo group
   - orLabel: label for atoms in an OR stereo group
   - andLabel: label for atoms in an AND stereo group
   - cipLabel: label for chiral atoms that aren't in a stereo group.
   - bondLabel: label for CIP stereochemistry on bonds

 If any label is empty, the corresponding annotations will not be added.

 The labels can contain the following placeholders:
   - {id} - the stereo group's index
   - {cip} - the atom or bond's CIP stereochemistry

 Note that CIP labels will only be added if CIP stereochemistry has been
 assigned to the molecule.
)DOC");

    m.def(
        "SimplifyEnhancedStereo", Chirality::simplifyEnhancedStereo, "mol"_a,
        "removeAffectedStereoGroups"_a = true,
        R"DOC(Simplifies the stereochemical representation of a molecule where all
specified stereocenters are in the same StereoGroup

  Arguments:
   - mol: molecule to modify
   - removeAffectedStereoGroups: if set then the affected StereoGroups will be removed

If all specified stereocenters are in the same AND or OR stereogroup, a
moleculeNote property will be set on the molecule with the value "AND
enantiomer" or "OR enantiomer". CIP labels, if present, are removed.
)DOC");

    m.def("_TestSetProps", testSetProps, "mol"_a);
    m.def("NeedsHs", MolOps::needsHs, "mol"_a,
          "returns whether or not the molecule needs to have Hs added");
    m.def(
        "CountAtomElec", MolOps::countAtomElec, "atom"_a,
        "returns the number of electrons available on an atom to donate for aromaticity");
    m.def("AtomHasConjugatedBond", MolOps::atomHasConjugatedBond, "atom"_a,
          "returns whether or not the atom is involved in a conjugated bond");

    nb::enum_<RDKit::SubsetMethod>(m, "SubsetMethod")
        .value("BONDS_BETWEEN_ATOMS", RDKit::SubsetMethod::BONDS_BETWEEN_ATOMS)
        .value("BONDS", RDKit::SubsetMethod::BONDS);

    nb::bind_vector<std::vector<bool>>(m, "BoolVector");

    nb::class_<RDKit::SubsetOptions>(m, "SubsetOptions")
        .def_rw("sanitize", &RDKit::SubsetOptions::sanitize,
                "Sanitize the resulting subset")
        .def_rw("clearComputedProps", &RDKit::SubsetOptions::clearComputedProps,
                "clear all computed props on the subsetted molecule")
        .def_rw("copyAsQuery", &RDKit::SubsetOptions::copyAsQuery,
                "Return the subset as a query")
        .def_rw("copyCoordinates", &RDKit::SubsetOptions::copyCoordinates,
                "Copy the active coordinates from the molecule")
        .def_rw("conformerIdx", &RDKit::SubsetOptions::conformerIdx,
                "What conformer idx to use for the coordinates default is -1")
        .def_rw("method", &RDKit::SubsetOptions::method,
                "Subsetting method to use");

    nb::bind_map<std::map<unsigned int, unsigned int>>(m, "UIntUIntMap");

    nb::class_<RDKit::SubsetInfo>(m, "SubsetInfo")
        .def_rw("atomMapping", &RDKit::SubsetInfo::atomMapping,
                "mapping from the original atom index to the subset atom index")
        .def_rw(
            "bondMapping", &RDKit::SubsetInfo::bondMapping,
            " mapping from the original bond index to the subset bond index");

    docString =
        "Extract a subgraph from an ROMol. Bonds, atoms, substance groups and \n\
stereo groups are only extracted to the subgraph if all participant entities \n\
are contained within the given atoms and bonds. \n\
\n\
ARGUMENTS:\n\
 - mol - starting mol \n\
 - atoms - indices atoms to extract \n\
 - bonds - indices bonds to extract \n\
 - subsetInfo - optional subsetInfo to record the atoms and bonds used \n\
 - options - optional subset options, note the method is ignored since all the atoms and bonds are sp \n\
\n";

    m.def("CopyMolSubset", copyMolSubsetHelper1, "mol"_a, "atomIndices"_a,
          "bondIndices"_a, "options"_a = SubsetOptions(), docString.c_str(),
          nb::rv_policy::take_ownership);
    m.def("CopyMolSubset", copyMolSubsetHelper2, "mol"_a, "atomIndices"_a,
          "bondIndices"_a, "subsetInfo"_a, "options"_a = SubsetOptions(),
          docString.c_str(), nb::rv_policy::take_ownership);

    docString =
        "Extract a subgraph from an ROMol. Bonds, atoms, substance groups and \n\
stereo groups are only extracted to the subgraph if all participant entities \n\
are contained within the given atoms and bonds. \n\
\n\
ARGUMENTS:\n\
 - mol - starting mol \n\
 - path - the indices of atoms or bonds to extract. If an index falls \n\
          outside of the acceptable indices, it is ignored.yes \n\
          Use SubsetMethod.BONDS to indicate a bond path and BONDS_BETWEEN_ATOMS \n\
          to indicate an atom path with any bond that includes atoms in the path. \n\
 - subsetInfo - optional subsetInfo to record the atoms and bonds used \n\
 - options - optional subset options, note the method is ignored since all the atoms and bonds are sp \n\
\n";

    m.def("CopyMolSubset", copyMolSubsetHelper3, "mol"_a, "path"_a,
          "options"_a = SubsetOptions(), "copy a subset of a molecule",
          nb::rv_policy::take_ownership);
    m.def("CopyMolSubset", copyMolSubsetHelper4, "mol"_a, "path"_a,
          "subsetInfo"_a, "options"_a = SubsetOptions(),
          "copy a subset of a molecule", nb::rv_policy::take_ownership);
#endif
  }
};

}  // namespace RDKit

void wrap_molops(nb::module_ &m) { RDKit::molops_wrapper::wrap(m); }
