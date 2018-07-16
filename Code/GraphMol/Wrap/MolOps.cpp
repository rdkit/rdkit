//
//  Copyright (C) 2003-2014 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define NO_IMPORT_ARRAY
#include "rdmolops.h"
#include <RDBoost/python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <string>
#include <math.h>

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <RDBoost/Wrap.h>
#include <RDBoost/python_streambuf.h>

#include <sstream>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
namespace python = boost::python;
using boost_adaptbx::python::streambuf;

namespace RDKit {
std::string molToSVG(const ROMol &mol, unsigned int width, unsigned int height,
                     python::object pyHighlightAtoms, bool kekulize,
                     unsigned int lineWidthMult, unsigned int fontSize,
                     bool includeAtomCircles, int confId) {
  RDUNUSED_PARAM(kekulize);
  std::unique_ptr<std::vector<int>> highlightAtoms =
      pythonObjectToVect(pyHighlightAtoms, static_cast<int>(mol.getNumAtoms()));
  std::stringstream outs;
  MolDraw2DSVG drawer(width, height, outs);
  drawer.setFontSize(fontSize / 24.);
  drawer.setLineWidth(drawer.lineWidth() * lineWidthMult);
  drawer.drawOptions().circleAtoms = includeAtomCircles;
  drawer.drawMolecule(mol, highlightAtoms.get(), nullptr, nullptr, confId);
  drawer.finishDrawing();
  return outs.str();
}
python::tuple fragmentOnSomeBondsHelper(const ROMol &mol,
                                        python::object pyBondIndices,
                                        unsigned int nToBreak, bool addDummies,
                                        python::object pyDummyLabels,
                                        python::object pyBondTypes,
                                        bool returnCutsPerAtom) {
  std::unique_ptr<std::vector<unsigned int>> bondIndices =
      pythonObjectToVect(pyBondIndices, mol.getNumBonds());
  if (!bondIndices.get()) throw_value_error("empty bond indices");

  std::vector<std::pair<unsigned int, unsigned int>> *dummyLabels = nullptr;
  if (pyDummyLabels) {
    unsigned int nVs =
        python::extract<unsigned int>(pyDummyLabels.attr("__len__")());
    dummyLabels = new std::vector<std::pair<unsigned int, unsigned int>>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      unsigned int v1 = python::extract<unsigned int>(pyDummyLabels[i][0]);
      unsigned int v2 = python::extract<unsigned int>(pyDummyLabels[i][1]);
      (*dummyLabels)[i] = std::make_pair(v1, v2);
    }
  }
  std::vector<Bond::BondType> *bondTypes = nullptr;
  if (pyBondTypes) {
    unsigned int nVs =
        python::extract<unsigned int>(pyBondTypes.attr("__len__")());
    if (nVs != bondIndices->size()) {
      throw_value_error("bondTypes shorter than bondIndices");
    }
    bondTypes = new std::vector<Bond::BondType>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      (*bondTypes)[i] = python::extract<Bond::BondType>(pyBondTypes[i]);
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
  python::list res;
  for (auto &frag : frags) {
    res.append(frag);
  }
  delete dummyLabels;
  delete bondTypes;
  if (cutsPerAtom) {
    python::list pyCutsPerAtom;
    for (auto &cut : *cutsPerAtom) {
      python::list localL;
      for (unsigned int j = 0; j < mol.getNumAtoms(); ++j) {
        localL.append(cut[j]);
      }
      pyCutsPerAtom.append(python::tuple(localL));
    }
    delete cutsPerAtom;
    python::list tres;
    tres.append(python::tuple(res));
    tres.append(python::tuple(pyCutsPerAtom));
    return python::tuple(tres);
  } else {
    return python::tuple(res);
  }
}

python::tuple getShortestPathHelper(const ROMol &mol, int aid1, int aid2) {
  if (aid1 < 0 || aid1 >= rdcast<int>(mol.getNumAtoms()) || aid2 < 0 ||
      aid2 >= rdcast<int>(mol.getNumAtoms())) {
    throw_value_error("bad atom index");
  }
  return static_cast<python::tuple>(MolOps::getShortestPath(mol, aid1, aid2));
}

ROMol *fragmentOnBondsHelper(const ROMol &mol, python::object pyBondIndices,
                             bool addDummies, python::object pyDummyLabels,
                             python::object pyBondTypes,
                             python::list pyCutsPerAtom) {
  std::unique_ptr<std::vector<unsigned int>> bondIndices =
      pythonObjectToVect(pyBondIndices, mol.getNumBonds());
  if (!bondIndices.get()) throw_value_error("empty bond indices");
  std::vector<std::pair<unsigned int, unsigned int>> *dummyLabels = nullptr;
  if (pyDummyLabels) {
    unsigned int nVs =
        python::extract<unsigned int>(pyDummyLabels.attr("__len__")());
    dummyLabels = new std::vector<std::pair<unsigned int, unsigned int>>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      unsigned int v1 = python::extract<unsigned int>(pyDummyLabels[i][0]);
      unsigned int v2 = python::extract<unsigned int>(pyDummyLabels[i][1]);
      (*dummyLabels)[i] = std::make_pair(v1, v2);
    }
  }
  std::vector<Bond::BondType> *bondTypes = nullptr;
  if (pyBondTypes) {
    unsigned int nVs =
        python::extract<unsigned int>(pyBondTypes.attr("__len__")());
    if (nVs != bondIndices->size()) {
      throw_value_error("bondTypes shorter than bondIndices");
    }
    bondTypes = new std::vector<Bond::BondType>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      (*bondTypes)[i] = python::extract<Bond::BondType>(pyBondTypes[i]);
    }
  }
  std::vector<unsigned int> *cutsPerAtom = nullptr;
  if (pyCutsPerAtom) {
    cutsPerAtom = new std::vector<unsigned int>;
    unsigned int nAts =
        python::extract<unsigned int>(pyCutsPerAtom.attr("__len__")());
    if (nAts < mol.getNumAtoms()) {
      throw_value_error("cutsPerAtom shorter than the number of atoms");
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

ROMol *renumberAtomsHelper(const ROMol &mol, python::object &pyNewOrder) {
  if (python::extract<unsigned int>(pyNewOrder.attr("__len__")()) <
      mol.getNumAtoms()) {
    throw_value_error("atomCounts shorter than the number of atoms");
  }
  std::unique_ptr<std::vector<unsigned int>> newOrder =
      pythonObjectToVect(pyNewOrder, mol.getNumAtoms());
  ROMol *res = MolOps::renumberAtoms(mol, *newOrder);
  return res;
}

namespace {
std::string getResidue(const ROMol &m, const Atom *at) {
  RDUNUSED_PARAM(m);
  if (at->getMonomerInfo()->getMonomerType() != AtomMonomerInfo::PDBRESIDUE)
    return "";
  return static_cast<const AtomPDBResidueInfo *>(at->getMonomerInfo())
      ->getResidueName();
}
std::string getChainId(const ROMol &m, const Atom *at) {
  RDUNUSED_PARAM(m);
  if (at->getMonomerInfo()->getMonomerType() != AtomMonomerInfo::PDBRESIDUE)
    return "";
  return static_cast<const AtomPDBResidueInfo *>(at->getMonomerInfo())
      ->getChainId();
}
}  // namespace
python::dict splitMolByPDBResidues(const ROMol &mol, python::object pyWhiteList,
                                   bool negateList) {
  std::vector<std::string> *whiteList = nullptr;
  if (pyWhiteList) {
    unsigned int nVs =
        python::extract<unsigned int>(pyWhiteList.attr("__len__")());
    whiteList = new std::vector<std::string>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      (*whiteList)[i] = python::extract<std::string>(pyWhiteList[i]);
    }
  }
  std::map<std::string, boost::shared_ptr<ROMol>> res =
      MolOps::getMolFragsWithQuery(mol, getResidue, false, whiteList,
                                   negateList);
  delete whiteList;

  python::dict pyres;
  for (std::map<std::string, boost::shared_ptr<ROMol>>::const_iterator iter =
           res.begin();
       iter != res.end(); ++iter) {
    pyres[iter->first] = iter->second;
  }
  return pyres;
}
python::dict splitMolByPDBChainId(const ROMol &mol, python::object pyWhiteList,
                                  bool negateList) {
  std::vector<std::string> *whiteList = nullptr;
  if (pyWhiteList) {
    unsigned int nVs =
        python::extract<unsigned int>(pyWhiteList.attr("__len__")());
    whiteList = new std::vector<std::string>(nVs);
    for (unsigned int i = 0; i < nVs; ++i) {
      (*whiteList)[i] = python::extract<std::string>(pyWhiteList[i]);
    }
  }
  std::map<std::string, boost::shared_ptr<ROMol>> res =
      MolOps::getMolFragsWithQuery(mol, getChainId, false, whiteList,
                                   negateList);
  delete whiteList;

  python::dict pyres;
  for (std::map<std::string, boost::shared_ptr<ROMol>>::const_iterator iter =
           res.begin();
       iter != res.end(); ++iter) {
    pyres[iter->first] = iter->second;
  }
  return pyres;
}

python::dict parseQueryDefFileHelper(python::object &input, bool standardize,
                                     std::string delimiter, std::string comment,
                                     unsigned int nameColumn,
                                     unsigned int smartsColumn) {
  python::extract<std::string> get_filename(input);
  std::map<std::string, ROMOL_SPTR> queryDefs;

  if (get_filename.check()) {
    parseQueryDefFile(get_filename(), queryDefs, standardize, delimiter,
                      comment, nameColumn, smartsColumn);
  } else {
    auto *sb = new streambuf(input);
    std::istream *istr = new streambuf::istream(*sb);
    parseQueryDefFile(istr, queryDefs, standardize, delimiter, comment,
                      nameColumn, smartsColumn);
    delete istr;
    delete sb;
  }

  python::dict res;
  for (std::map<std::string, ROMOL_SPTR>::const_iterator iter =
           queryDefs.begin();
       iter != queryDefs.end(); ++iter) {
    res[iter->first] = iter->second;
  }

  return res;
}

void addRecursiveQueriesHelper(ROMol &mol, python::dict replDict,
                               std::string propName) {
  std::map<std::string, ROMOL_SPTR> replacements;
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(replDict.keys().attr("__len__")());
       ++i) {
    ROMol *m = python::extract<ROMol *>(replDict.values()[i]);
    ROMOL_SPTR nm(new ROMol(*m));
    std::string k = python::extract<std::string>(replDict.keys()[i]);
    replacements[k] = nm;
  }
  addRecursiveQueries(mol, replacements, propName);
}

ROMol *addHs(const ROMol &orig, bool explicitOnly, bool addCoords,
             python::object onlyOnAtoms, bool addResidueInfo) {
  std::unique_ptr<std::vector<unsigned int>> onlyOn;
  if (onlyOnAtoms) {
    onlyOn = pythonObjectToVect(onlyOnAtoms, orig.getNumAtoms());
  }
  ROMol *res = MolOps::addHs(orig, explicitOnly, addCoords, onlyOn.get(),
                             addResidueInfo);
  return res;
}
int getSSSR(ROMol &mol) {
  VECT_INT_VECT rings;
  int nr = MolOps::findSSSR(mol, rings);
  return nr;
}

PyObject *replaceSubstructures(const ROMol &orig, const ROMol &query,
                               const ROMol &replacement,
                               bool replaceAll = false,
                               unsigned int replacementConnectionPoint = 0,
                               bool useChirality = false) {
  std::vector<ROMOL_SPTR> v =
      replaceSubstructs(orig, query, replacement, replaceAll,
                        replacementConnectionPoint, useChirality);
  PyObject *res = PyTuple_New(v.size());
  for (unsigned int i = 0; i < v.size(); ++i) {
    PyTuple_SetItem(res, i, python::converter::shared_ptr_to_python(v[i]));
  }
  return res;
}

void addRecursiveQuery(ROMol &mol, const ROMol &query, unsigned int atomIdx,
                       bool preserveExistingQuery) {
  if (atomIdx >= mol.getNumAtoms()) {
    throw_value_error("atom index exceeds mol.GetNumAtoms()");
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
MolOps::SanitizeFlags sanitizeMol(ROMol &mol, boost::uint64_t sanitizeOps,
                                  bool catchErrors) {
  RWMol &wmol = static_cast<RWMol &>(mol);
  unsigned int operationThatFailed;
  if (catchErrors) {
    try {
      MolOps::sanitizeMol(wmol, operationThatFailed, sanitizeOps);
    } catch (const MolSanitizeException &e) {
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
  RWMol *res = static_cast<RWMol *>(new ROMol(mol, false));
  return res;
}

ROMol *getNormal(const RWMol &mol) {
  ROMol *res = static_cast<ROMol *>(new RWMol(mol));
  return res;
}

void kekulizeMol(ROMol &mol, bool clearAromaticFlags = false) {
  RWMol &wmol = static_cast<RWMol &>(mol);
  MolOps::Kekulize(wmol, clearAromaticFlags);
}

void cleanupMol(ROMol &mol) {
  RWMol &rwmol = static_cast<RWMol &>(mol);
  MolOps::cleanUp(rwmol);
}

void setAromaticityMol(ROMol &mol, MolOps::AromaticityModel model) {
  RWMol &wmol = static_cast<RWMol &>(mol);
  MolOps::setAromaticity(wmol, model);
}

void setConjugationMol(ROMol &mol) {
  RWMol &wmol = static_cast<RWMol &>(mol);
  MolOps::setConjugation(wmol);
}

void assignRadicalsMol(ROMol &mol) {
  RWMol &wmol = static_cast<RWMol &>(mol);
  MolOps::assignRadicals(wmol);
}

void setHybridizationMol(ROMol &mol) {
  RWMol &wmol = static_cast<RWMol &>(mol);
  MolOps::setHybridization(wmol);
}

VECT_INT_VECT getSymmSSSR(ROMol &mol) {
  VECT_INT_VECT rings;
  MolOps::symmetrizeSSSR(mol, rings);
  return rings;
}
PyObject *getDistanceMatrix(ROMol &mol, bool useBO = false,
                            bool useAtomWts = false, bool force = false,
                            const char *prefix = nullptr) {
  int nats = mol.getNumAtoms();
  npy_intp dims[2];
  dims[0] = nats;
  dims[1] = nats;
  double *distMat;

  distMat = MolOps::getDistanceMat(mol, useBO, useAtomWts, force, prefix);

  PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);

  memcpy(PyArray_DATA(res), static_cast<void *>(distMat),
         nats * nats * sizeof(double));

  return PyArray_Return(res);
}
PyObject *get3DDistanceMatrix(ROMol &mol, int confId = -1,
                              bool useAtomWts = false, bool force = false,
                              const char *prefix = nullptr) {
  int nats = mol.getNumAtoms();
  npy_intp dims[2];
  dims[0] = nats;
  dims[1] = nats;
  double *distMat;

  distMat = MolOps::get3DDistanceMat(mol, confId, useAtomWts, force, prefix);

  PyArrayObject *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);

  memcpy(PyArray_DATA(res), static_cast<void *>(distMat),
         nats * nats * sizeof(double));

  if (prefix == nullptr || std::string(prefix) == "") {
    delete[] distMat;
  }
  return PyArray_Return(res);
}

PyObject *getAdjacencyMatrix(ROMol &mol, bool useBO = false, int emptyVal = 0,
                             bool force = false, const char *prefix = nullptr) {
  int nats = mol.getNumAtoms();
  npy_intp dims[2];
  dims[0] = nats;
  dims[1] = nats;

  double *tmpMat =
      MolOps::getAdjacencyMatrix(mol, useBO, emptyVal, force, prefix);

  PyArrayObject *res;
  if (useBO) {
    // if we're using valence, the results matrix is made up of doubles
    res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    memcpy(PyArray_DATA(res), static_cast<void *>(tmpMat),
           nats * nats * sizeof(double));
  } else {
    res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_INT);
    int *data = static_cast<int *>(PyArray_DATA(res));
    for (int i = 0; i < nats; i++) {
      for (int j = 0; j < nats; j++) {
        data[i * nats + j] = (int)round(tmpMat[i * nats + j]);
      }
    }
  }
  return PyArray_Return(res);
}

python::tuple GetMolFragsWithMapping(
    const ROMol &mol, bool asMols, bool sanitizeFrags,
    python::object frags = python::object(),
    python::object fragsMolAtomMapping = python::object()) {
  python::list res;

  if (!asMols) {
    VECT_INT_VECT fragsVec;
    MolOps::getMolFrags(mol, fragsVec);

    for (unsigned int i = 0; i < fragsVec.size(); ++i) {
      python::list tpl;
      for (unsigned int j = 0; j < fragsVec[i].size(); ++j) {
        tpl.append(fragsVec[i][j]);
      }
      res.append(python::tuple(tpl));
    }
  } else {
    std::vector<std::vector<int>> fragsMolAtomMappingVec;
    std::vector<int> fragsVec;
    std::vector<boost::shared_ptr<ROMol>> molFrags;
    python::list &fragsList = reinterpret_cast<python::list &>(frags);
    python::list &fragsMolAtomMappingList =
        reinterpret_cast<python::list &>(fragsMolAtomMapping);
    bool hasFrags = fragsList != python::object();
    bool hasFragsMolAtomMapping = fragsMolAtomMappingList != python::object();
    molFrags =
        hasFrags || hasFragsMolAtomMapping
            ? MolOps::getMolFrags(
                  mol, sanitizeFrags, hasFrags ? &fragsVec : NULL,
                  hasFragsMolAtomMapping ? &fragsMolAtomMappingVec : NULL)
            : MolOps::getMolFrags(mol, sanitizeFrags);
    if (hasFrags) {
      for (unsigned int i = 0; i < fragsVec.size(); ++i)
        fragsList.append(fragsVec[i]);
    }
    if (hasFragsMolAtomMapping) {
      for (unsigned int i = 0; i < fragsMolAtomMappingVec.size(); ++i) {
        python::list perFragMolAtomMappingTpl;
        for (unsigned int j = 0; j < fragsMolAtomMappingVec[i].size(); ++j)
          perFragMolAtomMappingTpl.append(fragsMolAtomMappingVec[i][j]);
        fragsMolAtomMappingList.append(python::tuple(perFragMolAtomMappingTpl));
      }
    }
    for (unsigned int i = 0; i < molFrags.size(); ++i) res.append(molFrags[i]);
  }
  return python::tuple(res);
}

python::tuple GetMolFrags(const ROMol &mol, bool asMols, bool sanitizeFrags) {
  return GetMolFragsWithMapping(mol, asMols, sanitizeFrags);
}

ExplicitBitVect *wrapLayeredFingerprint(
    const ROMol &mol, unsigned int layerFlags, unsigned int minPath,
    unsigned int maxPath, unsigned int fpSize, python::list atomCounts,
    ExplicitBitVect *includeOnlyBits, bool branchedPaths,
    python::object fromAtoms) {
  std::unique_ptr<std::vector<unsigned int>> lFromAtoms =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::vector<unsigned int> *atomCountsV = nullptr;
  if (atomCounts) {
    atomCountsV = new std::vector<unsigned int>;
    unsigned int nAts =
        python::extract<unsigned int>(atomCounts.attr("__len__")());
    if (nAts < mol.getNumAtoms()) {
      throw_value_error("atomCounts shorter than the number of atoms");
    }
    atomCountsV->resize(nAts);
    for (unsigned int i = 0; i < nAts; ++i) {
      (*atomCountsV)[i] = python::extract<unsigned int>(atomCounts[i]);
    }
  }

  ExplicitBitVect *res;
  res = RDKit::LayeredFingerprintMol(mol, layerFlags, minPath, maxPath, fpSize,
                                     atomCountsV, includeOnlyBits,
                                     branchedPaths, lFromAtoms.get());

  if (atomCountsV) {
    for (unsigned int i = 0; i < atomCountsV->size(); ++i) {
      atomCounts[i] = (*atomCountsV)[i];
    }
    delete atomCountsV;
  }

  return res;
}
ExplicitBitVect *wrapPatternFingerprint(const ROMol &mol, unsigned int fpSize,
                                        python::list atomCounts,
                                        ExplicitBitVect *includeOnlyBits) {
  std::vector<unsigned int> *atomCountsV = nullptr;
  if (atomCounts) {
    atomCountsV = new std::vector<unsigned int>;
    unsigned int nAts =
        python::extract<unsigned int>(atomCounts.attr("__len__")());
    if (nAts < mol.getNumAtoms()) {
      throw_value_error("atomCounts shorter than the number of atoms");
    }
    atomCountsV->resize(nAts);
    for (unsigned int i = 0; i < nAts; ++i) {
      (*atomCountsV)[i] = python::extract<unsigned int>(atomCounts[i]);
    }
  }

  ExplicitBitVect *res;
  res = RDKit::PatternFingerprintMol(mol, fpSize, atomCountsV, includeOnlyBits);

  if (atomCountsV) {
    for (unsigned int i = 0; i < atomCountsV->size(); ++i) {
      atomCounts[i] = (*atomCountsV)[i];
    }
    delete atomCountsV;
  }

  return res;
}

ExplicitBitVect *wrapRDKFingerprintMol(
    const ROMol &mol, unsigned int minPath, unsigned int maxPath,
    unsigned int fpSize, unsigned int nBitsPerHash, bool useHs,
    double tgtDensity, unsigned int minSize, bool branchedPaths,
    bool useBondOrder, python::object atomInvariants, python::object fromAtoms,
    python::object atomBits, python::object bitInfo) {
  std::unique_ptr<std::vector<unsigned int>> lAtomInvariants =
      pythonObjectToVect<unsigned int>(atomInvariants);
  std::unique_ptr<std::vector<unsigned int>> lFromAtoms =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::vector<std::vector<boost::uint32_t>> *lAtomBits = nullptr;
  std::map<boost::uint32_t, std::vector<std::vector<int>>> *lBitInfo = nullptr;
  // if(!(atomBits.is_none())){
  if (atomBits != python::object()) {
    lAtomBits =
        new std::vector<std::vector<boost::uint32_t>>(mol.getNumAtoms());
  }
  if (bitInfo != python::object()) {
    lBitInfo = new std::map<boost::uint32_t, std::vector<std::vector<int>>>;
  }
  ExplicitBitVect *res;
  res = RDKit::RDKFingerprintMol(mol, minPath, maxPath, fpSize, nBitsPerHash,
                                 useHs, tgtDensity, minSize, branchedPaths,
                                 useBondOrder, lAtomInvariants.get(),
                                 lFromAtoms.get(), lAtomBits, lBitInfo);

  if (lAtomBits) {
    python::list &pyl = static_cast<python::list &>(atomBits);
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      python::list tmp;
      BOOST_FOREACH (boost::uint32_t v, (*lAtomBits)[i]) { tmp.append(v); }
      pyl.append(tmp);
    }
    delete lAtomBits;
  }
  if (lBitInfo) {
    python::dict &pyd = static_cast<python::dict &>(bitInfo);
    for (auto &it : (*lBitInfo)) {
      python::list temp;
      std::vector<std::vector<int>>::iterator itset;
      for (itset = it.second.begin(); itset != it.second.end(); ++itset) {
        python::list temp2;
        for (int &i : *itset) {
          temp2.append(i);
        }
        temp.append(temp2);
      }
      if (!pyd.has_key(it.first)) {
        pyd[it.first] = temp;
      }
    }
    delete lBitInfo;
  }

  return res;
}

SparseIntVect<boost::uint64_t> *wrapUnfoldedRDKFingerprintMol(
    const ROMol &mol, unsigned int minPath, unsigned int maxPath, bool useHs,
    bool branchedPaths, bool useBondOrder, python::object atomInvariants,
    python::object fromAtoms, python::object atomBits, python::object bitInfo) {
  std::unique_ptr<std::vector<unsigned int>> lAtomInvariants =
      pythonObjectToVect<unsigned int>(atomInvariants);
  std::unique_ptr<std::vector<unsigned int>> lFromAtoms =
      pythonObjectToVect(fromAtoms, mol.getNumAtoms());
  std::vector<std::vector<boost::uint64_t>> *lAtomBits = nullptr;
  std::map<boost::uint64_t, std::vector<std::vector<int>>> *lBitInfo = nullptr;

  // if(!(atomBits.is_none())){
  if (atomBits != python::object()) {
    lAtomBits =
        new std::vector<std::vector<boost::uint64_t>>(mol.getNumAtoms());
  }
  if (bitInfo != python::object()) {
    lBitInfo = new std::map<boost::uint64_t, std::vector<std::vector<int>>>;
  }

  SparseIntVect<boost::uint64_t> *res;
  res = getUnfoldedRDKFingerprintMol(
      mol, minPath, maxPath, useHs, branchedPaths, useBondOrder,
      lAtomInvariants.get(), lFromAtoms.get(), lAtomBits, lBitInfo);

  if (lAtomBits) {
    python::list &pyl = static_cast<python::list &>(atomBits);
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      python::list tmp;
      BOOST_FOREACH (boost::uint64_t v, (*lAtomBits)[i]) { tmp.append(v); }
      pyl.append(tmp);
    }
    delete lAtomBits;
  }
  if (lBitInfo) {
    python::dict &pyd = static_cast<python::dict &>(bitInfo);
    for (auto &it : (*lBitInfo)) {
      python::list temp;
      std::vector<std::vector<int>>::iterator itset;
      for (itset = it.second.begin(); itset != it.second.end(); ++itset) {
        python::list temp2;
        for (int &i : *itset) {
          temp2.append(i);
        }
        temp.append(temp2);
      }
      if (!pyd.has_key(it.first)) {
        pyd[it.first] = temp;
      }
    }
    delete lBitInfo;
  }

  return res;
}

python::object findAllSubgraphsOfLengthsMtoNHelper(const ROMol &mol,
                                                   unsigned int lowerLen,
                                                   unsigned int upperLen,
                                                   bool useHs = false,
                                                   int rootedAtAtom = -1) {
  if (lowerLen > upperLen) {
    throw_value_error("lowerLen > upperLen");
  }

  INT_PATH_LIST_MAP oMap = findAllSubgraphsOfLengthsMtoN(
      mol, lowerLen, upperLen, useHs, rootedAtAtom);
  python::list res;
  for (unsigned int i = lowerLen; i <= upperLen; ++i) {
    python::list tmp;
    const PATH_LIST &pth = oMap[i];
    for (const auto &pthit : pth) {
      tmp.append(python::tuple(pthit));
    }
    res.append(tmp);
  }
  return python::tuple(res);
};

ROMol *pathToSubmolHelper(const ROMol &mol, python::object &path, bool useQuery,
                          python::object atomMap) {
  ROMol *result;
  PATH_TYPE pth;
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(path.attr("__len__")()); ++i) {
    pth.push_back(python::extract<unsigned int>(path[i]));
  }
  std::map<int, int> mapping;
  result = Subgraphs::pathToSubmol(mol, pth, useQuery, mapping);
  if (atomMap != python::object()) {
    // make sure the optional argument actually was a dictionary
    python::dict typecheck = python::extract<python::dict>(atomMap);
    atomMap.attr("clear")();
    for (std::map<int, int>::const_iterator mIt = mapping.begin();
         mIt != mapping.end(); ++mIt) {
      atomMap[mIt->first] = mIt->second;
    }
  }
  return result;
}

ROMol *adjustQueryPropertiesHelper(const ROMol &mol, python::object pyparams) {
  MolOps::AdjustQueryParameters params;
  if (pyparams != python::object()) {
    params = python::extract<MolOps::AdjustQueryParameters>(pyparams);
  }
  return MolOps::adjustQueryProperties(mol, &params);
}

ROMol *replaceCoreHelper(const ROMol &mol, const ROMol &core,
                         python::object match, bool replaceDummies,
                         bool labelByIndex, bool requireDummyMatch = false) {
  // convert input to MatchVect
  MatchVectType matchVect;

  unsigned int length = python::extract<unsigned int>(match.attr("__len__")());

  for (unsigned int i = 0; i < length; ++i) {
    int sz = 1;
    if (PyObject_HasAttrString(static_cast<python::object>(match[i]).ptr(),
                               "__len__")) {
      sz = python::extract<unsigned int>(match[i].attr("__len__")());
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
        v2 = python::extract<int>(match[i]);
        break;
      case 2:
        v1 = python::extract<int>(match[i][0]);
        v2 = python::extract<int>(match[i][1]);
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

struct molops_wrapper {
  static void wrap() {
    std::string docString;
    python::enum_<MolOps::SanitizeFlags>("SanitizeFlags")
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
        .value("SANITIZE_ADJUSTHS", MolOps::SANITIZE_ADJUSTHS)
        .value("SANITIZE_ALL", MolOps::SANITIZE_ALL)
        .export_values();
    ;

    // ------------------------------------------------------------------------
    docString =
        "Assign stereochemistry to bonds based on coordinates and a conformer.\n\
        DEPRECATED: Please use the version that takes a conformer ID instead.\n\
        \n\
  ARGUMENTS:\n\
  \n\
    - mol: the molecule to be modified\n\
    - conformer: Conformer providing the coordinates\n\
\n";
    python::def("DetectBondStereoChemistry", DetectBondStereoChemistry,
                (python::arg("mol"), python::arg("conformer")),
                docString.c_str());
    docString =
        "Assign stereochemistry to bonds based on coordinates.\n\
        \n\
  ARGUMENTS:\n\
  \n\
    - mol: the molecule to be modified\n\
    - confId: Conformer to use for the coordinates\n\
\n";
    python::def("DetectBondStereochemistry", MolOps::detectBondStereochemistry,
                (python::arg("mol"), python::arg("confId") = -1),
                docString.c_str());

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
    python::def(
        "SanitizeMol", sanitizeMol,
        (python::arg("mol"), python::arg("sanitizeOps") = MolOps::SANITIZE_ALL,
         python::arg("catchErrors") = false),
        docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Get the smallest set of simple rings for a molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use.\n\
\n\
  RETURNS: a sequence of sequences containing the rings found as atom ids\n\
         The length of this will be equal to NumBonds-NumAtoms+1 for single-fragment molecules.\n\
\n";
    python::def("GetSSSR", getSSSR, docString.c_str());

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
\n\
  RETURNS: a sequence of sequences containing the rings found as atom ids\n\
\n";
    python::def("GetSymmSSSR", getSymmSSSR, docString.c_str());

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
    python::def("FastFindRings", MolOps::fastFindRings, docString.c_str());

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
    python::def("AddHs", addHs,
                (python::arg("mol"), python::arg("explicitOnly") = false,
                 python::arg("addCoords") = false,
                 python::arg("onlyOnAtoms") = python::object(),
                 python::arg("addResidueInfo") = false),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());

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
\n";
    python::def("RemoveHs",
                (ROMol * (*)(const ROMol &, bool, bool, bool)) MolOps::removeHs,
                (python::arg("mol"), python::arg("implicitOnly") = false,
                 python::arg("updateExplicitCount") = false,
                 python::arg("sanitize") = true),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());

    python::def("MergeQueryHs",
                (ROMol * (*)(const ROMol &, bool)) & MolOps::mergeQueryHs,
                (python::arg("mol"), python::arg("mergeUnmappedOnly") = false),
                "merges hydrogens into their neighboring atoms as queries",
                python::return_value_policy<python::manage_new_object>());

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
    python::def(
        "DeleteSubstructs", deleteSubstructs,
        (python::arg("mol"), python::arg("query"),
         python::arg("onlyFrags") = false, python::arg("useChirality") = false),
        docString.c_str(),
        python::return_value_policy<python::manage_new_object>());
    docString = "Do a Murcko decomposition and return the scaffold";
    python::def("MurckoDecompose", MurckoDecompose, (python::arg("mol")),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());
    docString = "Combine the atoms from two molecules to produce a third";
    python::def("CombineMols", combineMols,
                (python::arg("mol1"), python::arg("mol2"),
                 python::arg("offset") = RDGeom::Point3D(0, 0, 0)),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());

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
\n\
  EXAMPLES:\n\
\n\
   The following examples substitute SMILES/SMARTS strings for molecules, you'd have\n\
   to actually use molecules:\n\
\n\
    - ReplaceSubstructs('CCOC','OC','NC') -> ('CCNC',)\n\
\n\
    - ReplaceSubstructs('COCCOC','OC','NC') -> ('COCCNC','CNCCOC')\n\
\n\
    - ReplaceSubstructs('COCCOC','OC','NC',True) -> ('CNCCNC',)\n\
\n\
    - ReplaceSubstructs('COCCOC','OC','CN',True,1) -> ('CNCCNC',)\n\
\n";
    python::def("ReplaceSubstructs", replaceSubstructures,
                (python::arg("mol"), python::arg("query"),
                 python::arg("replacement"), python::arg("replaceAll") = false,
                 python::arg("replacementConnectionPoint") = 0,
                 python::arg("useChirality") = false),
                docString.c_str());

    // ------------------------------------------------------------------------
    docString = "Adds named recursive queries to atoms\n";
    python::def(
        "MolAddRecursiveQueries", addRecursiveQueriesHelper,
        (python::arg("mol"), python::arg("queries"), python::arg("propName")),
        docString.c_str());

    docString = "reads query definitions from a simply formatted file\n";
    python::def(
        "ParseMolQueryDefFile", parseQueryDefFileHelper,
        (python::arg("fileobj"), python::arg("standardize") = true,
         python::arg("delimiter") = "\t", python::arg("comment") = "//",
         python::arg("nameColumn") = 0, python::arg("smartsColumn") = 1),
        docString.c_str());

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
    python::def("GetDistanceMatrix", getDistanceMatrix,
                (python::arg("mol"), python::arg("useBO") = false,
                 python::arg("useAtomWts") = false,
                 python::arg("force") = false, python::arg("prefix") = ""),
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
    python::def("Get3DDistanceMatrix", get3DDistanceMatrix,
                (python::arg("mol"), python::arg("confId") = -1,
                 python::arg("useAtomWts") = false,
                 python::arg("force") = false, python::arg("prefix") = ""),
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
    python::def("GetAdjacencyMatrix", getAdjacencyMatrix,
                (python::arg("mol"), python::arg("useBO") = false,
                 python::arg("emptyVal") = 0, python::arg("force") = false,
                 python::arg("prefix") = ""),
                docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Kekulizes the molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
    - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the \n\
      molecule will be marked non-aromatic following the kekulization.\n\
      Default value is 0.\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
    python::def("Kekulize", kekulizeMol,
                (python::arg("mol"), python::arg("clearAromaticFlags") = false),
                docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "cleans up certain common bad functionalities in the molecule\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n\
  NOTES:\n\
\n\
    - The molecule is modified in place.\n\
\n";
    python::def("Cleanup", cleanupMol, (python::arg("mol")), docString.c_str());

    python::enum_<MolOps::AromaticityModel>("AromaticityModel")
        .value("AROMATICITY_DEFAULT", MolOps::AROMATICITY_DEFAULT)
        .value("AROMATICITY_RDKIT", MolOps::AROMATICITY_RDKIT)
        .value("AROMATICITY_SIMPLE", MolOps::AROMATICITY_SIMPLE)
        .value("AROMATICITY_MDL", MolOps::AROMATICITY_MDL)
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
    python::def("SetAromaticity", setAromaticityMol,
                (python::arg("mol"),
                 python::arg("model") = MolOps::AROMATICITY_DEFAULT),
                docString.c_str());

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
    python::def("SetConjugation", setConjugationMol, (python::arg("mol")),
                docString.c_str());
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
    python::def("SetHybridization", setHybridizationMol, (python::arg("mol")),
                docString.c_str());
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
    python::def("AssignRadicals", assignRadicalsMol, (python::arg("mol")),
                docString.c_str());

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
    python::def(
        "FindAllSubgraphsOfLengthN", &findAllSubgraphsOfLengthN,
        (python::arg("mol"), python::arg("length"),
         python::arg("useHs") = false, python::arg("rootedAtAtom") = -1),
        docString.c_str());
    // ------------------------------------------------------------------------
    docString =
        "Finds all subgraphs of a particular length in a molecule\n\
  See documentation for FindAllSubgraphsOfLengthN for definitions\n\
\n";
    python::def(
        "FindAllSubgraphsOfLengthMToN", &findAllSubgraphsOfLengthsMtoNHelper,
        (python::arg("mol"), python::arg("min"), python::arg("max"),
         python::arg("useHs") = false, python::arg("rootedAtAtom") = -1),
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
    python::def("FindUniqueSubgraphsOfLengthN", &findUniqueSubgraphsOfLengthN,
                (python::arg("mol"), python::arg("length"),
                 python::arg("useHs") = false, python::arg("useBO") = true,
                 python::arg("rootedAtAtom") = -1),
                docString.c_str());

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
    python::def("FindAllPathsOfLengthN", &findAllPathsOfLengthN,
                (python::arg("mol"), python::arg("length"),
                 python::arg("useBonds") = true, python::arg("useHs") = false,
                 python::arg("rootedAtAtom") = -1),
                docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Finds the bonds within a certain radius of an atom in a molecule\n\
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
  RETURNS: a vector of bond IDs\n\
\n\
\n";
    python::def("FindAtomEnvironmentOfRadiusN", &findAtomEnvironmentOfRadiusN,
                (python::arg("mol"), python::arg("radius"),
                 python::arg("rootedAtAtom"), python::arg("useHs") = false),
                docString.c_str());

    python::def("PathToSubmol", pathToSubmolHelper,
                (python::arg("mol"), python::arg("path"),
                 python::arg("useQuery") = false,
                 python::arg("atomMap") = python::object()),
                "", python::return_value_policy<python::manage_new_object>());

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
    - frags: (optional, defaults to None) if this is provided as an empty list,\n\
      the result will be mol.GetNumAtoms() long on return and will contain the\n\
      fragment assignment for each Atom\n\
    - fragsMolAtomMapping: (optional, defaults to None) if this is provided as\n\
      an empty list, the result will be a a numFrags long list on return, and\n\
      each entry will contain the indices of the Atoms in that fragment:\n\
      [(0, 1, 2, 3), (4, 5)]\n\
\n\
  RETURNS: a tuple of tuples with IDs for the atoms in each fragment\n\
           or a tuple of molecules.\n\
\n";
    python::def("GetMolFrags", &GetMolFragsWithMapping,
                (python::arg("mol"), python::arg("asMols") = false,
                 python::arg("sanitizeFrags") = true,
                 python::arg("frags") = python::object(),
                 python::arg("fragsMolAtomMapping") = python::object()),
                docString.c_str());

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
    python::def(
        "SplitMolByPDBResidues", &splitMolByPDBResidues,
        (python::arg("mol"), python::arg("whiteList") = python::object(),
         python::arg("negateList") = false),
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
    python::def(
        "SplitMolByPDBChainId", &splitMolByPDBChainId,
        (python::arg("mol"), python::arg("whiteList") = python::object(),
         python::arg("negateList") = false),
        docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Returns the formal charge for the molecule.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
\n";
    python::def("GetFormalCharge", &MolOps::getFormalCharge, docString.c_str());

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
    python::def("GetShortestPath", getShortestPathHelper, docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Does the CIP stereochemistry assignment \n\
  for the molecule's atoms (R/S) and double bond (Z/E).\n\
  Chiral atoms will have a property '_CIPCode' indicating\n\
  their chiral code.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - cleanIt: (optional) if provided, atoms with a chiral specifier that aren't\n\
      actually chiral (e.g. atoms with duplicate substituents or only 2 substituents,\n\
      etc.) will have their chiral code set to CHI_UNSPECIFIED. Bonds with \n\
      STEREOCIS/STEREOTRANS specified that have duplicate substituents based upon the CIP \n\
      atom ranks will be marked STEREONONE. \n\
    - force: (optional) causes the calculation to be repeated, even if it has already\n\
      been done\n\
    - flagPossibleStereoCenters (optional)   set the _ChiralityPossible property on\n\
      atoms that are possible stereocenters\n\
\n";
    python::def("AssignStereochemistry", MolOps::assignStereochemistry,
                (python::arg("mol"), python::arg("cleanIt") = false,
                 python::arg("force") = false,
                 python::arg("flagPossibleStereoCenters") = false),
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
    python::def("FindPotentialStereoBonds", MolOps::findPotentialStereoBonds,
                (python::arg("mol"), python::arg("cleanIt") = false),
                docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Removes all stereochemistry info from the molecule.\n\
\n";
    python::def("RemoveStereochemistry", MolOps::removeStereochemistry,
                (python::arg("mol")), docString.c_str());

    // ------------------------------------------------------------------------
    docString =
        "Sets the chiral tags on a molecule's atoms based on \n\
  a 3D conformation.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule to use\n\
    - confId: the conformer id to use, -1 for the default \n\
    - replaceExistingTags: if True, existing stereochemistry information will be cleared \n\
                           before running the calculation. \n\
\n";
    python::def("AssignAtomChiralTagsFromStructure",
                MolOps::assignChiralTypesFrom3D,
                (python::arg("mol"), python::arg("confId") = -1,
                 python::arg("replaceExistingTags") = true),
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
    python::def(
        "RDKFingerprint", wrapRDKFingerprintMol,
        (python::arg("mol"), python::arg("minPath") = 1,
         python::arg("maxPath") = 7, python::arg("fpSize") = 2048,
         python::arg("nBitsPerHash") = 2, python::arg("useHs") = true,
         python::arg("tgtDensity") = 0.0, python::arg("minSize") = 128,
         python::arg("branchedPaths") = true,
         python::arg("useBondOrder") = true, python::arg("atomInvariants") = 0,
         python::arg("fromAtoms") = 0,
         python::arg("atomBits") = python::object(),
         python::arg("bitInfo") = python::object()),
        docString.c_str(),
        python::return_value_policy<python::manage_new_object>());
    python::scope().attr("_RDKFingerprint_version") =
        RDKit::RDKFingerprintMolVersion;

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

    python::def(
        "UnfoldedRDKFingerprintCountBased", wrapUnfoldedRDKFingerprintMol,
        (python::arg("mol"), python::arg("minPath") = 1,
         python::arg("maxPath") = 7, python::arg("useHs") = true,
         python::arg("branchedPaths") = true,
         python::arg("useBondOrder") = true, python::arg("atomInvariants") = 0,
         python::arg("fromAtoms") = 0,
         python::arg("atomBits") = python::object(),
         python::arg("bitInfo") = python::object()),
        docString.c_str(),
        python::return_value_policy<python::manage_new_object>());

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
    python::def(
        "LayeredFingerprint", wrapLayeredFingerprint,
        (python::arg("mol"), python::arg("layerFlags") = 0xFFFFFFFF,
         python::arg("minPath") = 1, python::arg("maxPath") = 7,
         python::arg("fpSize") = 2048,
         python::arg("atomCounts") = python::list(),
         python::arg("setOnlyBits") = (ExplicitBitVect *)nullptr,
         python::arg("branchedPaths") = true, python::arg("fromAtoms") = 0),
        docString.c_str(),
        python::return_value_policy<python::manage_new_object>());
    python::scope().attr("_LayeredFingerprint_version") =
        RDKit::LayeredFingerprintMolVersion;
    python::scope().attr("LayeredFingerprint_substructLayers") =
        RDKit::substructLayers;

    // ------------------------------------------------------------------------
    docString =
        "A fingerprint using SMARTS patterns \n\
\n\
  NOTE: This function is experimental. The API or results may change from\n\
    release to release.\n";
    python::def("PatternFingerprint", wrapPatternFingerprint,
                (python::arg("mol"), python::arg("fpSize") = 2048,
                 python::arg("atomCounts") = python::list(),
                 python::arg("setOnlyBits") = (ExplicitBitVect *)nullptr),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());

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
    python::def("WedgeMolBonds", WedgeMolBonds, docString.c_str());

    docString =
        "Set the wedging on an individual bond from a molecule.\n\
   The wedging scheme used is that from Mol files.\n\
\n\
  ARGUMENTS:\n\
\n\
    - bond: the bond to update\n\
    - atom ID: the atom from which to do the wedging\n\
    - conformer: the conformer to use to determine wedge direction\n\
\n\
\n";
    python::def("WedgeBond", WedgeBond, docString.c_str());

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
    python::def("ReplaceSidechains", replaceSidechains,
                (python::arg("mol"), python::arg("coreQuery"),
                 python::arg("useChirality") = false),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());

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
EXAMPLES:\n\
    >>> from rdkit.Chem import MolToSmiles, MolFromSmiles, ReplaceCore\n\
    >>> mol = MolFromSmiles('C1ONNCC1')\n\
    >>> core = MolFromSmiles('NN')\n\
\n\
    Note: Using isomericSmiles is necessary to see the labels.\n\
    >>> MolToSmiles(ReplaceCore(mol, core, mol.GetSubstructMatch(core)), isomericSmiles=True)\n\
    '[1*]OCCC[2*]'\n\
\n\
    Since NN is symmetric, we should actually get two matches here if we don't\n\
    uniquify the matches.\n\
    >>> [MolToSmiles(ReplaceCore(mol, core, match), isomericSmiles=True)\n\
    ...     for match in mol.GetSubstructMatches(core, uniquify=False)]\n\
    ['[1*]OCCC[2*]', '[1*]CCCO[2*]']\n\
\n\
";
    python::def("ReplaceCore", replaceCoreHelper,
                (python::arg("mol"), python::arg("core"),
                 python::arg("matches"), python::arg("replaceDummies") = true,
                 python::arg("labelByIndex") = false,
                 python::arg("requireDummyMatch") = false),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());
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
   in the core, like all good computer scientists expect, atoms indices start at 0.\n\
   >>> MolToSmiles(ReplaceCore(MolFromSmiles('CCN1CCC1'),MolFromSmiles('C1CCN1'),\n\
   ...                         labelByIndex=True),\n\
   ...   isomericSmiles=True)\n\
   '[3*]CC'\n\
\n\
   Non-core matches just return None\n\
   >>> ReplaceCore(MolFromSmiles('CCC1CC1'),MolFromSmiles('C1CCC1'))\n\
\n\
   The bond between atoms are considered part of the core and are removed as well\n\
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
   of the core or as part of the decoration (r-group)\n\
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
    python::def("ReplaceCore",
                (ROMol * (*)(const ROMol &, const ROMol &, bool, bool, bool,
                             bool)) replaceCore,
                (python::arg("mol"), python::arg("coreQuery"),
                 python::arg("replaceDummies") = true,
                 python::arg("labelByIndex") = false,
                 python::arg("requireDummyMatch") = false,
                 python::arg("useChirality") = false),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());

    docString = "Return a new molecule with all BRICS bonds broken";
    python::def("FragmentOnBRICSBonds", MolFragmenter::fragmentOnBRICSBonds,
                (python::arg("mol")), docString.c_str(),
                python::return_value_policy<python::manage_new_object>());

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
    python::def("FragmentOnBonds", fragmentOnBondsHelper,
                (python::arg("mol"), python::arg("bondIndices"),
                 python::arg("addDummies") = true,
                 python::arg("dummyLabels") = python::object(),
                 python::arg("bondTypes") = python::object(),
                 python::arg("cutsPerAtom") = python::list()),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());
    docString = "fragment on some bonds";
    python::def(
        "FragmentOnSomeBonds", fragmentOnSomeBondsHelper,
        (python::arg("mol"), python::arg("bondIndices"),
         python::arg("numToBreak") = 1, python::arg("addDummies") = true,
         python::arg("dummyLabels") = python::object(),
         python::arg("bondTypes") = python::object(),
         python::arg("returnCutsPerAtom") = false),
        docString.c_str());

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
    python::def(
        "AddRecursiveQuery", addRecursiveQuery,
        (python::arg("mol"), python::arg("query"), python::arg("atomIdx"),
         python::arg("preserveExistingQuery") = true),
        docString.c_str());

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
    python::def("RenumberAtoms", renumberAtomsHelper,
                (python::arg("mol"), python::arg("newOrder")),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());

    // ------------------------------------------------------------------------
    docString = "Returns svg for a molecule";
    python::def("MolToSVG", molToSVG,
                (python::arg("mol"), python::arg("width") = 300,
                 python::arg("height") = 300,
                 python::arg("highlightAtoms") = python::object(),
                 python::arg("kekulize") = true,
                 python::arg("lineWidthMult") = 1, python::arg("fontSize") = 12,
                 python::arg("includeAtomCircles") = true),
                docString.c_str());

    python::enum_<MolOps::AdjustQueryWhichFlags>("AdjustQueryWhichFlags")
        .value("ADJUST_IGNORENONE", MolOps::ADJUST_IGNORENONE)
        .value("ADJUST_IGNORECHAINS", MolOps::ADJUST_IGNORECHAINS)
        .value("ADJUST_IGNORERINGS", MolOps::ADJUST_IGNORERINGS)
        .value("ADJUST_IGNOREDUMMIES", MolOps::ADJUST_IGNOREDUMMIES)
        .value("ADJUST_IGNORENONDUMMIES", MolOps::ADJUST_IGNORENONDUMMIES)
        .value("ADJUST_IGNOREALL", MolOps::ADJUST_IGNOREALL)
        .export_values();

    docString =
        "Parameters controlling which components of the query atoms are adjusted.\n\
\n\
Attributes:\n\
  - adjustDegree: \n\
      modified atoms have an explicit-degree query added based on their degree in the query \n\
  - adjustHeavyDegree: \n\
      modified atoms have a heavy-atom-degree query added based on their degree in the query \n\
  - adjustDegreeFlags: \n\
      controls which atoms have a degree query added \n\
  - adjustRingCount: \n\
      modified atoms have a ring-count query added based on their ring count in the query \n\
  - adjustRingCountFlags: \n\
      controls which atoms have a ring-count query added \n\
  - makeDummiesQueries: \n\
      dummy atoms that do not have a specified isotope are converted to any-atom queries \n\
  - aromatizeIfPossible: \n\
      attempts aromaticity perception on the molecule \n\
  - makeBondsGeneric: \n\
      convert bonds to generic (any) bonds \n\
  - makeBondsGenericFlags: \n\
      controls which bonds are made generic \n\
  - makeAtomsGeneric: \n\
      convert atoms to generic (any) atoms \n\
  - makeAtomsGenericFlags: \n\
      controls which atoms are made generic \n\
  - adjustRingChain: \n\
      modified atoms have a ring-chain query added based on whether or not they are in a ring \n\
  - adjustRingChainFlags: \n\
      controls which atoms have a ring-chain query added \n\
\n\
A note on the flags controlling which atoms/bonds are modified: \n\
   These generally limit the set of atoms/bonds to be modified.\n\
   For example if ADJUST_RINGSONLY is set, then only atoms in rings will be modified.\n\
       ADJUST_IGNORENONE causes all atoms to be modified\n\
       ADJUST_SETALL sets all of the ADJUST flags\n\
   Some of the options obviously make no sense for bonds\n\
";
    python::class_<MolOps::AdjustQueryParameters>("AdjustQueryParameters",
                                                  docString.c_str())
        .def_readwrite("adjustDegree",
                       &MolOps::AdjustQueryParameters::adjustDegree)
        .def_readwrite("adjustDegreeFlags",
                       &MolOps::AdjustQueryParameters::adjustDegreeFlags)
        .def_readwrite("adjustHeavyDegree",
                       &MolOps::AdjustQueryParameters::adjustHeavyDegree)
        .def_readwrite("adjustHeavyDegreeFlags",
                       &MolOps::AdjustQueryParameters::adjustHeavyDegreeFlags)
        .def_readwrite("adjustRingCount",
                       &MolOps::AdjustQueryParameters::adjustRingCount)
        .def_readwrite("adjustRingCountFlags",
                       &MolOps::AdjustQueryParameters::adjustRingCountFlags)
        .def_readwrite("makeDummiesQueries",
                       &MolOps::AdjustQueryParameters::makeDummiesQueries)
        .def_readwrite("aromatizeIfPossible",
                       &MolOps::AdjustQueryParameters::aromatizeIfPossible)
        .def_readwrite("makeBondsGeneric",
                       &MolOps::AdjustQueryParameters::makeBondsGeneric)
        .def_readwrite("makeBondsGenericFlags",
                       &MolOps::AdjustQueryParameters::makeBondsGenericFlags)
        .def_readwrite("makeAtomsGeneric",
                       &MolOps::AdjustQueryParameters::makeAtomsGeneric)
        .def_readwrite("makeAtomsGenericFlags",
                       &MolOps::AdjustQueryParameters::makeAtomsGenericFlags)
        .def_readwrite("adjustRingChain",
                       &MolOps::AdjustQueryParameters::adjustRingChain)
        .def_readwrite("adjustRingChainFlags",
                       &MolOps::AdjustQueryParameters::adjustRingChainFlags);

    docString =
        "Returns a new molecule where the query properties of atoms have been "
        "modified.";
    python::def("AdjustQueryProperties", adjustQueryPropertiesHelper,
                (python::arg("mol"), python::arg("params") = python::object()),
                docString.c_str(),
                python::return_value_policy<python::manage_new_object>());
  };
};
}  // namespace RDKit

void wrap_molops() { RDKit::molops_wrapper::wrap(); }
