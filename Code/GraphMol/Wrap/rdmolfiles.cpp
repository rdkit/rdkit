//
//  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "rdmolops.h"

#include <RDGeneral/BoostStartInclude.h>
#include <RDBoost/python.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostStartInclude.h>

#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/Canon.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/FileParsers/SequenceWriters.h>
#include <GraphMol/FileParsers/PNGParser.h>
#include <GraphMol/MarvinParse/MarvinParser.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>

#include <RDBoost/Wrap.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SanitException.h>
#include <string.h>

namespace python = boost::python;
using namespace RDKit;

void rdBadFileExceptionTranslator(RDKit::BadFileException const &x) {
  std::ostringstream ss;
  ss << "File error: " << x.what();
  PyErr_SetString(PyExc_IOError, ss.str().c_str());
}
void rdFileParseExceptionTranslator(RDKit::FileParseException const &x) {
  std::ostringstream ss;
  ss << "File parsing error: " << x.what();
  PyErr_SetString(PyExc_RuntimeError, ss.str().c_str());
}

namespace RDKit {
std::string pyObjectToString(python::object input) {
  python::extract<std::string> ex(input);
  if (ex.check()) {
    return ex();
  }
  std::wstring ws = python::extract<std::wstring>(input);
  return std::string(ws.begin(), ws.end());
}

ROMol *MolFromSmiles(python::object ismiles, bool sanitize,
                     python::dict replDict) {
  std::map<std::string, std::string> replacements;
  const auto items = replDict.items();
  for (unsigned int i = 0; i < python::len(items); ++i) {
    const auto item = items[i];
    replacements[python::extract<std::string>(item[0])] =
        python::extract<std::string>(item[1]);
  }
  RWMol *newM;
  std::string smiles = pyObjectToString(ismiles);
  try {
    newM = SmilesToMol(smiles, 0, sanitize, &replacements);
  } catch (...) {
    newM = nullptr;
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromSmarts(python::object ismarts, bool mergeHs,
                     python::dict replDict) {
  std::map<std::string, std::string> replacements;
  const auto items = replDict.items();
  for (unsigned int i = 0; i < python::len(items); ++i) {
    const auto item = items[i];
    replacements[python::extract<std::string>(item[0])] =
        python::extract<std::string>(item[1]);
  }
  std::string smarts = pyObjectToString(ismarts);

  RWMol *newM;
  try {
    newM = SmartsToMol(smarts, 0, mergeHs, &replacements);
  } catch (...) {
    newM = nullptr;
  }
  return static_cast<ROMol *>(newM);
}
ROMol *MolFromTPLFile(const std::string &filename, bool sanitize = true,
                      bool skipFirstConf = false) {
  RWMol *newM;
  try {
    newM = TPLFileToMol(filename, sanitize, skipFirstConf);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (...) {
    newM = nullptr;
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromTPLBlock(python::object itplBlock, bool sanitize = true,
                       bool skipFirstConf = false) {
  std::istringstream inStream(pyObjectToString(itplBlock));
  unsigned int line = 0;
  RWMol *newM;
  try {
    newM = TPLDataStreamToMol(&inStream, line, sanitize, skipFirstConf);
  } catch (...) {
    newM = nullptr;
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromMolFileHelper(const std::string &molFilename, bool sanitize,
                            bool removeHs, bool strictParsing) {
  RWMol *newM = nullptr;
  try {
    newM = MolFileToMol(molFilename, sanitize, removeHs, strictParsing);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromMolBlock(python::object imolBlock, bool sanitize, bool removeHs,
                       bool strictParsing) {
  std::istringstream inStream(pyObjectToString(imolBlock));
  unsigned int line = 0;
  RWMol *newM = nullptr;
  try {
    newM =
        MolDataStreamToMol(inStream, line, sanitize, removeHs, strictParsing);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromMolFile(const std::string &molFilename, bool sanitize,
                      bool removeHs, bool strictParsing) {
  RWMol *newM = nullptr;
  try {
    newM = MolFileToMol(molFilename, sanitize, removeHs, strictParsing);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

RDKit::ROMol *MolFromSCSRBlock(const std::string &molBlock, bool sanitize,
                               bool removeHs, python::object pyParams) {
  RDKit::v2::FileParsers::MolFromSCSRParams scsrParams;
  if (pyParams) {
    scsrParams =
        python::extract<RDKit::v2::FileParsers::MolFromSCSRParams>(pyParams);
  }
  std::istringstream inStream(molBlock);
  unsigned int line = 0;
  try {
    RDKit::v2::FileParsers::MolFileParserParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    params.strictParsing = false;
    auto mol = RDKit::v2::FileParsers::MolFromSCSRDataStream(
        inStream, line, params, scsrParams);

    return static_cast<ROMol *>(mol.release());

  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(nullptr);
}

RDKit::ROMol *MolFromSCSRFile(const std::string &molFilename, bool sanitize,
                              bool removeHs, python::object pyParams) {
  RDKit::v2::FileParsers::MolFromSCSRParams scsrParams;
  if (pyParams) {
    scsrParams =
        python::extract<RDKit::v2::FileParsers::MolFromSCSRParams>(pyParams);
  }
  try {
    RDKit::v2::FileParsers::MolFileParserParams params;
    params.sanitize = sanitize;
    params.removeHs = removeHs;
    params.strictParsing = false;
    auto mol = RDKit::v2::FileParsers::MolFromSCSRFile(molFilename, params,
                                                       scsrParams);

    return static_cast<ROMol *>(mol.release());

  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(nullptr);
}

ROMol *MolFromMrvFile(const std::string &molFilename, bool sanitize,
                      bool removeHs) {
  RWMol *newM = nullptr;
  try {
    newM = MrvFileToMol(molFilename, sanitize, removeHs);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromMrvBlock(python::object imolBlock, bool sanitize, bool removeHs) {
  std::istringstream inStream(pyObjectToString(imolBlock));
  RWMol *newM = nullptr;
  try {
    newM = MrvDataStreamToMol(inStream, sanitize, removeHs);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromXYZFile(const char *xyzFilename) {
  RWMol *newM = nullptr;
  try {
    newM = XYZFileToMol(xyzFilename);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromXYZBlock(python::object ixyzBlock) {
  std::istringstream inStream(pyObjectToString(ixyzBlock));
  RWMol *newM = nullptr;
  try {
    newM = XYZDataStreamToMol(inStream);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromSVG(python::object imolBlock, bool sanitize, bool removeHs) {
  RWMol *res = nullptr;
  res = RDKitSVGToMol(pyObjectToString(imolBlock), sanitize, removeHs);
  return static_cast<ROMol *>(res);
}

ROMol *MolFromMol2File(const std::string &molFilename, bool sanitize = true,
                       bool removeHs = true, bool cleanupSubstructures = true) {
  RWMol *newM;
  try {
    newM = Mol2FileToMol(molFilename, sanitize, removeHs, Mol2Type::CORINA,
                         cleanupSubstructures);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (...) {
    newM = nullptr;
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromMol2Block(std::string mol2Block, bool sanitize = true,
                        bool removeHs = true,
                        bool cleanupSubstructures = true) {
  std::istringstream inStream(mol2Block);
  RWMol *newM;
  try {
    newM = Mol2DataStreamToMol(inStream, sanitize, removeHs, Mol2Type::CORINA,
                               cleanupSubstructures);
  } catch (...) {
    newM = nullptr;
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromPDBFile(const std::string &filename, bool sanitize, bool removeHs,
                      unsigned int flavor, bool proximityBonding) {
  RWMol *newM = nullptr;
  try {
    newM = PDBFileToMol(filename, sanitize, removeHs, flavor, proximityBonding);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromPDBBlock(python::object molBlock, bool sanitize, bool removeHs,
                       unsigned int flavor, bool proximityBonding) {
  std::istringstream inStream(pyObjectToString(molBlock));
  RWMol *newM = nullptr;
  try {
    newM = PDBDataStreamToMol(inStream, sanitize, removeHs, flavor,
                              proximityBonding);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromSequence(python::object seq, bool sanitize, int flavor) {
  RWMol *newM = nullptr;
  try {
    newM = SequenceToMol(pyObjectToString(seq), sanitize, flavor);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}
ROMol *MolFromFASTA(python::object seq, bool sanitize, int flavor) {
  RWMol *newM = nullptr;
  try {
    newM = FASTAToMol(pyObjectToString(seq), sanitize, flavor);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}
ROMol *MolFromHELM(python::object seq, bool sanitize) {
  RWMol *newM = nullptr;
  try {
    newM = HELMToMol(pyObjectToString(seq), sanitize);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

std::string molFragmentToSmarts(const ROMol &mol, python::object atomsToUse,
                                python::object bondsToUse,
                                bool doIsomericSmarts = true) {
  auto atomIndices =
      pythonObjectToVect(atomsToUse, static_cast<int>(mol.getNumAtoms()));
  if (!atomIndices) {
    throw_value_error("atomsToUse argument must be non-empty");
  }
  auto bondIndices =
      pythonObjectToVect(bondsToUse, static_cast<int>(mol.getNumBonds()));
  return RDKit::MolFragmentToSmarts(mol, *atomIndices, bondIndices.get(),
                                    doIsomericSmarts);
}

std::string molFragmentToCXSmarts(const ROMol &mol, python::object atomsToUse,
                                  python::object bondsToUse,
                                  bool doIsomericSmarts = true) {
  auto atomIndices =
      pythonObjectToVect(atomsToUse, static_cast<int>(mol.getNumAtoms()));
  if (!atomIndices) {
    throw_value_error("atomsToUse argument must be non-empty");
  }
  auto bondIndices =
      pythonObjectToVect(bondsToUse, static_cast<int>(mol.getNumBonds()));
  return RDKit::MolFragmentToCXSmarts(mol, *atomIndices, bondIndices.get(),
                                      doIsomericSmarts);
}

struct smilesfrag_gen {
  std::string operator()(const ROMol &mol, const SmilesWriteParams &ps,
                         const std::vector<int> &atomsToUse,
                         const std::vector<int> *bondsToUse,
                         const std::vector<std::string> *atomSymbols,
                         const std::vector<std::string> *bondSymbols) {
    return MolFragmentToSmiles(mol, ps, atomsToUse, bondsToUse, atomSymbols,
                               bondSymbols);
  }
};
struct cxsmilesfrag_gen {
  std::string operator()(const ROMol &mol, const SmilesWriteParams &ps,
                         const std::vector<int> &atomsToUse,
                         const std::vector<int> *bondsToUse,
                         const std::vector<std::string> *atomSymbols,
                         const std::vector<std::string> *bondSymbols) {
    return MolFragmentToCXSmiles(mol, ps, atomsToUse, bondsToUse, atomSymbols,
                                 bondSymbols);
  }
};

template <typename F>
std::string MolFragmentToSmilesHelper1(const ROMol &mol,
                                       const SmilesWriteParams &params,
                                       python::object atomsToUse,
                                       python::object bondsToUse,
                                       python::object atomSymbols,
                                       python::object bondSymbols) {
  auto avect =
      pythonObjectToVect(atomsToUse, static_cast<int>(mol.getNumAtoms()));
  if (!avect.get() || !(avect->size())) {
    throw_value_error("atomsToUse must not be empty");
  }
  auto bvect =
      pythonObjectToVect(bondsToUse, static_cast<int>(mol.getNumBonds()));
  std::unique_ptr<std::vector<std::string>> asymbols =
      pythonObjectToVect<std::string>(atomSymbols);
  std::unique_ptr<std::vector<std::string>> bsymbols =
      pythonObjectToVect<std::string>(bondSymbols);
  if (asymbols.get() && asymbols->size() != mol.getNumAtoms()) {
    throw_value_error("length of atom symbol list != number of atoms");
  }
  if (bsymbols.get() && bsymbols->size() != mol.getNumBonds()) {
    throw_value_error("length of bond symbol list != number of bonds");
  }

  std::string res = F()(mol, params, *avect.get(), bvect.get(), asymbols.get(),
                        bsymbols.get());
  return res;
}

template <typename F>
std::string MolFragmentToSmilesHelper2(
    const ROMol &mol, python::object atomsToUse, python::object bondsToUse,
    python::object atomSymbols, python::object bondSymbols,
    bool doIsomericSmiles, bool doKekule, int rootedAtAtom, bool canonical,
    bool allBondsExplicit, bool allHsExplicit) {
  SmilesWriteParams ps;
  ps.doIsomericSmiles = doIsomericSmiles;
  ps.doKekule = doKekule;
  ps.rootedAtAtom = rootedAtAtom;
  ps.canonical = canonical;
  ps.allBondsExplicit = allBondsExplicit;
  ps.allHsExplicit = allHsExplicit;
  return MolFragmentToSmilesHelper1<F>(mol, ps, atomsToUse, bondsToUse,
                                       atomSymbols, bondSymbols);
}
std::vector<unsigned int> CanonicalRankAtoms(
    const ROMol &mol, bool breakTies = true, bool includeChirality = true,
    bool includeIsotopes = true, bool includeAtomMaps = true,
    bool includeChiralPresence = false) {
  std::vector<unsigned int> ranks(mol.getNumAtoms());
  const bool includeStereoGroups = true;
  const bool useNonStereoRanks = false;

  Canon::rankMolAtoms(mol, ranks, breakTies, includeChirality, includeIsotopes,
                      includeAtomMaps, includeChiralPresence,
                      includeStereoGroups, useNonStereoRanks);
  return ranks;
}

std::vector<int> CanonicalRankAtomsInFragment(
    const ROMol &mol, python::object atomsToUse, python::object bondsToUse,
    python::object atomSymbols, bool breakTies = true,
    bool includeChirality = true, bool includeIsotopes = true,
    bool includeAtomMaps = true, bool includeChiralPresence = false) {
  std::unique_ptr<std::vector<int>> avect =
      pythonObjectToVect(atomsToUse, static_cast<int>(mol.getNumAtoms()));
  if (!avect.get() || !(avect->size())) {
    throw_value_error("atomsToUse must not be empty");
  }
  std::unique_ptr<std::vector<int>> bvect =
      pythonObjectToVect(bondsToUse, static_cast<int>(mol.getNumBonds()));
  std::unique_ptr<std::vector<std::string>> asymbols =
      pythonObjectToVect<std::string>(atomSymbols);
  if (asymbols.get() && asymbols->size() != mol.getNumAtoms()) {
    throw_value_error("length of atom symbol list != number of atoms");
  }

  boost::dynamic_bitset<> atoms(mol.getNumAtoms());
  for (size_t i = 0; i < avect->size(); ++i) {
    atoms[(*avect)[i]] = true;
  }

  boost::dynamic_bitset<> bonds(mol.getNumBonds());
  for (size_t i = 0; bvect.get() && i < bvect->size(); ++i) {
    bonds[(*bvect)[i]] = true;
  }

  std::vector<unsigned int> ranks(mol.getNumAtoms());
  Canon::rankFragmentAtoms(mol, ranks, atoms, bonds, asymbols.get(), breakTies,
                           includeChirality, includeIsotopes,
                           includeChiralPresence);

  std::vector<int> resRanks(mol.getNumAtoms());
  // set unused ranks to -1 for the Python interface
  for (size_t i = 0; i < atoms.size(); ++i) {
    if (!atoms[i]) {
      resRanks[i] = -1;
    } else {
      resRanks[i] = static_cast<int>(ranks[i]);
    }
  }

  return resRanks;
}

ROMol *MolFromSmilesHelper(python::object ismiles,
                           const v2::SmilesParse::SmilesParserParams &params) {
  std::string smiles = pyObjectToString(ismiles);

  try {
    return v2::SmilesParse::MolFromSmiles(smiles, params).release();
  } catch (...) {
    return nullptr;
  }
}

ROMol *MolFromSmartsHelper(python::object ismiles,
                           const v2::SmilesParse::SmartsParserParams &params) {
  std::string smiles = pyObjectToString(ismiles);

  try {
    return v2::SmilesParse::MolFromSmarts(smiles, params).release();
  } catch (...) {
    return nullptr;
  }
}

python::list MolToRandomSmilesHelper(const ROMol &mol, unsigned int numSmiles,
                                     unsigned int randomSeed,
                                     bool doIsomericSmiles, bool doKekule,
                                     bool allBondsExplicit,
                                     bool allHsExplicit) {
  auto res = MolToRandomSmilesVect(mol, numSmiles, randomSeed, doIsomericSmiles,
                                   doKekule, allBondsExplicit, allHsExplicit);
  python::list pyres;
  for (auto smi : res) {
    pyres.append(smi);
  }
  return pyres;
}

ROMol *MolFromPNGFile(const std::string &filename, python::object pyParams) {
  SmilesParserParams params;
  if (pyParams) {
    params = python::extract<SmilesParserParams>(pyParams);
  }
  ROMol *newM = nullptr;
  try {
    newM = PNGFileToMol(filename, params);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return newM;
}

ROMol *MolFromPNGString(python::object png, python::object pyParams) {
  SmilesParserParams params;
  if (pyParams) {
    params = python::extract<SmilesParserParams>(pyParams);
  }
  ROMol *newM = nullptr;
  try {
    newM = PNGStringToMol(pyObjectToString(png), params);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return newM;
}

python::object addMolToPNGFileHelperParams(const ROMol &mol,
                                           python::object fname,
                                           const PNGMetadataParams &params) {
  std::string cstr = python::extract<std::string>(fname);

  auto res = addMolToPNGFile(mol, cstr, params);

  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

python::object addMolToPNGFileHelper(const ROMol &mol, python::object fname,
                                     bool includePkl, bool includeSmiles,
                                     bool includeMol) {
  PNGMetadataParams params;
  params.includePkl = includePkl;
  params.includeSmiles = includeSmiles;
  params.includeMol = includeMol;
  return addMolToPNGFileHelperParams(mol, fname, params);
}

python::object addMolToPNGStringHelperParams(const ROMol &mol,
                                             python::object png,
                                             const PNGMetadataParams &params) {
  std::string cstr = python::extract<std::string>(png);

  auto res = addMolToPNGString(mol, cstr, params);

  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

python::object addMolToPNGStringHelper(const ROMol &mol, python::object png,
                                       bool includePkl, bool includeSmiles,
                                       bool includeMol) {
  PNGMetadataParams params;
  params.includePkl = includePkl;
  params.includeSmiles = includeSmiles;
  params.includeMol = includeMol;
  return addMolToPNGStringHelperParams(mol, png, params);
}

python::object addMetadataToPNGFileHelper(python::dict pymetadata,
                                          python::object fname) {
  std::string cstr = python::extract<std::string>(fname);

  std::vector<std::pair<std::string, std::string>> metadata;

  const auto items = pymetadata.items();
  for (unsigned int i = 0; i < python::len(items); ++i) {
    const auto item = items[i];
    std::string key = python::extract<std::string>(item[0]);
    std::string val = python::extract<std::string>(item[1]);
    metadata.push_back(std::make_pair(key, val));
  }

  auto res = addMetadataToPNGFile(cstr, metadata);

  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

python::object addMetadataToPNGStringHelper(python::dict pymetadata,
                                            python::object png) {
  std::string cstr = python::extract<std::string>(png);

  std::vector<std::pair<std::string, std::string>> metadata;
  const auto items = pymetadata.items();
  for (unsigned int i = 0; i < python::len(items); ++i) {
    const auto item = items[i];
    std::string key = python::extract<std::string>(item[0]);
    std::string val = python::extract<std::string>(item[1]);
    metadata.push_back(std::make_pair(key, val));
  }

  auto res = addMetadataToPNGString(cstr, metadata);

  python::object retval = python::object(
      python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

python::object MolsFromPNGFile(const std::string &filename,
                               const std::string &tag,
                               python::object pyParams) {
  SmilesParserParams params;
  if (pyParams) {
    params = python::extract<SmilesParserParams>(pyParams);
  }
  std::vector<std::unique_ptr<ROMol>> mols;
  try {
    mols = PNGFileToMols(filename, tag, params);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  python::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(mol.release());
    res.append(sptr);
  }
  return python::tuple(res);
}

python::tuple MolsFromPNGString(python::object png, const std::string &tag,
                                python::object pyParams) {
  SmilesParserParams params;
  if (pyParams) {
    params = python::extract<SmilesParserParams>(pyParams);
  }
  auto mols = PNGStringToMols(pyObjectToString(png), tag, params);
  python::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(mol.release());
    res.append(sptr);
  }
  return python::tuple(res);
}

python::object MolsFromCDXMLFile(const std::string &filename, bool sanitize,
                                 bool removeHs) {
  std::vector<std::unique_ptr<RWMol>> mols;
  try {
    mols = CDXMLFileToMols(filename, sanitize, removeHs);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  python::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return python::tuple(res);
}

python::tuple MolsFromCDXMLHelper(python::object cdxml,
                                  python::object pyParams) {
  RDKit::v2::CDXMLParser::CDXMLParserParams params;
  if (pyParams) {
    params =
        python::extract<RDKit::v2::CDXMLParser::CDXMLParserParams>(pyParams);
  }
  auto mols =
      RDKit::v2::CDXMLParser::MolsFromCDXML(pyObjectToString(cdxml), params);
  python::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return python::tuple(res);
}

python::object MolsFromCDXMLFileHelper(const std::string &filename,
                                       python::object pyParams) {
  RDKit::v2::CDXMLParser::CDXMLParserParams params(
      true, true, RDKit::v2::CDXMLParser::CDXMLFormat::Auto);
  if (pyParams) {
    params =
        python::extract<RDKit::v2::CDXMLParser::CDXMLParserParams>(pyParams);
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  try {
    mols = RDKit::v2::CDXMLParser::MolsFromCDXMLFile(filename, params);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  python::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return python::tuple(res);
}

python::tuple MolsFromCDXML(python::object cdxml, bool sanitize,
                            bool removeHs) {
  auto mols = CDXMLToMols(pyObjectToString(cdxml), sanitize, removeHs);
  python::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return python::tuple(res);
}
namespace {
python::object translateMetadataToList(
    const std::vector<std::pair<std::string, std::string>> &metadata) {
  python::list resAsList;
  for (const auto &[key, value] : metadata) {
    // keys are safe to extract:
    // but values may include binary, so we convert them directly to bytes:
    python::object val = python::object(python::handle<>(
        PyBytes_FromStringAndSize(value.c_str(), value.length())));
    resAsList.append(python::make_tuple(key, val));
  }
  return resAsList;
}
python::object translateMetadataToDict(
    const std::vector<std::pair<std::string, std::string>> &metadata) {
  python::dict resAsDict;
  for (const auto &[key, value] : metadata) {
    // keys are safe to extract:
    // but values may include binary, so we convert them directly to bytes:
    resAsDict[key] = python::object(python::handle<>(
        PyBytes_FromStringAndSize(value.c_str(), value.length())));
  }
  return resAsDict;
}

}  // namespace
python::object MetadataFromPNGFile(python::object fname, bool asList) {
  std::string cstr = python::extract<std::string>(fname);
  auto metadata = PNGFileToMetadata(cstr);
  if (asList) {
    return translateMetadataToList(metadata);
  }
  return translateMetadataToDict(metadata);
}

python::object MetadataFromPNGString(python::object png, bool asList) {
  std::string cstr = python::extract<std::string>(png);
  auto metadata = PNGStringToMetadata(cstr);
  if (asList) {
    return translateMetadataToList(metadata);
  }
  return translateMetadataToDict(metadata);
}

void CanonicalizeEnhancedStereo(ROMol &mol) {
  Canon::canonicalizeEnhancedStereo(mol);
}

std::string MolToV2KMolBlockHelper(const ROMol &mol, python::object pyParams,
                                   int confId) {
  MolWriterParams params;
  if (pyParams) {
    params = python::extract<MolWriterParams>(pyParams);
  }
  return MolToV2KMolBlock(mol, params, confId);
}

}  // namespace RDKit

// MolSupplier stuff
#ifdef SUPPORT_COMPRESSED_SUPPLIERS
void wrap_compressedsdsupplier();
#endif
void wrap_sdsupplier();
void wrap_forwardsdsupplier();
void wrap_tdtsupplier();
void wrap_smisupplier();
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
void wrap_maesupplier();
#endif

// mol writer stuff
void wrap_smiwriter();
void wrap_sdwriter();
void wrap_tdtwriter();
void wrap_pdbwriter();
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
void wrap_maewriter();
#endif

// MultithreadedMolSupplier stuff
void wrap_multiSmiSupplier();
void wrap_multiSDSupplier();

BOOST_PYTHON_MODULE(rdmolfiles) {
  std::string docString;

  python::scope().attr("__doc__") =
      "Module containing RDKit functionality for working with molecular file "
      "formats.";
  python::register_exception_translator<RDKit::BadFileException>(
      &rdBadFileExceptionTranslator);

  python::register_exception_translator<RDKit::FileParseException>(
      &rdFileParseExceptionTranslator);

  docString =
      "Construct a molecule from a TPL file.\n\n\
  ARGUMENTS:\n\
\n\
    - fileName: name of the file to read\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - skipFirstConf: (optional) skips reading the first conformer.\n\
      Defaults to False.\n\
      This should be set to True when reading TPLs written by \n\
      the CombiCode.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromTPLFile", RDKit::MolFromTPLFile,
              (python::arg("fileName"), python::arg("sanitize") = true,
               python::arg("skipFirstConf") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a TPL block.\n\n\
  ARGUMENTS:\n\
\n\
    - tplBlock: string containing the TPL block\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - skipFirstConf: (optional) skips reading the first conformer.\n\
      Defaults to False.\n\
      This should be set to True when reading TPLs written by \n\
      the CombiCode.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromTPLBlock", RDKit::MolFromTPLBlock,
              (python::arg("tplBlock"), python::arg("sanitize") = true,
               python::arg("skipFirstConf") = false),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a Mol file.\n\n\
  ARGUMENTS:\n\
\n\
    - fileName: name of the file to read\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to true.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
    - strictParsing: (optional) if this is false, the parser is more lax about.\n\
      correctness of the content.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def(
      "MolFromMolFile", RDKit::MolFromMolFileHelper,
      (python::arg("molFileName"), python::arg("sanitize") = true,
       python::arg("removeHs") = true, python::arg("strictParsing") = true),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a Mol block.\n\n\
  ARGUMENTS:\n\
\n\
    - molBlock: string containing the Mol block\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
    - strictParsing: (optional) if this is false, the parser is more lax about.\n\
      correctness of the content.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def(
      "MolFromMolBlock", RDKit::MolFromMolBlock,
      (python::arg("molBlock"), python::arg("sanitize") = true,
       python::arg("removeHs") = true, python::arg("strictParsing") = true),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a Marvin (Mrv) file.\n\n\
  ARGUMENTS:\n\
\n\
    - fileName: name of the file to read\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to true.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromMrvFile", RDKit::MolFromMrvFile,
              (python::arg("molFileName"), python::arg("sanitize") = true,
               python::arg("removeHs") = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a Marvin (mrv) block.\n\n\
  ARGUMENTS:\n\
\n\
    - molBlock: string containing the Marvin block\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromMrvBlock", RDKit::MolFromMrvBlock,
              (python::arg("mrvBlock"), python::arg("sanitize") = true,
               python::arg("removeHs") = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from an XYZ file.\n\n\
  ARGUMENTS:\n\
\n\
    - xyzFileName: name of the file to read\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromXYZFile", RDKit::MolFromXYZFile,
              (python::arg("xyzFileName")), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from an XYZ string.\n\n\
  ARGUMENTS:\n\
\n\
    - xyzBlock: the XYZ data to read\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromXYZBlock", RDKit::MolFromXYZBlock,
              (python::arg("xyzBlock")), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from an RDKit-generate SVG string.\n\n\
  ARGUMENTS:\n\
\n\
    - svg: string containing the SVG data (must include molecule metadata)\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n\
  NOTE: this functionality should be considered beta.\n\
\n";
  python::def("MolFromRDKitSVG", RDKit::MolFromSVG,
              (python::arg("svg"), python::arg("sanitize") = true,
               python::arg("removeHs") = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a Tripos Mol2 file.\n\n\
  NOTE:\n\
    The parser expects the atom-typing scheme used by Corina.\n\
    Atom types from Tripos' dbtranslate are less supported.\n\
    Other atom typing schemes are unlikely to work.\n\
\n\
  ARGUMENTS:\n                                  \
\n\
    - mol2FileName: name of the file to read\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to true.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
    - cleanupSubstructures: (optional) toggles standardizing some \n\
      substructures found in mol2 files.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromMol2File", RDKit::MolFromMol2File,
              (python::arg("mol2FileName"), python::arg("sanitize") = true,
               python::arg("removeHs") = true,
               python::arg("cleanupSubstructures") = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a Tripos Mol2 block.\n\n\
  NOTE:\n\
    The parser expects the atom-typing scheme used by Corina.\n\
    Atom types from Tripos' dbtranslate are less supported.\n\
    Other atom typing schemes are unlikely to work.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol2Block: string containing the Mol2 block\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
    - cleanupSubstructures: (optional) toggles standardizing some \n\
      substructures found in mol2 files.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromMol2Block", RDKit::MolFromMol2Block,
              (python::arg("mol2Block"), python::arg("sanitize") = true,
               python::arg("removeHs") = true,
               python::arg("cleanupSubstructures") = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a Mol file.\n\n\
  ARGUMENTS:\n\
\n\
    - molFileName: name of the file to read\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to true.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
    - strictParsing: (optional) if this is false, the parser is more lax about.\n\
      correctness of the content.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def(
      "MolFromMolFile", RDKit::MolFromMolFile,
      (python::arg("molFileName"), python::arg("sanitize") = true,
       python::arg("removeHs") = true, python::arg("strictParsing") = true),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a Mol block.\n\n\
  ARGUMENTS:\n\
\n\
    - molBlock: string containing the Mol block\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
    - strictParsing: (optional) if this is false, the parser is more lax about.\n\
      correctness of the content.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def(
      "MolFromMolBlock", RDKit::MolFromMolBlock,
      (python::arg("molBlock"), python::arg("sanitize") = true,
       python::arg("removeHs") = true, python::arg("strictParsing") = true),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  python::class_<RDKit::MolWriterParams, boost::noncopyable>(
      "MolWriterParams", "Parameters controlling Mol writing")
      .def_readwrite(
          "includeStereo", &RDKit::MolWriterParams::includeStereo,
          "toggles inclusion of stereochemistry information (default=True)")
      .def_readwrite(
          "kekulize", &RDKit::MolWriterParams::kekulize,
          "triggers kekulization of the molecule before it is written (default=True)")
      .def_readwrite(
          "forceV3000", &RDKit::MolWriterParams::forceV3000,
          "force generation a V3000 mol block (happens automatically with more than 999 atoms or bonds)(default=False)")
      .def_readwrite(
          "precision", &RDKit::MolWriterParams::precision,
          "precision of coordinates (only available in V3000)(default=false)")
      .def("__setattr__", &safeSetattr);

  python::class_<RDKit::v2::FileParsers::MolFromSCSRParams, boost::noncopyable>(
      "MolFromSCSRParams",
      "Parameters controlling conversion of an SCSRMol to a Mol")
      .def_readwrite(
          "includeLeavingGroups",
          &RDKit::v2::FileParsers::MolFromSCSRParams::includeLeavingGroups,
          "include leaving groups atoms if not substited at that position")
      .def_readwrite(
          "scsrTemplateNames",
          &RDKit::v2::FileParsers::MolFromSCSRParams::scsrTemplateNames,
          "If True, the first template name in the Sgroup is used as the Sgroup label")
      .def_readwrite(
          "scsrBaseHbondOptions",
          &RDKit::v2::FileParsers::MolFromSCSRParams::scsrBaseHbondOptions,
          "One of Ignore, UseSapAll(default) , UseSapOne, Auto")
      .def("__setattr__", &safeSetattr);

  docString =
      "Construct a molecule from an SCSR Mol block.\n\n\
      ARGUMENTS:\n\
    \n\
        - molBlock: string containing the SCSR Mol block\n\
    \n\
        - sanitize: (optional) toggles sanitization of the molecule.\n\
          Defaults to True.\n\
    \n\
        - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
          This only make sense when sanitization is done.\n\
          Defaults to true.\n\
    \n\
        - molFromSCSRParams : MolFromSCSRParams to control conversion\n\
    \n RETURNS :\n\
    \n a Mol object, None on failure.\n\
    \n ";
  python::def("MolFromSCSRBlock", RDKit::MolFromSCSRBlock,
              (python::arg("molBlock"), python::arg("sanitize") = true,
               python::arg("removeHs") = true,
               python::arg("molFromSCSRParams") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from an SCSR Mol block.\n\n\
      ARGUMENTS:\n\
    \n\
        - filename: string containing the SCSR filename\n\
    \n\
        - sanitize: (optional) toggles sanitization of the molecule.\n\
          Defaults to True.\n\
    \n\
        - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
          This only make sense when sanitization is done.\n\
          Defaults to true.\n\
    \n\
        - molFromSCSRParams : MolFromSCSRParams to control conversion\n\
    \n RETURNS :\n\
    \n a Mol object, None on failure.\n\
    \n ";
  python::def("MolFromSCSRFile", RDKit::MolFromSCSRFile,
              (python::arg("filename"), python::arg("sanitize") = true,
               python::arg("removeHs") = true,
               python::arg("molFromSCSRParams") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Returns a Mol block for a molecule\n\
  Arguments:\n\
    - mol: the molecule\n\
    - params: the MolWriterParams\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolToMolBlock",
      (std::string(*)(const ROMol &, const MolWriterParams &,
                      int))RDKit::MolToMolBlock,
      (python::arg("mol"), python::arg("params"), python::arg("confId") = -1),
      docString.c_str());

  docString =
      "Returns a Mol block for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - includeStereo: (optional) toggles inclusion of stereochemical\n\
      information in the output\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
    - kekulize: (optional) triggers kekulization of the molecule before it's written,\n\
      as suggested by the MDL spec.\n\
    - forceV3000 (optional) force generation a V3000 mol block (happens automatically with \n\
      more than 999 atoms or bonds)\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToMolBlock",
              (std::string(*)(const ROMol &, bool, int, bool,
                              bool))RDKit::MolToMolBlock,
              (python::arg("mol"), python::arg("includeStereo") = true,
               python::arg("confId") = -1, python::arg("kekulize") = true,
               python::arg("forceV3000") = false),
              docString.c_str());

  docString =
      "Returns a V3000 Mol block for a molecule\n\
   ARGUMENTS:\n\
\n \
     - mol: the molecule\n\
     - params: the MolWriterParams\n\
     - confId: (optional) selects which conformation to output (-1 = default)\n\
\n \
   RETURNS:\n\
\n \
     a string\n\
\n ";
  python::def(
      "MolToV3KMolBlock",
      (std::string(*)(const ROMol &, const MolWriterParams &,
                      int))RDKit::MolToV3KMolBlock,
      (python::arg("mol"), python::arg("params"), python::arg("confId") = -1),
      docString.c_str());

  docString =
      "Returns a V3000 Mol block for a molecule\n\
   ARGUMENTS:\n\
\n \
     - mol: the molecule\n\
     - includeStereo: (optional) toggles inclusion of stereochemical\n\
       information in the output\n\
     - confId: (optional) selects which conformation to output (-1 = default)\n\
     - kekulize: (optional) triggers kekulization of the molecule before it's written,\n\
       as suggested by the MDL spec.\n\
\n \
   RETURNS:\n\
\n \
     a string\n\
\n ";

  python::def(
      "MolToV3KMolBlock",
      (std::string(*)(const ROMol &, bool, int, bool))RDKit::MolToV3KMolBlock,
      (python::arg("mol"), python::arg("includeStereo") = true,
       python::arg("confId") = -1, python::arg("kekulize") = true),
      docString.c_str());

  docString =
      R"DOC(Returns a V2000 Mol block for a molecule
   ARGUMENTS:

     - mol: the molecule
     - params: the MolWriterParams
     - confId: (optional) selects which conformation to output (-1 = default)

   RETURNS:

     a string

   NOTE: this function throws a ValueError if the molecule has more than 999 atoms, bonds, or SGroups
)DOC";
  python::def("MolToV2KMolBlock", MolToV2KMolBlockHelper,
              (python::arg("mol"), python::arg("params") = python::object(),
               python::arg("confId") = -1),
              docString.c_str());

  docString =
      "Writes a Mol file for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - filename: the file to write to\n\
    - params: the MolWriterParams\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
\n";
  python::def("MolToMolFile",
              (void (*)(const ROMol &, const std::string &,
                        const MolWriterParams &, int))RDKit::MolToMolFile,
              (python::arg("mol"), python::arg("filename"),
               python::arg("params"), python::arg("confId") = -1),
              docString.c_str());

  docString =
      "Writes a Mol file for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - filename: the file to write to\n\
    - includeStereo: (optional) toggles inclusion of stereochemical\n\
      information in the output\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
    - kekulize: (optional) triggers kekulization of the molecule before it's written,\n\
      as suggested by the MDL spec.\n\
    - forceV3000 (optional) force generation a V3000 mol block (happens automatically with \n\
      more than 999 atoms or bonds)\n\
\n";
  python::def(
      "MolToMolFile",
      (void (*)(const ROMol &, const std::string &, bool, int, bool,
                bool))RDKit::MolToMolFile,
      (python::arg("mol"), python::arg("filename"),
       python::arg("includeStereo") = true, python::arg("confId") = -1,
       python::arg("kekulize") = true, python::arg("forceV3000") = false),
      docString.c_str());

  docString =
      "Writes a V3000 Mol file for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - filename: the file to write to\n\
    - params: the MolWriterParams\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
\n";
  python::def("MolToV3KMolFile",
              (void (*)(const ROMol &, const std::string &,
                        const MolWriterParams &, int))RDKit::MolToV3KMolFile,
              (python::arg("mol"), python::arg("filename"),
               python::arg("params") = true, python::arg("confId") = -1),
              docString.c_str());

  docString =
      "Writes a V3000 Mol file for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - filename: the file to write to\n\
    - includeStereo: (optional) toggles inclusion of stereochemical\n\
      information in the output\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
    - kekulize: (optional) triggers kekulization of the molecule before it's written,\n\
      as suggested by the MDL spec.\n\
\n";
  python::def("MolToV3KMolFile",
              (void (*)(const ROMol &, const std::string &, bool, int,
                        bool))RDKit::MolToV3KMolFile,
              (python::arg("mol"), python::arg("filename"),
               python::arg("includeStereo") = true, python::arg("confId") = -1,
               python::arg("kekulize") = true),
              docString.c_str());
  //

  docString =
      "Returns a Marvin (Mrv) Mol block for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - includeStereo: (optional) toggles inclusion of stereochemical\n\
      information in the output\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
    - kekulize: (optional) triggers kekulization of the molecule before it's written.\n\
    - prettyPrint: (optional) makes the output more human readable.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToMrvBlock",
              (std::string(*)(const ROMol &, bool, int, bool,
                              bool))RDKit::MolToMrvBlock,
              (python::arg("mol"), python::arg("includeStereo") = true,
               python::arg("confId") = -1, python::arg("kekulize") = true,
               python::arg("prettyPrint") = false),
              docString.c_str());

  docString =
      "Returns a Marvin (Mrv) Mol block for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - params: marvin write params\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolToMrvBlock",
      (std::string(*)(const ROMol &, const MrvWriterParams &,
                      int))RDKit::MolToMrvBlock,
      (python::arg("mol"), python::arg("params"), python::arg("confId") = -1),
      docString.c_str());

  docString =
      "Writes a Marvin (MRV) file for a molecule\n\
   ARGUMENTS:\n\
\n\
     - mol: the molecule\n\
     - filename: the file to write to\n\
     - includeStereo: (optional) toggles inclusion of stereochemical\n\
       information in the output\n\
     - confId: (optional) selects which conformation to output (-1 = default)\n\
     - kekulize: (optional) triggers kekulization of the molecule before it's written.\n\
     - prettyPrint: (optional) makes the output more human readable.\n\
\n";
  python::def(
      "MolToMrvFile",
      (void (*)(const ROMol &, const std::string &, bool, int, bool,
                bool))RDKit::MolToMrvFile,
      (python::arg("mol"), python::arg("filename"),
       python::arg("includeStereo") = true, python::arg("confId") = -1,
       python::arg("kekulize") = true, python::arg("prettyPrint") = false),
      docString.c_str());

  docString =
      "Writes a Marvin (MRV) file for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - filename: the file to write to\n\
    - params: marvin write params\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
\n";
  python::def("MolToMrvFile",
              (void (*)(const ROMol &, const std::string &,
                        const MrvWriterParams &, int))RDKit::MolToMrvFile,
              (python::arg("mol"), python::arg("filename"),
               python::arg("params"), python::arg("confId") = -1),
              docString.c_str());

  docString =
      "Writes a CML block for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - confId: (optional) selects which conformation to output\n\
    - kekulize: (optional) triggers kekulization of the molecule before it's written\n\
\n";
  python::def("MolToCMLBlock", RDKit::MolToCMLBlock,
              (python::arg{"mol"}, python::arg{"confId"} = -1,
               python::arg{"kekulize"} = true),
              docString.c_str());

  docString =
      "Writes a CML file for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - filename: the file to write to\n\
    - confId: (optional) selects which conformation to output\n\
    - kekulize: (optional) triggers kekulization of the molecule before it's written\n\
\n";
  python::def("MolToCMLFile", RDKit::MolToCMLFile,
              (python::arg{"mol"}, python::arg{"filename"},
               python::arg{"confId"} = -1, python::arg{"kekulize"} = true),
              docString.c_str());

  docString =
      "Returns a XYZ block for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
    - precision: precision of the coordinates\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToXYZBlock", RDKit::MolToXYZBlock,
              (python::arg("mol"), python::arg("confId") = -1,
               python::arg("precision") = 6),
              docString.c_str());

  docString =
      "Writes a XYZ file for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - filename: the file to write to\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
    - precision: precision of the coordinates\n\
\n";
  python::def("MolToXYZFile", RDKit::MolToXYZFile,
              (python::arg("mol"), python::arg("filename"),
               python::arg("confId") = -1, python::arg("precision") = 6),
              docString.c_str());

  //

  python::class_<RDKit::v2::SmilesParse::SmilesParserParams,
                 boost::noncopyable>("SmilesParserParams",
                                     "Parameters controlling SMILES Parsing")
      .def_readwrite("debugParse",
                     &RDKit::v2::SmilesParse::SmilesParserParams::debugParse,
                     "controls the amount of debugging information produced")
      .def_readwrite("parseName",
                     &RDKit::v2::SmilesParse::SmilesParserParams::parseName,
                     "controls whether or not the molecule name is also parsed")
      .def_readwrite(
          "allowCXSMILES",
          &RDKit::v2::SmilesParse::SmilesParserParams::allowCXSMILES,
          "controls whether or not the CXSMILES extensions are parsed")
      .def_readwrite(
          "strictCXSMILES",
          &RDKit::v2::SmilesParse::SmilesParserParams::strictCXSMILES,
          "controls whether or not problems in CXSMILES parsing "
          "causes molecule parsing to fail")
      .def_readwrite("sanitize",
                     &RDKit::v2::SmilesParse::SmilesParserParams::sanitize,
                     "controls whether or not the molecule is sanitized before "
                     "being returned")
      .def_readwrite("removeHs",
                     &RDKit::v2::SmilesParse::SmilesParserParams::removeHs,
                     "controls whether or not Hs are removed before the "
                     "molecule is returned")
      .def_readwrite("flagPossible",
                     &RDKit::v2::SmilesParse::SmilesParserParams::flagPossible,
                     "controls whether or not to flag possible chiral centers")
      .def("__setattr__", &safeSetattr);
  python::class_<RDKit::v2::SmilesParse::SmartsParserParams,
                 boost::noncopyable>("SmartsParserParams",
                                     "Parameters controlling SMARTS Parsing")
      .def_readwrite("debugParse",
                     &RDKit::v2::SmilesParse::SmartsParserParams::debugParse,
                     "controls the amount of debugging information produced")
      .def_readwrite("parseName",
                     &RDKit::v2::SmilesParse::SmartsParserParams::parseName,
                     "controls whether or not the molecule name is also parsed")
      .def_readwrite(
          "allowCXSMILES",
          &RDKit::v2::SmilesParse::SmartsParserParams::allowCXSMILES,
          "controls whether or not the CXSMILES extensions are parsed")
      .def_readwrite(
          "strictCXSMILES",
          &RDKit::v2::SmilesParse::SmartsParserParams::strictCXSMILES,
          "controls whether or not problems in CXSMILES parsing "
          "causes molecule parsing to fail")
      .def_readwrite(
          "mergeHs", &RDKit::v2::SmilesParse::SmartsParserParams::mergeHs,
          "toggles merging H atoms in the SMARTS into neighboring atoms")
      .def("__setattr__", &safeSetattr);

  docString =
      "Construct a molecule from a SMILES string.\n\n\
     ARGUMENTS:\n\
\n\
       - SMILES: the smiles string\n\
\n\
       - params: used to provide optional parameters for the SMILES parsing\n\
\n\
     RETURNS:\n\
\n\
       a Mol object, None on failure.\n\
\n";
  python::def("MolFromSmiles", MolFromSmilesHelper,
              (python::arg("SMILES"), python::arg("params")), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a SMILES string.\n\n\
  ARGUMENTS:\n\
\n\
    - SMILES: the smiles string\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - replacements: (optional) a dictionary of replacement strings (see below)\n\
      Defaults to {}.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n\
   The optional replacements dict can be used to do string substitution of abbreviations \n\
   in the input SMILES. The set of substitutions is repeatedly looped through until \n\
   the string no longer changes. It is the responsibility of the caller to make sure \n\
   that substitutions results in legal and sensible SMILES. \n\
\n\
   Examples of replacements: \n\
\n\
     CC{Q}C with {'{Q}':'OCCO'} -> CCOCCOC  \n\n\
     C{A}C{Q}C with {'{Q}':'OCCO', '{A}':'C1(CC1)'} -> CC1(CC1)COCCOC  \n\n\
     C{A}C{Q}C with {'{Q}':'{X}CC{X}', '{A}':'C1CC1', '{X}':'N'} -> CC1CC1CNCCNC  \n\n\
\n";
  python::def("MolFromSmiles", RDKit::MolFromSmiles,
              (python::arg("SMILES"), python::arg("sanitize") = true,
               python::arg("replacements") = python::dict()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString = "Construct an atom from a SMILES string";
  python::def("AtomFromSmiles", SmilesToAtom, python::arg("SMILES"),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString = "Construct a bond from a SMILES string";
  python::def("BondFromSmiles", SmilesToBond, python::arg("SMILES"),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a SMARTS string.\n\n\
  ARGUMENTS:\n\
\n\
    - SMARTS: the smarts string\n\
\n\
    - mergeHs: (optional) toggles the merging of explicit Hs in the query into the attached\n\
      atoms.  So, for example, 'C[H]' becomes '[C;!H0]'.\n\
      Defaults to 0.\n\
\n\
    - replacements: (optional) a dictionary of replacement strings (see below)\n\
      Defaults to {}. See the documentation for MolFromSmiles for an explanation.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromSmarts", RDKit::MolFromSmarts,
              (python::arg("SMARTS"), python::arg("mergeHs") = false,
               python::arg("replacements") = python::dict()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString = "Construct an atom from a SMARTS string";
  python::def("AtomFromSmarts", SmartsToAtom, python::arg("SMARTS"),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString = "Construct a bond from a SMARTS string";
  python::def("BondFromSmarts", SmartsToBond, python::arg("SMILES"),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a SMARTS string.\n\n\
     ARGUMENTS:\n\
\n\
       - SMARTS: the smarts string\n\
\n\
       - params: used to provide optional parameters for the SMARTS parsing\n\
\n\
     RETURNS:\n\
\n\
       a Mol object, None on failure.\n\
\n";
  python::def("MolFromSmarts", MolFromSmartsHelper,
              (python::arg("SMARTS"), python::arg("params")), docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  python::class_<RDKit::SmilesWriteParams, boost::noncopyable>(
      "SmilesWriteParams", "Parameters controlling SMILES writing")
      .def_readwrite("doIsomericSmiles",
                     &RDKit::SmilesWriteParams::doIsomericSmiles,
                     "include stereochemistry and isotope information")
      .def_readwrite(
          "doKekule", &RDKit::SmilesWriteParams::doKekule,
          "kekulize the molecule before generating the SMILES and output "
          "single/double bonds. NOTE that the output is not canonical and that "
          "this will thrown an exception if the molecule cannot be kekulized")
      .def_readwrite("canonical", &RDKit::SmilesWriteParams::canonical,
                     "generate canonical SMILES")
      .def_readwrite(
          "cleanStereo", &RDKit::SmilesWriteParams::cleanStereo,
          "chiral centers are removed if they have duplicate sidechains")
      .def_readwrite("allBondsExplicit",
                     &RDKit::SmilesWriteParams::allBondsExplicit,
                     "include symbols for all bonds")
      .def_readwrite("allHsExplicit", &RDKit::SmilesWriteParams::allHsExplicit,
                     "provide hydrogen counts for every atom")
      .def_readwrite(
          "doRandom", &RDKit::SmilesWriteParams::doRandom,
          "randomize the output order. The resulting SMILES is not canonical")
      .def_readwrite("rootedAtAtom", &RDKit::SmilesWriteParams::rootedAtAtom,
                     "make sure the SMILES starts at the specified atom. The "
                     "resulting SMILES is not canonical")
      .def_readwrite(
          "includeDativeBonds", &RDKit::SmilesWriteParams::includeDativeBonds,
          "include the RDKit extension for dative bonds. Otherwise dative bonds will be written as single bonds")
      .def_readwrite("ignoreAtomMapNumbers",
                     &RDKit::SmilesWriteParams::ignoreAtomMapNumbers,
                     "ignore atom map numbers when canonicalizing the molecule")
      .def("__setattr__", &safeSetattr);

  python::def("MolToSmiles",
              (std::string(*)(const ROMol &,
                              const SmilesWriteParams &))RDKit::MolToSmiles,
              (python::arg("mol"), python::arg("params")),
              "Returns the canonical SMILES string for a molecule");

  docString =
      "Returns the canonical SMILES string for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - isomericSmiles: (optional) include information about stereochemistry in\n\
      the SMILES.  Defaults to true.\n\
    - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in\n\
      the SMILES.  Defaults to false.\n\
    - rootedAtAtom: (optional) if non-negative, this forces the SMILES \n\
      to start at a particular atom. Defaults to -1.  If not -1, overrides\n\
      canonical setting.\n\
    - canonical: (optional) if false no attempt will be made to canonicalize\n\
      the molecule. Defaults to true.\n\
    - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
    - allHsExplicit: (optional) if true, all H counts will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
    - doRandom: (optional) if true, randomize the traversal of the molecule graph,\n\
      so we can generate random smiles. Defaults to false.  If true, overrides\n\
      canonical setting.\n\
    - ignoreAtomMapNumbers (optional) if true, ignores any atom map numbers when\n\
      canonicalizing the molecule \n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolToSmiles",
      (std::string(*)(const ROMol &, bool, bool, int, bool, bool, bool, bool,
                      bool))RDKit::MolToSmiles,
      (python::arg("mol"), python::arg("isomericSmiles") = true,
       python::arg("kekuleSmiles") = false, python::arg("rootedAtAtom") = -1,
       python::arg("canonical") = true, python::arg("allBondsExplicit") = false,
       python::arg("allHsExplicit") = false, python::arg("doRandom") = false,
       python::arg("ignoreAtomMapNumbers") = false),
      docString.c_str());

  docString =
      "Returns the canonical SMILES string for a fragment of a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - params: the SmilesWriteParams \n\
    - atomsToUse : a list of atoms to include in the fragment\n\
    - bondsToUse : (optional) a list of bonds to include in the fragment\n\
      if not provided, all bonds between the atoms provided\n\
      will be included.\n\
    - atomSymbols : (optional) a list with the symbols to use for the atoms\n\
      in the SMILES. This should have be mol.GetNumAtoms() long.\n\
    - bondSymbols : (optional) a list with the symbols to use for the bonds\n\
      in the SMILES. This should have be mol.GetNumBonds() long.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolFragmentToSmiles", MolFragmentToSmilesHelper1<smilesfrag_gen>,
              (python::arg("mol"), python::arg("params"),
               python::arg("atomsToUse"), python::arg("bondsToUse") = 0,
               python::arg("atomSymbols") = 0, python::arg("bondSymbols") = 0),
              docString.c_str());

  docString =
      "Returns the canonical SMILES string for a fragment of a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - atomsToUse : a list of atoms to include in the fragment\n\
    - bondsToUse : (optional) a list of bonds to include in the fragment\n\
      if not provided, all bonds between the atoms provided\n\
      will be included.\n\
    - atomSymbols : (optional) a list with the symbols to use for the atoms\n\
      in the SMILES. This should have be mol.GetNumAtoms() long.\n\
    - bondSymbols : (optional) a list with the symbols to use for the bonds\n\
      in the SMILES. This should have be mol.GetNumBonds() long.\n\
    - isomericSmiles: (optional) include information about stereochemistry in\n\
      the SMILES.  Defaults to true.\n\
    - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in\n\
      the SMILES.  Defaults to false.\n\
    - rootedAtAtom: (optional) if non-negative, this forces the SMILES \n\
      to start at a particular atom. Defaults to -1.  If not -1, over-rides\n\
      setting for canonical.\n\
    - canonical: (optional) if false no attempt will be made to canonicalize\n\
      the molecule. Defaults to true.\n\
    - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
    - allHsExplicit: (optional) if true, all H counts will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolFragmentToSmiles", MolFragmentToSmilesHelper2<smilesfrag_gen>,
      (python::arg("mol"), python::arg("atomsToUse"),
       python::arg("bondsToUse") = 0, python::arg("atomSymbols") = 0,
       python::arg("bondSymbols") = 0, python::arg("isomericSmiles") = true,
       python::arg("kekuleSmiles") = false, python::arg("rootedAtAtom") = -1,
       python::arg("canonical") = true, python::arg("allBondsExplicit") = false,
       python::arg("allHsExplicit") = false),
      docString.c_str());

  python::enum_<RDKit::SmilesWrite::CXSmilesFields>("CXSmilesFields")
      .value("CX_NONE", RDKit::SmilesWrite::CXSmilesFields::CX_NONE)
      .value("CX_ATOM_LABELS",
             RDKit::SmilesWrite::CXSmilesFields::CX_ATOM_LABELS)
      .value("CX_MOLFILE_VALUES",
             RDKit::SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES)
      .value("CX_COORDS", RDKit::SmilesWrite::CXSmilesFields::CX_COORDS)
      .value("CX_RADICALS", RDKit::SmilesWrite::CXSmilesFields::CX_RADICALS)
      .value("CX_ATOM_PROPS", RDKit::SmilesWrite::CXSmilesFields::CX_ATOM_PROPS)
      .value("CX_LINKNODES", RDKit::SmilesWrite::CXSmilesFields::CX_LINKNODES)
      .value("CX_ENHANCEDSTEREO",
             RDKit::SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO)
      .value("CX_SGROUPS", RDKit::SmilesWrite::CXSmilesFields::CX_SGROUPS)
      .value("CX_POLYMER", RDKit::SmilesWrite::CXSmilesFields::CX_POLYMER)
      .value("CX_BOND_CFG", RDKit::SmilesWrite::CXSmilesFields::CX_BOND_CFG)
      .value("CX_BOND_ATROPISOMER",
             RDKit::SmilesWrite::CXSmilesFields::CX_BOND_ATROPISOMER)
      .value("CX_COORDINATE_BONDS",
             RDKit::SmilesWrite::CXSmilesFields::CX_COORDINATE_BONDS)
      .value("CX_ZERO_BONDS", RDKit::SmilesWrite::CXSmilesFields::CX_ZERO_BONDS)
      .value("CX_ALL", RDKit::SmilesWrite::CXSmilesFields::CX_ALL)
      .value("CX_ALL_BUT_COORDS",
             RDKit::SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS);

  python::enum_<RDKit::v2::FileParsers::SCSRBaseHbondOptions>(
      "SCSRBaseHbondOptions")
      .value("Ignore", RDKit::v2::FileParsers::SCSRBaseHbondOptions::Ignore)
      .value("UseSapAll",
             RDKit::v2::FileParsers::SCSRBaseHbondOptions::UseSapAll)
      .value("UseSapOne",
             RDKit::v2::FileParsers::SCSRBaseHbondOptions::UseSapOne)
      .value("Auto", RDKit::v2::FileParsers::SCSRBaseHbondOptions::Auto);

  python::enum_<RDKit::v2::FileParsers::SCSRTemplateNames>("SCSRTemplateNames")
      .value("UseFirstName",
             RDKit::v2::FileParsers::SCSRTemplateNames::UseFirstName)
      .value("UseSecondName",
             RDKit::v2::FileParsers::SCSRTemplateNames::UseSecondName)
      .value("AsEntered", RDKit::v2::FileParsers::SCSRTemplateNames::AsEntered);

  python::enum_<RDKit::RestoreBondDirOption>("RestoreBondDirOption")
      .value("RestoreBondDirOptionClear",
             RDKit::RestoreBondDirOption::RestoreBondDirOptionClear)
      .value("RestoreBondDirOptionTrue",
             RDKit::RestoreBondDirOption::RestoreBondDirOptionTrue);

  python::def(
      "MolToCXSmiles",
      (std::string(*)(const ROMol &, const SmilesWriteParams &, std::uint32_t,
                      RestoreBondDirOption))RDKit::MolToCXSmiles,
      (python::arg("mol"), python::arg("params"),
       python::arg("flags") = RDKit::SmilesWrite::CXSmilesFields::CX_ALL,
       python::arg("restoreBondDirs") =
           RDKit::RestoreBondDirOption::RestoreBondDirOptionClear),
      "Returns the CXSMILES string for a molecule");

  docString =
      "Returns the CXSMILES string for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - isomericSmiles: (optional) include information about stereochemistry in\n\
      the SMILES.  Defaults to true.\n\
    - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in\n\
      the SMILES.  Defaults to false.\n\
    - rootedAtAtom: (optional) if non-negative, this forces the SMILES \n\
      to start at a particular atom. Defaults to -1.\n\
    - canonical: (optional) if false no attempt will be made to canonicalize\n\
      the molecule. Defaults to true.\n\
    - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
    - allHsExplicit: (optional) if true, all H counts will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
    - doRandom: (optional) if true, randomizes the traversal of the molecule graph,\n\
      so we can generate random smiles. Defaults to false.  If true, overrides\n\
      canonical setting.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolToCXSmiles",
      (std::string(*)(const ROMol &, bool, bool, int, bool, bool, bool,
                      bool))RDKit::MolToCXSmiles,
      (python::arg("mol"), python::arg("isomericSmiles") = true,
       python::arg("kekuleSmiles") = false, python::arg("rootedAtAtom") = -1,
       python::arg("canonical") = true, python::arg("allBondsExplicit") = false,
       python::arg("allHsExplicit") = false, python::arg("doRandom") = false),
      docString.c_str());

  docString =
      "Returns the CXSMILES string for a fragment of a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - params: the SmilesWriteParams \n\
    - atomsToUse : a list of atoms to include in the fragment\n\
    - bondsToUse : (optional) a list of bonds to include in the fragment\n\
      if not provided, all bonds between the atoms provided\n\
      will be included.\n\
    - atomSymbols : (optional) a list with the symbols to use for the atoms\n\
      in the SMILES. This should have be mol.GetNumAtoms() long.\n\
    - bondSymbols : (optional) a list with the symbols to use for the bonds\n\
      in the SMILES. This should have be mol.GetNumBonds() long.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolFragmentToCXSmiles",
              MolFragmentToSmilesHelper1<cxsmilesfrag_gen>,
              (python::arg("mol"), python::arg("params"),
               python::arg("atomsToUse"), python::arg("bondsToUse") = 0,
               python::arg("atomSymbols") = 0, python::arg("bondSymbols") = 0),
              docString.c_str());

  docString =
      "Returns the CXSMILES string for a fragment of a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - atomsToUse : a list of atoms to include in the fragment\n\
    - bondsToUse : (optional) a list of bonds to include in the fragment\n\
      if not provided, all bonds between the atoms provided\n\
      will be included.\n\
    - atomSymbols : (optional) a list with the symbols to use for the atoms\n\
      in the SMILES. This should have be mol.GetNumAtoms() long.\n\
    - bondSymbols : (optional) a list with the symbols to use for the bonds\n\
      in the SMILES. This should have be mol.GetNumBonds() long.\n\
    - isomericSmiles: (optional) include information about stereochemistry in\n\
      the SMILES.  Defaults to true.\n\
    - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in\n\
      the SMILES.  Defaults to false.\n\
    - rootedAtAtom: (optional) if non-negative, this forces the SMILES \n\
      to start at a particular atom. Defaults to -1.  If not -1, overrides\n\
      canonical setting.\n\
    - canonical: (optional) if false no attempt will be made to canonicalize\n\
      the molecule. Defaults to true.\n\
    - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
    - allHsExplicit: (optional) if true, all H counts will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolFragmentToCXSmiles", MolFragmentToSmilesHelper2<cxsmilesfrag_gen>,
      (python::arg("mol"), python::arg("atomsToUse"),
       python::arg("bondsToUse") = 0, python::arg("atomSymbols") = 0,
       python::arg("bondSymbols") = 0, python::arg("isomericSmiles") = true,
       python::arg("kekuleSmiles") = false, python::arg("rootedAtAtom") = -1,
       python::arg("canonical") = true, python::arg("allBondsExplicit") = false,
       python::arg("allHsExplicit") = false),
      docString.c_str());

  docString =
      "Returns a SMARTS string for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - isomericSmiles: (optional) include information about stereochemistry in\n\
      the SMARTS.  Defaults to true.\n\
    - rootedAtomAtom: (optional) the atom index to start the SMARTS from.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToSmarts",
              (std::string(*)(const ROMol &, bool, int))RDKit::MolToSmarts,
              (python::arg("mol"), python::arg("isomericSmiles") = true,
               python::arg("rootedAtAtom") = -1),
              docString.c_str());

  docString =
      "Returns a SMARTS string for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - params: SmilesWriteParams controlling the SMARTS generation\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToSmarts",
              (std::string(*)(const ROMol &,
                              const SmilesWriteParams &))RDKit::MolToSmarts,
              (python::arg("mol"), python::arg("params")), docString.c_str());

  docString =
      "Returns a SMARTS string for a fragment of a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - atomsToUse: indices of atoms to include in the SMARTS string\n\
    - bondsToUse: indices of bonds to include in the SMARTS string (optional)\n\
    - isomericSmarts: (optional) include information about stereochemistry in\n\
      the SMARTS.  Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolFragmentToSmarts", molFragmentToSmarts,
      (python::arg("mol"), python::arg("atomsToUse"),
       python::arg("bondsToUse") = 0, python::arg("isomericSmarts") = true),
      docString.c_str());

  docString =
      "Returns a SMARTS string for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - isomericSmiles: (optional) include information about stereochemistry in\n\
      the SMARTS.  Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToCXSmarts",
              (std::string(*)(const ROMol &, bool))RDKit::MolToCXSmarts,
              (python::arg("mol"), python::arg("isomericSmiles") = true),
              docString.c_str());

  docString =
      "Returns a SMARTS string for a fragment of a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - atomsToUse: indices of atoms to include in the SMARTS string\n\
    - bondsToUse: indices of bonds to include in the SMARTS string (optional)\n\
    - isomericSmarts: (optional) include information about stereochemistry in\n\
      the SMARTS.  Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolFragmentToCXSmarts", molFragmentToCXSmarts,
      (python::arg("mol"), python::arg("atomsToUse"),
       python::arg("bondsToUse") = 0, python::arg("isomericSmarts") = true),
      docString.c_str());

  docString =
      "Writes a molecule to a TPL file.\n\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - fileName: name of the file to write\n\
    - partialChargeProp: name of the property to use for partial charges\n\
      Defaults to '_GasteigerCharge'.\n\
    - writeFirstConfTwice: Defaults to False.\n\
      This should be set to True when writing TPLs to be read by \n\
      the CombiCode.\n\
\n";
  python::def("MolToTPLFile", RDKit::MolToTPLFile,
              (python::arg("mol"), python::arg("fileName"),
               python::arg("partialChargeProp") = "_GasteigerCharge",
               python::arg("writeFirstConfTwice") = false),
              docString.c_str());

  docString =
      "Returns the Tpl block for a molecule.\n\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - partialChargeProp: name of the property to use for partial charges\n\
      Defaults to '_GasteigerCharge'.\n\
    - writeFirstConfTwice: Defaults to False.\n\
      This should be set to True when writing TPLs to be read by \n\
      the CombiCode.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToTPLBlock", RDKit::MolToTPLText,
              (python::arg("mol"),
               python::arg("partialChargeProp") = "_GasteigerCharge",
               python::arg("writeFirstConfTwice") = false),
              docString.c_str());

  docString =
      "Construct a molecule from a PDB file.\n\n\
  ARGUMENTS:\n\
\n\
    - pdbFileName: name of the file to read\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to true.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
    - flavor: (optional) \n\
\n\
    - proximityBonding: (optional) toggles automatic proximity bonding\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromPDBFile", RDKit::MolFromPDBFile,
              (python::arg("pdbFileName"), python::arg("sanitize") = true,
               python::arg("removeHs") = true, python::arg("flavor") = 0,
               python::arg("proximityBonding") = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a PDB block.\n\n\
  ARGUMENTS:\n\
\n\
    - molBlock: string containing the PDB block\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
    - flavor: (optional) \n\
\n\
    - proximityBonding: (optional) toggles automatic proximity bonding\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromPDBBlock", RDKit::MolFromPDBBlock,
              (python::arg("molBlock"), python::arg("sanitize") = true,
               python::arg("removeHs") = true, python::arg("flavor") = 0,
               python::arg("proximityBonding") = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Returns a PDB block for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
    - flavor: (optional) \n\
            - flavor & 1 : Write MODEL/ENDMDL lines around each record \n\
            - flavor & 2 : Don't write any CONECT records \n\
            - flavor & 4 : Write CONECT records in both directions \n\
            - flavor & 8 : Don't use multiple CONECTs to encode bond order \n\
            - flavor & 16 : Write MASTER record \n\
            - flavor & 32 : Write TER record \n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToPDBBlock", RDKit::MolToPDBBlock,
              (python::arg("mol"), python::arg("confId") = -1,
               python::arg("flavor") = 0),
              docString.c_str());

  docString =
      "Writes a PDB file for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - filename: name of the file to write\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
    - flavor: (optional) \n\
            - flavor & 1 : Write MODEL/ENDMDL lines around each record \n\
            - flavor & 2 : Don't write any CONECT records \n\
            - flavor & 4 : Write CONECT records in both directions \n\
            - flavor & 8 : Don't use multiple CONECTs to encode bond order \n\
            - flavor & 16 : Write MASTER record \n\
            - flavor & 32 : Write TER record \n\
\n";
  python::def("MolToPDBFile", RDKit::MolToPDBFile,
              (python::arg("mol"), python::arg("filename"),
               python::arg("confId") = -1, python::arg("flavor") = 0),
              docString.c_str());

  docString =
      "Construct a molecule from a sequence string (currently supports standard amino acids, DNA and RNA bases).\n\n\
  ARGUMENTS:\n\
\n\
    - text: string containing the sequence\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - flavor: (optional)\n\
        - 0 Protein, L amino acids (default)\n\
        - 1 Protein, D amino acids\n\
        - 2 RNA, no cap\n\
        - 3 RNA, 5' cap\n\
        - 4 RNA, 3' cap\n\
        - 5 RNA, both caps\n\
        - 6 DNA, no cap\n\
        - 7 DNA, 5' cap\n\
        - 8 DNA, 3' cap\n\
        - 9 DNA, both caps\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromSequence", RDKit::MolFromSequence,
              (python::arg("text"), python::arg("sanitize") = true,
               python::arg("flavor") = 0),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Returns the sequence string for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
\n\
  NOTE: the molecule should contain monomer information in AtomMonomerInfo structures \n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToSequence", RDKit::MolToSequence, (python::arg("mol")),
              docString.c_str());

  docString =
      "Construct a molecule from a FASTA string (currently supports standard amino acids, DNA and RNA bases).\n\n\
  ARGUMENTS:\n\
\n\
    - text: string containing the FASTA\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
- flavor: (optional)\n\
    - 0 Protein, L amino acids (default)\n\
    - 1 Protein, D amino acids\n\
    - 2 RNA, no cap\n\
    - 3 RNA, 5' cap\n\
    - 4 RNA, 3' cap\n\
    - 5 RNA, both caps\n\
    - 6 DNA, no cap\n\
    - 7 DNA, 5' cap\n\
    - 8 DNA, 3' cap\n\
    - 9 DNA, both caps\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromFASTA", RDKit::MolFromFASTA,
              (python::arg("text"), python::arg("sanitize") = true,
               python::arg("flavor") = 0),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Returns the FASTA string for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
\n\
  NOTE: the molecule should contain monomer information in AtomMonomerInfo structures \n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToFASTA", RDKit::MolToFASTA, (python::arg("mol")),
              docString.c_str());

  docString =
      "Construct a molecule from a HELM string (currently supports standard amino acids, DNA and RNA bases).\n\n\
  ARGUMENTS:\n\
\n\
    - text: string containing the HELM\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromHELM", RDKit::MolFromHELM,
              (python::arg("text"), python::arg("sanitize") = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Returns the HELM string for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
\n\
  NOTE: the molecule should contain monomer information in AtomMonomerInfo structures \n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToHELM", RDKit::MolToHELM, (python::arg("mol")),
              docString.c_str());

  docString =
      "Returns the canonical atom ranking for each atom of a molecule fragment.\n\
  If breakTies is False, this returns the symmetry class for each atom.  The symmetry\n\
  class is used by the canonicalization routines to type each atom based on the whole\n\
  chemistry of the molecular graph.  Any atom with the same rank (symmetry class) is\n\
  indistinguishable.  For example:\n\
\n\
    >>> mol = MolFromSmiles('C1NCN1')\n\
    >>> list(CanonicalRankAtoms(mol, breakTies=False))\n\
    [0,1,0,1]\n\
\n\
  In this case the carbons have the same symmetry class and the nitrogens have the same\n\
  symmetry class.  From the perspective of the Molecular Graph, they are identical.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - breakTies: (optional) force breaking of ranked ties [default=True]\n\
    - includeChirality: (optional) use chiral information when computing rank [default=True]\n\
    - includeIsotopes: (optional) use isotope information when computing rank [default=True]\n\
    - includeAtomMaps: (optional) use atom map information when computing rank [default=True]\n\
    - includeChiralPresence: (optional) use information about whether or not chirality is specified when computing rank [default=False]\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("CanonicalRankAtoms", CanonicalRankAtoms,
              (python::arg("mol"), python::arg("breakTies") = true,
               python::arg("includeChirality") = true,
               python::arg("includeIsotopes") = true,
               python::arg("includeAtomMaps") = true,
               python::arg("includeChiralPresence") = false),
              docString.c_str());

  docString =
      "Returns the canonical atom ranking for each atom of a molecule fragment\n\
  See help(CanonicalRankAtoms) for more information.\n\
\n\
   >>> mol = MolFromSmiles('C1NCN1.C1NCN1')\n\
   >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(0,4), breakTies=False))\n\
   [4,6,4,6,-1,-1,-1,-1]\n\
   >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(4,8), breakTies=False))\n\
   [-1,-1,-1,-1,4,6,4,6]\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - atomsToUse : a list of atoms to include in the fragment\n\
    - bondsToUse : (optional) a list of bonds to include in the fragment\n\
      if not provided, no bonds will be used\n\
    - atomSymbols : (optional) a list with the symbols to use for the atoms\n\
      in the SMILES. This should have be mol.GetNumAtoms() long.\n\
    - breakTies: (optional) force breaking of ranked ties\n\
    - includeChirality: (optional) use chiral information when computing rank [default=True]\n\
    - includeIsotopes: (optional) use isotope information when computing rank [default=True]\n\
    - includeAtomMaps: (optional) use atom map information when computing rank [default=True]\n\
    - includeChiralPresence: (optional) use information about whether or not chirality is specified when computing rank [default=False]\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "CanonicalRankAtomsInFragment", CanonicalRankAtomsInFragment,
      (python::arg("mol"), python::arg("atomsToUse"),
       python::arg("bondsToUse") = 0, python::arg("atomSymbols") = 0,
       python::arg("breakTies") = true, python::arg("includeChirality") = true,
       python::arg("includeIsotopes") = true,
       python::arg("includeAtomMaps") = true,
       python::arg("includeChiralPresence") = false),
      docString.c_str());

  python::def("CanonicalizeEnhancedStereo", CanonicalizeEnhancedStereo,
              (python::arg("mol")));

  python::def(
      "CreateAtomIntPropertyList", FileParserUtils::createAtomIntPropertyList,
      (python::arg("mol"), python::arg("propName"),
       python::arg("missingValueMarker") = "", python::arg("lineSize") = 190),
      "creates a list property on the molecule from individual atom property "
      "values");
  python::def(
      "CreateAtomDoublePropertyList",
      FileParserUtils::createAtomDoublePropertyList,
      (python::arg("mol"), python::arg("propName"),
       python::arg("missingValueMarker") = "", python::arg("lineSize") = 190),
      "creates a list property on the molecule from individual atom property "
      "values");
  python::def(
      "CreateAtomBoolPropertyList", FileParserUtils::createAtomBoolPropertyList,
      (python::arg("mol"), python::arg("propName"),
       python::arg("missingValueMarker") = "", python::arg("lineSize") = 190),
      "creates a list property on the molecule from individual atom property "
      "values");
  python::def(
      "CreateAtomStringPropertyList",
      FileParserUtils::createAtomStringPropertyList,
      (python::arg("mol"), python::arg("propName"),
       python::arg("missingValueMarker") = "", python::arg("lineSize") = 190),
      "creates a list property on the molecule from individual atom property "
      "values");

  python::def(
      "CreateBondIntPropertyList", FileParserUtils::createBondIntPropertyList,
      (python::arg("mol"), python::arg("propName"),
       python::arg("missingValueMarker") = "", python::arg("lineSize") = 190),
      "creates a list property on the molecule from individual bond property "
      "values");
  python::def(
      "CreateBondDoublePropertyList",
      FileParserUtils::createBondDoublePropertyList,
      (python::arg("mol"), python::arg("propName"),
       python::arg("missingValueMarker") = "", python::arg("lineSize") = 190),
      "creates a list property on the molecule from individual bond property "
      "values");
  python::def(
      "CreateBondBoolPropertyList", FileParserUtils::createBondBoolPropertyList,
      (python::arg("mol"), python::arg("propName"),
       python::arg("missingValueMarker") = "", python::arg("lineSize") = 190),
      "creates a list property on the molecule from individual bond property "
      "values");
  python::def(
      "CreateBondStringPropertyList",
      FileParserUtils::createBondStringPropertyList,
      (python::arg("mol"), python::arg("propName"),
       python::arg("missingValueMarker") = "", python::arg("lineSize") = 190),
      "creates a list property on the molecule from individual bond property "
      "values");

  python::def(
      "MolToRandomSmilesVect", RDKit::MolToRandomSmilesHelper,
      (python::arg("mol"), python::arg("numSmiles"),
       python::arg("randomSeed") = 0, python::arg("isomericSmiles") = true,
       python::arg("kekuleSmiles") = false,
       python::arg("allBondsExplicit") = false,
       python::arg("allHsExplicit") = false),
      "returns a list of SMILES generated using the randomSmiles algorithm");

#ifdef RDK_USE_BOOST_IOSTREAMS
  python::class_<RDKit::PNGMetadataParams, boost::noncopyable>(
      "PNGMetadataParams",
      "Parameters controlling metadata included in PNG images")
      .def_readwrite("includePkl", &RDKit::PNGMetadataParams::includePkl,
                     "toggles inclusion of molecule pickle (default=True)")
      .def_readwrite("includeSmiles", &RDKit::PNGMetadataParams::includeSmiles,
                     "toggles inclusion of molecule CXSMILES (default=True)")
      .def_readwrite("includeMol", &RDKit::PNGMetadataParams::includeMol,
                     "toggles inclusion of molecule molblock (default=False)")
      .def_readwrite(
          "propertyFlags", &RDKit::PNGMetadataParams::propertyFlags,
          "choose properties to be included in the pickle (default=rdkit.Chem.rdchem.PropertyPickleOptions.NoProps)")
      .def_readwrite(
          "smilesWriteParams", &RDKit::PNGMetadataParams::smilesWriteParams,
          "choose SmilesWriteParams for the CXSMILES string (default=rdkit.Chem.rdmolfiles.SmilesWriteParams())")
      .def_readwrite(
          "cxSmilesFlags", &RDKit::PNGMetadataParams::cxSmilesFlags,
          "choose CXSMILES fields to be included in the CXSMILES string (default=rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL)")
      .def_readwrite(
          "restoreBondDirs", &RDKit::PNGMetadataParams::restoreBondDirs,
          "choose what to do with bond dirs in the CXSMILES string (default=rdkit.Chem.rdmolfiles.RestoreBondDirOption.RestoreBondDirOptionClear)")
      .def("__setattr__", &safeSetattr);

  docString =
      R"DOC(Construct a molecule from metadata in a PNG string.

     ARGUMENTS:

       - png: the PNG string

       - params: used to provide optional parameters for the metadata parsing

     RETURNS:
       a Mol object, None on failure.
  )DOC";
  python::def("MolFromPNGString", MolFromPNGString,
              (python::arg("png"), python::arg("params") = python::object()),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      R"DOC(Construct a molecule from metadata in a PNG file.

     ARGUMENTS:

       - filename: the PNG filename

       - params: used to provide optional parameters for the metadata parsing

     RETURNS:
       a Mol object, None on failure.)DOC";

  std::string cdxml_notes = R"DOC()DOC";

  python::def(
      "MolFromPNGFile", MolFromPNGFile,
      (python::arg("filename"), python::arg("params") = python::object()),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());

  python::def("MolsFromPNGString", MolsFromPNGString,
              (python::arg("png"), python::arg("tag") = PNGData::pklTag,
               python::arg("params") = python::object()),
              "returns a tuple of molecules constructed from the PNG string");
  python::def("MolsFromPNGFile", MolsFromPNGFile,
              (python::arg("filename"), python::arg("tag") = PNGData::pklTag,
               python::arg("params") = python::object()),
              "returns a tuple of molecules constructed from the PNG file");
#endif

  docString =
      R"DOC(Construct a molecule from a cdxml file.

     Note that the CDXML format is large and complex, the RDKit doesn't support
     full functionality, just the base ones required for molecule and
     reaction parsing.

     ARGUMENTS:

       - filename: the cdxml filename

       - sanitize: if True, sanitize the molecules [default True]

       - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

     RETURNS:
       an iterator of parsed Mol objects.)DOC";

  python::def("MolsFromCDXMLFile", MolsFromCDXMLFile,
              (python::arg("filename"), python::arg("sanitize") = true,
               python::arg("removeHs") = true),
              docString.c_str());

  docString =
      R"DOC(Construct a molecule from a cdxml string.

     Note that the CDXML format is large and complex, the RDKit doesn't support
     full functionality, just the base ones required for molecule and
     reaction parsing.

     ARGUMENTS:

       - cdxml: the cdxml string

       - sanitize: if True, sanitize the molecules [default True]

       - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

     RETURNS:
       an iterator of parsed Mol objects.)DOC";

  python::def("MolsFromCDXML", MolsFromCDXML,
              (python::arg("cdxml"), python::arg("sanitize") = true,
               python::arg("removeHs") = true),
              docString.c_str());

  python::enum_<RDKit::v2::CDXMLParser::CDXMLFormat>("CDXMLFormat")
      .value("CDXML", RDKit::v2::CDXMLParser::CDXMLFormat::CDXML)
      .value("CDX", RDKit::v2::CDXMLParser::CDXMLFormat::CDX)
      .value("Auto", RDKit::v2::CDXMLParser::CDXMLFormat::Auto);

  python::class_<RDKit::v2::CDXMLParser::CDXMLParserParams, boost::noncopyable>(
      "CDXMLParserParams",
      "Parameters controlling conversion of a CDXML document to molecules",
      python::init<>(python::args("self"), "Construct a default CDXMLFormat"))
      .def(python::init<bool, bool, RDKit::v2::CDXMLParser::CDXMLFormat>(
          python::args("self", "sanitize", "removeHs", "format")))
      .def_readwrite("sanitize",
                     &RDKit::v2::CDXMLParser::CDXMLParserParams::sanitize,
                     "controls whether or not the molecule is sanitized before "
                     "being returned")
      .def_readwrite("removeHs",
                     &RDKit::v2::CDXMLParser::CDXMLParserParams::removeHs,
                     "controls whether or not Hs are removed before the "
                     "molecule is returned")
      .def_readwrite(
          "format", &RDKit::v2::CDXMLParser::CDXMLParserParams::format,
          "ChemDraw format One of Auto, CDXML, CDX.  For data streams, Auto defaults to CDXML")
      .def("__setattr__", &safeSetattr);

  docString =
      R"DOC(Construct a molecule from a cdxml file.

     Note: that the CDXML format is large and complex, the RDKit doesn't support
     full functionality, just the base ones required for molecule and
     reaction parsing.

     Note: If the ChemDraw extensions are available,
        CDXMLFormat::Auto attempts to see if the input string is CDXML or CDX,
     If not, it defaults to CDXML

     ARGUMENTS:

       - filename: the cdxml filename

       - pyParams: CDXParserParams, see CDXParserParams for usage

     RETURNS:
       a tuple  of parsed Mol objects.)DOC";

  python::def("MolsFromCDXMLFile", MolsFromCDXMLFileHelper,
              (python::arg("filename"), python::arg("params")),
              docString.c_str());

  docString =
      R"DOC(Construct a molecule from a cdxml string.

     Note that the CDXML format is large and complex, the RDKit doesn't support
     full functionality, just the base ones required for molecule and
     reaction parsing.

     Note: in this function CDXMLFormat::Auto currently defaults to CDXML

     ARGUMENTS:

       - cdxml: the cdxml string

       - pyParams: CDXParserParams, see CDXParserParams for usage

     RETURNS:
       a tuple of parsed Mol objects.)DOC";

  python::def("MolsFromCDXML", MolsFromCDXMLHelper,
              (python::arg("cdxml"), python::arg("params")), docString.c_str());

  docString = "Returns true if the RDKit is built with ChemDraw CDX support";
  python::def("HasChemDrawCDXSupport",
              RDKit::v2::CDXMLParser::hasChemDrawCDXSupport, docString.c_str());

#ifdef RDK_USE_BOOST_IOSTREAMS
  docString =
      R"DOC(Adds molecular metadata to PNG data read from a file.

     ARGUMENTS:

       - mol: the molecule

       - filename: the PNG filename

       - includePkl: include the RDKit's internal binary format in the output

       - includeSmiles: include CXSmiles in the output

       - includeMol: include CTAB (Mol) in the output

     RETURNS:
       the updated PNG data)DOC";
  python::def(
      "MolMetadataToPNGFile", addMolToPNGFileHelper,
      (python::arg("mol"), python::arg("filename"),
       python::arg("includePkl") = true, python::arg("includeSmiles") = true,
       python::arg("includeMol") = false),
      docString.c_str());

  docString =
      R"DOC(Adds molecular metadata to PNG data read from a file.

     ARGUMENTS:

       - mol: the molecule

       - filename: the PNG filename

       - params: an instance of PNGMetadataParams

     RETURNS:
       the updated PNG data)DOC";
  python::def(
      "MolMetadataToPNGFile", addMolToPNGFileHelperParams,
      (python::arg("mol"), python::arg("filename"), python::arg("params")),
      docString.c_str());

  docString =
      R"DOC(Adds molecular metadata to a PNG string.

     ARGUMENTS:

       - mol: the molecule

       - png: the PNG string

       - includePkl: include the RDKit's internal binary format in the output

       - includeSmiles: include CXSmiles in the output

       - includeMol: include CTAB (Mol) in the output

     RETURNS:
       the updated PNG data)DOC";
  python::def(
      "MolMetadataToPNGString", addMolToPNGStringHelper,
      (python::arg("mol"), python::arg("png"), python::arg("includePkl") = true,
       python::arg("includeSmiles") = true, python::arg("includeMol") = false),
      docString.c_str());

  docString =
      R"DOC(Adds molecular metadata to a PNG string.

     ARGUMENTS:

       - mol: the molecule

       - png: the PNG string

       - params: an instance of PNGMetadataParams

     RETURNS:
       the updated PNG data)DOC";
  python::def("MolMetadataToPNGString", addMolToPNGStringHelperParams,
              (python::arg("mol"), python::arg("png"), python::arg("params")),
              docString.c_str());

  docString =
      R"DOC(Adds metadata to PNG data read from a file.

     ARGUMENTS:

       - metadata: dict with the metadata to be written
                   (keys and values should be strings)

       - filename: the PNG filename

     RETURNS:
       the updated PNG data)DOC";
  python::def("AddMetadataToPNGFile", addMetadataToPNGFileHelper,
              (python::arg("metadata"), python::arg("filename")),
              docString.c_str());

  docString =
      R"DOC(Adds metadata to a PNG string.

     ARGUMENTS:

       - metadata: dict with the metadata to be written
                   (keys and values should be strings)

       - png: the PNG string

     RETURNS:
       the updated PNG data)DOC";
  python::def("AddMetadataToPNGString", addMetadataToPNGStringHelper,
              (python::arg("metadata"), python::arg("png")), docString.c_str());

  python::def("MetadataFromPNGFile", MetadataFromPNGFile,
              (python::arg("filename"), python::arg("asList") = false),
              "Returns a dict with all metadata from the PNG file. Keys are "
              "strings, values are bytes. "
              "If asList is True, a list of (key, value) tuples is returned; "
              "this enables retrieving multiple values sharing the same key.");

  python::def("MetadataFromPNGString", MetadataFromPNGString,
              (python::arg("png"), python::arg("asList") = false),
              "Returns a dict with all metadata from the PNG string. Keys are "
              "strings, values are bytes. "
              "If asList is True, a list of (key, value) tuples is returned; "
              "this enables retrieving multiple values sharing the same key.");
#endif
/********************************************************
 * MolSupplier stuff
 *******************************************************/
#ifdef SUPPORT_COMPRESSED_SUPPLIERS
  wrap_compressedsdsupplier();
#endif
  wrap_sdsupplier();
  wrap_forwardsdsupplier();
  wrap_tdtsupplier();
  wrap_smisupplier();
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
  wrap_maesupplier();
#endif
  // wrap_pdbsupplier();

  /********************************************************
   * MolWriter stuff
   *******************************************************/
  wrap_smiwriter();
  wrap_sdwriter();
  wrap_tdtwriter();
  wrap_pdbwriter();
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
  wrap_maewriter();
#endif

#ifdef RDK_BUILD_THREADSAFE_SSS
  /********************************************************
   * MultithreadedMolWriter stuff
   *******************************************************/
  wrap_multiSmiSupplier();
  wrap_multiSDSupplier();
#endif
}
