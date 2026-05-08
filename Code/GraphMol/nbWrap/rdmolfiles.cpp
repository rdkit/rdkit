//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/dynamic_bitset.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/pair.h>

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

#include <RDBoost/Wrap_nb.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SanitException.h>
#include <string.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;

namespace RDKit {
std::string pyObjectToString(nb::object input) {
  if (nb::isinstance<nb::str>(input)) {
    return nb::cast<std::string>(input);
  }
  std::wstring ws = nb::cast<std::wstring>(input);
  return std::string(ws.begin(), ws.end());
}

ROMol *MolFromSmiles(nb::object ismiles, bool sanitize, nb::dict replDict) {
  std::map<std::string, std::string> replacements;
  const auto items = replDict.items();
  for (unsigned int i = 0; i < nb::len(items); ++i) {
    const auto item = items[i];
    replacements[nb::cast<std::string>(item[0])] =
        nb::cast<std::string>(item[1]);
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

ROMol *MolFromSmarts(nb::object ismarts, bool mergeHs, nb::dict replDict) {
  std::map<std::string, std::string> replacements;
  const auto items = replDict.items();
  for (unsigned int i = 0; i < nb::len(items); ++i) {
    const auto item = items[i];
    replacements[nb::cast<std::string>(item[0])] =
        nb::cast<std::string>(item[1]);
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
    throw nb::python_error();
  } catch (...) {
    newM = nullptr;
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromTPLBlock(nb::object itplBlock, bool sanitize = true,
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
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromMolBlock(nb::object imolBlock, bool sanitize, bool removeHs,
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
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

RDKit::ROMol *MolFromSCSRBlock(
    const std::string &molBlock, bool sanitize, bool removeHs,
    RDKit::v2::FileParsers::MolFromSCSRParams *pyparams) {
  RDKit::v2::FileParsers::MolFromSCSRParams scsrParams;
  if (pyparams) {
    scsrParams = *pyparams;
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

RDKit::ROMol *MolFromSCSRFile(
    const std::string &molFilename, bool sanitize, bool removeHs,
    RDKit::v2::FileParsers::MolFromSCSRParams *pyparams) {
  RDKit::v2::FileParsers::MolFromSCSRParams scsrParams;
  if (pyparams) {
    scsrParams = *pyparams;
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
    throw nb::python_error();
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
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromMrvBlock(nb::object imolBlock, bool sanitize, bool removeHs) {
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

ROMol *MolFromXYZBlock(nb::object ixyzBlock) {
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

ROMol *MolFromSVG(nb::object imolBlock, bool sanitize, bool removeHs) {
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
    throw nb::python_error();
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
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromPDBBlock(nb::object molBlock, bool sanitize, bool removeHs,
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

ROMol *MolFromSequence(nb::object seq, bool sanitize, int flavor) {
  RWMol *newM = nullptr;
  try {
    newM = SequenceToMol(pyObjectToString(seq), sanitize, flavor);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}
ROMol *MolFromFASTA(nb::object seq, bool sanitize, int flavor) {
  RWMol *newM = nullptr;
  try {
    newM = FASTAToMol(pyObjectToString(seq), sanitize, flavor);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}
ROMol *MolFromHELM(nb::object seq, bool sanitize) {
  RWMol *newM = nullptr;
  try {
    newM = HELMToMol(pyObjectToString(seq), sanitize);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

std::string molFragmentToSmarts(const ROMol &mol, nb::object atomsToUse,
                                nb::object bondsToUse,
                                bool doIsomericSmarts = true) {
  auto atomIndices =
      pythonObjectToVect(atomsToUse, static_cast<int>(mol.getNumAtoms()));
  if (!atomIndices) {
    throw ValueErrorException("atomsToUse argument must be non-empty");
  }
  auto bondIndices =
      pythonObjectToVect(bondsToUse, static_cast<int>(mol.getNumBonds()));
  return RDKit::MolFragmentToSmarts(mol, *atomIndices, bondIndices.get(),
                                    doIsomericSmarts);
}

std::string molFragmentToCXSmarts(const ROMol &mol, nb::object atomsToUse,
                                  nb::object bondsToUse,
                                  bool doIsomericSmarts = true) {
  auto atomIndices =
      pythonObjectToVect(atomsToUse, static_cast<int>(mol.getNumAtoms()));
  if (!atomIndices) {
    throw ValueErrorException("atomsToUse argument must be non-empty");
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
std::string MolFragmentToSmilesHelper1(
    const ROMol &mol, const SmilesWriteParams &params, nb::object atomsToUse,
    nb::object bondsToUse, nb::object atomSymbols, nb::object bondSymbols) {
  auto avect =
      pythonObjectToVect(atomsToUse, static_cast<int>(mol.getNumAtoms()));
  if (!avect.get() || !(avect->size())) {
    throw ValueErrorException("atomsToUse must not be empty");
  }
  auto bvect =
      pythonObjectToVect(bondsToUse, static_cast<int>(mol.getNumBonds()));
  std::unique_ptr<std::vector<std::string>> asymbols =
      pythonObjectToVect<std::string>(atomSymbols);
  std::unique_ptr<std::vector<std::string>> bsymbols =
      pythonObjectToVect<std::string>(bondSymbols);
  if (asymbols.get() && asymbols->size() != mol.getNumAtoms()) {
    throw ValueErrorException("length of atom symbol list != number of atoms");
  }
  if (bsymbols.get() && bsymbols->size() != mol.getNumBonds()) {
    throw ValueErrorException("length of bond symbol list != number of bonds");
  }

  std::string res = F()(mol, params, *avect.get(), bvect.get(), asymbols.get(),
                        bsymbols.get());
  return res;
}

template <typename F>
std::string MolFragmentToSmilesHelper2(
    const ROMol &mol, nb::object atomsToUse, nb::object bondsToUse,
    nb::object atomSymbols, nb::object bondSymbols, bool doIsomericSmiles,
    bool doKekule, int rootedAtAtom, bool canonical, bool allBondsExplicit,
    bool allHsExplicit) {
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
    const ROMol &mol, nb::object atomsToUse, nb::object bondsToUse,
    nb::object atomSymbols, bool breakTies = true, bool includeChirality = true,
    bool includeIsotopes = true, bool includeAtomMaps = true,
    bool includeChiralPresence = false) {
  std::unique_ptr<std::vector<int>> avect =
      pythonObjectToVect(atomsToUse, static_cast<int>(mol.getNumAtoms()));
  if (!avect.get() || !(avect->size())) {
    throw ValueErrorException("atomsToUse must not be empty");
  }
  std::unique_ptr<std::vector<int>> bvect =
      pythonObjectToVect(bondsToUse, static_cast<int>(mol.getNumBonds()));
  std::unique_ptr<std::vector<std::string>> asymbols =
      pythonObjectToVect<std::string>(atomSymbols);
  if (asymbols.get() && asymbols->size() != mol.getNumAtoms()) {
    throw ValueErrorException("length of atom symbol list != number of atoms");
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
                           includeChirality, includeIsotopes, includeAtomMaps,
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

ROMol *MolFromSmilesHelper(nb::object ismiles,
                           const SmilesParserParams &params) {
  std::string smiles = pyObjectToString(ismiles);

  try {
    return SmilesToMol(smiles, params);
  } catch (...) {
    return nullptr;
  }
}

ROMol *MolFromSmartsHelper(nb::object ismiles,
                           const SmartsParserParams &params) {
  std::string smiles = pyObjectToString(ismiles);

  try {
    return SmartsToMol(smiles, params);
  } catch (...) {
    return nullptr;
  }
}

nb::list MolToRandomSmilesHelper(const ROMol &mol, unsigned int numSmiles,
                                 unsigned int randomSeed, bool doIsomericSmiles,
                                 bool doKekule, bool allBondsExplicit,
                                 bool allHsExplicit) {
  auto res = MolToRandomSmilesVect(mol, numSmiles, randomSeed, doIsomericSmiles,
                                   doKekule, allBondsExplicit, allHsExplicit);
  nb::list pyres;
  for (auto smi : res) {
    pyres.append(smi);
  }
  return pyres;
}

ROMol *MolFromPNGFile(const std::string &filename,
                      SmilesParserParams *pyParams) {
  SmilesParserParams params;
  if (pyParams) {
    params = *pyParams;
  }
  ROMol *newM = nullptr;
  try {
    newM = PNGFileToMol(filename, params);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  return newM;
}

ROMol *MolFromPNGString(nb::object png, SmilesParserParams *pyParams) {
  SmilesParserParams params;
  if (pyParams) {
    params = *pyParams;
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

nb::object addMolToPNGFileHelperParams(const ROMol &mol, nb::object fname,
                                       const PNGMetadataParams &params) {
  std::string cstr = nb::cast<std::string>(fname);

  auto res = addMolToPNGFile(mol, cstr, params);

  nb::object retval = nb::object(nb::steal<nb::object>(
      PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

nb::object addMolToPNGFileHelper(const ROMol &mol, nb::object fname,
                                 bool includePkl, bool includeSmiles,
                                 bool includeMol) {
  PNGMetadataParams params;
  params.includePkl = includePkl;
  params.includeSmiles = includeSmiles;
  params.includeMol = includeMol;
  return addMolToPNGFileHelperParams(mol, fname, params);
}

nb::object addMolToPNGStringHelperParams(const ROMol &mol, nb::object png,
                                         const PNGMetadataParams &params) {
  std::string cstr = nb::cast<std::string>(png);

  auto res = addMolToPNGString(mol, cstr, params);

  nb::object retval = nb::object(nb::steal<nb::object>(
      PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

nb::object addMolToPNGStringHelper(const ROMol &mol, nb::object png,
                                   bool includePkl, bool includeSmiles,
                                   bool includeMol) {
  PNGMetadataParams params;
  params.includePkl = includePkl;
  params.includeSmiles = includeSmiles;
  params.includeMol = includeMol;
  return addMolToPNGStringHelperParams(mol, png, params);
}

nb::object addMetadataToPNGFileHelper(nb::dict pymetadata, nb::object fname) {
  std::string cstr = nb::cast<std::string>(fname);

  std::vector<std::pair<std::string, std::string>> metadata;

  const auto items = pymetadata.items();
  for (unsigned int i = 0; i < nb::len(items); ++i) {
    const auto item = items[i];
    std::string key = nb::cast<std::string>(item[0]);
    std::string val = nb::cast<std::string>(item[1]);
    metadata.push_back(std::make_pair(key, val));
  }

  auto res = addMetadataToPNGFile(cstr, metadata);

  nb::object retval = nb::object(nb::steal<nb::object>(
      PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

nb::object addMetadataToPNGStringHelper(nb::dict pymetadata, nb::object png) {
  std::string cstr = nb::cast<std::string>(png);

  std::vector<std::pair<std::string, std::string>> metadata;
  const auto items = pymetadata.items();
  for (unsigned int i = 0; i < nb::len(items); ++i) {
    const auto item = items[i];
    std::string key = nb::cast<std::string>(item[0]);
    std::string val = nb::cast<std::string>(item[1]);
    metadata.push_back(std::make_pair(key, val));
  }

  auto res = addMetadataToPNGString(cstr, metadata);

  nb::object retval = nb::object(nb::steal<nb::object>(
      PyBytes_FromStringAndSize(res.c_str(), res.length())));
  return retval;
}

nb::object MolsFromPNGFile(const std::string &filename, const std::string &tag,
                           SmilesParserParams *pyParams) {
  SmilesParserParams params;
  if (pyParams) {
    params = *pyParams;
  }
  std::vector<std::unique_ptr<ROMol>> mols;
  try {
    mols = PNGFileToMols(filename, tag, params);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  nb::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(mol.release());
    res.append(sptr);
  }
  return nb::tuple(res);
}

nb::tuple MolsFromPNGString(nb::object png, const std::string &tag,
                            SmilesParserParams *pyParams) {
  SmilesParserParams params;
  if (pyParams) {
    params = *pyParams;
  }
  auto mols = PNGStringToMols(pyObjectToString(png), tag, params);
  nb::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(mol.release());
    res.append(sptr);
  }
  return nb::tuple(res);
}

nb::object MolsFromCDXMLFile(const std::string &filename, bool sanitize,
                             bool removeHs) {
  std::vector<std::unique_ptr<RWMol>> mols;
  try {
    mols = CDXMLFileToMols(filename, sanitize, removeHs);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  nb::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return nb::tuple(res);
}

nb::tuple MolsFromCDXMLHelper(
    nb::object cdxml, RDKit::v2::CDXMLParser::CDXMLParserParams *pyParams) {
  RDKit::v2::CDXMLParser::CDXMLParserParams params;
  if (pyParams) {
    params = *pyParams;
  }
  auto mols =
      RDKit::v2::CDXMLParser::MolsFromCDXML(pyObjectToString(cdxml), params);
  nb::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return nb::tuple(res);
}

nb::object MolsFromCDXMLFileHelper(
    const std::string &filename,
    RDKit::v2::CDXMLParser::CDXMLParserParams *pyParams) {
  RDKit::v2::CDXMLParser::CDXMLParserParams params(
      true, true, RDKit::v2::CDXMLParser::CDXMLFormat::Auto);
  if (pyParams) {
    params = *pyParams;
  }
  std::vector<std::unique_ptr<RWMol>> mols;
  try {
    mols = RDKit::v2::CDXMLParser::MolsFromCDXMLFile(filename, params);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.what());
    throw nb::python_error();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.what() << std::endl;
  } catch (...) {
  }
  nb::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return nb::tuple(res);
}

nb::tuple MolsFromCDXML(nb::object cdxml, bool sanitize, bool removeHs) {
  auto mols = CDXMLToMols(pyObjectToString(cdxml), sanitize, removeHs);
  nb::list res;
  for (auto &mol : mols) {
    // take ownership of the data from the unique_ptr
    ROMOL_SPTR sptr(static_cast<ROMol *>(mol.release()));
    res.append(sptr);
  }
  return nb::tuple(res);
}
namespace {
nb::object translateMetadataToList(
    const std::vector<std::pair<std::string, std::string>> &metadata) {
  nb::list resAsList;
  for (const auto &[key, value] : metadata) {
    // keys are safe to extract:
    // but values may include binary, so we convert them directly to bytes:
    nb::object val = nb::object(nb::steal<nb::object>(
        PyBytes_FromStringAndSize(value.c_str(), value.length())));
    resAsList.append(nb::make_tuple(key, val));
  }
  return resAsList;
}
nb::object translateMetadataToDict(
    const std::vector<std::pair<std::string, std::string>> &metadata) {
  nb::dict resAsDict;
  for (const auto &[key, value] : metadata) {
    // keys are safe to extract:
    // but values may include binary, so we convert them directly to bytes:
    resAsDict[nb::str(key.c_str())] = nb::object(nb::steal<nb::object>(
        PyBytes_FromStringAndSize(value.c_str(), value.length())));
  }
  return resAsDict;
}

}  // namespace
nb::object MetadataFromPNGFile(nb::object fname, bool asList) {
  std::string cstr = nb::cast<std::string>(fname);
  auto metadata = PNGFileToMetadata(cstr);
  if (asList) {
    return translateMetadataToList(metadata);
  }
  return translateMetadataToDict(metadata);
}

nb::object MetadataFromPNGString(nb::object png, bool asList) {
  std::string cstr = nb::cast<std::string>(png);
  auto metadata = PNGStringToMetadata(cstr);
  if (asList) {
    return translateMetadataToList(metadata);
  }
  return translateMetadataToDict(metadata);
}

void CanonicalizeEnhancedStereo(ROMol &mol) {
  Canon::canonicalizeEnhancedStereo(mol);
}

std::string MolToV2KMolBlockHelper(const ROMol &mol,
                                   const MolWriterParams *params, int confId) {
  MolWriterParams localParams;
  if (!params) {
    localParams = *params;
  }
  return MolToV2KMolBlock(mol, localParams, confId);
}

}  // namespace RDKit

// MolSupplier stuff
#ifdef SUPPORT_COMPRESSED_SUPPLIERS
void wrap_compressedsdsupplier(nb::module_ &m);
#endif
void wrap_sdsupplier(nb::module_ &m);
void wrap_forwardsdsupplier(nb::module_ &m);
void wrap_tdtsupplier(nb::module_ &m);
void wrap_smisupplier(nb::module_ &m);
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
void wrap_maesupplier(nb::module_ &m);
#endif

// mol writer stuff
void wrap_smiwriter(nb::module_ &m);
void wrap_sdwriter(nb::module_ &m);
void wrap_tdtwriter(nb::module_ &m);
void wrap_pdbwriter(nb::module_ &m);
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
void wrap_maewriter(nb::module_ &m);
#endif

// MultithreadedMolSupplier stuff
void wrap_multiSmiSupplier(nb::module_ &m);
void wrap_multiSDSupplier(nb::module_ &m);

NB_MODULE(rdmolfiles, m) {
  std::string docString;

  m.doc() =
      "Module containing RDKit functionality for working with molecular file "
      "formats.";
  nb::exception<BadFileException>(m, "BadFileException", PyExc_IOError);
  nb::exception<FileParseException>(m, "FileParseException",
                                    PyExc_RuntimeError);

  docString =
      R"DOC(Construct a molecule from a TPL file.

    ARGUMENTS:

      - fileName: name of the file to read

      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.

      - skipFirstConf: (optional) skips reading the first conformer.
        Defaults to False.
        This should be set to True when reading TPLs written by 
        the CombiCode.

    RETURNS:

      a Mol object, None on failure.

)DOC";
  m.def("MolFromTPLFile", RDKit::MolFromTPLFile, "fileName"_a,
        "sanitize"_a = true, "skipFirstConf"_a = false, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a TPL block.

    ARGUMENTS:

      - tplBlock: string containing the TPL block

      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.

      - skipFirstConf: (optional) skips reading the first conformer.
        Defaults to False.
        This should be set to True when reading TPLs written by 
        the CombiCode.

    RETURNS:

      a Mol object, None on failure.

)DOC";
  m.def("MolFromTPLBlock", RDKit::MolFromTPLBlock, "tplBlock"_a,
        "sanitize"_a = true, "skipFirstConf"_a = false, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a Mol file.

    ARGUMENTS:

      - fileName: name of the file to read

      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to true.

      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.

      - strictParsing: (optional) if this is false, the parser is more lax about.
 correctness of the content.
 Defaults to true.

 RETURNS :

 a Mol object, None on failure.

 )DOC";
  m.def("MolFromMolFile", RDKit::MolFromMolFileHelper, "molFileName"_a,
        "sanitize"_a = true, "removeHs"_a = true, "strictParsing"_a = true,
        docString.c_str(), nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a Mol block.

    ARGUMENTS:

      - molBlock: string containing the Mol block

      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.

      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.

      - strictParsing: (optional) if this is false, the parser is more lax about.
 correctness of the content.
 Defaults to true.

 RETURNS :

 a Mol object, None on failure
          .

 )DOC";
  m.def("MolFromMolBlock", RDKit::MolFromMolBlock, "molBlock"_a,
        "sanitize"_a = true, "removeHs"_a = true, "strictParsing"_a = true,
        docString.c_str(), nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a Marvin (Mrv) file.

    ARGUMENTS:

      - fileName: name of the file to read

      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to true.

      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.

    RETURNS:

      a Mol object, None on failure.

)DOC";
  m.def("MolFromMrvFile", RDKit::MolFromMrvFile, "molFileName"_a,
        "sanitize"_a = true, "removeHs"_a = true, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a Marvin (mrv) block.

    ARGUMENTS:

      - molBlock: string containing the Marvin block

      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.

      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.

    RETURNS:

      a Mol object, None on failure.

)DOC";
  m.def("MolFromMrvBlock", RDKit::MolFromMrvBlock, "mrvBlock"_a,
        "sanitize"_a = true, "removeHs"_a = true, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from an XYZ file.
    ARGUMENTS:
  
      - xyzFileName: name of the file to read
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromXYZFile", RDKit::MolFromXYZFile, "xyzFileName"_a,
        docString.c_str(), nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from an XYZ string.
    ARGUMENTS:
  
      - xyzBlock: the XYZ data to read
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromXYZBlock", RDKit::MolFromXYZBlock, "xyzBlock"_a,
        docString.c_str(), nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from an RDKit-generate SVG string.
    ARGUMENTS:
  
      - svg: string containing the SVG data (must include molecule
      metadata)
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.
  
      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.
  
    RETURNS:
  
      a Mol object, None on failure.
  
    NOTE: this functionality should be considered beta.
  )DOC";
  m.def("MolFromRDKitSVG", RDKit::MolFromSVG, "svg"_a, "sanitize"_a = true,
        "removeHs"_a = true, docString.c_str(), nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a Tripos Mol2 file.
    NOTE:
      The parser expects the atom-typing scheme used by Corina.
      Atom types from Tripos' dbtranslate are less supported.
      Other atom typing schemes are unlikely to work.
  
    ARGUMENTS:                                  \
  
      - mol2FileName: name of the file to read
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to true.
  
      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.
  
      - cleanupSubstructures: (optional) toggles standardizing some 
        substructures found in mol2 files.
        Defaults to true.
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromMol2File", RDKit::MolFromMol2File, "mol2FileName"_a,
        "sanitize"_a = true, "removeHs"_a = true,
        "cleanupSubstructures"_a = true, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a Tripos Mol2 block.
    NOTE:
      The parser expects the atom-typing scheme used by Corina.
      Atom types from Tripos' dbtranslate are less supported.
      Other atom typing schemes are unlikely to work.
  
    ARGUMENTS:
  
      - mol2Block: string containing the Mol2 block
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.
  
      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.
  
      - cleanupSubstructures: (optional) toggles standardizing some
        substructures found in mol2 files.
        Defaults to true.
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromMol2Block", RDKit::MolFromMol2Block, "mol2Block"_a,
        "sanitize"_a = true, "removeHs"_a = true,
        "cleanupSubstructures"_a = true, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a Mol file.
    ARGUMENTS:
  
      - molFileName: name of the file to read
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to true.
  
      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.
  
      - strictParsing: (optional) if this is false, the parser is more lax
      about.
        correctness of the content.
        Defaults to true.
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromMolFile", RDKit::MolFromMolFile, "molFileName"_a,
        "sanitize"_a = true, "removeHs"_a = true, "strictParsing"_a = true,
        docString.c_str(), nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a Mol block.
    ARGUMENTS:
  
      - molBlock: string containing the Mol block
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.
  
      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.
  
      - strictParsing: (optional) if this is false, the parser is more lax
      about.
        correctness of the content.
        Defaults to true.
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromMolBlock", RDKit::MolFromMolBlock, "molBlock"_a,
        "sanitize"_a = true, "removeHs"_a = true, "strictParsing"_a = true,
        docString.c_str(), nb::rv_policy::take_ownership);

  nb::class_<RDKit::MolWriterParams>(m, "MolWriterParams",
                                     "Parameters controlling Mol writing")
      .def_rw("includeStereo", &RDKit::MolWriterParams::includeStereo,
              "toggles inclusion of stereochemistry information (default=True)")
      .def_rw(
          "kekulize", &RDKit::MolWriterParams::kekulize,
          "triggers kekulization of the molecule before it is written (default=True)")
      .def_rw(
          "forceV3000", &RDKit::MolWriterParams::forceV3000,
          "force generation a V3000 mol block (happens automatically with more than 999 atoms or bonds)(default=False)")
      .def_rw(
          "precision", &RDKit::MolWriterParams::precision,
          "precision of coordinates (only available in V3000)(default=false)")
      .def("__setattr__", &safeSetattr);

  // nb::class_<RDKit::v2::FileParsers::MolFromSCSRParams>(
  //     m, "MolFromSCSRParams",
  //     "Parameters controlling conversion of an SCSRMol to a Mol")
  //     .def_rw("includeLeavingGroups",
  //             &RDKit::v2::FileParsers::MolFromSCSRParams::includeLeavingGroups,
  //             "include leaving groups atoms if not substited at that
  //             position")
  //     .def_rw(
  //         "scsrTemplateNames",
  //         &RDKit::v2::FileParsers::MolFromSCSRParams::scsrTemplateNames,
  //         "If True, the first template name in the Sgroup is used as the
  //         Sgroup label")
  //     .def_rw("scsrBaseHbondOptions",
  //             &RDKit::v2::FileParsers::MolFromSCSRParams::scsrBaseHbondOptions,
  //             "One of Ignore, UseSapAll(default) , UseSapOne, Auto")
  //     .def("__setattr__", &safeSetattr);

  // docString =
  //     R"DOC(Construct a molecule from an SCSR Mol block.
  //       ARGUMENTS:

  //         - molBlock: string containing the SCSR Mol block

  //         - sanitize: (optional) toggles sanitization of the molecule.
  //           Defaults to True.

  //         - removeHs: (optional) toggles removing hydrogens from the
  //         molecule.
  //           This only make sense when sanitization is done.
  //           Defaults to true.

  //         - molFromSCSRParams : MolFromSCSRParams to control conversion
  //      RETURNS :
  //      a Mol object, None on failure.
  //      )DOC";
  // m.def("MolFromSCSRBlock", RDKit::MolFromSCSRBlock, "molBlock"_a,
  //       "sanitize"_a = true, "removeHs"_a = true,
  //       "molFromSCSRParams"_a = nb::none(), docString.c_str(),
  //       nb::rv_policy::take_ownership);

  // docString =
  //     R"DOC(Construct a molecule from an SCSR Mol block.
  //       ARGUMENTS:

  //         - filename: string containing the SCSR filename

  //         - sanitize: (optional) toggles sanitization of the molecule.
  //           Defaults to True.

  //         - removeHs: (optional) toggles removing hydrogens from the
  //         molecule.
  //           This only make sense when sanitization is done.
  //           Defaults to true.

  //         - molFromSCSRParams : MolFromSCSRParams to control conversion
  //      RETURNS :
  //      a Mol object, None on failure.
  //      )DOC";
  // m.def("MolFromSCSRFile", RDKit::MolFromSCSRFile, "filename"_a,
  //       "sanitize"_a = true, "removeHs"_a = true,
  //       "molFromSCSRParams"_a = nb::none(), docString.c_str(),
  //       nb::rv_policy::take_ownership);

  docString =
      R"DOC(Returns a Mol block for a molecule
    Arguments:
      - mol: the molecule
      - params: the MolWriterParams
      - confId: (optional) selects which conformation to output (-1 =
      default)

    RETURNS:

      a string
  )DOC";
  m.def("MolToMolBlock",
        (std::string (*)(const ROMol &, const MolWriterParams &,
                         int))RDKit::MolToMolBlock,
        "mol"_a, "params"_a, "confId"_a = -1, docString.c_str());

  docString =
      R"DOC(Returns a Mol block for a molecule
    ARGUMENTS:

      - mol: the molecule
      - includeStereo: (optional) toggles inclusion of stereochemical
        information in the output
      - confId: (optional) selects which conformation to output (-1 =
      default)
      - kekulize: (optional) triggers kekulization of the molecule before
      it's written,
        as suggested by the MDL spec.
      - forceV3000 (optional) force generation a V3000 mol block (happens
      automatically with
        more than 999 atoms or bonds)

    RETURNS:

      a string
  )DOC";
  m.def("MolToMolBlock",
        (std::string (*)(const ROMol &, bool, int, bool,
                         bool))RDKit::MolToMolBlock,
        "mol"_a, "includeStereo"_a = true, "confId"_a = -1, "kekulize"_a = true,
        "forceV3000"_a = false, docString.c_str());

  docString =
      R"DOC(Returns a V3000 Mol block for a molecule
     ARGUMENTS:
   \
       - mol: the molecule
       - params: the MolWriterParams
       - confId: (optional) selects which conformation to output (-1 =
       default)
   \
     RETURNS:
   \
       a string
   )DOC";
  m.def("MolToV3KMolBlock",
        (std::string (*)(const ROMol &, const MolWriterParams &,
                         int))RDKit::MolToV3KMolBlock,
        "mol"_a, "params"_a, "confId"_a = -1, docString.c_str());

  docString =
      R"DOC(Returns a V3000 Mol block for a molecule
     ARGUMENTS:
   \
       - mol: the molecule
       - includeStereo: (optional) toggles inclusion of stereochemical
         information in the output
       - confId: (optional) selects which conformation to output (-1 =
       default)
       - kekulize: (optional) triggers kekulization of the molecule before
       it's written,
         as suggested by the MDL spec.
   \
     RETURNS:
   \
       a string
   )DOC";
  m.def(
      "MolToV3KMolBlock",
      (std::string (*)(const ROMol &, bool, int, bool))RDKit::MolToV3KMolBlock,
      "mol"_a, "includeStereo"_a = true, "confId"_a = -1, "kekulize"_a = true,
      docString.c_str());

  docString =
      R"DOC(Returns a V2000 Mol block for a molecule
     ARGUMENTS:

       - mol: the molecule
       - params: the MolWriterParams
       - confId: (optional) selects which conformation to output (-1 =
       default)

     RETURNS:

       a string

     NOTE: this function throws a ValueError if the molecule has more than
     999 atoms, bonds, or SGroups
  )DOC";
  m.def("MolToV2KMolBlock", MolToV2KMolBlockHelper, "mol"_a,
        "params"_a = nb::none(), "confId"_a = -1, docString.c_str());

  docString =
      R"DOC(Writes a Mol file for a molecule
    ARGUMENTS:

      - mol: the molecule
      - filename: the file to write to
      - params: the MolWriterParams
      - confId: (optional) selects which conformation to output (-1 =
      default)
  )DOC";
  m.def("MolToMolFile",
        (void (*)(const ROMol &, const std::string &, const MolWriterParams &,
                  int))RDKit::MolToMolFile,
        "mol"_a, "filename"_a, "params"_a, "confId"_a = -1, docString.c_str());

  docString =
      R"DOC(Writes a Mol file for a molecule
    ARGUMENTS:

      - mol: the molecule
      - filename: the file to write to
      - includeStereo: (optional) toggles inclusion of stereochemical
        information in the output
      - confId: (optional) selects which conformation to output (-1 =
      default)
      - kekulize: (optional) triggers kekulization of the molecule before
      it's written,
        as suggested by the MDL spec.
      - forceV3000 (optional) force generation a V3000 mol block (happens
      automatically with
        more than 999 atoms or bonds)
  )DOC";
  m.def("MolToMolFile",
        (void (*)(const ROMol &, const std::string &, bool, int, bool,
                  bool))RDKit::MolToMolFile,
        "mol"_a, "filename"_a, "includeStereo"_a = true, "confId"_a = -1,
        "kekulize"_a = true, "forceV3000"_a = false, docString.c_str());

  docString =
      R"DOC(Writes a V3000 Mol file for a molecule
    ARGUMENTS:

      - mol: the molecule
      - filename: the file to write to
      - params: the MolWriterParams
      - confId: (optional) selects which conformation to output (-1 =
      default)
  )DOC";
  m.def("MolToV3KMolFile",
        (void (*)(const ROMol &, const std::string &, const MolWriterParams &,
                  int))RDKit::MolToV3KMolFile,
        "mol"_a, "filename"_a, "params"_a = true, "confId"_a = -1,
        docString.c_str());

  docString =
      R"DOC(Writes a V3000 Mol file for a molecule
    ARGUMENTS:

      - mol: the molecule
      - filename: the file to write to
      - includeStereo: (optional) toggles inclusion of stereochemical
        information in the output
      - confId: (optional) selects which conformation to output (-1 =
      default)
      - kekulize: (optional) triggers kekulization of the molecule before
      it's written,
        as suggested by the MDL spec.
  )DOC";
  m.def("MolToV3KMolFile",
        (void (*)(const ROMol &, const std::string &, bool, int,
                  bool))RDKit::MolToV3KMolFile,
        "mol"_a, "filename"_a, "includeStereo"_a = true, "confId"_a = -1,
        "kekulize"_a = true, docString.c_str());
  //

  docString =
      R"DOC(Returns a Marvin (Mrv) Mol block for a molecule
    ARGUMENTS:

      - mol: the molecule
      - includeStereo: (optional) toggles inclusion of stereochemical
        information in the output
      - confId: (optional) selects which conformation to output (-1 =
      default)
      - kekulize: (optional) triggers kekulization of the molecule before
      it's written.
      - prettyPrint: (optional) makes the output more human readable.

    RETURNS:

      a string
  )DOC";
  m.def("MolToMrvBlock",
        (std::string (*)(const ROMol &, bool, int, bool,
                         bool))RDKit::MolToMrvBlock,
        "mol"_a, "includeStereo"_a = true, "confId"_a = -1, "kekulize"_a = true,
        "prettyPrint"_a = false, docString.c_str());

  docString =
      R"DOC(Returns a Marvin (Mrv) Mol block for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - params: marvin write params
      - confId: (optional) selects which conformation to output (-1 =
      default)
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToMrvBlock",
        (std::string (*)(const ROMol &, const MrvWriterParams &,
                         int))RDKit::MolToMrvBlock,
        "mol"_a, "params"_a, "confId"_a = -1, docString.c_str());

  docString =
      R"DOC(Writes a Marvin (MRV) file for a molecule
     ARGUMENTS:
  
       - mol: the molecule
       - filename: the file to write to
       - includeStereo: (optional) toggles inclusion of stereochemical
         information in the output
       - confId: (optional) selects which conformation to output (-1 =
       default)
       - kekulize: (optional) triggers kekulization of the molecule before
       it's written.
       - prettyPrint: (optional) makes the output more human readable.
  )DOC";
  m.def("MolToMrvFile",
        (void (*)(const ROMol &, const std::string &, bool, int, bool,
                  bool))RDKit::MolToMrvFile,
        "mol"_a, "filename"_a, "includeStereo"_a = true, "confId"_a = -1,
        "kekulize"_a = true, "prettyPrint"_a = false, docString.c_str());

  docString =
      R"DOC(Writes a Marvin (MRV) file for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - filename: the file to write to
      - params: marvin write params
      - confId: (optional) selects which conformation to output (-1 =
      default)
  )DOC";
  m.def("MolToMrvFile",
        (void (*)(const ROMol &, const std::string &, const MrvWriterParams &,
                  int))RDKit::MolToMrvFile,
        "mol"_a, "filename"_a, "params"_a, "confId"_a = -1, docString.c_str());

  docString =
      R"DOC(Writes a CML block for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - confId: (optional) selects which conformation to output
      - kekulize: (optional) triggers kekulization of the molecule before
      it's written
  )DOC";
  m.def("MolToCMLBlock", RDKit::MolToCMLBlock, "mol"_a, "confId"_a = -1,
        "kekulize"_a = true, docString.c_str());

  docString =
      R"DOC(Writes a CML file for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - filename: the file to write to
      - confId: (optional) selects which conformation to output
      - kekulize: (optional) triggers kekulization of the molecule before
      it's written
  )DOC";
  m.def("MolToCMLFile", RDKit::MolToCMLFile, "mol"_a, "filename"_a,
        "confId"_a = -1, "kekulize"_a = true, docString.c_str());

  docString =
      R"DOC(Returns a XYZ block for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - confId: (optional) selects which conformation to output (-1 =
      default)
      - precision: precision of the coordinates
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToXYZBlock", RDKit::MolToXYZBlock, "mol"_a, "confId"_a = -1,
        "precision"_a = 6, docString.c_str());

  docString =
      R"DOC(Writes a XYZ file for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - filename: the file to write to
      - confId: (optional) selects which conformation to output (-1 =
      default)
      - precision: precision of the coordinates
  )DOC";
  m.def("MolToXYZFile", RDKit::MolToXYZFile, "mol"_a, "filename"_a,
        "confId"_a = -1, "precision"_a = 6, docString.c_str());

  nb::class_<RDKit::SmilesParserParams>(m, "SmilesParserParams",
                                        "Parameters controlling SMILES parsing")
      .def_rw("debugParse", &RDKit::SmilesParserParams::debugParse,
              "controls the amount of debugging information produced")
      .def_rw("parseName", &RDKit::SmilesParserParams::parseName,
              "controls whether or not the molecule name is also parsed")
      .def_rw("allowCXSMILES", &RDKit::SmilesParserParams::allowCXSMILES,
              "controls whether or not the CXSMILES extensions are parsed")
      .def_rw("strictCXSMILES", &RDKit::SmilesParserParams::strictCXSMILES,
              "controls whether or not problems in CXSMILES parsing "
              "causes molecule parsing to fail")
      .def_rw("sanitize", &RDKit::SmilesParserParams::sanitize,
              "controls whether or not the molecule is sanitized before "
              "being returned")
      .def_rw("removeHs", &RDKit::SmilesParserParams::removeHs,
              "controls whether or not Hs are removed before the "
              "molecule is returned")
      .def("__setattr__", &safeSetattr);
  nb::class_<RDKit::SmartsParserParams>(m, "SmartsParserParams",
                                        "Parameters controlling SMARTS parsing")
      .def_rw("debugParse", &RDKit::SmartsParserParams::debugParse,
              "controls the amount of debugging information produced")
      .def_rw("parseName", &RDKit::SmartsParserParams::parseName,
              "controls whether or not the molecule name is also parsed")
      .def_rw("allowCXSMILES", &RDKit::SmartsParserParams::allowCXSMILES,
              "controls whether or not the CXSMILES extensions are parsed")
      .def_rw("strictCXSMILES", &RDKit::SmartsParserParams::strictCXSMILES,
              "controls whether or not problems in CXSMILES parsing "
              "causes molecule parsing to fail")
      .def_rw("mergeHs", &RDKit::SmartsParserParams::mergeHs,
              "toggles merging H atoms in the SMARTS into neighboring atoms")
      .def("__setattr__", &safeSetattr);

  docString =
      R"DOC(Construct a molecule from a SMILES string.
       ARGUMENTS:
  
         - SMILES: the smiles string
  
         - params: used to provide optional parameters for the SMILES
         parsing
  
       RETURNS:
  
         a Mol object, None on failure.
  )DOC";
  m.def("MolFromSmiles", MolFromSmilesHelper, "SMILES"_a, "params"_a,
        docString.c_str(), nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a SMILES string.
    ARGUMENTS:
  
      - SMILES: the smiles string
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.
  
      - replacements: (optional) a dictionary of replacement strings (see
      below)
        Defaults to {}.
  
    RETURNS:
  
      a Mol object, None on failure.
  
     The optional replacements dict can be used to do string substitution
     of abbreviations in the input SMILES. The set of substitutions is
     repeatedly looped through until the string no longer changes. It is
     the responsibility of the caller to make sure that substitutions
     results in legal and sensible SMILES.
  
     Examples of replacements:
  
       CC{Q}C with {'{Q}':'OCCO'} -> CCOCCOC
       C{A}C{Q}C with {'{Q}':'OCCO', '{A}':'C1(CC1)'} -> CC1(CC1)COCCOC
       C{A}C{Q}C with {'{Q}':'{X}CC{X}', '{A}':'C1CC1', '{X}':'N'} ->
       CC1CC1CNCCNC
  )DOC";
  m.def("MolFromSmiles", RDKit::MolFromSmiles, "SMILES"_a, "sanitize"_a = true,
        "replacements"_a = nb::dict(), docString.c_str(),
        nb::rv_policy::take_ownership);

  docString = R"DOC(Construct an atom from a SMILES string)DOC";
  m.def("AtomFromSmiles", SmilesToAtom, "SMILES"_a, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString = R"DOC(Construct a bond from a SMILES string)DOC";
  m.def("BondFromSmiles", SmilesToBond, "SMILES"_a, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a SMARTS string.
    ARGUMENTS:
  
      - SMARTS: the smarts string
  
      - mergeHs: (optional) toggles the merging of explicit Hs in the query
      into the attached
        atoms.  So, for example, 'C[H]' becomes '[C;!H0]'.
        Defaults to 0.
  
      - replacements: (optional) a dictionary of replacement strings (see
      below)
        Defaults to {}. See the documentation for MolFromSmiles for an
        explanation.
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromSmarts", RDKit::MolFromSmarts, "SMARTS"_a, "mergeHs"_a = false,
        "replacements"_a = nb::dict(), docString.c_str(),
        nb::rv_policy::take_ownership);

  docString = R"DOC(Construct an atom from a SMARTS string)DOC";
  m.def("AtomFromSmarts", SmartsToAtom, "SMARTS"_a, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString = R"DOC(Construct a bond from a SMARTS string)DOC";
  m.def("BondFromSmarts", SmartsToBond, "SMARTS"_a, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a SMARTS string.
       ARGUMENTS:
  
         - SMARTS: the smarts string
  
         - params: used to provide optional parameters for the SMARTS
         parsing
  
       RETURNS:
  
         a Mol object, None on failure.
  )DOC";
  m.def("MolFromSmarts", MolFromSmartsHelper, "SMARTS"_a, "params"_a,
        docString.c_str(), nb::rv_policy::take_ownership);

  nb::class_<RDKit::SmilesWriteParams>(m, "SmilesWriteParams",
                                       "Parameters controlling SMILES writing")
      .def_rw("doIsomericSmiles", &RDKit::SmilesWriteParams::doIsomericSmiles,
              "include stereochemistry and isotope information")
      .def_rw("doKekule", &RDKit::SmilesWriteParams::doKekule,
              "kekulize the molecule before generating the SMILES and output "
              "single/double bonds. NOTE that the output is not canonical "
              "and that this will thrown an exception if the molecule "
              "cannot be kekulized")
      .def_rw("canonical", &RDKit::SmilesWriteParams::canonical,
              "generate canonical SMILES")
      .def_rw("cleanStereo", &RDKit::SmilesWriteParams::cleanStereo,
              "chiral centers are removed if they have duplicate "
              "sidechains")
      .def_rw("allBondsExplicit", &RDKit::SmilesWriteParams::allBondsExplicit,
              "include symbols for all bonds")
      .def_rw("allHsExplicit", &RDKit::SmilesWriteParams::allHsExplicit,
              "provide hydrogen counts for every atom")
      .def_rw("doRandom", &RDKit::SmilesWriteParams::doRandom,
              "randomize the output order. The resulting SMILES is not "
              "canonical")
      .def_rw("rootedAtAtom", &RDKit::SmilesWriteParams::rootedAtAtom,
              "make sure the SMILES starts at the specified atom. The "
              "resulting SMILES is not canonical")
      .def_rw("includeDativeBonds",
              &RDKit::SmilesWriteParams::includeDativeBonds,
              "include the "
              "RDKit extension for dative bonds. Otherwise dative bonds will "
              "be written as single bonds")
      .def_rw("ignoreAtomMapNumbers",
              &RDKit::SmilesWriteParams::ignoreAtomMapNumbers,
              "ignore atom map numbers when canonicalizing the molecule")
      .def("__setattr__", &safeSetattr);

  m.def("MolToSmiles",
        (std::string (*)(const ROMol &,
                         const SmilesWriteParams &))RDKit::MolToSmiles,
        "mol"_a, "params"_a,
        "Returns the canonical SMILES string for a molecule");

  docString =
      R"DOC(Returns the canonical SMILES string for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - isomericSmiles: (optional) include information about
      stereochemistry in
        the SMILES.  Defaults to true.
      - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds)
      in
        the SMILES.  Defaults to false.
      - rootedAtAtom: (optional) if non-negative, this forces the SMILES
        to start at a particular atom. Defaults to -1.  If not -1,
        overrides
        canonical setting.
      - canonical: (optional) if false no attempt will be made to
      canonicalize
        the molecule. Defaults to true.
      - allBondsExplicit: (optional) if true, all bond orders will be
      explicitly indicated
        in the output SMILES. Defaults to false.
      - allHsExplicit: (optional) if true, all H counts will be explicitly
      indicated
        in the output SMILES. Defaults to false.
      - doRandom: (optional) if true, randomize the traversal of the
      molecule graph,
        so we can generate random smiles. Defaults to false.  If true,
        overrides
        canonical setting.
      - ignoreAtomMapNumbers (optional) if true, ignores any atom map
      numbers when
        canonicalizing the molecule
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToSmiles",
        (std::string (*)(const ROMol &, bool, bool, int, bool, bool, bool, bool,
                         bool))RDKit::MolToSmiles,
        "mol"_a, "isomericSmiles"_a = true, "kekuleSmiles"_a = false,
        "rootedAtAtom"_a = -1, "canonical"_a = true,
        "allBondsExplicit"_a = false, "allHsExplicit"_a = false,
        "doRandom"_a = false, "ignoreAtomMapNumbers"_a = false,
        docString.c_str());

  docString =
      R"DOC(Returns the canonical SMILES string for a fragment of a
        molecule
    ARGUMENTS:
  
      - mol: the molecule
      - params: the SmilesWriteParams
      - atomsToUse : a list of atoms to include in the fragment
      - bondsToUse : (optional) a list of bonds to include in the fragment
        if not provided, all bonds between the atoms provided
        will be included.
      - atomSymbols : (optional) a list with the symbols to use for the
      atoms
        in the SMILES. This should have be mol.GetNumAtoms() long.
      - bondSymbols : (optional) a list with the symbols to use for the
      bonds
        in the SMILES. This should have be mol.GetNumBonds() long.
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolFragmentToSmiles", MolFragmentToSmilesHelper1<smilesfrag_gen>,
        "mol"_a, "params"_a, "atomsToUse"_a, "bondsToUse"_a = nb::none(),
        "atomSymbols"_a = nb::none(), "bondSymbols"_a = nb::none(),
        docString.c_str());

  docString =
      R"DOC(Returns the canonical SMILES string for a fragment of a
        molecule
    ARGUMENTS:
  
      - mol: the molecule
      - atomsToUse : a list of atoms to include in the fragment
      - bondsToUse : (optional) a list of bonds to include in the fragment
        if not provided, all bonds between the atoms provided
        will be included.
      - atomSymbols : (optional) a list with the symbols to use for the
      atoms
        in the SMILES. This should have be mol.GetNumAtoms() long.
      - bondSymbols : (optional) a list with the symbols to use for the
      bonds
        in the SMILES. This should have be mol.GetNumBonds() long.
      - isomericSmiles: (optional) include information about
      stereochemistry in
        the SMILES.  Defaults to true.
      - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds)
      in
        the SMILES.  Defaults to false.
      - rootedAtAtom: (optional) if non-negative, this forces the SMILES
        to start at a particular atom. Defaults to -1.  If not -1,
        over-rides
        setting for canonical.
      - canonical: (optional) if false no attempt will be made to
      canonicalize
        the molecule. Defaults to true.
      - allBondsExplicit: (optional) if true, all bond orders will be
      explicitly indicated
        in the output SMILES. Defaults to false.
      - allHsExplicit: (optional) if true, all H counts will be explicitly
      indicated
        in the output SMILES. Defaults to false.
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolFragmentToSmiles", MolFragmentToSmilesHelper2<smilesfrag_gen>,
        "mol"_a, "atomsToUse"_a, "bondsToUse"_a = nb::none(),
        "atomSymbols"_a = nb::none(), "bondSymbols"_a = nb::none(),
        "isomericSmiles"_a = true, "kekuleSmiles"_a = false,
        "rootedAtAtom"_a = -1, "canonical"_a = true,
        "allBondsExplicit"_a = false, "allHsExplicit"_a = false,
        docString.c_str());

  nb::enum_<RDKit::SmilesWrite::CXSmilesFields>(m, "CXSmilesFields")
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

  nb::enum_<RDKit::v2::FileParsers::SCSRBaseHbondOptions>(
      m, "SCSRBaseHbondOptions")
      .value("Ignore", RDKit::v2::FileParsers::SCSRBaseHbondOptions::Ignore)
      .value("UseSapAll",
             RDKit::v2::FileParsers::SCSRBaseHbondOptions::UseSapAll)
      .value("UseSapOne",
             RDKit::v2::FileParsers::SCSRBaseHbondOptions::UseSapOne)
      .value("Auto", RDKit::v2::FileParsers::SCSRBaseHbondOptions::Auto);

  nb::enum_<RDKit::v2::FileParsers::SCSRTemplateNames>(m, "SCSRTemplateNames")
      .value("UseFirstName",
             RDKit::v2::FileParsers::SCSRTemplateNames::UseFirstName)
      .value("UseSecondName",
             RDKit::v2::FileParsers::SCSRTemplateNames::UseSecondName)
      .value("AsEntered", RDKit::v2::FileParsers::SCSRTemplateNames::AsEntered);

  nb::enum_<RDKit::RestoreBondDirOption>(m, "RestoreBondDirOption")
      .value("RestoreBondDirOptionClear",
             RDKit::RestoreBondDirOption::RestoreBondDirOptionClear)
      .value("RestoreBondDirOptionTrue",
             RDKit::RestoreBondDirOption::RestoreBondDirOptionTrue);

  m.def(
      "MolToCXSmiles",
      (std::string (*)(const ROMol &, const SmilesWriteParams &, std::uint32_t,
                       RestoreBondDirOption))RDKit::MolToCXSmiles,
      "mol"_a, "params"_a,
      "flags"_a = RDKit::SmilesWrite::CXSmilesFields::CX_ALL,
      "restoreBondDirs"_a =
          RDKit::RestoreBondDirOption::RestoreBondDirOptionClear,
      "Returns the CXSMILES string for a molecule");

  docString =
      R"DOC(Returns the CXSMILES string for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - isomericSmiles: (optional) include information about
      stereochemistry in
        the SMILES.  Defaults to true.
      - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds)
      in
        the SMILES.  Defaults to false.
      - rootedAtAtom: (optional) if non-negative, this forces the SMILES
        to start at a particular atom. Defaults to -1.
      - canonical: (optional) if false no attempt will be made to
      canonicalize
        the molecule. Defaults to true.
      - allBondsExplicit: (optional) if true, all bond orders will be
      explicitly indicated
        in the output SMILES. Defaults to false.
      - allHsExplicit: (optional) if true, all H counts will be explicitly
      indicated
        in the output SMILES. Defaults to false.
      - doRandom: (optional) if true, randomizes the traversal of the
      molecule graph,
        so we can generate random smiles. Defaults to false.  If true,
        overrides
        canonical setting.
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToCXSmiles",
        (std::string (*)(const ROMol &, bool, bool, int, bool, bool, bool,
                         bool))RDKit::MolToCXSmiles,
        "mol"_a, "isomericSmiles"_a = true, "kekuleSmiles"_a = false,
        "rootedAtAtom"_a = -1, "canonical"_a = true,
        "allBondsExplicit"_a = false, "allHsExplicit"_a = false,
        "doRandom"_a = false, docString.c_str());

  docString =
      R"DOC(Returns the CXSMILES string for a fragment of a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - params: the SmilesWriteParams
      - atomsToUse : a list of atoms to include in the fragment
      - bondsToUse : (optional) a list of bonds to include in the fragment
        if not provided, all bonds between the atoms provided
        will be included.
      - atomSymbols : (optional) a list with the symbols to use for the
      atoms
        in the SMILES. This should have be mol.GetNumAtoms() long.
      - bondSymbols : (optional) a list with the symbols to use for the
      bonds
        in the SMILES. This should have be mol.GetNumBonds() long.
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolFragmentToCXSmiles", MolFragmentToSmilesHelper1<cxsmilesfrag_gen>,
        "mol"_a, "params"_a, "atomsToUse"_a, "bondsToUse"_a = nb::none(),
        "atomSymbols"_a = nb::none(), "bondSymbols"_a = nb::none(),
        docString.c_str());

  docString =
      R"DOC(Returns the CXSMILES string for a fragment of a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - atomsToUse : a list of atoms to include in the fragment
      - bondsToUse : (optional) a list of bonds to include in the fragment
        if not provided, all bonds between the atoms provided
        will be included.
      - atomSymbols : (optional) a list with the symbols to use for the
      atoms
        in the SMILES. This should have be mol.GetNumAtoms() long.
      - bondSymbols : (optional) a list with the symbols to use for the
      bonds
        in the SMILES. This should have be mol.GetNumBonds() long.
      - isomericSmiles: (optional) include information about
      stereochemistry in
        the SMILES.  Defaults to true.
      - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds)
      in
        the SMILES.  Defaults to false.
      - rootedAtAtom: (optional) if non-negative, this forces the SMILES
        to start at a particular atom. Defaults to -1.  If not -1,
        overrides
        canonical setting.
      - canonical: (optional) if false no attempt will be made to
      canonicalize
        the molecule. Defaults to true.
      - allBondsExplicit: (optional) if true, all bond orders will be
      explicitly indicated
        in the output SMILES. Defaults to false.
      - allHsExplicit: (optional) if true, all H counts will be explicitly
      indicated
        in the output SMILES. Defaults to false.
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolFragmentToCXSmiles", MolFragmentToSmilesHelper2<cxsmilesfrag_gen>,
        "mol"_a, "atomsToUse"_a, "bondsToUse"_a = nb::none(),
        "atomSymbols"_a = nb::none(), "bondSymbols"_a = nb::none(),
        "isomericSmiles"_a = true, "kekuleSmiles"_a = false,
        "rootedAtAtom"_a = -1, "canonical"_a = true,
        "allBondsExplicit"_a = false, "allHsExplicit"_a = false,
        docString.c_str());

  docString =
      R"DOC(Returns a SMARTS string for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - isomericSmiles: (optional) include information about
      stereochemistry in
        the SMARTS.  Defaults to true.
      - rootedAtomAtom: (optional) the atom index to start the SMARTS
      from.
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToSmarts",
        (std::string (*)(const ROMol &, bool, int))RDKit::MolToSmarts, "mol"_a,
        "isomericSmiles"_a = true, "rootedAtAtom"_a = -1, docString.c_str());

  docString =
      R"DOC(Returns a SMARTS string for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - params: SmilesWriteParams controlling the SMARTS generation
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToSmarts",
        (std::string (*)(const ROMol &,
                         const SmilesWriteParams &))RDKit::MolToSmarts,
        "mol"_a, "params"_a, docString.c_str());

  docString =
      R"DOC(Returns a SMARTS string for a fragment of a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - atomsToUse: indices of atoms to include in the SMARTS string
      - bondsToUse: indices of bonds to include in the SMARTS string
      (optional)
      - isomericSmarts: (optional) include information about
      stereochemistry in
        the SMARTS.  Defaults to true.
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolFragmentToSmarts", molFragmentToSmarts, "mol"_a, "atomsToUse"_a,
        "bondsToUse"_a = nb::none(), "isomericSmarts"_a = true,
        docString.c_str());

  docString =
      R"DOC(Returns a SMARTS string for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - isomericSmiles: (optional) include information about
      stereochemistry in
        the SMARTS.  Defaults to true.
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToCXSmarts",
        (std::string (*)(const ROMol &, bool))RDKit::MolToCXSmarts, "mol"_a,
        "isomericSmiles"_a = true, docString.c_str());

  docString =
      R"DOC(Returns a SMARTS string for a fragment of a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - atomsToUse: indices of atoms to include in the SMARTS string
      - bondsToUse: indices of bonds to include in the SMARTS string
      (optional)
      - isomericSmarts: (optional) include information about
      stereochemistry in
        the SMARTS.  Defaults to true.
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolFragmentToCXSmarts", molFragmentToCXSmarts, "mol"_a, "atomsToUse"_a,
        "bondsToUse"_a = nb::none(), "isomericSmarts"_a = true,
        docString.c_str());

  docString =
      R"DOC(Writes a molecule to a TPL file.
    ARGUMENTS:
  
      - mol: the molecule
      - fileName: name of the file to write
      - partialChargeProp: name of the property to use for partial charges
        Defaults to '_GasteigerCharge'.
      - writeFirstConfTwice: Defaults to False.
        This should be set to True when writing TPLs to be read by
        the CombiCode.
  )DOC";
  m.def("MolToTPLFile", RDKit::MolToTPLFile, "mol"_a, "fileName"_a,
        "partialChargeProp"_a = "_GasteigerCharge",
        "writeFirstConfTwice"_a = false, docString.c_str());

  docString =
      R"DOC(Returns the Tpl block for a molecule.
    ARGUMENTS:
  
      - mol: the molecule
      - partialChargeProp: name of the property to use for partial charges
        Defaults to '_GasteigerCharge'.
      - writeFirstConfTwice: Defaults to False.
        This should be set to True when writing TPLs to be read by
        the CombiCode.
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToTPLBlock", RDKit::MolToTPLText, "mol"_a,
        "partialChargeProp"_a = "_GasteigerCharge",
        "writeFirstConfTwice"_a = false, docString.c_str());

  docString =
      R"DOC(Construct a molecule from a PDB file.
    ARGUMENTS:
  
      - pdbFileName: name of the file to read
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to true.
  
      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.
  
      - flavor: (optional)
  
      - proximityBonding: (optional) toggles automatic proximity bonding
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromPDBFile", RDKit::MolFromPDBFile, "pdbFileName"_a,
        "sanitize"_a = true, "removeHs"_a = true, "flavor"_a = 0,
        "proximityBonding"_a = true, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from a PDB block.
    ARGUMENTS:
  
      - molBlock: string containing the PDB block
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.
  
      - removeHs: (optional) toggles removing hydrogens from the molecule.
        This only make sense when sanitization is done.
        Defaults to true.
  
      - flavor: (optional)
  
      - proximityBonding: (optional) toggles automatic proximity bonding
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromPDBBlock", RDKit::MolFromPDBBlock, "molBlock"_a,
        "sanitize"_a = true, "removeHs"_a = true, "flavor"_a = 0,
        "proximityBonding"_a = true, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Returns a PDB block for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - confId: (optional) selects which conformation to output (-1 =
      default)
      - flavor: (optional)
              - flavor & 1 : Write MODEL/ENDMDL lines around each record
              - flavor & 2 : Don't write any CONECT records
              - flavor & 4 : Write CONECT records in both directions
              - flavor & 8 : Don't use multiple CONECTs to encode bond
              order
  
              - flavor & 16 : Write MASTER record
              - flavor & 32 : Write TER record
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToPDBBlock", RDKit::MolToPDBBlock, "mol"_a, "confId"_a = -1,
        "flavor"_a = 0, docString.c_str());

  docString =
      R"DOC(Writes a PDB file for a molecule
    ARGUMENTS:
  
      - mol: the molecule
      - filename: name of the file to write
      - confId: (optional) selects which conformation to output (-1 =
      default)
      - flavor: (optional)
              - flavor & 1 : Write MODEL/ENDMDL lines around each record
              - flavor & 2 : Don't write any CONECT records
              - flavor & 4 : Write CONECT records in both directions
              - flavor & 8 : Don't use multiple CONECTs to encode bond
              order
  
              - flavor & 16 : Write MASTER record
              - flavor & 32 : Write TER record
  )DOC";
  m.def("MolToPDBFile", RDKit::MolToPDBFile, "mol"_a, "filename"_a,
        "confId"_a = -1, "flavor"_a = 0, docString.c_str());

  docString =
      R"DOC(Construct a molecule from a sequence string (currently
        supports standard amino acids, DNA and RNA bases).
    ARGUMENTS:
  
      - text: string containing the sequence
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.
  
      - flavor: (optional)
          - 0 Protein, L amino acids (default)
          - 1 Protein, D amino acids
          - 2 RNA, no cap
          - 3 RNA, 5' cap
          - 4 RNA, 3' cap
          - 5 RNA, both caps
          - 6 DNA, no cap
          - 7 DNA, 5' cap
          - 8 DNA, 3' cap
          - 9 DNA, both caps
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromSequence", RDKit::MolFromSequence, "text"_a,
        "sanitize"_a = true, "flavor"_a = 0, docString.c_str(),
        nb::rv_policy::take_ownership);

  docString =
      R"DOC(Returns the sequence string for a molecule
    ARGUMENTS:
  
      - mol: the molecule
  
    NOTE: the molecule should contain monomer information in
    AtomMonomerInfo structures
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToSequence", RDKit::MolToSequence, "mol"_a, docString.c_str());

  docString =
      R"DOC(Construct a molecule from a FASTA string (currently supports
        standard amino acids, DNA and RNA bases).
    ARGUMENTS:
  
      - text: string containing the FASTA
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to True.
  
  - flavor: (optional)
      - 0 Protein, L amino acids (default)
      - 1 Protein, D amino acids
      - 2 RNA, no cap
      - 3 RNA, 5' cap
      - 4 RNA, 3' cap
      - 5 RNA, both caps
      - 6 DNA, no cap
      - 7 DNA, 5' cap
      - 8 DNA, 3' cap
      - 9 DNA, both caps
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromFASTA", RDKit::MolFromFASTA, "text"_a, "sanitize"_a = true,
        "flavor"_a = 0, docString.c_str(), nb::rv_policy::take_ownership);

  docString =
      R"DOC(Returns the FASTA string for a molecule
    ARGUMENTS:
  
      - mol: the molecule
  
    NOTE: the molecule should contain monomer information in
    AtomMonomerInfo structures
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToFASTA", RDKit::MolToFASTA, "mol"_a, docString.c_str());

  docString =
      R"DOC(Construct a molecule from a HELM string (currently supports
        standard
        amino acids, DNA and RNA bases).
    ARGUMENTS:
  
      - text: string containing the HELM
  
      - sanitize: (optional) toggles sanitization of the molecule.
        Defaults to true.
  
    RETURNS:
  
      a Mol object, None on failure.
  )DOC";
  m.def("MolFromHELM", RDKit::MolFromHELM, "text"_a, "sanitize"_a = true,
        docString.c_str(), nb::rv_policy::take_ownership);

  docString =
      R"DOC(Returns the HELM string for a molecule
    ARGUMENTS:
  
      - mol: the molecule
  
    NOTE: the molecule should contain monomer information in
    AtomMonomerInfo structures
  
    RETURNS:
  
      a string
  )DOC";
  m.def("MolToHELM", RDKit::MolToHELM, "mol"_a, docString.c_str());

  docString =
      R"DOC(Returns the canonical atom ranking for each atom of a
        molecule fragment.
    If breakTies is False, this returns the symmetry class for each atom.
    The symmetry class is used by the canonicalization routines to type
    each atom based on the whole chemistry of the molecular graph.  Any
    atom with the same rank (symmetry class) is indistinguishable.  For
    example:
  
      >>> mol = MolFromSmiles('C1NCN1')
      >>> list(CanonicalRankAtoms(mol, breakTies=False))
      [0,1,0,1]
  
    In this case the carbons have the same symmetry class and the nitrogens
    have the same
    symmetry class.  From the perspective of the Molecular Graph, they are
    identical.
  
    ARGUMENTS:
  
      - mol: the molecule
      - breakTies: (optional) force breaking of ranked ties [default=True]
      - includeChirality: (optional) use chiral information when computing
      rank [default=True]
      - includeIsotopes: (optional) use isotope information when computing
      rank [default=True]
      - includeAtomMaps: (optional) use atom map information when computing
      rank [default=True]
      - includeChiralPresence: (optional) use information about whether or
      not chirality is specified when computing rank [default=False]
  
    RETURNS:
  
      a string
  )DOC";
  m.def("CanonicalRankAtoms", CanonicalRankAtoms, "mol"_a, "breakTies"_a = true,
        "includeChirality"_a = true, "includeIsotopes"_a = true,
        "includeAtomMaps"_a = true, "includeChiralPresence"_a = false,
        docString.c_str());

  docString =
      R"DOC(Returns the canonical atom ranking for each atom of a
        molecule fragment
    See help(CanonicalRankAtoms) for more information.
  
     >>> mol = MolFromSmiles('C1NCN1.C1NCN1')
     >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(0,4),
     breakTies=False))
     [4,6,4,6,-1,-1,-1,-1]
     >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(4,8),
     breakTies=False))
     [-1,-1,-1,-1,4,6,4,6]
  
    ARGUMENTS:
  
      - mol: the molecule
      - atomsToUse : a list of atoms to include in the fragment
      - bondsToUse : (optional) a list of bonds to include in the fragment
        if not provided, no bonds will be used
      - atomSymbols : (optional) a list with the symbols to use for the
      atoms
        in the SMILES. This should have be mol.GetNumAtoms() long.
      - breakTies: (optional) force breaking of ranked ties
      - includeChirality: (optional) use chiral information when computing
      rank [default=True]
      - includeIsotopes: (optional) use isotope information when computing
      rank [default=True]
      - includeAtomMaps: (optional) use atom map information when computing
      rank [default=True]
      - includeChiralPresence: (optional) use information about whether or
      not chirality is specified when computing rank [default=False]
  
    RETURNS:
  
      a string
  )DOC";
  m.def("CanonicalRankAtomsInFragment", CanonicalRankAtomsInFragment, "mol"_a,
        "atomsToUse"_a, "bondsToUse"_a = 0, "atomSymbols"_a = 0,
        "breakTies"_a = true, "includeChirality"_a = true,
        "includeIsotopes"_a = true, "includeAtomMaps"_a = true,
        "includeChiralPresence"_a = false, docString.c_str());

  m.def("CanonicalizeEnhancedStereo", CanonicalizeEnhancedStereo, "mol"_a);

  m.def(
      "CreateAtomIntPropertyList", FileParserUtils::createAtomIntPropertyList,
      "mol"_a, "propName"_a, "missingValueMarker"_a = "", "lineSize"_a = 190,
      "creates a list property on the molecule from individual atom property values");
  m.def(
      "CreateAtomDoublePropertyList",
      FileParserUtils::createAtomDoublePropertyList, "mol"_a, "propName"_a,
      "missingValueMarker"_a = "", "lineSize"_a = 190,
      "creates a list property on the molecule from individual atom property values");
  m.def(
      "CreateAtomBoolPropertyList", FileParserUtils::createAtomBoolPropertyList,
      "mol"_a, "propName"_a, "missingValueMarker"_a = "", "lineSize"_a = 190,
      "creates a list property on the molecule from individual atom property values");
  m.def(
      "CreateAtomStringPropertyList",
      FileParserUtils::createAtomStringPropertyList, "mol"_a, "propName"_a,
      "missingValueMarker"_a = "", "lineSize"_a = 190,
      "creates a list property on the molecule from individual atom property values");

  m.def(
      "CreateBondIntPropertyList", FileParserUtils::createBondIntPropertyList,
      "mol"_a, "propName"_a, "missingValueMarker"_a = "", "lineSize"_a = 190,
      "creates a list property on the molecule from individual bond property values");
  m.def(
      "CreateBondDoublePropertyList",
      FileParserUtils::createBondDoublePropertyList, "mol"_a, "propName"_a,
      "missingValueMarker"_a = "", "lineSize"_a = 190,
      "creates a list property on the molecule from individual bond property values");
  m.def(
      "CreateBondBoolPropertyList", FileParserUtils::createBondBoolPropertyList,
      "mol"_a, "propName"_a, "missingValueMarker"_a = "", "lineSize"_a = 190,
      "creates a list property on the molecule from individual bond property values");
  m.def(
      "CreateBondStringPropertyList",
      FileParserUtils::createBondStringPropertyList, "mol"_a, "propName"_a,
      "missingValueMarker"_a = "", "lineSize"_a = 190,
      "creates a list property on the molecule from individual bond property values");

  m.def("MolToRandomSmilesVect", RDKit::MolToRandomSmilesHelper, "mol"_a,
        "numSmiles"_a, "randomSeed"_a = 0, "isomericSmiles"_a = true,
        "kekuleSmiles"_a = false, "allBondsExplicit"_a = false,
        "allHsExplicit"_a = false,
        "returns a list of SMILES generated using the randomSmiles algorithm");

#ifdef RDK_USE_BOOST_IOSTREAMS
  nb::class_<RDKit::PNGMetadataParams>(
      m, "PNGMetadataParams",
      "Parameters controlling metadata included in PNG images")
      .def_rw("includePkl", &RDKit::PNGMetadataParams::includePkl,
              "toggles inclusion of molecule pickle (default=True)")
      .def_rw("includeSmiles", &RDKit::PNGMetadataParams::includeSmiles,
              "toggles inclusion of molecule CXSMILES (default=True)")
      .def_rw("includeMol", &RDKit::PNGMetadataParams::includeMol,
              "toggles inclusion of molecule molblock (default=False)")
      .def_rw(
          "propertyFlags", &RDKit::PNGMetadataParams::propertyFlags,
          "choose properties to be included in the pickle (default=rdkit.Chem.rdchem.PropertyPickleOptions.NoProps)")
      .def_rw(
          "smilesWriteParams", &RDKit::PNGMetadataParams::smilesWriteParams,
          "choose SmilesWriteParams for the CXSMILES string (default=rdkit.Chem.rdmolfiles.SmilesWriteParams())")
      .def_rw(
          "cxSmilesFlags", &RDKit::PNGMetadataParams::cxSmilesFlags,
          "choose CXSMILES fields to be included in the CXSMILES string (default=rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL)")
      .def_rw(
          "restoreBondDirs", &RDKit::PNGMetadataParams::restoreBondDirs,
          "choose what to do with bond dirs in the CXSMILES string (default=rdkit.Chem.rdmolfiles.RestoreBondDirOption.RestoreBondDirOptionClear)")
      .def("__setattr__", &safeSetattr);

  docString =
      R"DOC(Construct a molecule from metadata in a PNG string.

       ARGUMENTS:

         - png: the PNG string

         - params: used to provide optional parameters for the metadata
         parsing

       RETURNS:
         a Mol object, None on failure.
    )DOC";
  m.def("MolFromPNGString", MolFromPNGString, "png"_a, "params"_a = nb::none(),
        docString.c_str(), nb::rv_policy::take_ownership);

  docString =
      R"DOC(Construct a molecule from metadata in a PNG file.

       ARGUMENTS:

         - filename: the PNG filename

         - params: used to provide optional parameters for the metadata
         parsing

       RETURNS:
         a Mol object, None on failure.)DOC";

  std::string cdxml_notes = R"DOC()DOC";

  m.def("MolFromPNGFile", MolFromPNGFile, "filename"_a, "params"_a = nb::none(),
        docString.c_str(), nb::rv_policy::take_ownership);

  m.def("MolsFromPNGString", MolsFromPNGString, "png"_a,
        "tag"_a = PNGData::pklTag, "params"_a = nb::none(),
        "returns a tuple of molecules constructed from the PNG string");
  m.def("MolsFromPNGFile", MolsFromPNGFile, "filename"_a,
        "tag"_a = PNGData::pklTag, "params"_a = nb::none(),
        "returns a tuple of molecules constructed from the PNG file");
#endif

  docString =
      R"DOC(Construct a molecule from a cdxml file.

       Note that the CDXML format is large and complex, the RDKit doesn't
       support full functionality, just the base ones required for molecule
       and reaction parsing.

       ARGUMENTS:

         - filename: the cdxml filename

         - sanitize: if True, sanitize the molecules [default True]

         - removeHs: if True, convert explicit Hs into implicit Hs.
         [default True]

       RETURNS:
         an iterator of parsed Mol objects.)DOC";

  m.def("MolsFromCDXMLFile", MolsFromCDXMLFile, "filename"_a,
        "sanitize"_a = true, "removeHs"_a = true, docString.c_str());

  docString =
      R"DOC(Construct a molecule from a cdxml string.

       Note that the CDXML format is large and complex, the RDKit doesn't
       support full functionality, just the base ones required for molecule
       and reaction parsing.

       ARGUMENTS:

         - cdxml: the cdxml string

         - sanitize: if True, sanitize the molecules [default True]

         - removeHs: if True, convert explicit Hs into implicit Hs.
         [default True]

       RETURNS:
         an iterator of parsed Mol objects.)DOC";

  m.def("MolsFromCDXML", MolsFromCDXML, "cdxml"_a, "sanitize"_a = true,
        "removeHs"_a = true, docString.c_str());

  nb::enum_<RDKit::v2::CDXMLParser::CDXMLFormat>(m, "CDXMLFormat")
      .value("CDXML", RDKit::v2::CDXMLParser::CDXMLFormat::CDXML)
      .value("CDX", RDKit::v2::CDXMLParser::CDXMLFormat::CDX)
      .value("Auto", RDKit::v2::CDXMLParser::CDXMLFormat::Auto);

  nb::class_<RDKit::v2::CDXMLParser::CDXMLParserParams>(
      m, "CDXMLParserParams",
      "Parameters controlling conversion of a CDXML document to molecules")
      .def(nb::init<>(), "Construct a default CDXMLFormat")
      .def(nb::init<bool, bool, RDKit::v2::CDXMLParser::CDXMLFormat>(),
           "sanitize"_a, "removeHs"_a, "format"_a)
      .def_rw("sanitize", &RDKit::v2::CDXMLParser::CDXMLParserParams::sanitize,
              "controls whether or not the molecule is sanitized before "
              "being returned")
      .def_rw("removeHs", &RDKit::v2::CDXMLParser::CDXMLParserParams::removeHs,
              "controls whether or not Hs are removed before the "
              "molecule is returned")
      .def_rw(
          "format", &RDKit::v2::CDXMLParser::CDXMLParserParams::format,
          "ChemDraw format One of Auto, CDXML, CDX.  For data streams, Auto defaults to CDXML")
      .def("__setattr__", &safeSetattr);

  docString =
      R"DOC(Construct a molecule from a cdxml file.

       Note: that the CDXML format is large and complex, the RDKit doesn't
       support full functionality, just the base ones required for molecule
       and reaction parsing.

       Note: If the ChemDraw extensions are available,
          CDXMLFormat::Auto attempts to see if the input string is CDXML or
          CDX,
       If not, it defaults to CDXML

       ARGUMENTS:

         - filename: the cdxml filename

         - pyParams: CDXParserParams, see CDXParserParams for usage

       RETURNS:
         a tuple  of parsed Mol objects.)DOC";

  m.def("MolsFromCDXMLFile", MolsFromCDXMLFileHelper, "filename"_a, "params"_a,
        docString.c_str());

  docString =
      R"DOC(Construct a molecule from a cdxml string.

       Note that the CDXML format is large and complex, the RDKit doesn't
       support full functionality, just the base ones required for molecule
       and reaction parsing.

       Note: in this function CDXMLFormat::Auto currently defaults to CDXML

       ARGUMENTS:

         - cdxml: the cdxml string

         - pyParams: CDXParserParams, see CDXParserParams for usage

       RETURNS:
         a tuple of parsed Mol objects.)DOC";

  m.def("MolsFromCDXML", MolsFromCDXMLHelper, "cdxml"_a, "params"_a,
        docString.c_str());

  docString =
      R"DOC(Returns true if the RDKit is built with ChemDraw CDX
        support)DOC";
  m.def("HasChemDrawCDXSupport", RDKit::v2::CDXMLParser::hasChemDrawCDXSupport,
        docString.c_str());

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
  m.def("MolMetadataToPNGFile", addMolToPNGFileHelper, "mol"_a, "filename"_a,
        "includePkl"_a = true, "includeSmiles"_a = true, "includeMol"_a = false,
        docString.c_str());

  docString =
      R"DOC(Adds molecular metadata to PNG data read from a file.

     ARGUMENTS:

       - mol: the molecule

       - filename: the PNG filename

       - params: an instance of PNGMetadataParams

     RETURNS:
       the updated PNG data)DOC";
  m.def("MolMetadataToPNGFile", addMolToPNGFileHelperParams, "mol"_a,
        "filename"_a, "params"_a, docString.c_str());

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
  m.def("MolMetadataToPNGString", addMolToPNGStringHelper, "mol"_a, "png"_a,
        "includePkl"_a = true, "includeSmiles"_a = true, "includeMol"_a = false,
        docString.c_str());

  docString =
      R"DOC(Adds molecular metadata to a PNG string.

     ARGUMENTS:

       - mol: the molecule

       - png: the PNG string

       - params: an instance of PNGMetadataParams

     RETURNS:
       the updated PNG data)DOC";
  m.def("MolMetadataToPNGString", addMolToPNGStringHelperParams, "mol"_a,
        "png"_a, "params"_a, docString.c_str());

  docString =
      R"DOC(Adds metadata to PNG data read from a file.

     ARGUMENTS:

       - metadata: dict with the metadata to be written
                   (keys and values should be strings)

       - filename: the PNG filename

     RETURNS:
       the updated PNG data)DOC";
  m.def("AddMetadataToPNGFile", addMetadataToPNGFileHelper, "metadata"_a,
        "filename"_a, docString.c_str());

  docString =
      R"DOC(Adds metadata to a PNG string.

     ARGUMENTS:

       - metadata: dict with the metadata to be written
                   (keys and values should be strings)

       - png: the PNG string

     RETURNS:
       the updated PNG data)DOC";
  m.def("AddMetadataToPNGString", addMetadataToPNGStringHelper, "metadata"_a,
        "png"_a, docString.c_str());

  m.def("MetadataFromPNGFile", MetadataFromPNGFile, "filename"_a,
        "asList"_a = false,
        "Returns a dict with all metadata from the PNG file. Keys are "
        "strings, values are bytes. "
        "If asList is True, a list of (key, value) tuples is returned; "
        "this enables retrieving multiple values sharing the same key.");

  m.def("MetadataFromPNGString", MetadataFromPNGString, "png"_a,
        "asList"_a = false,
        "Returns a dict with all metadata from the PNG string. Keys are "
        "strings, values are bytes. "
        "If asList is True, a list of (key, value) tuples is returned; "
        "this enables retrieving multiple values sharing the same key.");
#endif
/********************************************************
 * MolSupplier stuff
 *******************************************************/
#ifdef SUPPORT_COMPRESSED_SUPPLIERS
  wrap_compressedsdsupplier(m);
#endif
  wrap_sdsupplier(m);
  wrap_forwardsdsupplier(m);
  wrap_tdtsupplier(m);
  wrap_smisupplier(m);
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
  wrap_maesupplier(m);
#endif

  //   /********************************************************
  //    * MolWriter stuff
  //    *******************************************************/
  wrap_smiwriter(m);
  wrap_sdwriter(m);
  wrap_tdtwriter(m);
  wrap_pdbwriter(m);
#ifdef RDK_BUILD_MAEPARSER_SUPPORT
  wrap_maewriter(m);
#endif

#ifdef RDK_BUILD_THREADSAFE_SSS
  /********************************************************
   * MultithreadedMolWriter stuff
   *******************************************************/
  wrap_multiSmiSupplier(m);
  wrap_multiSDSupplier(m);
#endif
}
