//
//  Copyright (C) 2003-2019 Greg Landrum and Rational Discovery LLC
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
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/FileParsers/SequenceParsers.h>
#include <GraphMol/FileParsers/SequenceWriters.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>

#include <RDBoost/Wrap.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SanitException.h>

namespace python = boost::python;
using namespace RDKit;

void rdSanitExceptionTranslator(RDKit::MolSanitizeException const &x) {
  std::ostringstream ss;
  ss << "Sanitization error: " << x.message();
  PyErr_SetString(PyExc_ValueError, ss.str().c_str());
}
void rdBadFileExceptionTranslator(RDKit::BadFileException const &x) {
  std::ostringstream ss;
  ss << "File error: " << x.message();
  PyErr_SetString(PyExc_IOError, ss.str().c_str());
}

namespace RDKit {
std::string pyObjectToString(python::object input) {
  python::extract<std::string> ex(input);
  if (ex.check()) return ex();
  std::wstring ws = python::extract<std::wstring>(input);
  return std::string(ws.begin(), ws.end());
}

ROMol *MolFromSmiles(python::object ismiles, bool sanitize,
                     python::dict replDict) {
  std::map<std::string, std::string> replacements;
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(replDict.keys().attr("__len__")());
       ++i) {
    replacements[python::extract<std::string>(replDict.keys()[i])] =
        python::extract<std::string>(replDict.values()[i]);
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
  for (unsigned int i = 0;
       i < python::extract<unsigned int>(replDict.keys().attr("__len__")());
       ++i) {
    replacements[python::extract<std::string>(replDict.keys()[i])] =
        python::extract<std::string>(replDict.values()[i]);
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
ROMol *MolFromTPLFile(const char *filename, bool sanitize = true,
                      bool skipFirstConf = false) {
  RWMol *newM;
  try {
    newM = TPLFileToMol(filename, sanitize, skipFirstConf);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.message());
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

ROMol *MolFromMolFile(const char *molFilename, bool sanitize, bool removeHs,
                      bool strictParsing) {
  RWMol *newM = nullptr;
  try {
    newM = MolFileToMol(molFilename, sanitize, removeHs, strictParsing);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.message());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.message() << std::endl;
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
    BOOST_LOG(rdWarningLog) << e.message() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromSVG(python::object imolBlock, bool sanitize, bool removeHs) {
  RWMol *res = nullptr;
  res = RDKitSVGToMol(pyObjectToString(imolBlock), sanitize, removeHs);
  return static_cast<ROMol *>(res);
}

ROMol *MolFromMol2File(const char *molFilename, bool sanitize = true,
                       bool removeHs = true, bool cleanupSubstructures = true) {
  RWMol *newM;
  try {
    newM = Mol2FileToMol(molFilename, sanitize, removeHs, CORINA,
                         cleanupSubstructures);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.message());
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
    newM = Mol2DataStreamToMol(inStream, sanitize, removeHs, CORINA,
                               cleanupSubstructures);
  } catch (...) {
    newM = nullptr;
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromPDBFile(const char *filename, bool sanitize, bool removeHs,
                      unsigned int flavor, bool proximityBonding) {
  RWMol *newM = nullptr;
  try {
    newM = PDBFileToMol(filename, sanitize, removeHs, flavor, proximityBonding);
  } catch (RDKit::BadFileException &e) {
    PyErr_SetString(PyExc_IOError, e.message());
    throw python::error_already_set();
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.message() << std::endl;
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
    BOOST_LOG(rdWarningLog) << e.message() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

ROMol *MolFromSequence(python::object seq, bool sanitize, int flavor) {
  RWMol *newM = nullptr;
  try {
    newM = SequenceToMol(pyObjectToString(seq), sanitize, flavor);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.message() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}
ROMol *MolFromFASTA(python::object seq, bool sanitize, int flavor) {
  RWMol *newM = nullptr;
  try {
    newM = FASTAToMol(pyObjectToString(seq), sanitize, flavor);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.message() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}
ROMol *MolFromHELM(python::object seq, bool sanitize) {
  RWMol *newM = nullptr;
  try {
    newM = HELMToMol(pyObjectToString(seq), sanitize);
  } catch (RDKit::FileParseException &e) {
    BOOST_LOG(rdWarningLog) << e.message() << std::endl;
  } catch (...) {
  }
  return static_cast<ROMol *>(newM);
}

struct smilesfrag_gen {
  std::string operator()(const ROMol &mol, const std::vector<int> &atomsToUse,
                         const std::vector<int> *bondsToUse,
                         const std::vector<std::string> *atomSymbols,
                         const std::vector<std::string> *bondSymbols,
                         bool doIsomericSmiles, bool doKekule, int rootedAtAtom,
                         bool canonical, bool allBondsExplicit,
                         bool allHsExplicit) {
    return MolFragmentToSmiles(
        mol, atomsToUse, bondsToUse, atomSymbols, bondSymbols, doIsomericSmiles,
        doKekule, rootedAtAtom, canonical, allBondsExplicit, allHsExplicit);
  }
};
struct cxsmilesfrag_gen {
  std::string operator()(const ROMol &mol, const std::vector<int> &atomsToUse,
                         const std::vector<int> *bondsToUse,
                         const std::vector<std::string> *atomSymbols,
                         const std::vector<std::string> *bondSymbols,
                         bool doIsomericSmiles, bool doKekule, int rootedAtAtom,
                         bool canonical, bool allBondsExplicit,
                         bool allHsExplicit) {
    return MolFragmentToCXSmiles(
        mol, atomsToUse, bondsToUse, atomSymbols, bondSymbols, doIsomericSmiles,
        doKekule, rootedAtAtom, canonical, allBondsExplicit, allHsExplicit);
  }
};

template <typename F>
std::string MolFragmentToSmilesHelper(
    const ROMol &mol, python::object atomsToUse, python::object bondsToUse,
    python::object atomSymbols, python::object bondSymbols,
    bool doIsomericSmiles, bool doKekule, int rootedAtAtom, bool canonical,
    bool allBondsExplicit, bool allHsExplicit) {
  std::unique_ptr<std::vector<int>> avect =
      pythonObjectToVect(atomsToUse, static_cast<int>(mol.getNumAtoms()));
  if (!avect.get() || !(avect->size())) {
    throw_value_error("atomsToUse must not be empty");
  }
  std::unique_ptr<std::vector<int>> bvect =
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

  std::string res =
      F()(mol, *avect.get(), bvect.get(), asymbols.get(), bsymbols.get(),
          doIsomericSmiles, doKekule, rootedAtAtom, canonical, allBondsExplicit,
          allHsExplicit);
  return res;
}

std::vector<unsigned int> CanonicalRankAtoms(const ROMol &mol,
                                             bool breakTies = true,
                                             bool includeChirality = true,
                                             bool includeIsotopes = true) {
  std::vector<unsigned int> ranks(mol.getNumAtoms());
  Canon::rankMolAtoms(mol, ranks, breakTies, includeChirality, includeIsotopes);
  return ranks;
}

std::vector<int> CanonicalRankAtomsInFragment(const ROMol &mol,
                                              python::object atomsToUse,
                                              python::object bondsToUse,
                                              python::object atomSymbols,
                                              python::object bondSymbols,
                                              bool breakTies = true)

{
  std::unique_ptr<std::vector<int>> avect =
      pythonObjectToVect(atomsToUse, static_cast<int>(mol.getNumAtoms()));
  if (!avect.get() || !(avect->size())) {
    throw_value_error("atomsToUse must not be empty");
  }
  std::unique_ptr<std::vector<int>> bvect =
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

  boost::dynamic_bitset<> atoms(mol.getNumAtoms());
  for (size_t i = 0; i < avect->size(); ++i) atoms[(*avect)[i]] = true;

  boost::dynamic_bitset<> bonds(mol.getNumBonds());
  for (size_t i = 0; bvect.get() && i < bvect->size(); ++i)
    bonds[(*bvect)[i]] = true;

  std::vector<unsigned int> ranks(mol.getNumAtoms());
  Canon::rankFragmentAtoms(mol, ranks, atoms, bonds, asymbols.get(),
                           bsymbols.get(), breakTies);

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
                           const SmilesParserParams &params) {
  std::string smiles = pyObjectToString(ismiles);

  return SmilesToMol(smiles, params);
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
#ifdef RDK_BUILD_COORDGEN_SUPPORT
void wrap_maesupplier();
#endif

// mol writer stuff
void wrap_smiwriter();
void wrap_sdwriter();
void wrap_tdtwriter();
void wrap_pdbwriter();

BOOST_PYTHON_MODULE(rdmolfiles) {
  std::string docString;

  python::scope().attr("__doc__") =
      "Module containing RDKit functionality for working with molecular file "
      "formats.";
  python::register_exception_translator<RDKit::BadFileException>(
      &rdBadFileExceptionTranslator);

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
              (python::arg("molBlock"), python::arg("sanitize") = true,
               python::arg("removeHs") = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "Construct a molecule from a Tripos Mol2 file.\n\n\
  NOTE:\n \
    The parser expects the atom-typing scheme used by Corina.\n\
    Atom types from Tripos' dbtranslate are less supported.\n\
    Other atom typing schemes are unlikely to work.\n\
\n\
  ARGUMENTS:\n                                  \
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
    - cleanupSubstructures: (optional) toggles standardizing some \n\
      substructures found in mol2 files.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";
  python::def("MolFromMol2File", RDKit::MolFromMol2File,
              (python::arg("molFileName"), python::arg("sanitize") = true,
               python::arg("removeHs") = true,
               python::arg("cleanupSubstructures") = true),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString =
      "Construct a molecule from a Tripos Mol2 block.\n\n\
  NOTE:\n \
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
              (python::arg("molBlock"), python::arg("sanitize") = true,
               python::arg("removeHs") = true,
               python::arg("cleanupSubstructures") = true),
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
  python::def("MolToMolBlock", RDKit::MolToMolBlock,
              (python::arg("mol"), python::arg("includeStereo") = true,
               python::arg("confId") = -1, python::arg("kekulize") = true,
               python::arg("forceV3000") = false),
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
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolToMolFile", RDKit::MolToMolFile,
      (python::arg("mol"), python::arg("filename"),
       python::arg("includeStereo") = true, python::arg("confId") = -1,
       python::arg("kekulize") = true, python::arg("forceV3000") = false),
      docString.c_str());

  python::class_<RDKit::SmilesParserParams, boost::noncopyable>(
      "SmilesParserParams", "Parameters controlling SMILES Parsing")
      .def_readwrite("maxIterations", &RDKit::SmilesParserParams::debugParse,
                     "controls the amount of debugging information produced")
      .def_readwrite("parseName", &RDKit::SmilesParserParams::parseName,
                     "controls whether or not the molecule name is also parsed")
      .def_readwrite(
          "allowCXSMILES", &RDKit::SmilesParserParams::allowCXSMILES,
          "controls whether or not the CXSMILES extensions are parsed")
      .def_readwrite("sanitize", &RDKit::SmilesParserParams::sanitize,
                     "controls whether or not the molecule is sanitized before "
                     "being returned")
      .def_readwrite("removeHs", &RDKit::SmilesParserParams::removeHs,
                     "controls whether or not Hs are removed before the "
                     "molecule is returned");
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
   the string no longer changes. It is the responsiblity of the caller to make sure \n\
   that substitutions results in legal and sensible SMILES. \n\
 \n\
   Examples of replacements: \n\
 \n\
     CC{Q}C with {'{Q}':'OCCO'} -> CCOCCOC  \n\
     C{A}C{Q}C with {'{Q}':'OCCO', '{A}':'C1(CC1)'} -> CC1(CC1)COCCOC  \n\
     C{A}C{Q}C with {'{Q}':'{X}CC{X}', '{A}':'C1CC1', '{X}':'N'} -> CC1CC1CCNCCNC  \n\
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
      "Returns the canonical SMILES string for a molecule\n\
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
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolToSmiles", RDKit::MolToSmiles,
      (python::arg("mol"), python::arg("isomericSmiles") = true,
       python::arg("kekuleSmiles") = false, python::arg("rootedAtAtom") = -1,
       python::arg("canonical") = true, python::arg("allBondsExplicit") = false,
       python::arg("allHsExplicit") = false, python::arg("doRandom") = false),
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
      to start at a particular atom. Defaults to -1.\n\
    - canonical: (optional) if false no attempt will be made to canonicalize\n\
      the molecule. Defaults to true.\n\
    - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
    - allHsExplicit: (optional) if true, all H counts will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
    - doRandom: (optional) if true, randomized the DFS transversal graph,\n\
      so we can generate random smiles. Defaults to false.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolFragmentToSmiles", MolFragmentToSmilesHelper<smilesfrag_gen>,
      (python::arg("mol"), python::arg("atomsToUse"),
       python::arg("bondsToUse") = 0, python::arg("atomSymbols") = 0,
       python::arg("bondSymbols") = 0, python::arg("isomericSmiles") = true,
       python::arg("kekuleSmiles") = false, python::arg("rootedAtAtom") = -1,
       python::arg("canonical") = true, python::arg("allBondsExplicit") = false,
       python::arg("allHsExplicit") = false),
      docString.c_str());

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
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolToCXSmiles", RDKit::MolToCXSmiles,
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
      to start at a particular atom. Defaults to -1.\n\
    - canonical: (optional) if false no attempt will be made to canonicalize\n\
      the molecule. Defaults to true.\n\
    - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
    - allHsExplicit: (optional) if true, all H counts will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
    - doRandom: (optional) if true, randomized the DFS transversal graph,\n\
      so we can generate random smiles. Defaults to false.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def(
      "MolFragmentToCXSmiles", MolFragmentToSmilesHelper<cxsmilesfrag_gen>,
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
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToSmarts", RDKit::MolToSmarts,
              (python::arg("mol"), python::arg("isomericSmiles") = true),
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
    - fileName: name of the file to read\n\
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
              (python::arg("molFileName"), python::arg("sanitize") = true,
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
            flavor & 1 : Write MODEL/ENDMDL lines around each record \n\
            flavor & 2 : Don't write any CONECT records \n\
            flavor & 4 : Write CONECT records in both directions \n\
            flavor & 8 : Don't use multiple CONECTs to encode bond order \n\
            flavor & 16 : Write MASTER record \n\
            flavor & 32 : Write TER record \n\
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
            flavor & 1 : Write MODEL/ENDMDL lines around each record \n\
            flavor & 2 : Don't write any CONECT records \n\
            flavor & 4 : Write CONECT records in both directions \n\
            flavor & 8 : Don't use multiple CONECTs to encode bond order \n\
            flavor & 16 : Write MASTER record \n\
            flavor & 32 : Write TER record \n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("MolToPDBFile", RDKit::MolToPDBFile,
              (python::arg("mol"), python::arg("filename"),
               python::arg("confId") = -1, python::arg("flavor") = 0),
              docString.c_str());

  docString =
      "Construct a molecule from a sequence string (currently only supports peptides).\n\n\
  ARGUMENTS:\n\
\n\
    - text: string containing the sequence\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
    - flavor: (optional)\n\
      0 Protein, L amino acids (default)\n\
      1 Protein, D amino acids\n\
      2 RNA, no cap\n\
      3 RNA, 5' cap\n\
      4 RNA, 3' cap\n\
      5 RNA, both caps\n\
      6 DNA, no cap\n\
      7 DNA, 5' cap\n\
      8 DNA, 3' cap\n\
      9 DNA, both caps\n\
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
      "Construct a molecule from a FASTA string (currently only supports peptides).\n\n\
  ARGUMENTS:\n\
\n\
    - text: string containing the FASTA\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to True.\n\
\n\
- flavor: (optional)\n\
  0 Protein, L amino acids (default)\n\
  1 Protein, D amino acids\n\
  2 RNA, no cap\n\
  3 RNA, 5' cap\n\
  4 RNA, 3' cap\n\
  5 RNA, both caps\n\
  6 DNA, no cap\n\
  7 DNA, 5' cap\n\
  8 DNA, 3' cap\n\
  9 DNA, both caps\n\
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
      "Construct a molecule from a HELM string (currently only supports peptides).\n\n\
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
  symmetry class.  From the perspective of the Molecular Graph, they are indentical.\n\
\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - breakTies: (optional) force breaking of ranked ties [default=True]\n\
    - includeChirality: (optional) use chiral information when computing rank [default=True]\n\
    - includeIsotopes: (optional) use isotope information when computing rank [default=True]\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("CanonicalRankAtoms", CanonicalRankAtoms,
              (python::arg("mol"), python::arg("breakTies") = true,
               python::arg("includeChirality") = true,
               python::arg("includeIsotopes") = true),
              docString.c_str());

  docString =
      "Returns the canonical atom ranking for each atom of a molecule fragment\n\
  See help(CanonicalRankAtoms) for more information.\n\
\n\
   >>> mol = MolFromSmiles('C1NCN1.C1NCN1')\n\
   >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(0,4), breakTies=False))\n\
   [0,1,0,1,-1,-1,-1,-1]\n\
   >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(4,8), breakTies=False))\n\
   [-1,-1,-1,-1,0,1,0,1]\n\
\n\
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
    - breakTies: (optional) force breaking of ranked ties\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";
  python::def("CanonicalRankAtomsInFragment", CanonicalRankAtomsInFragment,
              (python::arg("mol"), python::arg("atomsToUse"),
               python::arg("bondsToUse") = 0, python::arg("atomSymbols") = 0,
               python::arg("bondSymbols") = 0, python::arg("breakTies") = true),
              docString.c_str());

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
#ifdef RDK_BUILD_COORDGEN_SUPPORT
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
}
