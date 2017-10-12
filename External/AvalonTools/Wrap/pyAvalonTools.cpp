// $Id: pyAvalonTools.cpp 4504 2011-03-04 13:57:07Z landrgr1 $
//
//  Created by Greg Landrum, July 2008
//

#include <boost/python.hpp>

#include <GraphMol/GraphMol.h>
#include <DataStructs/ExplicitBitVect.h>
#include <RDBoost/Wrap.h>
#include <AvalonTools.h>
#include <boost/cstdint.hpp>

extern "C" {
#include "struchk.h"
}

namespace python = boost::python;

namespace {

RDKit::SparseIntVect<boost::uint32_t> *getAvalonCountFP(const RDKit::ROMol &mol,
                                                        unsigned int nBits,
                                                        bool isQuery,
                                                        unsigned int bitFlags) {
  auto *res = new RDKit::SparseIntVect<boost::uint32_t>(nBits);
  AvalonTools::getAvalonCountFP(mol, *res, nBits, isQuery, bitFlags);
  return res;
}
RDKit::SparseIntVect<boost::uint32_t> *getAvalonCountFP(const std::string &data,
                                                        bool isSmiles,
                                                        unsigned int nBits,
                                                        bool isQuery,
                                                        unsigned int bitFlags) {
  auto *res = new RDKit::SparseIntVect<boost::uint32_t>(nBits);
  AvalonTools::getAvalonCountFP(data, isSmiles, *res, nBits, isQuery, bitFlags);
  return res;
}

ExplicitBitVect *getAvalonFP(const RDKit::ROMol &mol, unsigned int nBits,
                             bool isQuery, bool resetVect,
                             unsigned int bitFlags) {
  auto *res = new ExplicitBitVect(nBits);
  AvalonTools::getAvalonFP(mol, *res, nBits, isQuery, resetVect, bitFlags);
  return res;
}
python::list getAvalonFPAsWords(const RDKit::ROMol &mol, unsigned int nBits,
                                bool isQuery, bool resetVect,
                                unsigned int bitFlags) {
  std::vector<boost::uint32_t> words;
  AvalonTools::getAvalonFP(mol, words, nBits, isQuery, resetVect, bitFlags);
  python::list res;
  for (std::vector<boost::uint32_t>::const_iterator ci = words.begin();
       ci != words.end(); ++ci) {
    res.append(static_cast<unsigned long>(*ci));
  }
  return res;
}

ExplicitBitVect *getAvalonFP(const std::string &data, bool isSmiles,
                             unsigned int nBits, bool isQuery, bool resetVect,
                             unsigned int bitFlags) {
  auto *res = new ExplicitBitVect(nBits);
  AvalonTools::getAvalonFP(data, isSmiles, *res, nBits, isQuery, resetVect,
                           bitFlags);
  return res;
}
python::list getAvalonFPAsWords(const std::string &data, bool isSmiles,
                                unsigned int nBits, bool isQuery,
                                bool resetVect, unsigned int bitFlags) {
  std::vector<boost::uint32_t> words;
  AvalonTools::getAvalonFP(data, isSmiles, words, nBits, isQuery, resetVect,
                           bitFlags);
  python::list res;
  for (std::vector<boost::uint32_t>::const_iterator ci = words.begin();
       ci != words.end(); ++ci) {
    res.append(static_cast<unsigned long>(*ci));
  }
  return res;
}

python::tuple CheckMolecule(const std::string &data, bool isSmiles) {
  int errs = 0;
  RDKit::ROMOL_SPTR rMol = AvalonTools::checkMol(errs, data, isSmiles);
  return python::make_tuple(errs, rMol);
  ;
}
python::tuple CheckMolecule(RDKit::ROMol &mol) {
  int errs = 0;
  RDKit::ROMOL_SPTR rMol = AvalonTools::checkMol(errs, mol);
  return python::make_tuple(errs, rMol);
  ;
}
python::tuple CheckMoleculeString(const std::string &data, bool isSmiles) {
  std::pair<std::string, int> res = AvalonTools::checkMolString(data, isSmiles);
  return python::make_tuple(res.second, res.first);
  ;
}

enum StruChkFlag {
  bad_molecule = BAD_MOLECULE,
  alias_conversion_failed = ALIAS_CONVERSION_FAILED,
  transformed = TRANSFORMED,
  fragments_found = FRAGMENTS_FOUND,
  either_warning = EITHER_WARNING,
  stereo_error = STEREO_ERROR,
  dubious_stereo_removed = DUBIOUS_STEREO_REMOVED,
  atom_clash = ATOM_CLASH,
  atom_check_failed = ATOM_CHECK_FAILED,
  size_check_failed = SIZE_CHECK_FAILED,
  recharged = RECHARGED,
  stereo_forced_bad = STEREO_FORCED_BAD,
  stereo_transformed = STEREO_TRANSFORMED,
  template_transformed = TEMPLATE_TRANSFORMED,
};

enum StruChkResult {
  success = 0,
  bad_set = BAD_SET,
  transformed_set = TRANSFORMED_SET,
};
}

BOOST_PYTHON_MODULE(pyAvalonTools) {
  python::scope().attr("__doc__") =
      "Module containing functionality from the Avalon toolkit.\n\
\n\
The functions currently exposed are:\n\
  - GetCanonSmiles()   : return the canonical smiles for a molecule\n\
  - GetAvalonFP()      : return the Avalon fingerprint for a molecule as\n\
                         an RDKit ExplicitBitVector\n\
  - GetAvalonCountFP()      : return the Avalon fingerprint for a molecule as\n\
                              an RDKit SparseIntVector\n\
  - Generate2DCoords() : use the Avalon coordinate generator to create\n\
                         a set of 2D coordinates for a molecule\n\
Each function can be called with either an RDKit molecule or some\n\
molecule data as text (e.g. a SMILES or an MDL mol block).\n\
\n\
See the individual docstrings for more information.\n\
";

  std::string docString = "returns canonical smiles for an RDKit molecule";
  python::def("GetCanonSmiles",
              (std::string(*)(RDKit::ROMol &, int))AvalonTools::getCanonSmiles,
              (python::arg("mol"), python::arg("flags") = -1),
              docString.c_str());
  docString =
      "Returns canonical smiles for some molecule data.\n\
If the isSmiles argument is true, the data is assumed to be SMILES, otherwise\n\
MDL mol data is assumed.";
  python::def("GetCanonSmiles",
              (std::string(*)(const std::string &, bool,
                              int))AvalonTools::getCanonSmiles,
              (python::arg("molData"), python::arg("isSmiles"),
               python::arg("flags") = -1),
              docString.c_str());
  docString = "returns the Avalon fingerprint for an RDKit molecule";
  python::def("GetAvalonFP",
              (ExplicitBitVect * (*)(const RDKit::ROMol &, unsigned int, bool,
                                     bool, unsigned int))getAvalonFP,
              (python::arg("mol"), python::arg("nBits") = 512,
               python::arg("isQuery") = false, python::arg("resetVect") = false,
               python::arg("bitFlags") = AvalonTools::avalonSimilarityBits),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString =
      "returns the Avalon fingerprint for some molecule data.\n\
If the isSmiles argument is true, the data is assumed to be SMILES, otherwise\n\
MDL mol data is assumed.";
  python::def("GetAvalonFP",
              (ExplicitBitVect * (*)(const std::string &, bool, unsigned int,
                                     bool, bool, unsigned int))getAvalonFP,
              (python::arg("molData"), python::arg("isSmiles"),
               python::arg("nBits") = 512, python::arg("isQuery") = false,
               python::arg("resetVect") = false,
               python::arg("bitFlags") = AvalonTools::avalonSimilarityBits),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString = "Generates 2d coordinates for an RDKit molecule";
  python::def("Generate2DCoords",
              (unsigned int (*)(RDKit::ROMol &, bool))AvalonTools::set2DCoords,
              (python::arg("mol"), python::arg("clearConfs") = true),
              docString.c_str());
  docString =
      "returns an MDL mol block with 2D coordinates for some molecule data.\n\
If the isSmiles argument is true, the data is assumed to be SMILES, otherwise\n\
MDL mol data is assumed.";
  python::def(
      "Generate2DCoords",
      (std::string(*)(const std::string &, bool))AvalonTools::set2DCoords,
      (python::arg("molData"), python::arg("isSmiles")), docString.c_str());

  docString =
      "returns the Avalon fingerprint for an RDKit molecule as a list of ints";
  python::def("GetAvalonFPAsWords",
              (python::list(*)(const RDKit::ROMol &, unsigned int, bool, bool,
                               unsigned int))getAvalonFPAsWords,
              (python::arg("mol"), python::arg("nBits") = 512,
               python::arg("isQuery") = false, python::arg("resetVect") = false,
               python::arg("bitFlags") = AvalonTools::avalonSimilarityBits),
              docString.c_str());

  docString = "returns the Avalon count fingerprint for an RDKit molecule";
  python::def("GetAvalonCountFP", (RDKit::SparseIntVect<boost::uint32_t> *
                                   (*)(const RDKit::ROMol &, unsigned int, bool,
                                       unsigned int))getAvalonCountFP,
              (python::arg("mol"), python::arg("nBits") = 512,
               python::arg("isQuery") = false,
               python::arg("bitFlags") = AvalonTools::avalonSimilarityBits),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());
  docString =
      "returns the Avalon count fingerprint for some molecule data.\n\
If the isSmiles argument is true, the data is assumed to be SMILES, otherwise\n\
MDL mol data is assumed.";
  python::def("GetAvalonCountFP", (RDKit::SparseIntVect<boost::uint32_t> *
                                   (*)(const std::string &, bool, unsigned int,
                                       bool, unsigned int))getAvalonCountFP,
              (python::arg("molData"), python::arg("isSmiles"),
               python::arg("nBits") = 512, python::arg("isQuery") = false,
               python::arg("bitFlags") = AvalonTools::avalonSimilarityBits),
              docString.c_str(),
              python::return_value_policy<python::manage_new_object>());

  docString =
      "returns the Avalon fingerprint for some molecule data as a list of ints.\n\
If the isSmiles argument is true, the data is assumed to be SMILES, otherwise\n\
MDL mol data is assumed.";
  python::def("GetAvalonFPAsWords",
              (python::list(*)(const std::string &, bool, unsigned int, bool,
                               bool, unsigned int))getAvalonFPAsWords,
              (python::arg("molData"), python::arg("isSmiles"),
               python::arg("nBits") = 512, python::arg("isQuery") = false,
               python::arg("resetVect") = false,
               python::arg("bitFlags") = AvalonTools::avalonSimilarityBits),
              docString.c_str());

  docString =
      "initializes the structure checker.\n\
The argument should contain option lines separated by embedded newlines.\
An empty string will be used if the argument is omitted.\
An non-zero error code is returned in case of failure.";
  python::def("InitializeCheckMol",
              (int (*)(const std::string &))AvalonTools::initCheckMol,
              (python::arg("options") = ""), docString.c_str());
  docString = "close open files used by molecule-checking functions.";
  python::def("CloseCheckMolFiles", AvalonTools::closeCheckMolFiles,
              docString.c_str());

  docString =
      "check a molecule passed in as a string.\n\
If the isSmiles argument is true, the string should represent the SMILES encoding\n\
of the molecule, otherwise it should be encoded as an MDL molfile.\n\
The first member of the return tuple contains the bit-encoded corrections made to the molecule.\n\
If possible, the molecule (corrected when appropriate) is returned as the second member of \n\
the return tuple. Otherwise, None is returned.";
  python::def("CheckMolecule",
              (python::tuple(*)(const std::string &, bool))CheckMolecule,
              (python::arg("molstring"), python::arg("isSmiles")),
              docString.c_str());

  docString =
      "check a molecule passed in as an RDKit molecule.\n\
The first member of the return tuple contains the bit-encoded corrections made to the molecule.\n\
If possible, the molecule (corrected when appropriate) is returned as the second member of \n\
the return tuple. Otherwise, None is returned.";
  python::def("CheckMolecule", (python::tuple(*)(RDKit::ROMol &))CheckMolecule,
              (python::arg("mol")), docString.c_str());

  docString =
      "check a molecule passed in as a string and returns the result as a string.\n\
If the isSmiles argument is true, the string should represent the SMILES encoding\n\
of the molecule, otherwise it should be encoded as an MDL molfile.\n\
The first member of the return tuple contains the bit-encoded corrections made to the molecule.\n\
If possible, a corrected CTAB for the molecule is returned as the second member of \n\
the return tuple.";
  python::def("CheckMoleculeString", CheckMoleculeString,
              (python::arg("molstring"), python::arg("isSmiles")),
              docString.c_str());

  python::def("GetCheckMolLog", AvalonTools::getCheckMolLog,
              "Returns the Struchk log for the last molecules processed.");
  
  python::scope().attr("avalonSSSBits") = AvalonTools::avalonSSSBits;
  python::scope().attr("avalonSimilarityBits") =
      AvalonTools::avalonSimilarityBits;

  python::enum_<StruChkFlag>("StruChkFlag")
      .value("bad_molecule", bad_molecule)
      .value("alias_conversion_failed", alias_conversion_failed)
      .value("transformed", transformed)
      .value("fragments_found", fragments_found)
      .value("either_warning", either_warning)
      .value("stereo_error", stereo_error)
      .value("dubious_stereo_removed", dubious_stereo_removed)
      .value("atom_clash", atom_clash)
      .value("atom_check_failed", atom_check_failed)
      .value("size_check_failed", size_check_failed)
      .value("recharged", recharged)
      .value("stereo_forced_bad", stereo_forced_bad)
      .value("stereo_transformed", stereo_transformed)
      .value("template_transformed", template_transformed);
  python::enum_<StruChkResult>("StruChkResult")
      .value("success", success)
      .value("bad_set", bad_set)
      .value("transformed_set", transformed_set);
}
