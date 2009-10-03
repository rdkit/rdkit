// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "rdmolops.h"
#include <fstream>
#include <boost/python.hpp>

#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/BadFileException.h>

#include <RDBoost/Wrap.h>
#include <RDBoost/Exceptions.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SanitException.h>

namespace python = boost::python;
using namespace RDKit;

void rdSanitExceptionTranslator(RDKit::MolSanitizeException const& x){
  std::ostringstream ss;
  ss << "Sanitization error: " << x.message();
  PyErr_SetString(PyExc_ValueError,ss.str().c_str());
}
void rdBadFileExceptionTranslator(RDKit::BadFileException const& x){
  std::ostringstream ss;
  ss << "File error: " << x.message();
  PyErr_SetString(PyExc_IOError,ss.str().c_str());
}


namespace RDKit{
  ROMol *MolFromSmiles(std::string smiles,bool sanitize=1){
    RWMol *newM;
    try {
      newM = SmilesToMol(smiles,0,sanitize);
    } catch (...) {
      newM=0;
    }
    return static_cast<ROMol *>(newM);
    if(!newM) return 0;
  }

  ROMol *MolFromSmarts(const char *smarts,bool mergeHs=false){
    RWMol *newM; 
    try {
      newM = SmartsToMol(smarts,0,mergeHs);
    } catch (...) {
      newM=0;
    }
    return static_cast<ROMol *>(newM);
  }
   
  ROMol *MolFromMolFile(const char *molFilename, bool sanitize=true, bool removeHs=true) {
    // sometimes problems with exception-passing between modules
    // causes the BadFileException to get lost. This is something
    // of a hack to route around that.
    std::ifstream inStream(molFilename);
    if (!inStream || (inStream.bad()) ) {
      std::ostringstream errout;
      errout << "Bad input file " << molFilename;
      PyErr_SetString(PyExc_IOError,errout.str().c_str());
      throw python::error_already_set();
    }
    
    RWMol *newM=0;
    try {
      newM = MolFileToMol(molFilename, sanitize,removeHs);
    } catch (const RDKit::BadFileException &e) {
      PyErr_SetString(PyExc_IOError,e.message());
      throw python::error_already_set();
    }
    catch (...) {
      newM=0;
    }
    return static_cast<ROMol *>(newM);
  }

  ROMol *MolFromMolBlock(std::string molBlock, bool sanitize=true, bool removeHs=true) {
    std::istringstream inStream(molBlock);
    unsigned int line = 0;
    RWMol *newM;
    try {
      newM = MolDataStreamToMol(inStream, line, sanitize, removeHs);
    } catch (...) {
      newM=0;
    }
    return static_cast<ROMol *>(newM);
  }
}

// MolSupplier stuff
#ifdef SUPPORT_COMPRESSED_SUPPLIERS
void wrap_compressedsdsupplier();
#endif
void wrap_sdsupplier();
void wrap_smisupplier();

// mol writer stuff
void wrap_smiwriter();
void wrap_sdwriter();


BOOST_PYTHON_MODULE(rdmolfiles)
{
  std::string docString;

  python::scope().attr("__doc__") =
    "Module containing RDKit functionality for working with molecular file formats."
    ;
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  python::register_exception_translator<RDKit::MolSanitizeException>(&rdSanitExceptionTranslator);

  docString="Construct a molecule from a Mol file.\n\n\
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
  python::def("MolFromMolFile", RDKit::MolFromMolFile,
	      (python::arg("molFileName"),
	       python::arg("sanitize")=true,
               python::arg("removeHs")=true),
	      docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Construct a molecule from a Mol block.\n\n\
  ARGUMENTS:\n\
\n\
    - molBlock: string containing the Mol block\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to 1.\n\
\n\
    - removeHs: (optional) toggles removing hydrogens from the molecule.\n\
      This only make sense when sanitization is done.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";  
  python::def("MolFromMolBlock", RDKit::MolFromMolBlock,
	      (python::arg("molBlock"),
	       python::arg("sanitize")=true,
	       python::arg("removeHs")=true),
	      docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Returns the a Mol block for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - includeStereo: (optional) toggles inclusion of stereochemical\n\
                     information in the output\n\
    - confId: (optional) selects which conformation to output (-1 = default)\n\
    - kekulize: (optional) triggers kekulization of the molecule before it's written,\n\
                as suggested by the MDL spec.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";  
  python::def("MolToMolBlock",RDKit::MolToMolBlock,
	      (python::arg("mol"),python::arg("includeStereo")=false,
	       python::arg("confId")=-1,python::arg("kekulize")=true),
	      docString.c_str());


  docString="Construct a molecule from a SMILES string.\n\n\
  ARGUMENTS:\n\
\n\
    - SMILES: the smiles string\n\
\n\
    - sanitize: (optional) toggles sanitization of the molecule.\n\
      Defaults to 1.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";  
  python::def("MolFromSmiles",RDKit::MolFromSmiles,
	      (python::arg("SMILES"),
	       python::arg("sanitize")=true),
	      docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Construct a molecule from a SMARTS string.\n\n\
  ARGUMENTS:\n\
\n\
    - SMARTS: the smarts string\n\
\n\
    - mergeHs: (optional) toggles the merging of explicit Hs in the query into the attached\n\
      atoms.  So, for example, 'C[H]' becomes '[C;!H0]'.\n\
      Defaults to 0.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";  
  python::def("MolFromSmarts",RDKit::MolFromSmarts,
	      (python::arg("SMARTS"),
	       python::arg("mergeHs")=false),
	      docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Returns the canonical SMILES string for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - isomericSmiles: (optional) include information about stereochemistry in\n\
      the SMILES.  Defaults to false.\n\
    - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in\n\
      the SMILES.  Defaults to false.\n\
    - rootedAtAtom: (optional) if non-negative, this forces the SMILES \n\
      to start at a particular atom. Defaults to -1.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";  
  python::def("MolToSmiles",RDKit::MolToSmiles,
	      (python::arg("mol"),
	       python::arg("isomericSmiles")=false,
	       python::arg("kekuleSmiles")=false,
	       python::arg("rootedAtAtom")=-1),
	      docString.c_str());

  docString="Returns a SMARTS string for a molecule\n\
  ARGUMENTS:\n\
\n\
    - mol: the molecule\n\
    - isomericSmarts: (optional) include information about stereochemistry in\n\
      the SMARTS.  Defaults to false.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";  
  python::def("MolToSmarts",RDKit::MolToSmarts,
	      (python::arg("mol"),python::arg("isomericSmiles")=false),
	      docString.c_str());
  
  /********************************************************
   * MolSupplier stuff
   *******************************************************/
#ifdef SUPPORT_COMPRESSED_SUPPLIERS
  wrap_compressedsdsupplier();
#endif
  wrap_sdsupplier();
  wrap_smisupplier();

  /********************************************************
   * MolWriter stuff
   *******************************************************/
  wrap_smiwriter();
  wrap_sdwriter();
}


