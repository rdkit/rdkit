// $Id$
//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "rdmolops.h"
#include <boost/python.hpp>

#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>

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
  ROMol *MolFromSmiles(std::string smiles,bool sanitize,
                       python::dict replDict){
    std::map<std::string,std::string> replacements;
    for(unsigned int i=0;i<python::extract<unsigned int>(replDict.keys().attr("__len__")());++i){
      replacements[python::extract<std::string>(replDict.keys()[i])]=python::extract<std::string>(replDict.values()[i]);
    }
    RWMol *newM;
    try {
      newM = SmilesToMol(smiles,0,sanitize,&replacements);
    } catch (...) {
      newM=0;
    }
    return static_cast<ROMol *>(newM);
  }

  ROMol *MolFromSmarts(const char *smarts,bool mergeHs,
                       python::dict replDict){
    std::map<std::string,std::string> replacements;
    for(unsigned int i=0;i<python::extract<unsigned int>(replDict.keys().attr("__len__")());++i){
      replacements[python::extract<std::string>(replDict.keys()[i])]=python::extract<std::string>(replDict.values()[i]);
    }
    RWMol *newM; 
    try {
      newM = SmartsToMol(smarts,0,mergeHs,&replacements);
    } catch (...) {
      newM=0;
    }
    return static_cast<ROMol *>(newM);
  }
   
  ROMol *MolFromTPLFile(const char *filename, bool sanitize=true,
			bool skipFirstConf=false ) {
    RWMol *newM;
    try {
      newM = TPLFileToMol(filename,sanitize,skipFirstConf);
    } catch (RDKit::BadFileException &e) {
      PyErr_SetString(PyExc_IOError,e.message());
      throw python::error_already_set();
    } catch (...) {
      newM=0;
    }
    return static_cast<ROMol *>(newM);
  }

  ROMol *MolFromTPLBlock(std::string tplBlock, bool sanitize=true,
			bool skipFirstConf=false ) {
    std::istringstream inStream(tplBlock);
    unsigned int line = 0;
    RWMol *newM;
    try {
      newM = TPLDataStreamToMol(&inStream,line,sanitize,skipFirstConf);
    } catch (...) {
      newM=0;
    }
    return static_cast<ROMol *>(newM);
  }

  ROMol *MolFromMolFile(const char *molFilename, bool sanitize, bool removeHs,bool strictParsing) {
    RWMol *newM=0;
    try {
      newM = MolFileToMol(molFilename, sanitize,removeHs,strictParsing);
    } catch (RDKit::BadFileException &e) {
      PyErr_SetString(PyExc_IOError,e.message());
      throw python::error_already_set();
    } catch (RDKit::FileParseException &e) {
      BOOST_LOG(rdWarningLog) << e.message() <<std::endl;
    } catch (...) {

    }
    return static_cast<ROMol *>(newM);
  }

  ROMol *MolFromMolBlock(std::string molBlock, bool sanitize, bool removeHs, bool strictParsing) {
    std::istringstream inStream(molBlock);
    unsigned int line = 0;
    RWMol *newM=0;
    try {
      newM = MolDataStreamToMol(inStream, line, sanitize, removeHs, strictParsing);
    }  catch (RDKit::FileParseException &e) {
      BOOST_LOG(rdWarningLog) << e.message() <<std::endl;
    } catch (...) {
    }
    return static_cast<ROMol *>(newM);
  }

  ROMol *MolFromMol2File(const char *molFilename, bool sanitize=true, bool removeHs=true) {
    RWMol *newM;
    try {
      newM = Mol2FileToMol(molFilename, sanitize,removeHs);
    } catch (RDKit::BadFileException &e) {
      PyErr_SetString(PyExc_IOError,e.message());
      throw python::error_already_set();
    } catch (...) {
      newM=0;
    }
    return static_cast<ROMol *>(newM);
  }

  ROMol *MolFromMol2Block(std::string mol2Block, bool sanitize=true, bool removeHs=true){
    std::istringstream inStream(mol2Block);
    RWMol *newM;
    try {
      newM = Mol2DataStreamToMol(inStream, sanitize, removeHs);
    } catch (...) {
      newM=0;
    }
    return static_cast<ROMol *>(newM);
  }

  std::string MolFragmentToSmilesHelper(const ROMol &mol,
                                        python::object atomsToUse,
                                        python::object bondsToUse,
                                        python::object atomSymbols,
                                        python::object bondSymbols,
                                        bool doIsomericSmiles,
                                        bool doKekule,
                                        int rootedAtAtom,
                                        bool canonical,
                                        bool allBondsExplicit
                                        ){
    std::vector<int> *avect=pythonObjectToVect(atomsToUse,static_cast<int>(mol.getNumAtoms()));
    if(!avect || !(avect->size())){
      throw_value_error("atomsToUse must not be empty");
    }
    std::vector<int> *bvect=pythonObjectToVect(bondsToUse,static_cast<int>(mol.getNumBonds()));
    std::vector<std::string> *asymbols=pythonObjectToVect<std::string>(atomSymbols);
    std::vector<std::string> *bsymbols=pythonObjectToVect<std::string>(bondSymbols);
    if(asymbols && asymbols->size()!=mol.getNumAtoms()){
      throw_value_error("length of atom symbol list != number of atoms");
    }
    if(bsymbols && bsymbols->size()!=mol.getNumBonds()){
      throw_value_error("length of bond symbol list != number of bonds");
    }
    
    std::string res=MolFragmentToSmiles(mol,*avect,bvect,asymbols,bsymbols,
                                        doIsomericSmiles,doKekule,rootedAtAtom,
                                        canonical,allBondsExplicit);
    if(avect) delete avect;
    if(bvect) delete bvect;
    if(asymbols) delete asymbols;
    if(bsymbols) delete bsymbols;
    return res;
  }

}

// MolSupplier stuff
#ifdef SUPPORT_COMPRESSED_SUPPLIERS
void wrap_compressedsdsupplier();
#endif
void wrap_sdsupplier();
void wrap_forwardsdsupplier();
void wrap_tdtsupplier();
void wrap_smisupplier();

// mol writer stuff
void wrap_smiwriter();
void wrap_sdwriter();
void wrap_tdtwriter();


BOOST_PYTHON_MODULE(rdmolfiles)
{
  std::string docString;

  python::scope().attr("__doc__") =
    "Module containing RDKit functionality for working with molecular file formats."
    ;
  python::register_exception_translator<IndexErrorException>(&translate_index_error);
  python::register_exception_translator<ValueErrorException>(&translate_value_error);
  python::register_exception_translator<RDKit::MolSanitizeException>(&rdSanitExceptionTranslator);
  python::register_exception_translator<RDKit::BadFileException>(&rdBadFileExceptionTranslator);


  docString="Construct a molecule from a TPL file.\n\n\
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
	      (python::arg("fileName"),
	       python::arg("sanitize")=true,
	       python::arg("skipFirstConf")=false),
	      docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Construct a molecule from a TPL block.\n\n\
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
	      (python::arg("tplBlock"),
	       python::arg("sanitize")=true,
	       python::arg("skipFirstConf")=false),
	      docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

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
    - strictParsing: (optional) if this is false, the parser is more lax about.\n\
      correctness of the content.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";  
  python::def("MolFromMolFile", RDKit::MolFromMolFile,
	      (python::arg("molFileName"),
	       python::arg("sanitize")=true,
               python::arg("removeHs")=true,
               python::arg("strictParsing")=true),
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
    - strictParsing: (optional) if this is false, the parser is more lax about.\n\
      correctness of the content.\n\
      Defaults to true.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";  
  python::def("MolFromMolBlock", RDKit::MolFromMolBlock,
	      (python::arg("molBlock"),
	       python::arg("sanitize")=true,
	       python::arg("removeHs")=true,
               python::arg("strictParsing")=true),
	      docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());

  docString="Construct a molecule from a Tripos Mol2 file.\n\n\
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
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";  
  python::def("MolFromMol2File", RDKit::MolFromMol2File,
	      (python::arg("molFileName"),
	       python::arg("sanitize")=true,
               python::arg("removeHs")=true),
	      docString.c_str(),
	      python::return_value_policy<python::manage_new_object>());
  docString="Construct a molecule from a Tripos Mol2 block.\n\n\
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
  python::def("MolFromMol2Block", RDKit::MolFromMol2Block,
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
  python::def("MolFromSmiles",RDKit::MolFromSmiles,
	      (python::arg("SMILES"),
	       python::arg("sanitize")=true,
               python::arg("replacements")=python::dict()),
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
    - replacements: (optional) a dictionary of replacement strings (see below)\n\
      Defaults to {}. See the documentation for MolFromSmiles for an explanation.\n\
\n\
  RETURNS:\n\
\n\
    a Mol object, None on failure.\n\
\n";  
  python::def("MolFromSmarts",RDKit::MolFromSmarts,
	      (python::arg("SMARTS"),
	       python::arg("mergeHs")=false,
               python::arg("replacements")=python::dict()),
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
    - canonical: (optional) if false no attempt will be made to canonicalize\n\
      the molecule. Defaults to true.\n\
    - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";  
  python::def("MolToSmiles",RDKit::MolToSmiles,
	      (python::arg("mol"),
	       python::arg("isomericSmiles")=false,
	       python::arg("kekuleSmiles")=false,
	       python::arg("rootedAtAtom")=-1,
	       python::arg("canonical")=true,
               python::arg("allBondsExplicit")=false),
	      docString.c_str());

  docString="Returns the canonical SMILES string for a fragment of a molecule\n\
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
      the SMILES.  Defaults to false.\n\
    - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in\n\
      the SMILES.  Defaults to false.\n\
    - rootedAtAtom: (optional) if non-negative, this forces the SMILES \n\
      to start at a particular atom. Defaults to -1.\n\
    - canonical: (optional) if false no attempt will be made to canonicalize\n\
      the molecule. Defaults to true.\n\
    - allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated\n\
      in the output SMILES. Defaults to false.\n\
\n\
  RETURNS:\n\
\n\
    a string\n\
\n";  
  python::def("MolFragmentToSmiles",MolFragmentToSmilesHelper,
	      (python::arg("mol"),
               python::arg("atomsToUse"),
               python::arg("bondsToUse")=0,
               python::arg("atomSymbols")=0,
               python::arg("bondSymbols")=0,
	       python::arg("isomericSmiles")=false,
	       python::arg("kekuleSmiles")=false,
	       python::arg("rootedAtAtom")=-1,
	       python::arg("canonical")=true,
               python::arg("allBondsExplicit")=false),
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

  docString="Writes a molecule to a TPL file.\n\n\
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
	      (python::arg("mol"),
	       python::arg("fileName"),
	       python::arg("partialChargeProp")="_GasteigerCharge",
	       python::arg("writeFirstConfTwice")=false),
	      docString.c_str());

  docString="Returns the Tpl block for a molecule.\n\n\
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
	       python::arg("partialChargeProp")="_GasteigerCharge",
	       python::arg("writeFirstConfTwice")=false),
	      docString.c_str());


  

  
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

  /********************************************************
   * MolWriter stuff
   *******************************************************/
  wrap_smiwriter();
  wrap_sdwriter();
  wrap_tdtwriter();


  

}


