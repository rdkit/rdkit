//
//  Copyright (C) 2010-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_FILEPARSERUTILS_H
#define RD_FILEPARSERUTILS_H

#include <string>
#include <iostream>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
class RWMol;
class Conformer;

namespace FileParserUtils {
template <typename T>
T stripSpacesAndCast(const std::string &input, bool acceptSpaces = false) {
  std::string trimmed = boost::trim_copy(input);
  if (acceptSpaces && trimmed == "") {
    return 0;
  } else {
    return boost::lexical_cast<T>(trimmed);
  }
}
RDKIT_FILEPARSERS_EXPORT int toInt(const std::string &input,
                                   bool acceptSpaces = false);
RDKIT_FILEPARSERS_EXPORT double toDouble(const std::string &input,
                                         bool acceptSpaces = true);

// reads a line from an MDL v3K CTAB
RDKIT_FILEPARSERS_EXPORT std::string getV3000Line(std::istream *inStream,
                                                  unsigned int &line);

// nAtoms and nBonds are ignored on input, set on output
RDKIT_FILEPARSERS_EXPORT bool ParseV3000CTAB(
    std::istream *inStream, unsigned int &line, RWMol *mol, Conformer *&conf,
    bool &chiralityPossible, unsigned int &nAtoms, unsigned int &nBonds,
    bool strictParsing = true, bool expectMEND = true);

// nAtoms and nBonds are used
RDKIT_FILEPARSERS_EXPORT bool ParseV2000CTAB(
    std::istream *inStream, unsigned int &line, RWMol *mol, Conformer *&conf,
    bool &chiralityPossible, unsigned int &nAtoms, unsigned int &nBonds,
    bool strictParsing = true);

RDKIT_FILEPARSERS_EXPORT Atom *replaceAtomWithQueryAtom(RWMol *mol, Atom *atom);

template <typename T>
void applyMolListPropsToAtoms(ROMol &mol, const std::string &prefix,
                              const std::string &missingValueMarker = "n/a") {
  for (auto pn : mol.getPropList()) {
    if (pn.find(prefix) == 0 && pn.size() > prefix.size()) {
      std::string strVect = mol.getProp<std::string>(pn);
      std::vector<std::string> tokens;
      boost::split(tokens, strVect, boost::is_any_of(" \t\n"),
                   boost::token_compress_on);
      if (tokens.size() < mol.getNumAtoms()) {
        BOOST_LOG(rdWarningLog)
            << "Property list " << pn << " too short, only " << tokens.size()
            << " elements found. Ignoring it." << std::endl;
        continue;
      }
      std::string mv = missingValueMarker;
      size_t first_token = 0;
      if (tokens.size() == mol.getNumAtoms() + 1 && tokens[0].front() == '[' &&
          tokens[0].back() == ']') {
        mv = tokens[0].substr(1);
        mv.pop_back();
        first_token = 1;
      }
      std::string atompn = pn.substr(prefix.size());
      for (size_t i = first_token; i < tokens.size(); ++i) {
        if (tokens[i] != mv) {
          T apv;
          try {
            apv = boost::lexical_cast<T>(tokens[i]);
          } catch (const boost::bad_lexical_cast &) {
            BOOST_LOG(rdWarningLog)
                << "Value " << tokens[i] << " for property " << pn
                << " of atom " << i << " can not be parsed. Ignoring it."
                << std::endl;
          }
          unsigned int atomid = i - first_token;
          mol.getAtomWithIdx(atomid)->setProp(atompn, apv);
        }
      }
    }
  }
}
inline void processMolPropertyLists(
    ROMol &mol, const std::string &missingValueMarker = "n/a") {
  applyMolListPropsToAtoms<std::string>(mol, "atom.prop.", missingValueMarker);
  applyMolListPropsToAtoms<std::int64_t>(mol, "atom.iprop.",
                                         missingValueMarker);
  applyMolListPropsToAtoms<double>(mol, "atom.dprop.", missingValueMarker);
  applyMolListPropsToAtoms<bool>(mol, "atom.bprop.", missingValueMarker);
}

template <typename T>
std::string getAtomPropertyList(ROMol &mol, const std::string &atomPropName, 
      std::string missingValueMarker = "", unsigned int lineSize=190 ){
    std::string res;
    std::string propVal;
    if(!missingValueMarker.empty()){
      propVal += boost::str(boost::format("[%s] ")%missingValueMarker);
    } else {
      missingValueMarker = "n/a";
    }
    for(const auto & atom : mol.atoms()) {
      std::string apVal=missingValueMarker;
      if(atom->hasProp(atomPropName)) {
        T tVal = atom->getProp<T>(atomPropName);
        apVal = boost::lexical_cast<std::string>(tVal);
        // seems like this should work, but it doesn't: atom->getProp(atomPropName,apVal);        
      }
      if(propVal.length() + apVal.length() +1 >= lineSize){
         // remove trailing space:
         propVal.pop_back();
         res += propVal+"\n";
         propVal = "";
      }
      propVal += apVal + " ";
    }
    if(!propVal.empty()) {
      // remove the trailing space:
      propVal.pop_back();
      res += propVal;
    }
    return res;
}
inline void createAtomIntPropertyList(ROMol &mol, const std::string &atomPropName, const std::string &missingValueMarker="", unsigned int lineSize=190){
  std::string molPropName = "atom.iprop."+atomPropName;
  mol.setProp(molPropName,getAtomPropertyList<boost::int64_t>(mol,atomPropName,missingValueMarker,lineSize));
}
inline void createAtomDoublePropertyList(ROMol &mol, const std::string &atomPropName, const std::string &missingValueMarker="", unsigned int lineSize=190){
  std::string molPropName = "atom.dprop."+atomPropName;
  mol.setProp(molPropName,getAtomPropertyList<double>(mol,atomPropName,missingValueMarker,lineSize));
}
inline void createAtomBoolPropertyList(ROMol &mol, const std::string &atomPropName, const std::string &missingValueMarker="", unsigned int lineSize=190){
  std::string molPropName = "atom.bprop."+atomPropName;
  mol.setProp(molPropName,getAtomPropertyList<bool>(mol,atomPropName,missingValueMarker,lineSize));
}
inline void createAtomStringPropertyList(ROMol &mol, const std::string &atomPropName, const std::string &missingValueMarker="", unsigned int lineSize=190){
  std::string molPropName = "atom.prop."+atomPropName;
  mol.setProp(molPropName,getAtomPropertyList<std::string>(mol,atomPropName,missingValueMarker,lineSize));
}




}  // namespace FileParserUtils
}  // namespace RDKit

#endif
