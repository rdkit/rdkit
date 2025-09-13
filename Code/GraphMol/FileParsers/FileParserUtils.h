//
//  Copyright (C) 2010-2022 Greg Landrum and other RDKit contributors
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
#include "FileParsers.h"
#include <string_view>

namespace RDKit {
class RWMol;
class Conformer;

namespace FileParserUtils {
RDKIT_FILEPARSERS_EXPORT inline std::string_view strip(
    std::string_view orig, std::string stripChars = " \t\r\n") {
  std::string_view res = orig;
  auto start = res.find_first_not_of(stripChars);
  if (start != std::string_view::npos) {
    auto end = res.find_last_not_of(stripChars) + 1;
    res = res.substr(start, end - start);
  } else {
    res = "";
  }
  return res;
}

template <typename T>
T stripSpacesAndCast(std::string_view input, bool acceptSpaces = false) {
  auto trimmed = strip(input, " ");
  if (acceptSpaces && trimmed.empty()) {
    return 0;
  } else {
    return boost::lexical_cast<T>(trimmed);
  }
}
template <typename T>
T stripSpacesAndCast(const std::string &input, bool acceptSpaces = false) {
  return stripSpacesAndCast<T>(std::string_view(input.c_str()), acceptSpaces);
}
RDKIT_FILEPARSERS_EXPORT int toInt(const std::string &input,
                                   bool acceptSpaces = true);
RDKIT_FILEPARSERS_EXPORT unsigned int toUnsigned(const std::string &input,
                                                 bool acceptSpaces = true);
RDKIT_FILEPARSERS_EXPORT double toDouble(const std::string &input,
                                         bool acceptSpaces = true);
RDKIT_FILEPARSERS_EXPORT int toInt(const std::string_view input,
                                   bool acceptSpaces = true);
RDKIT_FILEPARSERS_EXPORT unsigned int toUnsigned(std::string_view input,
                                                 bool acceptSpaces = true);
RDKIT_FILEPARSERS_EXPORT double toDouble(const std::string_view input,
                                         bool acceptSpaces = true);

// gets a V3000 CTAB for a molecule
RDKIT_FILEPARSERS_EXPORT std::string getV3000CTAB(
    const ROMol &tmol, const boost::dynamic_bitset<> &wasAromatic,
    int confId = -1, unsigned int precision = 6);
//! \overload
inline std::string getV3000CTAB(const ROMol &tmol, int confId = -1,
                                unsigned int precision = 6) {
  boost::dynamic_bitset<> wasAromatic(tmol.getNumBonds());
  return getV3000CTAB(tmol, wasAromatic, confId, precision);
};
// reads a line from an MDL v3K CTAB
RDKIT_FILEPARSERS_EXPORT std::string getV3000Line(std::istream *inStream,
                                                  unsigned int &line);

// nAtoms and nBonds are ignored on input, set on output
RDKIT_FILEPARSERS_EXPORT bool ParseV3000CTAB(
    std::istream *inStream, unsigned int &line, RWMol *mol, Conformer *&conf,
    bool &chiralityPossible, unsigned int &nAtoms, unsigned int &nBonds,
    bool strictParsing = true, bool expectMEND = true,
    bool expectMacroAtoms = false);

// nAtoms and nBonds are used
RDKIT_FILEPARSERS_EXPORT bool ParseV2000CTAB(
    std::istream *inStream, unsigned int &line, RWMol *mol, Conformer *&conf,
    bool &chiralityPossible, unsigned int &nAtoms, unsigned int &nBonds,
    bool strictParsing = true);

//! finishes up the processing (sanitization, etc.) of a molecule read from
//! CTAB
RDKIT_FILEPARSERS_EXPORT void finishMolProcessing(
    RWMol *res, bool chiralityPossible,
    const v2::FileParsers::MolFileParserParams &ps);
//! \overload
inline void finishMolProcessing(RWMol *res, bool chiralityPossible,
                                bool sanitize, bool removeHs) {
  v2::FileParsers::MolFileParserParams ps;
  ps.sanitize = sanitize;
  ps.removeHs = removeHs;
  finishMolProcessing(res, chiralityPossible, ps);
}

//! Deprecated, please use QueryOps::replaceAtomWithQueryAtom instead
RDKIT_FILEPARSERS_EXPORT Atom *replaceAtomWithQueryAtom(RWMol *mol, Atom *atom);

//! applies a particular property to the atoms as an atom property list
template <typename T, typename U>
void applyMolListProp(ROMol &mol, const std::string &pn,
                      const std::string &prefix,
                      const std::string &missingValueMarker, size_t nItems,
                      U getter) {
  std::string itempn = pn.substr(prefix.size());
  std::string strVect = mol.getProp<std::string>(pn);
  std::vector<std::string> tokens;
  boost::split(tokens, strVect, boost::is_any_of(" \t\n"),
               boost::token_compress_on);
  if (tokens.size() < nItems) {
    BOOST_LOG(rdWarningLog)
        << "Property list " << pn << " too short, only " << tokens.size()
        << " elements found. Ignoring it." << std::endl;
    return;
  }
  std::string mv = missingValueMarker;
  size_t first_token = 0;
  if (tokens.size() == nItems + 1 && tokens[0].front() == '[' &&
      tokens[0].back() == ']') {
    mv = std::string(tokens[0].begin() + 1, tokens[0].end() - 1);
    first_token = 1;
  }
  if (mv.empty()) {
    BOOST_LOG(rdWarningLog) << "Missing value marker for property " << pn
                            << " is empty." << std::endl;
  }
  for (size_t i = first_token; i < tokens.size(); ++i) {
    if (tokens[i] != mv) {
      unsigned int itemid = i - first_token;
      try {
        T apv = boost::lexical_cast<T>(tokens[i]);
        getter(itemid)->setProp(itempn, apv);
      } catch (const boost::bad_lexical_cast &) {
        BOOST_LOG(rdWarningLog)
            << "Value " << tokens[i] << " for property " << pn << " of item "
            << itemid << " can not be parsed. Ignoring it." << std::endl;
      }
    }
  }
}

//! applies a particular property to the atoms as an atom property list
template <typename T>
[[deprecated("use applyMolListProp instead")]]
void applyMolListPropToAtoms(ROMol &mol, const std::string &pn,
                             const std::string &prefix,
                             const std::string &missingValueMarker = "n/a") {
  auto getter = [&mol](size_t which) { return mol.getAtomWithIdx(which); };
  applyMolListProp<T>(mol, pn, prefix, missingValueMarker, mol.getNumAtoms(),
                      getter);
}

template <typename T, typename U>
void applyMolListProps(ROMol &mol, const std::string &prefix, size_t nItems,
                       U getter, const std::string missingValueMarker = "n/a") {
  for (auto pn : mol.getPropList()) {
    if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
      applyMolListProp<T>(mol, pn, prefix, missingValueMarker, nItems, getter);
    }
  }
}

//! applies all properties matching a particular prefix as an atom property
//! list
template <typename T>
[[deprecated("use applyMolListProps instead")]]
void applyMolListPropsToAtoms(ROMol &mol, const std::string &prefix,
                              const std::string missingValueMarker = "n/a") {
  auto getter = [&mol](size_t which) { return mol.getAtomWithIdx(which); };
  auto nItems = mol.getNumAtoms();
  for (auto pn : mol.getPropList()) {
    if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
      applyMolListProp<T>(mol, pn, prefix, missingValueMarker, nItems, getter);
    }
  }
}

static const std::string atomPropPrefix = "atom.";
static const std::string bondPropPrefix = "bond.";

//! if the property name matches our rules for atom property lists, we'll
//! apply it to the atoms
inline void processMolPropertyList(
    ROMol &mol, const std::string pn,
    const std::string &missingValueMarker = "n/a") {
  auto propSetter = [&](const std::string &propPrefix, auto getter,
                        size_t nItems) {
    std::string prefix = propPrefix + "prop.";
    if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
      applyMolListProp<std::string>(mol, pn, prefix, missingValueMarker, nItems,
                                    getter);
    } else {
      prefix = propPrefix + "iprop.";
      if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
        applyMolListProp<std::int64_t>(mol, pn, prefix, missingValueMarker,
                                       nItems, getter);
      } else {
        prefix = propPrefix + "dprop.";
        if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
          applyMolListProp<double>(mol, pn, prefix, missingValueMarker, nItems,
                                   getter);
        } else {
          prefix = propPrefix + "bprop.";
          if (pn.find(prefix) == 0 && pn.length() > prefix.length()) {
            applyMolListProp<bool>(mol, pn, prefix, missingValueMarker, nItems,
                                   getter);
          }
        }
      }
    }
  };

  if (pn.find(atomPropPrefix) == 0 && pn.length() > atomPropPrefix.length()) {
    propSetter(
        atomPropPrefix,
        [&mol](size_t which) { return mol.getAtomWithIdx(which); },
        mol.getNumAtoms());
  } else if (pn.find(bondPropPrefix) == 0 &&
             pn.length() > bondPropPrefix.length()) {
    propSetter(
        bondPropPrefix,
        [&mol](size_t which) { return mol.getBondWithIdx(which); },
        mol.getNumBonds());
  }
}
//! loops over all properties and applies the ones that match the rules for
//! atom property lists to the atoms and bonds
inline void processMolPropertyLists(
    ROMol &mol, const std::string &missingValueMarker = "n/a") {
  for (auto pn : mol.getPropList()) {
    processMolPropertyList(mol, pn, missingValueMarker);
  }
}
template <typename T>
std::string getAtomPropertyList(ROMol &mol, const std::string &atomPropName,
                                std::string missingValueMarker = "",
                                unsigned int lineSize = 190) {
  std::string res;
  std::string propVal;
  if (!missingValueMarker.empty()) {
    propVal += boost::str(boost::format("[%s] ") % missingValueMarker);
  } else {
    missingValueMarker = "n/a";
  }
  for (const auto &atom : mol.atoms()) {
    std::string apVal = missingValueMarker;
    if (atom->hasProp(atomPropName)) {
      T tVal = atom->getProp<T>(atomPropName);
      apVal = boost::lexical_cast<std::string>(tVal);
      // seems like this should work, but it doesn't:
      // atom->getProp(atomPropName,apVal);
    }
    if (propVal.length() + apVal.length() + 1 >= lineSize) {
      // remove trailing space:
      propVal.pop_back();
      res += propVal + "\n";
      propVal = "";
    }
    propVal += apVal + " ";
  }
  if (!propVal.empty()) {
    // remove the trailing space:
    propVal.pop_back();
    res += propVal;
  }
  return res;
}
inline void createAtomIntPropertyList(
    ROMol &mol, const std::string &atomPropName,
    const std::string &missingValueMarker = "", unsigned int lineSize = 190) {
  std::string molPropName = "atom.iprop." + atomPropName;
  mol.setProp(molPropName,
              getAtomPropertyList<boost::int64_t>(
                  mol, atomPropName, missingValueMarker, lineSize));
}
inline void createAtomDoublePropertyList(
    ROMol &mol, const std::string &atomPropName,
    const std::string &missingValueMarker = "", unsigned int lineSize = 190) {
  std::string molPropName = "atom.dprop." + atomPropName;
  mol.setProp(molPropName,
              getAtomPropertyList<double>(mol, atomPropName, missingValueMarker,
                                          lineSize));
}
inline void createAtomBoolPropertyList(
    ROMol &mol, const std::string &atomPropName,
    const std::string &missingValueMarker = "", unsigned int lineSize = 190) {
  std::string molPropName = "atom.bprop." + atomPropName;
  mol.setProp(molPropName,
              getAtomPropertyList<bool>(mol, atomPropName, missingValueMarker,
                                        lineSize));
}
inline void createAtomStringPropertyList(
    ROMol &mol, const std::string &atomPropName,
    const std::string &missingValueMarker = "", unsigned int lineSize = 190) {
  std::string molPropName = "atom.prop." + atomPropName;
  mol.setProp(molPropName,
              getAtomPropertyList<std::string>(mol, atomPropName,
                                               missingValueMarker, lineSize));
}

RDKIT_FILEPARSERS_EXPORT void moveAdditionalPropertiesToSGroups(RWMol &mol);

}  // namespace FileParserUtils
}  // namespace RDKit

#endif
