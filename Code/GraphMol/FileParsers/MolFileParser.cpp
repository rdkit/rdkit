//
//  Copyright (C) 2002-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include "FileParsers.h"
#include "FileParserUtils.h"
#include "MolSGroupParsing.h"
#include "MolFileStereochem.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/Sgroup.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/RDLog.h>

#include <fstream>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <typeinfo>
#include <exception>
#ifdef RDKIT_USE_BOOST_REGEX
#include <boost/regex.hpp>
using boost::regex;
using boost::regex_match;
using boost::smatch;
#else
#include <regex>
using std::regex;
using std::regex_match;
using std::smatch;
#endif
#include <sstream>
#include <locale>
#include <stdlib.h>
#include <cstdio>

using namespace RDKit::SGroupParsing;

namespace RDKit {

namespace FileParserUtils {

int toInt(const std::string &input, bool acceptSpaces) {
  int res = 0;
  // don't need to worry about locale stuff here because
  // we're not going to have delimiters
  res = strtol(input.c_str(), nullptr, 10);
  if (!res && !acceptSpaces && input[0] == ' ') {
    std::string trimmed = boost::trim_copy(input);
    if (trimmed.length() == 0) throw boost::bad_lexical_cast();
  }
  return res;
}

double toDouble(const std::string &input, bool acceptSpaces) {
  double res = atof(input.c_str());
  if (res == 0.0 && !acceptSpaces && input[0] == ' ') {
    std::string trimmed = boost::trim_copy(input);
    if (trimmed.length() == 0) throw boost::bad_lexical_cast();
  }
  return res;
}

std::string getV3000Line(std::istream *inStream, unsigned int &line) {
  // FIX: technically V3K blocks are case-insensitive. We should really be
  // up-casing everything here.
  PRECONDITION(inStream, "bad stream");
  std::string res, tempStr;

  ++line;
  tempStr = getLine(inStream);
  if (tempStr.size() < 7 || tempStr.substr(0, 7) != "M  V30 ") {
    std::ostringstream errout;
    errout << "Line " << line << " does not start with 'M  V30 '" << std::endl;
    throw FileParseException(errout.str());
  }
  // FIX: do we need to handle trailing whitespace after a -?
  while (tempStr[tempStr.length() - 1] == '-') {
    // continuation character, append what we read:
    res += tempStr.substr(7, tempStr.length() - 8);
    // and then read another line:
    ++line;
    tempStr = getLine(inStream);
    if (tempStr.size() < 7 || tempStr.substr(0, 7) != "M  V30 ") {
      std::ostringstream errout;
      errout << "Line " << line << " does not start with 'M  V30 '"
             << std::endl;
      throw FileParseException(errout.str());
    }
  }
  res += tempStr.substr(7, tempStr.length() - 7);

  return res;
}

Atom *replaceAtomWithQueryAtom(RWMol *mol, Atom *atom) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(atom, "bad atom");
  if (atom->hasQuery()) return atom;

  QueryAtom qa(*atom);
  unsigned int idx = atom->getIdx();

  if (atom->getFormalCharge() != 0) {
    qa.expandQuery(makeAtomFormalChargeQuery(atom->getFormalCharge()));
  }
  if (atom->hasProp(common_properties::_hasMassQuery)) {
    qa.expandQuery(makeAtomMassQuery(static_cast<int>(atom->getMass())));
  }
  mol->replaceAtom(idx, &qa);
  return mol->getAtomWithIdx(idx);
}
}  // namespace FileParserUtils
using RDKit::FileParserUtils::getV3000Line;

namespace {

void completeQueryAndChildren(ATOM_EQUALS_QUERY *query, Atom *tgt,
                              int magicVal) {
  PRECONDITION(query, "no query");
  PRECONDITION(tgt, "no atom");
  if (query->getVal() == magicVal) {
    int tgtVal = query->getDataFunc()(tgt);
    query->setVal(tgtVal);
  }
  QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI childIt;
  for (childIt = query->beginChildren(); childIt != query->endChildren();
       ++childIt) {
    completeQueryAndChildren((ATOM_EQUALS_QUERY *)(childIt->get()), tgt,
                             magicVal);
  }
}
void completeMolQueries(RWMol *mol, int magicVal = 0xDEADBEEF) {
  for (ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms();
       ++ai) {
    if ((*ai)->hasQuery()) {
      ATOM_EQUALS_QUERY *query =
          static_cast<ATOM_EQUALS_QUERY *>((*ai)->getQuery());
      completeQueryAndChildren(query, *ai, magicVal);
    }
  }
}

bool startsWith(const std::string &haystack, const char *needle, size_t size) {
  return haystack.compare(0u, size, needle, size) == 0;
}

//! parse a collection block to find enhanced stereo groups
std::string parseEnhancedStereo(std::istream *inStream, unsigned int &line,
                                RWMol *mol) {
  // Lines like (absolute, relative, racemic):
  // M  V30 MDLV30/STEABS ATOMS=(2 2 3)
  // M  V30 MDLV30/STEREL1 ATOMS=(1 12)
  // M  V30 MDLV30/STERAC1 ATOMS=(1 12)
  const regex stereo_label(
      R"regex(MDLV30/STE(...)[0-9]* +ATOMS=\(([0-9]+) +(.*)\))regex");

  smatch match;
  std::vector<StereoGroup> groups;

  // Read the collection until the end
  auto tempStr = getV3000Line(inStream, line);
  boost::to_upper(tempStr);
  while (!startsWith(tempStr, "END", 3)) {
    // If this line in the collection is part of a stereo group
    if (regex_match(tempStr, match, stereo_label)) {
      StereoGroupType grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;

      if (match[1] == "ABS") {
        grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
      } else if (match[1] == "REL") {
        grouptype = RDKit::StereoGroupType::STEREO_OR;
      } else if (match[1] == "RAC") {
        grouptype = RDKit::StereoGroupType::STEREO_AND;
      } else {
        std::ostringstream errout;
        errout << "Unrecognized stereogroup type : '" << tempStr << "' on line"
               << line;
        throw FileParseException(errout.str());
      }

      const unsigned int count = FileParserUtils::toInt(match[2], true);
      std::vector<Atom *> atoms;
      std::stringstream ss(match[3]);
      unsigned int index;
      for (size_t i = 0; i < count; ++i) {
        ss >> index;
        // atoms are 1 indexed in molfiles
        atoms.push_back(mol->getAtomWithIdx(index - 1));
      }
      groups.emplace_back(grouptype, std::move(atoms));
    } else {
      // skip collection types we don't know how to read. Only one documented
      // is MDLV30/HILITE
      BOOST_LOG(rdWarningLog) << "Skipping unrecognized collection type at "
                                 "line "
                              << line << ": " << tempStr << std::endl;
    }
    tempStr = getV3000Line(inStream, line);
  }

  if (!groups.empty()) {
    mol->setStereoGroups(std::move(groups));
  }
  tempStr = getV3000Line(inStream, line);
  return tempStr;
}

//*************************************
//
// Every effort has been made to adhere to MDL's standard
// for mol files
//
//*************************************

void ParseOldAtomList(RWMol *mol, const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  unsigned int idx;
  try {
    idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(0, 3)) -
          1;
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(0, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }

  URANGE_CHECK(idx, mol->getNumAtoms());
  QueryAtom a(*(mol->getAtomWithIdx(idx)));

  auto *q = new ATOM_OR_QUERY;
  q->setDescription("AtomOr");

  switch (text[4]) {
    case 'T':
      q->setNegation(true);
      break;
    case 'F':
      q->setNegation(false);
      break;
    default:
      std::ostringstream errout;
      errout << "Unrecognized atom-list query modifier: '" << text[14]
             << "' on line " << line;
      throw FileParseException(errout.str());
  }

  int nQueries;
  try {
    nQueries = FileParserUtils::toInt(text.substr(9, 1));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(9, 1) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }

  RANGE_CHECK(0, nQueries, 5);
  for (int i = 0; i < nQueries; i++) {
    int pos = 11 + i * 4;
    int atNum;
    try {
      atNum = FileParserUtils::toInt(text.substr(pos, 3));
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(pos, 3) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
    RANGE_CHECK(0, atNum, 200);  // goofy!
    q->addChild(
        QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(atNum)));
    if (!i) a.setAtomicNum(atNum);
  }

  a.setQuery(q);
  a.setProp(common_properties::_MolFileAtomQuery, 1);

  mol->replaceAtom(idx, &a);
}

void ParseChargeLine(RWMol *mol, const std::string &text, bool firstCall,
                     unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  CHG"), "bad charge line");

  // if this line is specified all the atom other than those specified
  // here should carry a charge of 0; but we should only do this once:
  if (firstCall) {
    for (ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms();
         ++ai) {
      (*ai)->setFormalCharge(0);
    }
  }

  int ie, nent;
  try {
    nent = FileParserUtils::toInt(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  int spos = 9;
  for (ie = 0; ie < nent; ie++) {
    int aid, chg;
    try {
      aid = FileParserUtils::toInt(text.substr(spos, 4));
      spos += 4;
      chg = FileParserUtils::toInt(text.substr(spos, 4));
      spos += 4;
      mol->getAtomWithIdx(aid - 1)->setFormalCharge(chg);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParseRadicalLine(RWMol *mol, const std::string &text, bool firstCall,
                      unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  RAD"), "bad charge line");

  // if this line is specified all the atom other than those specified
  // here should carry a charge of 0; but we should only do this once:
  if (firstCall) {
    for (ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms();
         ++ai) {
      (*ai)->setFormalCharge(0);
    }
  }

  int ie, nent;
  try {
    nent = FileParserUtils::toInt(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  int spos = 9;
  for (ie = 0; ie < nent; ie++) {
    int aid, rad;
    std::ostringstream errout;

    try {
      aid = FileParserUtils::toInt(text.substr(spos, 4));
      spos += 4;
      rad = FileParserUtils::toInt(text.substr(spos, 4));
      spos += 4;

      switch (rad) {
        case 1:
          mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(2);
          break;
        case 2:
          mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(1);
          break;
        case 3:
          mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(2);
          break;
        default:
          errout << "Unrecognized radical value " << rad << " for atom "
                 << aid - 1 << " on line " << line << std::endl;
          throw FileParseException(errout.str());
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParsePXALine(RWMol *mol, const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  PXA", "bad PXA line");
  unsigned int pos = 7;
  try {
    unsigned int atIdx =
        FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(pos, 3));
    pos += 3;
    mol->getAtomWithIdx(atIdx - 1)->setProp(
        "_MolFile_PXA", text.substr(pos, text.length() - pos));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(pos, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
}

void ParseIsotopeLine(RWMol *mol, const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  ISO"), "bad isotope line");

  unsigned int nent;
  try {
    nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  unsigned int spos = 9;
  for (unsigned int ie = 0; ie < nent; ie++) {
    unsigned int aid;
    try {
      aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      Atom *atom = mol->getAtomWithIdx(aid - 1);
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        int isotope = FileParserUtils::toInt(text.substr(spos, 4));
        if (isotope < 0) {
          BOOST_LOG(rdErrorLog)
              << " atom " << aid
              << " has a negative isotope value. line:  " << line << std::endl;
        } else {
          atom->setIsotope(isotope);
          spos += 4;
        }
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParseSubstitutionCountLine(RWMol *mol, const std::string &text,
                                unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  SUB"), "bad SUB line");

  unsigned int nent;
  try {
    nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  unsigned int spos = 9;
  for (unsigned int ie = 0; ie < nent; ie++) {
    unsigned int aid;
    int count;
    try {
      aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      Atom *atom = mol->getAtomWithIdx(aid - 1);
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        count = FileParserUtils::toInt(text.substr(spos, 4));
        if (count == 0) continue;
        ATOM_EQUALS_QUERY *q = makeAtomExplicitDegreeQuery(0);
        switch (count) {
          case -1:
            q->setVal(0);
            break;
          case -2:
            q->setVal(atom->getDegree());
            break;
          case 1:
          case 2:
          case 3:
          case 4:
          case 5:
            q->setVal(count);
            break;
          case 6:
            BOOST_LOG(rdWarningLog) << " atom degree query with value 6 found. "
                                       "This will not match degree >6. The MDL "
                                       "spec says it should.  line: "
                                    << line;
            q->setVal(6);
            break;
          default:
            std::ostringstream errout;
            errout << "Value " << count
                   << " is not supported as a degree query. line: " << line;
            throw FileParseException(errout.str());
        }
        if (!atom->hasQuery()) {
          atom = FileParserUtils::replaceAtomWithQueryAtom(mol, atom);
        }
        atom->expandQuery(q, Queries::COMPOSITE_AND);
        spos += 4;
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParseUnsaturationLine(RWMol *mol, const std::string &text,
                           unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  UNS"), "bad UNS line");

  unsigned int nent;
  try {
    nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  unsigned int spos = 9;
  for (unsigned int ie = 0; ie < nent; ie++) {
    unsigned int aid;
    int count;
    try {
      aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      Atom *atom = mol->getAtomWithIdx(aid - 1);
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        count = FileParserUtils::toInt(text.substr(spos, 4));
        if (count == 0) {
          continue;
        } else if (count == 1) {
          ATOM_EQUALS_QUERY *q = makeAtomUnsaturatedQuery();
          if (!atom->hasQuery()) {
            atom = FileParserUtils::replaceAtomWithQueryAtom(mol, atom);
          }
          atom->expandQuery(q, Queries::COMPOSITE_AND);
        } else {
          std::ostringstream errout;
          errout << "Value " << count
                 << " is not supported as an unsaturation "
                    "query (only 0 and 1 are allowed). "
                    "line: "
                 << line;
          throw FileParseException(errout.str());
        }
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParseRingBondCountLine(RWMol *mol, const std::string &text,
                            unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  RBC"), "bad RBC line");

  unsigned int nent;
  try {
    nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  unsigned int spos = 9;
  for (unsigned int ie = 0; ie < nent; ie++) {
    unsigned int aid;
    int count;
    try {
      aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      Atom *atom = mol->getAtomWithIdx(aid - 1);
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        count = FileParserUtils::toInt(text.substr(spos, 4));
        if (count == 0) continue;
        ATOM_EQUALS_QUERY *q = makeAtomRingBondCountQuery(0);
        switch (count) {
          case -1:
            q->setVal(0);
            break;
          case -2:
            q->setVal(0xDEADBEEF);
            mol->setProp(common_properties::_NeedsQueryScan, 1);
            break;
          case 1:
          case 2:
          case 3:
            q->setVal(count);
            break;
          case 4:
            delete q;
            q = static_cast<ATOM_EQUALS_QUERY *>(new ATOM_LESSEQUAL_QUERY);
            q->setVal(4);
            q->setDescription("AtomRingBondCount");
            q->setDataFunc(queryAtomRingBondCount);
            break;
          default:
            std::ostringstream errout;
            errout << "Value " << count
                   << " is not supported as a ring-bond count query. line: "
                   << line;
            throw FileParseException(errout.str());
        }
        if (!atom->hasQuery()) {
          atom = FileParserUtils::replaceAtomWithQueryAtom(mol, atom);
        }
        atom->expandQuery(q, Queries::COMPOSITE_AND);
        spos += 4;
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParseZCHLine(RWMol *mol, const std::string &text, unsigned int line) {
  // part of Alex Clark's ZBO proposal
  // from JCIM 51:3149-57 (2011)
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  ZCH"), "bad ZCH line");

  unsigned int nent;
  try {
    nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  unsigned int spos = 9;
  for (unsigned int ie = 0; ie < nent; ie++) {
    unsigned int aid = 0;
    int val = 0;
    try {
      aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        val = FileParserUtils::stripSpacesAndCast<int>(text.substr(spos, 4));
      }
      if (!aid || aid > mol->getNumAtoms() || aid == 0) {
        std::ostringstream errout;
        errout << "Bad ZCH specification on line " << line;
        throw FileParseException(errout.str());
      }
      spos += 4;
      --aid;
      Atom *atom = mol->getAtomWithIdx(aid);
      if (!atom) {
        std::ostringstream errout;
        errout << "Atom " << aid << " from ZCH specification on line " << line
               << " not found";
        throw FileParseException(errout.str());
      } else {
        atom->setFormalCharge(val);
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParseHYDLine(RWMol *mol, const std::string &text, unsigned int line) {
  // part of Alex Clark's ZBO proposal
  // from JCIM 51:3149-57 (2011)
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  HYD"), "bad HYD line");

  unsigned int nent;
  try {
    nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  unsigned int spos = 9;
  for (unsigned int ie = 0; ie < nent; ie++) {
    unsigned int aid = 0;
    int val = -1;
    try {
      aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        val = FileParserUtils::stripSpacesAndCast<int>(text.substr(spos, 4));
      }
      if (!aid || aid > mol->getNumAtoms() || aid == 0) {
        std::ostringstream errout;
        errout << "Bad HYD specification on line " << line;
        throw FileParseException(errout.str());
      }
      spos += 4;
      --aid;
      Atom *atom = mol->getAtomWithIdx(aid);
      if (!atom) {
        std::ostringstream errout;
        errout << "Atom " << aid << " from HYD specification on line " << line
               << " not found";
        throw FileParseException(errout.str());
      } else {
        if (val >= 0) {
          atom->setProp("_ZBO_H", true);
          atom->setNumExplicitHs(val);
        }
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParseZBOLine(RWMol *mol, const std::string &text, unsigned int line) {
  // part of Alex Clark's ZBO proposal
  // from JCIM 51:3149-57 (2011)
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  ZBO"), "bad ZBO line");

  unsigned int nent;
  try {
    nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  unsigned int spos = 9;
  for (unsigned int ie = 0; ie < nent; ie++) {
    unsigned int bid = 0;
    unsigned int order = 0;
    try {
      bid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        order = FileParserUtils::stripSpacesAndCast<unsigned int>(
            text.substr(spos, 4));
      }
      if (!bid || bid > mol->getNumBonds() || bid == 0) {
        std::ostringstream errout;
        errout << "Bad ZBO specification on line " << line;
        throw FileParseException(errout.str());
      }
      spos += 4;
      --bid;
      Bond *bnd = mol->getBondWithIdx(bid);
      if (!bnd) {
        std::ostringstream errout;
        errout << "Bond " << bid << " from ZBO specification on line " << line
               << " not found";
        throw FileParseException(errout.str());
      } else {
        if (order == 0) {
          bnd->setBondType(Bond::ZERO);
        } else {
          bnd->setBondType(static_cast<Bond::BondType>(order));
        }
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParseMarvinSmartsLine(RWMol *mol, const std::string &text,
                           unsigned int line) {
  const unsigned int atomNumStart = 10;
  const unsigned int smartsStart = 15;
  const unsigned int SMA = 7;
  // M  MRV SMA   1 [*;A]
  // 01234567890123456789
  //           1111111111
  if (text.substr(SMA, 3) != "SMA") {
    return;
  }

  unsigned int idx;
  std::string idxTxt = text.substr(atomNumStart, smartsStart - atomNumStart);
  try {
    idx = FileParserUtils::stripSpacesAndCast<unsigned int>(idxTxt) - 1;
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << idxTxt << "' to an atom index on line "
           << line;
    throw FileParseException(errout.str());
  }

  URANGE_CHECK(idx, mol->getNumAtoms());
  // Should we check the validity of the marvin line here?  Should we
  // automatically
  //   Add these as recursive smarts?  I tend to think so...
  std::string sma = text.substr(smartsStart);
  Atom *at = mol->getAtomWithIdx(idx);
  at->setProp(common_properties::MRV_SMA, sma);
  RWMol *m = 0;
  try {
    m = SmartsToMol(sma);
  } catch (...) {
    // Is this every used?
  }

  if (m) {
    QueryAtom::QUERYATOM_QUERY *query = new RecursiveStructureQuery(m);
    if (!at->hasQuery()) {
      QueryAtom qAt(*at);
      int oidx = at->getIdx();
      mol->replaceAtom(oidx, &qAt);
      at = mol->getAtomWithIdx(oidx);
    }
    at->expandQuery(query, Queries::COMPOSITE_AND);
    at->setProp(common_properties::_MolFileAtomQuery, 1);
  } else {
    std::ostringstream errout;
    errout << "Cannot parse smarts: '" << sma << "' on line " << line;
    throw FileParseException(errout.str());
  }
}

void ParseNewAtomList(RWMol *mol, const std::string &text, unsigned int line) {
  if (text.size() < 15) {
    std::ostringstream errout;
    errout << "Atom list line too short: '" << text << "'";
    throw FileParseException(errout.str());
  }
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  ALS"),
               "bad atom list line");

  unsigned int idx;
  try {
    idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(7, 3)) -
          1;
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(7, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  URANGE_CHECK(idx, mol->getNumAtoms());
  QueryAtom *a = nullptr;

  int nQueries;
  try {
    nQueries = FileParserUtils::toInt(text.substr(10, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(10, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }

  if (!nQueries) {
    BOOST_LOG(rdWarningLog) << "Empty atom list: '" << text << "' on line "
                            << line << "." << std::endl;
    return;
  }

  if (nQueries < 0) {
    std::ostringstream errout;
    errout << "negative length atom list: '" << text << "' on line " << line
           << "." << std::endl;
    throw FileParseException(errout.str());
  }
  for (unsigned int i = 0; i < static_cast<unsigned int>(nQueries); i++) {
    unsigned int pos = 16 + i * 4;
    if (text.size() < pos + 4) {
      std::ostringstream errout;
      errout << "Atom list line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }

    std::string atSymb = text.substr(pos, 4);
    atSymb.erase(atSymb.find(' '), atSymb.size());
    int atNum = PeriodicTable::getTable()->getAtomicNumber(atSymb);
    if (!i) {
      a = new QueryAtom(*(mol->getAtomWithIdx(idx)));
      // replace the query:
      a->setAtomicNum(atNum);
      a->setQuery(makeAtomNumQuery(atNum));
    } else {
      a->expandQuery(makeAtomNumQuery(atNum), Queries::COMPOSITE_OR, true);
    }
  }
  ASSERT_INVARIANT(a, "no atom built");
  a->setProp(common_properties::_MolFileAtomQuery, 1);
  switch (text[14]) {
    case 'T':
      a->getQuery()->setNegation(true);
      break;
    case 'F':
      a->getQuery()->setNegation(false);
      break;
    default:
      std::ostringstream errout;
      errout << "Unrecognized atom-list query modifier: '" << text[14]
             << "' on line " << line;
      delete a;
      throw FileParseException(errout.str());
  }

  mol->replaceAtom(idx, a);
  delete a;
}

void ParseV3000RGroups(RWMol *mol, Atom *&atom, const std::string &text,
                       unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(atom, "bad atom");
  if (text[0] != '(' || text[text.size() - 1] != ')') {
    std::ostringstream errout;
    errout << "Bad RGROUPS specification '" << text << "' on line " << line
           << ". Missing parens.";
    throw FileParseException(errout.str());
  }
  std::vector<std::string> splitToken;
  std::string resid = text.substr(1, text.size() - 2);
  boost::split(splitToken, resid, boost::is_any_of(std::string(" ")));
  if (splitToken.size() < 1) {
    std::ostringstream errout;
    errout << "Bad RGROUPS specification '" << text << "' on line " << line
           << ". Missing values.";
    throw FileParseException(errout.str());
  }
  unsigned int nRs;
  try {
    nRs = FileParserUtils::stripSpacesAndCast<unsigned int>(splitToken[0]);
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << splitToken[0] << "' to int on line" << line;
    throw FileParseException(errout.str());
  }
  if (splitToken.size() < nRs + 1) {
    std::ostringstream errout;
    errout << "Bad RGROUPS specification '" << text << "' on line " << line
           << ". Not enough values.";
    throw FileParseException(errout.str());
  }
  for (unsigned int i = 0; i < nRs; ++i) {
    unsigned int rLabel;
    try {
      rLabel =
          FileParserUtils::stripSpacesAndCast<unsigned int>(splitToken[i + 1]);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << splitToken[i + 1] << "' to int on line"
             << line;
      throw FileParseException(errout.str());
    }
    atom = FileParserUtils::replaceAtomWithQueryAtom(mol, atom);
    atom->setProp(common_properties::_MolFileRLabel, rLabel);
    std::string dLabel = "R" + std::to_string(rLabel);
    atom->setProp(common_properties::dummyLabel, dLabel);
    atom->setIsotope(rLabel);
    atom->setQuery(makeAtomNullQuery());
  }
}

void ParseRGroupLabels(RWMol *mol, const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  RGP"),
               "bad R group label line");

  int nLabels;
  try {
    nLabels = FileParserUtils::toInt(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }

  for (int i = 0; i < nLabels; i++) {
    int pos = 10 + i * 8;
    unsigned int atIdx;
    try {
      atIdx = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(pos, 3));
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(pos, 3) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
    unsigned int rLabel;
    try {
      rLabel = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(pos + 4, 3));
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(pos + 4, 3)
             << "' to int on line " << line;
      throw FileParseException(errout.str());
    }
    atIdx -= 1;
    if (atIdx > mol->getNumAtoms()) {
      std::ostringstream errout;
      errout << "Attempt to set R group label on nonexistent atom " << atIdx
             << " on line " << line;
      throw FileParseException(errout.str());
    }
    QueryAtom qatom(*(mol->getAtomWithIdx(atIdx)));
    qatom.setProp(common_properties::_MolFileRLabel, rLabel);

    // set the dummy label so that this is shown correctly
    // in other pieces of the code :
    // (this was sf.net issue 3316600)
    std::string dLabel = "R" + std::to_string(rLabel);
    qatom.setProp(common_properties::dummyLabel, dLabel);

    // the CTFile spec (June 2005 version) technically only allows
    // R labels up to 32. Since there are three digits, we'll accept
    // anything: so long as it's positive and less than 1000:
    if (rLabel > 0 && rLabel < 999) {
      qatom.setIsotope(rLabel);
    }
    qatom.setQuery(makeAtomNullQuery());
    mol->replaceAtom(atIdx, &qatom);
  }
}

void ParseAtomAlias(RWMol *mol, std::string text, const std::string &nextLine,
                    unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 2) == std::string("A "), "bad atom alias line");

  unsigned int idx;
  try {
    idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(3, 3)) -
          1;
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(3, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  URANGE_CHECK(idx, mol->getNumAtoms());
  Atom *at = mol->getAtomWithIdx(idx);
  at->setProp(common_properties::molFileAlias, nextLine);
}

void ParseAtomValue(RWMol *mol, std::string text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 2) == std::string("V "), "bad atom value line");

  unsigned int idx;
  try {
    idx = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(3, 3)) -
          1;
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(3, 3) << "' to int on line"
           << line;
    throw FileParseException(errout.str());
  }
  URANGE_CHECK(idx, mol->getNumAtoms());
  Atom *at = mol->getAtomWithIdx(idx);
  at->setProp(common_properties::molFileValue,
              text.substr(7, text.length() - 7));
}

Atom *ParseMolFileAtomLine(const std::string text, RDGeom::Point3D &pos,
                           unsigned int line) {
  std::string symb;
  int massDiff, chg, hCount;

  if (text.size() < 34) {
    std::ostringstream errout;
    errout << "Atom line too short: '" << text << "' on line " << line;
    throw FileParseException(errout.str());
  }

  try {
    pos.x = FileParserUtils::toDouble(text.substr(0, 10));
    pos.y = FileParserUtils::toDouble(text.substr(10, 10));
    pos.z = FileParserUtils::toDouble(text.substr(20, 10));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot process coordinates on line " << line;
    throw FileParseException(errout.str());
  }
  symb = text.substr(31, 3);
  boost::trim(symb);

  // REVIEW: should we handle missing fields at the end of the line?
  massDiff = 0;
  if (text.size() >= 36 && text.substr(34, 2) != " 0") {
    try {
      massDiff = FileParserUtils::toInt(text.substr(34, 2), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(34, 2) << "' to into on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
  chg = 0;
  if (text.size() >= 39 && text.substr(36, 3) != "  0") {
    try {
      chg = FileParserUtils::toInt(text.substr(36, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(36, 3) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
  hCount = 0;
  if (text.size() >= 45 && text.substr(42, 3) != "  0") {
    try {
      hCount = FileParserUtils::toInt(text.substr(42, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(42, 3) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
  }
  auto *res = new Atom;
  if (symb == "L" || symb == "A" || symb == "Q" || symb == "*" ||
      symb == "LP" || symb == "R" || symb == "R#" ||
      (symb[0] == 'R' && symb >= "R0" && symb <= "R99")) {
    if (symb == "A" || symb == "Q" || symb == "*") {
      auto *query = new QueryAtom(0);
      if (symb == "*") {
        // according to the MDL spec, these match anything
        query->setQuery(makeAtomNullQuery());
      } else if (symb == "Q") {
        query->setQuery(makeQAtomQuery());
      } else if (symb == "A") {
        query->setQuery(makeAAtomQuery());
      }
      delete res;
      res = query;
      // queries have no implicit Hs:
      res->setNoImplicit(true);
    } else {
      res->setAtomicNum(0);
    }
    if (massDiff == 0 && symb[0] == 'R') {
      if (symb.length() > 1) {
        std::string rlabel = "";
        rlabel = symb.substr(1, symb.length() - 1);
        int rnumber;
        try {
          rnumber = boost::lexical_cast<int>(rlabel);
        } catch (boost::bad_lexical_cast &) {
          rnumber = -1;
        }
        if (rnumber >= 0) res->setIsotope(rnumber);
      }
    }
  } else if (symb == "D") {  // mol blocks support "D" and "T" as shorthand...
                             // handle that.
    res->setAtomicNum(1);
    res->setIsotope(2);
  } else if (symb == "T") {  // mol blocks support "D" and "T" as shorthand...
                             // handle that.
    res->setAtomicNum(1);
    res->setIsotope(3);
  } else {
    if (symb.size() == 2 && symb[1] >= 'A' && symb[1] <= 'Z')
      symb[1] = static_cast<char>(tolower(symb[1]));
    try {
      res->setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(symb));
    } catch (const Invar::Invariant &e) {
      delete res;
      throw FileParseException(e.getMessage());
    }
  }

  // res->setPos(pX,pY,pZ);
  if (chg != 0) res->setFormalCharge(4 - chg);

  // FIX: this does not appear to be correct
  if (hCount == 1) {
    res->setNoImplicit(true);
  }

  if (massDiff != 0) {
    int defIso =
        PeriodicTable::getTable()->getMostCommonIsotope(res->getAtomicNum());
    int dIso = defIso + massDiff;
    if (dIso < 0) {
      BOOST_LOG(rdWarningLog)
          << " atom " << res->getIdx()
          << " has a negative isotope offset. line:  " << line << std::endl;
    }
    res->setIsotope(dIso);
    res->setProp(common_properties::_hasMassQuery, true);
  }

  if (text.size() >= 42 && text.substr(39, 3) != "  0") {
    int parity = 0;
    try {
      parity = FileParserUtils::toInt(text.substr(39, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(39, 3) << "' to int on line "
             << line;
      delete res;
      throw FileParseException(errout.str());
    }
    res->setProp(common_properties::molParity, parity);
  }

  if (text.size() >= 48 && text.substr(45, 3) != "  0") {
    int stereoCare = 0;
    try {
      stereoCare = FileParserUtils::toInt(text.substr(45, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(45, 3) << "' to int on line "
             << line;
      delete res;
      throw FileParseException(errout.str());
    }
    res->setProp("molStereoCare", stereoCare);
  }
  if (text.size() >= 51 && text.substr(48, 3) != "  0") {
    int totValence = 0;
    try {
      totValence = FileParserUtils::toInt(text.substr(48, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(48, 3) << "' to int on line "
             << line;
      delete res;
      throw FileParseException(errout.str());
    }
    if (totValence != 0) {
      // only set if it's a non-default value
      res->setProp(common_properties::molTotValence, totValence);
    }
  }
  if (text.size() >= 57 && text.substr(54, 3) != "  0") {
    int rxnRole = 0;
    try {
      rxnRole = FileParserUtils::toInt(text.substr(54, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(54, 3) << "' to int on line "
             << line;
      delete res;
      throw FileParseException(errout.str());
    }
    if (rxnRole != 0) {
      // only set if it's a non-default value
      res->setProp(common_properties::molRxnRole, rxnRole);
    }
  }
  if (text.size() >= 60 && text.substr(57, 3) != "  0") {
    int rxnComponent = 0;
    try {
      rxnComponent = FileParserUtils::toInt(text.substr(57, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(57, 3) << "' to int on line "
             << line;
      delete res;
      throw FileParseException(errout.str());
    }
    if (rxnComponent != 0) {
      // only set if it's a non-default value
      res->setProp(common_properties::molRxnComponent, rxnComponent);
    }
  }
  if (text.size() >= 63 && text.substr(60, 3) != "  0") {
    int atomMapNumber = 0;
    try {
      atomMapNumber = FileParserUtils::toInt(text.substr(60, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(60, 3) << "' to int on line "
             << line;
      delete res;
      throw FileParseException(errout.str());
    }
    res->setProp(common_properties::molAtomMapNumber, atomMapNumber);
  }
  if (text.size() >= 66 && text.substr(63, 3) != "  0") {
    int inversionFlag = 0;
    try {
      inversionFlag = FileParserUtils::toInt(text.substr(63, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(63, 3) << "' to int on line "
             << line;
      delete res;
      throw FileParseException(errout.str());
    }
    res->setProp(common_properties::molInversionFlag, inversionFlag);
  }
  if (text.size() >= 69 && text.substr(66, 3) != "  0") {
    int exactChangeFlag = 0;
    try {
      exactChangeFlag = FileParserUtils::toInt(text.substr(66, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(66, 3) << "' to int on line "
             << line;
      delete res;
      throw FileParseException(errout.str());
    }
    res->setProp("molExactChangeFlag", exactChangeFlag);
  }
  return res;
}

Bond *ParseMolFileBondLine(const std::string &text, unsigned int line) {
  unsigned int idx1, idx2, bType, stereo;
  int spos = 0;

  if (text.size() < 9) {
    std::ostringstream errout;
    errout << "Bond line too short: '" << text << "' on line " << line;
    throw FileParseException(errout.str());
  }

  try {
    idx1 = FileParserUtils::toInt(text.substr(spos, 3));
    spos += 3;
    idx2 = FileParserUtils::toInt(text.substr(spos, 3));
    spos += 3;
    bType = FileParserUtils::toInt(text.substr(spos, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(spos, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }

  // adjust the numbering
  idx1--;
  idx2--;

  Bond::BondType type;
  Bond *res = nullptr;
  switch (bType) {
    case 1:
      type = Bond::SINGLE;
      res = new Bond;
      break;
    case 2:
      type = Bond::DOUBLE;
      res = new Bond;
      break;
    case 3:
      type = Bond::TRIPLE;
      res = new Bond;
      break;
    case 4:
      type = Bond::AROMATIC;
      res = new Bond;
      break;
    case 0:
      type = Bond::UNSPECIFIED;
      res = new Bond;
      BOOST_LOG(rdWarningLog)
          << "bond with order 0 found on line " << line
          << ". This is not part of the MDL specification." << std::endl;
      break;
    default:
      type = Bond::UNSPECIFIED;
      // it's a query bond of some type
      res = new QueryBond;
      if (bType == 8) {
        BOND_NULL_QUERY *q;
        q = makeBondNullQuery();
        res->setQuery(q);
      } else if (bType == 6) {
        res->setQuery(makeSingleOrAromaticBondQuery());
        res->setProp(common_properties::_MolFileBondQuery, 1);
      } else if (bType == 5 || bType == 7) {
        BOND_OR_QUERY *q;
        q = new BOND_OR_QUERY;
        if (bType == 5) {
          // single or double
          q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
              makeBondOrderEqualsQuery(Bond::SINGLE)));
          q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
              makeBondOrderEqualsQuery(Bond::DOUBLE)));
          q->setDescription("BondOr");
          res->setProp(common_properties::_MolFileBondQuery, 1);
        } else if (bType == 7) {
          // double or aromatic
          q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
              makeBondOrderEqualsQuery(Bond::DOUBLE)));
          q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
              makeBondOrderEqualsQuery(Bond::AROMATIC)));
          q->setDescription("BondOr");
        }
        res->setQuery(q);
        res->setProp(common_properties::_MolFileBondQuery, 1);
      } else {
        BOND_NULL_QUERY *q;
        q = makeBondNullQuery();
        res->setQuery(q);
        BOOST_LOG(rdWarningLog)
            << "unrecognized query bond type, " << bType << ", found on line "
            << line << ". Using an \"any\" query." << std::endl;
      }
      break;
  }
  res->setBeginAtomIdx(idx1);
  res->setEndAtomIdx(idx2);
  res->setBondType(type);
  res->setProp(common_properties::_MolFileBondType, bType);

  if (text.size() >= 12 && text.substr(9, 3) != "  0") {
    try {
      stereo = FileParserUtils::toInt(text.substr(9, 3));
      switch (stereo) {
        case 0:
          res->setBondDir(Bond::NONE);
          break;
        case 1:
          res->setBondDir(Bond::BEGINWEDGE);
          break;
        case 6:
          res->setBondDir(Bond::BEGINDASH);
          break;
        case 3:  // "either" double bond
          res->setBondDir(Bond::EITHERDOUBLE);
          res->setStereo(Bond::STEREOANY);
          break;
        case 4:  // "either" single bond
          res->setBondDir(Bond::UNKNOWN);
          break;
      }
      res->setProp(common_properties::_MolFileBondStereo, stereo);
    } catch (boost::bad_lexical_cast &) {
      ;
    }
  }
  if (text.size() >= 18 && text.substr(15, 3) != "  0") {
    try {
      int topology = FileParserUtils::toInt(text.substr(15, 3));
      if (topology) {
        if (!res->hasQuery()) {
          auto *qBond = new QueryBond(*res);
          delete res;
          res = qBond;
        }
        BOND_EQUALS_QUERY *q = makeBondIsInRingQuery();
        switch (topology) {
          case 1:
            break;
          case 2:
            q->setNegation(true);
            break;
          default:
            std::ostringstream errout;
            errout << "Unrecognized bond topology specifier: " << topology
                   << " on line " << line;
            throw FileParseException(errout.str());
        }
        res->expandQuery(q);
      }
    } catch (boost::bad_lexical_cast &) {
      ;
    }
  }
  if (text.size() >= 21 && text.substr(18, 3) != "  0") {
    try {
      int reactStatus = FileParserUtils::toInt(text.substr(18, 3));
      res->setProp("molReactStatus", reactStatus);
    } catch (boost::bad_lexical_cast &) {
      ;
    }
  }
  return res;
}

void ParseMolBlockAtoms(std::istream *inStream, unsigned int &line,
                        unsigned int nAtoms, RWMol *mol, Conformer *conf) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(conf, "bad conformer");
  for (unsigned int i = 1; i <= nAtoms; ++i) {
    ++line;
    std::string tempStr = getLine(inStream);
    if (inStream->eof()) {
      throw FileParseException("EOF hit while reading atoms");
    }
    RDGeom::Point3D pos;
    Atom *atom = ParseMolFileAtomLine(tempStr, pos, line);
    unsigned int aid = mol->addAtom(atom, false, true);
    conf->setAtomPos(aid, pos);
    mol->setAtomBookmark(atom, i);
  }
}

// returns whether or not any sign of chirality was detected
void ParseMolBlockBonds(std::istream *inStream, unsigned int &line,
                        unsigned int nBonds, RWMol *mol,
                        bool &chiralityPossible) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(mol, "bad molecule");
  for (unsigned int i = 1; i <= nBonds; ++i) {
    ++line;
    std::string tempStr = getLine(inStream);
    if (inStream->eof()) {
      throw FileParseException("EOF hit while reading bonds");
    }
    Bond *bond = ParseMolFileBondLine(tempStr, line);
    // if we got an aromatic bond set the flag on the bond and the connected
    // atoms
    if (bond->getBondType() == Bond::AROMATIC) {
      bond->setIsAromatic(true);
      mol->getAtomWithIdx(bond->getBeginAtomIdx())->setIsAromatic(true);
      mol->getAtomWithIdx(bond->getEndAtomIdx())->setIsAromatic(true);
    }
    // if the bond might have chirality info associated with it, set a flag:
    if (bond->getBondDir() != Bond::NONE &&
        bond->getBondDir() != Bond::UNKNOWN) {
      chiralityPossible = true;
    }
    mol->addBond(bond, true);
    mol->setBondBookmark(bond, i);
  }
}

bool ParseMolBlockProperties(std::istream *inStream, unsigned int &line,
                             RWMol *mol, bool strictParsing) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(mol, "bad molecule");
  // older mol files can have an atom list block here
  std::string tempStr = getLine(inStream);
  ++line;
  // there is apparently some software out there that puts a
  // blank line in mol blocks before the "M  END". If we aren't
  // doing strict parsing, deal with that here.
  if (!tempStr.size()) {
    if (!strictParsing) {
      tempStr = getLine(inStream);
      ++line;
    } else {
      std::ostringstream errout;
      errout << "Problems encountered parsing Mol data, unexpected blank line "
                "found at line "
             << line;
      throw FileParseException(errout.str());
    }
  } else {
    if (tempStr[0] != 'M' && tempStr[0] != 'A' && tempStr[0] != 'V' &&
        tempStr[0] != 'G' && tempStr[0] != 'S') {
      ParseOldAtomList(mol, tempStr, line);
    }
  }

  IDX_TO_SGROUP_MAP sGroupMap;
  IDX_TO_STR_VECT_MAP dataFieldsMap;
  bool fileComplete = false;
  bool firstChargeLine = true;
  unsigned int SCDcounter = 0;
  unsigned int lastDataSGroup = 0;
  std::ostringstream currentDataField;
  std::string lineBeg = tempStr.substr(0, 6);
  while (!inStream->eof() && lineBeg != "M  END" &&
         tempStr.substr(0, 4) != "$$$$") {
    if (tempStr[0] == 'A') {
      line++;
      std::string nextLine = getLine(inStream);
      if (lineBeg != "M  END") {
        ParseAtomAlias(mol, tempStr, nextLine, line);
      }
    } else if (tempStr[0] == 'G') {
      BOOST_LOG(rdWarningLog)
          << " deprecated group abbreviation ignored on line " << line
          << std::endl;
      // we need to skip the next line, which holds the abbreviation:
      line++;
      tempStr = getLine(inStream);
    } else if (tempStr[0] == 'V') {
      ParseAtomValue(mol, tempStr, line);
    } else if (lineBeg == "S  SKP") {
      int nToSkip = FileParserUtils::toInt(tempStr.substr(6, 3));
      if (nToSkip < 0) {
        std::ostringstream errout;
        errout << "negative skip value " << nToSkip << " on line " << line;
        throw FileParseException(errout.str());
      }
      for (unsigned int i = 0; i < static_cast<unsigned int>(nToSkip); ++i) {
        ++line;
        tempStr = getLine(inStream);
      }
    } else if (lineBeg == "M  ALS")
      ParseNewAtomList(mol, tempStr, line);
    else if (lineBeg == "M  ISO")
      ParseIsotopeLine(mol, tempStr, line);
    else if (lineBeg == "M  RGP")
      ParseRGroupLabels(mol, tempStr, line);
    else if (lineBeg == "M  RBC")
      ParseRingBondCountLine(mol, tempStr, line);
    else if (lineBeg == "M  SUB")
      ParseSubstitutionCountLine(mol, tempStr, line);
    else if (lineBeg == "M  UNS")
      ParseUnsaturationLine(mol, tempStr, line);
    else if (lineBeg == "M  CHG") {
      ParseChargeLine(mol, tempStr, firstChargeLine, line);
      firstChargeLine = false;
    } else if (lineBeg == "M  RAD") {
      ParseRadicalLine(mol, tempStr, firstChargeLine, line);
      firstChargeLine = false;
    } else if (lineBeg == "M  PXA") {
      ParsePXALine(mol, tempStr, line);

      /* SGroup parsing start */
    } else if (lineBeg == "M  STY") {
      ParseSGroupV2000STYLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SST") {
      ParseSGroupV2000SSTLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SLB") {
      ParseSGroupV2000SLBLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SCN") {
      ParseSGroupV2000SCNLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SDS") {
      ParseSGroupV2000SDSLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SAL" || lineBeg == "M  SBL" ||
               lineBeg == "M  SPA") {
      ParseSGroupV2000VectorDataLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SMT") {
      ParseSGroupV2000SMTLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SDI") {
      ParseSGroupV2000SDILine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  CRS") {
      std::ostringstream errout;
      errout << "Unsupported SGroup subtype '" << lineBeg << "' on line "
             << line;
      throw FileParseException(errout.str());
    } else if (lineBeg == "M  SBV") {
      ParseSGroupV2000SBVLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SDT") {
      ParseSGroupV2000SDTLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SDD") {
      ParseSGroupV2000SDDLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SCD" || lineBeg == "M  SED") {
      ParseSGroupV2000SCDSEDLine(sGroupMap, dataFieldsMap, mol, tempStr, line,
                                 strictParsing, SCDcounter, lastDataSGroup,
                                 currentDataField);
    } else if (lineBeg == "M  SPL") {
      ParseSGroupV2000SPLLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SNC") {
      ParseSGroupV2000SNCLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SAP") {
      ParseSGroupV2000SAPLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SCL") {
      ParseSGroupV2000SCLLine(sGroupMap, mol, tempStr, line);
    } else if (lineBeg == "M  SBT") {
      ParseSGroupV2000SBTLine(sGroupMap, mol, tempStr, line);

      /* SGroup parsing end */
    } else if (lineBeg == "M  ZBO")
      ParseZBOLine(mol, tempStr, line);
    else if (lineBeg == "M  ZCH") {
      ParseZCHLine(mol, tempStr, line);
    } else if (lineBeg == "M  HYD") {
      ParseHYDLine(mol, tempStr, line);
    } else if (lineBeg == "M  MRV") {
      ParseMarvinSmartsLine(mol, tempStr, line);
    }
    line++;
    tempStr = getLine(inStream);
    lineBeg = tempStr.substr(0, 6);
  }
  if (tempStr[0] == 'M' && tempStr.substr(0, 6) == "M  END") {
    // All went well, make final updates to SGroups, and add them to Mol
    for (const auto &sgroup : sGroupMap) {
      sgroup.second.setProp("DATAFIELDS", dataFieldsMap[sgroup.first]);
      addSGroup(*mol, sgroup.second);
    }

    fileComplete = true;
  }
  return fileComplete;
}

Atom *ParseV3000AtomSymbol(std::string token, unsigned int &line) {
  bool negate = false;
  boost::trim(token);
  std::string cpy = token;
  boost::to_upper(cpy);
  if (cpy.size() > 3 && cpy.substr(0, 3) == "NOT") {
    negate = true;
    token = token.substr(3, token.size() - 3);
    boost::trim(token);
  }

  Atom *res = nullptr;
  if (token[0] == '[') {
    // atom list:
    if (token[token.length() - 1] != ']') {
      std::ostringstream errout;
      errout << "Bad atom token '" << token << "' on line: " << line;
      throw FileParseException(errout.str());
    }
    token = token.substr(1, token.size() - 2);

    std::vector<std::string> splitToken;
    boost::split(splitToken, token, boost::is_any_of(","));

    for (std::vector<std::string>::const_iterator stIt = splitToken.begin();
         stIt != splitToken.end(); ++stIt) {
      std::string atSymb = boost::trim_copy(*stIt);
      if (atSymb == "") continue;
      if (atSymb.size() == 2 && atSymb[1] >= 'A' && atSymb[1] <= 'Z')
        atSymb[1] = static_cast<char>(tolower(atSymb[1]));

      int atNum = PeriodicTable::getTable()->getAtomicNumber(atSymb);
      if (!res) {
        res = new QueryAtom(atNum);
      } else {
        res->expandQuery(makeAtomNumQuery(atNum), Queries::COMPOSITE_OR, true);
      }
    }
    res->getQuery()->setNegation(negate);
  } else {
    if (negate) {
      std::ostringstream errout;
      errout << "NOT tokens only supported for atom lists. line " << line;
      throw FileParseException(errout.str());
    }
    // it's a normal CTAB atom symbol:
    // NOTE: "R" and "R0"-"R99" are not in the v3K CTAB spec, but we're going to
    // support them anyway
    if (token == "R" || (token[0] == 'R' && token >= "R0" && token <= "R99") ||
        token == "R#" || token == "A" || token == "Q" || token == "*") {
      if (token == "A" || token == "Q" || token == "*") {
        res = new QueryAtom(0);
        if (token == "*") {
          // according to the MDL spec, these match anything
          res->setQuery(makeAtomNullQuery());
        } else if (token == "Q") {
          auto *q = new ATOM_OR_QUERY;
          q->setDescription("AtomOr");
          q->setNegation(true);
          q->addChild(
              QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(6)));
          q->addChild(
              QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(1)));
          res->setQuery(q);
        } else if (token == "A") {
          res->setQuery(makeAtomNumQuery(1));
          res->getQuery()->setNegation(true);
        }
        // queries have no implicit Hs:
        res->setNoImplicit(true);
      } else {
        res = new Atom(1);
        res->setAtomicNum(0);
      }
      if (token[0] == 'R' && token >= "R0" && token <= "R99") {
        std::string rlabel = "";
        rlabel = token.substr(1, token.length() - 1);
        int rnumber;
        try {
          rnumber = boost::lexical_cast<int>(rlabel);
        } catch (boost::bad_lexical_cast &) {
          rnumber = -1;
        }
        if (rnumber >= 0) res->setIsotope(rnumber);
      }
    } else if (token == "D") {  // mol blocks support "D" and "T" as
                                // shorthand... handle that.
      res = new Atom(1);
      res->setIsotope(2);
    } else if (token == "T") {  // mol blocks support "D" and "T" as
                                // shorthand... handle that.
      res = new Atom(1);
      res->setIsotope(3);
    } else {
      if (token.size() == 2 && token[1] >= 'A' && token[1] <= 'Z')
        token[1] = static_cast<char>(tolower(token[1]));

      res = new Atom(PeriodicTable::getTable()->getAtomicNumber(token));
    }
  }

  POSTCONDITION(res, "no atom built");
  return res;
}

bool splitAssignToken(const std::string &token, std::string &prop,
                      std::string &val) {
  std::vector<std::string> splitToken;
  boost::split(splitToken, token, boost::is_any_of("="));
  if (splitToken.size() != 2) {
    return false;
  }
  prop = splitToken[0];
  boost::to_upper(prop);
  val = splitToken[1];
  return true;
}

template <class T>
void ParseV3000AtomProps(RWMol *mol, Atom *&atom, typename T::iterator &token,
                         const T &tokens, unsigned int &line) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(atom, "bad atom");
  std::ostringstream errout;
  while (token != tokens.end()) {
    std::string prop, val;
    if (!splitAssignToken(*token, prop, val)) {
      errout << "Invalid atom property: '" << *token << "' for atom "
             << atom->getIdx() + 1 << " on line " << line << std::endl;
      throw FileParseException(errout.str());
    }

    if (prop == "CHG") {
      int charge = FileParserUtils::toInt(val);
      if (!atom->hasQuery()) {
        atom->setFormalCharge(charge);
      } else {
        atom->expandQuery(makeAtomFormalChargeQuery(charge));
      }
    } else if (prop == "RAD") {
      // FIX handle queries here
      switch (FileParserUtils::toInt(val)) {
        case 0:
          break;
        case 1:
          atom->setNumRadicalElectrons(2);
          break;
        case 2:
          atom->setNumRadicalElectrons(1);
          break;
        case 3:
          atom->setNumRadicalElectrons(2);
          break;
        default:
          errout << "Unrecognized RAD value " << val << " for atom "
                 << atom->getIdx() + 1 << " on line " << line << std::endl;
          throw FileParseException(errout.str());
      }
    } else if (prop == "MASS") {
      // the documentation for V3000 CTABs says that this should contain the
      // "absolute atomic weight" (whatever that means).
      // Online examples seem to have integer (isotope) values and Marvin won't
      // even read something that has a float.
      // We'll go with the int
      int v;
      double dv;
      try {
        v = FileParserUtils::toInt(val);
      } catch (boost::bad_lexical_cast &) {
        try {
          dv = FileParserUtils::toDouble(val);
          v = static_cast<int>(floor(dv));
        } catch (boost::bad_lexical_cast &) {
          v = -1;
        }
      }
      if (v < 0) {
        errout << "Bad value for MASS :" << val << " for atom "
               << atom->getIdx() + 1 << " on line " << line << std::endl;
        throw FileParseException(errout.str());
      } else {
        if (!atom->hasQuery()) {
          atom->setIsotope(v);
        } else {
          atom->expandQuery(makeAtomIsotopeQuery(v));
        }
      }
    } else if (prop == "CFG") {
      int cfg = FileParserUtils::toInt(val);
      switch (cfg) {
        case 0:
          break;
        case 1:
        case 2:
        case 3:
          atom->setProp(common_properties::molParity, cfg);
          break;
        default:
          errout << "Unrecognized CFG value : " << val << " for atom "
                 << atom->getIdx() + 1 << " on line " << line << std::endl;
          throw FileParseException(errout.str());
      }
    } else if (prop == "HCOUNT") {
      if (val != "0") {
        int hcount = FileParserUtils::toInt(val);
        if (!atom->hasQuery()) {
          atom = FileParserUtils::replaceAtomWithQueryAtom(mol, atom);
        }
        if (hcount == -1) hcount = 0;
        atom->expandQuery(makeAtomHCountQuery(hcount));
      }
    } else if (prop == "UNSAT") {
      if (val == "1") {
        if (!atom->hasQuery()) {
          atom = FileParserUtils::replaceAtomWithQueryAtom(mol, atom);
        }
        atom->expandQuery(makeAtomUnsaturatedQuery());
      }
    } else if (prop == "RBCNT") {
      if (val != "0") {
        int rbcount = FileParserUtils::toInt(val);
        if (!atom->hasQuery()) {
          atom = FileParserUtils::replaceAtomWithQueryAtom(mol, atom);
        }
        if (rbcount == -1) rbcount = 0;
        atom->expandQuery(makeAtomRingBondCountQuery(rbcount));
      }
    } else if (prop == "VAL") {
      if (val != "0") {
        int totval = FileParserUtils::toInt(val);
        atom->setProp(common_properties::molTotValence, totval);
      }
    } else if (prop == "RGROUPS") {
      ParseV3000RGroups(mol, atom, val, line);
      // FIX
    }
    ++token;
  }
}

void tokenizeV3000Line(std::string line, std::vector<std::string> &tokens) {
  tokens.clear();
  bool inQuotes = false, inParens = false;
  unsigned int start = 0;
  unsigned int pos = 0;
  while (pos < line.size()) {
    if (line[pos] == ' ' || line[pos] == '\t') {
      if (start == pos) {
        ++start;
        ++pos;
      } else if (!inQuotes && !inParens) {
        tokens.push_back(line.substr(start, pos - start));
        ++pos;
        start = pos;
      } else {
        ++pos;
      }
    } else if (line[pos] == ')' && inParens) {
      tokens.push_back(line.substr(start, pos - start + 1));
      inParens = false;
      ++pos;
      start = pos;
    } else if (line[pos] == '(' && !inQuotes) {
      inParens = true;
      ++pos;
    } else if (line[pos] == '"' && !inParens) {
      if (pos + 1 < line.size() && line[pos + 1] == '"') {
        pos += 2;
      } else if (inQuotes) {
        // don't push on the quotes themselves
        tokens.push_back(line.substr(start + 1, pos - start - 1));
        ++pos;
        start = pos;
        inQuotes = false;
      } else {
        ++pos;
        inQuotes = true;
      }
    } else {
      ++pos;
    }
  }
  if (start != pos) {
    tokens.push_back(line.substr(start, line.size() - start));
  }
#if 0
      std::cerr<<"tokens: ";
      std::copy(tokens.begin(),tokens.end(),std::ostream_iterator<std::string>(std::cerr,"|"));
      std::cerr<<std::endl;
#endif
}

void ParseV3000AtomBlock(std::istream *inStream, unsigned int &line,
                         unsigned int nAtoms, RWMol *mol, Conformer *conf) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(nAtoms > 0, "bad atom count");
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(conf, "bad conformer");
  std::string tempStr;
  std::vector<std::string> splitLine;

  tempStr = getV3000Line(inStream, line);
  if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN ATOM") {
    std::ostringstream errout;
    errout << "BEGIN ATOM line not found on line " << line;
    throw FileParseException(errout.str());
  }
  for (unsigned int i = 0; i < nAtoms; ++i) {
    tempStr = getV3000Line(inStream, line);
    std::string trimmed = boost::trim_copy(tempStr);

    std::vector<std::string> tokens;
    std::vector<std::string>::iterator token;

    tokenizeV3000Line(trimmed, tokens);
    token = tokens.begin();

    if (token == tokens.end()) {
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line" << line;
      throw FileParseException(errout.str());
    }
    unsigned int molIdx = atoi(token->c_str());

    // start with the symbol:
    ++token;
    if (token == tokens.end()) {
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line " << line;
      throw FileParseException(errout.str());
    }
    Atom *atom = ParseV3000AtomSymbol(*token, line);

    // now the position;
    RDGeom::Point3D pos;
    ++token;
    if (token == tokens.end()) {
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line " << line;
      throw FileParseException(errout.str());
    }
    pos.x = atof(token->c_str());
    ++token;
    if (token == tokens.end()) {
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line " << line;
      throw FileParseException(errout.str());
    }
    pos.y = atof(token->c_str());
    ++token;
    if (token == tokens.end()) {
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line " << line;
      throw FileParseException(errout.str());
    }
    pos.z = atof(token->c_str());
    // the map number:
    ++token;
    if (token == tokens.end()) {
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line " << line;
      throw FileParseException(errout.str());
    }
    int mapNum = atoi(token->c_str());
    if (mapNum > 0) {
      atom->setProp(common_properties::molAtomMapNumber, mapNum);
    }
    ++token;

    unsigned int aid = mol->addAtom(atom, false, true);

    // additional properties this may change the atom,
    // so be careful with it:
    ParseV3000AtomProps(mol, atom, token, tokens, line);

    mol->setAtomBookmark(atom, molIdx);
    conf->setAtomPos(aid, pos);
  }
  tempStr = getV3000Line(inStream, line);
  if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END ATOM") {
    std::ostringstream errout;
    errout << "END ATOM line not found on line " << line;
    throw FileParseException(errout.str());
  }

  bool nonzeroZ = hasNonZeroZCoords(*conf);
  if (mol->hasProp(common_properties::_3DConf)) {
    conf->set3D(true);
    mol->clearProp(common_properties::_3DConf);
    if (!nonzeroZ) {
      BOOST_LOG(rdWarningLog)
          << "Warning: molecule is tagged as 3D, but all Z coords are zero"
          << std::endl;
    }
  } else {
    conf->set3D(nonzeroZ);
  }
}
void ParseV3000BondBlock(std::istream *inStream, unsigned int &line,
                         unsigned int nBonds, RWMol *mol,
                         bool &chiralityPossible) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(nBonds > 0, "bad bond count");
  PRECONDITION(mol, "bad molecule");

  auto tempStr = getV3000Line(inStream, line);
  if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN BOND") {
    throw FileParseException("BEGIN BOND line not found");
  }
  for (unsigned int i = 0; i < nBonds; ++i) {
    tempStr = boost::trim_copy(getV3000Line(inStream, line));
    std::vector<std::string> splitLine;
    tokenizeV3000Line(tempStr, splitLine);
    if (splitLine.size() < 4) {
      std::ostringstream errout;
      errout << "bond line " << line << " is too short";
      throw FileParseException(errout.str());
    }
    Bond *bond;
    unsigned int bondIdx = atoi(splitLine[0].c_str());
    unsigned int bType = atoi(splitLine[1].c_str());
    unsigned int a1Idx = atoi(splitLine[2].c_str());
    unsigned int a2Idx = atoi(splitLine[3].c_str());

    switch (bType) {
      case 1:
        bond = new Bond(Bond::SINGLE);
        break;
      case 2:
        bond = new Bond(Bond::DOUBLE);
        break;
      case 3:
        bond = new Bond(Bond::TRIPLE);
        break;
      case 4:
        bond = new Bond(Bond::AROMATIC);
        bond->setIsAromatic(true);
        break;
      case 9:
        bond = new Bond(Bond::DATIVE);
        break;
      case 0:
        bond = new Bond(Bond::UNSPECIFIED);
        BOOST_LOG(rdWarningLog)
            << "bond with order 0 found on line " << line
            << ". This is not part of the MDL specification." << std::endl;
        break;
      default:
        // it's a query bond of some type
        bond = new QueryBond;
        if (bType == 8) {
          BOND_NULL_QUERY *q;
          q = makeBondNullQuery();
          bond->setQuery(q);
        } else if (bType == 6) {
          bond->setQuery(makeSingleOrAromaticBondQuery());
          bond->setProp(common_properties::_MolFileBondQuery, 1);
        } else if (bType == 5 || bType == 7) {
          BOND_OR_QUERY *q;
          q = new BOND_OR_QUERY;
          if (bType == 5) {
            // single or double
            q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
                makeBondOrderEqualsQuery(Bond::SINGLE)));
            q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
                makeBondOrderEqualsQuery(Bond::DOUBLE)));
            q->setDescription("BondOr");
          } else if (bType == 7) {
            // double or aromatic
            q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
                makeBondOrderEqualsQuery(Bond::DOUBLE)));
            q->addChild(QueryBond::QUERYBOND_QUERY::CHILD_TYPE(
                makeBondOrderEqualsQuery(Bond::AROMATIC)));
            q->setDescription("BondOr");
          }
          bond->setProp(common_properties::_MolFileBondQuery, 1);
          bond->setQuery(q);
        } else {
          BOND_NULL_QUERY *q;
          q = makeBondNullQuery();
          bond->setQuery(q);
          BOOST_LOG(rdWarningLog)
              << "unrecognized query bond type, " << bType << ", found on line "
              << line << ". Using an \"any\" query." << std::endl;
        }
        break;
    }
    bond->setProp(common_properties::_MolFileBondType, bType);

    // additional bond properties:
    unsigned int lPos = 4;
    std::ostringstream errout;
    while (lPos < splitLine.size()) {
      std::string prop, val;
      if (!splitAssignToken(splitLine[lPos], prop, val)) {
        errout << "bad bond property '" << splitLine[lPos] << "' on line "
               << line;
        throw FileParseException(errout.str());
      }
      if (prop == "CFG") {
        unsigned int cfg = atoi(val.c_str());
        switch (cfg) {
          case 0:
            break;
          case 1:
            bond->setBondDir(Bond::BEGINWEDGE);
            chiralityPossible = true;
            break;
          case 2:
            if (bType == 1)
              bond->setBondDir(Bond::UNKNOWN);
            else if (bType == 2) {
              bond->setBondDir(Bond::EITHERDOUBLE);
              bond->setStereo(Bond::STEREOANY);
            }
            break;
          case 3:
            bond->setBondDir(Bond::BEGINDASH);
            chiralityPossible = true;
            break;
          default:
            errout << "bad bond CFG " << val << "' on line " << line;
            throw FileParseException(errout.str());
        }
        bond->setProp(common_properties::_MolFileBondCfg, cfg);
      } else if (prop == "TOPO") {
        if (val != "0") {
          if (!bond->hasQuery()) {
            auto *qBond = new QueryBond(*bond);
            delete bond;
            bond = qBond;
          }
          BOND_EQUALS_QUERY *q = makeBondIsInRingQuery();
          if (val == "1") {
            // nothing
          } else if (val == "2") {
            q->setNegation(true);
          } else {
            errout << "bad bond TOPO " << val << "' on line " << line;
            throw FileParseException(errout.str());
          }
          bond->expandQuery(q);
        }
      } else if (prop == "RXCTR") {
        int reactStatus = FileParserUtils::toInt(val);
        bond->setProp("molReactStatus", reactStatus);
      } else if (prop == "STBOX") {
      } else if (prop == "ENDPTS") {
        bond->setProp(common_properties::_MolFileBondEndPts, val);
      } else if (prop == "ATTACH") {
        bond->setProp(common_properties::_MolFileBondAttach, val);
      }
      ++lPos;
    }

    bond->setBeginAtomIdx(mol->getAtomWithBookmark(a1Idx)->getIdx());
    bond->setEndAtomIdx(mol->getAtomWithBookmark(a2Idx)->getIdx());
    mol->addBond(bond, true);
    if (bond->getIsAromatic()) {
      mol->getAtomWithIdx(bond->getBeginAtomIdx())->setIsAromatic(true);
      mol->getAtomWithIdx(bond->getEndAtomIdx())->setIsAromatic(true);
    }
    mol->setBondBookmark(bond, bondIdx);
  }
  tempStr = getV3000Line(inStream, line);
  if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END BOND") {
    std::ostringstream errout;
    errout << "END BOND line not found at line " << line;
    throw FileParseException(errout.str());
  }
}

void ProcessMolProps(RWMol *mol) {
  PRECONDITION(mol, "no molecule");
  for (RWMol::AtomIterator atomIt = mol->beginAtoms();
       atomIt != mol->endAtoms(); ++atomIt) {
    Atom *atom = *atomIt;
    int totV;
    if (atom->getPropIfPresent(common_properties::molTotValence, totV) &&
        !atom->hasProp("_ZBO_H")) {
      if (totV == 0) continue;
      atom->setNoImplicit(true);
      if (totV == 15     // V2000
          || totV == -1  // v3000
      ) {
        atom->setNumExplicitHs(0);
      } else {
        if (atom->getExplicitValence() > totV) {
          BOOST_LOG(rdWarningLog)
              << "atom " << atom->getIdx() << " has specified valence (" << totV
              << ") smaller than the drawn valence "
              << atom->getExplicitValence() << "." << std::endl;
          atom->setNumExplicitHs(0);
        } else {
          atom->setNumExplicitHs(totV - atom->getExplicitValence());
        }
      }
    }
  }
}

}  // namespace
namespace FileParserUtils {
bool ParseV3000CTAB(std::istream *inStream, unsigned int &line, RWMol *mol,
                    Conformer *&conf, bool &chiralityPossible,
                    unsigned int &nAtoms, unsigned int &nBonds,
                    bool strictParsing, bool expectMEND) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(mol, "bad molecule");

  std::string tempStr;
  std::vector<std::string> splitLine;

  bool fileComplete = false;

  tempStr = getV3000Line(inStream, line);
  boost::to_upper(tempStr);
  if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN CTAB") {
    std::ostringstream errout;
    errout << "BEGIN CTAB line not found on line " << line;
    throw FileParseException(errout.str());
  }

  tempStr = getV3000Line(inStream, line);
  boost::to_upper(tempStr);
  if (tempStr.size() < 8 || tempStr.substr(0, 7) != "COUNTS ") {
    std::ostringstream errout;
    errout << "Bad counts line : '" << tempStr << "' on line " << line;
    throw FileParseException(errout.str());
  }
  std::string trimmed =
      boost::trim_copy(tempStr.substr(7, tempStr.length() - 7));
  boost::split(splitLine, trimmed, boost::is_any_of(" \t"),
               boost::token_compress_on);
  if (splitLine.size() < 2) {
    std::ostringstream errout;
    errout << "Bad counts line : '" << tempStr << "' on line " << line;
    throw FileParseException(errout.str());
  }

  nAtoms = FileParserUtils::toInt(splitLine[0]);
  nBonds = FileParserUtils::toInt(splitLine[1]);
  if (!nAtoms) {
    throw FileParseException("molecule has no atoms");
  }
  conf = new Conformer(nAtoms);

  unsigned int nSgroups = 0, n3DConstraints = 0, chiralFlag = 0;
  (void)chiralFlag;  // needs to be read
  if (splitLine.size() > 2) nSgroups = FileParserUtils::toInt(splitLine[2]);
  if (splitLine.size() > 3)
    n3DConstraints = FileParserUtils::toInt(splitLine[3]);
  if (splitLine.size() > 4) chiralFlag = FileParserUtils::toInt(splitLine[4]);

  ParseV3000AtomBlock(inStream, line, nAtoms, mol, conf);
  if (nBonds) {
    ParseV3000BondBlock(inStream, line, nBonds, mol, chiralityPossible);
  }

  if (nSgroups) {
    ParseV3000SGroupsBlock(inStream, line, nSgroups, mol, strictParsing);
  }

  if (n3DConstraints) {
    BOOST_LOG(rdWarningLog)
        << "3D constraint information in mol block igored at line " << line
        << std::endl;
    tempStr = getV3000Line(inStream, line);
    boost::to_upper(tempStr);
    if (tempStr.length() < 11 || tempStr.substr(0, 11) != "BEGIN OBJ3D") {
      std::ostringstream errout;
      errout << "BEGIN OBJ3D line not found on line " << line;
      throw FileParseException(errout.str());
    }
    for (unsigned int i = 0; i < n3DConstraints; ++i)
      tempStr = getV3000Line(inStream, line);
    tempStr = getV3000Line(inStream, line);
    boost::to_upper(tempStr);
    if (tempStr.length() < 9 || tempStr.substr(0, 9) != "END OBJ3D") {
      std::ostringstream errout;
      errout << "END OBJ3D line not found on line " << line;
      throw FileParseException(errout.str());
    }
  }

  tempStr = getV3000Line(inStream, line);
  // do link nodes:
  boost::to_upper(tempStr);
  while (tempStr.length() > 8 && tempStr.substr(0, 8) == "LINKNODE") {
    tempStr = getV3000Line(inStream, line);
    boost::to_upper(tempStr);
  }

  while (tempStr.length() > 5 && tempStr.substr(0, 5) == "BEGIN") {
    if (tempStr.length() > 15 && tempStr.substr(6, 10) == "COLLECTION") {
      tempStr = parseEnhancedStereo(inStream, line, mol);
    } else {
      // skip blocks we don't know how to read
      BOOST_LOG(rdWarningLog)
          << "skipping block at line " << line << ": '" << tempStr <<"'"<< std::endl;
      while (tempStr.length() < 3 || tempStr.substr(0, 3) != "END") {
        tempStr = getV3000Line(inStream, line);
      }
      tempStr = getV3000Line(inStream, line);
    }
  }

  boost::to_upper(tempStr);
  if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END CTAB") {
    throw FileParseException("END CTAB line not found");
  }

  if (expectMEND) {
    tempStr = getLine(inStream);
    ++line;
    if (tempStr[0] == 'M' && tempStr.substr(0, 6) == "M  END") {
      fileComplete = true;
    }
  } else {
    fileComplete = true;
  }

  mol->addConformer(conf, true);
  conf = nullptr;

  return fileComplete;
}

bool ParseV2000CTAB(std::istream *inStream, unsigned int &line, RWMol *mol,
                    Conformer *&conf, bool &chiralityPossible,
                    unsigned int &nAtoms, unsigned int &nBonds,
                    bool strictParsing) {
  RDUNUSED_PARAM(strictParsing);
  conf = new Conformer(nAtoms);
  if (nAtoms == 0) {
    conf->set3D(false);
  } else {
    ParseMolBlockAtoms(inStream, line, nAtoms, mol, conf);

    bool nonzeroZ = hasNonZeroZCoords(*conf);
    if (mol->hasProp(common_properties::_3DConf)) {
      conf->set3D(true);
      mol->clearProp(common_properties::_3DConf);
      if (!nonzeroZ) {
        BOOST_LOG(rdWarningLog)
            << "Warning: molecule is tagged as 3D, but all Z coords are zero"
            << std::endl;
      }
    } else {
      conf->set3D(nonzeroZ);
    }
  }
  mol->addConformer(conf, true);
  conf = nullptr;

  ParseMolBlockBonds(inStream, line, nBonds, mol, chiralityPossible);

  bool fileComplete =
      ParseMolBlockProperties(inStream, line, mol, strictParsing);
  return fileComplete;
}

}  // namespace FileParserUtils

//------------------------------------------------
//
//  Read a molecule from a stream
//
//------------------------------------------------
RWMol *MolDataStreamToMol(std::istream *inStream, unsigned int &line,
                          bool sanitize, bool removeHs, bool strictParsing) {
  PRECONDITION(inStream, "no stream");
  std::string tempStr;
  bool fileComplete = false;
  bool chiralityPossible = false;
  Utils::LocaleSwitcher ls;
  // mol name
  line++;
  tempStr = getLine(inStream);
  if (inStream->eof()) {
    return nullptr;
  }
  auto *res = new RWMol();
  res->setProp(common_properties::_Name, tempStr);

  // info
  line++;
  tempStr = getLine(inStream);
  res->setProp("_MolFileInfo", tempStr);
  if (tempStr.length() >= 22) {
    std::string dimLabel = tempStr.substr(20, 2);
    // Unless labelled as 3D we assume 2D
    if (dimLabel == "3d" || dimLabel == "3D") {
      res->setProp(common_properties::_3DConf, 1);
    }
  }
  // comments
  line++;
  tempStr = getLine(inStream);
  res->setProp("_MolFileComments", tempStr);

  unsigned int nAtoms = 0, nBonds = 0, nLists = 0, chiralFlag = 0, nsText = 0,
               nRxnComponents = 0;
  int nReactants = 0, nProducts = 0, nIntermediates = 0;
  (void)nLists;  // read from the file but unused
  (void)nsText;
  (void)nRxnComponents;
  (void)nReactants;
  (void)nProducts;
  (void)nIntermediates;
  // counts line, this is where we really get started
  line++;
  tempStr = getLine(inStream);

  if (tempStr.size() < 6) {
    if (res) {
      delete res;
      res = nullptr;
    }
    std::ostringstream errout;
    errout << "Counts line too short: '" << tempStr << "' on line" << line;
    throw FileParseException(errout.str());
  }

  unsigned int spos = 0;
  // this needs to go into a try block because if the lexical_cast throws an
  // exception we want to catch and delete mol before leaving this function
  try {
    nAtoms = FileParserUtils::toInt(tempStr.substr(spos, 3), true);
    spos = 3;
    nBonds = FileParserUtils::toInt(tempStr.substr(spos, 3), true);
    spos = 6;
  } catch (boost::bad_lexical_cast &) {
    if (res) {
      delete res;
      res = nullptr;
    }
    std::ostringstream errout;
    errout << "Cannot convert '" << tempStr.substr(spos, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  try {
    spos = 6;
    if (tempStr.size() >= 9)
      nLists = FileParserUtils::toInt(tempStr.substr(spos, 3), true);

    spos = 12;
    if (tempStr.size() >= spos + 3)
      chiralFlag = FileParserUtils::toInt(tempStr.substr(spos, 3), true);

    spos = 15;
    if (tempStr.size() >= spos + 3)
      nsText = FileParserUtils::toInt(tempStr.substr(spos, 3), true);

    spos = 18;
    if (tempStr.size() >= spos + 3)
      nRxnComponents = FileParserUtils::toInt(tempStr.substr(spos, 3), true);

    spos = 21;
    if (tempStr.size() >= spos + 3)
      nReactants = FileParserUtils::toInt(tempStr.substr(spos, 3), true);

    spos = 24;
    if (tempStr.size() >= spos + 3)
      nProducts = FileParserUtils::toInt(tempStr.substr(spos, 3), true);

    spos = 27;
    if (tempStr.size() >= spos + 3)
      nIntermediates = FileParserUtils::toInt(tempStr.substr(spos, 3), true);

  } catch (boost::bad_lexical_cast &) {
    // some SD files (such as some from NCI) lack all the extra information
    // on the header line, so ignore problems parsing there.
  }

  unsigned int ctabVersion = 2000;
  if (tempStr.size() > 35) {
    if (tempStr.size() < 39 || tempStr[34] != 'V') {
      std::ostringstream errout;
      errout << "CTAB version string invalid at line " << line;
      if (strictParsing) {
        delete res;
        res = nullptr;
        throw FileParseException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
      }
    } else if (tempStr.substr(34, 5) == "V3000") {
      ctabVersion = 3000;
    } else if (tempStr.substr(34, 5) != "V2000") {
      std::ostringstream errout;
      errout << "Unsupported CTAB version: '" << tempStr.substr(34, 5)
             << "' at line " << line;
      if (strictParsing) {
        delete res;
        res = nullptr;
        throw FileParseException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
      }
    }
  }

  if (chiralFlag) {
    res->setProp(common_properties::_MolFileChiralFlag, chiralFlag);
  }

  Conformer *conf = nullptr;
  try {
    if (ctabVersion == 2000) {
      fileComplete = FileParserUtils::ParseV2000CTAB(inStream, line, res, conf,
                                                     chiralityPossible, nAtoms,
                                                     nBonds, strictParsing);
    } else {
      if (nAtoms != 0 || nBonds != 0) {
        std::ostringstream errout;
        errout << "V3000 mol blocks should have 0s in the initial counts line. "
                  "(line: "
               << line << ")";
        if (strictParsing) {
          delete res;
          res = nullptr;
          throw FileParseException(errout.str());
        } else {
          BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
        }
      }
      fileComplete = FileParserUtils::ParseV3000CTAB(inStream, line, res, conf,
                                                     chiralityPossible, nAtoms,
                                                     nBonds, strictParsing);
    }
  } catch (MolFileUnhandledFeatureException &e) {
    // unhandled mol file feature, just delete the result
    delete res;
    delete conf;
    res = nullptr;
    conf = nullptr;
    BOOST_LOG(rdErrorLog) << " Unhandled CTAB feature: '" << e.message()
                          << "'. Molecule skipped." << std::endl;

    if (!inStream->eof()) tempStr = getLine(inStream);
    ++line;
    while (!inStream->eof() && tempStr.substr(0, 6) != "M  END" &&
           tempStr.substr(0, 4) != "$$$$") {
      tempStr = getLine(inStream);
      ++line;
    }
    if (!inStream->eof() || tempStr.substr(0, 6) == "M  END" ||
        tempStr.substr(0, 4) == "$$$$")
      fileComplete = true;
    else
      fileComplete = false;
  } catch (FileParseException &e) {
    // catch our exceptions and throw them back after cleanup
    delete res;
    delete conf;
    res = nullptr;
    conf = nullptr;
    throw e;
  }

  if (!fileComplete) {
    delete res;
    delete conf;
    res = nullptr;
    conf = nullptr;
    std::ostringstream errout;
    errout
        << "Problems encountered parsing Mol data, M  END missing around line "
        << line;
    throw FileParseException(errout.str());
  }

  if (res) {
    res->clearAllAtomBookmarks();
    res->clearAllBondBookmarks();

    // calculate explicit valence on each atom:
    for (RWMol::AtomIterator atomIt = res->beginAtoms();
         atomIt != res->endAtoms(); ++atomIt) {
      (*atomIt)->calcExplicitValence(false);
    }

    // postprocess mol file flags
    ProcessMolProps(res);

    // update the chirality and stereo-chemistry
    //
    // NOTE: we detect the stereochemistry before sanitizing/removing
    // hydrogens because the removal of H atoms may actually remove
    // the wedged bond from the molecule.  This wipes out the only
    // sign that chirality ever existed and makes us sad... so first
    // perceive chirality, then remove the Hs and sanitize.
    //
    const Conformer &conf = res->getConformer();
    if (chiralityPossible) {
      if (!conf.is3D()) {
        DetectAtomStereoChemistry(*res, &conf);
      } else {
        res->updatePropertyCache(false);
        MolOps::assignChiralTypesFrom3D(*res, conf.getId(), true);
      }
    }

    if (sanitize) {
      try {
        if (removeHs) {
          MolOps::removeHs(*res, false, false);
        } else {
          MolOps::sanitizeMol(*res);
        }
        // now that atom stereochem has been perceived, the wedging
        // information is no longer needed, so we clear
        // single bond dir flags:
        ClearSingleBondDirFlags(*res);

        // unlike DetectAtomStereoChemistry we call detectBondStereochemistry
        // here after sanitization because we need the ring information:
        MolOps::detectBondStereochemistry(*res);
      } catch (...) {
        delete res;
        res = nullptr;
        throw;
      }
      MolOps::assignStereochemistry(*res, true, true, true);
    } else {
      // we still need to do something about double bond stereochemistry
      // (was github issue 337)
      // now that atom stereochem has been perceived, the wedging
      // information is no longer needed, so we clear
      // single bond dir flags:
      ClearSingleBondDirFlags(*res);
      MolOps::detectBondStereochemistry(*res);
    }

    if (res->hasProp(common_properties::_NeedsQueryScan)) {
      res->clearProp(common_properties::_NeedsQueryScan);
      completeMolQueries(res);
    }
  }
  return res;
}

RWMol *MolDataStreamToMol(std::istream &inStream, unsigned int &line,
                          bool sanitize, bool removeHs, bool strictParsing) {
  return MolDataStreamToMol(&inStream, line, sanitize, removeHs, strictParsing);
}
//------------------------------------------------
//
//  Read a molecule from a string
//
//------------------------------------------------
RWMol *MolBlockToMol(const std::string &molBlock, bool sanitize, bool removeHs,
                     bool strictParsing) {
  std::istringstream inStream(molBlock);
  unsigned int line = 0;
  return MolDataStreamToMol(inStream, line, sanitize, removeHs, strictParsing);
}

//------------------------------------------------
//
//  Read a molecule from a file
//
//------------------------------------------------
RWMol *MolFileToMol(const std::string &fName, bool sanitize, bool removeHs,
                    bool strictParsing) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  RWMol *res = nullptr;
  if (!inStream.eof()) {
    unsigned int line = 0;
    res = MolDataStreamToMol(inStream, line, sanitize, removeHs, strictParsing);
  }
  return res;
}
}  // namespace RDKit
