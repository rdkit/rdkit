//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
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
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include "FileParsers.h"
#include "FileParserUtils.h"
#include "MolSGroupParsing.h"
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/StereoGroup.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/Atropisomers.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/GenericGroups/GenericGroups.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/Chirality.h>

#include <fstream>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <typeinfo>
#include <exception>
#include <charconv>
#include <regex>
#include <sstream>
#include <locale>
#include <cstdlib>
#include <cstdio>
#include <string_view>

using namespace RDKit::SGroupParsing;
using std::regex;
using std::regex_match;
using std::smatch;

namespace RDKit {

namespace FileParserUtils {

int toInt(const std::string_view input, bool acceptSpaces) {
  // don't need to worry about locale stuff here because
  // we're not going to have delimiters

  // sanity check on the input since strtol doesn't do it for us:
  const char *txt = input.data();
  for (size_t i = 0u; i < input.size() && *txt != '\x00'; ++i) {
    if ((*txt >= '0' && *txt <= '9') || (acceptSpaces && *txt == ' ') ||
        *txt == '+' || *txt == '-') {
      ++txt;
    } else {
      throw boost::bad_lexical_cast();
    }
  }
  // remove leading spaces
  txt = input.data();
  unsigned int sz = input.size();
  if (acceptSpaces) {
    while (*txt == ' ') {
      ++txt;
      --sz;
      // have we run off the end of the view?
      if (sz < 1U) {
        return 0;
      }
    }
  }
  int res = 0;
  std::from_chars(txt, txt + sz, res);

  return res;
}
int toInt(const std::string &input, bool acceptSpaces) {
  return toInt(std::string_view(input.c_str()), acceptSpaces);
}
unsigned int toUnsigned(const std::string_view input, bool acceptSpaces) {
  // don't need to worry about locale stuff here because
  // we're not going to have delimiters

  // sanity check on the input since strtol doesn't do it for us:
  const char *txt = input.data();
  for (size_t i = 0u; i < input.size() && *txt != '\x00'; ++i) {
    if ((*txt >= '0' && *txt <= '9') || (acceptSpaces && *txt == ' ') ||
        *txt == '+') {
      ++txt;
    } else {
      throw boost::bad_lexical_cast();
    }
  }
  // remove leading spaces
  txt = input.data();
  unsigned int sz = input.size();
  if (acceptSpaces) {
    while (*txt == ' ') {
      ++txt;
      --sz;
      // have we run off the end of the view?
      if (sz < 1U) {
        return 0;
      }
    }
  }
  unsigned int res = 0;
  std::from_chars(txt, txt + sz, res);
  return res;
}
unsigned int toUnsigned(const std::string &input, bool acceptSpaces) {
  return toUnsigned(std::string_view(input.c_str()), acceptSpaces);
}
double toDouble(const std::string_view input, bool acceptSpaces) {
  // sanity check on the input since strtol doesn't do it for us:
  const char *txt = input.data();
  for (size_t i = 0u; i < input.size() && *txt != '\x00'; ++i) {
    // check for ',' and '.' because locale
    if ((*txt >= '0' && *txt <= '9') || (acceptSpaces && *txt == ' ') ||
        *txt == '+' || *txt == '-' || *txt == ',' || *txt == '.') {
      ++txt;
    } else {
      throw boost::bad_lexical_cast();
    }
  }
  // unfortunately from_chars() with doubles didn't work on g++ until v11.1
  // and the status with clang is hard to figure out... we remain old-school
  // remove leading spaces
  double res = atof(input.data());
  return res;
}
double toDouble(const std::string &input, bool acceptSpaces) {
  return toDouble(std::string_view(input.c_str()), acceptSpaces);
}
std::string getV3000Line(std::istream *inStream, unsigned int &line) {
  // FIX: technically V3K blocks are case-insensitive. We should really be
  // up-casing everything here.
  PRECONDITION(inStream, "bad stream");
  std::string res;
  ++line;
  auto inl = getLine(inStream);
  std::string_view tempStr = inl;
  if (tempStr.size() < 7 || tempStr.substr(0, 7) != "M  V30 ") {
    std::ostringstream errout;
    errout << "Line " << line << " does not start with 'M  V30 '" << std::endl;
    throw FileParseException(errout.str());
  }
  // FIX: do we need to handle trailing whitespace after a -?
  while (tempStr.back() == '-') {
    // continuation character, append what we read:
    res += tempStr.substr(7, tempStr.length() - 8);
    // and then read another line:
    ++line;
    inl = getLine(inStream);
    tempStr = inl;
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
  return QueryOps::replaceAtomWithQueryAtom(mol, atom);
}
}  // namespace FileParserUtils
using RDKit::FileParserUtils::getV3000Line;

namespace {

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
      R"regex(MDLV30/STE(...)([0-9]*) +ATOMS=\(([0-9]+) +(.*)\) *)regex");

  smatch match;
  std::vector<StereoGroup> groups;

  // Read the collection until the end
  auto tempStr = getV3000Line(inStream, line);
  boost::to_upper(tempStr);
  while (!startsWith(tempStr, "END", 3)) {
    // If this line in the collection is part of a stereo group
    if (regex_match(tempStr, match, stereo_label)) {
      StereoGroupType grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
      unsigned groupid = 0;

      if (match[1] == "ABS") {
        grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
      } else if (match[1] == "REL") {
        grouptype = RDKit::StereoGroupType::STEREO_OR;
        groupid = FileParserUtils::toUnsigned(match[2], true);
      } else if (match[1] == "RAC") {
        grouptype = RDKit::StereoGroupType::STEREO_AND;
        groupid = FileParserUtils::toUnsigned(match[2], true);
      } else {
        std::ostringstream errout;
        errout << "Unrecognized stereogroup type : '" << tempStr << "' on line"
               << line;
        throw FileParseException(errout.str());
      }

      const unsigned int count = FileParserUtils::toUnsigned(match[3], true);
      std::vector<Atom *> atoms;
      std::stringstream ss(match[4]);
      unsigned int index;
      for (size_t i = 0; i < count; ++i) {
        ss >> index;
        // atoms are 1 indexed in molfiles
        atoms.push_back(mol->getAtomWithIdx(index - 1));
      }
      std::vector<Bond *> newBonds;
      groups.emplace_back(grouptype, std::move(atoms), std::move(newBonds),
                          groupid);
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

void ParseOldAtomList(RWMol *mol, const std::string_view &text,
                      unsigned int line) {
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
      delete q;
      std::ostringstream errout;
      errout << "Unrecognized atom-list query modifier: '" << text[4]
             << "' on line " << line;
      throw FileParseException(errout.str());
  }

  int nQueries;
  try {
    nQueries = FileParserUtils::toInt(text.substr(9, 1));
  } catch (const std::out_of_range &) {
    delete q;
    std::ostringstream errout;
    errout << "Cannot convert position 9 of '" << text << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  } catch (boost::bad_lexical_cast &) {
    delete q;
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
    } catch (const std::out_of_range &) {
      delete q;
      std::ostringstream errout;
      errout << "Cannot convert position " << pos << " of '" << text
             << "' to int on line " << line;
      throw FileParseException(errout.str());
    } catch (boost::bad_lexical_cast &) {
      delete q;
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(pos, 3) << "' to int on line "
             << line;
      throw FileParseException(errout.str());
    }
    RANGE_CHECK(0, atNum, 200);  // goofy!
    q->addChild(
        QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(atNum)));
    if (!i) {
      a.setAtomicNum(atNum);
    }
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
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
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
        case 0:
          // This shouldn't be required, but let's make sure.
          mol->getAtomWithIdx(aid - 1)->setNumRadicalElectrons(0);
          break;
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
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParsePXALine(RWMol *mol, const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  PXA", "bad PXA line");
  unsigned int pos = 7;
  try {
    auto atIdx =
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
        }
      }
      spos += 4;
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
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
    int count = 0;
    try {
      aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      Atom *atom = mol->getAtomWithIdx(aid - 1);
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        count = FileParserUtils::toInt(text.substr(spos, 4));
      }
      spos += 4;
      if (count == 0) {
        continue;
      }
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
        atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
      }
      atom->expandQuery(q, Queries::COMPOSITE_AND);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
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
    int count = 0;
    try {
      aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      Atom *atom = mol->getAtomWithIdx(aid - 1);
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        count = FileParserUtils::toInt(text.substr(spos, 4));
      }
      spos += 4;
      if (count == 0) {
        continue;
      } else if (count == 1) {
        ATOM_EQUALS_QUERY *q = makeAtomUnsaturatedQuery();
        if (!atom->hasQuery()) {
          atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
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
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
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
    int count = 0;
    try {
      aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      Atom *atom = mol->getAtomWithIdx(aid - 1);
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        count = FileParserUtils::toInt(text.substr(spos, 4));
      }
      spos += 4;
      if (count == 0) {
        continue;
      }
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
        atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
      }
      atom->expandQuery(q, Queries::COMPOSITE_AND);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
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
      if (!aid || aid > mol->getNumAtoms()) {
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
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
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
      if (!aid || aid > mol->getNumAtoms()) {
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
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
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
      if (!bid || bid > mol->getNumBonds()) {
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
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
      throw FileParseException(errout.str());
    }
  }
}

void ParseMarvinSmartsLine(RWMol *mol, const std::string &text,
                           unsigned int line) {
  const unsigned int atomNumStart = 10;
  const unsigned int smartsStart = 15;
  // M  MRV SMA   1 [*;A]
  // 01234567890123456789
  //           1111111111
  if (text.substr(0, 10) != "M  MRV SMA") {
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
  RWMol *m = nullptr;
  try {
    m = SmartsToMol(sma);
  } catch (...) {
    // Is this ever used?
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

void ParseAttachPointLine(RWMol *mol, const std::string &text,
                          unsigned int line, bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  APO"), "bad APO line");

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
      if (!aid || aid > mol->getNumAtoms()) {
        std::ostringstream errout;
        errout << "Bad APO specification on line " << line;
        throw FileParseException(errout.str());
      }
      spos += 4;
      --aid;
      Atom *atom = mol->getAtomWithIdx(aid);
      if (!atom) {
        std::ostringstream errout;
        errout << "Atom " << aid << " from APO specification on line " << line
               << " not found";
        throw FileParseException(errout.str());
      } else {
        if (val < 0 || val > 3) {
          std::ostringstream errout;
          errout << "Value " << val << " from APO specification on line "
                 << line << " is invalid";
          throw FileParseException(errout.str());
        } else if (val) {
          if (val == 3) {
            // this is -1 in v3k mol blocks, so use that:
            val = -1;
          }
          if (atom->hasProp(common_properties::molAttachPoint)) {
            std::ostringstream errout;
            errout << "Multiple ATTCHPT values for atom " << atom->getIdx() + 1
                   << " on line " << line;
            if (strictParsing) {
              throw FileParseException(errout.str());
            } else {
              BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
            }
          } else {
            atom->setProp(common_properties::molAttachPoint, val);
          }
        }
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
      throw FileParseException(errout.str());
    }
  }
}

// the format differs between V2000 and V3000, so we have to do a bit of
// translation here
void ParseLinkNodeLine(RWMol *mol, const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == std::string("M  LIN"), "bad LIN line");

  unsigned int nent;
  try {
    nent = FileParserUtils::stripSpacesAndCast<unsigned int>(text.substr(6, 3));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(6, 3) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  std::string propVal = "";
  unsigned int spos = 9;
  for (unsigned int ie = 0; ie < nent; ie++) {
    try {
      auto aid = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      if (!aid || aid > mol->getNumAtoms()) {
        std::ostringstream errout;
        errout << "LIN specification has bad atom idx on line " << line;
        throw FileParseException(errout.str());
      }
      spos += 4;

      if (text.size() < spos + 4 || text.substr(spos, 4) == "    ") {
        std::ostringstream errout;
        errout << "LIN specification missing repeat count on line " << line;
        throw FileParseException(errout.str());
      }
      auto repeatCount = FileParserUtils::stripSpacesAndCast<unsigned int>(
          text.substr(spos, 4));
      spos += 4;
      if (repeatCount < 2) {
        std::ostringstream errout;
        errout << "LIN specification: repeat count must be >=2 on line "
               << line;
        throw FileParseException(errout.str());
      }
      unsigned int substB = 0;
      unsigned int substC = 0;
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        substB = FileParserUtils::stripSpacesAndCast<unsigned int>(
            text.substr(spos, 4));
      }
      spos += 4;
      if (text.size() >= spos + 4 && text.substr(spos, 4) != "    ") {
        substC = FileParserUtils::stripSpacesAndCast<unsigned int>(
            text.substr(spos, 4));
      }
      spos += 4;

      if (!substB || substB > mol->getNumAtoms() ||
          substC > mol->getNumAtoms()) {
        std::ostringstream errout;
        errout << "LIN specification has bad substituent idx on line " << line;
        throw FileParseException(errout.str());
      }

      boost::format formatter;
      if (substC) {
        formatter = boost::format("1 %1% 2 %2% %3% %2% %4%") % repeatCount %
                    aid % substB % substC;
      } else {
        formatter = boost::format("1 %1% 1 %2% %3%") % repeatCount % aid %
                    substB % substC;
      }
      if (!propVal.empty()) {
        propVal += "|";
      }
      propVal += formatter.str();
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(spos, 4)
             << "' to int on line " << line;
      throw FileParseException(errout.str());
    }
    mol->setProp(common_properties::molFileLinkNodes, propVal);
  }
}

// Recursively populates queryVect with COMPOSITE_AND queries
// present in the input query. If the logic of the input query
// is more complex, it returns nullptr and empty set.
// The returned ptr should only be checked for not being null
// and not used for any other purposes, as the actual result is
// the queryVect
const QueryAtom::QUERYATOM_QUERY *getAndQueries(
    const QueryAtom::QUERYATOM_QUERY *q,
    std::vector<const QueryAtom::QUERYATOM_QUERY *> &queryVect) {
  if (q) {
    auto qOrig = q;
    for (auto cq = qOrig->beginChildren(); cq != qOrig->endChildren(); ++cq) {
      if (q == qOrig && q->getDescription() != "AtomAnd") {
        q = nullptr;
        break;
      }
      q = getAndQueries(cq->get(), queryVect);
    }
    if (q == qOrig) {
      queryVect.push_back(q);
    }
  }
  if (!q) {
    queryVect.clear();
  }
  return q;
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
  QueryAtom *a = nullptr;
  QueryAtom *qaOrig = nullptr;
  QueryAtom::QUERYATOM_QUERY *qOrig = nullptr;
  Atom *aOrig = mol->getAtomWithIdx(idx);
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
      if (aOrig->hasQuery()) {
        qaOrig = dynamic_cast<QueryAtom *>(aOrig);
        if (qaOrig) {
          qOrig = qaOrig->getQuery();
        }
      }
      a = new QueryAtom(*aOrig);
      a->setAtomicNum(atNum);
      if (!qOrig) {
        qOrig = a->getQuery()->copy();
      }
      a->setQuery(makeAtomNumQuery(atNum));
    } else {
      a->expandQuery(makeAtomNumQuery(atNum), Queries::COMPOSITE_OR, true);
      // For COMPOSITE_OR query atoms, reset atomic num to 0 such that they are
      // exported as "*" in SMILES
      a->setAtomicNum(0);
    }
  }
  ASSERT_INVARIANT(a, "no atom built");
  if (qOrig) {
    std::vector<const QueryAtom::QUERYATOM_QUERY *> queryVect;
    if (getAndQueries(qOrig, queryVect)) {
      for (const auto &q : queryVect) {
        if (q->getDescription() != "AtomAtomicNum") {
          a->expandQuery(q->copy(), Queries::COMPOSITE_AND, true);
        }
      }
    }
    if (!qaOrig) {
      delete qOrig;
    }
  }
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

void ParseV3000RGroups(RWMol *mol, Atom *&atom, std::string_view text,
                       unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(atom, "bad atom");
  if (text[0] != '(' || text.back() != ')') {
    std::ostringstream errout;
    errout << "Bad RGROUPS specification '" << text << "' on line " << line
           << ". Missing parens.";
    throw FileParseException(errout.str());
  }
  std::vector<std::string> splitToken;
  std::string resid = std::string(text.substr(1, text.size() - 2));
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
    atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
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

namespace {
void setRGPProps(const std::string_view symb, Atom *res) {
  PRECONDITION(res, "bad atom pointer");
  // set the dummy label so that this is shown correctly
  // in other pieces of the code :
  std::string symbc(symb);
  res->setProp(common_properties::dummyLabel, symbc);
}

void lookupAtomicNumber(Atom *res, const std::string &symb,
                        bool strictParsing) {
  try {
    res->setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(symb));
  } catch (const Invar::Invariant &e) {
    if (strictParsing || symb.empty()) {
      throw FileParseException(e.what());
    } else {
      res->setAtomicNum(0);
      res->setProp(common_properties::dummyLabel, symb);
    }
  }
}

}  // namespace

Atom *ParseMolFileAtomLine(const std::string_view text, RDGeom::Point3D &pos,
                           unsigned int line, bool strictParsing) {
  std::string symb;
  int massDiff, chg, hCount;

  if ((strictParsing && text.size() < 34) || text.size() < 32) {
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
      errout << "Cannot convert '" << text.substr(34, 2) << "' to int on line "
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
  std::unique_ptr<Atom> res(new Atom);
  bool isComplexQueryName =
      std::find(complexQueries.begin(), complexQueries.end(), symb) !=
      complexQueries.end();
  if (isComplexQueryName || symb == "L" || symb == "*" || symb == "LP" ||
      symb == "R" || symb == "R#" ||
      (symb[0] == 'R' && symb >= "R0" && symb <= "R99")) {
    if (isComplexQueryName || symb == "*" || symb == "R") {
      auto *query = new QueryAtom(0);
      if (symb == "*" || symb == "R") {
        // according to the MDL spec, these match anything
        query->setQuery(makeAtomNullQuery());
      } else if (isComplexQueryName) {
        convertComplexNameToQuery(query, symb);
      }
      res.reset(query);
      // queries have no implicit Hs:
      res->setNoImplicit(true);
    } else {
      res->setAtomicNum(0);
    }
    if (massDiff == 0 && symb[0] == 'R') {
      if (symb.length() > 1 && symb >= "R0" && symb <= "R99") {
        std::string rlabel = "";
        rlabel = symb.substr(1, symb.length() - 1);
        int rnumber;
        try {
          rnumber = boost::lexical_cast<int>(rlabel);
        } catch (boost::bad_lexical_cast &) {
          rnumber = -1;
        }
        if (rnumber >= 0) {
          res->setIsotope(rnumber);
        }
      }
    }
    if (symb[0] == 'R') {
      // we used to skip R# here because that really should be handled by an
      // RGP spec, but that turned out to not be permissive enough... <sigh>
      setRGPProps(symb, res.get());
    }
  } else if (symb == "D") {  // mol blocks support "D" and "T" as shorthand...
                             // handle that.
    res->setAtomicNum(1);
    res->setIsotope(2);
  } else if (symb == "T") {  // mol blocks support "D" and "T" as shorthand...
                             // handle that.
    res->setAtomicNum(1);
    res->setIsotope(3);
  } else if (symb == "Pol" || symb == "Mod") {
    res->setAtomicNum(0);
    res->setProp(common_properties::dummyLabel, symb);
  } else if (GenericGroups::genericMatchers.find(symb) !=
             GenericGroups::genericMatchers.end()) {
    res.reset(new QueryAtom(0));
    res->setProp(common_properties::atomLabel, std::string(symb));
  } else {
    if (symb.size() == 2 && symb[1] >= 'A' && symb[1] <= 'Z') {
      symb[1] = static_cast<char>(tolower(symb[1]));
    }
    lookupAtomicNumber(res.get(), symb, strictParsing);
  }

  // res->setPos(pX,pY,pZ);
  if (chg != 0) {
    res->setFormalCharge(4 - chg);
  }

  if (hCount >= 1) {
    if (!res->hasQuery()) {
      auto qatom = new QueryAtom(*res);
      res.reset(qatom);
    }
    res->setNoImplicit(true);
    if (hCount > 1) {
      ATOM_EQUALS_QUERY *oq = makeAtomImplicitHCountQuery(hCount - 1);
      auto nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(
          hCount - 1, oq->getDataFunc(),
          std::string("less_") + oq->getDescription());
      res->expandQuery(nq);
      delete oq;
    } else {
      res->expandQuery(makeAtomImplicitHCountQuery(0));
    }
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
  }

  if (text.size() >= 42 && text.substr(39, 3) != "  0") {
    int parity = 0;
    try {
      parity = FileParserUtils::toInt(text.substr(39, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(39, 3) << "' to int on line "
             << line;
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
      throw FileParseException(errout.str());
    }
    res->setProp(common_properties::molStereoCare, stereoCare);
  }
  if (text.size() >= 51 && text.substr(48, 3) != "  0") {
    int totValence = 0;
    try {
      totValence = FileParserUtils::toInt(text.substr(48, 3), true);
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert '" << text.substr(48, 3) << "' to int on line "
             << line;
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
      throw FileParseException(errout.str());
    }
    res->setProp("molExactChangeFlag", exactChangeFlag);
  }
  return res.release();
}

Bond *ParseMolFileBondLine(const std::string_view text, unsigned int line) {
  unsigned int idx1, idx2, bType, stereo;
  int spos = 0;

  if (text.size() < 9) {
    std::ostringstream errout;
    errout << "Bond line too short: '" << text << "' on line " << line;
    throw FileParseException(errout.str());
  }

  try {
    idx1 = FileParserUtils::toUnsigned(text.substr(spos, 3));
    spos += 3;
    idx2 = FileParserUtils::toUnsigned(text.substr(spos, 3));
    spos += 3;
    bType = FileParserUtils::toUnsigned(text.substr(spos, 3));
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
    case 9:
      type = Bond::DATIVE;
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
      } else if (bType == 5) {
        res->setQuery(makeSingleOrDoubleBondQuery());
        res->setProp(common_properties::_MolFileBondQuery, 1);
      } else if (bType == 6) {
        res->setQuery(makeSingleOrAromaticBondQuery());
        res->setProp(common_properties::_MolFileBondQuery, 1);
      } else if (bType == 7) {
        res->setQuery(makeDoubleOrAromaticBondQuery());
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
      stereo = FileParserUtils::toUnsigned(text.substr(9, 3));
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
}  // namespace

void ParseMolBlockAtoms(std::istream *inStream, unsigned int &line,
                        unsigned int nAtoms, RWMol *mol, Conformer *conf,
                        bool strictParsing) {
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
    Atom *atom = ParseMolFileAtomLine(tempStr, pos, line, strictParsing);
    unsigned int aid = mol->addAtom(atom, false, true);
    conf->setAtomPos(aid, pos);
    mol->setAtomBookmark(atom, i);
  }
}

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
    }
    // if the bond might have chirality info associated with it, set a flag:
    if (bond->getBondDir() != Bond::NONE &&
        bond->getBondDir() != Bond::UNKNOWN) {
      chiralityPossible = true;
    }
    // v2k has no way to set stereoCare on bonds, so set the property if both
    // the beginning and end atoms have it set:
    int care1 = 0;
    int care2 = 0;
    if (!bond->hasProp(common_properties::molStereoCare) &&
        mol->getAtomWithIdx(bond->getBeginAtomIdx())
            ->getPropIfPresent(common_properties::molStereoCare, care1) &&
        mol->getAtomWithIdx(bond->getEndAtomIdx())
            ->getPropIfPresent(common_properties::molStereoCare, care2)) {
      if (care1 && care2) {
        bond->setProp(common_properties::molStereoCare, 1);
      }
    }
    mol->addBond(bond, true);
    mol->setBondBookmark(bond, i);
  }
}

bool checkAttachmentPointsAreValid(
    const RWMol *mol, std::pair<const int, SubstanceGroup> &sgroup) {
  bool res = true;
  int nAtoms = static_cast<int>(mol->getNumAtoms());
  std::vector<SubstanceGroup::AttachPoint> &attachPoints =
      sgroup.second.getAttachPoints();
  for (auto &attachPoint : attachPoints) {
    if (attachPoint.lvIdx == nAtoms) {
      const std::vector<unsigned int> &bonds = sgroup.second.getBonds();
      if (bonds.size() == 1) {
        const auto bond = mol->getBondWithIdx(bonds.front());
        if (bond->getBeginAtomIdx() == attachPoint.aIdx ||
            bond->getEndAtomIdx() == attachPoint.aIdx) {
          attachPoint.lvIdx = bond->getOtherAtomIdx(attachPoint.aIdx);
        }
      }
    }
    if (attachPoint.lvIdx == nAtoms) {
      BOOST_LOG(rdWarningLog)
          << "Could not infer missing lvIdx on malformed SAP line for SGroup "
          << sgroup.first << std::endl;
      res = false;
    }
  }
  return res;
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
      ParseOldAtomList(mol, std::string_view(tempStr.c_str()), line);
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
  while (!inStream->eof() && !inStream->fail() && lineBeg != "M  END" &&
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
    } else if (lineBeg == "M  ALS") {
      ParseNewAtomList(mol, tempStr, line);
    } else if (lineBeg == "M  ISO") {
      ParseIsotopeLine(mol, tempStr, line);
    } else if (lineBeg == "M  RGP") {
      ParseRGroupLabels(mol, tempStr, line);
    } else if (lineBeg == "M  RBC") {
      ParseRingBondCountLine(mol, tempStr, line);
    } else if (lineBeg == "M  SUB") {
      ParseSubstitutionCountLine(mol, tempStr, line);
    } else if (lineBeg == "M  UNS") {
      ParseUnsaturationLine(mol, tempStr, line);
    } else if (lineBeg == "M  CHG") {
      ParseChargeLine(mol, tempStr, firstChargeLine, line);
      firstChargeLine = false;
    } else if (lineBeg == "M  RAD") {
      ParseRadicalLine(mol, tempStr, firstChargeLine, line);
      firstChargeLine = false;
    } else if (lineBeg == "M  PXA") {
      ParsePXALine(mol, tempStr, line);

      /* SGroup parsing start */
    } else if (lineBeg == "M  STY") {
      ParseSGroupV2000STYLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SST") {
      ParseSGroupV2000SSTLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SLB") {
      ParseSGroupV2000SLBLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SCN") {
      ParseSGroupV2000SCNLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SDS") {
      ParseSGroupV2000SDSLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SAL" || lineBeg == "M  SBL" ||
               lineBeg == "M  SPA") {
      ParseSGroupV2000VectorDataLine(sGroupMap, mol, tempStr, line,
                                     strictParsing);
    } else if (lineBeg == "M  SMT") {
      ParseSGroupV2000SMTLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SDI") {
      ParseSGroupV2000SDILine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  CRS") {
      std::ostringstream errout;
      errout << "Unsupported SGroup subtype '" << lineBeg << "' on line "
             << line;
      throw FileParseException(errout.str());
    } else if (lineBeg == "M  SBV") {
      ParseSGroupV2000SBVLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SDT") {
      ParseSGroupV2000SDTLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SDD") {
      ParseSGroupV2000SDDLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SCD" || lineBeg == "M  SED") {
      ParseSGroupV2000SCDSEDLine(sGroupMap, dataFieldsMap, mol, tempStr, line,
                                 strictParsing, SCDcounter, lastDataSGroup,
                                 currentDataField);
    } else if (lineBeg == "M  SPL") {
      ParseSGroupV2000SPLLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SNC") {
      ParseSGroupV2000SNCLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SAP") {
      ParseSGroupV2000SAPLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SCL") {
      ParseSGroupV2000SCLLine(sGroupMap, mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  SBT") {
      ParseSGroupV2000SBTLine(sGroupMap, mol, tempStr, line, strictParsing);

      /* SGroup parsing end */
    } else if (lineBeg == "M  ZBO") {
      ParseZBOLine(mol, tempStr, line);
    } else if (lineBeg == "M  ZCH") {
      ParseZCHLine(mol, tempStr, line);
    } else if (lineBeg == "M  HYD") {
      ParseHYDLine(mol, tempStr, line);
    } else if (lineBeg == "M  MRV") {
      ParseMarvinSmartsLine(mol, tempStr, line);
    } else if (lineBeg == "M  APO") {
      ParseAttachPointLine(mol, tempStr, line, strictParsing);
    } else if (lineBeg == "M  LIN") {
      ParseLinkNodeLine(mol, tempStr, line);
    }
    line++;
    tempStr = getLine(inStream);
    lineBeg = tempStr.substr(0, 6);
  }
  if (tempStr[0] == 'M' && tempStr.substr(0, 6) == "M  END") {
    // All went well, make final updates to SGroups, and add them to Mol
    for (auto &sgroup : sGroupMap) {
      if (sgroup.second.getIsValid()) {
        sgroup.second.setProp("DATAFIELDS", dataFieldsMap[sgroup.first]);
        sgroup.second.setIsValid(checkAttachmentPointsAreValid(mol, sgroup));
      }
      if (sgroup.second.getIsValid()) {
        addSubstanceGroup(*mol, sgroup.second);
      } else {
        std::ostringstream errout;
        errout << "SGroup " << sgroup.first << " is invalid";
        if (strictParsing) {
          throw FileParseException(errout.str());
        } else {
          BOOST_LOG(rdWarningLog)
              << errout.str() << " and will be ignored" << std::endl;
        }
      }
    }

    fileComplete = true;
  }
  return fileComplete;
}

Atom *ParseV3000AtomSymbol(std::string_view token, unsigned int &line,
                           bool strictParsing) {
  bool negate = false;
  token = FileParserUtils::strip(token);
  if (token.size() > 3 && (token[0] == 'N' || token[0] == 'n') &&
      (token[1] == 'O' || token[1] == 'o') &&
      (token[2] == 'T' || token[2] == 't')) {
    negate = true;
    token = token.substr(3, token.size() - 3);
    token = FileParserUtils::strip(token);
  }

  std::unique_ptr<Atom> res;
  if (token[0] == '[') {
    // atom list:
    if (token.back() != ']') {
      std::ostringstream errout;
      errout << "Bad atom token '" << token << "' on line: " << line;
      throw FileParseException(errout.str());
    }
    token = token.substr(1, token.size() - 2);

    std::vector<std::string> splitToken;
    boost::split(splitToken, token, boost::is_any_of(","));

    for (std::vector<std::string>::const_iterator stIt = splitToken.begin();
         stIt != splitToken.end(); ++stIt) {
      std::string_view stoken = *stIt;
      std::string atSymb(FileParserUtils::strip(stoken));
      if (atSymb.empty()) {
        continue;
      }
      if (atSymb.size() == 2 && atSymb[1] >= 'A' && atSymb[1] <= 'Z') {
        atSymb[1] = static_cast<char>(tolower(atSymb[1]));
      }

      int atNum = PeriodicTable::getTable()->getAtomicNumber(atSymb);
      if (!res) {
        res.reset(new QueryAtom(atNum));
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
    bool isComplexQueryName =
        std::find(complexQueries.begin(), complexQueries.end(), token) !=
        complexQueries.end();
    if (isComplexQueryName || token == "R" ||
        (token[0] == 'R' && token >= "R0" && token <= "R99") || token == "R#" ||
        token == "*") {
      if (isComplexQueryName || token == "*") {
        res.reset(new QueryAtom(0));
        if (token == "*") {
          // according to the MDL spec, these match anything
          res->setQuery(makeAtomNullQuery());
        } else if (isComplexQueryName) {
          convertComplexNameToQuery(res.get(), token);
        }
        // queries have no implicit Hs:
        res->setNoImplicit(true);
      } else {
        res.reset(new Atom(1));
        res->setAtomicNum(0);
      }
      if (token[0] == 'R' && token >= "R0" && token <= "R99") {
        auto rlabel = token.substr(1, token.length() - 1);
        int rnumber;
        try {
          rnumber = boost::lexical_cast<int>(rlabel);
        } catch (boost::bad_lexical_cast &) {
          rnumber = -1;
        }
        if (rnumber >= 0) {
          res->setIsotope(rnumber);
        }
      }
      if (token[0] == 'R') {
        // we used to skip R# here because that really should be handled by an
        // RGP spec, but that turned out to not be permissive enough... <sigh>
        setRGPProps(token, res.get());
      }
    } else if (token == "D") {  // mol blocks support "D" and "T" as
                                // shorthand... handle that.
      res.reset(new Atom(1));
      res->setIsotope(2);
    } else if (token == "T") {  // mol blocks support "D" and "T" as
                                // shorthand... handle that.
      res.reset(new Atom(1));
      res->setIsotope(3);
    } else if (token == "Pol" || token == "Mod") {
      res.reset(new Atom(0));
      res->setProp(common_properties::dummyLabel, std::string(token));
    } else if (GenericGroups::genericMatchers.find(std::string(token)) !=
               GenericGroups::genericMatchers.end()) {
      res.reset(new QueryAtom(0));
      res->setProp(common_properties::atomLabel, std::string(token));
    } else {
      std::string tcopy(token);
      if (token.size() == 2 && token[1] >= 'A' && token[1] <= 'Z') {
        tcopy[1] = static_cast<char>(tolower(token[1]));
      }
      res.reset(new Atom(0));
      lookupAtomicNumber(res.get(), tcopy, strictParsing);
    }
  }

  POSTCONDITION(res, "no atom built");
  return res.release();
}

bool splitAssignToken(std::string_view token, std::string &prop,
                      std::string_view &val) {
  auto equalsLoc = token.find("=");
  if (equalsLoc == token.npos || equalsLoc != token.rfind("=")) {
    return false;
  }
  prop = token.substr(0, equalsLoc);
  boost::to_upper(prop);
  val = token.substr(equalsLoc + 1);
  return true;
}

template <class T>
void ParseV3000AtomProps(RWMol *mol, Atom *&atom, typename T::iterator &token,
                         const T &tokens, unsigned int &line,
                         bool strictParsing) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(atom, "bad atom");
  std::ostringstream errout;
  while (token != tokens.end()) {
    std::string prop;
    std::string_view val;
    if (!splitAssignToken(*token, prop, val)) {
      errout << "Invalid atom property: '" << *token << "' for atom "
             << atom->getIdx() + 1 << " on line " << line << std::endl;
      throw FileParseException(errout.str());
    }

    if (prop == "CHG") {
      auto charge = FileParserUtils::toInt(val);
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
      // Online examples seem to have integer (isotope) values and Marvin
      // won't even read something that has a float. We'll go with the int
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
      auto cfg = FileParserUtils::toInt(val);
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
        auto hcount = FileParserUtils::toInt(val);
        if (!atom->hasQuery()) {
          atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
        }
        if (hcount == -1) {
          hcount = 0;
        }
        if (hcount > 0) {
          ATOM_EQUALS_QUERY *oq = makeAtomImplicitHCountQuery(hcount);
          auto nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(
              hcount, oq->getDataFunc(),
              std::string("less_") + oq->getDescription());
          atom->expandQuery(nq);
          delete oq;
        } else {
          atom->expandQuery(makeAtomImplicitHCountQuery(0));
        }
      }
    } else if (prop == "UNSAT") {
      if (val == "1") {
        if (!atom->hasQuery()) {
          atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
        }
        atom->expandQuery(makeAtomUnsaturatedQuery());
      }
    } else if (prop == "RBCNT") {
      if (val != "0") {
        auto rbcount = FileParserUtils::toInt(val);
        if (!atom->hasQuery()) {
          atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
        }
        atom->setProp(common_properties::molRingBondCount, rbcount);
        if (rbcount == -1) {
          rbcount = 0;
        } else if (rbcount == -2) {
          // Ring bonds can only be counted during post processing
          mol->setProp(common_properties::_NeedsQueryScan, 1);
          rbcount = 0xDEADBEEF;
        } else if (rbcount > 4) {
          rbcount = 4;
        }
        atom->expandQuery(makeAtomRingBondCountQuery(rbcount));
      }
    } else if (prop == "VAL") {
      if (val != "0") {
        auto totval = FileParserUtils::toInt(val);
        atom->setProp(common_properties::molTotValence, totval);
      }
    } else if (prop == "RGROUPS") {
      ParseV3000RGroups(mol, atom, val, line);
      // FIX
    } else if (prop == "STBOX") {
      if (val != "0") {
        auto ival = FileParserUtils::toInt(val);
        atom->setProp(common_properties::molStereoCare, ival);
      }
    } else if (prop == "SUBST") {
      if (val != "0") {
        auto ival = FileParserUtils::toInt(val);
        atom->setProp(common_properties::molSubstCount, ival);
      }
    } else if (prop == "EXACHG") {
      if (val != "0") {
        auto ival = FileParserUtils::toInt(val);
        atom->setProp(common_properties::molRxnExactChange, ival);
      }
    } else if (prop == "INVRET") {
      if (val != "0") {
        auto ival = FileParserUtils::toInt(val);
        atom->setProp(common_properties::molInversionFlag, ival);
      }
    } else if (prop == "ATTCHPT") {
      if (val != "0") {
        auto ival = FileParserUtils::toInt(val);
        if (atom->hasProp(common_properties::molAttachPoint)) {
          errout << "Multiple ATTCHPT values for atom " << atom->getIdx() + 1
                 << " on line " << line;
          if (strictParsing) {
            throw FileParseException(errout.str());
          } else {
            BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
            errout.str(std::string());
          }
        } else {
          atom->setProp(common_properties::molAttachPoint, ival);
        }
      }
    } else if (prop == "ATTCHORD") {
      if (val != "0") {
        auto ival = FileParserUtils::toInt(val);
        atom->setProp(common_properties::molAttachOrder, ival);
      }
    } else if (prop == "CLASS") {
      atom->setProp(common_properties::molAtomClass, std::string(val));
    } else if (prop == "SEQID") {
      if (val != "0") {
        auto ival = FileParserUtils::toInt(val);
        atom->setProp(common_properties::molAtomSeqId, ival);
      }
    }
    ++token;
  }
}

void tokenizeV3000Line(std::string_view line,
                       std::vector<std::string_view> &tokens) {
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

bool calculate3dFlag(const RWMol &mol, const Conformer &conf,
                     bool chiralityPossible) {
  int marked3d = 0;
  if (mol.getPropIfPresent(common_properties::_3DConf, marked3d)) {
    mol.clearProp(common_properties::_3DConf);
  }

  bool nonzeroZ = hasNonZeroZCoords(conf);

  if (!nonzeroZ && marked3d == 1) {
    // If we have no Z coordinates, mark the structure 2D if we see any
    // 2D stereo markers, or stay as 3D if
    if (chiralityPossible) {
      BOOST_LOG(rdWarningLog)
          << "Warning: molecule is tagged as 3D, but all Z coords are zero and 2D stereo "
             "markers have been found, marking the mol as 2D."
          << std::endl;
      return false;
    }
    return true;
  } else if (marked3d == 0 && nonzeroZ) {
    BOOST_LOG(rdWarningLog)
        << "Warning: molecule is tagged as 2D, but at least one Z coordinate is not zero. "
           "Marking the mol as 3D."
        << std::endl;
    return true;
  }

  return nonzeroZ;
}

void ParseV3000AtomBlock(std::istream *inStream, unsigned int &line,
                         unsigned int nAtoms, RWMol *mol, Conformer *conf,
                         bool strictParsing) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(nAtoms > 0, "bad atom count");
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(conf, "bad conformer");
  std::vector<std::string> splitLine;

  auto inl = getV3000Line(inStream, line);
  std::string_view tempStr = inl;
  if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN ATOM") {
    std::ostringstream errout;
    errout << "BEGIN ATOM line not found on line " << line;
    throw FileParseException(errout.str());
  }
  for (unsigned int i = 0; i < nAtoms; ++i) {
    inl = getV3000Line(inStream, line);
    tempStr = inl;
    auto trimmed = FileParserUtils::strip(tempStr);

    std::vector<std::string_view> tokens;
    std::vector<std::string_view>::iterator token;

    tokenizeV3000Line(trimmed, tokens);
    token = tokens.begin();

    if (token == tokens.end()) {
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line" << line;
      throw FileParseException(errout.str());
    }
    unsigned int molIdx = 0;
    std::from_chars(token->data(), token->data() + token->size(), molIdx);

    // start with the symbol:
    ++token;
    if (token == tokens.end()) {
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line " << line;
      throw FileParseException(errout.str());
    }
    Atom *atom = ParseV3000AtomSymbol(*token, line, strictParsing);

    // now the position;
    RDGeom::Point3D pos;
    ++token;
    if (token == tokens.end()) {
      delete atom;
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line " << line;
      throw FileParseException(errout.str());
    }

    pos.x = atof(std::string(*token).c_str());
    ++token;
    if (token == tokens.end()) {
      delete atom;
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line " << line;
      throw FileParseException(errout.str());
    }
    pos.y = atof(std::string(*token).c_str());
    ++token;
    if (token == tokens.end()) {
      delete atom;
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line " << line;
      throw FileParseException(errout.str());
    }
    pos.z = atof(std::string(*token).c_str());
    // the map number:
    ++token;
    if (token == tokens.end()) {
      delete atom;
      std::ostringstream errout;
      errout << "Bad atom line : '" << tempStr << "' on line " << line;
      throw FileParseException(errout.str());
    }
    int mapNum = atoi(std::string(*token).c_str());
    if (mapNum > 0) {
      atom->setProp(common_properties::molAtomMapNumber, mapNum);
    }
    ++token;

    unsigned int aid = mol->addAtom(atom, false, true);

    // additional properties this may change the atom,
    // so be careful with it:
    ParseV3000AtomProps(mol, atom, token, tokens, line, strictParsing);

    mol->setAtomBookmark(atom, molIdx);
    conf->setAtomPos(aid, pos);
  }
  inl = getV3000Line(inStream, line);
  tempStr = inl;
  if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END ATOM") {
    std::ostringstream errout;
    errout << "END ATOM line not found on line " << line;
    throw FileParseException(errout.str());
  }
}

void ParseV3000BondBlock(std::istream *inStream, unsigned int &line,
                         unsigned int nBonds, RWMol *mol,
                         bool &chiralityPossible) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(nBonds > 0, "bad bond count");
  PRECONDITION(mol, "bad molecule");

  auto inl = getV3000Line(inStream, line);
  std::string_view tempStr = inl;
  if (tempStr.length() < 10 || tempStr.substr(0, 10) != "BEGIN BOND") {
    throw FileParseException("BEGIN BOND line not found");
  }
  for (unsigned int i = 0; i < nBonds; ++i) {
    inl = getV3000Line(inStream, line);
    tempStr = inl;
    tempStr = FileParserUtils::strip(tempStr);
    std::vector<std::string_view> splitLine;
    tokenizeV3000Line(tempStr, splitLine);
    if (splitLine.size() < 4) {
      std::ostringstream errout;
      errout << "bond line " << line << " is too short";
      throw FileParseException(errout.str());
    }
    Bond *bond;
    unsigned int bondIdx = 0;
    std::from_chars(splitLine[0].data(),
                    splitLine[0].data() + splitLine[0].size(), bondIdx);
    unsigned int bType = 0;
    std::from_chars(splitLine[1].data(),
                    splitLine[1].data() + splitLine[1].size(), bType);
    unsigned int a1Idx = 0;
    std::from_chars(splitLine[2].data(),
                    splitLine[2].data() + splitLine[2].size(), a1Idx);
    unsigned int a2Idx = 0;
    std::from_chars(splitLine[3].data(),
                    splitLine[3].data() + splitLine[3].size(), a2Idx);

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
      case 10:
        bond = new Bond(Bond::HYDROGEN);
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
        } else if (bType == 5) {
          bond->setQuery(makeSingleOrDoubleBondQuery());
          bond->setProp(common_properties::_MolFileBondQuery, 1);
        } else if (bType == 6) {
          bond->setQuery(makeSingleOrAromaticBondQuery());
          bond->setProp(common_properties::_MolFileBondQuery, 1);
        } else if (bType == 7) {
          bond->setQuery(makeDoubleOrAromaticBondQuery());
          bond->setProp(common_properties::_MolFileBondQuery, 1);
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
      std::string prop;
      std::string_view val;
      if (!splitAssignToken(splitLine[lPos], prop, val)) {
        errout << "bad bond property '" << splitLine[lPos] << "' on line "
               << line;
        throw FileParseException(errout.str());
      }
      if (prop == "CFG") {
        unsigned int cfg = 0;
        std::from_chars(val.data(), val.data() + val.size(), cfg);
        switch (cfg) {
          case 0:
            break;
          case 1:
            bond->setBondDir(Bond::BEGINWEDGE);
            chiralityPossible = true;
            break;
          case 2:
            if (bType == 1) {
              bond->setBondDir(Bond::UNKNOWN);
            } else if (bType == 2) {
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
        bond->setProp(common_properties::molReactStatus, reactStatus);
      } else if (prop == "STBOX") {
        bond->setProp(common_properties::molStereoCare, std::string(val));
      } else if (prop == "ENDPTS") {
        bond->setProp(common_properties::_MolFileBondEndPts, std::string(val));
      } else if (prop == "ATTACH") {
        bond->setProp(common_properties::_MolFileBondAttach, std::string(val));
      }
      ++lPos;
    }

    bond->setBeginAtomIdx(mol->getAtomWithBookmark(a1Idx)->getIdx());
    bond->setEndAtomIdx(mol->getAtomWithBookmark(a2Idx)->getIdx());
    mol->addBond(bond, true);
    mol->setBondBookmark(bond, bondIdx);

    // set the stereoCare property on the bond if it's not set already and
    // both the beginning and end atoms have it set:
    int care1 = 0;
    int care2 = 0;
    if (!bond->hasProp(common_properties::molStereoCare) &&
        mol->getAtomWithIdx(bond->getBeginAtomIdx())
            ->getPropIfPresent(common_properties::molStereoCare, care1) &&
        mol->getAtomWithIdx(bond->getEndAtomIdx())
            ->getPropIfPresent(common_properties::molStereoCare, care2)) {
      if (care1 == care2) {
        bond->setProp(common_properties::molStereoCare, care1);
      }
    }
  }
  inl = getV3000Line(inStream, line);
  tempStr = inl;
  if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END BOND") {
    std::ostringstream errout;
    errout << "END BOND line not found at line " << line;
    throw FileParseException(errout.str());
  }
}
// The documentation about MRV_COORDINATE_BOND_TYPE in
// https://docs.chemaxon.com/display/docs/chemaxon-specific-information-in-mdl-mol-files.md
// seems to be wrong: it says the only data field in this group contains the
// index for the coordinate atom. But behavior in Marvin Sketch seems to
// indicate that it references the bond index instead (see
// https://github.com/rdkit/rdkit/issues/4473)

void processMrvCoordinateBond(RWMol &mol, const SubstanceGroup &sg) {
  std::vector<std::string> dataFields;
  if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    if (dataFields.empty()) {
      BOOST_LOG(rdWarningLog)
          << "ignoring MRV_COORDINATE_BOND_TYPE SGroup without data fields."
          << std::endl;
      return;
    }

    auto coordinate_bond_idx =
        FileParserUtils::toUnsigned(dataFields[0], true) - 1;

    if (dataFields.size() > 1) {
      BOOST_LOG(rdWarningLog) << "ignoring extra data fields in "
                                 "MRV_COORDINATE_BOND_TYPE SGroup for bond "
                              << coordinate_bond_idx << '.' << std::endl;
    }

    Bond *old_bond = nullptr;
    try {
      old_bond = mol.getBondWithIdx(coordinate_bond_idx);
    } catch (const Invar::Invariant &) {
      BOOST_LOG(rdWarningLog)
          << "molecule does not contain a bond matching the "
             "MRV_COORDINATE_BOND_TYPE SGroup for bond "
          << coordinate_bond_idx << ", ignoring." << std::endl;
      return;
    }

    if (!old_bond || old_bond->getBondType() != Bond::BondType::UNSPECIFIED) {
      BOOST_LOG(rdWarningLog)
          << "MRV_COORDINATE_BOND_TYPE SGroup with value "
          << coordinate_bond_idx
          << " does not reference a query bond, ignoring." << std::endl;
      return;
    }

    Bond new_bond(Bond::BondType::DATIVE);
    auto preserveProps = true;
    auto keepSGroups = true;
    mol.replaceBond(coordinate_bond_idx, &new_bond, preserveProps, keepSGroups);
  }
}

void processSMARTSQ(RWMol &mol, const SubstanceGroup &sg) {
  std::string field;
  if (sg.getPropIfPresent("QUERYOP", field) && field != "=") {
    BOOST_LOG(rdWarningLog) << "unrecognized QUERYOP '" << field
                            << "' for SMARTSQ. Query ignored." << std::endl;
    return;
  }
  std::vector<std::string> dataFields;
  if (!sg.getPropIfPresent("DATAFIELDS", dataFields) || dataFields.empty()) {
    BOOST_LOG(rdWarningLog)
        << "empty FIELDDATA for SMARTSQ. Query ignored." << std::endl;
    return;
  }
  if (dataFields.size() > 1) {
    BOOST_LOG(rdWarningLog)
        << "multiple FIELDDATA values for SMARTSQ. Taking the first."
        << std::endl;
  }
  const std::string &sma = dataFields[0];
  if (sma.empty()) {
    BOOST_LOG(rdWarningLog)
        << "Skipping empty SMARTS value for SMARTSQ." << std::endl;
    return;
  }

  for (auto aidx : sg.getAtoms()) {
    auto at = mol.getAtomWithIdx(aidx);

    std::unique_ptr<RWMol> m;
    try {
      m.reset(SmartsToMol(sma));
    } catch (...) {
      // Is this ever used?
    }

    if (!m || !m->getNumAtoms()) {
      BOOST_LOG(rdWarningLog)
          << "SMARTS for SMARTSQ '" << sma
          << "' could not be parsed or has no atoms. Ignoring it." << std::endl;
      return;
    }

    if (!at->hasQuery()) {
      QueryAtom qAt(*at);
      int oidx = at->getIdx();
      mol.replaceAtom(oidx, &qAt);
      at = mol.getAtomWithIdx(oidx);
    }
    QueryAtom::QUERYATOM_QUERY *query = nullptr;
    if (m->getNumAtoms() == 1) {
      query = m->getAtomWithIdx(0)->getQuery()->copy();
    } else {
      query = new RecursiveStructureQuery(m.release());
    }
    at->setQuery(query);
    at->setProp(common_properties::MRV_SMA, sma);
    at->setProp(common_properties::_MolFileAtomQuery, 1);
  }
}

void processMrvImplicitH(RWMol &mol, const SubstanceGroup &sg) {
  std::vector<std::string> dataFields;
  if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    for (const auto &df : dataFields) {
      if (df.substr(0, 6) == "IMPL_H") {
        auto val = FileParserUtils::toInt(df.substr(6));
        for (auto atIdx : sg.getAtoms()) {
          if (atIdx < mol.getNumAtoms()) {
            // if the atom has aromatic bonds to it, then set the explicit
            // value, otherwise skip it.
            auto atom = mol.getAtomWithIdx(atIdx);
            bool hasAromaticBonds = false;
            for (auto bndI :
                 boost::make_iterator_range(mol.getAtomBonds(atom))) {
              auto bnd = (mol)[bndI];
              if (bnd->getIsAromatic() ||
                  bnd->getBondType() == Bond::AROMATIC) {
                hasAromaticBonds = true;
                break;
              }
            }
            if (hasAromaticBonds) {
              atom->setNumExplicitHs(val);
            } else {
              BOOST_LOG(rdWarningLog)
                  << "MRV_IMPLICIT_H SGroup on atom without aromatic "
                     "bonds, "
                  << atIdx << ", ignored." << std::endl;
            }
          } else {
            BOOST_LOG(rdWarningLog)
                << "bad atom index, " << atIdx
                << ", found in MRV_IMPLICIT_H SGroup. Ignoring it."
                << std::endl;
          }
        }
      }
    }
  }
}

void processZBO(RWMol &mol, const SubstanceGroup &sg) {
  for (auto bidx : sg.getBonds()) {
    auto bond = mol.getBondWithIdx(bidx);
    bond->setBondType(Bond::BondType::ZERO);
  }
}

void processZCH(RWMol &mol, const SubstanceGroup &sg) {
  RDUNUSED_PARAM(mol);
  std::vector<std::string> dataFields;
  if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    if (dataFields.empty()) {
      BOOST_LOG(rdWarningLog)
          << "ignoring ZCHG SGroup without data fields." << std::endl;
      return;
    }
    for (const auto &df : dataFields) {
      std::string trimmed = boost::trim_copy(df);
      std::vector<std::string> splitLine;
      boost::split(splitLine, trimmed, boost::is_any_of(";"),
                   boost::token_compress_off);
      const auto &aids = sg.getAtoms();
      if (splitLine.size() < aids.size()) {
        BOOST_LOG(rdWarningLog)
            << "DATAFIELDS in ZCH SGroup is shorter than the number of atoms in the SGroup. Ignoring it."
            << std::endl;
        continue;
      }
      for (auto i = 0u; i < aids.size(); ++i) {
        auto aid = aids[i];
        auto atom = mol.getAtomWithIdx(aid);
        auto val = 0;
        if (!splitLine[i].empty()) {
          val = FileParserUtils::toInt(splitLine[i]);
        }
        atom->setFormalCharge(val);
      }
    }
  }
}
void processHYD(RWMol &mol, const SubstanceGroup &sg) {
  std::vector<std::string> dataFields;
  if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    if (dataFields.empty()) {
      BOOST_LOG(rdWarningLog)
          << "ignoring HYD SGroup without data fields." << std::endl;
      return;
    }
    for (const auto &df : dataFields) {
      std::string trimmed = boost::trim_copy(df);
      std::vector<std::string> splitLine;
      boost::split(splitLine, trimmed, boost::is_any_of(";"),
                   boost::token_compress_off);
      const auto &aids = sg.getAtoms();
      if (splitLine.size() < aids.size()) {
        BOOST_LOG(rdWarningLog)
            << "DATAFIELDS in HYD SGroup is shorter than the number of atoms in the SGroup. Ignoring it."
            << std::endl;
        continue;
      }
      for (auto i = 0u; i < aids.size(); ++i) {
        auto aid = aids[i];
        auto atom = mol.getAtomWithIdx(aid);
        auto val = 0;
        if (!splitLine[i].empty()) {
          val = FileParserUtils::toInt(splitLine[i]);
        }
        atom->setProp("_ZBO_H", true);
        atom->setNumExplicitHs(val);
      }
    }
  }
}

// process (and remove) SGroups which modify the structure
// and which we can unambiguously apply
void processSGroups(RWMol *mol) {
  std::vector<unsigned int> sgsToRemove;
  unsigned int sgIdx = 0;
  for (auto &sg : getSubstanceGroups(*mol)) {
    if (sg.getProp<std::string>("TYPE") == "DAT") {
      std::string field;
      if (sg.getPropIfPresent("FIELDNAME", field)) {
        if (field == "MRV_COORDINATE_BOND_TYPE") {
          // V2000 support for coordinate bonds
          processMrvCoordinateBond(*mol, sg);
          sgsToRemove.push_back(sgIdx);
          continue;
        } else if (field == "MRV_IMPLICIT_H") {
          // CXN extension to specify implicit Hs, used for aromatic rings
          processMrvImplicitH(*mol, sg);
          sgsToRemove.push_back(sgIdx);
          continue;
        } else if (field == "ZBO") {
          // RDKit extension for zero-order bonds
          processZBO(*mol, sg);
          sgsToRemove.push_back(sgIdx);
          continue;
        } else if (field == "ZCH") {
          // RDKit extension for charge on atoms involved in zero-order bonds
          processZCH(*mol, sg);
          sgsToRemove.push_back(sgIdx);
          continue;
        } else if (field == "HYD") {
          // RDKit extension for hydrogen-count on atoms involved in
          // zero-order bonds
          processHYD(*mol, sg);
          sgsToRemove.push_back(sgIdx);
          continue;
        }
      }
      if (sg.getPropIfPresent("QUERYTYPE", field) &&
          (field == "SMARTSQ" || field == "SQ")) {
        processSMARTSQ(*mol, sg);
        sgsToRemove.push_back(sgIdx);
        continue;
      }
    }
    ++sgIdx;
  }
  // now remove the S groups we processed, we saved indices so do this in
  // backwards
  auto &sgs = getSubstanceGroups(*mol);
  for (auto it = sgsToRemove.rbegin(); it != sgsToRemove.rend(); ++it) {
    sgs.erase(sgs.begin() + *it);
  }
}

void ProcessMolProps(RWMol *mol) {
  PRECONDITION(mol, "no molecule");
  // we have to loop the ugly way because we may need to actually replace an
  // atom
  for (unsigned int aidx = 0; aidx < mol->getNumAtoms(); ++aidx) {
    auto atom = mol->getAtomWithIdx(aidx);
    int ival = 0;
    if (atom->getPropIfPresent(common_properties::molSubstCount, ival) &&
        ival != 0) {
      if (!atom->hasQuery()) {
        atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
      }
      bool gtQuery = false;
      if (ival == -1) {
        ival = 0;
      } else if (ival == -2) {
        // as drawn
        ival = atom->getDegree();
      } else if (ival >= 6) {
        // 6 or more
        gtQuery = true;
      }
      if (!gtQuery) {
        atom->expandQuery(makeAtomExplicitDegreeQuery(ival));
      } else {
        // create a temp query the normal way so that we can be sure to get
        // the description right
        std::unique_ptr<ATOM_EQUALS_QUERY> tmp{
            makeAtomExplicitDegreeQuery(ival)};
        atom->expandQuery(makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(
            ival, tmp->getDataFunc(),
            std::string("less_") + tmp->getDescription()));
      }
    }
    if (atom->getPropIfPresent(common_properties::molTotValence, ival) &&
        ival != 0 && !atom->hasProp("_ZBO_H")) {
      atom->setNoImplicit(true);
      if (ival == 15     // V2000
          || ival == -1  // v3000
      ) {
        atom->setNumExplicitHs(0);
      } else {
        if (atom->getExplicitValence() > ival) {
          BOOST_LOG(rdWarningLog)
              << "atom " << atom->getIdx() << " has specified valence (" << ival
              << ") smaller than the drawn valence "
              << atom->getExplicitValence() << "." << std::endl;
          atom->setNumExplicitHs(0);
        } else {
          atom->setNumExplicitHs(ival - atom->getExplicitValence());
        }
      }
    }
    atom->clearProp(common_properties::molTotValence);
  }
  processSGroups(mol);
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

  nAtoms = FileParserUtils::toUnsigned(splitLine[0]);
  nBonds = FileParserUtils::toUnsigned(splitLine[1]);
  conf = new Conformer(nAtoms);

  unsigned int nSgroups = 0, n3DConstraints = 0, chiralFlag = 0;

  if (splitLine.size() > 2) {
    nSgroups = FileParserUtils::toUnsigned(splitLine[2]);
  }
  if (splitLine.size() > 3) {
    n3DConstraints = FileParserUtils::toUnsigned(splitLine[3]);
  }
  if (splitLine.size() > 4) {
    chiralFlag = FileParserUtils::toUnsigned(splitLine[4]);
  }

  mol->setProp(common_properties::_MolFileChiralFlag, chiralFlag);

  if (nAtoms) {
    ParseV3000AtomBlock(inStream, line, nAtoms, mol, conf, strictParsing);
  }
  if (nBonds) {
    ParseV3000BondBlock(inStream, line, nBonds, mol, chiralityPossible);
  }

  tempStr = getV3000Line(inStream, line);
  // do link nodes:
  boost::to_upper(tempStr);
  while (tempStr.length() > 8 && tempStr.substr(0, 8) == "LINKNODE") {
    boost::to_upper(tempStr);
    // if the line has nothing on it we just ignore it
    if (tempStr.size() > 9) {
      std::string existing = "";
      if (mol->getPropIfPresent(common_properties::molFileLinkNodes,
                                existing)) {
        existing += "|";
      }
      existing += tempStr.substr(9);  // skip the "LINKNODE "
      mol->setProp(common_properties::molFileLinkNodes, existing);
    }
    tempStr = getV3000Line(inStream, line);
  }

  if (nSgroups) {
    boost::to_upper(tempStr);
    if (tempStr.length() < 12 || tempStr.substr(0, 12) != "BEGIN SGROUP") {
      std::ostringstream errout;
      errout << "BEGIN SGROUP line not found on line " << line;
      if (strictParsing) {
        throw FileParseException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
      }
    } else {
      tempStr =
          ParseV3000SGroupsBlock(inStream, line, nSgroups, mol, strictParsing);
      boost::to_upper(tempStr);
      if (tempStr.length() < 10 || tempStr.substr(0, 10) != "END SGROUP") {
        std::ostringstream errout;
        errout << "END SGROUP line not found on line " << line;
        if (strictParsing) {
          throw FileParseException(errout.str());
        } else {
          BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
        }
      } else {
        tempStr = getV3000Line(inStream, line);
      }
    }
  }

  while (tempStr.length() > 5 && tempStr.substr(0, 5) == "BEGIN") {
    if (tempStr.length() > 15 && tempStr.substr(6, 10) == "COLLECTION") {
      tempStr = parseEnhancedStereo(inStream, line, mol);
    } else {
      // skip blocks we don't know how to read
      BOOST_LOG(rdWarningLog) << "skipping block at line " << line << ": '"
                              << tempStr << "'" << std::endl;
      while (tempStr.length() < 3 || tempStr.substr(0, 3) != "END") {
        tempStr = getV3000Line(inStream, line);
      }
      tempStr = getV3000Line(inStream, line);
    }
  }

  if (n3DConstraints) {
    BOOST_LOG(rdWarningLog)
        << "3D constraint information in mol block ignored at line " << line
        << std::endl;
    boost::to_upper(tempStr);
    if (tempStr.length() < 11 || tempStr.substr(0, 11) != "BEGIN OBJ3D") {
      std::ostringstream errout;
      errout << "BEGIN OBJ3D line not found on line " << line;
      if (strictParsing) {
        throw FileParseException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
      }
    }
    for (unsigned int i = 0; i < n3DConstraints; ++i) {
      tempStr = getV3000Line(inStream, line);
    }
    tempStr = getV3000Line(inStream, line);
    boost::to_upper(tempStr);
    if (tempStr.length() < 9 || tempStr.substr(0, 9) != "END OBJ3D") {
      std::ostringstream errout;
      errout << "END OBJ3D line not found on line " << line;
      if (strictParsing) {
        throw FileParseException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
      }
    } else {
      tempStr = getV3000Line(inStream, line);
    }
  }

  boost::to_upper(tempStr);
  if (tempStr.length() < 8 || tempStr.substr(0, 8) != "END CTAB") {
    if (strictParsing) {
      throw FileParseException("END CTAB line not found");
    } else {
      BOOST_LOG(rdWarningLog) << "END CTAB line not found." << std::endl;
    }
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

  auto is3d = calculate3dFlag(*mol, *conf, chiralityPossible);
  conf->set3D(is3d);
  mol->addConformer(conf, true);
  conf = nullptr;

  return fileComplete;
}

bool ParseV2000CTAB(std::istream *inStream, unsigned int &line, RWMol *mol,
                    Conformer *&conf, bool &chiralityPossible,
                    unsigned int &nAtoms, unsigned int &nBonds,
                    bool strictParsing) {
  conf = new Conformer(nAtoms);

  if (nAtoms == 0) {
    conf->set3D(false);
  } else {
    ParseMolBlockAtoms(inStream, line, nAtoms, mol, conf, strictParsing);
  }
  ParseMolBlockBonds(inStream, line, nBonds, mol, chiralityPossible);

  auto is3d = calculate3dFlag(*mol, *conf, chiralityPossible);
  conf->set3D(is3d);
  mol->addConformer(conf, true);
  conf = nullptr;

  bool fileComplete =
      ParseMolBlockProperties(inStream, line, mol, strictParsing);
  return fileComplete;
}

void finishMolProcessing(
    RWMol *res, bool chiralityPossible,
    const RDKit::v2::FileParsers::MolFileParserParams &params) {
  if (!res) {
    return;
  }
  res->clearAllAtomBookmarks();
  res->clearAllBondBookmarks();

  if (params.expandAttachmentPoints) {
    MolOps::expandAttachmentPoints(*res);
  }

  // calculate explicit valence on each atom:
  for (auto atom : res->atoms()) {
    atom->calcExplicitValence(false);
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
  if (chiralityPossible || conf.is3D()) {
    if (!conf.is3D()) {
      DetectAtomStereoChemistry(*res, &conf);
    } else {
      res->updatePropertyCache(false);
      MolOps::assignChiralTypesFrom3D(*res, conf.getId(), true);
    }
  }

  Atropisomers::detectAtropisomerChirality(*res, &conf);

  // now that atom stereochem has been perceived, the wedging
  // information is no longer needed, so we clear
  // single bond dir flags:
  MolOps::clearSingleBondDirFlags(*res);

  if (params.sanitize) {
    if (params.removeHs) {
      // Bond stereo detection must happen before H removal, or
      // else we might be removing stereogenic H atoms in double
      // bonds (e.g. imines). But before we run stereo detection,
      // we need to run mol cleanup so don't have trouble with
      // e.g. nitro groups. Sadly, this a;; means we will find
      // run both cleanup and ring finding twice (a fast find
      // rings in bond stereo detection, and another in
      // sanitization's SSSR symmetrization).
      unsigned int failedOp = 0;
      MolOps::sanitizeMol(*res, failedOp, MolOps::SANITIZE_CLEANUP);
      MolOps::detectBondStereochemistry(*res);
      MolOps::removeHs(*res, false, false);
    } else {
      MolOps::sanitizeMol(*res);
      MolOps::detectBondStereochemistry(*res);
    }

    MolOps::assignStereochemistry(*res, true, true, true);
  } else {
    MolOps::detectBondStereochemistry(*res);
  }

  if (res->hasProp(common_properties::_NeedsQueryScan)) {
    res->clearProp(common_properties::_NeedsQueryScan);
    QueryOps::completeMolQueries(res);
  }
}
}  // namespace FileParserUtils

namespace v2 {
namespace FileParsers {
//------------------------------------------------
//
//  Read a molecule from a stream
//
//------------------------------------------------
std::unique_ptr<RWMol> MolFromMolDataStream(std::istream &inStream,
                                            unsigned int &line,
                                            const MolFileParserParams &params) {
  std::string tempStr;
  bool fileComplete = false;
  bool chiralityPossible = false;
  Utils::LocaleSwitcher ls;
  // mol name
  line++;
  tempStr = getLine(inStream);
  if (inStream.eof()) {
    return nullptr;
  }
  auto res = std::make_unique<RWMol>();
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
      res = nullptr;
    }
    std::ostringstream errout;
    errout << "Counts line too short: '" << tempStr << "' on line" << line;
    throw FileParseException(errout.str());
  }

  unsigned int spos = 0;
  // this needs to go into a try block because if the lexical_cast throws an
  // exception we want to catch throw a different exception
  try {
    nAtoms = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    spos = 3;
    nBonds = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    spos = 6;
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << tempStr.substr(spos, 3)
           << "' to unsigned int on line " << line;
    throw FileParseException(errout.str());
  }
  try {
    spos = 6;
    if (tempStr.size() >= 9) {
      nLists = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 12;
    if (tempStr.size() >= spos + 3) {
      chiralFlag = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 15;
    if (tempStr.size() >= spos + 3) {
      nsText = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 18;
    if (tempStr.size() >= spos + 3) {
      nRxnComponents =
          FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 21;
    if (tempStr.size() >= spos + 3) {
      nReactants = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 24;
    if (tempStr.size() >= spos + 3) {
      nProducts = FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

    spos = 27;
    if (tempStr.size() >= spos + 3) {
      nIntermediates =
          FileParserUtils::toUnsigned(tempStr.substr(spos, 3), true);
    }

  } catch (boost::bad_lexical_cast &) {
    // some SD files (such as some from NCI) lack all the extra information
    // on the header line, so ignore problems parsing there.
  }

  unsigned int ctabVersion = 2000;
  if (tempStr.size() > 35) {
    if (tempStr.size() < 39 || tempStr[34] != 'V') {
      std::ostringstream errout;
      errout << "CTAB version string invalid at line " << line;
      if (params.strictParsing) {
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
      if (params.strictParsing) {
        throw FileParseException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
      }
    }
  }

  res->setProp(common_properties::_MolFileChiralFlag, chiralFlag);

  Conformer *conf = nullptr;
  try {
    if (ctabVersion == 2000) {
      fileComplete = FileParserUtils::ParseV2000CTAB(
          &inStream, line, res.get(), conf, chiralityPossible, nAtoms, nBonds,
          params.strictParsing);
    } else {
      if (nAtoms != 0 || nBonds != 0) {
        std::ostringstream errout;
        errout << "V3000 mol blocks should have 0s in the initial counts line. "
                  "(line: "
               << line << ")";
        if (params.strictParsing) {
          throw FileParseException(errout.str());
        } else {
          BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
        }
      }
      fileComplete = FileParserUtils::ParseV3000CTAB(
          &inStream, line, res.get(), conf, chiralityPossible, nAtoms, nBonds,
          params.strictParsing);
    }
  } catch (MolFileUnhandledFeatureException &e) {
    // unhandled mol file feature, show an error
    res.reset();
    delete conf;
    conf = nullptr;
    BOOST_LOG(rdErrorLog) << " Unhandled CTAB feature: '" << e.what()
                          << "'. Molecule skipped." << std::endl;

    if (!inStream.eof()) {
      tempStr = getLine(inStream);
    }
    ++line;
    while (!inStream.eof() && !inStream.fail() &&
           tempStr.substr(0, 6) != "M  END" && tempStr.substr(0, 4) != "$$$$") {
      tempStr = getLine(inStream);
      ++line;
    }
    fileComplete = !inStream.eof() || tempStr.substr(0, 6) == "M  END" ||
                   tempStr.substr(0, 4) == "$$$$";
  } catch (FileParseException &e) {
    // catch our exceptions and throw them back after cleanup
    delete conf;
    conf = nullptr;
    throw e;
  }

  if (!fileComplete) {
    delete conf;
    conf = nullptr;
    std::ostringstream errout;
    errout
        << "Problems encountered parsing Mol data, M  END missing around line "
        << line;
    throw FileParseException(errout.str());
  }

  if (res) {
    FileParserUtils::finishMolProcessing(res.get(), chiralityPossible, params);
  }
  return res;
}

//------------------------------------------------
//
//  Read a molecule from a string
//
//------------------------------------------------
std::unique_ptr<RWMol> MolFromMolBlock(const std::string &molBlock,
                                       const MolFileParserParams &params) {
  std::istringstream inStream(molBlock);
  unsigned int line = 0;
  return MolFromMolDataStream(inStream, line, params);
}

//------------------------------------------------
//
//  Read a molecule from a file
//
//------------------------------------------------
std::unique_ptr<RWMol> MolFromMolFile(const std::string &fName,
                                      const MolFileParserParams &params) {
  std::ifstream inStream(fName.c_str());
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }
  if (!inStream.eof()) {
    unsigned int line = 0;
    return MolFromMolDataStream(inStream, line, params);
  } else {
    return std::unique_ptr<RWMol>();
  }
}
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
