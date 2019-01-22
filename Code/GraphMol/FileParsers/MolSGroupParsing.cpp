//
//  Copyright (C) 2002-2018 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FileParsers.h"
#include "FileParserUtils.h"
#include "MolSGroupParsing.h"

#include <RDGeneral/FileParseException.h>

namespace RDKit {
namespace SGroupParsing {

/* ------------------ V2000 Utils  ------------------ */

unsigned int ParseSGroupIntField(const std::string &text, unsigned int line,
                                 unsigned int &pos, bool isFieldCounter) {
  ++pos;  // Account for separation space
  unsigned int fieldValue;
  size_t len = 3 - isFieldCounter;  // field counters are smaller
  try {
    fieldValue = FileParserUtils::toInt(text.substr(pos, len));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(pos, len) << "' to int on line "
           << line;
    throw FileParseException(errout.str());
  }
  pos += len;
  return fieldValue;
}

double ParseSGroupDoubleField(const std::string &text, unsigned int line,
                              unsigned int &pos) {
  size_t len = 10;
  double fieldValue;
  try {
    fieldValue = FileParserUtils::toDouble(text.substr(pos, len));
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Cannot convert '" << text.substr(pos, len)
           << "' to double on line " << line;
    throw FileParseException(errout.str());
  }
  pos += len;
  return fieldValue;
}

void ParseSGroupV2000STYLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  STY", "bad STY line");

  unsigned int pos = 6;
  unsigned int nent = ParseSGroupIntField(text, line, pos, true);

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup STY line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }

    int nbr = ParseSGroupIntField(text, line, pos);

    std::string typ = text.substr(pos + 1, 3);
    if (SGroupChecks::isValidType(typ)) {
      sGroupMap.emplace(nbr, SGroup(mol, typ));
    } else {
      std::ostringstream errout;
      errout << "S group " << typ << " on line " << line;
      throw MolFileUnhandledFeatureException(errout.str());
    }
    pos += 4;
  }
}

void ParseSGroupV2000VectorDataLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                                    const std::string &text,
                                    unsigned int line) {
  PRECONDITION(mol, "bad mol");

  std::string typ = text.substr(3, 3);

  void (SGroup::*sGroupAddIndexedElement)(const int) = nullptr;

  if (typ == "SAL") {
    sGroupAddIndexedElement = &SGroup::addAtomWithBookmark;
  } else if (typ == "SBL") {
    sGroupAddIndexedElement = &SGroup::addBondWithBookmark;
  } else if (typ == "SPA") {
    sGroupAddIndexedElement = &SGroup::addParentAtomWithBookmark;
  } else {
    std::ostringstream errout;
    errout << "Unsupported SGroup line '" << typ
           << "' passed to Vector Data parser ";
    throw FileParseException(errout.str());
  }

  unsigned int pos = 6;
  unsigned int sgIdx = ParseSGroupIntField(text, line, pos);
  unsigned int nent = ParseSGroupIntField(text, line, pos, true);

  for (unsigned int i = 0; i < nent; ++i) {
    if (text.size() < pos + 4) {
      std::ostringstream errout;
      errout << "SGroup line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }

    int nbr = ParseSGroupIntField(text, line, pos);
    (sGroupMap.at(sgIdx).*sGroupAddIndexedElement)(nbr);
  }
}

void ParseSGroupV2000SDILine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SDI", "bad SDI line");

  unsigned int pos = 6;
  unsigned int sgIdx = ParseSGroupIntField(text, line, pos);

  unsigned int nCoords = ParseSGroupIntField(text, line, pos, true);
  if (nCoords != 4) {
    std::ostringstream errout;
    errout << "Unexpected number of coordinates for SDI on line " << line;
    throw FileParseException(errout.str());
  }

  SGroup::Bracket bracket;
  for (unsigned int i = 0; i < 2; ++i) {
    double x = ParseSGroupDoubleField(text, line, pos);
    double y = ParseSGroupDoubleField(text, line, pos);
    double z = 0.;
    bracket[i] = RDGeom::Point3D(x, y, z);
  }
  bracket[2] = RDGeom::Point3D(0., 0., 0.);
  sGroupMap.at(sgIdx).addBracket(bracket);
}

void ParseSGroupV2000SSTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SST", "bad SST line");

  unsigned int pos = 6;
  unsigned int nent = ParseSGroupIntField(text, line, pos, true);

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SST line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }

    unsigned int sgIdx = ParseSGroupIntField(text, line, pos);

    std::string subType = text.substr(++pos, 3);

    if (!SGroupChecks::isValidSubType(subType)) {
      std::ostringstream errout;
      errout << "Unsupported SGroup subtype '" << subType << "' on line "
             << line;
      throw FileParseException(errout.str());
    }

    sGroupMap.at(sgIdx).setProp("SUBTYPE", subType);
    pos += 3;
  }
}

void ParseSGroupV2000SMTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SMT", "bad SMT line");

  unsigned int pos = 6;
  unsigned int sgIdx = ParseSGroupIntField(text, line, pos);
  auto &sgroup = sGroupMap.at(sgIdx);
  ++pos;

  std::string label = text.substr(pos, text.length() - pos);

  if (sgroup.getProp<std::string>("TYPE") ==
      "MUL") {  // Case of multiple groups
    sgroup.setProp("MULT", label);

  } else {  // Case of abbreviation groups, but we might not have seen a SCL
            // line yet
    sgroup.setProp("LABEL", label);
  }
}

void ParseSGroupV2000SLBLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SLB", "bad SLB line");

  unsigned int pos = 6;
  unsigned int nent = ParseSGroupIntField(text, line, pos, true);

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SLB line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }

    unsigned int sgIdx = ParseSGroupIntField(text, line, pos);
    unsigned int id = ParseSGroupIntField(text, line, pos);

    if (id != 0 && !SGroupChecks::isSGroupIdFree(*mol, id)) {
      std::ostringstream errout;
      errout << "SGroup ID '" << id
             << "' is assigned to more than onge SGroup, on line " << line;
      throw FileParseException(errout.str());
    }

    sGroupMap.at(sgIdx).setProp<unsigned int>("ID", id);
  }
}

void ParseSGroupV2000SCNLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SCN", "bad SCN line");

  unsigned int pos = 6;
  unsigned int nent = ParseSGroupIntField(text, line, pos, true);

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SCN line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }

    unsigned int sgIdx = ParseSGroupIntField(text, line, pos);

    std::string connect = text.substr(++pos, 2);

    if (!SGroupChecks::isValidConnectType(connect)) {
      std::ostringstream errout;
      errout << "Unsupported SGroup connection type '" << connect
             << "' on line " << line;
      throw FileParseException(errout.str());
    }

    sGroupMap.at(sgIdx).setProp("CONNECT", connect);
    pos += 3;
  }
}

void ParseSGroupV2000SDSLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 10) == "M  SDS EXP", "bad SDS line");

  unsigned int pos = 10;
  unsigned int nent = ParseSGroupIntField(text, line, pos, true);

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 4) {
      std::ostringstream errout;
      errout << "SGroup SDS line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }
    unsigned int sgIdx = ParseSGroupIntField(text, line, pos);

    sGroupMap.at(sgIdx).setProp("ESTATE", "E");
  }
}

void ParseSGroupV2000SBVLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SBV", "bad SBV line");

  unsigned int pos = 6;

  unsigned int sgIdx = ParseSGroupIntField(text, line, pos);
  auto &sgroup = sGroupMap.at(sgIdx);

  unsigned int bondMark = ParseSGroupIntField(text, line, pos);
  Bond *bond = mol->getUniqueBondWithBookmark(bondMark);

  RDGeom::Point3D vector;
  if (sgroup.getProp<std::string>("TYPE") == "SUP") {
    vector.x = ParseSGroupDoubleField(text, line, pos);
    vector.y = ParseSGroupDoubleField(text, line, pos);
    vector.z = 0.;
  }

  sgroup.addCState(bond->getIdx(), vector);
}

void ParseSGroupV2000SDTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SDT", "bad SDT line");

  unsigned int pos = 6;
  unsigned int sgIdx = ParseSGroupIntField(text, line, pos);
  auto &sgroup = sGroupMap.at(sgIdx);

  std::string fieldName = text.substr(++pos, 30);
  boost::trim_right(fieldName);
  pos += 30;
  std::string fieldType = text.substr(pos, 2);
  boost::trim_right(fieldType);
  pos += 2;
  std::string fieldInfo = text.substr(pos, 20);
  boost::trim_right(fieldInfo);
  pos += 20;
  std::string queryType = text.substr(pos, 2);
  boost::trim_right(queryType);
  pos += 2;
  std::string queryOp = text.substr(pos, text.length() - pos);
  boost::trim_right(queryOp);

  sgroup.setProp("FIELDNAME", fieldName);
  sgroup.setProp("FIELDTYPE", fieldType);
  sgroup.setProp("FIELDINFO", fieldInfo);
  sgroup.setProp("QUERYTYPE", queryType);
  sgroup.setProp("QUERYOP", queryOp);
}

void ParseSGroupV2000SDDLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SDD", "bad SDD line");

  unsigned int pos = 6;
  unsigned int sgIdx = ParseSGroupIntField(text, line, pos);

  // Store the rest of the line as is.
  ++pos;
  sGroupMap.at(sgIdx).setProp("FIELDDISP",
                              text.substr(pos, text.length() - pos));
}

void ParseSGroupV2000SCDSEDLine(IDX_TO_SGROUP_MAP &sGroupMap,
                                IDX_TO_STR_VECT_MAP &dataFieldsMap, RWMol *mol,
                                const std::string &text, unsigned int line,
                                bool strictParsing, unsigned int &counter,
                                unsigned int &lastDataSGroup,
                                std::ostringstream &currentDataField) {
  PRECONDITION(mol, "bad mol");

  unsigned int pos = 3;
  std::string type = text.substr(pos, 3);
  pos += 3;

  unsigned int sgIdx = ParseSGroupIntField(text, line, pos);

  if (lastDataSGroup != 0 && lastDataSGroup != sgIdx) {
    std::ostringstream errout;
    errout << "Found a Data Field not matching the the SGroup of the last Data "
              "Field at line "
           << line;
    throw FileParseException(errout.str());
  } else if (lastDataSGroup == 0 && type == "SCD") {
    lastDataSGroup = sgIdx;
  } else if (type == "SED") {
    lastDataSGroup = 0;
  }

  auto &sgroup = sGroupMap.at(sgIdx);

  // this group must already have seen an SDT line
  if (!sgroup.hasProp("FIELDNAME")) {
    std::ostringstream errout;
    errout << "Found a SCD line without a previous SDT specification at line "
           << line;
    throw FileParseException(errout.str());
  }

  if (strictParsing) {
    if (type == "SCD" && counter > 2) {
      std::ostringstream errout;
      errout << "Found too many consecutive SCD lines, (#" << (counter + 1)
             << " at line " << line << ") for SGroup " << sgIdx;
      throw FileParseException(errout.str());
    }
  }

  currentDataField << text.substr(++pos, 69);

  if (type == "SED") {
    std::string trimmedData = boost::trim_right_copy(currentDataField.str());
    dataFieldsMap[sgIdx].push_back(trimmedData.substr(0, 200));
    currentDataField.str("");
    counter = 0;
  } else {
    ++counter;
  }
}

void ParseSGroupV2000SPLLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SPL", "bad SPL line");

  unsigned int pos = 6;
  unsigned int nent = ParseSGroupIntField(text, line, pos, true);

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SPL line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }

    unsigned int sgIdx = ParseSGroupIntField(text, line, pos);
    unsigned int parentIdx = ParseSGroupIntField(text, line, pos);

    // both parentIdx and size are offset by +1, so it's fine
    if (parentIdx > sGroupMap.size()) {
      std::ostringstream errout;
      errout << "SGroup SPL line contains wrong parent SGroup (" << parentIdx
             << ") for SGroup " << sgIdx << "on line " << line;
      throw FileParseException(errout.str());
    }

    sGroupMap.at(sgIdx).setProp<unsigned int>("PARENT", parentIdx - 1);
  }
}

void ParseSGroupV2000SNCLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SNC", "bad SNC line");

  unsigned int pos = 6;
  unsigned int nent = ParseSGroupIntField(text, line, pos, true);

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SNC line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }

    unsigned int sgIdx = ParseSGroupIntField(text, line, pos);

    unsigned int compno = ParseSGroupIntField(text, line, pos);
    if (compno > 256u) {
      std::ostringstream errout;
      errout << "SGroup SNC value over 256: '" << compno << "' on line "
             << line;
      throw FileParseException(errout.str());
    }
    sGroupMap.at(sgIdx).setProp<unsigned int>("COMPNO", compno);
  }
}

// PXA does not have a V3000 equivalent
void ParseSGroupV2000PXALine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  PXA", "bad PXA line");

  unsigned int pos = 6;

  unsigned int sgIdx = ParseSGroupIntField(text, line, pos);

  ++pos;
  sGroupMap.at(sgIdx).setProp("PXA", text.substr(pos, text.length() - pos));
}

void ParseSGroupV2000SAPLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SAP", "bad SAP line");

  unsigned int pos = 6;

  unsigned int sgIdx = ParseSGroupIntField(text, line, pos);

  unsigned int nent = ParseSGroupIntField(text, line, pos, true);

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SAP line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }

    unsigned int aIdxMark = ParseSGroupIntField(text, line, pos);
    unsigned int lvIdxMark = ParseSGroupIntField(text, line, pos);

    unsigned int aIdx = mol->getAtomWithBookmark(aIdxMark)->getIdx();
    int lvIdx = -1;

    if (lvIdxMark != 0) {
      lvIdx = mol->getAtomWithBookmark(lvIdxMark)->getIdx();
    }

    sGroupMap.at(sgIdx).addAttachPoint(aIdx, lvIdx, text.substr(++pos, 2));
  }
}

void ParseSGroupV2000SCLLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SCL", "bad SCL line");

  unsigned int pos = 6;

  unsigned int sgIdx = ParseSGroupIntField(text, line, pos);

  ++pos;
  sGroupMap.at(sgIdx).setProp("CLASS", text.substr(pos, text.length() - pos));
}

void ParseSGroupV2000SBTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SBT", "bad SBT line");

  unsigned int pos = 6;
  unsigned int nent = ParseSGroupIntField(text, line, pos, true);

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SBT line too short: '" << text << "' on line " << line;
      throw FileParseException(errout.str());
    }

    unsigned int sgIdx = ParseSGroupIntField(text, line, pos);
    unsigned int bracketType = ParseSGroupIntField(text, line, pos);

    auto &sgroup = sGroupMap.at(sgIdx);

    if (bracketType == 0) {
      sgroup.setProp("BRKTYP", "BRACKET");
    } else if (bracketType == 1) {
      sgroup.setProp("BRKTYP", "PAREN");
    } else {
      std::ostringstream errout;
      errout << "Invalid SBT value '" << bracketType << "' on line " << line;
      throw FileParseException(errout.str());
    }
  }
}

/* ------------------ V3000 Utils  ------------------ */

template <class T>
std::vector<T> ParseV3000Array(std::stringstream &stream) {
  stream.get();  // discard parentheses

  unsigned int count;
  stream >> count;
  std::vector<T> values;
  values.reserve(count);
  T value;
  for (unsigned i = 0; i < count; ++i) {
    stream >> value;
    values.push_back(value);
  }
  stream.get();  // discard parentheses
  return values;
}

void ParseV3000CStateLabel(RWMol *mol, SGroup &sgroup,
                           std::stringstream &stream, unsigned int line) {
  stream.get();  // discard parentheses

  unsigned int count;
  unsigned int bondMark;
  stream >> count >> bondMark;

  std::string type = sgroup.getProp<std::string>("TYPE");

  if ((type != "SUP" && count != 1) || (type == "SUP" && count != 4)) {
    std::ostringstream errout;
    errout << "Unexpected number of fields for CSTATE field on line " << line;
    throw FileParseException(errout.str());
  }

  Bond *bond = mol->getUniqueBondWithBookmark(bondMark);

  RDGeom::Point3D vector;
  if (type == "SUP") {
    stream >> vector.x >> vector.y >> vector.z;
  }
  sgroup.addCState(bond->getIdx(), vector);

  stream.get();  // discard final parentheses
}

void ParseV3000SAPLabel(RWMol *mol, SGroup &sgroup, std::stringstream &stream) {
  stream.get();  // discard parentheses

  unsigned int count;
  unsigned int aIdxMark;
  std::string lvIdxStr;  // In V3000 this may be a string
  std::string sapIdStr;
  stream >> count >> aIdxMark >> lvIdxStr >> sapIdStr;

  // remove final parentheses that gets parsed into sapIdStr
  sapIdStr.pop_back();

  unsigned int aIdx = mol->getAtomWithBookmark(aIdxMark)->getIdx();
  int lvIdx = -1;

  boost::to_upper(lvIdxStr);
  if (lvIdxStr == "AIDX") {
    lvIdx = aIdx;
  } else {
    unsigned int lvIdxTmp = FileParserUtils::toInt(lvIdxStr);
    if (lvIdxTmp > 0) {
      lvIdx = mol->getAtomWithBookmark(lvIdxTmp)->getIdx();
    }
  }

  sgroup.addAttachPoint(aIdx, lvIdx, sapIdStr);
}

std::string ParseV3000StringPropLabel(std::stringstream &stream) {
  std::string strValue;

  // TODO: this should be improved to be able to handle escaped quotes

  auto nextChar = stream.peek();
  if (nextChar == '"') {
    stream.get();
    std::getline(stream, strValue, '"');
  } else if (nextChar == '\'') {
    std::getline(stream, strValue, '\'');
  } else {
    stream >> strValue;
  }

  boost::trim_right(strValue);
  return strValue;
}

void ParseV3000ParseLabel(const std::string &label,
                          std::stringstream &lineStream, STR_VECT &dataFields,
                          unsigned int &line, SGroup &sgroup, size_t nSgroups,
                          RWMol *mol, bool &strictParsing) {
  // TODO: remove this once we find out how to handle XBHEAD & XBCORR
  if (label == "XBHEAD" || label == "XBCORR") {
    std::ostringstream errout;
    errout << "XBHEAD or XBCORR labels (found on line " << line
           << ") are not yet supported";
    throw FileParseException(errout.str());
  }

  if (label == "ATOMS") {
    for (auto atomIdx : ParseV3000Array<unsigned int>(lineStream)) {
      sgroup.addAtomWithBookmark(atomIdx);
    }
  } else if (label == "PATOMS") {
    for (auto patomIdx : ParseV3000Array<unsigned int>(lineStream)) {
      sgroup.addParentAtomWithBookmark(patomIdx);
    }
  } else if (label == "CBONDS" || label == "XBONDS") {
    for (auto bondIdx : ParseV3000Array<unsigned int>(lineStream)) {
      sgroup.addBondWithBookmark(bondIdx);
    }
  } else if (label == "BRKXYZ") {
    auto coords = ParseV3000Array<double>(lineStream);
    if (coords.size() != 9) {
      std::ostringstream errout;
      errout << "Unexpected number of coordinates for BRKXYZ on line " << line;
      throw FileParseException(errout.str());
    }

    SGroup::Bracket bracket;
    for (unsigned int i = 0; i < 3; ++i) {
      bracket[i] = RDGeom::Point3D(*(coords.begin() + (3 * i)),
                                   *(coords.begin() + (3 * i) + 1),
                                   *(coords.begin() + (3 * i) + 2));
    }
    sgroup.addBracket(bracket);
  } else if (label == "CSTATE") {
    ParseV3000CStateLabel(mol, sgroup, lineStream, line);
  } else if (label == "SAP") {
    ParseV3000SAPLabel(mol, sgroup, lineStream);
  } else if (label == "PARENT") {
    // Store relationship until all SGroups have been read
    unsigned int parentIdx;
    lineStream >> parentIdx;

    // both parentIdx and nSGroups are offset by +1, so it's fine
    if (parentIdx > nSgroups) {
      std::ostringstream errout;
      errout << "Wrong parent SGroup '" << parentIdx << "' on line " << line;
      throw FileParseException(errout.str());
    }

    sgroup.setProp<unsigned int>("PARENT", parentIdx - 1);
  } else if (label == "COMPNO") {
    unsigned int compno;
    lineStream >> compno;
    if (compno > 256u) {
      std::ostringstream errout;
      errout << "SGroup SNC value over 256: '" << compno << "' on line "
             << line;
      throw FileParseException(errout.str());
    }
    sgroup.setProp<unsigned int>("COMPNO", compno);
  } else if (label == "FIELDDATA") {
    auto strValue = ParseV3000StringPropLabel(lineStream);
    if (strictParsing) {
      strValue = strValue.substr(0, 200);
    }
    dataFields.push_back(strValue);

  } else {
    // Parse string props
    auto strValue = ParseV3000StringPropLabel(lineStream);

    if (label == "SUBTYPE" && !SGroupChecks::isValidSubType(strValue)) {
      std::ostringstream errout;
      errout << "Unsupported SGroup subtype '" << strValue << "' on line "
             << line;
      throw FileParseException(errout.str());
    } else if (label == "CONNECT" &&
               !SGroupChecks::isValidConnectType(strValue)) {
      std::ostringstream errout;
      errout << "Unsupported SGroup connection type '" << strValue
             << "' on line " << line;
      throw FileParseException(errout.str());
    }

    sgroup.setProp(label, strValue);
  }
}

void ParseV3000SGroupsBlock(std::istream *inStream, unsigned int &line,
                            unsigned int nSgroups, RWMol *mol,
                            bool &strictParsing) {
  std::string tempStr = FileParserUtils::getV3000Line(inStream, line);
  boost::to_upper(tempStr);
  if (tempStr.length() < 12 || tempStr.substr(0, 12) != "BEGIN SGROUP") {
    std::ostringstream errout;
    errout << "BEGIN SGROUP line not found on line " << line;
    throw FileParseException(errout.str());
  }

  unsigned int defaultLineNum = 0;
  std::string defaultString;

  // SGroups may be written in unsorted ID order, according to spec, so we will
  // temporarily store them in a map before adding them to the mol
  IDX_TO_SGROUP_MAP sGroupMap;

  std::unordered_map<std::string, std::stringstream> defaultLabels;

  tempStr = FileParserUtils::getV3000Line(inStream, line);

  // Store defaults
  if (tempStr.substr(0, 7) == "DEFAULT") {
    defaultString = tempStr.substr(7);
    defaultLineNum = line;
    boost::trim_right(defaultString);
    tempStr = FileParserUtils::getV3000Line(inStream, line);
    boost::trim_right(tempStr);
  }

  for (unsigned int si = 0; si < nSgroups; ++si) {
    unsigned int sequenceId;
    unsigned int externalId;
    std::string type;

    std::stringstream lineStream(tempStr);
    lineStream >> sequenceId;
    lineStream >> type;
    lineStream >> externalId;

    std::set<std::string> parsedLabels;

    if (strictParsing && !SGroupChecks::isValidType(type)) {
      std::ostringstream errout;
      errout << "Unsupported SGroup type '" << type << "' on line " << line;
      throw MolFileUnhandledFeatureException(errout.str());
    }

    SGroup sgroup(mol, type);
    STR_VECT dataFields;

    if (externalId > 0) {
      if (!SGroupChecks::isSGroupIdFree(*mol, externalId)) {
        std::ostringstream errout;
        errout << "Existing SGroup ID '" << externalId
               << "' assigned to a second SGroup on line " << line;
        throw FileParseException(errout.str());
      }

      sgroup.setProp<unsigned int>("ID", externalId);
    }

    while (!lineStream.eof()) {
      char spacer;
      std::string label;

      lineStream.get(spacer);
      if (lineStream.gcount() == 0) {
        continue;
      } else if (spacer != ' ') {
        std::ostringstream errout;
        errout << "Found character '" << spacer
               << "' when expecting a separator (space) on line " << line;
        throw FileParseException(errout.str());
      }

      std::getline(lineStream, label, '=');
      ParseV3000ParseLabel(label, lineStream, dataFields, line, sgroup,
                           nSgroups, mol, strictParsing);
      parsedLabels.insert(label);
    }

    // Process defaults
    lineStream.clear();
    lineStream.str(defaultString);
    while (!lineStream.eof()) {
      char spacer;
      std::string label;

      lineStream.get(spacer);
      if (lineStream.gcount() == 0) {
        continue;
      } else if (spacer != ' ') {
        std::ostringstream errout;
        errout << "Found character '" << spacer
               << "' when expecting a separator (space) in DEFAULTS on line "
               << defaultLineNum;
        throw FileParseException(errout.str());
      }

      std::getline(lineStream, label, '=');
      if (std::find(parsedLabels.begin(), parsedLabels.end(), label) ==
          parsedLabels.end()) {
        ParseV3000ParseLabel(label, lineStream, dataFields, defaultLineNum,
                             sgroup, nSgroups, mol, strictParsing);
      } else {
        spacer = lineStream.peek();
        if (spacer == ' ') {
          std::ostringstream errout;
          errout << "Found unexpected whitespace at DEFAULT label " << label;
          throw FileParseException(errout.str());
        } else if (spacer == '(') {
          std::getline(lineStream, label, ')');
          lineStream.get(spacer);
        } else if (spacer == '"') {
          lineStream.get(spacer);
          std::getline(lineStream, label, '"');
        } else {
          std::getline(lineStream, label, ' ');
          lineStream.putback(' ');
        }
      }
    }

    sgroup.setProp("DATAFIELDS", dataFields);
    sGroupMap.emplace(sequenceId, sgroup);

    tempStr = FileParserUtils::getV3000Line(inStream, line);
    boost::trim_right(tempStr);
  }

  boost::to_upper(tempStr);
  if (tempStr.length() < 10 || tempStr.substr(0, 10) != "END SGROUP") {
    std::ostringstream errout;
    errout << "END SGROUP line not found on line " << line;
    throw FileParseException(errout.str());
  }

  // SGroups successfully parsed, now add them to the molecule sorted by index
  for (const auto &sg : sGroupMap) {
    addSGroup(*mol, sg.second);
  }
}

}  // namespace SGroupParsing
}  // namespace RDKit
