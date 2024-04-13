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
  } catch (const std::out_of_range &) {
    std::ostringstream errout;
    errout << "SGroup line too short: '" << text << "' on line " << line;
    throw FileParseException(errout.str());
  }
  pos += len;
  return fieldValue;
}

unsigned int ParseSGroupIntField(bool &ok, bool strictParsing,
                                 const std::string &text, unsigned int line,
                                 unsigned int &pos, bool isFieldCounter) {
  ok = true;
  unsigned int res = 0;
  try {
    res = ParseSGroupIntField(text, line, pos, isFieldCounter);
  } catch (const std::exception &e) {
    if (strictParsing) {
      throw;
    } else {
      ok = false;
      BOOST_LOG(rdWarningLog) << e.what() << std::endl;
    }
  }
  return res;
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
  } catch (const std::out_of_range &) {
    std::ostringstream errout;
    errout << "SGroup line too short: '" << text << "' on line " << line;
    throw FileParseException(errout.str());
  }
  pos += len;
  return fieldValue;
}

double ParseSGroupDoubleField(bool &ok, bool strictParsing,
                              const std::string &text, unsigned int line,
                              unsigned int &pos) {
  ok = true;
  double res = 0.;
  try {
    res = ParseSGroupDoubleField(text, line, pos);
  } catch (const std::exception &e) {
    if (strictParsing) {
      throw;
    } else {
      ok = false;
      BOOST_LOG(rdWarningLog) << e.what() << std::endl;
    }
  }
  return res;
}

SubstanceGroup *FindSgIdx(IDX_TO_SGROUP_MAP &sGroupMap, int sgIdx,
                          unsigned int line) {
  auto sgIt = sGroupMap.find(sgIdx);
  if (sgIt == sGroupMap.end()) {
    BOOST_LOG(rdWarningLog) << "SGroup " << sgIdx << " referenced on line "
                            << line << " not found." << std::endl;
    return nullptr;
  }
  return &sgIt->second;
}

void ParseSGroupV2000STYLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  STY", "bad STY line");

  unsigned int pos = 6;
  bool ok;
  unsigned int nent =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    return;
  }

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup STY line too short: '" << text << "' on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      return;
    }

    unsigned int sequenceId =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      return;
    }

    std::string typ = text.substr(pos + 1, 3);
    if (SubstanceGroupChecks::isValidType(typ)) {
      auto sgroup = SubstanceGroup(mol, typ);
      sgroup.setProp<unsigned int>("index", sequenceId);
      sGroupMap.emplace(sequenceId, sgroup);
    } else {
      std::ostringstream errout;
      errout << "S group " << typ << " on line " << line;
      SGroupWarnOrThrow<MolFileUnhandledFeatureException>(strictParsing,
                                                          errout.str());
    }
    pos += 4;
  }
}

void ParseSGroupV2000VectorDataLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                                    const std::string &text, unsigned int line,
                                    bool strictParsing) {
  PRECONDITION(mol, "bad mol");

  std::string typ = text.substr(3, 3);

  void (SubstanceGroup::*sGroupAddIndexedElement)(const int) = nullptr;

  if (typ == "SAL") {
    sGroupAddIndexedElement = &SubstanceGroup::addAtomWithBookmark;
  } else if (typ == "SBL") {
    sGroupAddIndexedElement = &SubstanceGroup::addBondWithBookmark;
  } else if (typ == "SPA") {
    sGroupAddIndexedElement = &SubstanceGroup::addParentAtomWithBookmark;
  } else {
    std::ostringstream errout;
    errout << "Unsupported SGroup line '" << typ
           << "' passed to Vector Data parser ";
    throw FileParseException(errout.str());
  }

  unsigned int pos = 6;
  bool ok;
  unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
  if (!ok) {
    return;
  }
  SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
  if (!sgroup) {
    return;
  }
  unsigned int nent =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    sgroup->setIsValid(false);
    return;
  }

  for (unsigned int i = 0; i < nent; ++i) {
    if (text.size() < pos + 4) {
      std::ostringstream errout;
      errout << "SGroup line too short: '" << text << "' on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      sgroup->setIsValid(false);
      return;
    }
    unsigned int nbr = ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      sgroup->setIsValid(false);
      return;
    }
    try {
      (sgroup->*sGroupAddIndexedElement)(nbr);
    } catch (const std::exception &e) {
      SGroupWarnOrThrow<>(strictParsing, e.what());
      sgroup->setIsValid(false);
      return;
    }
  }
}

void ParseSGroupV2000SDILine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SDI", "bad SDI line");

  unsigned int pos = 6;
  bool ok;
  unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
  if (!ok) {
    return;
  }
  SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
  if (!sgroup) {
    return;
  }

  unsigned int nCoords =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    sgroup->setIsValid(false);
    return;
  }
  if (nCoords != 4) {
    std::ostringstream errout;
    errout << "Unexpected number of coordinates for SDI on line " << line;
    SGroupWarnOrThrow<>(strictParsing, errout.str());
    sgroup->setIsValid(false);
    return;
  }

  SubstanceGroup::Bracket bracket;
  for (unsigned int i = 0; i < 2; ++i) {
    double x = ParseSGroupDoubleField(ok, strictParsing, text, line, pos);
    if (!ok) {
      sgroup->setIsValid(false);
      return;
    }
    double y = ParseSGroupDoubleField(ok, strictParsing, text, line, pos);
    if (!ok) {
      sgroup->setIsValid(false);
      return;
    }
    double z = 0.;
    bracket[i] = RDGeom::Point3D(x, y, z);
  }
  bracket[2] = RDGeom::Point3D(0., 0., 0.);
  try {
    sgroup->addBracket(bracket);
  } catch (const std::exception &e) {
    SGroupWarnOrThrow<>(strictParsing, e.what());
    sgroup->setIsValid(false);
    return;
  }
}

void ParseSGroupV2000SSTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SST", "bad SST line");

  unsigned int pos = 6;
  bool ok;
  unsigned int nent =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    return;
  }

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SST line too short: '" << text << "' on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      return;
    }

    unsigned int sgIdx =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      return;
    }
    SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
    if (!sgroup) {
      return;
    };

    std::string subType = text.substr(++pos, 3);

    if (!SubstanceGroupChecks::isValidSubType(subType)) {
      std::ostringstream errout;
      errout << "Unsupported SGroup subtype '" << subType << "' on line "
             << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      sgroup->setIsValid(false);
      return;
    }

    sgroup->setProp("SUBTYPE", subType);
    pos += 3;
  }
}

void ParseSGroupV2000SMTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SMT", "bad SMT line");

  unsigned int pos = 6;
  bool ok;
  unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
  if (!ok) {
    return;
  }
  SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
  if (!sgroup) {
    return;
  }
  ++pos;

  if (pos >= text.length()) {
    std::ostringstream errout;
    errout << "SGroup line too short: '" << text << "' on line " << line;
    SGroupWarnOrThrow<>(strictParsing, errout.str());
    sgroup->setIsValid(false);
    return;
  }
  std::string label = text.substr(pos, text.length() - pos);

  if (sgroup->getProp<std::string>("TYPE") ==
      "MUL") {  // Case of multiple groups
    sgroup->setProp("MULT", label);

  } else {  // Case of abbreviation groups, but we might not have seen a SCL
            // line yet
    sgroup->setProp("LABEL", label);
  }
}

void ParseSGroupV2000SLBLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SLB", "bad SLB line");

  unsigned int pos = 6;
  bool ok;
  unsigned int nent =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    return;
  }

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SLB line too short: '" << text << "' on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      return;
    }

    unsigned int sgIdx =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      return;
    }
    SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
    if (!sgroup) {
      return;
    }
    unsigned int id = ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      sgroup->setIsValid(false);
      return;
    }
    if (id != 0 && !SubstanceGroupChecks::isSubstanceGroupIdFree(*mol, id)) {
      std::ostringstream errout;
      errout << "SGroup ID '" << id
             << "' is assigned to more than one SGroup, on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      sgroup->setIsValid(false);
      return;
    }

    sgroup->setProp<unsigned int>("ID", id);
  }
}

void ParseSGroupV2000SCNLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SCN", "bad SCN line");

  unsigned int pos = 6;
  bool ok;
  unsigned int nent =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    return;
  }

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 7) {
      std::ostringstream errout;
      errout << "SGroup SCN line too short: '" << text << "' on line " << line;
      errout << "\n needed: " << pos + 7 << " found: " << text.size();
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      return;
    }

    unsigned int sgIdx =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      return;
    }
    SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
    if (!sgroup) {
      return;
    }

    std::string connect = text.substr(++pos, 2);

    if (!SubstanceGroupChecks::isValidConnectType(connect)) {
      std::ostringstream errout;
      errout << "Unsupported SGroup connection type '" << connect
             << "' on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      sgroup->setIsValid(false);
      return;
    }

    sgroup->setProp("CONNECT", connect);
    pos += 3;
  }
}

void ParseSGroupV2000SDSLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 10) == "M  SDS EXP", "bad SDS line");

  unsigned int pos = 10;
  bool ok;
  unsigned int nent =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    return;
  }

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 4) {
      std::ostringstream errout;
      errout << "SGroup SDS line too short: '" << text << "' on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      return;
    }
    unsigned int sgIdx =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      return;
    }
    SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
    if (!sgroup) {
      return;
    }

    sgroup->setProp("ESTATE", "E");
  }
}

void ParseSGroupV2000SBVLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SBV", "bad SBV line");

  unsigned int pos = 6;
  bool ok;
  unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
  if (!ok) {
    return;
  }
  SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
  if (!sgroup) {
    return;
  }

  unsigned int bondMark =
      ParseSGroupIntField(ok, strictParsing, text, line, pos);
  if (!ok) {
    sgroup->setIsValid(false);
    return;
  }
  Bond *bond = mol->getUniqueBondWithBookmark(bondMark);

  RDGeom::Point3D vector;
  if (sgroup->getProp<std::string>("TYPE") == "SUP") {
    vector.x = ParseSGroupDoubleField(ok, strictParsing, text, line, pos);
    if (!ok) {
      sgroup->setIsValid(false);
      return;
    }
    vector.y = ParseSGroupDoubleField(ok, strictParsing, text, line, pos);
    if (!ok) {
      sgroup->setIsValid(false);
      return;
    }
    vector.z = 0.;
  }

  try {
    sgroup->addCState(bond->getIdx(), vector);
  } catch (const std::exception &e) {
    SGroupWarnOrThrow<>(strictParsing, e.what());
    sgroup->setIsValid(false);
    return;
  }
}

void ParseSGroupV2000SDTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SDT", "bad SDT line");

  unsigned int pos = 6;
  bool ok;
  unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
  if (!ok) {
    return;
  }
  SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
  if (!sgroup) {
    return;
  }

  std::string fieldName;
  std::string fieldType;
  std::string fieldInfo;
  std::string queryType;
  std::string queryOp;

  try {
    fieldName = text.substr(++pos, 30);
    boost::trim_right(fieldName);
    pos += 30;
    fieldType = text.substr(pos, 2);
    boost::trim_right(fieldType);
    pos += 2;
    fieldInfo = text.substr(pos, 20);
    boost::trim_right(fieldInfo);
    pos += 20;
    queryType = text.substr(pos, 2);
    boost::trim_right(queryType);
    pos += 2;
    queryOp = text.substr(pos, text.length() - pos);
    boost::trim_right(queryOp);
  } catch (const std::out_of_range &) {
    // all kinds of wild things out there... this insulates us from them without
    // making the code super complicated
  }

  // only add entries for the remaining properties if they aren't blank
  if (!fieldName.empty()) {
    sgroup->setProp("FIELDNAME", fieldName);
  }
  if (!fieldType.empty()) {
    sgroup->setProp("FIELDTYPE", fieldType);
  }
  if (!fieldInfo.empty()) {
    sgroup->setProp("FIELDINFO", fieldInfo);
  }
  if (!queryType.empty()) {
    sgroup->setProp("QUERYTYPE", queryType);
  }
  if (!queryOp.empty()) {
    sgroup->setProp("QUERYOP", queryOp);
  }
}

void ParseSGroupV2000SDDLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SDD", "bad SDD line");

  unsigned int pos = 6;
  bool ok;
  unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
  if (!ok) {
    return;
  }
  SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
  if (!sgroup) {
    return;
  }

  // Store the rest of the line as is.
  ++pos;
  if (pos < text.length()) {
    sgroup->setProp("FIELDDISP", text.substr(pos, text.length() - pos));
  }
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

  bool ok;
  unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
  if (!ok) {
    return;
  }
  SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
  if (!sgroup) {
    return;
  }

  if (lastDataSGroup != 0 && lastDataSGroup != sgIdx) {
    std::ostringstream errout;
    errout << "Found a Data Field not matching the SGroup of the last Data "
              "Field at line "
           << line;
    SGroupWarnOrThrow<>(strictParsing, errout.str());
    sgroup->setIsValid(false);
    return;
  } else if (lastDataSGroup == 0 && type == "SCD") {
    lastDataSGroup = sgIdx;
  } else if (type == "SED") {
    lastDataSGroup = 0;
  }

  // have we already seen an SDT line?
  if (!sgroup->hasProp("FIELDNAME")) {
    // one can read the docs and draw the conclusion that this is mandatory,
    // but it's also possible to interpret them the other way, and we know
    // that there are CTABs out there with empty fieldnames in SDT lines,
    // so let's just issue a warning and accept it.
    BOOST_LOG(rdWarningLog)
        << "Found a SCD/SED line with missing/empty SDT specification at line "
        << line << std::endl;
  }

  if (strictParsing) {
    if (type == "SCD" && counter > 2) {
      std::ostringstream errout;
      errout << "Found too many consecutive SCD lines, (#" << (counter + 1)
             << " at line " << line << ") for SGroup " << sgIdx;
      throw FileParseException(errout.str());
    }
  }

  if (pos + 1 < text.length()) {
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
}

void ParseSGroupV2000SPLLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SPL", "bad SPL line");

  unsigned int pos = 6;
  bool ok;
  unsigned int nent =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    return;
  }

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SPL line too short: '" << text << "' on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      return;
    }

    unsigned int sgIdx =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      return;
    }
    SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
    if (!sgroup) {
      return;
    }
    unsigned int parentIdx = ParseSGroupIntField(text, line, pos);

    sgroup->setProp<unsigned int>("PARENT", parentIdx);
  }
}

void ParseSGroupV2000SNCLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SNC", "bad SNC line");

  unsigned int pos = 6;
  bool ok;
  unsigned int nent =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    return;
  }

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SNC line too short: '" << text << "' on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      return;
    }

    unsigned int sgIdx =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      return;
    }
    SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
    if (!sgroup) {
      return;
    }

    unsigned int compno =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      sgroup->setIsValid(false);
      return;
    }
    if (compno > 256u) {
      std::ostringstream errout;
      errout << "SGroup SNC value over 256: '" << compno << "' on line "
             << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      sgroup->setIsValid(false);
      return;
    }
    sgroup->setProp<unsigned int>("COMPNO", compno);
  }
}

void ParseSGroupV2000SAPLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SAP", "bad SAP line");

  unsigned int pos = 6;
  bool ok;
  unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
  if (!ok) {
    return;
  }
  SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
  if (!sgroup) {
    return;
  }

  unsigned int nent =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    sgroup->setIsValid(false);
    return;
  }

  for (unsigned int ie = 0; ie < nent; ++ie) {
    int lvIdx = -1;
    if (text.size() < pos + 11) {
      std::ostringstream errout;
      errout << "SGroup SAP line too short: '" << text << "' on line " << line;
      if (strictParsing) {
        throw FileParseException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
        if (text.size() < pos + 4) {
          sgroup->setIsValid(false);
          return;
        }
        lvIdx = mol->getNumAtoms();
      }
    }

    std::string id = "  ";
    unsigned int aIdxMark =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      sgroup->setIsValid(false);
      return;
    }
    unsigned int aIdx = mol->getAtomWithBookmark(aIdxMark)->getIdx();

    if (lvIdx == -1) {
      unsigned int lvIdxMark =
          ParseSGroupIntField(ok, strictParsing, text, line, pos);
      if (!ok) {
        sgroup->setIsValid(false);
        return;
      }
      if (lvIdxMark != 0) {
        lvIdx = mol->getAtomWithBookmark(lvIdxMark)->getIdx();
      }
      if (text.size() >= pos + 3) {
        id = text.substr(pos + 1, 2);
        pos += 3;
      }
    }

    try {
      sgroup->addAttachPoint(aIdx, lvIdx, id);
    } catch (const std::exception &e) {
      SGroupWarnOrThrow<>(strictParsing, e.what());
      sgroup->setIsValid(false);
      return;
    }
  }
}

void ParseSGroupV2000SCLLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SCL", "bad SCL line");

  unsigned int pos = 6;
  bool ok;
  unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
  if (!ok) {
    return;
  }
  SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
  if (!sgroup) {
    return;
  }
  if (pos + 1 >= text.length()) {
    std::ostringstream errout;
    errout << "SGroup SCL line too short: '" << text << "' on line " << line;
    SGroupWarnOrThrow<>(strictParsing, errout.str());
    sgroup->setIsValid(false);
    return;
  }

  ++pos;
  sgroup->setProp("CLASS", text.substr(pos, text.length() - pos));
}

void ParseSGroupV2000SBTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(text.substr(0, 6) == "M  SBT", "bad SBT line");

  unsigned int pos = 6;
  bool ok;
  unsigned int nent =
      ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
  if (!ok) {
    return;
  }

  for (unsigned int ie = 0; ie < nent; ++ie) {
    if (text.size() < pos + 8) {
      std::ostringstream errout;
      errout << "SGroup SBT line too short: '" << text << "' on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      return;
    }

    unsigned int sgIdx =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      return;
    }
    SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
    if (!sgroup) {
      return;
    }
    unsigned int bracketType =
        ParseSGroupIntField(ok, strictParsing, text, line, pos);
    if (!ok) {
      sgroup->setIsValid(false);
      return;
    }

    if (bracketType == 0) {
      sgroup->setProp("BRKTYP", "BRACKET");
    } else if (bracketType == 1) {
      sgroup->setProp("BRKTYP", "PAREN");
    } else {
      std::ostringstream errout;
      errout << "Invalid SBT value '" << bracketType << "' on line " << line;
      SGroupWarnOrThrow<>(strictParsing, errout.str());
      sgroup->setIsValid(false);
      return;
    }
  }
}

/* ------------------ V3000 Utils  ------------------ */

template <class T>
std::vector<T> ParseV3000Array(std::stringstream &stream, int maxV,
                               bool strictParsing) {
  auto paren = stream.get();  // discard parentheses
  if (paren != '(') {
    BOOST_LOG(rdWarningLog)
        << "WARNING: first character of V3000 array is not '('" << std::endl;
  }

  unsigned int count = 0;
  stream >> count;
  std::vector<T> values;
  if (maxV >= 0 && count > static_cast<unsigned int>(maxV)) {
    SGroupWarnOrThrow(strictParsing, "invalid count value");
    return values;
  }

  values.reserve(count);
  T value;
  for (unsigned i = 0; i < count; ++i) {
    stream >> value;
    values.push_back(value);
  }
  paren = stream.get();  // discard parentheses
  if (paren != ')') {
    BOOST_LOG(rdWarningLog)
        << "WARNING: final character of V3000 array is not ')'" << std::endl;
  }
  return values;
}

// force instantiation of the versions of this that we use
template std::vector<unsigned int> ParseV3000Array(std::stringstream &stream,
                                                   int, bool);
template std::vector<int> ParseV3000Array(std::stringstream &stream, int, bool);

void ParseV3000CStateLabel(RWMol *mol, SubstanceGroup &sgroup,
                           std::stringstream &stream, unsigned int line,
                           bool strictParsing) {
  stream.get();  // discard parentheses

  unsigned int count;
  unsigned int bondMark;
  stream >> count >> bondMark;

  std::string type = sgroup.getProp<std::string>("TYPE");

  if ((type != "SUP" && count != 1) || (type == "SUP" && count != 4)) {
    std::ostringstream errout;
    errout << "Unexpected number of fields for CSTATE field on line " << line;
    SGroupWarnOrThrow<>(strictParsing, errout.str());
    sgroup.setIsValid(false);
    return;
  }

  Bond *bond = mol->getUniqueBondWithBookmark(bondMark);

  RDGeom::Point3D vector;
  if (type == "SUP") {
    stream >> vector.x >> vector.y >> vector.z;
  }
  try {
    sgroup.addCState(bond->getIdx(), vector);
  } catch (const std::exception &e) {
    SGroupWarnOrThrow<>(strictParsing, e.what());
    sgroup.setIsValid(false);
    return;
  }

  stream.get();  // discard final parentheses
}

void ParseV3000SAPLabel(RWMol *mol, SubstanceGroup &sgroup,
                        std::stringstream &stream, bool strictParsing) {
  stream.get();  // discard parentheses

  unsigned int count = 0;
  unsigned int aIdxMark = 0;
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

  try {
    sgroup.addAttachPoint(aIdx, lvIdx, sapIdStr);
  } catch (const std::exception &e) {
    SGroupWarnOrThrow<>(strictParsing, e.what());
    sgroup.setIsValid(false);
    return;
  }
}

std::string ParseV3000StringPropLabel(std::stringstream &stream) {
  std::string strValue;

  auto nextChar = stream.peek();
  if (nextChar == ' ') {
    // empty value, we peeked at the next field's separator
    return strValue;
  } else if (nextChar == '"') {
    // skip the opening quote:
    stream.get();

    // this is a bit gross because it's legal to include a \" in a value,
    // but the way that's done is by doubling it. So
    // FIELDINFO=""""
    // should assign the value \" to FIELDINFO
    char chr;
    while (stream.get(chr)) {
      if (chr == '"') {
        nextChar = stream.peek();

        // if the next element in the stream is a \" then we have a quoted \".
        // Otherwise we're done
        if (nextChar != '"') {
          break;
        } else {
          // skip the second \"
          stream.get();
        }
      }
      strValue += chr;
    }
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
                          unsigned int line, SubstanceGroup &sgroup, size_t,
                          RWMol *mol, bool strictParsing) {
  PRECONDITION(mol, "bad mol");
  // TODO: we could handle these in a more structured way
  try {
    if (label == "XBHEAD" || label == "XBCORR") {
      std::vector<unsigned int> bvect = ParseV3000Array<unsigned int>(
          lineStream, mol->getNumBonds(), strictParsing);
      std::transform(bvect.begin(), bvect.end(), bvect.begin(),
                     [](unsigned int v) -> unsigned int { return v - 1; });
      sgroup.setProp(label, bvect);
    } else if (label == "ATOMS") {
      for (auto atomIdx : ParseV3000Array<unsigned int>(
               lineStream, mol->getNumAtoms(), strictParsing)) {
        sgroup.addAtomWithBookmark(atomIdx);
      }
    } else if (label == "PATOMS") {
      for (auto patomIdx : ParseV3000Array<unsigned int>(
               lineStream, mol->getNumAtoms(), strictParsing)) {
        sgroup.addParentAtomWithBookmark(patomIdx);
      }
    } else if (label == "CBONDS" || label == "XBONDS") {
      for (auto bondIdx : ParseV3000Array<unsigned int>(
               lineStream, mol->getNumBonds(), strictParsing)) {
        sgroup.addBondWithBookmark(bondIdx);
      }
    } else if (label == "BRKXYZ") {
      auto coords = ParseV3000Array<double>(lineStream, 9, strictParsing);
      if (coords.size() != 9) {
        std::ostringstream errout;
        errout << "Unexpected number of coordinates for BRKXYZ on line "
               << line;
        throw FileParseException(errout.str());
      }

      SubstanceGroup::Bracket bracket;
      for (unsigned int i = 0; i < 3; ++i) {
        bracket[i] = RDGeom::Point3D(*(coords.begin() + (3 * i)),
                                     *(coords.begin() + (3 * i) + 1),
                                     *(coords.begin() + (3 * i) + 2));
      }
      sgroup.addBracket(bracket);
    } else if (label == "CSTATE") {
      ParseV3000CStateLabel(mol, sgroup, lineStream, line, strictParsing);
    } else if (label == "SAP") {
      ParseV3000SAPLabel(mol, sgroup, lineStream, strictParsing);
    } else if (label == "PARENT") {
      // Store relationship until all SGroups have been read
      unsigned int parentIdx;
      lineStream >> parentIdx;
      sgroup.setProp<unsigned int>("PARENT", parentIdx);
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

      if (label == "SUBTYPE" &&
          !SubstanceGroupChecks::isValidSubType(strValue)) {
        std::ostringstream errout;
        errout << "Unsupported SGroup subtype '" << strValue << "' on line "
               << line;
        throw FileParseException(errout.str());
      } else if (label == "CONNECT" &&
                 !SubstanceGroupChecks::isValidConnectType(strValue)) {
        std::ostringstream errout;
        errout << "Unsupported SGroup connection type '" << strValue
               << "' on line " << line;
        throw FileParseException(errout.str());
      }

      sgroup.setProp(label, strValue);
    }
  } catch (const std::exception &e) {
    SGroupWarnOrThrow<>(strictParsing, e.what());
    sgroup.setIsValid(false);
    return;
  }
}

std::string ParseV3000SGroupsBlock(std::istream *inStream, unsigned int line,
                                   unsigned int nSgroups, RWMol *mol,
                                   bool strictParsing) {
  PRECONDITION(inStream, "no stream");
  PRECONDITION(mol, "no molecule");
  unsigned int defaultLineNum = 0;
  std::string defaultString;

  // SGroups may be written in unsorted ID order, according to spec, so we will
  // temporarily store them in a map before adding them to the mol
  IDX_TO_SGROUP_MAP sGroupMap;

  std::unordered_map<std::string, std::stringstream> defaultLabels;

  auto tempStr = FileParserUtils::getV3000Line(inStream, line);

  // Store defaults
  if (tempStr.substr(0, 7) == "DEFAULT" && tempStr.length() > 8) {
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
    if (strictParsing && !SubstanceGroupChecks::isValidType(type)) {
      std::ostringstream errout;
      errout << "Unsupported SGroup type '" << type << "' on line " << line;
      throw MolFileUnhandledFeatureException(errout.str());
    }

    SubstanceGroup sgroup(mol, type);
    STR_VECT dataFields;

    sgroup.setProp<unsigned int>("index", sequenceId);
    if (externalId > 0) {
      if (!SubstanceGroupChecks::isSubstanceGroupIdFree(*mol, externalId)) {
        std::ostringstream errout;
        errout << "Existing SGroup ID '" << externalId
               << "' assigned to a second SGroup on line " << line;
        if (strictParsing) {
          throw FileParseException(errout.str());
        } else {
          BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
          sgroup.setIsValid(false);
        }
      }

      sgroup.setProp<unsigned int>("ID", externalId);
    }

    while (sgroup.getIsValid() && !lineStream.eof() && !lineStream.fail()) {
      char spacer;
      std::string label;

      lineStream.get(spacer);
      if (lineStream.gcount() == 0) {
        continue;
      } else if (spacer != ' ') {
        std::ostringstream errout;
        errout << "Found character '" << spacer
               << "' when expecting a separator (space) on line " << line;
        if (strictParsing) {
          throw FileParseException(errout.str());
        } else {
          BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
          sgroup.setIsValid(false);
          continue;
        }
      }

      std::getline(lineStream, label, '=');
      ParseV3000ParseLabel(label, lineStream, dataFields, line, sgroup,
                           nSgroups, mol, strictParsing);
      parsedLabels.insert(label);
    }

    // Process defaults
    lineStream.clear();
    lineStream.str(defaultString);
    while (sgroup.getIsValid() && !lineStream.eof() && !lineStream.fail()) {
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
        if (strictParsing) {
          throw FileParseException(errout.str());
        } else {
          BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
          sgroup.setIsValid(false);
          continue;
        }
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
          if (strictParsing) {
            throw FileParseException(errout.str());
          } else {
            BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
            sgroup.setIsValid(false);
            continue;
          }
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

  if (sGroupMap.size() != nSgroups) {
    std::ostringstream errout;
    errout << "Found " << sGroupMap.size() << " SGroups when " << nSgroups
           << " were expected." << std::endl;
    if (strictParsing) {
      throw FileParseException(errout.str());
    } else {
      BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    }
  }
  // SGroups successfully parsed, now add them to the molecule
  for (const auto &sg : sGroupMap) {
    if (sg.second.getIsValid()) {
      addSubstanceGroup(*mol, sg.second);
    } else {
      BOOST_LOG(rdWarningLog) << "SGroup " << sg.first
                              << " is invalid and will be ignored" << std::endl;
    }
  }
  return tempStr;
}

}  // namespace SGroupParsing
}  // namespace RDKit
