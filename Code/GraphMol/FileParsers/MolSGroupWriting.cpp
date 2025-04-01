//
//  Copyright (C) 2002-2018 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string/join.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include "FileParsers.h"
#include "MolSGroupWriting.h"

namespace RDKit {
namespace SGroupWriting {

/* ------------------ V2000 Utils  ------------------ */

std::string BuildV2000STYLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;

  const auto &sgroups = getSubstanceGroups(mol);
  for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    temp << FormatV2000IntField(1 + (sg - sgroups.begin()))
         << FormatV2000StringField(sg->getProp<std::string>("TYPE"), 3, true,
                                   true);
    if (++count == 8) {
      ret << "M  STY" << FormatV2000NumEntriesField(8) << temp.str() << "\n";
      temp.str("");
      count = 0;
    }
  }
  if (count) {
    ret << "M  STY" << FormatV2000NumEntriesField(count) << temp.str() << "\n";
  }

  return ret.str();
}

std::string BuildV2000StringPropLines(const unsigned int entriesPerLine,
                                      const ROMol &mol,
                                      const std::string &propName,
                                      const std::string &propCode,
                                      const unsigned int fieldWitdh) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  const auto &sgroups = getSubstanceGroups(mol);
  for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    std::string propValue;
    // Write field only if defined
    if (sg->getPropIfPresent(propName, propValue)) {
      temp << FormatV2000IntField(1 + (sg - sgroups.begin()))
           << FormatV2000StringField(propValue, fieldWitdh, true, true);
      if (++count == entriesPerLine) {
        ret << "M  " << propCode << FormatV2000NumEntriesField(entriesPerLine)
            << temp.str() << "\n";
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  " << propCode << FormatV2000NumEntriesField(count) << temp.str()
        << "\n";
  }

  return ret.str();
}

std::string BuildV2000SLBLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  const auto &sgroups = getSubstanceGroups(mol);
  for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    unsigned int id;
    // Write value if assigned, else 0
    if (sg->getPropIfPresent("ID", id)) {
      temp << FormatV2000IntField(1 + (sg - sgroups.begin()))
           << FormatV2000IntField(id);
      if (++count == 8) {
        ret << "M  SLB" << FormatV2000NumEntriesField(8) << temp.str() << "\n";
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SLB" << FormatV2000NumEntriesField(count) << temp.str() << "\n";
  }

  return ret.str();
}

std::string BuildV2000SDSLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  const auto &sgroups = getSubstanceGroups(mol);
  for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    // Write field only if defined
    std::string eState;
    if (sg->getPropIfPresent("ESTATE", eState) && eState == "E") {
      temp << FormatV2000IntField(1 + (sg - sgroups.begin()));
      if (++count == 15) {
        ret << "M  SDS EXP" << FormatV2000NumEntriesField(15) << temp.str()
            << "\n";
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SDS EXP" << FormatV2000NumEntriesField(count) << temp.str()
        << "\n";
  }

  return ret.str();
}

std::string BuildV2000SPLLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  const auto &sgroups = getSubstanceGroups(mol);
  for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    // Write field only if a parent is defined
    unsigned int parentIdx = -1;
    if (sg->getPropIfPresent("PARENT", parentIdx)) {
      temp << FormatV2000IntField(1 + (sg - sgroups.begin()))
           << FormatV2000IntField(parentIdx);
      if (++count == 8) {
        ret << "M  SPL" << FormatV2000NumEntriesField(8) << temp.str() << "\n";
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SPL" << FormatV2000NumEntriesField(count) << temp.str() << "\n";
  }

  return ret.str();
}

std::string BuildV2000SNCLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  const auto &sgroups = getSubstanceGroups(mol);
  for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    unsigned int compno;
    // Write field only if compno is set
    if (sg->getPropIfPresent("COMPNO", compno)) {
      temp << FormatV2000IntField(1 + (sg - sgroups.begin()))
           << FormatV2000IntField(compno);
      if (++count == 8) {
        ret << "M  SNC" << FormatV2000NumEntriesField(8) << temp.str() << "\n";
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SNC" << FormatV2000NumEntriesField(count) << temp.str() << "\n";
  }

  return ret.str();
}

std::string BuildV2000SBTLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  const auto &sgroups = getSubstanceGroups(mol);
  for (auto sg = sgroups.begin(); sg != sgroups.end(); ++sg) {
    std::string bracketType;
    if (sg->getPropIfPresent("BRKTYP", bracketType)) {
      unsigned int idx = 1 + (sg - sgroups.begin());
      if (bracketType == "BRACKET") {
        temp << FormatV2000IntField(idx) << FormatV2000IntField(0);
      } else if (bracketType == "PAREN") {
        temp << FormatV2000IntField(idx) << FormatV2000IntField(1);
      } else {
        std::ostringstream errout;
        errout << "Invalid BRKTYP value '" << bracketType << "' for SGroup "
               << idx;
        throw SubstanceGroupException(errout.str());
      }
      if (++count == 8) {
        ret << "M  SBT" << FormatV2000NumEntriesField(8) << temp.str() << "\n";
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SBT" << FormatV2000NumEntriesField(count) << temp.str() << "\n";
  }

  return ret.str();
}

template <class T>
std::string BuildV2000IdxVectorDataLines(const unsigned int entriesPerLine,
                                         const unsigned int sGroupId,
                                         const std::string &code,
                                         const T &idxVector) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (const auto &idx : idxVector) {
    temp << FormatV2000IntField(1 + idx);
    if (++count == entriesPerLine) {
      ret << "M  " << code << FormatV2000IntField(sGroupId)
          << FormatV2000NumEntriesField(entriesPerLine) << temp.str() << "\n";
      temp.str("");
      count = 0;
    }
  }
  if (count) {
    ret << "M  " << code << FormatV2000IntField(sGroupId)
        << FormatV2000NumEntriesField(count) << temp.str() << "\n";
  }
  return ret.str();
}

std::string BuildV2000SMTLine(const int idx, const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  std::string smtValue;
  if ((sgroup.getProp<std::string>("TYPE") == "MUL" &&
       sgroup.getPropIfPresent("MULT", smtValue)) ||
      sgroup.getPropIfPresent("LABEL", smtValue)) {
    ret << "M  SMT" << FormatV2000IntField(idx)
        << FormatV2000StringField(smtValue, 69, false, true) << "\n";
  }
  return ret.str();
}

std::string BuildV2000SDILine(const int idx, const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  const std::vector<SubstanceGroup::Bracket> brackets = sgroup.getBrackets();

  for (const auto &bracket : brackets) {
    ret << "M  SDI" << FormatV2000IntField(idx)
        << FormatV2000NumEntriesField(4);

    for (unsigned int iPoint = 0; iPoint < 2; ++iPoint) {
      ret << FormatV2000DoubleField(bracket.at(iPoint).x);
      ret << FormatV2000DoubleField(bracket.at(iPoint).y);
    }
    ret << "\n";
  }

  return ret.str();
}

std::string BuildV2000SBVLine(const int idx, const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  for (const auto &cstate : sgroup.getCStates()) {
    ret << "M  SBV" << FormatV2000IntField(idx)
        << FormatV2000IntField(cstate.bondIdx + 1);
    if (sgroup.getProp<std::string>("TYPE") == "SUP") {
      ret << FormatV2000DoubleField(cstate.vector.x);
      ret << FormatV2000DoubleField(cstate.vector.y);
    }
    ret << "\n";
  }

  return ret.str();
}

std::string BuildV2000SDTLine(const int idx, const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  std::string sdtValue;
  if (sgroup.getPropIfPresent("FIELDNAME", sdtValue)) {
    ret << "M  SDT" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sdtValue, 30, true, true);

    if (sgroup.getPropIfPresent("FIELDTYPE", sdtValue)) {
      ret << FormatV2000StringField(sdtValue, 2, true, false);
    } else {
      ret << " T";
    }

    if (sgroup.getPropIfPresent("FIELDINFO", sdtValue)) {
      ret << FormatV2000StringField(sdtValue, 20, true, false);
    }

    if (sgroup.getPropIfPresent("QUERYTYPE", sdtValue)) {
      ret << FormatV2000StringField(sdtValue, 2, true, false);
    }
    if (sgroup.getPropIfPresent("QUERYOP", sdtValue)) {
      ret << FormatV2000StringField(sdtValue, 15, true, false);
    }

    ret << "\n";
  }
  return ret.str();
}

std::string BuildV2000SDDLine(const int idx, const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  std::string sddValue;
  if (sgroup.getPropIfPresent("FIELDDISP", sddValue)) {
    ret << "M  SDD" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sddValue, 69, false, true);
    ret << "\n";
  }

  return ret.str();
}

std::string BuildV2000SCDSEDLines(const int idx, const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  STR_VECT dataFields;
  if (sgroup.getPropIfPresent("DATAFIELDS", dataFields)) {
    for (const auto &data : dataFields) {
      unsigned int length = data.size();
      if (length > 200) {
        std::ostringstream errout;
        errout << "Data field '" << data << "' in SGroup " << idx
               << " is longer than limit of 200 characters.";
        throw SubstanceGroupException(errout.str());
      }
      unsigned int start = 0;
      unsigned int end = 69;
      for (; length > end; start += 69, end += 69) {
        std::string dataChunk = data.substr(start, end);
        ret << "M  SCD" << FormatV2000IntField(idx)
            << FormatV2000StringField(dataChunk, 69, true, true) << "\n";
      }
      std::string dataChunk = data.substr(start, end);
      ret << "M  SED" << FormatV2000IntField(idx)
          << FormatV2000StringField(dataChunk, 69, false, true) << "\n";
    }
  }

  return ret.str();
}

std::string BuildV2000SAPLines(const int idx, const SubstanceGroup &sgroup) {
  std::ostringstream ret;
  std::ostringstream temp;

  const std::vector<SubstanceGroup::AttachPoint> saps =
      sgroup.getAttachPoints();

  unsigned int count = 0;
  unsigned int entriesPerLine = 6;
  for (const auto &sap : saps) {
    temp << FormatV2000IntField(1 + sap.aIdx);

    // lvIdx == -1 will turn into 0, which is right (see spec)
    temp << FormatV2000IntField(1 + sap.lvIdx);

    temp << FormatV2000StringField(sap.id, 2, false, true);
    if (++count == entriesPerLine) {
      ret << "M  SAP" << FormatV2000IntField(idx)
          << FormatV2000IntField(entriesPerLine) << temp.str() << "\n";
      temp.str("");
      count = 0;
    }
  }
  if (count) {
    ret << "M  SAP" << FormatV2000IntField(idx)
        << FormatV2000NumEntriesField(count) << temp.str() << "\n";
  }
  return ret.str();
}

std::string BuildV2000SCLLine(const int idx, const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  std::string sclValue;
  if (sgroup.getPropIfPresent("CLASS", sclValue)) {
    ret << "M  SCL" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sclValue, 69, false, true);
    ret << "\n";
  }

  return ret.str();
}

const std::string GetMolFileSGroupInfo(const RWMol &mol) {
  std::ostringstream ret;

  // multiple group per line properties
  ret << BuildV2000STYLines(mol);
  ret << BuildV2000SLBLines(mol);
  ret << BuildV2000StringPropLines(8, mol, "SUBTYPE", "SST", 3);
  ret << BuildV2000StringPropLines(8, mol, "CONNECT", "SCN", 3);
  ret << BuildV2000SDSLines(mol);
  ret << BuildV2000SPLLines(mol);
  ret << BuildV2000SNCLines(mol);
  ret << BuildV2000SBTLines(mol);

  // single group per line properties
  unsigned int idx = 0;
  for (const auto &sgroup : getSubstanceGroups(mol)) {
    ++idx;
    ret << BuildV2000IdxVectorDataLines(15, idx, "SAL", sgroup.getAtoms());

    // SPA line must always come somewhere after the SAL line.
    ret << BuildV2000IdxVectorDataLines(15, idx, "SPA",
                                        sgroup.getParentAtoms());

    ret << BuildV2000IdxVectorDataLines(15, idx, "SBL", sgroup.getBonds());
    // Write CRS line -- CRS still not supported
    ret << BuildV2000SDILine(idx, sgroup);

    ret << BuildV2000SMTLine(idx, sgroup);
    ret << BuildV2000SBVLine(idx, sgroup);

    ret << BuildV2000SDTLine(idx, sgroup);
    ret << BuildV2000SDDLine(idx, sgroup);
    // SCD/SED must come after SDT
    ret << BuildV2000SCDSEDLines(idx, sgroup);

    ret << BuildV2000SAPLines(idx, sgroup);
    ret << BuildV2000SCLLine(idx, sgroup);
  }
  return ret.str();
}

/* ------------------ V3000 Utils  ------------------ */

template <class T>
std::string BuildV3000IdxVectorDataBlock(const std::string &key,
                                         const std::vector<T> &dataVector) {
  return BuildV3000IdxVectorDataBlock(key, dataVector.begin(),
                                      dataVector.end());
}

template <class Iterator>
std::string BuildV3000IdxVectorDataBlock(const std::string &key,
                                         const Iterator &dataVectorBegin,
                                         const Iterator &dataVectorEnd) {
  std::ostringstream ret;
  size_t size = dataVectorEnd - dataVectorBegin;
  if (size) {
    ret << ' ' << key << "=(" << size;
    for (auto itr = dataVectorBegin; itr < dataVectorEnd; ++itr) {
      ret << ' ' << 1 + *itr;
    }
    ret << ')';
  }
  return ret.str();
}

std::string BuildV3000BondsBlock(const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  auto isXBond = [&sgroup](unsigned int bondIdx) {
    return SubstanceGroup::BondType::XBOND == sgroup.getBondType(bondIdx);
  };

  auto bonds = sgroup.getBonds();
  auto first_cbond = std::stable_partition(bonds.begin(), bonds.end(), isXBond);

  ret << BuildV3000IdxVectorDataBlock("XBONDS", bonds.begin(), first_cbond);
  ret << BuildV3000IdxVectorDataBlock("CBONDS", first_cbond, bonds.end());

  if (sgroup.hasProp("XBHEAD")) {
    auto v = sgroup.getProp<std::vector<unsigned int>>("XBHEAD");
    ret << BuildV3000IdxVectorDataBlock("XBHEAD", v.begin(), v.end());
  }
  if (sgroup.hasProp("XBCORR")) {
    auto v = sgroup.getProp<std::vector<unsigned int>>("XBCORR");
    ret << BuildV3000IdxVectorDataBlock("XBCORR", v.begin(), v.end());
  }

  return ret.str();
}

std::string FormatV3000StringPropertyBlock(const std::string &prop,
                                           const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  std::string propValue;
  if (sgroup.getPropIfPresent(prop, propValue)) {
    if (!propValue.empty()) {
      ret << ' ' << prop << '=';
      // CTAB spec says: "Strings that contain blank spaces or start with left
      // parenthesis or double quote, must be surrounded by double quotes A
      // double quote can be entered literally by doubling it."
      // However, BIOVIA Draw 2020 doesn't correctly parse values like
      // foo"" or foo(bar) but does fine with "foo""" and "foo(bar)"
      // and both BIOVIA Draw and Marvin Sketch happily ignore the theoretically
      // extra quotes.
      bool needsQuotes = propValue.find(' ') != std::string::npos ||
                         propValue.find('"') != std::string::npos ||
                         propValue.find('(') != std::string::npos;
      if (needsQuotes) {
        ret << "\"";
      }

      for (auto chr : propValue) {
        ret << chr;
        // double quotes need to be doubled on output:
        if (chr == '"') {
          ret << chr;
        }
      }

      if (needsQuotes) {
        ret << "\"";
      }
    }
  }

  return ret.str();
}

std::string FormatV3000ParentBlock(const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  unsigned int parentIdx = -1;
  if (sgroup.getPropIfPresent("PARENT", parentIdx)) {
    ret << " PARENT=" << parentIdx;
  }

  return ret.str();
}

std::string FormatV3000CompNoBlock(const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  unsigned int compno;

  if (sgroup.getPropIfPresent("COMPNO", compno)) {
    ret << " COMPNO=" << compno;
  }

  return ret.str();
}

std::string FormatV3000BracketBlock(
    const std::vector<SubstanceGroup::Bracket> brackets) {
  std::ostringstream ret;

  for (const auto &bracket : brackets) {
    ret << " BRKXYZ=(9";
    for (unsigned int iPoint = 0; iPoint < 2; ++iPoint) {
      ret << ' ' << FormatV3000DoubleField(bracket.at(iPoint).x);
      ret << ' ' << FormatV3000DoubleField(bracket.at(iPoint).y);
      ret << " 0";  // z coordinate is 0 by format specification
    }
    ret << " 0 0 0";  // 3rd point is 0 by format specification
    ret << ")";
  }

  return ret.str();
}

std::string FormatV3000FieldDataBlock(const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  STR_VECT dataFields;
  if (sgroup.getPropIfPresent("DATAFIELDS", dataFields)) {
    for (const auto &data : dataFields) {
      ret << " FIELDDATA=\"" << data << "\"";
    }
  }

  return ret.str();
}

std::string FormatV3000CStateBlock(const SubstanceGroup &sgroup) {
  std::ostringstream ret;

  for (const auto &cstate : sgroup.getCStates()) {
    unsigned int xbondIdx = 1 + cstate.bondIdx;
    ret << " CSTATE=(";
    if (sgroup.getProp<std::string>("TYPE") == "SUP") {
      ret << "4 " << xbondIdx;
      ret << ' ' << FormatV3000DoubleField(cstate.vector.x);
      ret << ' ' << FormatV3000DoubleField(cstate.vector.y);
      ret << " 0";
    } else {
      ret << "1 " << xbondIdx;
    }
    ret << ")";
  }

  return ret.str();
}

std::string FormatV3000AttachPointBlock(
    const std::vector<SubstanceGroup::AttachPoint> &saps) {
  std::ostringstream ret;

  for (const auto &sap : saps) {
    ret << " SAP=(3 " << (1 + sap.aIdx);

    if (sap.lvIdx != -1 && sap.aIdx == static_cast<unsigned int>(sap.lvIdx)) {
      ret << " aidx";
    } else {
      // lvIdx == -1 will turn into 0, which is right (see spec)
      ret << ' ' << (1 + sap.lvIdx);
    }

    ret << ' ' << sap.id << ")";
  }

  return ret.str();
}

namespace {
void addBlockToSGroupString(std::string block, std::string &currentLine,
                            std::ostringstream &os) {
  if (block.empty()) {
    return;
  }
  if (currentLine.length() + block.length() < 78) {
    currentLine += block;
  } else {
    os << currentLine << " -\n";
    unsigned int length = block.size();
    unsigned int start = 0;
    while (length - start >= 73) {
      os << "M  V30";
      if (start) {
        os << ' ';
      }
      os << block.substr(start, 72);
      start += 72;
      if (start < length) {
        // need to write more, so add another "-"
        os << "-\n";
      }
    }
    if (start < length) {
      currentLine =
          "M  V30" + std::string(start ? " " : "") + block.substr(start, 73);
    }
  }
}
}  // namespace

//! Write a SGroup line. The order of the labels is defined in the spec.
const std::string GetV3000MolFileSGroupLines(const unsigned int idx,
                                             const SubstanceGroup &sgroup) {
  std::ostringstream os;

  unsigned int id = 0;
  sgroup.getPropIfPresent("ID", id);

  std::string currLine = (boost::format("M  V30 %d %s %d") % idx %
                          sgroup.getProp<std::string>("TYPE") % id)
                             .str();
  addBlockToSGroupString(
      BuildV3000IdxVectorDataBlock("ATOMS", sgroup.getAtoms()), currLine, os);
  // also writes XBHEAD and XBCORR
  addBlockToSGroupString(BuildV3000BondsBlock(sgroup), currLine, os);
  addBlockToSGroupString(
      BuildV3000IdxVectorDataBlock("PATOMS", sgroup.getParentAtoms()), currLine,
      os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("SUBTYPE", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("MULT", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("CONNECT", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000ParentBlock(sgroup), currLine, os);
  addBlockToSGroupString(FormatV3000CompNoBlock(sgroup), currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("LABEL", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000BracketBlock(sgroup.getBrackets()),
                         currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("ESTATE", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000CStateBlock(sgroup), currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("FIELDNAME", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("FIELDINFO", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("FIELDDISP", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("QUERYTYPE", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("QUERYOP", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000FieldDataBlock(sgroup), currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("CLASS", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000AttachPointBlock(sgroup.getAttachPoints()),
                         currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("BRKTYP", sgroup),
                         currLine, os);
  addBlockToSGroupString(FormatV3000StringPropertyBlock("SEQID", sgroup),
                         currLine, os);
  std::string res;
  if (!currLine.empty() && currLine != "M  V30") {
    os << currLine << "\n";
    res = os.str();
  } else {
    // there's an extra " -" at the end that we need to remove
    res = os.str();
    res = res.substr(0, res.length() - 3) + "\n";
  }
  return res;
}

}  // namespace SGroupWriting
}  // namespace RDKit
