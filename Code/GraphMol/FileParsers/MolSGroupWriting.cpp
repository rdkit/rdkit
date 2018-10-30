//
//  Copyright (C) 2002-2018 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/algorithm/string/join.hpp>

#include "FileParsers.h"
#include "MolSGroupWriting.h"

namespace RDKit {
namespace SGroupWriting {

/* ------------------ V2000 Utils  ------------------ */

std::string BuildV2000STYLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = 1 + (sGroupItr - mol.beginSGroups());
    temp << FormatV2000IntField(idx)
         << FormatV2000StringField((*sGroupItr)->getType(), 3, true, true);
    if (++count == 8) {
      ret << "M  STY" << FormatV2000NumEntriesField(8) << temp.str()
          << std::endl;
      temp.str("");
      count = 0;
    }
  }
  if (count) {
    ret << "M  STY" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
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
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    if ((*sGroupItr)->hasProp(propName)) {  // Write field only if defined
      temp << FormatV2000IntField(idx)
           << FormatV2000StringField((*sGroupItr)->getProp(propName),
                                     fieldWitdh, true, true);
      if (++count == entriesPerLine) {
        ret << "M  " << propCode << FormatV2000NumEntriesField(entriesPerLine)
            << temp.str() << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  " << propCode << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SLBLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    unsigned int id = (*sGroupItr)->getId();
    if (id > 0) {  // Write field only if specific id was assigned
      temp << FormatV2000IntField(idx) << FormatV2000IntField(id);
      if (++count == 8) {
        ret << "M  SLB" << FormatV2000NumEntriesField(8) << temp.str()
            << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SLB" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SDSLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    // Write field only if defined
    if ((*sGroupItr)->hasProp("ESTATE") &&
        (*sGroupItr)->getProp("ESTATE") == "E") {
      temp << FormatV2000IntField(idx);
      if (++count == 15) {
        ret << "M  SDS EXP" << FormatV2000NumEntriesField(15) << temp.str()
            << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SDS EXP" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SPLLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;
  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    auto parent = (*sGroupItr)->getParent();
    if (parent) {  // Write field only if a parent is defined
      unsigned int parentIdx = 1 + parent->getIndexInMol();
      temp << FormatV2000IntField(idx) << FormatV2000IntField(parentIdx);
      if (++count == 8) {
        ret << "M  SPL" << FormatV2000NumEntriesField(8) << temp.str()
            << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SPL" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SNCLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;
  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    unsigned int compno = (*sGroupItr)->getCompNo();
    if (compno > 0) {  // Write field only if compno is set
      temp << FormatV2000IntField(idx) << FormatV2000IntField(compno);
      if (++count == 8) {
        ret << "M  SNC" << FormatV2000NumEntriesField(8) << temp.str()
            << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SNC" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SBTLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    if ((*sGroupItr)->hasProp("BRKTYP")) {
      std::string bracketType = (*sGroupItr)->getProp("BRKTYP");
      if (bracketType == "BRACKET") {
        temp << FormatV2000IntField(idx) << FormatV2000IntField(0);
      } else if (bracketType == "PAREN") {
        temp << FormatV2000IntField(idx) << FormatV2000IntField(1);
      } else {
        std::ostringstream errout;
        errout << "Invalid BRKTYP value '" << bracketType << "' for SGroup "
               << idx;
        throw SGroupException(errout.str());
      }
      if (++count == 8) {
        ret << "M  SBT" << FormatV2000NumEntriesField(8) << temp.str()
            << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SBT" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

template <class T>
std::string BuildV2000IdxVectorDataLines(const unsigned int entriesPerLine,
                                         const unsigned int sGroupId,
                                         const std::string &code,
                                         const T &dataVector) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (const auto &element : dataVector) {
    temp << FormatV2000IntField(1 + element->getIdx());
    if (++count == entriesPerLine) {
      ret << "M  " << code << FormatV2000IntField(sGroupId)
          << FormatV2000NumEntriesField(entriesPerLine) << temp.str()
          << std::endl;
      temp.str("");
      count = 0;
    }
  }
  if (count) {
    ret << "M  " << code << FormatV2000IntField(sGroupId)
        << FormatV2000NumEntriesField(count) << temp.str() << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SMTLine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  if (sgroup->getType() == "MUL" && sgroup->hasProp("MULT")) {
    ret << "M  SMT" << FormatV2000IntField(idx)
        << FormatV2000StringField(sgroup->getProp("MULT"), 69, false, true)
        << std::endl;
  } else if (sgroup->hasProp("LABEL")) {
    ret << "M  SMT" << FormatV2000IntField(idx)
        << FormatV2000StringField(sgroup->getProp("LABEL"), 69, false, true)
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SDILine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  const std::vector<SGroup::Bracket> brackets = sgroup->getBrackets();

  for (const auto &bracket : brackets) {
    ret << "M  SDI" << FormatV2000IntField(idx)
        << FormatV2000NumEntriesField(4);

    for (unsigned int iPoint = 0; iPoint < 2; ++iPoint) {
      ret << FormatV2000DoubleField(bracket.at(iPoint).x);
      ret << FormatV2000DoubleField(bracket.at(iPoint).y);
    }
    ret << std::endl;
  }

  return ret.str();
}

std::string BuildV2000SBVLine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  for (const auto &cstate : sgroup->getCStates()) {
    ret << "M  SBV" << FormatV2000IntField(idx)
        << FormatV2000IntField(cstate.bond->getIdx() + 1);
    if (cstate.vector) {
      ret << FormatV2000DoubleField(cstate.vector->x);
      ret << FormatV2000DoubleField(cstate.vector->y);
    }
    ret << std::endl;
  }

  return ret.str();
}

std::string BuildV2000SDTLine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  if (sgroup->hasProp("FIELDNAME")) {
    ret << "M  SDT" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sgroup->getProp("FIELDNAME"), 30, true, true);

    if (sgroup->hasProp("FIELDTYPE")) {
      ret << FormatV2000StringField(sgroup->getProp("FIELDTYPE"), 2, true,
                                    false);
    } else {
      ret << " T";
    }

    if (sgroup->hasProp("FIELDINFO")) {
      ret << FormatV2000StringField(sgroup->getProp("FIELDINFO"), 20, true,
                                    false);
    }

    if (sgroup->hasProp("QUERYTYPE")) {
      ret << FormatV2000StringField(sgroup->getProp("QUERYTYPE"), 2, true,
                                    false);
    }
    if (sgroup->hasProp("QUERYOP")) {
      ret << FormatV2000StringField(sgroup->getProp("QUERYOP"), 15, true,
                                    false);
    }

    ret << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SDDLine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  if (sgroup->hasProp("FIELDDISP")) {
    ret << "M  SDD" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sgroup->getProp("FIELDDISP"), 69, false,
                                  true);
    ret << std::endl;
  }

  return ret.str();
}

std::string BuildV2000SCDSEDLines(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  for (const auto &data : sgroup->getDataFields()) {
    unsigned int length = data.size();
    if (length > 200) {
      std::ostringstream errout;
      errout << "Data field '" << data << "' in SGroup " << sgroup->getId()
             << " is longer than limit of 200 characters.";
      throw SGroupException(errout.str());
    }
    unsigned int start = 0;
    unsigned int end = 69;
    for (; length > end; start += 69, end += 69) {
      std::string dataChunk = data.substr(start, end);
      ret << "M  SCD" << FormatV2000IntField(idx)
          << FormatV2000StringField(dataChunk, 69, true, true) << std::endl;
    }
    std::string dataChunk = data.substr(start, end);
    ret << "M  SED" << FormatV2000IntField(idx)
        << FormatV2000StringField(dataChunk, 69, false, true) << std::endl;
  }

  return ret.str();
}

std::string BuildV2000PXALine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  if (sgroup->hasProp("PXA")) {
    ret << "M  PXA" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sgroup->getProp("PXA"), 69, false, true);
    ret << std::endl;
  }

  return ret.str();
}

std::string BuildV2000SAPLines(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;
  std::ostringstream temp;

  const std::vector<SGroup::AttachPoint> saps = sgroup->getAttachPoints();

  unsigned int count = 0;
  unsigned int entriesPerLine = 6;
  for (const auto &sap : saps) {
    temp << FormatV2000IntField(1 + sap.aAtom->getIdx());
    if (sap.lvAtom != nullptr) {
      temp << FormatV2000IntField(1 + sap.lvAtom->getIdx());
    } else {
      temp << FormatV2000IntField(0);
    }
    temp << FormatV2000StringField(sap.id, 2, false, true);
    if (++count == entriesPerLine) {
      ret << "M  SAP" << FormatV2000IntField(idx)
          << FormatV2000IntField(entriesPerLine) << temp.str() << std::endl;
      temp.str("");
      count = 0;
    }
  }
  if (count) {
    ret << "M  SAP" << FormatV2000IntField(idx)
        << FormatV2000NumEntriesField(count) << temp.str() << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SCLLine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  if (sgroup->hasProp("CLASS")) {
    ret << "M  SCL" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sgroup->getProp("CLASS"), 69, false, true);
    ret << std::endl;
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
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = 1 + (sGroupItr - mol.beginSGroups());
    ret << BuildV2000IdxVectorDataLines(15, idx, "SAL",
                                        (*sGroupItr)->getAtoms());
    ret << BuildV2000IdxVectorDataLines(15, idx, "SPA",
                                        (*sGroupItr)->getPAtoms());
    ret << BuildV2000IdxVectorDataLines(15, idx, "SBL",
                                        (*sGroupItr)->getBonds());
    // Write CRS line -- CRS still not supported
    ret << BuildV2000SDILine(idx, sGroupItr->get());

    ret << BuildV2000SMTLine(idx, sGroupItr->get());
    ret << BuildV2000SBVLine(idx, sGroupItr->get());

    ret << BuildV2000SDTLine(idx, sGroupItr->get());
    ret << BuildV2000SDDLine(idx, sGroupItr->get());
    // SCD/SED must come after SDT
    ret << BuildV2000SCDSEDLines(idx, sGroupItr->get());

    ret << BuildV2000PXALine(idx, sGroupItr->get());
    ret << BuildV2000SAPLines(idx, sGroupItr->get());
    ret << BuildV2000SCLLine(idx, sGroupItr->get());
  }

  return ret.str();
}

/* ------------------ V3000 Utils  ------------------ */

template <class T>
std::string BuildV3000IdxVectorDataBlock(const std::string &key,
                                         const std::vector<T *> &dataVector) {
  return BuildV3000IdxVectorDataBlock(key, dataVector.begin(),
                                      dataVector.end());
}

template <class Iterator>
std::string BuildV3000IdxVectorDataBlock(const std::string &key,
                                         const Iterator &dataVectorBegin,
                                         const Iterator &dataVectorEnd) {
  using T = typename std::iterator_traits<Iterator>::value_type;

  std::ostringstream ret;

  unsigned int size = dataVectorEnd - dataVectorBegin;

  if (size) {
    auto getOutputIdx = [](T element) -> std::string {
      return std::to_string(1 + element->getIdx());
    };

    std::vector<std::string> tempStr(1 + size);
    tempStr[0] = std::to_string(size);

    std::transform(dataVectorBegin, dataVectorEnd, tempStr.begin() + 1,
                   getOutputIdx);

    ret << ' ' << key << "=(" << boost::algorithm::join(tempStr, " ") << ')';
  }

  return ret.str();
}

std::string BuildV3000BondsBlock(const SGROUP_SPTR &sgroup) {
  std::ostringstream ret;

  auto isXBond = [&sgroup](Bond *bond) {
    return SGroup::BondType::XBOND == sgroup->getBondType(bond);
  };

  auto bonds = sgroup->getBonds();
  auto first_cbond = std::stable_partition(bonds.begin(), bonds.end(), isXBond);

  ret << BuildV3000IdxVectorDataBlock("XBONDS", bonds.begin(), first_cbond);
  ret << BuildV3000IdxVectorDataBlock("CBONDS", first_cbond, bonds.end());

  return ret.str();
}

std::string FormatV3000StringPropertyBlock(const std::string &prop,
                                           const SGROUP_SPTR &sgroup) {
  std::ostringstream ret;

  if (sgroup->hasProp(prop)) {
    ret << ' ' << prop << '=';
    std::string value = sgroup->getProp(prop);
    bool hasSpaces = (value.end() != find(value.begin(), value.end(), ' '));

    if (hasSpaces || value.empty()) {
      ret << "\"";
    }

    ret << value;

    if (hasSpaces || value.empty()) {
      ret << "\"";
    }
  }

  return ret.str();
}

std::string FormatV3000ParentBlock(const SGROUP_SPTR &sgroup) {
  std::ostringstream ret;

  auto parent = sgroup->getParent();

  if (parent) {
    unsigned int parentIdx = 1 + parent->getIndexInMol();
    ret << " PARENT=" << parentIdx;
  }

  return ret.str();
}

std::string FormatV3000CompNoBlock(const SGROUP_SPTR &sgroup) {
  std::ostringstream ret;

  unsigned int compno = sgroup->getCompNo();

  if (compno > 0) {
    ret << " COMPNO=" << compno;
  }

  return ret.str();
}

std::string FormatV3000BracketBlock(
    const std::vector<SGroup::Bracket> brackets) {
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

std::string FormatV3000FieldDataBlock(const SGROUP_SPTR &sgroup) {
  std::ostringstream ret;

  for (const auto &data : sgroup->getDataFields()) {
    ret << " FIELDDATA=\"" << data << "\"";
  }

  return ret.str();
}

std::string FormatV3000CStateBlock(const std::vector<SGroup::CState> &cstates) {
  std::ostringstream ret;

  for (const auto &cstate : cstates) {
    unsigned int xbondIdx = 1 + cstate.bond->getIdx();
    ret << " CSTATE=(";
    if (cstate.vector) {
      ret << "4 " << xbondIdx;
      ret << ' ' << FormatV3000DoubleField(cstate.vector->x);
      ret << ' ' << FormatV3000DoubleField(cstate.vector->y);
      ret << " 0";
    } else {
      ret << "1 " << xbondIdx;
    }
    ret << ")";
  }

  return ret.str();
}

std::string FormatV3000AttachPointBlock(
    const std::vector<SGroup::AttachPoint> &saps) {
  std::ostringstream ret;

  for (const auto &sap : saps) {
    ret << " SAP=(3 " << (1 + sap.aAtom->getIdx());

    if (sap.aAtom == sap.lvAtom) {
      ret << " aidx";
    } else if (sap.lvAtom) {
      ret << ' ' << (1 + sap.lvAtom->getIdx());
    } else {
      ret << " 0";
    }

    ret << ' ' << sap.id << ")";
  }

  return ret.str();
}

const std::string GetV3000MolFileSGroupLines(const unsigned int idx,
                                             const SGROUP_SPTR &sgroup) {
  std::ostringstream os;

  os << ' ' << idx;
  os << ' ' << sgroup->getType();
  os << ' ' << sgroup->getId();

  os << BuildV3000IdxVectorDataBlock("ATOMS", sgroup->getAtoms());
  os << BuildV3000BondsBlock(sgroup);
  os << BuildV3000IdxVectorDataBlock("PATOMS", sgroup->getPAtoms());
  os << FormatV3000StringPropertyBlock("SUBTYPE", sgroup);
  os << FormatV3000StringPropertyBlock("MULT", sgroup);
  os << FormatV3000StringPropertyBlock("CONNECT", sgroup);
  os << FormatV3000ParentBlock(sgroup);
  os << FormatV3000CompNoBlock(sgroup);
  // XBHEAD -> part of V2000 CRS, not supported yet
  // XBCORR -> part of V2000 CRS, not supported yet
  os << FormatV3000StringPropertyBlock("LABEL", sgroup);
  os << FormatV3000BracketBlock(sgroup->getBrackets());
  os << FormatV3000StringPropertyBlock("ESTATE", sgroup);
  os << FormatV3000CStateBlock(sgroup->getCStates());
  os << FormatV3000StringPropertyBlock("FIELDNAME", sgroup);
  os << FormatV3000StringPropertyBlock("FIELDINFO", sgroup);
  os << FormatV3000StringPropertyBlock("FIELDDISP", sgroup);
  os << FormatV3000StringPropertyBlock("QUERYTYPE", sgroup);
  os << FormatV3000StringPropertyBlock("QUERYOP", sgroup);
  os << FormatV3000FieldDataBlock(sgroup);
  os << FormatV3000StringPropertyBlock("CLASS", sgroup);
  os << FormatV3000AttachPointBlock(sgroup->getAttachPoints());
  os << FormatV3000StringPropertyBlock("BRKTYP", sgroup);
  os << FormatV3000StringPropertyBlock("SEQID", sgroup);

  std::string sGroupBlock = os.str();
  unsigned int length = sGroupBlock.size();
  os.str("");

  unsigned int start = 0;
  while (length - start > 73) {
    os << "M  V30 " << sGroupBlock.substr(start, 72) << '-' << std::endl;
    start += 72;
  }
  os << "M  V30 " << sGroupBlock.substr(start, 73) << std::endl;

  return os.str();
}

}  // namespace SGroupWriting
}  // namespace RDKit
