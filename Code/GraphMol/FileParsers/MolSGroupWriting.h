//
//  Copyright (C) 2002-2021 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#pragma once

#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>
#include <GraphMol/SubstanceGroup.h>

namespace RDKit {
namespace SGroupWriting {
typedef std::unordered_map<int, SubstanceGroup> IDX_TO_SGROUP_MAP;

/* ------------------ Inlined Formatters ------------------ */

inline std::string FormatV2000IntField(int value) {
  auto fmt = boost::format(" %3d") % value;
  return fmt.str();
}

inline std::string FormatV2000NumEntriesField(int value) {
  auto fmt = boost::format(" %2d") % value;
  return fmt.str();
}

inline std::string FormatV2000DoubleField(double value) {
  auto fmt = boost::format("%10.4f") % value;
  return fmt.str();
}

inline std::string FormatV2000StringField(const std::string &value,
                                          unsigned int fieldSize, bool pad,
                                          bool addSeparator) {
  std::ostringstream os;
  if (addSeparator) {
    os << ' ';
  }
  if (value.size() >= fieldSize) {
    os << value.substr(0, fieldSize);
  } else if (pad) {
    os << std::setw(fieldSize) << std::left << value;
  } else {
    os << value;
  }
  return os.str();
}

inline std::string FormatV3000DoubleField(double value) {
  return boost::trim_copy(FormatV2000DoubleField(value));
}

/* ------------------ V2000 Utils  ------------------ */

std::string BuildV2000STYLines(const ROMol &mol);

std::string BuildV2000StringPropLines(const unsigned int entriesPerLine,
                                      const ROMol &mol,
                                      const std::string &propName,
                                      const std::string &propCode,
                                      const unsigned int fieldWitdh);

std::string BuildV2000SLBLines(const ROMol &mol);

std::string BuildV2000SDSLines(const ROMol &mol);

std::string BuildV2000SPLLines(const ROMol &mol);

std::string BuildV2000SNCLines(const ROMol &mol);

std::string BuildV2000SBTLines(const ROMol &mol);

template <class T>
std::string BuildV2000IdxVectorDataLines(const unsigned int entriesPerLine,
                                         const unsigned int sGroupId,
                                         const std::string &code,
                                         const T &dataVector);

std::string BuildV2000SMTLine(const int idx, const SubstanceGroup *sgroup);

std::string BuildV2000SDILine(const int idx, const SubstanceGroup *sgroup);

std::string BuildV2000SBVLine(const int idx, const SubstanceGroup *sgroup);

std::string BuildV2000SDTLine(const int idx, const SubstanceGroup *sgroup);

std::string BuildV2000SDDLine(const int idx, const SubstanceGroup *sgroup);

std::string BuildV2000SCDSEDLines(const int idx, const SubstanceGroup *sgroup);

std::string BuildV2000SAPLines(const int idx, const SubstanceGroup *sgroup);

std::string BuildV2000SCLLine(const int idx, const SubstanceGroup *sgroup);
const std::string GetMolFileSGroupInfo(const RWMol &mol);

/* ------------------ V3000 Utils  ------------------ */

template <class T>
std::string BuildV3000IdxVectorDataBlock(const std::string &key,
                                         const std::vector<T *> &dataVector);

template <class Iterator>
std::string BuildV3000IdxVectorDataBlock(const std::string &key,
                                         const Iterator &dataVectorBegin,
                                         const Iterator &dataVectorEnd);

/* Classify bonds between XBONDS and CBOfindP work on a copy of
 * bonds vector to prevent reordering of original vector */
std::string BuildV3000BondsBlock(const SubstanceGroup &sgroup);

std::string FormatV3000StringPropertyBlock(const std::string &prop,
                                           const SubstanceGroup &sgroup);

std::string FormatV3000ParentBlock(const SubstanceGroup &sgroup);

std::string FormatV3000CompNoBlock(const SubstanceGroup &sgroup);

std::string FormatV3000BracketBlock(
    const std::vector<SubstanceGroup::Bracket> brackets);

std::string FormatV3000CStateBlock(
    const std::vector<SubstanceGroup::CState> &cstates);

const std::string GetV3000MolFileSGroupLines(const unsigned int idx,
                                             const SubstanceGroup &sgroup);
}  // namespace SGroupWriting
}  // namespace RDKit
