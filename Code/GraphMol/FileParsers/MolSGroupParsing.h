//
//  Copyright (C) 2002-2018 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>

#pragma once
#include <GraphMol/SubstanceGroup.h>
#include <RDGeneral/FileParseException.h>
#include <sstream>

namespace RDKit {

namespace SGroupParsing {
typedef std::map<int, SubstanceGroup> IDX_TO_SGROUP_MAP;
typedef std::map<int, STR_VECT> IDX_TO_STR_VECT_MAP;

/* ------------------ V2000 Utils  ------------------ */

unsigned int ParseSGroupIntField(const std::string &text, unsigned int line,
                                 unsigned int &pos,
                                 bool isFieldCounter = false);

unsigned int ParseSGroupIntField(bool &ok, bool strictParsing,
                                 const std::string &text, unsigned int line,
                                 unsigned int &pos,
                                 bool isFieldCounter = false);

double ParseSGroupDoubleField(const std::string &text, unsigned int line,
                              unsigned int &pos);

double ParseSGroupDoubleField(bool &ok, bool strictParsing,
                              const std::string &text, unsigned int line,
                              unsigned int &pos);

SubstanceGroup *FindSgIdx(IDX_TO_SGROUP_MAP &sGroupMap, int sgIdx,
                          unsigned int line);

template <class Exc = FileParseException>
void SGroupWarnOrThrow(bool strictParsing, const std::string &msg) {
  if (strictParsing) {
    throw Exc(msg);
  } else {
    BOOST_LOG(rdWarningLog) << msg << std::endl;
  }
}

void ParseSGroupV2000STYLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000VectorDataLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                                    const std::string &text, unsigned int line,
                                    bool strictParsing = true);

void ParseSGroupV2000SDILine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SSTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SMTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SLBLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SCNLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SDSLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SBVLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SDTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SDDLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SCDSEDLine(IDX_TO_SGROUP_MAP &sGroupMap,
                                IDX_TO_STR_VECT_MAP &dataFieldsMap, RWMol *mol,
                                const std::string &text, unsigned int line,
                                bool strictParsing, unsigned int &counter,
                                unsigned int &lastDataSGroup,
                                std::ostringstream &currentDataField);

void ParseSGroupV2000SPLLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SNCLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

//! if the SAP line is malformed and has no lvIdx and no id,
//! lvIdx is set to mol->getNumAtoms() and id is set to "  "
//! the user is responsible for replacing lvIdx with the correct
//! index: if d_bonds.size() == 1, and one of the bond atom indices
//! is aIdx, the other can be safely assigned to lvIdx
void ParseSGroupV2000SAPLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SCLLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

void ParseSGroupV2000SBTLine(IDX_TO_SGROUP_MAP &sGroupMap, RWMol *mol,
                             const std::string &text, unsigned int line,
                             bool strictParsing = true);

/* ------------------ V3000 Utils  ------------------ */

template <class T>
RDKIT_FILEPARSERS_EXPORT std::vector<T> ParseV3000Array(
    std::stringstream &stream, int maxV = -1, bool strictParsing = false);
#if defined(_MSC_VER) && defined(RDKIT_DYN_LINK)
template RDKIT_FILEPARSERS_EXPORT std::vector<int> ParseV3000Array(
    std::stringstream &, int, bool);
template RDKIT_FILEPARSERS_EXPORT std::vector<unsigned int> ParseV3000Array(
    std::stringstream &, int, bool);
#endif
template <class T>
std::vector<T> ParseV3000Array(const std::string &s, int maxV = -1,
                               bool strictParsing = false) {
  std::stringstream stream(s);
  return ParseV3000Array<T>(stream, maxV, strictParsing);
}

void ParseV3000CStateLabel(RWMol *mol, SubstanceGroup &sgroup,
                           std::stringstream &stream, unsigned int line,
                           bool strictParsing = true);

void ParseV3000SAPLabel(RWMol *mol, SubstanceGroup &sgroup,
                        std::stringstream &stream, bool strictParsing = true);

std::string ParseV3000StringPropLabel(std::stringstream &stream);

// returns the last line read in the SGroups block
std::string ParseV3000SGroupsBlock(std::istream *inStream, unsigned int line,
                                   unsigned int nSgroups, RWMol *mol,
                                   bool strictParsing);

}  // namespace SGroupParsing
}  // namespace RDKit
