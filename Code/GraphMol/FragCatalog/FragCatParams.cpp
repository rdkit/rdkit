// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "FragCatParams.h"
#include "FragCatalogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <sstream>
#include <cstdint>

namespace RDKit {
using std::int32_t;
using std::uint32_t;

FragCatParams::FragCatParams(unsigned int lLen, unsigned int uLen,
                             const std::string &fgroupFile, double tol) {
  d_funcGroups.clear();
  d_typeStr = "Fragment Catalog Parameters";
  CHECK_INVARIANT(lLen <= uLen,
                  "The upper length for fragments must be >= lower length");
  d_lowerFragLen = lLen;
  d_upperFragLen = uLen;
  d_tolerance = tol;
  d_funcGroups = readFuncGroups(fgroupFile);
}

FragCatParams::FragCatParams(const FragCatParams &other) {
  d_funcGroups.clear();
  // copy consttructor
  d_typeStr = other.getTypeStr();
  d_lowerFragLen = other.getLowerFragLength();
  d_upperFragLen = other.getUpperFragLength();
  d_tolerance = other.getTolerance();

  // std::cout << "In param copier\n";
  const MOL_SPTR_VECT &ofgrps = other.getFuncGroups();
  // const MOL_PTR_VECT &ofgrps = other.getFuncGroups();
  MOL_SPTR_VECT::const_iterator fgi;
  // MOL_PTR_VECT_CI fgi;
  for (fgi = ofgrps.begin(); fgi != ofgrps.end(); fgi++) {
    auto *nmol = new ROMol(*(fgi->get()));
    // ROMol *nmol = new ROMol(*(*fgi));
    d_funcGroups.push_back(ROMOL_SPTR(nmol));
    // d_funcGroups.push_back(nmol);
  }
}

FragCatParams::FragCatParams(const std::string &pickle) {
  d_typeStr = "Fragment Catalog Parameters";
  this->initFromString(pickle);
}

FragCatParams::~FragCatParams() {}

const MOL_SPTR_VECT &FragCatParams::getFuncGroups() const {
  return d_funcGroups;
}

const ROMol *FragCatParams::getFuncGroup(unsigned int fid) const {
  URANGE_CHECK(fid, d_funcGroups.size());
  // return d_funcGroups[fid];
  return d_funcGroups[fid].get();
}

void FragCatParams::toStream(std::ostream &ss) const {
  ss << d_lowerFragLen << " " << d_upperFragLen << " " << d_tolerance << "\n";
  ss << d_funcGroups.size() << "\n";
  for (const auto &d_funcGroup : d_funcGroups) {
    std::string text;
    d_funcGroup->getProp(common_properties::_Name, text);
    ss << text;
    ss << "\t";
    d_funcGroup->getProp(common_properties::_fragSMARTS, text);
    ss << text;
    ss << "\n";
  }
}
std::string FragCatParams::Serialize() const {
  std::stringstream ss;
  toStream(ss);
  return ss.str();
}

void FragCatParams::initFromString(const std::string &text) {
  std::stringstream ss(text);
  initFromStream(ss);
}

void FragCatParams::initFromStream(std::istream &ss) {
  ss >> d_lowerFragLen;
  ss >> d_upperFragLen;
  ss >> d_tolerance;
  int nGroups;
  ss >> nGroups;
  // std::cout << "Reading " << nGroups << " Groups" << std::endl;

  d_funcGroups = readFuncGroups(ss, nGroups);
}
}  // namespace RDKit
