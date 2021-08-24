//
//  Copyright (C) 2018-2021 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "FragmentCatalogParams.h"
#include "FragmentCatalogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <sstream>

namespace RDKit {
namespace MolStandardize {

#include "fragmentPatterns.in"

FragmentCatalogParams::FragmentCatalogParams(const std::string &fgroupFile) {
  d_funcGroups.clear();
  if (fgroupFile.empty()) {
    d_funcGroups = readFuncGroups(defaults::defaultFragments);
  } else {
    d_funcGroups = readFuncGroups(fgroupFile);
  }
}

FragmentCatalogParams::FragmentCatalogParams(std::istream &fgroupStream) {
  d_funcGroups.clear();
  d_funcGroups = readFuncGroups(fgroupStream);
}

FragmentCatalogParams::FragmentCatalogParams(
    const std::vector<std::pair<std::string, std::string>> &data) {
  d_funcGroups.clear();
  d_funcGroups = readFuncGroups(data);
};

FragmentCatalogParams::FragmentCatalogParams(
    const FragmentCatalogParams &other) {
  d_typeStr = other.d_typeStr;
  d_funcGroups.clear();

  const std::vector<std::shared_ptr<ROMol>> &ofgrps = other.getFuncGroups();
  for (auto &fgi : ofgrps) {
    std::shared_ptr<ROMol> nmol(new ROMol(*fgi));
    d_funcGroups.push_back(nmol);
  }
}

FragmentCatalogParams::~FragmentCatalogParams() {}

const std::vector<std::shared_ptr<ROMol>>
    &FragmentCatalogParams::getFuncGroups() const {
  return d_funcGroups;
}

const ROMol *FragmentCatalogParams::getFuncGroup(unsigned int fid) const {
  URANGE_CHECK(fid, d_funcGroups.size());
  return d_funcGroups[fid].get();
}

void FragmentCatalogParams::toStream(std::ostream &ss) const {
  ss << d_funcGroups.size() << "\n";
  //	for (const auto &d_funcGroup : d_funcGroups) {
  //		std::string text;
  //		d_funcGroup->getProp(common_properties::_Name, text);
  //		ss << text;
  //		ss << "\t";
  //		d_funcGroup->getProp(common_properties::_fragSMARTS, text);
  //		ss << text;
  //		ss << "\n";
  //	}
}

std::string FragmentCatalogParams::Serialize() const {
  std::stringstream ss;
  toStream(ss);
  return ss.str();
}

void FragmentCatalogParams::initFromStream(std::istream &) {
  UNDER_CONSTRUCTION("not implemented");
}

void FragmentCatalogParams::initFromString(const std::string &) {
  UNDER_CONSTRUCTION("not implemented");
}

}  // namespace MolStandardize
}  // namespace RDKit
