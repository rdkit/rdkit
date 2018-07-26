#include "FragmentCatalogParams.h"
#include "FragmentCatalogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <sstream>

namespace RDKit {
namespace MolStandardize {

FragmentCatalogParams::FragmentCatalogParams(const std::string &fgroupFile) {
  d_funcGroups.clear();
  d_funcGroups = readFuncGroups(fgroupFile);
}

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
  // return d_funcGroups[fid];
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

void FragmentCatalogParams::initFromStream(std::istream &ss) {}

void FragmentCatalogParams::initFromString(const std::string &text) {}

}  // namespace MolStandardize
}  // namespace RDKit
