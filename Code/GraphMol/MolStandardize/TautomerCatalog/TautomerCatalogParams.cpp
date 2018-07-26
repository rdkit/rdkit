#include "TautomerCatalogParams.h"
#include "TautomerCatalogUtils.h"
#include <GraphMol/RDKitBase.h>
#include <sstream>

namespace RDKit {
namespace MolStandardize {

TautomerCatalogParams::TautomerCatalogParams(const std::string &tautomerFile) {
  d_transforms.clear();
  d_transforms = readTautomers(tautomerFile);
}

TautomerCatalogParams::TautomerCatalogParams(
    const TautomerCatalogParams &other) {
  d_typeStr = other.d_typeStr;
  d_transforms.clear();

  const std::vector<TautomerTransform> &transforms = other.getTransforms();
  for (const auto &transform : transforms) {
    d_transforms.push_back(transform);
  }
}

TautomerCatalogParams::~TautomerCatalogParams() {}

const std::vector<TautomerTransform> &TautomerCatalogParams::getTransforms()
    const {
  return d_transforms;
}

const TautomerTransform TautomerCatalogParams::getTransform(
    unsigned int fid) const {
  URANGE_CHECK(fid, d_transforms.size());
  return d_transforms[fid];  //.get();
}

void TautomerCatalogParams::toStream(std::ostream &ss) const {
  ss << d_transforms.size() << "\n";
}

std::string TautomerCatalogParams::Serialize() const {
  std::stringstream ss;
  toStream(ss);
  return ss.str();
}

void TautomerCatalogParams::initFromStream(std::istream &ss) {}

void TautomerCatalogParams::initFromString(const std::string &text) {}

}  // namespace MolStandardize
}  // namespace RDKit
