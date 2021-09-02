//
//  Copyright (C) 2006-2021 Greg Landrum and other RDKit contributors
//
//
#include "MolCatalogParams.h"
#include <sstream>
#include <RDGeneral/Invariant.h>

namespace RDKit {

MolCatalogParams::MolCatalogParams(const std::string &pickle) {
  d_typeStr = "MolCatalog Parameters";
  this->initFromString(pickle);
}

MolCatalogParams::~MolCatalogParams() {}

void MolCatalogParams::toStream(std::ostream &) const {
  // at the moment this is a no-op
}
std::string MolCatalogParams::Serialize() const {
  std::stringstream ss;
  toStream(ss);
  return ss.str();
}

void MolCatalogParams::initFromString(const std::string &text) {
  std::stringstream ss(text);
  initFromStream(ss);
}

void MolCatalogParams::initFromStream(std::istream &) {
  // at the moment this is a no-op
}
}  // namespace RDKit
