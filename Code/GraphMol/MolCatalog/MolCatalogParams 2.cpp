// $Id$
//
//  Copyright (C) 2006 Greg Landrum
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

void MolCatalogParams::toStream(std::ostream &ss) const {
  RDUNUSED_PARAM(ss);
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

void MolCatalogParams::initFromStream(std::istream &ss) {
  RDUNUSED_PARAM(ss);
  // at the moment this is a no-op
}
}
