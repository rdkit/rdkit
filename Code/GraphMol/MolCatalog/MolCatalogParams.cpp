// $Id$
//
//  Copyright (C) 2006 Greg Landrum
//
//
#include "MolCatalogParams.h"
#include <sstream>

namespace RDKit{

  MolCatalogParams::MolCatalogParams(const std::string &pickle) {
    d_typeStr = "MolCatalog Parameters";
    this->initFromString(pickle);
  }

  MolCatalogParams::~MolCatalogParams() {
  }

  void MolCatalogParams::toStream(std::ostream &ss) const {
    // at the moment this is a no-op
  }
  std::string MolCatalogParams::Serialize() const {
    std::stringstream ss;
    toStream(ss);
    return ss.str();
  }

  void MolCatalogParams::initFromString(const std::string &text){
    std::stringstream ss(text);
    initFromStream(ss);
  }
  
  void MolCatalogParams::initFromStream(std::istream &ss){
    // at the moment this is a no-op
  }
}
