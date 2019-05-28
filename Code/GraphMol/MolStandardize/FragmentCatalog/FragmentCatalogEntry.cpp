//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "FragmentCatalogEntry.h"

#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/StreamOps.h>
#include <GraphMol/MolPickler.h>
#include <iostream>
#include <fstream>

namespace RDKit {
namespace MolStandardize {

// FragmentCatalogEntry::FragmentCatalogEntry(const ROMol *omol, const PATH_TYPE
// &path) { 	PRECONDITION(omol, "bad mol");
//}
//
// FragmentCatalogEntry::FragmentCatalogEntry(const std::string &pickle) {
//	dp_props = new Dict();
//	this->initFromString(pickle);
//}

void FragmentCatalogEntry::toStream(std::ostream &ss) const {
  MolPickler::pickleMol(*dp_mol, ss);

  std::int32_t tmpInt;
  tmpInt = getBitId();
  streamWrite(ss, tmpInt);

  tmpInt = d_descrip.size();
  streamWrite(ss, tmpInt);
  ss.write(d_descrip.c_str(), tmpInt * sizeof(char));
}

std::string FragmentCatalogEntry::Serialize() const {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);
  toStream(ss);
  return ss.str();
}

void FragmentCatalogEntry::initFromStream(std::istream &ss) {
  // the molecule:
  dp_mol = new ROMol();
  MolPickler::molFromPickle(ss, *dp_mol);

  std::int32_t tmpInt;
  // the bitId:
  streamRead(ss, tmpInt);
  setBitId(tmpInt);

  // the description:
  streamRead(ss, tmpInt);
  char *tmpText = new char[tmpInt + 1];
  ss.read(tmpText, tmpInt * sizeof(char));
  tmpText[tmpInt] = 0;
  d_descrip = tmpText;
  delete[] tmpText;
}

void FragmentCatalogEntry::initFromString(const std::string &text) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);
  // initialize the stream:
  ss.write(text.c_str(), text.length());
  // now start reading out values:
  initFromStream(ss);
}

}  // namespace MolStandardize
}  // namespace RDKit
