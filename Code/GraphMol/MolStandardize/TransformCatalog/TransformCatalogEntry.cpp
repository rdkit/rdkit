//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "TransformCatalogEntry.h"

#include <RDGeneral/types.h>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/StreamOps.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <iostream>
#include <fstream>

namespace RDKit {
namespace MolStandardize {

void TransformCatalogEntry::toStream(std::ostream &ss) const {
  ReactionPickler::pickleReaction(*dp_transform, ss);

  std::int32_t tmpInt;
  tmpInt = getBitId();
  streamWrite(ss, tmpInt);

  tmpInt = d_descrip.size();
  streamWrite(ss, tmpInt);
  ss.write(d_descrip.c_str(), tmpInt * sizeof(char));
}

std::string TransformCatalogEntry::Serialize() const {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);
  toStream(ss);
  return ss.str();
}

void TransformCatalogEntry::initFromStream(std::istream &ss) {
  // the reaction:
  dp_transform = new ChemicalReaction();
  ReactionPickler::reactionFromPickle(ss, *dp_transform);

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

void TransformCatalogEntry::initFromString(const std::string &text) {
  std::stringstream ss(std::ios_base::binary | std::ios_base::out |
                       std::ios_base::in);
  // initialize the stream:
  ss.write(text.c_str(), text.length());
  // now start reading out values:
  initFromStream(ss);
}

}  // namespace MolStandardize
}  // namespace RDKit
