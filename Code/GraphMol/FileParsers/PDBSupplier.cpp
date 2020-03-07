//
//  Copyright (C) 2013 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <fstream>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>

namespace RDKit {

PDBMolSupplier::PDBMolSupplier(std::istream *inStream, bool takeOwnership,
                               bool sanitize, bool removeHs,
                               unsigned int flavor, bool proximityBonding) {
  dp_inStream = inStream;
  df_owner = takeOwnership;
  df_sanitize = sanitize;
  df_removeHs = removeHs;
  d_flavor = flavor;
  df_proximityBonding = proximityBonding;
}

PDBMolSupplier::PDBMolSupplier(const std::string &fileName, bool sanitize,
                               bool removeHs, unsigned int flavor,
                               bool proximityBonding) {
  dp_inStream = openAndCheckStream(fileName);
  df_owner = true;
  df_sanitize = sanitize;
  df_removeHs = removeHs;
  d_flavor = flavor;
  df_proximityBonding = proximityBonding;
}

void PDBMolSupplier::init() {}
void PDBMolSupplier::reset() {}

ROMol *PDBMolSupplier::next() {
  return (ROMol *)PDBDataStreamToMol(dp_inStream, df_sanitize, df_removeHs,
                                     d_flavor, df_proximityBonding);
}

bool PDBMolSupplier::atEnd() {
  if (dp_inStream->eof()) {
    return true;
  }
  int ch = dp_inStream->peek();
  return ch == -1;
}
}  // namespace RDKit
