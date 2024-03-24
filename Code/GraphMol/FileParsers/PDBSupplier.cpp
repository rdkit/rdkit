//
//  Copyright (C) 2013-2024 Greg Landrum, NextMove Software,
//    and other RDKit contributors
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

namespace v2 {
namespace FileParsers {
PDBMolSupplier::PDBMolSupplier(std::istream *inStream, bool takeOwnership,
                               const PDBParserParams &params) {
  dp_inStream = inStream;
  df_owner = takeOwnership;
  d_params = params;
}

PDBMolSupplier::PDBMolSupplier(const std::string &fileName,
                               const PDBParserParams &params) {
  dp_inStream = openAndCheckStream(fileName);
  df_owner = true;
  d_params = params;
}

void PDBMolSupplier::init() {}
void PDBMolSupplier::reset() {}

std::unique_ptr<RWMol> PDBMolSupplier::next() {
  return MolFromPDBDataStream(*dp_inStream, d_params);
}

bool PDBMolSupplier::atEnd() {
  if (dp_inStream->eof()) {
    return true;
  }
  int ch = dp_inStream->peek();
  return ch == -1;
}
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
