// $Id$
//
//  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <iostream>
#include <fstream>
#include <memory>

#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/StreamOps.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include "ConformerParser.h"

#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
using namespace RDKit;
using namespace RDKit::ConformerParser;

void test1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test ConformerParser." << std::endl;

  std::unique_ptr<ROMol> mol{SmilesToMol("CCC")};
  std::vector<std::vector<double>> coords;
  std::string rdbase = getenv("RDBASE");
  std::string fName = rdbase + "/Code/GraphMol/test_data/water_coords_bad.trx";
  bool ok = false;
  try {
    readAmberTrajectory(fName, coords, mol->getNumAtoms());
  } catch (ValueErrorException &) {
    ok = true;
  }
  TEST_ASSERT(ok);

  fName = rdbase + "/Code/GraphMol/test_data/water_coords_bad2.trx";
  ok = false;
  try {
    readAmberTrajectory(fName, coords, mol->getNumAtoms());
  } catch (ValueErrorException &) {
    // std::cout << e.what() << std::endl;
    ok = true;
  }
  TEST_ASSERT(ok);

  fName = rdbase + "/Code/GraphMol/test_data/water_coords.trx";
  readAmberTrajectory(fName, coords, mol->getNumAtoms());
  TEST_ASSERT(coords.size() == 1);
  TEST_ASSERT(coords[0].size() == 9);
  INT_VECT res = addConformersFromList(*mol, coords);
  TEST_ASSERT((*mol).getNumConformers() == 1);
  TEST_ASSERT((*mol).getConformer().getNumAtoms() == 3);
  TEST_ASSERT((*mol).getConformer(0).getAtomPos(0).x == 0.1941767);

  coords.push_back(coords[0]);
  mol->clearConformers();
  res = addConformersFromList(*mol, coords);
  TEST_ASSERT(mol->getNumConformers() == 2);
  mol->clearConformers();
  res = addConformersFromList(*mol, coords, 1);
  TEST_ASSERT(mol->getNumConformers() == 1);

  coords.resize(0);
  mol->clearConformers();
  fName = rdbase + "/Code/GraphMol/test_data/water_coords2.trx";
  readAmberTrajectory(fName, coords, mol->getNumAtoms());
  TEST_ASSERT(coords.size() == 2);
  TEST_ASSERT(coords[1].size() == 9);
  res = addConformersFromList(*mol, coords);
  TEST_ASSERT((*mol).getNumConformers() == 2);

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main() {
  RDLog::InitLogs();
#if 1
  test1();
#endif
}
