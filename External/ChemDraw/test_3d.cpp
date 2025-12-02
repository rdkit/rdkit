//
//  Copyright (c) 2025 Glysade Inc and other RDkit contributors
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

#include "chemdraw.h"
#include <catch2/catch_all.hpp>
#include "RDGeneral/test.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <RDGeneral/FileParseException.h>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SmilesParse/CanonicalizeStereoGroups.h>

using namespace RDKit;
using namespace RDKit::v2;

TEST_CASE("Round TRIP") {
  std::string path =
      std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/";
  std::string code_path = std::string(getenv("RDBASE"));

  // Eventually this catch test is to see if round tripping mol 3d -> chemdraw
  // returns
  //  reasonable coords, however chemdraw seems to forget about the original
  //  scale
  //   and converts to pixel drawing coords, so this test is kind of meaningless
  SECTION("3D structs") {
    auto fname =
        code_path + "/Code/GraphMol/FileParsers/test_data/Issue3514824.mol";
    auto mol = MolFileToMol(fname);
    REQUIRE(mol);
    auto &conf = mol->getConformer(0);
    for (auto bond : mol->bonds()) {
      auto p1 = conf.getAtomPos(bond->getBeginAtomIdx());
      auto p2 = conf.getAtomPos(bond->getEndAtomIdx());
      auto length = (p1 - p2).length();
      std::cerr << bond->getIdx() << " : " << length << std::endl;
      ;
    }
    std::cerr << "----------" << std::endl;
    {
      auto fname2 =
          code_path + "/Code/GraphMol/FileParsers/test_data/Issue3514824.cdxml";
      auto mols = MolsFromChemDrawFile(fname2);
      auto &conf2 = mols[0]->getConformer(0);
      for (auto bond : mols[0]->bonds()) {
        auto p1 = conf2.getAtomPos(bond->getBeginAtomIdx());
        auto p2 = conf2.getAtomPos(bond->getEndAtomIdx());
        auto length = (p1 - p2).length();
        std::cerr << bond->getIdx() << " : " << length << std::endl;
      }
    }
    std::cerr << "----------" << std::endl;
    {
      auto cdx = MolToChemDrawBlock(*mol);
      auto mols = MolsFromChemDrawBlock(cdx);
      auto &conf2 = mols[0]->getConformer(0);
      for (auto bond : mols[0]->bonds()) {
        auto p1 = conf2.getAtomPos(bond->getBeginAtomIdx());
        auto p2 = conf2.getAtomPos(bond->getEndAtomIdx());
        auto length = (p1 - p2).length();
        std::cerr << bond->getIdx() << " : " << length << std::endl;
      }
    }
    delete mol;
  }
}
