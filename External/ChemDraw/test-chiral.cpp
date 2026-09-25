//
//  Copyright (c) 2024 Glysade Inc and other RDkit contributors
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

#include <filesystem>

using namespace RDKit;
using namespace RDKit::v2;

TEST_CASE("Geometry") {
  std::string path =
      std::string(getenv("RDBASE")) + "/External/ChemDraw/test_data/";
  SECTION("R/S Tetrahedral") {
    //_sleep(10 * 1000);

    {
      auto fname = path + "geometry-tetrahedral.cdxml";
      auto mols = MolsFromChemDrawFile(fname);
      REQUIRE(mols.size());  // [C@H]1(C2)[C@@H]2C1
      auto mol = "[C@H]1(C2)[C@@H]2C1"_smiles;
      auto smi = MolToSmiles(*mol);
      REQUIRE(smi == MolToSmiles(*mols[0]));
    }
    {
      auto fname = path + "geometry-tetrahedral-2.cdxml";
      auto mols = MolsFromChemDrawFile(fname);
      REQUIRE(mols.size());
      auto mol = "[C@H]1(C2)[C@@H]2C1"_smiles;
      auto smi = MolToSmiles(*mol);
      REQUIRE(smi == MolToSmiles(*mols[0]));
    }

    {
      auto fname = path + "geometry-tetrahedral-3.cdxml";
      auto mols = MolsFromChemDrawFile(fname);
      REQUIRE(mols.size());
      auto mol = "C1CC[C@H]2CCCC[C@@H]2C1"_smiles;
      auto smi = MolToSmiles(*mol);
      REQUIRE(smi == MolToSmiles(*mols[0]));
    }

    /* this one we still get wrong...
    {
      auto fname = path + "geometry-tetrahedral-4.cdxml";
      auto mols = MolsFromChemDrawFile(fname);
      REQUIRE(mols.size());
      auto mol =
    "CC(S[C@@H]1CC2=C([H])C(CC[C@]2(C)[C@@]3([H])CC([H])([H])[C@]4(C)[C@](OC5=O)(CC5([H])[H])CC[C@@]4([H])[C@]13[H])=O)=O"_smiles;
      auto smi = MolToSmiles(*mol);
      std::cerr << "** " << smi << std::endl;
      REQUIRE(smi == MolToSmiles(*mols[0]));
    }
    */
  }
}
