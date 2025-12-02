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
#include "chemdrawreaction.h"
#include <catch2/catch_all.hpp>
#include "RDGeneral/test.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <filesystem>
using namespace RDKit;
using namespace RDKit::v2;

TEST_CASE("CDXML Parser") {
  std::string cdxmlbase =
      std::string(getenv("RDBASE")) + "/Code/GraphMol/test_data/CDXML/";
  SECTION("CDXML REACTION") {
    auto fname = cdxmlbase + "rxn2.cdxml";
    std::vector<std::string> expected = {
        "Cl[c:1]1[cH:4][cH:3][cH:2][cH:6][cH:5]1",
        "OC(O)B[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1",
        "[cH:1]1[cH:4][cH:3][cH:2][c:6](-[c:7]2[cH:8][cH:9][cH:10][cH:11][cH:12]2)[cH:5]1"};

    auto rxns = ChemDrawFileToChemicalReactions(fname);
    CHECK(rxns.size() == 1);
    unsigned int i = 0;
    int count = 0;
    for (auto &mol : rxns[0]->getReactants()) {
      CHECK(mol->getProp<unsigned int>("CDX_SCHEME_ID") == 397);
      CHECK(mol->getProp<unsigned int>("CDX_STEP_ID") == 398);
      CHECK(mol->getProp<unsigned int>("CDX_REAGENT_ID") == i++);
      CHECK(MolToSmiles(*mol) == expected[count++]);
    }
    i = 0;
    for (auto &mol : rxns[0]->getProducts()) {
      CHECK(mol->getProp<unsigned int>("CDX_SCHEME_ID") == 397);
      CHECK(mol->getProp<unsigned int>("CDX_STEP_ID") == 398);
      CHECK(mol->getProp<unsigned int>("CDX_PRODUCT_ID") == i++);
      CHECK(MolToSmiles(*mol) == expected[count++]);
    }

    auto smarts = ChemicalReactionToRxnSmarts(*rxns[0]);
    CHECK(
        smarts ==
        "[#6&D2:2]1:[#6&D2:3]:[#6&D2:4]:[#6&D3:1](:[#6&D2:5]:[#6&D2:6]:1)-[#17&D1].[#6&D3](-[#5&D2]-[#6&D3:7]1:[#6&D2:8]:[#6&D2:9]:[#6&D2:10]:[#6&D2:11]:[#6&D2:12]:1)(-[#8&D1])-[#8&D1]>>[#6&D2:1]1:[#6&D2:5]:[#6&D3:6](:[#6&D2:2]:[#6&D2:3]:[#6&D2:4]:1)-[#6&D3:7]1:[#6&D2:8]:[#6&D2:9]:[#6&D2:10]:[#6&D2:11]:[#6&D2:12]:1");
  }

  SECTION("Github #7528 CDXML Grouped Agents in Reactions") {
    // The failing case had fragments grouped with labels, ensure the grouped
    // cersion and the ungrouped versions have the same results
    auto fname = cdxmlbase + "github7467-grouped-fragments.cdxml";
    auto rxns = ChemDrawFileToChemicalReactions(fname);
    CHECK(rxns.size() == 1);
    fname = cdxmlbase + "github7467-ungrouped-fragments.cdxml";
    auto rxns2 = ChemDrawFileToChemicalReactions(fname);

    CHECK(ChemicalReactionToRxnSmarts(*rxns[0]) ==
          ChemicalReactionToRxnSmarts(*rxns2[0]));

    // Check to see if our understanding of grouped reagents in reactions is
    // correct
    fname = cdxmlbase + "reaction-with-grouped-templates.cdxml";
    auto rxns3 = ChemDrawFileToChemicalReactions(fname);
    CHECK(rxns3.size() == 1);
    std::string rxnb = R"RXN($RXN

      Mrv2004  062120241319

  2  0
$MOL

  Mrv2004 06212413192D

  5  5  0  0  0  0            999 V2000
    2.6221   -4.6475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6221   -5.4725    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4070   -5.7274    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8918   -5.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4070   -4.3926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  1  1  0  0  0  0
M  END
$MOL

  Mrv2004 06212413192D

 11 11  0  0  0  0            999 V2000
    6.9305   -4.5100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9305   -5.3350    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6450   -5.7475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3594   -5.3350    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3594   -4.5100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6450   -4.0975    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6171   -4.4825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6171   -5.3075    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.4020   -5.5624    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.8868   -4.8950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.4020   -4.2276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  6  1  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  7  8  2  0  0  0  0
 11  7  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  2  0  0  0  0
 10 11  1  0  0  0  0
M  END
)RXN";
    std::unique_ptr<ChemicalReaction> rxn_mb{RxnBlockToChemicalReaction(rxnb)};
    // CDXMLToReaction is sanitized by default, this might be a mistake...
    unsigned int failed;
    RxnOps::sanitizeRxn(
        *rxn_mb, failed,
        RxnOps::SANITIZE_ADJUST_REACTANTS | RxnOps::SANITIZE_ADJUST_PRODUCTS,
        RxnOps::MatchOnlyAtRgroupsAdjustParams());

    CHECK(rxns3[0]->getNumReactantTemplates() ==
          rxn_mb->getNumReactantTemplates());
    CHECK(ChemicalReactionToRxnSmarts(*rxns3[0]) ==
          ChemicalReactionToRxnSmarts(*rxn_mb));
  }
}
