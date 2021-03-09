//
//  Copyright (C) 2021 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include "catch.hpp"

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include "Embedder.h"
#include <tuple>

using namespace RDKit;

TEST_CASE("Torsions not found in fused macrocycles", "[macrocycles]") {
  RDLog::InitLogs();
  SECTION("reported") {
    // this is 6VY8 from the PDB
    auto mol1 =
        "CC[C@H](C)[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)[C@@H]2CCCN2C(=O)[C@H]([C@@H](C)CC)NC(=O)[C@H](CO)NC(=O)[C@H](CCCC[NH3+])NC(=O)[C@H]([C@@H](C)O)NC(O)[C@@H]2CN3NNC[C@H]3C[C@H](NC1=O)C(O)N[C@@H](Cc1ccccc1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC(=O)O)C(=O)NCC(=O)N[C@@H](CCCNC(N)=[NH2+])C(=O)N2"_smiles;
    REQUIRE(mol1);
    MolOps::addHs(*mol1);
    ForceFields::CrystalFF::CrystalFFDetails details;
    bool useExpTorsions = true;
    bool useSmallRingTorsions = false;
    bool useMacrocycleTorsions = true;
    bool useBasicKnowledge = true;
    unsigned int version = 2;
    bool verbose = true;
    std::stringstream sstrm;
    rdInfoLog->SetTee(sstrm);
    ForceFields::CrystalFF::getExperimentalTorsions(
        *mol1, details, useExpTorsions, useSmallRingTorsions,
        useMacrocycleTorsions, useBasicKnowledge, version, verbose);
    rdInfoLog->ClearTee();
    auto txt = sstrm.str();
    CHECK(txt.find("{9-}") != std::string::npos);
  }
  SECTION("edges") {
    std::vector<std::tuple<std::string, bool, unsigned int>> tests{
        {"O=C1CNC(=O)C2CCC(N1)NC(=O)CNC2=O", true, 15},  // 9-9
        {"O=C1NC2CCC(C(=O)N1)C(=O)NCC(=O)N2", true, 4},  // 9-8
        {"O=C1NC2CCC(C(=O)N1)C(=O)NC(=O)N2", false, 0},  // 8-8
        {"O=C1CC(=O)NC2NC(=O)CC(=O)NC(N1)NC(=O)CC(=O)N2", true,
         18}};  // 12-12-12
    for (const auto &tpl : tests) {
      std::unique_ptr<RWMol> m{SmilesToMol(std::get<0>(tpl))};
      REQUIRE(m);
      MolOps::addHs(*m);
      ForceFields::CrystalFF::CrystalFFDetails details;
      bool useExpTorsions = true;
      bool useSmallRingTorsions = false;
      bool useMacrocycleTorsions = true;
      bool useBasicKnowledge = true;
      unsigned int version = 2;
      bool verbose = true;
      std::stringstream sstrm;
      rdInfoLog->SetTee(sstrm);
      std::cerr << "-----------" << std::endl;
      ForceFields::CrystalFF::getExperimentalTorsions(
          *m, details, useExpTorsions, useSmallRingTorsions,
          useMacrocycleTorsions, useBasicKnowledge, version, verbose);
      rdInfoLog->ClearTee();
      auto txt = sstrm.str();
      if (std::get<1>(tpl)) {
        CHECK(txt.find("{9-}") != std::string::npos);
      } else {
        CHECK(txt.find("{9-}") == std::string::npos);
      }
      CHECK(details.expTorsionAngles.size() == std::get<2>(tpl));
    }
  }
}
