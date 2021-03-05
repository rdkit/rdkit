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

using namespace RDKit;

TEST_CASE("Torsions not found in fused macrocycles", "[macrocycles]") {
  RDLog::InitLogs();
  SECTION("build ff") {
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
}