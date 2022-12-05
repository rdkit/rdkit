//
//  Copyright (C) 2022 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "catch.hpp"

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/test.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/Fingerprints/TopologicalTorsionGenerator.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>

using namespace RDKit;

TEST_CASE("includeRedundantEnvironments") {
  auto mol = "CC(=O)O"_smiles;
  REQUIRE(mol);
  SECTION("basics") {
    {
      std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen{
          MorganFingerprint::getMorganGenerator<std::uint32_t>(2)};
      REQUIRE(fpgen);
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp{
          fpgen->getCountFingerprint(*mol)};
      REQUIRE(fp);
      CHECK(fp->getTotalVal() == 8);
    }
    {
      // turn on inclusion of redundant bits
      unsigned int radius = 2;
      bool countSimulation = false;
      bool includeChirality = false;
      bool useBondTypes = true;
      bool onlyNonzeroInvariants = false;
      bool includeRedundantEnvironments = true;
      std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen{
          MorganFingerprint::getMorganGenerator<std::uint32_t>(
              2, countSimulation, includeChirality, useBondTypes,
              onlyNonzeroInvariants, includeRedundantEnvironments)};
      REQUIRE(fpgen);
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp{
          fpgen->getCountFingerprint(*mol)};
      REQUIRE(fp);
      CHECK(fp->getTotalVal() == 12);
    }
  }
}