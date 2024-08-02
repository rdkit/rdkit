//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/test_fixtures.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("github #7556: chiral sulfur in conjugated rings") {
  SECTION("as reported") {
    auto m = "CC1=CC(Cl)=CC2=C1N=[S@](C)N=C2N"_smiles;
    REQUIRE(m);
    CHECK(!m->getBondBetweenAtoms(8, 9)->getIsConjugated());
    CHECK(!m->getBondBetweenAtoms(9, 11)->getIsConjugated());
    REQUIRE(m->getAtomWithIdx(9)->getChiralTag() != Atom::CHI_UNSPECIFIED);
  }
}