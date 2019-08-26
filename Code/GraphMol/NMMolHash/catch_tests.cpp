//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include "molhash.h"

#include <iostream>
#include <fstream>

using namespace RDKit;

TEST_CASE("Basic MolHash","[molhash]") {
  SECTION("basics") {
    auto om = "C1CCCC(O)C1c1ccnc(OC)c1"_smiles;
    REQUIRE(om);
    {
      std::unique_ptr<RWMol> m(new RWMol(*om));
      auto hsh = MolHash(m.get(),HashFunction::AnonymousGraph);  
      CHECK(hsh == "***1****(*2*****2*)*1");
    }
    {
        std::unique_ptr<RWMol> m(new RWMol(*om));
        auto hsh = MolHash(m.get(),HashFunction::ElementGraph);
      CHECK(hsh == "COC1CC(C2CCCCC2O)CCN1");
    }
    {
        std::unique_ptr<RWMol> m(new RWMol(*om));
        auto hsh = MolHash(m.get(),HashFunction::CanonicalSmiles);
      CHECK(hsh == "COc1cc(C2CCCCC2O)ccn1");
    }
    {
        std::unique_ptr<RWMol> m(new RWMol(*om));
        auto hsh = MolHash(m.get(),HashFunction::MurckoScaffold);  
      CHECK(hsh == "c1cc(C2CCCCC2)ccn1");
    }
    {
        std::unique_ptr<RWMol> m(new RWMol(*om));
        auto hsh = MolHash(m.get(),HashFunction::ExtendedMurcko);
      CHECK(hsh == "*c1cc(C2CCCCC2*)ccn1");
    }
    {
        std::unique_ptr<RWMol> m(new RWMol(*om));
        auto hsh = MolHash(m.get(),HashFunction::MolFormula);
      CHECK(hsh == "C12H17NO2");
    }
    {
        std::unique_ptr<RWMol> m(new RWMol(*om));
        auto hsh = MolHash(m.get(),HashFunction::AtomBondCounts);
      CHECK(hsh == "15,16");
    }
    {
        std::unique_ptr<RWMol> m(new RWMol(*om));
        auto hsh = MolHash(m.get(),HashFunction::DegreeVector);
      CHECK(hsh == "0,4,9,2");
    }
    {
        std::unique_ptr<RWMol> m(new RWMol(*om));
        auto hsh = MolHash(m.get(),HashFunction::Mesomer);
      CHECK(hsh == "CO[C]1[CH][C](C2CCCCC2O)[CH][CH][N]1_0");
    }   
  }
  SECTION("tautomers") {
      auto om = "C(CC1=NNC=C1)C1=CNC=N1"_smiles;
      REQUIRE(om);
      {
          std::unique_ptr<RWMol> m(new RWMol(*om));
          auto hsh = MolHash(m.get(), HashFunction::HetAtomTautomer);
          CHECK(hsh == "[CH]1[CH][C](CC[C]2[CH][N][CH][N]2)[N][N]1_2_0");
      }
      {
          std::unique_ptr<RWMol> m(new RWMol(*om));
          auto hsh = MolHash(m.get(), HashFunction::HetAtomProtomer);
          CHECK(hsh == "[CH]1[CH][C](CC[C]2[CH][N][CH][N]2)[N][N]1_2");
      }
  }

}