//
//  Copyright (C) 2022-2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include "AlignMolecules.h"
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

using namespace RDKit;

TEST_CASE("symmetric functional groups") {
  SECTION("basics") {
    auto m1 =
        "CCC(=O)[O-] "
        "|(-1.11,0.08,-0.29;0.08,-0.18,0.58;1.34,0.03,-0.16;1.74,1.22,-0.32;2.06,-1.04,-0.66)|"_smiles;
    REQUIRE(m1);

    // swap the bond orders to the Os
    RWMol m2(*m1);
    m2.getAtomWithIdx(3)->setFormalCharge(-1);
    m2.getAtomWithIdx(4)->setFormalCharge(0);
    m2.getBondBetweenAtoms(2, 3)->setBondType(Bond::BondType::SINGLE);
    m2.getBondBetweenAtoms(2, 4)->setBondType(Bond::BondType::DOUBLE);

    {
      auto rmsd = MolAlign::getBestRMS(m2, *m1);
      CHECK(rmsd == Catch::Approx(0.0).margin(1e-3));
    }
    {
      // previous behavior
      int probeId = -1;
      int refId = -1;
      std::vector<MatchVectType> mp;
      int maxMatches = 1e6;
      bool symmetrize = false;
      auto rmsd = MolAlign::getBestRMS(m2, *m1, probeId, refId, mp, maxMatches,
                                       symmetrize);
      CHECK(rmsd == Catch::Approx(0.747).margin(1e-3));
    }
  }
  SECTION("terminal sulfate1") {
    auto m1 =
        "CS(=O)(=O)[O-] |(-0.93,-0.06,-0.04;0.82,0.07,0.13;1.27,-0.04,1.54;1.21,1.40,-0.48;1.53,-1.11,-0.82)|"_smiles;
    REQUIRE(m1);

    // swap the bond orders to the Os
    RWMol m2(*m1);
    m2.getAtomWithIdx(2)->setFormalCharge(-1);
    m2.getAtomWithIdx(4)->setFormalCharge(0);
    m2.getBondBetweenAtoms(1, 2)->setBondType(Bond::BondType::SINGLE);
    m2.getBondBetweenAtoms(1, 4)->setBondType(Bond::BondType::DOUBLE);

    {
      auto rmsd = MolAlign::getBestRMS(m2, *m1);
      CHECK(rmsd == Catch::Approx(0.0).margin(1e-3));
    }
    {
      // previous behavior
      int probeId = -1;
      int refId = -1;
      std::vector<MatchVectType> mp;
      int maxMatches = 1e6;
      bool symmetrize = false;
      auto rmsd = MolAlign::getBestRMS(m2, *m1, probeId, refId, mp, maxMatches,
                                       symmetrize);
      CHECK(rmsd == Catch::Approx(0.097).margin(1e-3));
    }
  }
  SECTION("terminal sulfate2") {
    auto m1 =
        "CS(=O)(=O)[O-] |(-0.93,-0.06,-0.04;0.82,0.07,0.13;1.27,-0.04,1.54;1.21,1.40,-0.48;1.53,-1.11,-0.82)|"_smiles;
    REQUIRE(m1);

    // swap the bond orders to the Os
    RWMol m2(*m1);
    m2.getAtomWithIdx(3)->setFormalCharge(-1);
    m2.getAtomWithIdx(4)->setFormalCharge(0);
    m2.getBondBetweenAtoms(1, 3)->setBondType(Bond::BondType::SINGLE);
    m2.getBondBetweenAtoms(1, 4)->setBondType(Bond::BondType::DOUBLE);

    {
      auto rmsd = MolAlign::getBestRMS(m2, *m1);
      CHECK(rmsd == Catch::Approx(0.0).margin(1e-3));
    }
    {
      // previous behavior
      int probeId = -1;
      int refId = -1;
      std::vector<MatchVectType> mp;
      int maxMatches = 1e6;
      bool symmetrize = false;
      auto rmsd = MolAlign::getBestRMS(m2, *m1, probeId, refId, mp, maxMatches,
                                       symmetrize);
      CHECK(rmsd == Catch::Approx(0.097).margin(1e-3));
    }
  }
}

#ifdef RDK_BUILD_THREADSAFE_SSS
TEST_CASE("multithreaded getBestRMS") {
  SECTION("basics") {
    // has 288 self matches
    auto m1 =
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)F |(-1.17097,1.42189,1.14513;-1.54917,0.262724,0.549205;-1.70317,-0.764739,1.49745;-2.82875,0.445186,-0.0104401;-0.695326,-0.20819,-0.58675;-1.32875,-1.40402,-1.02194;-0.794122,0.733556,-1.61075;0.695194,-0.600926,-0.295382;1.26585,-1.00432,-1.52316;0.671971,-1.69988,0.563393;1.62838,0.438987,0.231506;1.1944,0.938004,1.42313;1.6862,1.50141,-0.682411;2.92826,-0.0596757,0.321018)|"_smiles;
    REQUIRE(m1);
    auto m2 =
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)F |(-1.4374,1.69863,0.0454955;-1.63421,0.384221,0.267923;-1.76372,0.215075,1.6529;-2.89568,0.031803,-0.240695;-0.665291,-0.57679,-0.303489;-1.12767,-1.85176,0.0867448;-0.706055,-0.620221,-1.69542;0.736642,-0.522726,0.120365;0.785974,-0.733214,1.50576;1.36074,-1.67057,-0.401737;1.57617,0.630117,-0.22085;1.06853,1.80237,0.30467;1.8928,0.75494,-1.55471;2.80916,0.458118,0.433039)|"_smiles;
    REQUIRE(m2);
    auto probeId = 0;
    auto refId = 0;
    std::vector<MatchVectType> mp;
    int maxMatches = 1e6;
    bool symmetrize = true;
    RDNumeric::DoubleVector *weights = nullptr;
    int numThreads = 1;
    auto ref = MolAlign::getBestRMS(*m2, *m1, probeId, refId, mp, maxMatches,
                                    symmetrize, weights, numThreads);
    numThreads = 4;
    auto mt_val = MolAlign::getBestRMS(*m2, *m1, probeId, refId, mp, maxMatches,
                                       symmetrize, weights, numThreads);
    CHECK(ref == Catch::Approx(mt_val).epsilon(0.00001));
  }
  SECTION("more symmetry") {
    std::string rdbase = getenv("RDBASE");
    std::string fname1 =
        rdbase + "/Code/GraphMol/MolAlign/test_data/symmetric.mol";
    std::unique_ptr<ROMol> m1{MolFileToMol(fname1)};
    REQUIRE(m1);

    auto probeId = 0;
    auto refId = 0;
    std::vector<MatchVectType> mp;
    int maxMatches = 1e6;
    bool symmetrize = true;
    RDNumeric::DoubleVector *weights = nullptr;
    {
      int numThreads = 1;
      auto start = std::chrono::high_resolution_clock::now();
      auto ref = MolAlign::getBestRMS(*m1, *m1, probeId, refId, mp, maxMatches,
                                      symmetrize, weights, numThreads);
      auto finish = std::chrono::high_resolution_clock::now();
      std::cerr << (finish - start).count() << std::endl;
      CHECK(ref == Catch::Approx(0.0).epsilon(0.00001));
    }
    {
      int numThreads = 4;
      auto start = std::chrono::high_resolution_clock::now();
      auto ref = MolAlign::getBestRMS(*m1, *m1, probeId, refId, mp, maxMatches,
                                      symmetrize, weights, numThreads);
      auto finish = std::chrono::high_resolution_clock::now();
      std::cerr << (finish - start).count() << std::endl;
      CHECK(ref == Catch::Approx(0.0).epsilon(0.00001));
    }
  }
}
#endif

TEST_CASE("getAllConformerBestRMS") {
  std::string rdbase = getenv("RDBASE");
  std::string fname1 =
      rdbase + "/Code/GraphMol/MolAlign/test_data/symmetric.confs.sdf";
  SDMolSupplier suppl(fname1);
  std::unique_ptr<ROMol> mol{suppl[0]};
  REQUIRE(mol);
  for (auto i = 1u; i < suppl.length(); ++i) {
    std::unique_ptr<ROMol> nm{suppl[i]};
    REQUIRE(nm);
    mol->addConformer(new Conformer(nm->getConformer()), true);
  }
  // CHECK(mol->getNumConformers() == 10);
  SECTION("basics") {
    auto nconfs = mol->getNumConformers();
    std::vector<double> rmsds;
    {
      auto start = std::chrono::high_resolution_clock::now();
      rmsds = MolAlign::getAllConformerBestRMS(*mol);
      CHECK(rmsds.size() == (nconfs * (nconfs - 1)) / 2);
      auto finish = std::chrono::high_resolution_clock::now();
      std::cerr << (finish - start).count() << std::endl;

      ROMol refMol(*mol);
      ROMol prbMol(*mol);
      auto refVal = MolAlign::getBestRMS(refMol, prbMol, 1, 0);
      CHECK(rmsds[0] == Catch::Approx(refVal).epsilon(0.00001));
    }
    std::vector<double> mtrmsds;
    {
      auto start = std::chrono::high_resolution_clock::now();
      int numThreads = 4;
      mtrmsds = MolAlign::getAllConformerBestRMS(*mol, numThreads);
      CHECK(mtrmsds.size() == (nconfs * (nconfs - 1)) / 2);
      auto finish = std::chrono::high_resolution_clock::now();
      std::cerr << (finish - start).count() << std::endl;
    }
    for (auto i = 0u; i < rmsds.size(); ++i) {
      CHECK(rmsds[i] == Catch::Approx(mtrmsds[i]).epsilon(0.00001));
    }
  }
}