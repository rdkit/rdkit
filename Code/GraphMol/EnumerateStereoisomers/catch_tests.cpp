//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <ranges>
#include <unordered_set>

#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/EnumerateStereoisomers/EnumerateStereoisomers.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::EnumerateStereoisomers;

TEST_CASE("Simple test") {
  auto m1 = "BrC=CC1OC(C2)(F)C2(Cl)C1"_smiles;
  REQUIRE(m1);
  {
    StereoisomerEnumerator enu(*m1);
    CHECK(enu.getStereoisomerCount() == 16);
    const static std::unordered_set<std::string> expected{
        R"(F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2)",
        R"(F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2)",
        R"(F[C@@]12C[C@@]1(Cl)C[C@H](/C=C/Br)O2)",
        R"(F[C@@]12C[C@@]1(Cl)C[C@H](/C=C\Br)O2)",
        R"(F[C@@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2)",
        R"(F[C@@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2)",
        R"(F[C@@]12C[C@]1(Cl)C[C@H](/C=C/Br)O2)",
        R"(F[C@@]12C[C@]1(Cl)C[C@H](/C=C\Br)O2)",
        R"(F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2)",
        R"(F[C@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2)",
        R"(F[C@]12C[C@@]1(Cl)C[C@H](/C=C/Br)O2)",
        R"(F[C@]12C[C@@]1(Cl)C[C@H](/C=C\Br)O2)",
        R"(F[C@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2)",
        R"(F[C@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2)",
        R"(F[C@]12C[C@]1(Cl)C[C@H](/C=C/Br)O2)",
        R"(F[C@]12C[C@]1(Cl)C[C@H](/C=C\Br)O2)",
    };
    std::unordered_set<std::string> got;
    std::cout << "[";
    while (auto isomer = enu.next()) {
      std::cout << "r\"" << MolToSmiles(*isomer) << "\"," << std::endl;
      got.insert(MolToSmiles(*isomer));
    }
    std::cout << "]" << std::endl;
    CHECK(got == expected);
  }
}

TEST_CASE("Embedding") {
  auto m1 = "BrC=CC1OC(C2)(F)C2(Cl)C1"_smiles;
  REQUIRE(m1);

  StereoEnumerationOptions opts;
  opts.tryEmbedding = true;
  StereoisomerEnumerator enu1(*m1, opts);
  const static std::unordered_set<std::string> expected{
      R"(F[C@@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2)",
      R"(F[C@@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2)",
      R"(F[C@@]12C[C@]1(Cl)C[C@H](/C=C/Br)O2)",
      R"(F[C@@]12C[C@]1(Cl)C[C@H](/C=C\Br)O2)",
      R"(F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2)",
      R"(F[C@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2)",
      R"(F[C@]12C[C@@]1(Cl)C[C@H](/C=C/Br)O2)",
      R"(F[C@]12C[C@@]1(Cl)C[C@H](/C=C\Br)O2)",
  };
  std::unordered_set<std::string> got;
  std::cout << "[";
  while (auto isomer = enu1.next()) {
    std::cout << "r\"" << MolToSmiles(*isomer) << "\"," << std::endl;
    got.insert(MolToSmiles(*isomer));
  }
  std::cout << "]" << std::endl;
  CHECK(got == expected);

  // Check we get the right number even with maxIsomers.
  opts.maxIsomers = 8;
  StereoisomerEnumerator enu2(*m1, opts);
  got.clear();
  while (auto isomer = enu2.next()) {
    got.insert(MolToSmiles(*isomer));
  }
  CHECK(got == expected);

  // Check no infinite loop when maxIsomers greater than
  // possible.
  opts.maxIsomers = 1024;
  StereoisomerEnumerator enu3(*m1, opts);
  got.clear();
  while (auto isomer = enu3.next()) {
    got.insert(MolToSmiles(*isomer));
  }
  CHECK(got == expected);
}

TEST_CASE("Unique") {
  auto m1 = "FC(Cl)C=CC=CC(F)Cl"_smiles;
  REQUIRE(m1);
  StereoEnumerationOptions opts;
  opts.unique = true;
  StereoisomerEnumerator enu(*m1, opts);
  const static std::unordered_set<std::string> expected{
      R"(F[C@@H](Cl)/C=C/C=C/[C@@H](F)Cl)",
      R"(F[C@@H](Cl)/C=C\C=C/[C@@H](F)Cl)",
      R"(F[C@@H](Cl)/C=C\C=C\[C@@H](F)Cl)",
      R"(F[C@H](Cl)/C=C/C=C/[C@@H](F)Cl)",
      R"(F[C@H](Cl)/C=C/C=C/[C@H](F)Cl)",
      R"(F[C@H](Cl)/C=C/C=C\[C@@H](F)Cl)",
      R"(F[C@H](Cl)/C=C\C=C/[C@@H](F)Cl)",
      R"(F[C@H](Cl)/C=C\C=C/[C@H](F)Cl)",
      R"(F[C@H](Cl)/C=C\C=C\[C@@H](F)Cl)",
      R"(F[C@H](Cl)/C=C\C=C\[C@H](F)Cl)",
  };
  std::unordered_set<std::string> got;
  std::cout << "[";
  while (auto isomer = enu.next()) {
    std::cout << "r\"" << MolToSmiles(*isomer) << "\"," << std::endl;
    got.insert(MolToSmiles(*isomer));
  }
  std::cout << "]" << std::endl;
  CHECK(got == expected);
}

TEST_CASE("Unassigned") {
  auto m1 = "C/C(F)=C/[C@@H](C)Cl"_smiles;
  REQUIRE(m1);
  StereoEnumerationOptions opts;
  StereoisomerEnumerator enu1(*m1, opts);
  CHECK(enu1.getStereoisomerCount() == 1);
  std::unordered_set<std::string> expected{"C/C(F)=C/[C@@H](C)Cl"};
  std::unordered_set<std::string> got;
  while (auto isomer = enu1.next()) {
    got.insert(MolToSmiles(*isomer));
    CHECK(!isomer->hasProp("_MolFileChiralFlag"));
  }
  CHECK(got == expected);

  // Enumerate bond stereo only
  auto m4 = "CC(F)=C[C@@H](C)Cl"_smiles;
  REQUIRE(m4);
  StereoisomerEnumerator enu4(*m4, opts);
  CHECK(enu4.getStereoisomerCount() == 2);
  got.clear();
  while (auto isomer = enu4.next()) {
    std::string prop;
    CHECK(isomer->getPropIfPresent<std::string>("_MolFileChiralFlag", prop));
    CHECK(prop == "1");
    std::cout << "r\"" << MolToSmiles(*isomer) << "\"," << std::endl;
    got.insert(MolToSmiles(*isomer));
  }
  expected = std::unordered_set<std::string>{
      R"(C/C(F)=C/[C@@H](C)Cl)",
      R"(C/C(F)=C\[C@@H](C)Cl)",
  };
  CHECK(got == expected);

  auto m2 = "BrC=C[C@H]1OC(C2)(F)C2(Cl)C1"_smiles;
  REQUIRE(m2);
  StereoisomerEnumerator enu2(*m2, opts);
  CHECK(enu2.getStereoisomerCount() == 8);
  expected = std::unordered_set<std::string>{
      R"(F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2)",
      R"(F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2)",
      R"(F[C@@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2)",
      R"(F[C@@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2)",
      R"(F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2)",
      R"(F[C@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2)",
      R"(F[C@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2)",
      R"(F[C@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2)",
  };
  std::cout << "[";
  got.clear();
  while (auto isomer = enu2.next()) {
    std::cout << "r\"" << MolToSmiles(*isomer) << "\"," << std::endl;
    got.insert(MolToSmiles(*isomer));
  }
  std::cout << "]" << std::endl;
  CHECK(got == expected);

  opts.onlyUnassigned = false;
  StereoisomerEnumerator enu3(*m2, opts);
  CHECK(enu3.getStereoisomerCount() == 16);
}

TEST_CASE("Subset") {
  std::string smi("Br");
  for (int i = 0; i < 20; ++i) {
    smi += "[CH](Cl)";
  }
  smi += "F";
  auto m1 = v2::SmilesParse::MolFromSmiles(smi);
  StereoEnumerationOptions opts;
  opts.randomSeed = 1066;
  StereoisomerEnumerator enu1(*m1, opts);
  CHECK(enu1.getStereoisomerCount() == 1048576);
  const static std::unordered_set<std::string> expected1{
      R"(F[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)Br)",
      R"(F[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)Br)",
      R"(F[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)Br)",
      R"(F[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)Br)",
      R"(F[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)Br)",
  };
  std::unordered_set<std::string> got;
  std::cout << "[";
  for (int i = 0; i < 5; ++i) {
    auto isomer = enu1.next();
    std::cout << "r\"" << MolToSmiles(*isomer) << "\"," << std::endl;
    got.insert(MolToSmiles(*isomer));
  }
  std::cout << "]" << std::endl;
  CHECK(got == expected1);

  opts.maxIsomers = 3;
  StereoisomerEnumerator enu2(*m1, opts);
  const static std::unordered_set<std::string> expected2{
      R"(F[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)Br)",
      R"(F[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)Br)",
      R"(F[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)Br)"};
  got.clear();
  std::cout << "[";
  while (auto isomer = enu2.next()) {
    std::cout << "r\"" << MolToSmiles(*isomer) << "\"," << std::endl;
    got.insert(MolToSmiles(*isomer));
  }
  std::cout << "]" << std::endl;
  CHECK(got == expected2);
}

TEST_CASE("EnhancedStereo") {
  const char *rdbase = getenv("RDBASE");
  REQUIRE(rdbase);
  std::string fname = std::string(rdbase) +
                      "/Code/GraphMol/FileParsers/test_data/two_centers_or.mol";
  auto m1 = v2::FileParsers::MolFromMolFile(fname);
  REQUIRE(m1);
  StereoisomerEnumerator enu1(*m1);
  const static std::unordered_set<std::string> expected{
      R"(C[C@H]([C@@H](C)F)[C@@H](C)Br)", R"(C[C@@H]([C@H](C)F)[C@@H](C)Br)"};
  std::unordered_set<std::string> got;
  while (auto isomer = enu1.next()) {
    std::string prop;
    CHECK(isomer->getPropIfPresent<std::string>("_MolFileChiralFlag", prop));
    CHECK(prop == "1");
    std::cout << "r\"" << MolToSmiles(*isomer) << "\"," << std::endl;
    got.insert(MolToSmiles(*isomer));
  }
  CHECK(got == expected);

  // Check the original was in there
  auto m2 = RWMol(*m1);
  m2.setStereoGroups(std::vector<StereoGroup>());
  auto origSmi = MolToSmiles(m2);
  auto it = got.find(origSmi);
  CHECK(it != got.end());

  // Make sure there are no extras if enumerating all
  StereoEnumerationOptions opts;
  opts.unique = false;
  opts.onlyUnassigned = false;
  got.clear();
  StereoisomerEnumerator enu2(*m1, opts);
  while (auto isomer = enu2.next()) {
    std::string prop;
    CHECK(isomer->getPropIfPresent<std::string>("_MolFileChiralFlag", prop));
    CHECK(prop == "1");
    std::cout << "r\"" << MolToSmiles(*isomer) << "\"," << std::endl;
    CHECK(got.insert(MolToSmiles(*isomer)).second);
  }
  CHECK(got.size() == 8);
}

TEST_CASE("Issue 2890") {
  auto m1 = "CC=CC"_smiles;
  m1->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOANY);
  StereoisomerEnumerator enu1(*m1);
  std::unordered_set<std::string> got;
  while (auto isomer = enu1.next()) {
    std::string prop;
    CHECK(isomer->getPropIfPresent<std::string>("_MolFileChiralFlag", prop));
    CHECK(prop == "1");
    got.insert(MolToSmiles(*isomer));
  }
  CHECK(got.size() == 2);
}

TEST_CASE("Issue 3231") {
  // EnumerateStereochemistry should clear CIP labels.  The original test
  // uses FindMolChiralCenters but that's Python only, too.  There are only
  // so may yaks that can be shaved in one go!
  auto getCIPLabels =
      [](ROMol &mol) -> std::vector<std::pair<unsigned int, std::string>> {
    MolOps::assignStereochemistry(mol);
    CIPLabeler::assignCIPLabels(mol);
    std::vector<std::pair<unsigned int, std::string>> result;
    for (auto &si : Chirality::findPotentialStereo(mol, true, true)) {
      if (si.type == Chirality::StereoType::Atom_Tetrahedral) {
        if (si.specified == Chirality::StereoSpecified::Unspecified) {
          result.push_back(std::make_pair(si.centeredOn, "?"));
        } else {
          auto atom = mol.getAtomWithIdx(si.centeredOn);
          auto cipLabel =
              atom->getProp<std::string>(common_properties::_CIPCode);
          result.push_back(std::make_pair(si.centeredOn, cipLabel));
        }
      }
    }
    return result;
  };

  auto m1 =
      "C[C@H](OC1=C(N)N=CC(C2=CN(C3C[C@H](C)NCC3)N=C2)=C1)C4=C(Cl)C=CC(F)=C4Cl"_smiles;
  REQUIRE(m1);
  auto origLabels = getCIPLabels(*m1);
  CHECK(origLabels == std::vector<std::pair<unsigned int, std::string>>{
                          {1, "S"}, {12, "?"}, {14, "S"}});
  StereoEnumerationOptions opts;
  opts.onlyUnassigned = false;
  opts.maxIsomers = 20;
  StereoisomerEnumerator enu(*m1, opts);
  std::vector<std::vector<std::pair<unsigned int, std::string>>> chiCents;
  while (auto isomer = enu.next()) {
    chiCents.push_back(getCIPLabels(*isomer));
  }
  std::ranges::sort(chiCents);
  chiCents.erase(std::unique(chiCents.begin(), chiCents.end()), chiCents.end());
  for (const auto &cc : chiCents) {
    for (const auto &label : cc) {
      std::cout << label.first << " " << label.second << " :: ";
    }
    std::cout << std::endl;
  }
  std::vector<std::vector<std::pair<unsigned int, std::string>>> expCents{
      {{1, "R"}, {12, "R"}, {14, "R"}}, {{1, "R"}, {12, "R"}, {14, "S"}},
      {{1, "R"}, {12, "S"}, {14, "R"}}, {{1, "R"}, {12, "S"}, {14, "S"}},
      {{1, "S"}, {12, "R"}, {14, "R"}}, {{1, "S"}, {12, "R"}, {14, "S"}},
      {{1, "S"}, {12, "S"}, {14, "R"}}, {{1, "S"}, {12, "S"}, {14, "S"}}};
  CHECK(chiCents.size() == expCents.size());
}

TEST_CASE("Issue 3505") {
  auto m1 = "CCC(C)Br"_smiles;
  REQUIRE(m1);
  StereoisomerEnumerator enu1(*m1);
  int count = 0;
  while (auto isomer = enu1.next()) {
    std::string prop;
    CHECK(isomer->getPropIfPresent<std::string>("_MolFileChiralFlag", prop));
    CHECK(prop == "1");
    auto at = isomer->getAtomWithIdx(2);
    CHECK((at->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CW ||
           at->getChiralTag() == Atom::ChiralType::CHI_TETRAHEDRAL_CCW));
    CHECK(at->hasProp("_ChiralityPossible"));
    ++count;
  }
  CHECK(count == 2);
}

TEST_CASE("Either or Double Stereo") {
  const char *rdbase = getenv("RDBASE");
  REQUIRE(rdbase);
  std::string fname = std::string(rdbase) +
                      "/Code/GraphMol/FileParsers/test_data/simple_either.mol";
  auto m1 = v2::FileParsers::MolFromMolFile(fname);
  REQUIRE(m1);
  StereoisomerEnumerator enu1(*m1);
  std::unordered_set<std::string> got;
  while (auto isomer = enu1.next()) {
    std::string prop;
    CHECK(isomer->getPropIfPresent<std::string>("_MolFileChiralFlag", prop));
    CHECK(prop == "1");
    got.insert(MolToSmiles(*isomer));
  }
  CHECK(got == std::unordered_set<std::string>{R"(C/C=C/C)", R"(C/C=C\C)"});
}

TEST_CASE("Embed many chirals") {
  std::string smi("C1");
  for (int i = 0; i < 40; ++i) {
    smi += "C(Cl)(Br)";
  }
  smi += "C1";
  auto m1 = v2::SmilesParse::MolFromSmiles(smi);
  REQUIRE(m1);
  StereoEnumerationOptions opts;
  opts.maxIsomers = 2;
  opts.tryEmbedding = true;
  StereoisomerEnumerator enu1(*m1, opts);
  std::unordered_set<std::string> got;
  while (auto isomer = enu1.next()) {
    got.insert(MolToSmiles(*isomer));
  }
  CHECK(got.size() == 2);
}

TEST_CASE("Issue 6405") {
  auto m1 = "O[C@H](Br)[C@H](F)C |&1:1,3|"_smiles;
  REQUIRE(m1);
  StereoisomerEnumerator enu1(*m1);
  while (auto isomer = enu1.next()) {
    std::string prop;
    CHECK(isomer->getPropIfPresent<std::string>("_MolFileChiralFlag", prop));
    CHECK(prop == "1");
    CHECK(isomer->getStereoGroups().empty());
  }
}

TEST_CASE("Issue 7516") {
  auto m1 = "CC1CC(C)C1"_smiles;
  REQUIRE(m1);
  StereoisomerEnumerator enu1(*m1);
  std::unordered_set<std::string> got;
  while (auto isomer = enu1.next()) {
    got.insert(MolToSmiles(*isomer));
  }
  CHECK(got == std::unordered_set<std::string>{R"(C[C@H]1C[C@@H](C)C1)",
                                               R"(C[C@H]1C[C@H](C)C1)"});

  auto m2 = "COC(=O)C1CC(NC(N)=O)C1"_smiles;
  REQUIRE(m2);
  StereoisomerEnumerator enu2(*m2);
  got.clear();
  while (auto isomer = enu2.next()) {
    got.insert(MolToSmiles(*isomer));
  }
  CHECK(got.size() == 2);

  auto m3 = "O=C(NC1CC2[NH+](C(C1)CC2)Cc3ccccc3)N"_smiles;
  REQUIRE(m2);
  StereoisomerEnumerator enu3(*m3);
  got.clear();
  while (auto isomer = enu3.next()) {
    got.insert(MolToSmiles(*isomer));
  }
  CHECK(got.size() == 8);
}

TEST_CASE("Issue 7608") {
  auto m1 = "N=C(N1CCC1)NCCOc2ccc(C(F)(F)F)cc2"_smiles;
  REQUIRE(m1);
  StereoisomerEnumerator enu1(*m1);
  std::unordered_set<std::string> got;
  while (auto isomer = enu1.next()) {
    got.insert(MolToSmiles(*isomer));
  }
  CHECK(got.size() == 1);
}

TEST_CASE("Atropisomers") {
  auto m1 = "CC1=CC=CC(I)=C1N1C(C)=CC=C1Br |wD:7.7,wU:8.9|"_smiles;
  std::cout << "Starting " << MolToCXSmiles(*m1) << std::endl;
  REQUIRE(m1);
  StereoEnumerationOptions opts;
  opts.onlyUnassigned = false;
  StereoisomerEnumerator enu1(*m1, opts);
  std::unordered_set<std::string> got;
  while (auto isomer = enu1.next()) {
    std::cout << "Isomer : "
              << MolToCXSmiles(*isomer, SmilesWriteParams(),
                               SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS)
              << std::endl;
    got.insert(MolToCXSmiles(*isomer, SmilesWriteParams(),
                             SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS));
  }
  CHECK(got.size() == 2);

  auto m2 =
      "C1=C(C2=CCCCCCCCCCC2)CCC1 |(-1.71835,-1.22732,;-1.61655,-0.232318,;-0.751953,0.270082,;-0.754753,1.27008,;0.109847,1.77248,;0.977247,1.27488,;1.84185,1.77728,;2.70925,1.27968,;2.71205,0.279682,;1.84745,-0.222718,;1.85025,-1.22272,;0.985647,-1.72512,;0.118247,-1.22752,;0.115447,-0.227518,;-2.53135,0.171882,;-3.19835,-0.573118,;-2.69595,-1.43772,),wU:1.14|"_smiles;
  REQUIRE(m1);
  StereoisomerEnumerator enu2(*m2, opts);
  got.clear();
  while (auto isomer = enu2.next()) {
    std::cout << "Isomer : "
              << MolToCXSmiles(*isomer, SmilesWriteParams(),
                               SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS)
              << std::endl;
    got.insert(MolToCXSmiles(*isomer, SmilesWriteParams(),
                             SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS));
  }
  for (auto &g : got) {
    std::cout << "r\"" << g << "\"," << std::endl;
  }
  CHECK(got.size() == 4);

  {
    auto m = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 20 21 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -0.024300 -1.051200 0.000000 0
M  V30 2 C 0.690200 -1.463700 0.000000 0
M  V30 3 O 0.690200 -2.288700 0.000000 0
M  V30 4 N 1.404500 -1.051200 0.000000 0
M  V30 5 C 2.119100 -1.463700 0.000000 0
M  V30 6 C 2.119100 -2.288700 0.000000 0
M  V30 7 C 1.404500 -2.701200 0.000000 0
M  V30 8 C 2.833500 -2.701200 0.000000 0
M  V30 9 C 3.547900 -2.288700 0.000000 0
M  V30 10 N 3.547900 -1.463700 0.000000 0
M  V30 11 C 2.833600 -1.051200 0.000000 0
M  V30 12 C 3.363800 -0.419100 0.000000 0
M  V30 13 C 4.176300 -0.562300 0.000000 0
M  V30 14 C 3.081600 0.356100 0.000000 0
M  V30 15 C 1.404500 -0.226000 0.000000 0
M  V30 16 N 2.155300 0.191900 0.000000 0
M  V30 17 C 2.161900 1.051300 0.000000 0
M  V30 18 C 1.417800 1.480900 0.000000 0
M  V30 19 C 0.676800 1.045500 0.000000 0
M  V30 20 C 0.690200 0.186400 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 4 2
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 1 6 8
M  V30 8 2 8 9
M  V30 9 1 9 10
M  V30 10 2 10 11
M  V30 11 1 11 12
M  V30 12 1 12 13
M  V30 13 1 12 14
M  V30 14 1 4 15
M  V30 15 2 15 16
M  V30 16 1 16 17
M  V30 17 2 17 18
M  V30 18 1 18 19
M  V30 19 2 19 20
M  V30 20 1 5 11 CFG=1
M  V30 21 1 20 15
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    StereoEnumerationOptions opts;
    opts.onlyUnassigned = false;
    StereoisomerEnumerator enu1(*m, opts);
    std::unordered_set<std::string> got;
    while (auto isomer = enu1.next()) {
      std::cout << "Isomer : "
                << MolToCXSmiles(*isomer, SmilesWriteParams(),
                                 SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS)
                << std::endl;
      got.insert(MolToCXSmiles(*isomer, SmilesWriteParams(),
                               SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS));
    }
    for (auto &g : got) {
      std::cout << "r\"" << g << "\"," << std::endl;
    }

    CHECK(got.size() == 2);
  }
}