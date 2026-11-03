//
//  Copyright (c) 2019-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
///
#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/test_fixtures.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <algorithm>
#include <regex>

using namespace RDKit;
using std::unique_ptr;

TEST_CASE("Github #1039") {
  const auto useLegacy = GENERATE(true, false);
  UseLegacyStereoPerceptionFixture fx(useLegacy);
  SECTION("double bond") {
    auto m1 = "C/C=C/C=C/C"_smiles;
    REQUIRE(m1);
    std::vector<unsigned int> bonds = {2};
    std::unique_ptr<ROMol> pieces(MolFragmenter::fragmentOnBonds(*m1, bonds));
    REQUIRE(pieces);
    CHECK(pieces->getNumAtoms() == 8);
    REQUIRE(pieces->getBondBetweenAtoms(3, 6));
    REQUIRE(pieces->getBondBetweenAtoms(2, 7));
    CHECK(pieces->getBondBetweenAtoms(3, 6)->getBondType() == Bond::SINGLE);
    CHECK(pieces->getBondBetweenAtoms(3, 6)->getBondDir() == Bond::ENDUPRIGHT);
    CHECK(pieces->getBondBetweenAtoms(2, 7)->getBondType() == Bond::SINGLE);
    CHECK(pieces->getBondBetweenAtoms(2, 7)->getBondDir() == Bond::ENDUPRIGHT);
    CHECK(MolToSmiles(*pieces) == "[2*]/C=C/C.[3*]/C=C/C");
  }
  SECTION("atomic stereo") {
    auto m1 = "C(C)(F)(Cl)O"_smiles;
    REQUIRE(m1);
    m1->getBondWithIdx(0)->setBondDir(Bond::BEGINWEDGE);
    std::vector<unsigned int> bonds = {0};
    std::unique_ptr<ROMol> pieces(MolFragmenter::fragmentOnBonds(*m1, bonds));
    REQUIRE(pieces);
    CHECK(pieces->getNumAtoms() == 7);
    REQUIRE(pieces->getBondBetweenAtoms(0, 6));
    REQUIRE(pieces->getBondBetweenAtoms(1, 5));
    CHECK(pieces->getBondBetweenAtoms(0, 6)->getBondDir() == Bond::BEGINWEDGE);
    CHECK(pieces->getBondBetweenAtoms(1, 5)->getBondDir() == Bond::NONE);
    // no actual stereo in the SMILES here since we haven't assigned it (need a
    // conformer to do that using wedging)
    CHECK(MolToSmiles(*pieces) == "*C.[1*]C(O)(F)Cl");
  }
  SECTION("bond stereo") {
    auto m = "O/C=N/C=C"_smiles;
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{{1, 1}};
    std::vector<unsigned int> bonds{0};
    auto resa = RDKit::MolFragmenter::fragmentOnBonds(*m, bonds);
    CHECK(MolToSmiles(*resa) == "*/C=N/C=C.[1*]O");
    // make sure we still have stereo atoms
    std::vector<std::vector<int>> expected_stereo_atoms{
        {5, 3},  // 5 is the new dummy atom, it was 0 before
        {},     {}, {}, {},
    };
    std::vector<std::vector<int>> received_stereo;
    for (auto *bond : resa->bonds()) {
      received_stereo.push_back(bond->getStereoAtoms());
    }
    CHECK(received_stereo == expected_stereo_atoms);
    delete resa;
  }
  {  // break non stereo atom bond
    auto m = "C/C(O)=N/C=C"_smiles;
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{{1, 1}};
    std::vector<unsigned int> bonds{0};
    auto resa = RDKit::MolFragmenter::fragmentOnBonds(*m, bonds);
    CHECK(MolToSmiles(*resa) == "*/C(O)=N/C=C.[1*]C");
    // make sure we still have stereo atoms
    std::vector<std::vector<int>> expected_stereo_atoms;
    if (useLegacy) {
      expected_stereo_atoms = {{}, {2, 4}, {}, {}, {}, {}};
    } else {
      expected_stereo_atoms = {{}, {6, 4}, {}, {}, {}, {}};
    }
    std::vector<std::vector<int>> received_stereo;
    for (auto *bond : resa->bonds()) {
      received_stereo.push_back(bond->getStereoAtoms());
    }
    CHECK(received_stereo == expected_stereo_atoms);
    delete resa;
  }
  {  // bond stereo should only be removed when deleting the double bond with
     // E/Z
    auto m = "O/C=N/C=C"_smiles;
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{{1, 1}};
    std::vector<std::string> expected = {
        "*/C=N/C=C.[1*]O",
        "[1*]=NC=C.[2*]=CO",  // bond stereo gone
        "[2*]C=C.[3*]/N=C/O", "[3*]=C.[4*]=C/N=C/O"};
    for (unsigned int i = 0; i < m->getNumBonds(); ++i) {
      std::vector<unsigned int> bonds{i};
      auto resa = RDKit::MolFragmenter::fragmentOnBonds(*m, bonds);
      auto smiles = MolToSmiles(*resa);
      CHECK(smiles == expected[i]);
      delete resa;
    }
  }
  {  // bond stereo should only be removed when deleting the double bond with
     // E/Z
    // chiral stereo should stay
    auto m = "O/C=N/[C@H](I)F"_smiles;
    std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{{1, 1}};
    std::vector<std::string> expected = {
        "*/C=N/[C@@H](F)I.[1*]O",
        "[1*]=N[C@@H](F)I.[2*]=CO",  // bond stereo gone
        "[2*][C@@H](F)I.[3*]/N=C/O", "[3*]I.[4*][C@H](F)/N=C/O",
        "[3*]F.[5*][C@@H](I)/N=C/O"};
    for (unsigned int i = 0; i < m->getNumBonds(); ++i) {
      std::vector<unsigned int> bonds{i};
      auto resa = RDKit::MolFragmenter::fragmentOnBonds(*m, bonds);
      auto smiles = MolToSmiles(*resa);
      CHECK(smiles == expected[i]);
      delete resa;
    }
  }
}

TEST_CASE("molzip") {
  const auto useLegacy = GENERATE(true, false);
  CAPTURE(useLegacy);
  UseLegacyStereoPerceptionFixture fx(useLegacy);
  SECTION("basic tests") {
    {
      auto a = "C[*:1]"_smiles;
      auto b = "N[*:1]"_smiles;
      auto mol = molzip(*a, *b);
      CHECK(MolToSmiles(*mol) == "CN");
    }
    {
      // 0 isotopes aren't mapped
      auto a = "C[*]"_smiles;
      auto b = "N[*]"_smiles;
      MolzipParams p;
      auto mol = molzip(*a, *b, p);
      CHECK(MolToSmiles(*mol) == "*C.*N");
    }
    {
      auto a = "C[1*]"_smiles;
      auto b = "N[1*]"_smiles;
      MolzipParams p;
      p.label = MolzipLabel::Isotope;
      auto mol = molzip(*a, *b, p);
      CHECK(MolToSmiles(*mol) == "CN");
    }
    {
      auto a = "[C@H](Br)([*:1])F"_smiles;
      auto b = "[*:1]N"_smiles;
      auto mol = molzip(*a, *b);
      CHECK(MolToSmiles(*mol) == "N[C@@H](F)Br");
    }
    {
      auto a = "[C@H]([*:1])(Br)F"_smiles;
      auto b = "[*:1]N"_smiles;
      auto mol = molzip(*a, *b);
      CHECK(MolToSmiles(*mol) == "N[C@H](F)Br");
    }

    {
      auto a = "[C@H]([*:1])(F)([*:2])"_smiles;
      auto b = "[*:1]N.[*:2]I"_smiles;
      auto mol = molzip(*a, *b);
      REQUIRE(MolToSmiles(*mol) == "N[C@@H](F)I");
    }

    {
      auto a = "[C@H]([Xe])(F)([V])"_smiles;
      auto b = "[Xe]N.[V]I"_smiles;
      MolzipParams params;
      params.label = MolzipLabel::AtomType;
      params.atomSymbols = {"Xe", "V"};
      auto mol = molzip(*a, *b, params);
      CHECK(MolToSmiles(*mol) == "N[C@@H](F)I");
    }

    {
      auto m = "OOO[C@](F)(I)N"_smiles;
      std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{{1, 1},
                                                                     {2, 2}};
      for (unsigned int i = 0; i < m->getNumBonds(); ++i) {
        for (unsigned int j = 0; j < m->getNumBonds(); ++j) {
          if (i != j) {
            std::vector<unsigned int> bonds{i, j};
            auto resa = RDKit::MolFragmenter::fragmentOnBonds(*m, bonds);
            MolzipParams p;
            p.label = MolzipLabel::FragmentOnBonds;
            CHECK(MolToSmiles(*molzip(*resa, p)) == MolToSmiles(*m));
            delete resa;
            // Now try using atom labels
            auto res = RDKit::MolFragmenter::fragmentOnBonds(*m, bonds, true,
                                                             &dummyLabels);
            for (auto *atom : res->atoms()) {
              if (atom->getIsotope()) {
                atom->setAtomMapNum(atom->getIsotope());
              }
            }
            CHECK(MolToSmiles(*molzip(*res)) == MolToSmiles(*m));
            delete res;
          }
        }
      }
    }
  }
  SECTION("use atom property as label") {
    auto a = "[C@H]([*])(F)([*])"_smiles;
    auto b = "[*]N.[*]I"_smiles;
    a->getAtomWithIdx(1)->setProp<unsigned int>("foo", 1);
    a->getAtomWithIdx(3)->setProp<unsigned int>("foo", 2);
    b->getAtomWithIdx(0)->setProp<unsigned int>("foo", 1);
    b->getAtomWithIdx(2)->setProp<unsigned int>("foo", 2);
    MolzipParams p;
    p.label = MolzipLabel::AtomProperty;
    p.atomProperty = "foo";
    auto mol = molzip(*a, *b, p);
    // chirality is "lost" here because [C@H]([*])(F)([*]) is considered
    // achiral
    CHECK(MolToSmiles(*mol) == "NC(F)I");
  }

  SECTION("test bond stereo") {
    {
      auto a = "F/C=C/[*:1]"_smiles;
      auto b = "[*:1]F"_smiles;
      auto mol = molzip(*a, *b);
      CHECK(MolToSmiles(*mol) == "F/C=C/F");
    }

    {
      auto a = "C=C/N=C/[*:1]"_smiles;
      auto b = "O[*:1]"_smiles;
      auto mol = molzip(*a, *b);
      CHECK(MolToSmiles(*mol) == "C=C/N=C/O");
    }

    {
      auto a = "C=C[*:1].O/C=N/[*:1]"_smiles;
      auto mol = molzip(*a);
      CHECK(MolToSmiles(*mol) == "C=C/N=C/O");
    }
    {  // test single mol isotope labels
      auto a = "C=C[1*].O/C=N/[1*]"_smiles;
      MolzipParams p;
      p.label = MolzipLabel::Isotope;
      auto mol = molzip(*a, p);
      CHECK(MolToSmiles(*mol) == "C=C/N=C/O");
    }
    {
      // double bondd stereo not handled
      // auto m =  "O/C=N/C=C/F"_smiles;
      auto m = "O/C=N/C=C"_smiles;
      std::vector<std::pair<unsigned int, unsigned int>> dummyLabels{{1, 1}};
      for (unsigned int i = 0; i < m->getNumBonds(); ++i) {
        std::vector<unsigned int> bonds{i};
        {
          std::unique_ptr<ROMol> resa{
              RDKit::MolFragmenter::fragmentOnBonds(*m, bonds)};
          auto smiles = MolToSmiles(*resa);
          if (std::count(smiles.begin(), smiles.end(), '/') != 2) {
            continue;  // we removed bond stereo in fragment to bonds!
          }
          MolzipParams p;
          p.label = MolzipLabel::FragmentOnBonds;
          CHECK(MolToSmiles(*molzip(*resa, p)) == MolToSmiles(*m));
        }
        {
          // Now try using atom labels
          std::unique_ptr<ROMol> res{RDKit::MolFragmenter::fragmentOnBonds(
              *m, bonds, true, &dummyLabels)};
          auto smiles = MolToSmiles(*res);

          if (std::count(smiles.begin(), smiles.end(), '/') != 2) {
            continue;  // we removed bond stereo in fragment to bonds!
          }
          for (auto *atom : res->atoms()) {
            if (atom->getIsotope()) {
              atom->setAtomMapNum(atom->getIsotope());
            }
          }

          CHECK(MolToSmiles(*molzip(*res)) == MolToSmiles(*m));
        }
      }
    }
  }
  SECTION("unzippable molecules") {
    auto a = "C[*:1]"_smiles;
    auto b = "N[*:2]"_smiles;
    auto mol = molzip(*a, *b);
    CHECK(MolToSmiles(*mol) == "C[*:1].N[*:2]");
  }
  {
    auto a = "[*:2]OC[*:1]"_smiles;
    auto b = "N[*:1]"_smiles;
    auto mol = molzip(*a, *b);
    CHECK(MolToSmiles(*mol) == "NCO[*:2]");
  }
  {
    auto a = "[*:1]OC[*:1]"_smiles;
    auto b = "N[*:1]"_smiles;
    bool caught = false;
    try {
      auto mol = molzip(*a, *b);
    } catch (Invar::Invariant &e) {
      CHECK(e.toUserString().find(
                "molzip: bond info already exists for end atom with label:1") !=
            std::string::npos);
      caught = true;
    }
    CHECK(caught == true);
  }

  {
    // check to see we can make a non sanitizable zipped mol
    MolzipParams p;
    p.enforceValenceRules = false;
    auto a = "CC(=[*:1])N"_smiles;
    auto b = "[*:1]-N=C"_smiles;
    auto mol = molzip(*a, *b, p);
  }

  {
    // check atom property zipping
    auto a = "C=C*.O/C=N/*"_smiles;
    a->getAtomWithIdx(2)->setProp<unsigned int>("fuse", 1);
    a->getAtomWithIdx(6)->setProp<unsigned int>("fuse", 1);
    MolzipParams p;
    p.atomProperty = "fuse";
    p.label = MolzipLabel::AtomProperty;
    auto mol = molzip(*a, p);
    CHECK(MolToSmiles(*mol) == "C=C/N=C/O");
  }

  SECTION("MolZip saves bonddir") {
    {  // a-* b<*
      auto a = "CC[*:1]"_smiles;
      auto b = "N[*:1]"_smiles;
      b->getBondWithIdx(0)->setBondDir(Bond::BondDir::BEGINWEDGE);
      auto mol = molzip(*a, *b);
      CHECK(MolToSmiles(*mol) == "CCN");
      CHECK(mol->getBondWithIdx(1)->getBondDir() == Bond::BondDir::BEGINWEDGE);
      CHECK(mol->getBondWithIdx(1)->getBeginAtom()->getIdx() == 2);
    }
    {  // a>* b-*
      auto a = "[*:1]CC"_smiles;
      auto b = "N[*:1]"_smiles;
      a->getBondWithIdx(0)->setBondDir(Bond::BondDir::BEGINWEDGE);
      auto mol = molzip(*a, *b);
      CHECK(MolToSmiles(*mol) == "CCN");
      CHECK(mol->getBondWithIdx(1)->getBondDir() == Bond::BondDir::BEGINWEDGE);
      CHECK(mol->getBondWithIdx(1)->getBeginAtom()->getIdx() == 2);
    }
    {  // a<* b-*
      auto a = "CC[*:1]"_smiles;
      auto b = "N[*:1]"_smiles;
      a->getBondWithIdx(1)->setBondDir(Bond::BondDir::BEGINWEDGE);
      auto mol = molzip(*a, *b);
      CHECK(MolToSmiles(*mol) == "CCN");
      CHECK(mol->getBondWithIdx(1)->getBondDir() == Bond::BondDir::BEGINWEDGE);
      CHECK(mol->getBondWithIdx(1)->getBeginAtom()->getIdx() == 1);
    }
    {  // a>* b-*
      auto a = "[*:1]CC"_smiles;
      auto b = "N[*:1]"_smiles;
      a->getBondWithIdx(0)->setBondDir(Bond::BondDir::BEGINWEDGE);
      auto mol = molzip(*a, *b);
      CHECK(MolToSmiles(*mol) == "CCN");
      CHECK(mol->getBondWithIdx(1)->getBondDir() == Bond::BondDir::BEGINWEDGE);
      CHECK(mol->getBondWithIdx(1)->getBeginAtom()->getIdx() == 2);
    }
    {  // a-* b<*
      auto a = "CC[*:1]"_smiles;
      auto b = "[*:1]N"_smiles;
      b->getBondWithIdx(0)->setBondDir(Bond::BondDir::BEGINWEDGE);
      auto mol = molzip(*a, *b);
      CHECK(MolToSmiles(*mol) == "CCN");
      CHECK(mol->getBondWithIdx(1)->getBondDir() == Bond::BondDir::BEGINWEDGE);
      CHECK(mol->getBondWithIdx(1)->getBeginAtom()->getIdx() == 1);
    }
  }
}

TEST_CASE(
    "Github4825: ReplaceCore should set stereo on ring bonds when it breaks "
    "rings") {
  SECTION("basics") {
    auto m = "C1C=CCC2=C1C=CC=N2"_smiles;
    REQUIRE(m);
    auto core = "c1ncccc1"_smiles;
    REQUIRE(core);
    std::unique_ptr<ROMol> res{replaceCore(*m, *core)};
    REQUIRE(res);
    auto mb = MolToV3KMolBlock(*res);
    CHECK(mb.find("CFG=2") == std::string::npos);
  }
  SECTION("adjacent") {
    auto m = "C1CC2=C(C=C1)C=CC=N2"_smiles;
    REQUIRE(m);
    auto core = "c1ncccc1"_smiles;
    REQUIRE(core);
    std::unique_ptr<ROMol> res{replaceCore(*m, *core)};
    REQUIRE(res);
    auto mb = MolToV3KMolBlock(*res);
    CHECK(mb.find("CFG=2") == std::string::npos);
  }
  SECTION("don't do larger rings") {
    auto m = "C1C=CCCCC2=C1C=CC=N2"_smiles;
    REQUIRE(m);
    auto core = "c1ncccc1"_smiles;
    REQUIRE(core);
    std::unique_ptr<ROMol> res{replaceCore(*m, *core)};
    REQUIRE(res);
    auto mb = MolToV3KMolBlock(*res);
    CHECK(mb.find("CFG=2") != std::string::npos);
  }
}

TEST_CASE(
    "Github5334: ReplaceCore should set stereo on ring bonds when it breaks "
    "rings") {
  SECTION("segfault") {
    auto a = "C([*:1])[*:2].[C@@H](Cl)([1*:1])[2*:2]"_smiles;
    bool caught = false;
    try {
      auto mol = molzip(*a);
      CHECK(false);
    } catch (Invar::Invariant &e) {
      CHECK(
          e.toUserString().find(
              "molzip: zipped Bond already exists, perhaps labels are duplicated") !=
          std::string::npos);
      caught = true;
    }
    CHECK(caught);
  }
}

TEST_CASE(
    "ReplaceCore handles chiral center with multiple bonds from core to chiral center") {
  auto structure = "C1CSCN[C@@]12(NCCCO2)"_smiles;
  auto core = "NCSCC"_smarts;
  std::unique_ptr<ROMol> res{replaceCore(*structure, *core, true, true)};
  REQUIRE(res);
  auto resultSmiles = MolToSmiles(*res);
  REQUIRE(resultSmiles == "*[C@]1([4*])NCCCO1");
}

TEST_CASE("Molzip with 2D coordinates", "[molzip]") {
  std::vector<std::string> frags = {"c1nc([1*:1])c([2*:2])c([3*:3])n1",
                                    "CC(C)(C#N)c1ccc([1*:1])cc1",
                                    "Nc1ccc(C#C[2*:2])cn1", "OCC[3*:3]"};
  std::vector<ROMOL_SPTR> mols = {
      ROMOL_SPTR(SmilesToMol(frags[0])), ROMOL_SPTR(SmilesToMol(frags[1])),
      ROMOL_SPTR(SmilesToMol(frags[2])), ROMOL_SPTR(SmilesToMol(frags[3]))};
  MolzipParams params;
  params.generateCoordinates = true;
  const auto zippedMol = molzip(mols, params);
  for (auto &mol : mols) {
    for (auto &atom : mol->atoms()) {
      atom->setIsotope(0);
      atom->setAtomMapNum(0);
    }
  }
  const auto &zippedConformer = zippedMol->getConformer();
  for (size_t i = 0; i < frags.size(); i++) {
    const auto sma =
        std::regex_replace(frags[i], std::regex(R"(\[\d\*\:\d\])"), "*");
    const auto query = SmartsToMol(sma);
    const auto &mol = mols[i];
    const auto &molConformer = mol->getConformer();
    const auto molMatches = SubstructMatch(*mol, *query);
    const auto zippedMatches = SubstructMatch(*zippedMol, *query);
    REQUIRE(molMatches.size() == 1);
    REQUIRE(zippedMatches.size() == 1);
    const auto molMatch = molMatches[0];
    const auto zippedMatch = zippedMatches[0];
    for (size_t j = 0; j < molMatch.size(); j++) {
      const auto &position1 = molConformer.getAtomPos(molMatch[i].second);
      const auto &position2 = zippedConformer.getAtomPos(zippedMatch[i].second);
      CHECK(position1.x == position2.x);
      CHECK(position1.y == position2.y);
      CHECK(position1.z == position2.z);
    }
    delete query;
  }
}

TEST_CASE("Molzip with split rings from rgroup", "[molzip]") {
  std::vector<std::string> frags = {"[*:1]CC[*:2]", "[*:1]NN[*:2]",
                                    "[*:1]NN[*:2]"};
  std::vector<ROMOL_SPTR> mols;
  for (auto smi : frags) {
    mols.push_back(ROMOL_SPTR(SmilesToMol(smi)));
  }
  const auto zippedMol = molzip(mols);
  CHECK(MolToSmiles(*zippedMol) == "C1CNN1");
}

TEST_CASE("Github #6034: FragmentOnBonds may create unexpected radicals") {
  auto m = "C[C@H](Cl)c1ccccc1"_smiles;

  REQUIRE(m);
  REQUIRE(m->getNumAtoms() == 9);

  REQUIRE(m->getAtomWithIdx(1)->getNoImplicit() == true);

  bool add_dummies = false;
  std::vector<unsigned int> bonds = {0};
  std::vector<std::pair<unsigned, unsigned>> *dummyLabels = nullptr;
  const std::vector<Bond::BondType> *bondTypes = nullptr;
  std::vector<unsigned> nCutsPerAtom(m->getNumAtoms(), 0);
  std::unique_ptr<ROMol> pieces(MolFragmenter::fragmentOnBonds(
      *m, bonds, add_dummies, dummyLabels, bondTypes, &nCutsPerAtom));
  REQUIRE(pieces);
  REQUIRE(nCutsPerAtom == std::vector<unsigned>{1, 1, 0, 0, 0, 0, 0, 0, 0});

  for (auto at : pieces->atoms()) {
    INFO("atom " + std::to_string(at->getIdx()));
    if (at->getAtomicNum() == 6) {
      CHECK(at->getNoImplicit() == (at->getIdx() == 1));
      CHECK(at->getTotalValence() == 4);
    }
  }
}

TEST_CASE(
    "GitHub #7327: SaltRemover may clear computed properties even if no atoms are removed",
    "[bug]") {
  auto m = "C=CC=O"_smiles;
  REQUIRE(m);
  REQUIRE(m->getRingInfo()->isSymmSssr());

  auto q = "[O,N]"_smarts;
  REQUIRE(q);

  bool onlyFrags = true;
  std::unique_ptr<ROMol> m2{deleteSubstructs(*m, *q, onlyFrags)};
  REQUIRE(m2);
  CHECK(m2->getNumAtoms() == m->getNumAtoms());  // No atoms removed
  CHECK(m2->getRingInfo()->isSymmSssr());
}

TEST_CASE(
    "Github #8288: molzip add linker bond functionality (fixes memory issue)",
    "[feature,bug]") {
  auto m = "[*:1][*:2].C[*:1].S[*:2]"_smiles;
  auto res = molzip(*m);
  CHECK(MolToSmiles(*res) == "CS");
}

TEST_CASE(
    "Github #8389: FragmentOnBonds does not preserve bond stereo when useLegacy=False") {
  auto useLegacy = GENERATE(true, false);
  CAPTURE(useLegacy);
  UseLegacyStereoPerceptionFixture fx(useLegacy);
  SECTION("as reported ") {
    auto m = "CC/C=C/N"_smiles;
    REQUIRE(m);
    std::unique_ptr<ROMol> res(
        MolFragmenter::fragmentOnBonds(*m, std::vector<unsigned int>{1}));
    REQUIRE(res);
    CHECK(res->getBondWithIdx(1)->getBondType() == Bond::BondType::DOUBLE);
    if (useLegacy) {
      CHECK(res->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOE);
    } else {
      CHECK(res->getBondWithIdx(1)->getStereo() ==
            Bond::BondStereo::STEREOTRANS);
    }

    CHECK(MolToSmiles(*res) == "[1*]/C=C/N.[2*]CC");
  }
}

TEST_CASE(
    "Github #9009: fragmentOnBonds should transfer queries on deleted bonds") {
  // Basic test
  auto m1 = "c1ccccc1-C1CCCCC1"_smarts;
  REQUIRE(m1);
  auto bond1 = m1->getBondWithIdx(1);
  REQUIRE(bond1);
  REQUIRE(bond1->hasQuery());
  std::unique_ptr<ROMol> splitM1(
      MolFragmenter::fragmentOnBonds(*m1, std::vector<unsigned int>{1, 5}));
  REQUIRE(splitM1);
  auto nb1 = splitM1->getBondWithIdx(11);
  REQUIRE(nb1);
  CHECK(nb1->hasQuery());
  CHECK(nb1->getQuery()->getDescription() == "SingleOrAromaticBond");

  auto nb2 = splitM1->getBondWithIdx(12);
  REQUIRE(nb2);
  CHECK(nb2->hasQuery());
  CHECK(nb2->getQuery()->getDescription() == "SingleOrAromaticBond");

  auto nb3 = splitM1->getBondWithIdx(13);
  REQUIRE(nb3);
  CHECK(nb3->hasQuery());
  CHECK(nb3->getQuery()->getDescription() == "BondOrder");
  auto nb4 = splitM1->getBondWithIdx(14);
  REQUIRE(nb4);
  CHECK(nb4->hasQuery());
  CHECK(nb4->getQuery()->getDescription() == "BondOrder");

  // Check it does nothing if there isn't a query
  auto m2 = "c1ccccc1"_smiles;
  REQUIRE(m2);
  auto bond2 = m2->getBondWithIdx(1);
  REQUIRE(bond2);
  REQUIRE(!bond2->hasQuery());
  std::unique_ptr<ROMol> splitM2(
      MolFragmenter::fragmentOnBonds(*m2, std::vector<unsigned int>{1}));
  REQUIRE(splitM2);
  auto nb5 = splitM2->getBondWithIdx(5);
  REQUIRE(nb5);
  CHECK(!nb5->hasQuery());
  auto nb6 = splitM2->getBondWithIdx(6);
  REQUIRE(nb6);
  CHECK(!nb6->hasQuery());

  // Check it does nothing if bond types are specified
  auto m3 = "c1ccccc1"_smiles;
  REQUIRE(m3);
  auto bond3 = m3->getBondWithIdx(1);
  REQUIRE(bond3);
  REQUIRE(!bond3->hasQuery());
  auto dummyLabels =
      std::vector<std::pair<unsigned int, unsigned int>>{{25, 26}};
  auto newTypes = std::vector<Bond::BondType>{Bond::DOUBLE};
  std::unique_ptr<ROMol> splitM3(MolFragmenter::fragmentOnBonds(
      *m3, std::vector<unsigned int>{1}, true, &dummyLabels, &newTypes));
  REQUIRE(splitM3);
  auto nb7 = splitM3->getBondWithIdx(5);
  REQUIRE(nb7);
  CHECK(!nb7->hasQuery());
  auto nb8 = splitM3->getBondWithIdx(6);
  REQUIRE(nb8);
  CHECK(!nb8->hasQuery());
}

TEST_CASE("fragmentOnBonds should not segfault on duplicate bondIndices") {
  // Basic test
  auto mol = "NCCO"_smiles;
  REQUIRE(mol);
  REQUIRE(mol->getBondWithIdx(0));
  REQUIRE(mol->getBondWithIdx(2));
  std::unique_ptr<ROMol> splitMol;
  CHECK_THROWS_AS(splitMol.reset(MolFragmenter::fragmentOnBonds(
                      *mol, std::vector<unsigned int>{0, 2, 2})),
                  ValueErrorException);
  REQUIRE(!splitMol);
}

TEST_CASE("align fragments") {
  SECTION("basics") {
    auto m = "CC[1*].O[1*] |(1,0,0;2,0,0;3.5,0,0;0,1,0;0,2.1,0)|"_smiles;
    REQUIRE(m);
    REQUIRE(m->getNumConformers() == 1);
    MolzipParams params;
    params.alignCoordinates = true;
    params.label = MolzipLabel::Isotope;
    auto res = molzip(*m, params);
    REQUIRE(res);
    REQUIRE(res->getNumConformers() == 1);
    CHECK(MolToXYZBlock(*res) == R"XYZ(3

C      1.000000    0.000000    0.000000
C      2.000000    0.000000    0.000000
O      3.100000    0.000000    0.000000
)XYZ");
  }
  SECTION("2D multiple attachment points") {
    auto m =
        "[2*]CCC[1*].[2*]O[1*] |(2,1.2;1,1;1,0;2,0;3,0;\
0,0;0,1;0,2.1)|"_smiles;
    REQUIRE(m);
    REQUIRE(m->getNumConformers() == 1);
    MolzipParams params;
    params.alignCoordinates = true;
    params.label = MolzipLabel::Isotope;
    auto res = molzip(*m, params);
    REQUIRE(res);
    REQUIRE(res->getNumConformers() == 1);
    CHECK(MolToXYZBlock(*res) == R"XYZ(4

C      1.000000    1.000000    0.000000
C      1.000000    0.000000    0.000000
C      2.000000    0.000000    0.000000
O      3.100000    0.000000    0.000000
)XYZ");
  }
  SECTION("3D, more real") {
    auto m1 =
        "[6*]C1CCCCC1 |(-1.08968,0.169713,-0.700677;0.254816,-0.328544,-0.482737;0.56879,-1.01027,0.899618;1.09147,0.21675,1.40903;2.37628,0.102667,0.427592;2.29478,-0.272363,-1.13357;1.24328,0.638789,-1.31851)|"_smiles;
    REQUIRE(m1);
    auto conf1 = m1->getConformer();
    auto m2 =
        "[6*]C1CCC1 |(3.18447,1.72863,2.51726;1.96576,1.94984,2.29932;1.43402,2.74418,3.09468;0.109542,1.58059,3.19096;0.655465,0.761837,2.3866)|"_smiles;
    REQUIRE(m2);
    auto conf2 = m2->getConformer();
    auto origv = conf2.getAtomPos(0) - conf2.getAtomPos(1);

    MolzipParams params;
    params.alignCoordinates = true;
    params.label = MolzipLabel::Isotope;
    auto res = molzip(*m1, *m2, params);
    REQUIRE(res);
    REQUIRE(res->getNumConformers() == 1);
    auto resconf = res->getConformer();
    // check that the first fragment's coords are unchanged
    CHECK_THAT(resconf.getAtomPos(0).x,
               Catch::Matchers::WithinAbs(conf1.getAtomPos(1).x, 1e-4));
    CHECK_THAT(resconf.getAtomPos(0).y,
               Catch::Matchers::WithinAbs(conf1.getAtomPos(1).y, 1e-4));
    CHECK_THAT(resconf.getAtomPos(0).z,
               Catch::Matchers::WithinAbs(conf1.getAtomPos(1).z, 1e-4));
    // make sure the bond length is right:
    auto v = resconf.getAtomPos(0) - resconf.getAtomPos(6);
    CHECK_THAT(v.length(), Catch::Matchers::WithinAbs(origv.length(), 1e-2));
  }
  SECTION("on top") {
    auto m3 =
        "[1*]c1ccccc1.[1*]c1ncccc1 |(2.57143,0,;1.07143,0,;0.321429,-1.29904,;-1.17857,-1.29904,;-1.92857,0,;-1.17857,1.29904,;0.321429,1.29904,;12.5714,0,;11.0714,0,;10.3214,1.29904,;8.82143,1.29904,;8.07143,0,;8.82143,-1.29904,;10.3214,-1.29904,)|"_smiles;
    MolzipParams params;
    params.alignCoordinates = true;
    params.label = MolzipLabel::Isotope;
    auto res = molzip(*m3, params);
    INFO(MolToCXSmiles(*res));
    auto resconf = res->getConformer();
    // check that we don't have overlapping atoms:
    for (unsigned int i = 0; i < res->getNumAtoms(); ++i) {
      auto api = resconf.getAtomPos(i);
      for (unsigned int j = i + 1; j < res->getNumAtoms(); ++j) {
        auto apj = resconf.getAtomPos(j);
        CHECK((api - apj).length() > 1.0);
      }
    }
  }
  SECTION("first fragment should not move") {
    auto m3 =
        "[1*]c1ccccc1.[1*]c1ncccc1 |(2.57143,0,;1.07143,0,;0.321429,-1.29904,;-1.17857,-1.29904,;-1.92857,0,;-1.17857,1.29904,;0.321429,1.29904,;12.5714,0,;11.0714,0,;10.3214,1.29904,;8.82143,1.29904,;8.07143,0,;8.82143,-1.29904,;10.3214,-1.29904,)|"_smiles;
    MolzipParams params;
    params.alignCoordinates = true;
    params.label = MolzipLabel::Isotope;
    auto res = molzip(*m3, params);
    INFO(MolToCXSmiles(*res));
    auto resconf = res->getConformer();
    auto molconf = m3->getConformer();
    for (unsigned int i = 0; i < 6; ++i) {
      auto rap = resconf.getAtomPos(i);
      auto map = molconf.getAtomPos(i + 1);
      CHECK_THAT(rap.x, Catch::Matchers::WithinAbs(map.x, 1e-4));
      CHECK_THAT(rap.y, Catch::Matchers::WithinAbs(map.y, 1e-4));
      CHECK_THAT(rap.z, Catch::Matchers::WithinAbs(map.z, 1e-4));
    }
  }
}

TEST_CASE("github #9468") {
  SECTION("as reported") {
    auto m = "CCCCC(=O)N1CCCC1"_smiles;
    REQUIRE(m);
    std::vector<unsigned int> bonds{{0, 1, 2, 3}};
    bool addDummies = false;
    std::unique_ptr<ROMol> res{
        MolFragmenter::fragmentOnBonds(*m, bonds, addDummies)};
    REQUIRE(res);
    CHECK(res->getNumAtoms() == m->getNumAtoms());
    CHECK(res->getAtomWithIdx(0)->getTotalNumHs() == 4);
    CHECK(res->getAtomWithIdx(1)->getTotalNumHs() == 4);
    CHECK(res->getAtomWithIdx(2)->getTotalNumHs() == 4);
    CHECK(res->getAtomWithIdx(3)->getTotalNumHs() == 4);
  }
  SECTION("aromaticity") {
    auto m = "Cn1cccc1"_smiles;
    REQUIRE(m);
    std::vector<unsigned int> bonds = {0};
    bool addDummies = false;
    std::unique_ptr<ROMol> res{
        MolFragmenter::fragmentOnBonds(*m, bonds, addDummies)};
    REQUIRE(res);
    CHECK(res->getNumAtoms() == m->getNumAtoms());
    CHECK(res->getAtomWithIdx(0)->getTotalNumHs() == 4);
    CHECK(res->getAtomWithIdx(1)->getTotalNumHs() == 1);
    CHECK(res->getAtomWithIdx(1)->getNumExplicitHs() == 1);
  }
  SECTION("no changes to H counts if dummies are added") {
    auto m = "CCCCC(=O)N1CCCC1"_smiles;
    REQUIRE(m);
    std::vector<unsigned int> bonds{{0, 1, 2, 3}};
    bool addDummies = true;
    std::unique_ptr<ROMol> res{
        MolFragmenter::fragmentOnBonds(*m, bonds, addDummies)};
    REQUIRE(res);
    CHECK(res->getNumAtoms() > m->getNumAtoms());
    for (unsigned int i = 0; i < 4; ++i) {
      CHECK(res->getAtomWithIdx(i)->getTotalNumHs() ==
            m->getAtomWithIdx(i)->getTotalNumHs());
    }
    // make sure dummies have implicit props too
    for (unsigned int i = m->getNumAtoms(); i < res->getNumAtoms(); ++i) {
      CHECK(res->getAtomWithIdx(i)->getTotalNumHs() == 0);
    }
  }
  SECTION("test case from Rachel Walker") {
    auto mol = "c1ccccc1CC"_smiles;
    REQUIRE(mol);
    const auto *coreBond = mol->getBondBetweenAtoms(5, 6);
    const auto *sidechainBond = mol->getBondBetweenAtoms(6, 7);
    REQUIRE(coreBond);
    REQUIRE(sidechainBond);
    bool addDummies = false;
    std::unique_ptr<ROMol> splitMol(MolFragmenter::fragmentOnBonds(
        *mol, {coreBond->getIdx(), sidechainBond->getIdx()}, addDummies));
    REQUIRE(splitMol);
    CHECK(splitMol->getAtomWithIdx(5)->getTotalNumHs() == 1);
    CHECK(splitMol->getAtomWithIdx(6)->getTotalNumHs() == 4);
    CHECK(splitMol->getAtomWithIdx(7)->getTotalNumHs() == 4);
  }
  SECTION("breaking double bonds") {
    auto mol = "C=CC"_smiles;
    REQUIRE(mol);
    std::vector<unsigned int> bonds{0};
    {
      bool addDummies = false;
      std::unique_ptr<ROMol> splitMol(
          MolFragmenter::fragmentOnBonds(*mol, bonds, addDummies));
      REQUIRE(splitMol);
      CHECK(splitMol->getAtomWithIdx(0)->getTotalNumHs() == 4);
      CHECK(splitMol->getAtomWithIdx(1)->getTotalNumHs() == 3);
    }
    {
      bool addDummies = true;
      std::unique_ptr<ROMol> splitMol(
          MolFragmenter::fragmentOnBonds(*mol, bonds, addDummies));
      REQUIRE(splitMol);
      CHECK(splitMol->getAtomWithIdx(0)->getTotalNumHs() == 2);
      CHECK(splitMol->getAtomWithIdx(1)->getTotalNumHs() == 1);
    }
  }
  SECTION("breaking multiple bonds with NoImplicit") {
    // The brackets are totally redundant here, but they are required
    // to trigger the issue
    auto mol = "C[C](C)(C)C"_smiles;
    REQUIRE(mol);
    std::vector<unsigned int> bonds{0, 1};
    {
      bool addDummies = false;
      std::unique_ptr<ROMol> splitMol(
          MolFragmenter::fragmentOnBonds(*mol, bonds, addDummies));
      REQUIRE(splitMol);

      auto atom = splitMol->getAtomWithIdx(1);
      CHECK(atom->getTotalNumHs() == 2);
      CHECK(atom->getNoImplicit() == true);
    }
  }
}

TEST_CASE("MurckoDecompose", "[MurckoDecompose]") {
  const std::vector<std::pair<std::string, std::string>> data = {
      {"C1CC1CCC1CC(C)C1", "C1CC1CCC1CCC1"},
      {"CC1CCC1", "C1CCC1"},
      {"NCNCC2CC2C1CC1O", "C1CC1C1CC1"},
      {"OC2C(C)C21C(N)C1C", "C2CC12CC1"},  // Spiro
      {"C1CC1C(=O)OC", "C1CC1"},           // Carbonyl outside scaffold
      {"C1CC1C=C", "C1CC1"},               // Double bond outside scaffold
      {"C1CC1C=CC1CC1C=CNNCO", "C1CC1C=CC1CC1"},  // Double bond in scaffold
      {"CC1CC1C(N)C1C(N)C1", "C1CC1CC1CC1"},
      {"C1CC1S(=O)C1CC1C=CNNCO", "C1CC1S(=O)C1CC1"},  // S=O group in scaffold
      {"O=SCNC1CC1S(=O)C1CC1C=CNNCO",
       "C1CC1S(=O)C1CC1"},  // S=O group outside scaffold
      {"C1CC1S(=O)(=O)C1CC1C=CNNCO",
       "C1CC1S(=O)(=O)C1CC1"},  // SO2 group in scaffold
      {"O=S(CNCNC)(=O)CNC1CC1S(=O)(=O)C1CC1C=CNNCO",
       "C1CC1S(=O)(=O)C1CC1"},                    // SO2 group outside scaffold
      {"C1CC1C=NO", "C1CC1"},                     // Hydroxamide
      {"C1CC1C(C(C)C)=NC1CC1", "C1CC1C=NC1CC1"},  // Hydroxamide
      {"C1CC1C#N", "C1CC1"},                      // Cyano group
      {"C1CC1C#CNC", "C1CC1"},                    // Acetylene group
      {"O=C1N(C)C(=O)N1C#CNC", "O=C1NC(=O)N1"},   // Acetylene group
      {"[O-][N+](=O)c1cc(ccc1Cl)NS(=O)(=O)Cc2ccccc2",
       "c1ccccc1NS(=O)(=O)Cc2ccccc2"},
      {"Cn1cccc1", "c1ccc[nH]1"},
      {"C1CC1[CH](C)C1CC1", "C1CC1CC1CC1"},
      {"CC(C)c1c(Cl)cc(F)c(Nc2ccc(C)cc2CC(=O)O)c1F", "c1ccc(Nc2ccccc2)cc1"},
      {"C1CC1C[C@](Cl)(F)C1CCC1", "C1CC1CCC1CCC1"},
      {"C1CC1.C1CC1", "C1CC1.C1CC1"},  // this was a crash at one point
      {"CCC", ""},
  };
  for (const auto &[smi, tgt] : data) {
    INFO(smi);
    std::unique_ptr<ROMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    std::unique_ptr<ROMol> nMol(MurckoDecompose(*mol));
    REQUIRE(nMol);
    if (tgt.empty()) {
      CHECK(nMol->getNumAtoms() == 0);
    } else {
      REQUIRE(nMol->getNumAtoms());
      MolOps::sanitizeMol(static_cast<RWMol &>(*nMol));
      std::unique_ptr<ROMol> tgtMol(SmilesToMol(tgt));
      REQUIRE(tgtMol);
      CHECK(MolToSmiles(*nMol) == MolToSmiles(*tgtMol));
    }
  }
}

TEST_CASE("Issue 3537675: MurckoDecompose with substituted non-carbon "
          "aromatic atoms",
          "[MurckoDecompose]") {
  for (const auto smi : {"c1cccp1C", "c1cccn1C"}) {
    INFO(smi);
    std::unique_ptr<ROMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    std::unique_ptr<ROMol> nMol(MurckoDecompose(*mol));
    REQUIRE(nMol);
    CHECK(nMol->getAtomWithIdx(4)->getNumExplicitHs());
  }
}

TEST_CASE("sf.net issue 275: MurckoDecompose with chiral atoms",
          "[MurckoDecompose]") {
  auto mol = "CCCCC[C@H]1CC[C@H](C(=O)O)CC1"_smiles;
  REQUIRE(mol);
  std::unique_ptr<ROMol> nMol(MurckoDecompose(*mol));
  REQUIRE(nMol);
  CHECK(MolToSmiles(*nMol, true) == "C1CCCCC1");
  CHECK(MolToSmiles(*nMol, false) == "C1CCCCC1");
}

TEST_CASE("MurckoDecompose keeps atom order and handles edge cases",
          "[MurckoDecompose]") {
  // atom map numbers record the input index of each atom, so the expected
  // non-canonical SMILES checks the atoms, their order, and their H counts
  const std::vector<std::pair<std::string, std::string>> data = {
      {"CCCCC", ""},
      {"[Na+].[Cl-]", ""},
      {"O=C1CCCC1", "[O:1]=[C:2]1[CH2:3][CH2:4][CH2:5][CH2:6]1"},
      {"CC(=O)C1CC1", "[CH2:4]1[CH2:5][CH2:6]1"},
      {"C1CC1CC(=O)CC1CC1",
       "[CH2:1]1[CH2:2][CH:3]1[CH2:4][C:5](=[O:6])[CH2:7][CH:8]1[CH2:9][CH2:"
       "10]1"},
      {"C1CCC2(CC1)CCCC2",
       "[CH2:1]1[CH2:2][CH2:3][C:4]2([CH2:5][CH2:6]1)[CH2:7][CH2:8][CH2:9]["
       "CH2:10]2"},
      {"C1CC11CC1", "[CH2:1]1[CH2:2][C:3]12[CH2:4][CH2:5]2"},
      {"c1cc[nH]c1CC", "[cH:1]1[cH:2][cH:3][nH:4][cH:5]1"},
      {"c1cc[n+](C)cc1", "[cH:1]1[cH:2][cH:3][nH+:4][cH:6][cH:7]1"},
      {"CC[C@H](N)c1ccccc1", "[cH:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1"},
      {"c1ccccc1C.CCCC.C1CC1CC",
       "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1.[CH2:12]1[CH2:13][CH2:14]1"},
  };
  SmilesWriteParams ps;
  ps.canonical = false;
  ps.allHsExplicit = true;
  for (const auto &[smi, expected] : data) {
    for (bool explicitHs : {false, true}) {
      INFO(smi << " explicitHs: " << explicitHs);
      std::unique_ptr<RWMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      for (auto atom : mol->atoms()) {
        atom->setAtomMapNum(atom->getIdx() + 1);
      }
      if (explicitHs) {
        MolOps::addHs(*mol);
      }
      std::unique_ptr<ROMol> res(MurckoDecompose(*mol));
      REQUIRE(res);
      CHECK(MolToSmiles(*res, ps) == expected);
    }
  }
}

TEST_CASE("MurckoDecompose on a large molecule", "[MurckoDecompose]") {
  // this used to scale as O(N^3) in the number of atoms
  const unsigned int nUnits = 200;
  std::string smi = "c1ccccc1";
  for (unsigned int i = 0; i < nUnits; ++i) {
    smi += "CC(CCO)Cc1ccc(cc1)";
  }
  smi += "C";
  std::unique_ptr<ROMol> mol(SmilesToMol(smi));
  REQUIRE(mol);
  std::unique_ptr<ROMol> res(MurckoDecompose(*mol));
  REQUIRE(res);
  // each unit keeps its three-atom linker and its ring, the CCO side chain
  // and the terminal methyl go
  CHECK(res->getNumAtoms() == 6 + 9 * nUnits);
  CHECK(res->getNumBonds() == 6 + 10 * nUnits);
  CHECK(!mol->hasProp(common_properties::DistanceMatrix_Paths));
}
