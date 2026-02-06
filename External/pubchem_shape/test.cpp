#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/RWMol.h>

#include "PubChemShape.hpp"
#include "GraphMol/SmilesParse/SmilesWrite.h"

using namespace RDKit;

constexpr double test_opt_param = 0.5;
constexpr unsigned int test_max_preiters = 3;
constexpr unsigned int test_max_postiters = 16;
constexpr bool test_use_colors = true;

TEST_CASE("basic alignment") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";

  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/test1.sdf");
  auto ref = suppl[0];
  REQUIRE(ref);
  auto probe = suppl[1];
  REQUIRE(probe);
  SECTION("basics") {
    RWMol cp(*probe);
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] =
        AlignMolecule(*ref, cp, matrix, -1, -1, test_use_colors, test_opt_param,
                      test_max_preiters, test_max_postiters);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.773, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.303, 0.005));
    for (auto i = 0u; i < probe->getNumAtoms(); ++i) {
      CHECK(probe->getConformer().getAtomPos(i).x !=
            cp.getConformer().getAtomPos(i).x);
    }
  }
  SECTION("from shape") {
    auto ref_shape = PrepareConformer(*ref);
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] =
        AlignMolecule(ref_shape, *probe, matrix, -1, test_use_colors,
                      test_opt_param, test_max_preiters, test_max_postiters);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.773, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.303, 0.005));
  }
  SECTION("RDKit features") {
    ref->clearProp("PUBCHEM_PHARMACOPHORE_FEATURES");
    probe->clearProp("PUBCHEM_PHARMACOPHORE_FEATURES");
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] =
        AlignMolecule(*ref, *probe, matrix, -1, -1, test_use_colors,
                      test_opt_param, test_max_preiters, test_max_postiters);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.773, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.231, 0.005));
  }
  SECTION("no colors") {
    std::cout << "no oclors" << std::endl;
    std::vector<float> matrix(12, 0.0);
    int refConfId = -1;
    int prbConfId = -1;
    bool useColors = false;
    auto [nbr_st, nbr_ct] =
        AlignMolecule(*ref, *probe, matrix, refConfId, prbConfId, useColors,
                      test_opt_param, test_max_preiters, test_max_postiters);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.773, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.000, 0.005));
  }
  SECTION("we are re-entrant") {
    RWMol cp(*probe);
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] =
        AlignMolecule(*ref, cp, matrix, -1, -1, test_use_colors, test_opt_param,
                      test_max_preiters, test_max_postiters);
    RWMol cp2(cp);
    auto [nbr_st2, nbr_ct2] =
        AlignMolecule(*ref, cp2, matrix, -1, -1, test_use_colors,
                      test_opt_param, test_max_preiters, test_max_postiters);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(nbr_st2, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(nbr_ct2, 0.005));
    for (auto i = 0u; i < probe->getNumAtoms(); ++i) {
      CHECK_THAT(
          cp.getConformer().getAtomPos(i).x,
          Catch::Matchers::WithinAbs(cp2.getConformer().getAtomPos(i).x, 0.05));
      CHECK_THAT(
          cp.getConformer().getAtomPos(i).y,
          Catch::Matchers::WithinAbs(cp2.getConformer().getAtomPos(i).y, 0.05));
      CHECK_THAT(
          cp.getConformer().getAtomPos(i).z,
          Catch::Matchers::WithinAbs(cp2.getConformer().getAtomPos(i).z, 0.05));
    }
  }
}

TEST_CASE("bulk") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";
  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/bulk.pubchem.sdf");
  auto ref = suppl[0];
  REQUIRE(ref);
  for (auto i = 1u; i < suppl.length(); ++i) {
    auto probe = suppl[1];
    REQUIRE(probe);
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] =
        AlignMolecule(*ref, *probe, matrix, -1, -1, test_use_colors,
                      test_opt_param, test_max_preiters, test_max_postiters);
    CHECK_THAT(nbr_st,
               Catch::Matchers::WithinAbs(
                   probe->getProp<float>("shape_align_shape_tanimoto"), 0.005));
    CHECK_THAT(nbr_ct,
               Catch::Matchers::WithinAbs(
                   probe->getProp<float>("shape_align_color_tanimoto"), 0.005));
  }
}

TEST_CASE("handling molecules with Hs") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";

  v2::FileParsers::MolFileParserParams params;
  params.removeHs = false;
  auto suppl =
      v2::FileParsers::SDMolSupplier(dirName + "/align_with_hs.sdf", params);
  auto ref = suppl[0];
  REQUIRE(ref);
  auto probe = suppl[1];
  REQUIRE(probe);
  ShapeInputOptions mol1Opts, mol2Opts;
  SECTION("basics") {
    RWMol cp(*probe);
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] =
        AlignMolecule(*ref, cp, matrix, -1, -1, test_use_colors, test_opt_param,
                      test_max_preiters, test_max_postiters);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.840, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.753, 0.005));
    for (auto i = 0u; i < cp.getNumAtoms(); ++i) {
      // the failure mode here was that Hs had HUGE coordinates
      auto pos = cp.getConformer().getAtomPos(i);
      CHECK((pos.x > -10 && pos.x < 10));
    }
  }
}

TEST_CASE("re-entrant") {
  SECTION("basics") {
    auto ref = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.333279 6.212121 0.00000 0
M  V30 2 C -5.393939 5.151461 0.100000 0
M  V30 3 C -6.454600 6.212121 0.000000 0
M  V30 4 C -5.393939 7.272781 0.100000 0
M  V30 5 C -2.833279 6.212121 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 1
M  V30 5 1 5 1
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
)CTAB"_ctab;

    auto probe = R"CTAB(
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.333279 8.212121 0.000000 0
M  V30 2 C -5.393939 7.151461 0.100000 0
M  V30 3 C -6.454600 8.212121 0.000000 0
M  V30 4 C -5.393939 9.272781 0.100000 0
M  V30 5 C -2.833279 8.212121 0.000000 0
M  V30 6 C -1.333279 8.212121 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 1
M  V30 5 1 5 1
M  V30 6 1 5 6
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
)CTAB"_ctab;

    MolOps::removeAllHs(*ref);
    MolOps::removeAllHs(*probe);

    RWMol cp(*probe);
    std::vector<float> matrix(12, 0.0);
    bool useColors = true;
    int refConfId = -1;
    int prbConfId = -1;
    double opt_param = 1.0;
    unsigned int preiters = 100u;
    unsigned int postiters = 100u;
    auto [nbr_st, nbr_ct] =
        AlignMolecule(*ref, cp, matrix, refConfId, prbConfId, useColors,
                      opt_param, preiters, postiters);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.923, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.882, 0.005));

    auto rmsd = MolAlign::CalcRMS(*ref, cp);
    CHECK_THAT(rmsd, Catch::Matchers::WithinAbs(0.255, 0.005));

    RWMol cp2(cp);
    auto [nbr_st2, nbr_ct2] =
        AlignMolecule(*ref, cp2, matrix, refConfId, prbConfId, useColors,
                      opt_param, preiters, postiters);
    CHECK_THAT(nbr_st2, Catch::Matchers::WithinAbs(nbr_st2, 0.005));
    CHECK_THAT(nbr_ct2, Catch::Matchers::WithinAbs(nbr_ct2, 0.005));

    auto rmsd2 = MolAlign::CalcRMS(cp, cp2);
    CHECK_THAT(rmsd2, Catch::Matchers::WithinAbs(0.0, 0.005));
  }
  SECTION("real example") {
    std::string dirName = getenv("RDBASE");
    dirName += "/External/pubchem_shape/test_data";

    v2::FileParsers::MolFileParserParams params;
    params.removeHs = false;
    auto suppl =
        v2::FileParsers::SDMolSupplier(dirName + "/P17612.sdf", params);
    auto ref = suppl[0];
    REQUIRE(ref);
    auto probe = suppl[1];
    REQUIRE(probe);

    MolOps::removeAllHs(*ref);
    MolOps::removeAllHs(*probe);

    RWMol cp(*probe);
    std::vector<float> matrix(12, 0.0);
    bool useColors = true;
    int refConfId = -1;
    int prbConfId = -1;
    double opt_param = 1.0;
    unsigned int preiters = 100u;
    unsigned int postiters = 100u;
    auto [nbr_st, nbr_ct] =
        AlignMolecule(*ref, cp, matrix, refConfId, prbConfId, useColors,
                      opt_param, preiters, postiters);

    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.483, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.046, 0.005));
    RWMol cp2(cp);
    auto [nbr_st2, nbr_ct2] =
        AlignMolecule(*ref, cp2, matrix, refConfId, prbConfId, useColors,
                      opt_param, preiters, postiters);

    CHECK_THAT(nbr_st2, Catch::Matchers::WithinAbs(nbr_st, 0.005));
    CHECK_THAT(nbr_ct2, Catch::Matchers::WithinAbs(nbr_ct, 0.005));

    auto rmsd = MolAlign::CalcRMS(cp, cp2);
    CHECK_THAT(rmsd, Catch::Matchers::WithinAbs(0.007, 0.005));
  }
}

TEST_CASE("Github #8096") {
  SECTION("as reported") {
    auto m1 =
        "[H]c1c([H])c([H])c([2H])c([H])c1[H] |(1.55967,1.91617,0.0546381;0.885536,1.07172,0.030849;1.38172,-0.23747,0.0274262;2.44539,-0.439501,0.0483424;0.470206,-1.27516,-0.00361916;0.856925,-2.30002,-0.00633525;-0.896665,-1.07227,-0.0310991;-1.60071,-1.87642,-0.0551085;-1.36315,0.22877,-0.0271173;-2.43593,0.379132,-0.0487835;-0.479018,1.29083,0.00359778;-0.823965,2.31421,0.00720933)|"_smiles;
    REQUIRE(m1);
    auto m2 =
        "[H]c1c([H])c([H])c([H])c([H])c1[H] |(-2.06264,-0.844763,-0.0261403;-1.04035,-0.481453,-0.0114878;-0.00743655,-1.41861,-0.0137121;-0.215455,-2.47997,-0.0295909;1.29853,-0.949412,0.00507497;2.12524,-1.65277,0.00390664;1.58501,0.395878,0.0254188;2.61997,0.704365,0.0394811;0.550242,1.31385,0.0273741;0.783172,2.37039,0.0434262;-0.763786,0.88847,0.00908113;-1.60557,1.58532,0.0100194)|"_smiles;
    REQUIRE(m2);
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] = AlignMolecule(*m1, *m2, matrix);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(1.0, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(1.0, 0.005));
  }
}

#ifdef RDK_USE_BOOST_SERIALIZATION
TEST_CASE("Serialization") {
  auto m1 =
      "[H]c1c([H])c([H])c([H])c([H])c1[H] |(-2.06264,-0.844763,-0.0261403;-1.04035,-0.481453,-0.0114878;-0.00743655,-1.41861,-0.0137121;-0.215455,-2.47997,-0.0295909;1.29853,-0.949412,0.00507497;2.12524,-1.65277,0.00390664;1.58501,0.395878,0.0254188;2.61997,0.704365,0.0394811;0.550242,1.31385,0.0273741;0.783172,2.37039,0.0434262;-0.763786,0.88847,0.00908113;-1.60557,1.58532,0.0100194)|"_smiles;
  REQUIRE(m1);
  auto shape = PrepareConformer(*m1);
  auto istr = shape.toString();

  ShapeInput shape2(istr);
  CHECK(shape2.coord == shape.coord);
  CHECK(shape2.alpha_vector == shape.alpha_vector);
  CHECK(shape2.atom_type_vector == shape.atom_type_vector);
  CHECK(shape2.volumeAtomIndexVector == shape.volumeAtomIndexVector);
  CHECK(shape2.colorAtomType2IndexVectorMap ==
        shape.colorAtomType2IndexVectorMap);
  CHECK(shape2.shift == shape.shift);
  CHECK_THAT(shape2.sov, Catch::Matchers::WithinAbs(253.764, 0.005));
  CHECK_THAT(shape2.sof, Catch::Matchers::WithinAbs(5.074, 0.005));
}
#endif

TEST_CASE("d2CutOff set") {
  // Previously, shape1.sov and shape2.sov were slightly different because
  // useCutOff was only being set in AlignMolecule, not PrepareConformer.
  // Thus aligning a different molecule meant that PrepareConformer gave
  // different results before and after the alignment.
  auto m1 =
      "c1ccc(-c2ccccc2)cc1 |(-3.26053,-0.0841607,-0.741909;-2.93383,0.123873,0.593407;-1.60713,0.377277,0.917966;-0.644758,0.654885,-0.0378428;0.743308,0.219134,0.168663;1.82376,1.0395,-0.0112769;3.01462,0.695405,0.613858;3.18783,-0.589771,1.09649;2.15761,-1.50458,1.01949;0.988307,-1.1313,0.385783;-1.1048,0.797771,-1.34022;-2.39754,0.435801,-1.69921)|"_smiles;
  REQUIRE(m1);
  auto shape1 = PrepareConformer(*m1, -1);
  std::vector<float> matrix(12, 0.0);
  auto m2 =
      "c1ccc(-c2ccccc2)cc1 |(-3.26053,-0.0841607,-0.741909;-2.93383,0.123873,0.593407;-1.60713,0.377277,0.917966;-0.644758,0.654885,-0.0378428;0.743308,0.219134,0.168663;1.82376,1.0395,-0.0112769;3.01462,0.695405,0.613858;3.18783,-0.589771,1.09649;2.15761,-1.50458,1.01949;0.988307,-1.1313,0.385783;-1.1048,0.797771,-1.34022;-2.39754,0.435801,-1.69921)|"_smiles;
  REQUIRE(m2);
  AlignMolecule(*m2, *m2, matrix);
  auto shape2 = PrepareConformer(*m1, -1);
  CHECK(shape1.sov == shape2.sov);
}

TEST_CASE("Overlay onto shape bug (Github8462)") {
  auto m1 =
      "c1ccc(-c2ccccc2)cc1 |(-3.26053,-0.0841607,-0.741909;-2.93383,0.123873,0.593407;-1.60713,0.377277,0.917966;-0.644758,0.654885,-0.0378428;0.743308,0.219134,0.168663;1.82376,1.0395,-0.0112769;3.01462,0.695405,0.613858;3.18783,-0.589771,1.09649;2.15761,-1.50458,1.01949;0.988307,-1.1313,0.385783;-1.1048,0.797771,-1.34022;-2.39754,0.435801,-1.69921)|"_smiles;
  REQUIRE(m1);
  ROMol m2(*m1);
  for (auto a : m2.atoms()) {
    auto &pos = m2.getConformer().getAtomPos(a->getIdx());
    pos.x += 3.0;
    pos.y += 2.0;
  }
  ROMol m3(m2);

  std::vector<float> matrix(12, 0.0);
  auto [st, ct] = AlignMolecule(*m1, m2, matrix);
  CHECK_THAT(st, Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(ct, Catch::Matchers::WithinAbs(1.0, 0.005));
  for (unsigned int i = 0; i < m1->getNumAtoms(); ++i) {
    auto pos1 = m1->getConformer().getAtomPos(i);
    auto pos2 = m2.getConformer().getAtomPos(i);
    CHECK_THAT((pos1 - pos2).length(), Catch::Matchers::WithinAbs(0.0, 0.005));
  }

  auto s1 = PrepareConformer(*m1, -1);
  auto [st1, ct1] = AlignMolecule(s1, m3, matrix);
  CHECK_THAT(st1, Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(ct1, Catch::Matchers::WithinAbs(1.0, 0.005));
  for (unsigned int i = 0; i < m3.getNumAtoms(); ++i) {
    RDGeom::Point3D pos1(s1.coord[3 * i], s1.coord[3 * i + 1],
                         s1.coord[3 * i + 2]);
    auto pos2 = m3.getConformer().getAtomPos(i);
    CHECK_THAT((pos1 - pos2).length(), Catch::Matchers::WithinAbs(0.0, 0.005));
  }
}

TEST_CASE("Shape subset") {
  auto m1 =
      "c1ccc(-c2ccccc2)cc1 |(-3.26053,-0.0841607,-0.741909;-2.93383,0.123873,0.593407;-1.60713,0.377277,0.917966;-0.644758,0.654885,-0.0378428;0.743308,0.219134,0.168663;1.82376,1.0395,-0.0112769;3.01462,0.695405,0.613858;3.18783,-0.589771,1.09649;2.15761,-1.50458,1.01949;0.988307,-1.1313,0.385783;-1.1048,0.797771,-1.34022;-2.39754,0.435801,-1.69921)|"_smiles;
  REQUIRE(m1);
  ShapeInputOptions shapeOpts;
  shapeOpts.atomSubset = std::vector<unsigned int>{0, 1, 2, 3, 10, 11};
  auto partShape = PrepareConformer(*m1, -1, shapeOpts);
  CHECK(partShape.coord.size() == 21);
  CHECK_THAT(partShape.sov, Catch::Matchers::WithinAbs(253.929, 0.005));
  CHECK_THAT(partShape.sof, Catch::Matchers::WithinAbs(5.074, 0.005));

  shapeOpts.atomSubset.clear();
  auto wholeShape = PrepareConformer(*m1, -1, shapeOpts);
  CHECK(wholeShape.coord.size() == 42);
  CHECK_THAT(wholeShape.sov, Catch::Matchers::WithinAbs(542.04, 0.005));
  CHECK_THAT(wholeShape.sof, Catch::Matchers::WithinAbs(10.148, 0.005));
}

TEST_CASE("Dummy radii") {
  auto m1 =
      "[Xe]c1ccccc1 |(0.392086,-2.22477,0.190651;0.232269,-1.38667,0.118385;-1.06274,-0.918982,0.0342466;-1.26098,0.446053,-0.0811879;-0.244035,1.36265,-0.11691;1.05134,0.875929,-0.031248;1.28797,-0.499563,0.0864097),atomProp:0.dummyLabel.*|"_smiles;
  auto shape1 = PrepareConformer(*m1, -1);
  CHECK(shape1.coord.size() == 24);

  // The dummy radius defaults to 2.16, the same as Xe, so these shapes should
  // come out the same.
  auto m2 =
      "*c1ccccc1 |(0.392086,-2.22477,0.190651;0.232269,-1.38667,0.118385;-1.06274,-0.918982,0.0342466;-1.26098,0.446053,-0.0811879;-0.244035,1.36265,-0.11691;1.05134,0.875929,-0.031248;1.28797,-0.499563,0.0864097),atomProp:0.dummyLabel.*|"_smiles;
  ShapeInputOptions shapeOpts;
  shapeOpts.includeDummies = true;
  auto shape2 = PrepareConformer(*m2, -1, shapeOpts);
  CHECK(shape2.coord.size() == 24);
  CHECK(shape1.sov == shape2.sov);

  shapeOpts.dummyRadius = 2.5;
  auto shape3 = PrepareConformer(*m2, -1, shapeOpts);
  CHECK_THAT(shape3.sov, Catch::Matchers::WithinAbs(427.925, 0.005));

  shapeOpts.includeDummies = false;
  auto shape4 = PrepareConformer(*m2, -1, shapeOpts);
  CHECK(shape4.coord.size() == 21);
  CHECK(shape4.sov < shape2.sov);
  CHECK_THAT(shape4.sov, Catch::Matchers::WithinAbs(254.578, 0.005));
}

TEST_CASE("Non-standard radii") {
  auto m1 =
      "[Xe]c1ccccc1 |(0.392086,-2.22477,0.190651;0.232269,-1.38667,0.118385;-1.06274,-0.918982,0.0342466;-1.26098,0.446053,-0.0811879;-0.244035,1.36265,-0.11691;1.05134,0.875929,-0.031248;1.28797,-0.499563,0.0864097),atomProp:0.dummyLabel.*|"_smiles;
  auto shape1 = PrepareConformer(*m1, -1);
  CHECK(shape1.coord.size() == 24);
  CHECK_THAT(shape1.sov, Catch::Matchers::WithinAbs(376.434, 0.005));

  ShapeInputOptions shapeOpts;
  // Benzene derivative with atom 4 with an N radius.
  shapeOpts.atomRadii =
      std::vector<std::pair<unsigned int, double>>{{0, 2.5}, {4, 1.55}};
  auto shape2 = PrepareConformer(*m1, -1, shapeOpts);
  CHECK_THAT(shape2.sov, Catch::Matchers::WithinAbs(412.666, 0.005));

  // Corresponding pyridine derivative.
  auto m2 =
      "[Xe]c1ccncc1 |(0.392086,-2.22477,0.190651;0.232269,-1.38667,0.118385;-1.06274,-0.918982,0.0342466;-1.26098,0.446053,-0.0811879;-0.244035,1.36265,-0.11691;1.05134,0.875929,-0.031248;1.28797,-0.499563,0.0864097),atomProp:0.dummyLabel.*|"_smiles;
  auto shape3 = PrepareConformer(*m2, -1, shapeOpts);
  CHECK(shape3.sov == shape2.sov);
}

TEST_CASE("Shape-Shape alignment") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";

  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/test1.sdf");
  auto ref = suppl[0];
  REQUIRE(ref);
  auto probe = suppl[1];
  REQUIRE(probe);
  RWMol probeCP(*probe);

  auto refShape = PrepareConformer(*ref, -1);
  auto probeShape = PrepareConformer(*probe, -1);

  std::vector<float> matrix(12, 0.0);
  auto [mol_st, mol_ct] = AlignMolecule(*ref, *probe, matrix, -1, -1);
  auto [shape_st, shape_ct] = AlignShape(refShape, probeShape, matrix);

  // Check that the same results are achieved when overlaying the probe
  // molecule onto the reference molecule and the probe shape onto the reference
  // shape. The shapes are transformed to their centroid and principal axes, so
  // the final coords will not match.  Just check the tanimotos.
  CHECK_THAT(shape_st, Catch::Matchers::WithinAbs(mol_st, 0.001));
  CHECK_THAT(shape_ct, Catch::Matchers::WithinAbs(mol_ct, 0.001));
}

TEST_CASE("Atoms excluded from Color features") {
  auto m1 =
      "Nc1ccccc1 |(0.392086,-2.22477,0.190651;0.232269,-1.38667,0.118385;-1.06274,-0.918982,0.0342466;-1.26098,0.446053,-0.0811879;-0.244035,1.36265,-0.11691;1.05134,0.875929,-0.031248;1.28797,-0.499563,0.0864097),atomProp:0.dummyLabel.*|"_smiles;
  auto shape1 = PrepareConformer(*m1, -1);
  // The aniline N comes out as a donor, acceptor and basic so gets 3
  // features.
  CHECK(shape1.coord.size() == 33);
  ShapeInputOptions opts;
  opts.notColorAtoms = std::vector<unsigned int>{0};
  auto shape2 = PrepareConformer(*m1, -1, opts);
  CHECK(shape2.coord.size() == 24);
}

TEST_CASE("Hs not properly transformed when hcount = feature count") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";

  SECTION("as reported") {
    v2::FileParsers::MolFileParserParams ps;
    ps.removeHs = false;
    auto mol1 =
        v2::FileParsers::MolFromMolFile(dirName + "/hcount_ex1_1.mol", ps);
    REQUIRE(mol1);
    auto mol2 =
        v2::FileParsers::MolFromMolFile(dirName + "/hcount_ex1_2.mol", ps);
    REQUIRE(mol2);

    {
      RWMol cp(*mol2);
      std::vector<float> matrix(12, 0.0);
      auto [nbr_st, nbr_ct] =
          AlignMolecule(*mol1, cp, matrix, -1, -1, true, 1.0, 30, 30);
      CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.911, 0.005));
      CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.555, 0.005));

      // the bug led to H atoms in stupid positions, so we can detect it by just
      // looking at bond lengths to Hs:
      for (auto i = cp.getNumHeavyAtoms(); i < cp.getNumAtoms(); ++i) {
        INFO("checking atom " << i);
        auto at = cp.getAtomWithIdx(i);
        for (auto nbr : cp.atomNeighbors(at)) {
          auto dist = (cp.getConformer().getAtomPos(i) -
                       cp.getConformer().getAtomPos(nbr->getIdx()))
                          .length();
          CHECK(dist < 1.2);  // should be a bond to H
        }
      }

      MolToMolFile(*mol1, "m1_out.mol");
      MolToMolFile(cp, "m2_out.mol");
    }
  }
}

TEST_CASE("custom feature points") {
  auto m1 =
      "O=CC=O |(-1.75978,0.148897,0;-0.621382,-0.394324,0;0.624061,0.3656,.1;1.7571,-0.120174,.1)|"_smiles;
  SECTION("using shapes") {
    auto shape1 = PrepareConformer(*m1, -1);
    // each carbonyl O gets one feature:
    CHECK(shape1.coord.size() == 18);
    ShapeInputOptions opts2;
    opts2.customFeatures = {{1, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0},
                            {2, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0}};
    auto shape2 = PrepareConformer(*m1, -1, opts2);
    CHECK(shape2.coord.size() == 18);

    {
      // confirm that we don't add the features if useColors is false
      ShapeInputOptions topts;
      topts.customFeatures = {
          {1, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0},
          {2, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0}};
      topts.useColors = false;
      auto tshape = PrepareConformer(*m1, -1, topts);
      CHECK(tshape.coord.size() == 12);
    }

    // we'll swap the features on the second shape so that the alignment has to
    // be inverted
    ShapeInputOptions opts3;
    opts3.customFeatures = {{2, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0},
                            {1, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0}};

    auto m2 = ROMol(*m1);
    auto shape3 = PrepareConformer(m2, -1, opts3);
    CHECK(shape3.coord.size() == 18);
    double opt_param = 0.5;
    std::vector<float> matrix(12, 0.0);
    auto [st, ct] = AlignShape(shape2, shape3, matrix, opt_param);
    CHECK_THAT(st, Catch::Matchers::WithinAbs(1.000, 0.001));
    CHECK_THAT(ct, Catch::Matchers::WithinAbs(0.997, 0.001));
    CHECK(shape3.coord[0] > 0);      // x coord of first atom
    CHECK(shape3.coord[3 * 3] < 0);  // x coord of fourth atom

    auto conf = m2.getConformer(-1);
    TransformConformer(shape2.shift, shape2.inertialRot, matrix, shape3, conf);
    CHECK(conf.getAtomPos(0).x > 0);
    CHECK(conf.getAtomPos(3).x < 0);
  }
  SECTION("using molecules") {
    ShapeInputOptions opts2;
    opts2.customFeatures = {{1, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0},
                            {2, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0}};

    auto m2 = ROMol(*m1);
    // we'll swap the features on the second shape so that the alignment has to
    // be inverted
    ShapeInputOptions opts3;
    opts3.customFeatures = {{2, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0},
                            {1, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0}};

    double opt_param = 0.5;
    std::vector<float> matrix(12, 0.0);
    auto [st, ct] =
        AlignMolecule(*m1, m2, matrix, opts2, opts3, -2, -2, opt_param);
    CHECK_THAT(st, Catch::Matchers::WithinAbs(1.000, 0.001));
    CHECK_THAT(ct, Catch::Matchers::WithinAbs(0.997, 0.001));
    auto conf = m2.getConformer(-1);
    CHECK(conf.getAtomPos(0).x > 0);
    CHECK(conf.getAtomPos(3).x < 0);
  }
}

TEST_CASE("Score No Overlay") {
  // These are 2 ligands used by Andy Grant and Co in their original paper
  // https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(19961115)17:14%3C1653::AID-JCC7%3E3.0.CO;2-K
  // Ligands as extracted from PDB, with a bit of munging to get them as
  // SMILES strings (downloaded the Ideal ligand structures from RCSB
  // as SDFs and transferred the corresponding atom coords from 3tmn and 1tmn).
  auto pdb_trp_3tmn =
      R"(N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2 |(37.935,40.394,-3.825;39.119,39.593,-4.13;38.758,38.486,-5.101;37.526,38.337,-5.395;39.716,37.852,-5.605;39.883,39.108,-2.906;39.086,38.098,-2.209;38.093,38.363,-1.34;37.565,37.179,-0.881;38.201,36.136,-1.471;39.193,36.684,-2.308;40.015,35.812,-3.036;39.846,34.441,-2.913;38.844,33.933,-2.075;38.015,34.752,-1.333),wU:1.0|)"_smiles;
  REQUIRE(pdb_trp_3tmn);
  auto pdb_0zn_1tmn =
      R"([C@H](CCc1ccccc1)(C(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2)CC(C)C |(35.672,41.482,-5.722;34.516,40.842,-6.512;34.843,39.355,-6.7;33.819,38.475,-7.45;33.825,38.414,-8.838;32.951,37.553,-9.53;32.064,36.747,-8.81;32.096,36.799,-7.402;32.985,37.656,-6.73;35.934,42.778,-6.452;36.833,42.858,-7.316;35.175,43.735,-6.275;35.516,41.561,-4.218;36.707,42.096,-3.513;38.055,41.449,-3.859;39.11,42.138,-3.959;37.975,40.129,-3.983;39.134,39.277,-4.298;38.825,38.04,-5.133;37.649,37.934,-5.605;39.788,37.369,-5.652;39.985,38.945,-3.037;39.221,37.953,-2.164;37.934,37.961,-1.823;37.579,36.695,-1.314;38.63,35.975,-1.286;39.736,36.771,-1.642;41.052,36.341,-1.48;41.213,35.042,-0.964;40.095,34.215,-0.69;38.765,34.665,-0.855;36.506,41.966,-2.002;37.6,42.757,-1.31;37.546,44.225,-1.728;37.408,42.58,0.19),wD:0.0,wU:17.21,13.33|)"_smiles;
  REQUIRE(pdb_0zn_1tmn);
  ShapeInputOptions mol1Opts, mol2Opts;
  {
    auto [singleShape, singleColor] =
        ScoreMolecule(*pdb_trp_3tmn, *pdb_trp_3tmn, mol1Opts, mol2Opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(1.0, 0.001));
  }
  {
    auto pdb_trp_3tmn_cp =
        R"(N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2 |(37.935,40.394,-3.825;39.119,39.593,-4.13;38.758,38.486,-5.101;37.526,38.337,-5.395;39.716,37.852,-5.605;39.883,39.108,-2.906;39.086,38.098,-2.209;38.093,38.363,-1.34;37.565,37.179,-0.881;38.201,36.136,-1.471;39.193,36.684,-2.308;40.015,35.812,-3.036;39.846,34.441,-2.913;38.844,33.933,-2.075;38.015,34.752,-1.333),wU:1.0|)"_smiles;
    RDGeom::Point3D trans{100.0, 100.0, 100.0};
    RDGeom::Transform3D transform_3d;
    transform_3d.SetTranslation(trans);
    MolTransforms::transformConformer(pdb_trp_3tmn_cp->getConformer(),
                                      transform_3d);

    auto [singleShape, singleColor] =
        ScoreMolecule(*pdb_trp_3tmn, *pdb_trp_3tmn_cp, mol1Opts, mol2Opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.0, 0.001));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.0, 0.001));
  }
  {
    auto [singleShape, singleColor] =
        ScoreMolecule(*pdb_0zn_1tmn, *pdb_0zn_1tmn, mol1Opts, mol2Opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(1.0, 0.001));
  }
  {
    auto [singleShape, singleColor] =
        ScoreMolecule(*pdb_trp_3tmn, *pdb_0zn_1tmn, mol1Opts, mol2Opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.3376, 0.0005));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.2674, 0.0005));
  }
  {
    ShapeInputOptions opts;
    opts.normalize = false;
    auto shape = PrepareConformer(*pdb_trp_3tmn, -1, opts);
    auto [singleShape, singleColor] =
        ScoreMolecule(shape, *pdb_0zn_1tmn, mol1Opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.3376, 0.0005));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.2674, 0.0005));
  }
  {
    ShapeInputOptions opts;
    opts.normalize = false;
    auto shape1 = PrepareConformer(*pdb_trp_3tmn, -1, opts);
    auto shape2 = PrepareConformer(*pdb_0zn_1tmn, -1, opts);
    auto [singleShape, singleColor] = ScoreShape(shape1, shape2, true);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.3376, 0.0005));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.2674, 0.0005));
  }
}

TEST_CASE("Iressa onto Tagrisso") {
  // Conformations from PubChem produced by Omega. Iressa rotated and translated
  // by a random amount.  PubChem puts them both in their inertial frame which
  // makes things too easy.
  auto tagrisso =
      "C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C |(-0.9161,3.8415,-2.9811;0.1848,3.1933,-2.588;0.1064,1.7789,-2.12;-0.9619,1.1797,-2.0847;1.3654,1.2872,-1.7553;1.6841,0.0144,-1.273;0.6638,-0.9235,-1.1146;0.9578,-2.1997,-0.6343;-0.0813,-3.1358,-0.4783;-1.4556,-2.9979,-0.1847;-2.1716,-4.1359,-0.1085;-3.4803,-3.9673,0.173;-4.0689,-2.7353,0.3728;-3.2269,-1.647,0.2676;-3.7311,-0.317,0.4568;-5.0275,0.0291,0.153;-5.1887,1.3569,0.4454;-6.4231,2.0889,0.2595;-4.0141,1.8811,0.9361;-3.7121,3.1796,1.3615;-2.4139,3.4249,1.8179;-1.4588,2.4106,1.8467;-1.7752,1.1164,1.4181;-3.0776,0.8453,0.9537;-1.9103,-1.7423,-0.011;2.2723,-2.5382,-0.3127;2.58,-3.7798,0.1575;2.539,-3.9651,1.571;3.2927,-1.6003,-0.4713;2.9986,-0.324,-0.9514;4.0475,0.61,-1.1047;4.7738,0.6956,-2.3546;4.4021,1.497,-0.0162;5.4401,0.8254,0.8736;5.8294,1.7155,1.9601;4.8213,1.7057,3.0218;7.1361,1.3324,2.4981)|"_smiles;
  REQUIRE(tagrisso);
  auto iressa =
      "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1 |(11.4672,-0.467948,5.63989;12.0133,0.532631,6.49693;11.2039,1.5801,6.81985;11.2014,2.71958,6.00975;10.3926,3.81652,6.29699;10.4038,4.90395,5.50623;9.58889,5.91871,5.85946;8.76443,5.96486,6.91838;8.77814,4.86059,7.68868;7.92337,4.86224,8.81914;7.44878,5.8925,9.64622;8.22182,7.03851,9.85619;7.75051,8.06265,10.6777;6.50441,7.94546,11.2936;6.06567,8.93802,12.0809;5.72932,6.80403,11.0875;4.19056,6.65372,11.8447;6.20047,5.78015,10.2656;9.57161,3.74547,7.43436;9.56851,2.60328,8.25407;10.3868,1.52933,7.93769;10.3797,0.419365,8.74203;11.3064,0.402096,9.81907;10.7104,-0.399165,10.9685;9.40938,0.22121,11.4678;9.64205,1.59049,11.9223;8.38006,2.22985,12.3199;8.64991,3.65266,12.8011;9.56883,3.64192,13.8942;10.8103,3.05101,13.5078;10.5931,1.61394,13.0425)|"_smiles;
  REQUIRE(iressa);
  std::vector<float> matrix(12, 0.0);
  auto sims =
      AlignMolecule(*tagrisso, *iressa, matrix, -1, -1, true, 0.5, 10, 30);
  CHECK_THAT(sims.first, Catch::Matchers::WithinAbs(0.582, 0.005));
  CHECK_THAT(sims.second, Catch::Matchers::WithinAbs(0.092, 0.005));
}

TEST_CASE("PCLOBSTER") {
  std::string lobster_file =
      "/home/dave/Projects/Lobster/LOBSTER_112024/all_ligands.sdf";
  auto suppl = SDMolSupplier(lobster_file);
  std::vector<std::shared_ptr<ROMol>> mols;
  while (!suppl.atEnd()) {
    auto mol = suppl.next();
    mols.emplace_back(mol);
  }
  std::cout << "Number of mols " << mols.size() << std::endl;
  double sum_st = 0.0, sum_ct = 0.0;
  int num = 0;
  std::mt19937 e2(1);
  std::uniform_real_distribution<double> unif(0, 1);
  std::vector<float> matrix(12, 0.0);
  for (size_t i = 1; i < mols.size(); i++) {
    for (size_t j = 0; j < i; j++) {
      if (unif(e2) > 0.001) {
        continue;
      }
      auto [st, ct] = AlignMolecule(*mols[i], *mols[j], matrix);
      sum_st += st;
      sum_ct += ct;
      ++num;
      if (!(num % 1000)) {
        std::cout << num << "  " << i << "  " << j << std::endl;
      }
    }
  }
  std::cout << "Mean st of " << num << " : " << sum_st / num << std::endl;
  std::cout << "Mean ct of " << num << " : " << sum_ct / num << std::endl;
}
