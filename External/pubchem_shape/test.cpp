#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/RWMol.h>

#include "PubChemShape.hpp"

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
      CHECK_THAT(cp.getConformer().getAtomPos(i).x,
                 Catch::Matchers::WithinAbs(cp2.getConformer().getAtomPos(i).x,
                                            0.005));
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
  SECTION("basics") {
    RWMol cp(*probe);
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] =
        AlignMolecule(*ref, cp, matrix, -1, -1, test_use_colors, test_opt_param,
                      test_max_preiters, test_max_postiters);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.837, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.694, 0.005));
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
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.501, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.107, 0.005));
    RWMol cp2(cp);
    auto [nbr_st2, nbr_ct2] =
        AlignMolecule(*ref, cp2, matrix, refConfId, prbConfId, useColors,
                      opt_param, preiters, postiters);
    CHECK_THAT(nbr_st2, Catch::Matchers::WithinAbs(nbr_st2, 0.005));
    CHECK_THAT(nbr_ct2, Catch::Matchers::WithinAbs(nbr_ct2, 0.005));

    auto rmsd = MolAlign::CalcRMS(cp, cp2);
    CHECK_THAT(rmsd, Catch::Matchers::WithinAbs(0.017, 0.005));
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
  // molecule onto the reference and the probe shape onto the reference shape

  CHECK_THAT(shape_st, Catch::Matchers::WithinAbs(mol_st, 0.001));
  CHECK_THAT(shape_ct, Catch::Matchers::WithinAbs(mol_ct, 0.001));
  for (unsigned int i = 0; i < probe->getNumAtoms(); i++) {
    const auto &pos = probe->getConformer().getAtomPos(i);
    CHECK_THAT(pos.x, Catch::Matchers::WithinAbs(
                          probeShape.coord[3 * i] - refShape.shift[0], 0.001));
    CHECK_THAT(pos.y,
               Catch::Matchers::WithinAbs(
                   probeShape.coord[3 * i + 1] - refShape.shift[1], 0.001));
    CHECK_THAT(pos.z,
               Catch::Matchers::WithinAbs(
                   probeShape.coord[3 * i + 2] - refShape.shift[2], 0.001));
  }

  // Also check the TransformConformer function
  TransformConformer(refShape.shift, matrix, probeShape,
                     probeCP.getConformer(-1));
  for (unsigned int i = 0; i < probe->getNumAtoms(); i++) {
    const auto &pos = probe->getConformer().getAtomPos(i);
    const auto &poscp = probeCP.getConformer().getAtomPos(i);
    CHECK_THAT(pos.x, Catch::Matchers::WithinAbs(poscp.x, 0.001));
    CHECK_THAT(pos.y, Catch::Matchers::WithinAbs(poscp.y, 0.001));
    CHECK_THAT(pos.z, Catch::Matchers::WithinAbs(poscp.z, 0.001));
  }
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
    CHECK_THAT(st, Catch::Matchers::WithinAbs(0.997, 0.001));
    CHECK_THAT(ct, Catch::Matchers::WithinAbs(0.978, 0.001));
    CHECK(shape3.coord[0] > 0);      // x coord of first atom
    CHECK(shape3.coord[3 * 3] < 0);  // x coord of fourth atom

    auto conf = m2.getConformer(-1);
    TransformConformer(shape2.shift, matrix, shape3, conf);
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
    CHECK_THAT(st, Catch::Matchers::WithinAbs(0.997, 0.001));
    CHECK_THAT(ct, Catch::Matchers::WithinAbs(0.978, 0.001));
    auto conf = m2.getConformer(-1);
    CHECK(conf.getAtomPos(0).x > 0);
    CHECK(conf.getAtomPos(3).x < 0);
  }
}
