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
constexpr double test_use_colors = true;

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
