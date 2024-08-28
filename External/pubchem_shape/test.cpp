#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/RWMol.h>

#include "PubChemShape.hpp"

using namespace RDKit;

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
    auto [nbr_st, nbr_ct] = AlignMolecule(*ref, cp, matrix);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.773, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.303, 0.005));
    for (auto i = 0; i < probe->getNumAtoms(); ++i) {
      CHECK(probe->getConformer().getAtomPos(i).x !=
            cp.getConformer().getAtomPos(i).x);
    }
  }
  SECTION("from shape") {
    auto ref_shape = PrepareConformer(*ref);
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] = AlignMolecule(ref_shape, *probe, matrix);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.773, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.303, 0.005));
  }
  SECTION("RDKit features") {
    ref->clearProp("PUBCHEM_PHARMACOPHORE_FEATURES");
    probe->clearProp("PUBCHEM_PHARMACOPHORE_FEATURES");
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] = AlignMolecule(*ref, *probe, matrix);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.773, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.231, 0.005));
  }
  SECTION("no colors") {
    std::vector<float> matrix(12, 0.0);
    int refConfId = -1;
    int prbConfId = -1;
    bool useColors = false;
    auto [nbr_st, nbr_ct] =
        AlignMolecule(*ref, *probe, matrix, refConfId, prbConfId, useColors);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.773, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.000, 0.005));
  }
}

TEST_CASE("bulk") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";
  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/bulk.pubchem.sdf");
  auto ref = suppl[0];
  REQUIRE(ref);
  for (auto i = 1; i < suppl.length(); ++i) {
    auto probe = suppl[1];
    REQUIRE(probe);
    std::vector<float> matrix(12, 0.0);
    auto [nbr_st, nbr_ct] = AlignMolecule(*ref, *probe, matrix);
    CHECK_THAT(nbr_st,
               Catch::Matchers::WithinAbs(
                   probe->getProp<float>("shape_align_shape_tanimoto"), 0.005));
    CHECK_THAT(nbr_ct,
               Catch::Matchers::WithinAbs(
                   probe->getProp<float>("shape_align_color_tanimoto"), 0.005));
  }
}
