//
//  Copyright (C) 2021-2025 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//
// Tests for the Roshambo2-based shape alignment code.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "Roshambo2Shape.hpp"
#include "ShapeInput.h"
#include "GraphMol/FileParsers/MolWriters.h"
#include "GraphMol/SmilesParse/SmilesWrite.h"

using namespace RDKit;

TEST_CASE("basic alignment") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";

  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/test1.sdf");
  auto ref = suppl[0];
  REQUIRE(ref);
  auto probe = suppl[1];
  REQUIRE(probe);
  ShapeAlign::ShapeOverlayOptions overlayOpts;
  overlayOpts.d_useColors = true;
  SECTION("setup") {
    auto refShape = ShapeAlign::ShapeInput(*ref, -1, overlayOpts);
    CHECK_THAT(refShape.getSelfOverlapVol(),
               Catch::Matchers::WithinAbs(591.073, 0.005));
    CHECK_THAT(refShape.getSelfOverlapColor(),
               Catch::Matchers::WithinAbs(16.114, 0.005));

    auto probeShape = ShapeAlign::ShapeInput(*probe, -1, overlayOpts);
    CHECK_THAT(probeShape.getSelfOverlapVol(),
               Catch::Matchers::WithinAbs(751.067, 0.005));
    CHECK_THAT(probeShape.getSelfOverlapColor(),
               Catch::Matchers::WithinAbs(30.026, 0.005));
  }

  SECTION("with color") {
    ROMol cp(*probe);
    const auto &[tsc, tcc] =
        ShapeAlign::AlignMolecule(*ref, cp, nullptr, overlayOpts);
    CHECK_THAT(tsc, Catch::Matchers::WithinAbs(0.702, 0.005));
    CHECK_THAT(tcc, Catch::Matchers::WithinAbs(0.556, 0.005));
  }

  SECTION("without color") {
    ROMol cp(*probe);
    overlayOpts.d_useColors = false;
    const auto &[tsnc, tcnc] =
        ShapeAlign::AlignMolecule(*ref, cp, nullptr, overlayOpts);
    CHECK_THAT(tsnc, Catch::Matchers::WithinAbs(0.761, 0.005));
    CHECK_THAT(tcnc, Catch::Matchers::WithinAbs(0.0, 0.005));
  }

  SECTION("collect transform") {
    overlayOpts.d_useColors = true;
    ROMol cp(*probe);
    RDGeom::Transform3D xform;
    const auto &[tsnc, tcnc] =
        ShapeAlign::AlignMolecule(*ref, cp, &xform, overlayOpts);
    CHECK_THAT(tsnc, Catch::Matchers::WithinAbs(0.702, 0.005));
    CHECK_THAT(tcnc, Catch::Matchers::WithinAbs(0.556, 0.005));
    // Check a few values from the transform, just to be sure
    CHECK_THAT(xform.getValUnchecked(0, 0),
               Catch::Matchers::WithinAbs(-0.864, 0.005));
    CHECK_THAT(xform.getValUnchecked(1, 1),
               Catch::Matchers::WithinAbs(-0.804, 0.005));
    CHECK_THAT(xform.getValUnchecked(2, 2),
               Catch::Matchers::WithinAbs(0.686, 0.005));
    CHECK_THAT(xform.getValUnchecked(3, 3),
               Catch::Matchers::WithinAbs(1.0, 0.005));
  }
}

TEST_CASE("bulk") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";
  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/bulk.pubchem.sdf");
  auto ref = suppl[0];
  REQUIRE(ref);
  std::string testout = dirName + "/bulk.pubchem_out.sdf";
  auto writer = SDWriter(testout);
  writer.write(*ref);
  for (auto i = 1u; i < suppl.length(); ++i) {
    auto probe = suppl[1];
    REQUIRE(probe);
    auto [nbr_st, nbr_ct] = ShapeAlign::AlignMolecule(*ref, *probe);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.811, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.525, 0.005));
    writer.write(*probe);
  }
  writer.close();
}

TEST_CASE("shape alignment") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";

  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/test1.sdf");
  auto ref = suppl[0];
  REQUIRE(ref);
  auto probe = suppl[1];
  REQUIRE(probe);
  auto refShape = ShapeAlign::ShapeInput(*ref, -1);
  auto probeShape = ShapeAlign::ShapeInput(*probe, -1);

  auto [st, ct] = ShapeAlign::AlignShape(refShape, probeShape, nullptr);
  CHECK_THAT(st, Catch::Matchers::WithinAbs(0.702, 0.005));
  CHECK_THAT(ct, Catch::Matchers::WithinAbs(0.556, 0.005));

  // And pre-normalizing the shapes
  refShape.normalizeCoords();
  auto [st1, ct1] = ShapeAlign::AlignShape(refShape, probeShape, nullptr);
  CHECK_THAT(st1, Catch::Matchers::WithinAbs(0.702, 0.005));
  CHECK_THAT(ct1, Catch::Matchers::WithinAbs(0.556, 0.005));
  probeShape.normalizeCoords();
  auto [st2, ct2] = ShapeAlign::AlignShape(refShape, probeShape, nullptr);
  CHECK_THAT(st2, Catch::Matchers::WithinAbs(0.702, 0.005));
  CHECK_THAT(ct2, Catch::Matchers::WithinAbs(0.556, 0.005));
}

TEST_CASE("Overlay onto shape bug (Github8462)") {
  auto m1 =
      R"(c1ccc(-c2ccccc2)cc1 |(-3.26053,-0.0841607,-0.741909;-2.93383,0.123873,0.593407;-1.60713,0.377277,0.917966;-0.644758,0.654885,-0.0378428;0.743308,0.219134,0.168663;1.82376,1.0395,-0.0112769;3.01462,0.695405,0.613858;3.18783,-0.589771,1.09649;2.15761,-1.50458,1.01949;0.988307,-1.1313,0.385783;-1.1048,0.797771,-1.34022;-2.39754,0.435801,-1.69921)|)"_smiles;
  REQUIRE(m1);
  ROMol m2(*m1);
  for (auto a : m2.atoms()) {
    auto &pos = m2.getConformer().getAtomPos(a->getIdx());
    pos.x += 3.0;
    pos.y += 2.0;
  }
  ROMol m3(m2);

  auto [st, ct] = ShapeAlign::AlignMolecule(*m1, m2);
  CHECK_THAT(st, Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(ct, Catch::Matchers::WithinAbs(1.0, 0.005));
  for (unsigned int i = 0; i < m1->getNumAtoms(); ++i) {
    auto pos1 = m1->getConformer().getAtomPos(i);
    auto pos2 = m2.getConformer().getAtomPos(i);
    CHECK_THAT((pos1 - pos2).length(), Catch::Matchers::WithinAbs(0.0, 0.005));
  }

  // Create the shape without normalization to mimic an arbitrary shape.
  ShapeAlign::ShapeOverlayOptions opts;
  opts.d_normalize = false;
  auto s1 = ShapeAlign::ShapeInput(*m1, -1, opts);
  auto [st1, ct1] = AlignMolecule(s1, m3);
  CHECK_THAT(st1, Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(ct1, Catch::Matchers::WithinAbs(1.0, 0.005));
  for (unsigned int i = 0; i < m3.getNumAtoms(); ++i) {
    RDGeom::Point3D pos1(s1.getCoords()[4 * i], s1.getCoords()[4 * i + 1],
                         s1.getCoords()[4 * i + 2]);
    auto pos2 = m3.getConformer().getAtomPos(i);
    CHECK_THAT((pos1 - pos2).length(), Catch::Matchers::WithinAbs(0.0, 0.01));
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
    RDGeom::Transform3D xform;
    auto [nbr_st, nbr_ct] = ShapeAlign::AlignMolecule(*ref, cp);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.828, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.994, 0.005));
    for (auto i = 0u; i < cp.getNumAtoms(); ++i) {
      // the failure mode here was that Hs had HUGE coordinates
      auto pos = cp.getConformer().getAtomPos(i);
      CHECK((pos.x > -10 && pos.x < 10));
    }
    // Check the fixed calculation which will not be identical but should be
    // very close.
    auto news = ShapeAlign::ScoreMolecule(*ref, cp);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(news.first, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(news.second, 0.005));
  }
}

TEST_CASE("Github #8096") {
  SECTION("as reported") {
    auto m1 =
        R"([H]c1c([H])c([H])c([2H])c([H])c1[H] |(1.55967,1.91617,0.0546381;0.885536,1.07172,0.030849;1.38172,-0.23747,0.0274262;2.44539,-0.439501,0.0483424;0.470206,-1.27516,-0.00361916;0.856925,-2.30002,-0.00633525;-0.896665,-1.07227,-0.0310991;-1.60071,-1.87642,-0.0551085;-1.36315,0.22877,-0.0271173;-2.43593,0.379132,-0.0487835;-0.479018,1.29083,0.00359778;-0.823965,2.31421,0.00720933)|)"_smiles;
    REQUIRE(m1);
    auto m2 =
        R"([H]c1c([H])c([H])c([H])c([H])c1[H] |(-2.06264,-0.844763,-0.0261403;-1.04035,-0.481453,-0.0114878;-0.00743655,-1.41861,-0.0137121;-0.215455,-2.47997,-0.0295909;1.29853,-0.949412,0.00507497;2.12524,-1.65277,0.00390664;1.58501,0.395878,0.0254188;2.61997,0.704365,0.0394811;0.550242,1.31385,0.0273741;0.783172,2.37039,0.0434262;-0.763786,0.88847,0.00908113;-1.60557,1.58532,0.0100194)|)"_smiles;
    REQUIRE(m2);
    auto [nbr_st, nbr_ct] = ShapeAlign::AlignMolecule(*m1, *m2);
    CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.994, 0.005));
    CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(1.0, 0.005));
  }
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
      auto [nbr_st, nbr_ct] = ShapeAlign::AlignMolecule(*mol1, cp);
      CHECK_THAT(nbr_st, Catch::Matchers::WithinAbs(0.901, 0.005));
      CHECK_THAT(nbr_ct, Catch::Matchers::WithinAbs(0.860, 0.005));
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
    }
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
  ShapeAlign::ShapeOverlayOptions mol1Opts;
  {
    auto [singleShape, singleColor] =
        ScoreMolecule(*pdb_trp_3tmn, *pdb_trp_3tmn, mol1Opts);
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
        ScoreMolecule(*pdb_trp_3tmn, *pdb_trp_3tmn_cp, mol1Opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.0, 0.001));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.0, 0.001));
  }
  {
    auto [singleShape, singleColor] =
        ScoreMolecule(*pdb_0zn_1tmn, *pdb_0zn_1tmn, mol1Opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(1.0, 0.001));
  }
  {
    auto [singleShape, singleColor] =
        ScoreMolecule(*pdb_trp_3tmn, *pdb_0zn_1tmn, mol1Opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.349, 0.005));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.342, 0.005));
  }
  {
    ShapeAlign::ShapeOverlayOptions opts;
    opts.d_normalize = false;
    auto shape = ShapeAlign::ShapeInput(*pdb_trp_3tmn, -1, opts);
    auto [singleShape, singleColor] =
        ScoreMolecule(shape, *pdb_0zn_1tmn, mol1Opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.349, 0.005));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.342, 0.005));
  }
  {
    ShapeAlign::ShapeOverlayOptions opts;
    opts.d_normalize = false;
    auto shape1 = ShapeAlign::ShapeInput(*pdb_trp_3tmn, -1, opts);
    auto shape2 = ShapeAlign::ShapeInput(*pdb_0zn_1tmn, -1, opts);
    auto [singleShape, singleColor] = ShapeAlign::ScoreShape(shape1, shape2);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.349, 0.005));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.342, 0.005));
  }
}

TEST_CASE("Iressa onto Tagrisso") {
  // Conformations from PubChem produced by Omega. Iressa rotated and translated
  // by a random amount.  PubChem puts them both in their inertial frame which
  // makes things too easy.
  auto tagrisso =
      R"(C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C |(-0.9161,3.8415,-2.9811;0.1848,3.1933,-2.588;0.1064,1.7789,-2.12;-0.9619,1.1797,-2.0847;1.3654,1.2872,-1.7553;1.6841,0.0144,-1.273;0.6638,-0.9235,-1.1146;0.9578,-2.1997,-0.6343;-0.0813,-3.1358,-0.4783;-1.4556,-2.9979,-0.1847;-2.1716,-4.1359,-0.1085;-3.4803,-3.9673,0.173;-4.0689,-2.7353,0.3728;-3.2269,-1.647,0.2676;-3.7311,-0.317,0.4568;-5.0275,0.0291,0.153;-5.1887,1.3569,0.4454;-6.4231,2.0889,0.2595;-4.0141,1.8811,0.9361;-3.7121,3.1796,1.3615;-2.4139,3.4249,1.8179;-1.4588,2.4106,1.8467;-1.7752,1.1164,1.4181;-3.0776,0.8453,0.9537;-1.9103,-1.7423,-0.011;2.2723,-2.5382,-0.3127;2.58,-3.7798,0.1575;2.539,-3.9651,1.571;3.2927,-1.6003,-0.4713;2.9986,-0.324,-0.9514;4.0475,0.61,-1.1047;4.7738,0.6956,-2.3546;4.4021,1.497,-0.0162;5.4401,0.8254,0.8736;5.8294,1.7155,1.9601;4.8213,1.7057,3.0218;7.1361,1.3324,2.4981)|)"_smiles;
  REQUIRE(tagrisso);
  auto iressa =
      R"(COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1 |(11.4672,-0.467948,5.63989;12.0133,0.532631,6.49693;11.2039,1.5801,6.81985;11.2014,2.71958,6.00975;10.3926,3.81652,6.29699;10.4038,4.90395,5.50623;9.58889,5.91871,5.85946;8.76443,5.96486,6.91838;8.77814,4.86059,7.68868;7.92337,4.86224,8.81914;7.44878,5.8925,9.64622;8.22182,7.03851,9.85619;7.75051,8.06265,10.6777;6.50441,7.94546,11.2936;6.06567,8.93802,12.0809;5.72932,6.80403,11.0875;4.19056,6.65372,11.8447;6.20047,5.78015,10.2656;9.57161,3.74547,7.43436;9.56851,2.60328,8.25407;10.3868,1.52933,7.93769;10.3797,0.419365,8.74203;11.3064,0.402096,9.81907;10.7104,-0.399165,10.9685;9.40938,0.22121,11.4678;9.64205,1.59049,11.9223;8.38006,2.22985,12.3199;8.64991,3.65266,12.8011;9.56883,3.64192,13.8942;10.8103,3.05101,13.5078;10.5931,1.61394,13.0425)|)"_smiles;
  REQUIRE(iressa);
  auto sims = ShapeAlign::AlignMolecule(*tagrisso, *iressa);
  CHECK_THAT(sims.first, Catch::Matchers::WithinAbs(0.571, 0.005));
  CHECK_THAT(sims.second, Catch::Matchers::WithinAbs(0.360, 0.005));
}

TEST_CASE("Optimise in place") {
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
  auto initScores = ShapeAlign::ScoreMolecule(*pdb_trp_3tmn, *pdb_0zn_1tmn);
  CHECK_THAT(initScores.first, Catch::Matchers::WithinAbs(0.349, 0.001));
  CHECK_THAT(initScores.second, Catch::Matchers::WithinAbs(0.342, 0.001));
  {
    ShapeAlign::ShapeOverlayOptions opts;
    opts.d_mode = ShapeAlign::StartMode::AS_IS;
    opts.d_normalize = false;
    ROMol cp(*pdb_0zn_1tmn);
    auto [singleShape, singleColor] =
        ShapeAlign::AlignMolecule(*pdb_trp_3tmn, cp, nullptr, opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.404, 0.001));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.354, 0.001));
  }
  {
    // With default settings, it does a poor job.
    ShapeAlign::ShapeOverlayOptions opts;
    ROMol cp(*pdb_0zn_1tmn);
    auto [singleShape, singleColor] =
        ShapeAlign::AlignMolecule(*pdb_trp_3tmn, cp, nullptr, opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.288, 0.001));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.383, 0.001));
  }
  {
    // And with shapes the same
    ShapeAlign::ShapeOverlayOptions opts;
    opts.d_mode = ShapeAlign::StartMode::AS_IS;
    opts.d_normalize = false;
    auto refShape = ShapeAlign::ShapeInput(*pdb_trp_3tmn, -1, opts);
    ROMol cp(*pdb_0zn_1tmn);
    auto [singleShape, singleColor] =
        ShapeAlign::AlignMolecule(refShape, cp, nullptr, opts);
    CHECK_THAT(singleShape, Catch::Matchers::WithinAbs(0.404, 0.001));
    CHECK_THAT(singleColor, Catch::Matchers::WithinAbs(0.354, 0.001));
  }
}