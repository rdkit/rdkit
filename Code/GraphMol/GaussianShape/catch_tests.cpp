//
//  Copyright (C) 2026 David Cosgrove and other RDKit contributors
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

#include <random>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

bool checkMolsHaveRoughlySameCoords(const ROMol &m1, const ROMol &m2,
                                    double margin = 0.005) {
  for (unsigned int i = 0; i < m1.getNumAtoms(); ++i) {
    auto pos1 = m1.getConformer().getAtomPos(i);
    auto pos2 = m2.getConformer().getAtomPos(i);
    if ((pos1 - pos2).length() > margin) {
      // So the error is printed in a relevant place.
      std::cout << i << " : " << m1.getAtomWithIdx(i)->getAtomicNum()
                << " :: " << (pos1 - pos2).length() << std::endl;
      CHECK_THAT((pos1 - pos2).length(),
                 Catch::Matchers::WithinAbs(0.0, margin));
      return false;
    }
  }
  return true;
}

TEST_CASE("basic alignment") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";

  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/test1.sdf");
  auto refT = suppl[0];
  auto ref = v2::SmilesParse::MolFromSmiles(MolToCXSmiles(*refT));
  REQUIRE(ref);
  auto probeT = suppl[1];
  auto probe = v2::SmilesParse::MolFromSmiles(MolToCXSmiles(*probeT));
  REQUIRE(probe);
  GaussianShape::ShapeOverlayOptions overlayOpts;
  overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
  overlayOpts.startMode = GaussianShape::StartMode::ROTATE_180;
  overlayOpts.nSteps = 50;
  GaussianShape::ShapeInputOptions shapeOpts;
  SECTION("setup") {
    auto refShape = GaussianShape::ShapeInput(*ref, -1, shapeOpts);
    CHECK_THAT(refShape.getShapeVolume(),
               Catch::Matchers::WithinAbs(591.057, 0.005));
    CHECK_THAT(refShape.getColorVolume(),
               Catch::Matchers::WithinAbs(31.935, 0.005));

    auto probeShape = GaussianShape::ShapeInput(*probe, -1, shapeOpts);
    CHECK_THAT(probeShape.getShapeVolume(),
               Catch::Matchers::WithinAbs(751.013, 0.005));
    CHECK_THAT(probeShape.getColorVolume(),
               Catch::Matchers::WithinAbs(42.530, 0.005));
  }
  SECTION("shape only") {
    ROMol cp(*probe);
    overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_ONLY;
    overlayOpts.startMode = GaussianShape::StartMode::ROTATE_180;

    GaussianShape::ShapeInputOptions tShapeOpts;
    tShapeOpts.useColors = false;
    const auto scores = GaussianShape::AlignMolecule(
        *ref, cp, tShapeOpts, tShapeOpts, nullptr, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.760, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.760, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.0, 0.005));
    // Check that a re-score gives the same answer.
    auto rescores = GaussianShape::ScoreMolecule(*ref, cp, shapeOpts, shapeOpts,
                                                 overlayOpts);
    CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
    CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
    CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
  }
  SECTION("shape plus color score") {
    overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
    overlayOpts.startMode = GaussianShape::StartMode::ROTATE_180;
    ROMol cp(*probe);
    const auto scores = GaussianShape::AlignMolecule(
        *ref, cp, shapeOpts, shapeOpts, nullptr, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.494, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.760, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.236, 0.005));
    // Check that a re-score gives the same answer.
    const auto rescores = GaussianShape::ScoreMolecule(*ref, cp, shapeOpts,
                                                       shapeOpts, overlayOpts);
    CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
    CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
    CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
  }
  SECTION("shape and color") {
    overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR;
    overlayOpts.startMode = GaussianShape::StartMode::ROTATE_180;
    ROMol cp(*probe);
    const auto scores = GaussianShape::AlignMolecule(
        *ref, cp, shapeOpts, shapeOpts, nullptr, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.502, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.731, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.273, 0.005));
    const auto rescores = GaussianShape::ScoreMolecule(*ref, cp, shapeOpts,
                                                       shapeOpts, overlayOpts);
    CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
    CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
    CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
  }
  SECTION("collect transform") {
    overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
    overlayOpts.startMode = GaussianShape::StartMode::ROTATE_180;
    ROMol cp(*probe);
    RDGeom::Transform3D xform;
    const auto scores = GaussianShape::AlignMolecule(
        *ref, cp, shapeOpts, shapeOpts, &xform, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.494, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.760, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.236, 0.005));
    // Check a few values from the transform, just to be sure
    CHECK_THAT(xform.getValUnchecked(0, 0),
               Catch::Matchers::WithinAbs(-0.879, 0.005));
    CHECK_THAT(xform.getValUnchecked(1, 1),
               Catch::Matchers::WithinAbs(-0.818, 0.005));
    CHECK_THAT(xform.getValUnchecked(2, 2),
               Catch::Matchers::WithinAbs(0.811, 0.005));
    CHECK_THAT(xform.getValUnchecked(3, 3),
               Catch::Matchers::WithinAbs(1.0, 0.005));
  }
  SECTION("shape plus color score a la pubchem") {
    overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
    overlayOpts.startMode = GaussianShape::StartMode::A_LA_PUBCHEM;
    GaussianShape::ShapeInputOptions shapeOpts2;
    for (const auto acr : std::vector{true, false}) {
      shapeOpts2.allCarbonRadii = acr;
      ROMol cp(*probe);
      const auto scores = GaussianShape::AlignMolecule(
          *ref, cp, shapeOpts2, shapeOpts2, nullptr, overlayOpts);
      if (acr) {
        CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.498, 0.005));
        CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.758, 0.005));
        CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.237, 0.005));
      } else {
        CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.497, 0.005));
        CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.761, 0.005));
        CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.234, 0.005));
      }
      // Check that a re-score gives the same answer.
      const auto rescores = GaussianShape::ScoreMolecule(
          *ref, cp, shapeOpts2, shapeOpts2, overlayOpts);
      CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
      CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
      CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
    }
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
  GaussianShape::ShapeOverlayOptions overlayOpts;
  overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
  overlayOpts.startMode = GaussianShape::StartMode::A_LA_PUBCHEM;
  GaussianShape::ShapeInputOptions shapeOpts;
  for (auto i = 1u; i < suppl.length(); ++i) {
    auto probe = suppl[1];
    REQUIRE(probe);
    auto scores = GaussianShape::AlignMolecule(*ref, *probe, shapeOpts,
                                               shapeOpts, nullptr, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.534, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.818, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.249, 0.005));
    const auto rescores = GaussianShape::ScoreMolecule(*ref, *probe);
    CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
    CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
    CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));

    writer.write(*probe);
    break;
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
  auto refShape = GaussianShape::ShapeInput(*ref, -1);
  auto probeShape = GaussianShape::ShapeInput(*probe, -1);

  const auto ovProbe =
      "FC1(F)C[C@H](C(O)=O)N(Cc2ocnc2)C1 |(-13.8851,-5.77603,4.50901;-13.6126,-6.61984,3.48692;-12.5692,-7.38618,3.85415;-13.3182,-5.84097,2.22308;-14.6769,-5.73477,1.52809;-15.1634,-4.32173,1.4896;-14.7717,-3.71164,0.344426;-15.8263,-3.78355,2.36382;-15.6188,-6.489,2.36457;-16.6719,-7.14028,1.5981;-17.7523,-6.17224,1.19507;-18.3875,-5.49676,2.19276;-19.2867,-4.7195,1.5291;-19.2789,-4.83963,0.2207;-18.2959,-5.77046,0.000232054;-14.831,-7.44357,3.14755),wU:4.4|"_smiles;
  RDGeom::Transform3D xform;
  auto scores = GaussianShape::AlignShape(refShape, probeShape, &xform);
  CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.494, 0.005));
  CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.760, 0.005));
  CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.237, 0.005));
  SmilesWriteParams params;
  params.canonical = false;
  // The input structure being from an SDF doesn't have the atoms in an order
  // that will make a SMILES string so bounce it through one for comparison.
  auto probeCp1 = v2::SmilesParse::MolFromSmiles(MolToCXSmiles(*probe, params));
  MolTransforms::transformConformer(probeCp1->getConformer(), xform);
  CHECK(checkMolsHaveRoughlySameCoords(*ovProbe, *probeCp1));
  // And pre-normalizing the shapes
  refShape.normalizeCoords();
  probeShape.normalizeCoords();
  RDGeom::Transform3D xform1;
  auto scores1 = GaussianShape::AlignShape(refShape, probeShape, &xform1);
  CHECK_THAT(scores1[0], Catch::Matchers::WithinAbs(0.498, 0.005));
  CHECK_THAT(scores1[1], Catch::Matchers::WithinAbs(0.760, 0.005));
  CHECK_THAT(scores1[2], Catch::Matchers::WithinAbs(0.237, 0.005));
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

  auto scores = GaussianShape::AlignMolecule(*m1, m2);
  CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK(checkMolsHaveRoughlySameCoords(*m1, m2, 0.005));

  // Create the shape without normalization to mimic an arbitrary shape.
  auto s1 = GaussianShape::ShapeInput(*m1, -1);
  auto scores1 = AlignMolecule(s1, m3);
  CHECK_THAT(scores1[0], Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(scores1[1], Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(scores1[2], Catch::Matchers::WithinAbs(1.0, 0.005));
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
    auto scores = GaussianShape::AlignMolecule(*ref, cp);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.707, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.834, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.581, 0.005));
    for (auto i = 0u; i < cp.getNumAtoms(); ++i) {
      // the failure mode here was that Hs had HUGE coordinates
      auto pos = cp.getConformer().getAtomPos(i);
      CHECK((pos.x > -10 && pos.x < 10));
    }
    // Check the rescore
    auto rescores = GaussianShape::ScoreMolecule(*ref, cp);
    CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
    CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
    CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
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
    auto scores = GaussianShape::AlignMolecule(*m1, *m2);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(1.0, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.0, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(1.0, 0.005));
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
      auto scores = GaussianShape::AlignMolecule(*mol1, cp);
      CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.744, 0.005));
      CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.918, 0.005));
      CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.570, 0.005));
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
  GaussianShape::ShapeInputOptions shapeOpts;
  {
    auto scores = ScoreMolecule(*pdb_trp_3tmn, *pdb_trp_3tmn, shapeOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(1.0, 0.001));
  }
  {
    auto pdb_trp_3tmn_cp =
        R"(N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2 |(37.935,40.394,-3.825;39.119,39.593,-4.13;38.758,38.486,-5.101;37.526,38.337,-5.395;39.716,37.852,-5.605;39.883,39.108,-2.906;39.086,38.098,-2.209;38.093,38.363,-1.34;37.565,37.179,-0.881;38.201,36.136,-1.471;39.193,36.684,-2.308;40.015,35.812,-3.036;39.846,34.441,-2.913;38.844,33.933,-2.075;38.015,34.752,-1.333),wU:1.0|)"_smiles;
    RDGeom::Point3D trans{100.0, 100.0, 100.0};
    RDGeom::Transform3D transform_3d;
    transform_3d.SetTranslation(trans);
    MolTransforms::transformConformer(pdb_trp_3tmn_cp->getConformer(),
                                      transform_3d);

    auto scores = ScoreMolecule(*pdb_trp_3tmn, *pdb_trp_3tmn_cp, shapeOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.0, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.0, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.0, 0.001));
  }
  {
    auto scores = ScoreMolecule(*pdb_0zn_1tmn, *pdb_0zn_1tmn, shapeOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(1.0, 0.001));
  }
  {
    auto scores = ScoreMolecule(*pdb_trp_3tmn, *pdb_0zn_1tmn, shapeOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.307, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.349, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.265, 0.005));
  }
  {
    auto shape = GaussianShape::ShapeInput(*pdb_trp_3tmn, -1, shapeOpts);
    auto scores = ScoreMolecule(shape, *pdb_0zn_1tmn, shapeOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.307, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.349, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.265, 0.005));
  }
  {
    auto shape1 = GaussianShape::ShapeInput(*pdb_trp_3tmn, -1, shapeOpts);
    auto shape2 = GaussianShape::ShapeInput(*pdb_0zn_1tmn, -1, shapeOpts);
    auto scores = GaussianShape::ScoreShape(shape1, shape2);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.307, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.349, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.265, 0.005));
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
  GaussianShape::ShapeOverlayOptions opts;
  opts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
  opts.startMode = GaussianShape::StartMode::A_LA_PUBCHEM;
  opts.nSteps = 100;
  GaussianShape::ShapeInputOptions shapeOpts;
  shapeOpts.allCarbonRadii = false;
  auto scores = GaussianShape::AlignMolecule(*tagrisso, *iressa, shapeOpts,
                                             shapeOpts, nullptr, opts);
  CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.332, 0.005));
  CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.569, 0.005));
  CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.095, 0.005));

  auto rescores = GaussianShape::ScoreMolecule(*tagrisso, *iressa, shapeOpts,
                                               shapeOpts, opts);
  CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(rescores[0], 0.005));
  CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(rescores[1], 0.005));
  CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(rescores[2], 0.005));
  auto aligned_iressa =
      "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1 |(3.32561,-4.49721,0.300927;2.53581,-4.01343,-0.78355;1.35773,-3.4056,-0.468489;0.193742,-4.17423,-0.373537;-1.03475,-3.60055,-0.0547434;-2.14373,-4.35592,0.0319024;-3.28917,-3.71723,0.345957;-3.44662,-2.4047,0.581874;-2.31707,-1.67806,0.486599;-2.43316,-0.286277,0.72753;-3.48408,0.624088,0.533795;-4.40883,0.422247,-0.495334;-5.45345,1.32641,-0.688028;-5.57796,2.43691,0.146911;-6.58691,3.29722,-0.0515475;-4.65785,2.64294,1.17504;-4.80639,4.01148,2.20886;-3.61355,1.73849,1.36799;-1.07209,-2.21475,0.170776;0.0939354,-1.43492,0.0770612;1.30021,-2.03981,-0.242431;2.43295,-1.27322,-0.333386;2.73336,-0.721848,-1.60808;3.45719,0.603454,-1.41299;2.59548,1.60491,-0.650515;1.34882,1.8536,-1.37136;0.470419,2.74595,-0.602204;-0.818024,2.99985,-1.37951;-0.524608,3.5688,-2.65607;0.30286,2.69369,-3.42367;1.61814,2.43113,-2.69569)|"_smiles;
  REQUIRE(aligned_iressa);
  checkMolsHaveRoughlySameCoords(*iressa, *aligned_iressa);
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
  // This is the overlay produced by the first test below, to make sure we
  // haven't broken anything.
  auto ov_pdb_0zn_1tmn =
      R"(CC(C)C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O |(38.4182,43.8068,-0.910588;38.2709,42.2972,-0.731075;38.1304,41.9045,0.733229;37.0364,41.7932,-1.45451;37.1781,42.1406,-2.93763;35.8859,41.9023,-3.62687;35.9519,42.0509,-5.13225;34.6741,41.7191,-5.92414;34.7628,40.2524,-6.36478;33.5811,39.6611,-7.16439;33.5076,39.8293,-8.54156;32.4798,39.2259,-9.29221;31.5187,38.4455,-8.64246;31.6293,38.2599,-7.2498;32.6704,38.8605,-6.52059;36.3695,43.3994,-5.66929;37.2254,43.4931,-6.57465;35.7737,44.4114,-5.28982;38.3937,41.376,-3.47909;39.5342,41.9164,-3.54877;38.1091,40.1191,-3.80052;39.1087,39.1759,-4.32928;39.9624,38.5242,-3.20191;39.1024,37.5196,-2.43949;37.8503,37.6527,-2.00683;37.3343,36.383,-1.67591;38.2646,35.5274,-1.83827;38.2218,34.1585,-1.6274;39.4752,33.5046,-1.63354;40.6899,34.1995,-1.85784;40.7011,35.5753,-2.15222;39.4587,36.2072,-2.14728;38.5746,38.1494,-5.32118;37.3737,38.2894,-5.71533;39.3977,37.4441,-6.00821),wD:6.6,wU:4.3,21.22|)"_smiles;
  auto initScores = GaussianShape::ScoreMolecule(*pdb_trp_3tmn, *pdb_0zn_1tmn);
  CHECK_THAT(initScores[0], Catch::Matchers::WithinAbs(0.307, 0.001));
  CHECK_THAT(initScores[1], Catch::Matchers::WithinAbs(0.349, 0.001));
  CHECK_THAT(initScores[2], Catch::Matchers::WithinAbs(0.265, 0.001));
  // The PDB atom order isn't canonical, so bounce in and out of SMILES
  // to make it easier to check.
  auto canon_probe =
      v2::SmilesParse::MolFromSmiles(MolToCXSmiles(*pdb_0zn_1tmn));
  {
    // This should just tweak the input overlay.
    GaussianShape::ShapeOverlayOptions opts;
    opts.startMode = GaussianShape::StartMode::ROTATE_0;
    opts.normalize = false;
    GaussianShape::ShapeInputOptions shapeOpts;
    ROMol cp(*canon_probe);
    RDGeom::Transform3D xform;
    auto scores = GaussianShape::AlignMolecule(*pdb_trp_3tmn, cp, shapeOpts,
                                               shapeOpts, &xform, opts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.322, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.396, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.247, 0.001));
    CHECK(checkMolsHaveRoughlySameCoords(cp, *ov_pdb_0zn_1tmn));
  }
  {
    // With default settings, it does a poor job.
    GaussianShape::ShapeOverlayOptions opts;
    GaussianShape::ShapeInputOptions shapeOpts;
    ROMol cp(*canon_probe);
    auto scores = GaussianShape::AlignMolecule(*pdb_trp_3tmn, cp, shapeOpts,
                                               shapeOpts, nullptr, opts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.199, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.365, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.033, 0.001));
  }
  {
    // And with reference shape the same
    GaussianShape::ShapeOverlayOptions opts;
    opts.startMode = GaussianShape::StartMode::ROTATE_0;
    opts.normalize = false;
    GaussianShape::ShapeInputOptions shapeOpts;
    auto refShape = GaussianShape::ShapeInput(*pdb_trp_3tmn, -1, shapeOpts);
    ROMol cp(*canon_probe);
    RDGeom::Transform3D xform;
    auto scores =
        GaussianShape::AlignMolecule(refShape, cp, shapeOpts, &xform, opts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.322, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.396, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.247, 0.001));
    CHECK(checkMolsHaveRoughlySameCoords(cp, *ov_pdb_0zn_1tmn));
    MolTransforms::transformConformer(cp.getConformer(), xform);
    ROMol cp1(*canon_probe);
    MolTransforms::transformConformer(cp1.getConformer(), xform);
    CHECK(checkMolsHaveRoughlySameCoords(cp1, *ov_pdb_0zn_1tmn));
  }
}

TEST_CASE("Fragment Mode") {
  // On the PDB overlay.
  auto pdb_trp_3tmn =
      R"(N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2 |(37.935,40.394,-3.825;39.119,39.593,-4.13;38.758,38.486,-5.101;37.526,38.337,-5.395;39.716,37.852,-5.605;39.883,39.108,-2.906;39.086,38.098,-2.209;38.093,38.363,-1.34;37.565,37.179,-0.881;38.201,36.136,-1.471;39.193,36.684,-2.308;40.015,35.812,-3.036;39.846,34.441,-2.913;38.844,33.933,-2.075;38.015,34.752,-1.333),wU:1.0|)"_smiles;
  REQUIRE(pdb_trp_3tmn);
  auto pdb_0zn_1tmn =
      R"([C@H](CCc1ccccc1)(C(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2)CC(C)C |(35.672,41.482,-5.722;34.516,40.842,-6.512;34.843,39.355,-6.7;33.819,38.475,-7.45;33.825,38.414,-8.838;32.951,37.553,-9.53;32.064,36.747,-8.81;32.096,36.799,-7.402;32.985,37.656,-6.73;35.934,42.778,-6.452;36.833,42.858,-7.316;35.175,43.735,-6.275;35.516,41.561,-4.218;36.707,42.096,-3.513;38.055,41.449,-3.859;39.11,42.138,-3.959;37.975,40.129,-3.983;39.134,39.277,-4.298;38.825,38.04,-5.133;37.649,37.934,-5.605;39.788,37.369,-5.652;39.985,38.945,-3.037;39.221,37.953,-2.164;37.934,37.961,-1.823;37.579,36.695,-1.314;38.63,35.975,-1.286;39.736,36.771,-1.642;41.052,36.341,-1.48;41.213,35.042,-0.964;40.095,34.215,-0.69;38.765,34.665,-0.855;36.506,41.966,-2.002;37.6,42.757,-1.31;37.546,44.225,-1.728;37.408,42.58,0.19),wD:0.0,wU:17.21,13.33|)"_smiles;
  REQUIRE(pdb_0zn_1tmn);
  GaussianShape::ShapeOverlayOptions opts;
  opts.nSteps = 100;
  opts.startMode = GaussianShape::StartMode::ROTATE_180_FRAGMENT;
  opts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
  auto probeShape = GaussianShape::ShapeInput(*pdb_trp_3tmn, -1);
  auto refShape = GaussianShape::ShapeInput(*pdb_0zn_1tmn, -1);
  RDGeom::Transform3D xform;
  // Use the smaller molecule as the probe
  auto scores = GaussianShape::AlignShape(refShape, probeShape, &xform, opts);
  // These are close to the values above for starting from the xtal structures.
  CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.285, 0.001));
  CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.406, 0.001));
  CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.164, 0.001));
}

TEST_CASE("custom feature points") {
  auto m1 =
      "O=CC=O |(-1.75978,0.148897,0;-0.621382,-0.394324,0;0.624061,0.3656,.1;1.7571,-0.120174,.1)|"_smiles;
  SECTION("using shapes") {
    auto shape1 = GaussianShape::ShapeInput(*m1, -1);
    // each carbonyl O gets one feature:
    CHECK(shape1.getCoords().size() == 24);
    GaussianShape::ShapeInputOptions opts2;
    opts2.customFeatures = GaussianShape::CustomFeatures{
        {1, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0},
        {2, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0}};
    auto shape2 = GaussianShape::ShapeInput(*m1, -1, opts2);
    CHECK(shape2.getCoords().size() == 24);

    {
      // confirm that we don't add the features if not requested.
      GaussianShape::ShapeInputOptions topts;
      topts.customFeatures = GaussianShape::CustomFeatures{
          {1, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0},
          {2, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0}};
      topts.useColors = false;
      auto tshape = GaussianShape::ShapeInput(*m1, -1, topts);
      CHECK(tshape.getCoords().size() == 16);
    }

    // we'll swap the features on the second shape so that the alignment has to
    // be inverted
    GaussianShape::ShapeInputOptions opts3;
    opts3.customFeatures = GaussianShape::CustomFeatures{
        {2, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0},
        {1, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0}};

    auto m2 = ROMol(*m1);
    auto shape3 = GaussianShape::ShapeInput(m2, -1, opts3);
    CHECK(shape3.getCoords().size() == 24);
    GaussianShape::ShapeOverlayOptions overlayOpts;
    overlayOpts.optParam = 0.5;
    RDGeom::Transform3D xform;
    auto scores = AlignShape(shape2, shape3, &xform, overlayOpts);

    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(1.000, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.000, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.999, 0.001));
    CHECK(shape3.getCoords()[0] > 0);      // x coord of first atom
    CHECK(shape3.getCoords()[3 * 4] < 0);  // x coord of fourth atom

    auto conf = m2.getConformer(-1);
    MolTransforms::transformConformer(conf, xform);
    CHECK(conf.getAtomPos(0).x > 0);
    CHECK(conf.getAtomPos(3).x < 0);
  }
  SECTION("using molecules") {
    GaussianShape::ShapeInputOptions opts2;
    opts2.customFeatures = GaussianShape::CustomFeatures{
        {1, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0},
        {2, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0}};

    auto m2 = ROMol(*m1);
    // we'll swap the features on the second shape so that the alignment has to
    // be inverted
    GaussianShape::ShapeInputOptions opts3;
    opts3.customFeatures = GaussianShape::CustomFeatures{
        {2, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0},
        {1, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0}};

    GaussianShape::ShapeOverlayOptions overlayOpts;
    overlayOpts.optParam = 0.5;
    std::vector<float> matrix(12, 0.0);
    auto scores = AlignMolecule(*m1, m2, opts2, opts3, nullptr, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(1.000, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.000, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.999, 0.001));
    auto conf = m2.getConformer(-1);
    CHECK(conf.getAtomPos(0).x > 0);
    CHECK(conf.getAtomPos(3).x < 0);
  }
}

TEST_CASE("Non-standard radii") {
  auto m1 =
      "[Xe]c1ccccc1 |(0.392086,-2.22477,0.190651;0.232269,-1.38667,0.118385;-1.06274,-0.918982,0.0342466;-1.26098,0.446053,-0.0811879;-0.244035,1.36265,-0.11691;1.05134,0.875929,-0.031248;1.28797,-0.499563,0.0864097),atomProp:0.dummyLabel.*|"_smiles;
  GaussianShape::ShapeInputOptions shapeOpts;
  shapeOpts.useColors = false;
  shapeOpts.allCarbonRadii = false;
  auto shape1 = GaussianShape::ShapeInput(*m1, -1, shapeOpts);

  CHECK(shape1.getCoords().size() == 28);
  CHECK_THAT(shape1.getShapeVolume(),
             Catch::Matchers::WithinAbs(387.396, 0.005));
  // mol1 with atom 4 with an N radius and a bigger Xe.
  shapeOpts.atomRadii =
      std::vector<std::pair<unsigned int, double>>{{0, 2.5}, {4, 1.55}};
  shapeOpts.allCarbonRadii = false;
  auto shape2 = GaussianShape::ShapeInput(*m1, -1, shapeOpts);
  CHECK_THAT(shape2.getShapeVolume(),
             Catch::Matchers::WithinAbs(425.051, 0.005));

  // Corresponding pyridine derivative.
  auto m2 =
      "[Xe]c1ccncc1 |(0.392086,-2.22477,0.190651;0.232269,-1.38667,0.118385;-1.06274,-0.918982,0.0342466;-1.26098,0.446053,-0.0811879;-0.244035,1.36265,-0.11691;1.05134,0.875929,-0.031248;1.28797,-0.499563,0.0864097),atomProp:0.dummyLabel.*|"_smiles;
  auto shape3 = GaussianShape::ShapeInput(*m2, -1, shapeOpts);
  CHECK(shape3.getShapeVolume() == shape2.getShapeVolume());
}

TEST_CASE("Shape subset") {
  auto m1 =
      "c1ccc(-c2ccccc2)cc1 |(-3.26053,-0.0841607,-0.741909;-2.93383,0.123873,0.593407;-1.60713,0.377277,0.917966;-0.644758,0.654885,-0.0378428;0.743308,0.219134,0.168663;1.82376,1.0395,-0.0112769;3.01462,0.695405,0.613858;3.18783,-0.589771,1.09649;2.15761,-1.50458,1.01949;0.988307,-1.1313,0.385783;-1.1048,0.797771,-1.34022;-2.39754,0.435801,-1.69921)|"_smiles;
  REQUIRE(m1);
  GaussianShape::ShapeInputOptions shapeOpts;
  shapeOpts.atomSubset = std::vector<unsigned int>{0, 1, 2, 3, 10, 11};
  auto partShape = GaussianShape::ShapeInput(*m1, -1, shapeOpts);
  CHECK(partShape.getCoords().size() == 28);
  CHECK_THAT(partShape.getShapeVolume(),
             Catch::Matchers::WithinAbs(229.437, 0.005));
  CHECK_THAT(partShape.getColorVolume(),
             Catch::Matchers::WithinAbs(5.316, 0.005));

  shapeOpts.atomSubset.clear();
  auto wholeShape = GaussianShape::ShapeInput(*m1, -1, shapeOpts);
  CHECK(wholeShape.getCoords().size() == 56);
  CHECK_THAT(wholeShape.getShapeVolume(),
             Catch::Matchers::WithinAbs(556.266, 0.005));
  CHECK_THAT(wholeShape.getColorVolume(),
             Catch::Matchers::WithinAbs(10.631, 0.005));
}

#if 1
// Using this for benchmarking at the moment.  It should probably come out
// later.
TEST_CASE("LOBSTER") {
  std::string lobster_file =
      "/Users/david/Projects/Lobster/LOBSTER_112024/all_ligands.sdf";
  auto suppl = SDMolSupplier(lobster_file);
  std::vector<std::shared_ptr<ROMol>> mols;
  while (!suppl.atEnd()) {
    auto mol = suppl.next();
    mols.emplace_back(mol);
  }
  std::cout << "Number of mols " << mols.size() << std::endl;
  double sum_st = 0.0, sum_ct = 0.0, sum_comb = 0.0;
  int num = 0;
  std::mt19937 e2(1);
  std::uniform_real_distribution<double> unif(0, 1);
  GaussianShape::ShapeInputOptions opts;
  GaussianShape::ShapeOverlayOptions overlayOpts;
  overlayOpts.startMode = GaussianShape::StartMode::A_LA_PUBCHEM;
  for (const auto acr : std::vector<bool>{true}) {
    opts.allCarbonRadii = acr;
    for (size_t i = 1; i < mols.size(); i++) {
      for (size_t j = 0; j < i; j++) {
        if (unif(e2) > 0.001) {
          continue;
        }
        auto scores = GaussianShape::AlignMolecule(*mols[i], *mols[j], opts,
                                                   opts, nullptr, overlayOpts);
        sum_comb += scores[0];
        sum_st += scores[1];
        sum_ct += scores[2];
        ++num;
        if (!(num % 1000)) {
          std::cout << num << "  " << i << "  " << j << std::endl;
        }
      }
    }
    std::cout << "Mean combo of " << num << " : " << sum_comb / num
              << std::endl;
    std::cout << "Mean st of " << num << " : " << sum_st / num << std::endl;
    std::cout << "Mean ct of " << num << " : " << sum_ct / num << std::endl;
    if (acr) {
      CHECK_THAT(sum_st / num, Catch::Matchers::WithinAbs(0.554, 0.005));
      CHECK_THAT(sum_ct / num, Catch::Matchers::WithinAbs(0.094, 0.005));
    } else {
      CHECK_THAT(sum_st / num, Catch::Matchers::WithinAbs(0.540, 0.005));
      CHECK_THAT(sum_ct / num, Catch::Matchers::WithinAbs(0.098, 0.005));
    }
  }
}
#endif