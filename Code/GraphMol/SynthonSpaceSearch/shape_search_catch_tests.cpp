//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <cstdio>

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::SynthonSpaceSearch;
using namespace RDKit::RascalMCES;

const char *rdbase = getenv("RDBASE");

void prepareMolecule(RWMol *mol) {
  MolOps::addHs(*mol);
  auto dgParams = DGeomHelpers::ETKDGv3;
  // dgParams.pruneRmsThresh = 1.0;
  dgParams.randomSeed = 1;
  DGeomHelpers::EmbedMultipleConfs(*mol, 100, dgParams);
  MolOps::removeHs(*mol);
}

std::map<std::string, std::unique_ptr<ROMol>> loadLibrary(
    const std::string inFilename) {
  v2::FileParsers::SmilesMolSupplierParams params;
  params.titleLine = false;
  v2::FileParsers::SmilesMolSupplier suppl(inFilename, params);
  std::map<std::string, std::unique_ptr<ROMol>> mols;
  while (!suppl.atEnd()) {
    auto mol = suppl.next();
    if (mol) {
      prepareMolecule(mol.get());
      std::string molName = mol->getProp<std::string>(common_properties::_Name);
      mols.insert(std::make_pair(
          molName,
          std::unique_ptr<ROMol>(static_cast<ROMol *>(mol.release()))));
    }
  }
  return mols;
};

TEST_CASE("Shape Small tests") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
#if 1
  // These are the source files for the shape databases.  Useful to keep
  // around in case the databases ever need updating.
  std::vector<std::string> libNames{
      fullRoot + "amide_space.txt",
      fullRoot + "triazole_space.txt",
      fullRoot + "urea_space.txt",
  };
  std::vector<std::string> enumLibNames{
      fullRoot + "amide_space_enum.smi",
      fullRoot + "triazole_space_enum.smi",
      fullRoot + "urea_space_enum.smi",
  };
#endif
  std::vector<std::string> dbNames{
      fullRoot + "amide_space_shapes.spc",
      fullRoot + "triazole_space_shapes.spc",
      fullRoot + "urea_space_shapes.spc",
  };

  // The search of the enumerated libraries give 4, 8, 4 hits
  // respectively.  The first of these queries will have
  // conformers generated, the second and third will just use
  // the one given.
  std::vector<std::string> querySmis{
      "O=C(c1ccccc1)N1CCCC1 |(0.0443291,-1.81486,-1.76886;0.0506321,-0.858174,-0.921491;1.37975,-0.430412,-0.483603;2.18964,-1.35506,0.144714;3.47088,-1.00454,0.585539;3.93803,0.297573,0.388032;3.1267,1.22739,-0.242406;1.85597,0.849751,-0.670531;-1.14837,-0.261434,-0.446583;-1.26073,0.836916,0.520219;-2.73583,1.04666,0.696614;-3.34033,-0.283345,0.290893;-2.46516,-0.679843,-0.874401)|",
      "C[C@H]1CCN(c2nnc(CO)n2C2CCCC2)C1 |(3.88187,-1.9608,1.02401;3.40947,-0.556473,0.685633;3.49763,-0.278493,-0.787477;2.18424,0.313597,-1.21035;1.39217,0.480256,0.0110759;0.36454,1.42337,0.278193;0.656593,2.66561,0.766052;-0.477075,3.33073,0.923413;-1.48168,2.52297,0.540741;-2.93641,2.8664,0.551747;-3.33354,3.43671,-0.657732;-0.965924,1.32594,0.134935;-1.71735,0.224133,-0.333621;-1.11549,-0.643316,-1.37025;-2.26431,-1.65079,-1.54776;-2.58926,-1.96617,-0.0975393;-2.00229,-0.836476,0.740437;1.90383,-0.551243,0.906815),wD:1.0|",
      "C[C@@H]1C[C@H](NC(=O)NC2COC2)CN(C(=O)c2nccnc2F)C1 |(-0.346914,-0.986206,-4.28744;-0.686863,-0.0357247,-3.13265;0.429505,-0.1946,-2.14134;0.21099,0.659676,-0.907145;1.06526,0.0812473,0.104663;2.29297,0.75201,0.373712;2.5837,1.80373,-0.246936;3.23325,0.27478,1.33777;4.47647,0.99197,1.57602;4.94347,1.01294,2.99117;5.59284,-0.21541,2.82613;5.71049,0.107583,1.47157;-1.19052,0.721766,-0.417623;-2.14086,-0.0964663,-1.12904;-3.25312,-0.745252,-0.540367;-3.95877,-1.43507,-1.38825;-3.7227,-0.763533,0.803759;-4.9481,-0.204581,1.05395;-5.52492,-0.18107,2.24654;-4.86554,-0.748644,3.30451;-3.63585,-1.32731,3.13734;-3.08234,-1.32608,1.89052;-1.84839,-1.9059,1.69441;-2.02702,-0.329978,-2.57998),wD:1.0,wU:3.3|"};
  // The synthon search gives 1 hit for the urea space, where the
  // brute-force search gives 4 because the fragment similarities fall
  // below the threshold.  For example, comparing [2*]c1nccnc1F from
  // the query with synthon N#CCc(cncc1)c1[2*] (689988332-107515102)
  // when the dummy atoms are aligned, which they should be for a
  // good synthon match, the feature score is low because the nitrogen
  // acceptors don't align.  In the full molecule overlay, that is
  // compensated for by other things.
  std::vector<size_t> expNumHits{3, 8, 1};
  std::vector<std::vector<double>> expScores{
      {0.962, 0.951, 0.916},
      {0.955, 0.934, 0.923, 0.906, 0.862, 0.854, 0.832, 0.816},
      {0.614}};
  std::vector<std::vector<std::string>> expNames{
      {"1-1;2-1;amide-1", "1-2;2-1;amide-1", "1-3;2-1;amide-1"},
      {"1-1;2-1;3-1;triazole-1", "1-1;2-1;3-2;triazole-1",
       "1-1;2-2;3-1;triazole-1", "1-1;2-2;3-2;triazole-1",
       "1-2;2-1;3-1;triazole-1", "1-2;2-1;3-2;triazole-1",
       "1-2;2-2;3-1;triazole-1", "1-2;2-2;3-2;triazole-1"},
      {"277310376-742385846;182115391-684092275;487354835-896308859;urea-3"}};
  ShapeBuildParams shapeBuildOptions;
  shapeBuildOptions.numConfs = 100;
  shapeBuildOptions.rmsThreshold = 0.5;
  shapeBuildOptions.numThreads = 1;
  shapeBuildOptions.shapeSimThreshold = 0.95;
  shapeBuildOptions.randomSeed = 0xdac;

  for (size_t i = 0; i < dbNames.size(); i++) {
    // if (i != 0) {
    // continue;
    // }
    std::cout << i << " : " << dbNames[i] << std::endl;
    SynthonSpace synthonspace;
#if 1
    // In case the databases ever need updating.
    bool cancelled = false;
    synthonspace.readTextFile(libNames[i], cancelled);
    synthonspace.buildSynthonShapes(cancelled, shapeBuildOptions);
    synthonspace.writeDBFile(dbNames[i]);
#endif
    synthonspace.readDBFile(dbNames[i]);
    SynthonSpaceSearchParams params;
    params.similarityCutoff = 0.8;
    params.fragSimilarityAdjuster = 0.2;
    params.approxSimilarityAdjuster = 0.2;
    params.numConformers = shapeBuildOptions.numConfs;
    params.numThreads = shapeBuildOptions.numThreads;
    params.confRMSThreshold = shapeBuildOptions.rmsThreshold;
    params.timeOut = 0;
    params.randomSeed = shapeBuildOptions.randomSeed;
    params.bestHit = true;
    auto queryMol = v2::SmilesParse::MolFromSmiles(querySmis[i]);
    auto results = synthonspace.shapeSearch(*queryMol, params);
    unsigned int j = 0;
    CHECK(expNumHits[i] == results.getHitMolecules().size());
    std::cout << "Number of hits : " << results.getHitMolecules().size()
              << std::endl;
    for (const auto &mol : results.getHitMolecules()) {
      std::cout << "Hit : "
                << mol->getProp<std::string>(common_properties::_Name) << " : "
                << mol->getProp<double>("Similarity") << " : "
                << MolToCXSmiles(*mol) << std::endl;
      auto scores = GaussianShape::ScoreMolecule(*queryMol, *mol);
      std::cout << "check scores : " << scores[0] << ", " << scores[1] << ", "
                << scores[2] << std::endl;
    }
    for (const auto &mol : results.getHitMolecules()) {
      CHECK(mol->getProp<std::string>(common_properties::_Name) ==
            expNames[i][j]);
      CHECK_THAT(mol->getProp<double>("Similarity"),
                 Catch::Matchers::WithinAbs(expScores[i][j++], 0.001));
      auto scores = GaussianShape::ScoreMolecule(*queryMol, *mol);
      CHECK_THAT(mol->getProp<double>("Similarity"),
                 Catch::Matchers::WithinAbs(scores[0], 0.001));
    }
#if 0
    // Leave this in for now, in case we need to check brute force search
    // in future.
    auto mols = loadLibrary(enumLibNames[i]);
    prepareMolecule(queryMol.get());
    // RDKit::SDWriter sdw2(enumOutputNames[i]);
    std::vector<float> matrix(12, 0.0);
    unsigned int numHits = 0;
    for (auto &[smiles, mol] : mols) {
      bool foundHit = false;
      for (unsigned int i = 0; i < queryMol->getNumConformers(); ++i) {
        for (unsigned int j = 0; j < mol->getNumConformers(); ++j) {
          auto [st, ct] = AlignMolecule(*queryMol, *mol, matrix, i, j);
          if (st + ct > params.similarityCutoff) {
            std::cout << mol->getProp<std::string>(common_properties::_Name)
                      << " hit at " << st << ", " << ct << " : " << st + ct
                      << " for " << i << ", " << j << std::endl;
            ++numHits;
            foundHit = true;
            // sdw2.write(*mol);
            break;
          }
        }
        if (foundHit) {
          break;
        }
      }
    }
#endif
  }
}

TEST_CASE("Shape DB Writer") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/doebner_miller_space.txt";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);
  CHECK(synthonspace.getNumReactions() == 1);
  ShapeBuildParams shapeBuildParams;
  synthonspace.buildSynthonShapes(cancelled, shapeBuildParams);

  auto spaceName = std::tmpnam(nullptr);
  synthonspace.writeDBFile(spaceName);

  SynthonSpace newsynthonspace;
  newsynthonspace.readDBFile(spaceName);
  CHECK(newsynthonspace.getNumReactions() == 1);
  std::shared_ptr<SynthonSet> irxn;
  CHECK_NOTHROW(irxn = newsynthonspace.getReaction("doebner-miller-quinoline"));

  const auto &orxn = synthonspace.getReaction("doebner-miller-quinoline");
  for (size_t i = 0; i < irxn->getSynthons().size(); ++i) {
    REQUIRE(irxn->getSynthons()[i].size() == orxn->getSynthons()[i].size());
    for (size_t j = 0; j < irxn->getSynthons()[i].size(); ++j) {
      REQUIRE(irxn->getSynthons()[i][j]
                  .second->getShapes()
                  ->getShapes()
                  .getNumShapes() == orxn->getSynthons()[i][j]
                                         .second->getShapes()
                                         ->getShapes()
                                         .getNumShapes());
      for (size_t k = 0; k < irxn->getSynthons()[i][j]
                                 .second->getShapes()
                                 ->getShapes()
                                 .getNumShapes();
           ++k) {
        const auto ishape = irxn->getSynthons()[i][j].second->getShapes().get();
        const auto oshape = orxn->getSynthons()[i][j].second->getShapes().get();
        CHECK_THAT(fabs(ishape->getShapes().getShapeVolume(k) -
                        oshape->getShapes().getShapeVolume(k)),
                   Catch::Matchers::WithinAbs(0.0, 1.0e-6));
      }
    }
  }
}

TEST_CASE("Unspecified Stereo") {
  auto m1 = "C[C@H](Cl)CCOOC(Cl)F"_smiles;
  REQUIRE(m1);
  CHECK(details::hasUnspecifiedStereo(*m1) == true);
  CHECK(details::countChiralAtoms(*m1) == 2);
  auto m2 = "C[C@H](Cl)CCOO[C@@H](Cl)F"_smiles;
  REQUIRE(m2);
  CHECK(details::hasUnspecifiedStereo(*m2) == false);
  CHECK(details::countChiralAtoms(*m2) == 2);

  auto m3 = "C[C@H](Cl)CCOO[C@@](Cl)(F)CC=CC"_smiles;
  REQUIRE(m3);
  CHECK(details::hasUnspecifiedStereo(*m3) == true);

  auto m4 = R"(C[C@H](Cl)CCOO[C@@](Cl)(F)C\C=C/C)"_smiles;
  REQUIRE(m4);
  CHECK(details::hasUnspecifiedStereo(*m4) == false);

  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string spaceName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/small_freedom_shapes.spc";
  SynthonSpace space;
  space.readDBFile(spaceName);

  ShapeBuildParams shapeOptions;
  shapeOptions.randomSeed = 1;
  bool cancelled = false;
  space.buildSynthonShapes(cancelled, shapeOptions);

  SynthonSpaceSearchParams params;
  params.similarityCutoff = 1.6;
  params.enumerateUnspecifiedStereo = false;

  // This should bale with no results because there's unspecified
  // stereochem.
  auto results = space.shapeSearch(*m1, params);
  CHECK(results.getHitMolecules().empty());
}

TEST_CASE("Bad elements") {
  // This bit of FreedomSpace 2024-09 threw an exception originally
  // due to the Pd atom not being recognised by the shape builder.
  std::istringstream iss(R"(SMILES	synton_id	synton#	reaction_id
Cc1sc2ncnc(NN[U])c2c1C	100000182719	1	20a	2024-09
ClC(Cl)=C(Cl)c1cccc(N[U])c1	100016989384	1	20a	2024-09
Cl[Pd]c1cccc(-c2ccccc2)c1N[U]	114813040699	1	20a	2024-09
C=CCC1(S(=O)(=O)[U])CC1	100000101377	2	20a	2024-09
C=Cc1ccc(S(=O)(=O)[U])cc1	100000034458	2	20a	2024-09
)");
  SynthonSpace space;
  bool cancelled = false;
  space.readStream(iss, cancelled);
  ShapeBuildParams shapeOptions;
  shapeOptions.randomSeed = 1;
  CHECK_NOTHROW(space.buildSynthonShapes(cancelled, shapeOptions));
}

TEST_CASE("Shape Best Hit Found") {
  // Make sure that even if the threshold isn't met, there's a 3D hit returned
  // of the best that was possible.
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string spaceName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/small_freedom_shapes.spc";

  SynthonSpace space;
  space.readDBFile(spaceName);
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  params.numThreads = -1;
  params.numConformers = 200;
  params.confRMSThreshold = 0.25;
  params.randomSeed = 0xdac;

  params.fragSimilarityAdjuster = 0.1;
  params.approxSimilarityAdjuster = 0.1;
  params.similarityCutoff = 1.0;
  params.timeOut = 0;

  auto mol =
      "CCC(C(=O)NCc1ccco1)N(Cc1sccc1C)C(C)C |(1.19967,-2.26511,-1.8853;-0.0674677,-1.53728,-1.54329;-0.0395195,-0.735921,-0.2957;0.978688,0.316536,-0.356493;0.610079,1.54481,-0.391485;2.37979,0.114038,-0.377756;3.31574,1.24811,-0.443092;4.723,0.774447,-0.458657;5.58509,0.534718,0.592748;6.76773,0.108395,-0.00740365;6.56938,0.109237,-1.38929;5.34763,0.509833,-1.60401;-1.33892,-0.260032,0.0922388;-2.08764,0.476264,-0.832263;-3.39845,0.92709,-0.383151;-4.99177,0.223473,-0.828611;-5.94415,1.34626,0.133161;-5.00918,2.1798,0.729651;-3.7235,1.96002,0.462253;-2.60368,2.81164,1.05297;-1.99138,-1.01276,1.09008;-1.10607,-0.965532,2.34165;-2.26411,-2.45544,0.780694)|"_smiles;
  REQUIRE(mol);
  auto results = space.shapeSearch(*mol, params);
  CHECK(results.getHitMolecules().empty());
  auto &bestHit = results.getBestHit();
  CHECK(bestHit);
  CHECK_NOTHROW(bestHit->getProp<std::string>(common_properties::_Name));
  CHECK(bestHit->getProp<double>("Similarity") < 1.0);
  CHECK(bestHit->getConformer().is3D());
}

TEST_CASE("Small REAL test") {
  std::string dbFile =
      "/Users/david/Projects/SynthonSpaceTests/REAL/2024-09_RID-4-Cozchemix/random_real_0_confs.spc";
  SynthonSpace space;
  space.readDBFile(dbFile, 13);
  std::cout << space.getNumReactions() << " reactions for "
            << space.getNumSynthons() << " giving " << space.getNumProducts()
            << " products." << std::endl;
  SynthonSpaceSearchParams params;
  params.similarityCutoff = 0.75;
  params.numThreads = 13;
  params.timeOut = 0;
  params.randomSeed = 1;
  params.useProgressBar = 60;
  auto queryMol = v2::FileParsers::MolFromMolFile(
      "/Users/david/Projects/SynthonSpaceTests/FreedomSpace/esomeprazole.sdf");
  auto results = space.shapeSearch(*queryMol, params);
  std::cout << "Number of hits : " << results.getHitMolecules().size()
            << std::endl;
  for (const auto &mol : results.getHitMolecules()) {
    std::cout << MolToCXSmiles(*mol) << std::endl;
  }
  auto &bestHit = results.getBestHit();
  CHECK(bestHit);
  std::cout << "Best hit : " << MolToCXSmiles(*bestHit) << std::endl;
  std::cout << bestHit->getProp<std::string>(common_properties::_Name);
}

TEST_CASE("Small m_1458bbb test") {
  std::string dbFile =
      "/Users/david/Projects/SynthonSpaceTests/REAL/2024-09_RID-4-Cozchemix/m_1458bbb_confs.spc";
  SynthonSpace space;
  space.readDBFile(dbFile, 13);
  std::cout << space.getNumReactions() << " reactions for "
            << space.getNumSynthons() << " giving " << space.getNumProducts()
            << " products." << std::endl;
  SynthonSpaceSearchParams params;
  params.similarityCutoff = 0.75;
  params.numThreads = 13;
  params.timeOut = 0;
  params.randomSeed = 1;
  params.useProgressBar = 60;
  params.confRMSThreshold = 0.5;
  params.numConformers = 100;
  auto queryMol =
      "C[C@H](OC(=O)[C@H]1CC[C@@H](CN(C)C)O1)c1ccc(S(=O)(=O)F)cc1 |(2.56475,-3.10298,0.955074;2.28607,-1.62018,0.834903;0.929836,-1.35524,0.537459;0.117132,-0.687892,1.44187;0.639023,-0.335941,2.52603;-1.28757,-0.397477,1.16468;-2.1698,-1.60577,1.00963;-3.55218,-0.939503,1.06929;-3.20651,0.442004,1.62104;-4.12741,0.889789,2.71027;-5.47952,1.12966,2.32503;-6.24538,-0.0137629,1.94604;-6.13257,1.76958,3.47883;-1.92681,0.272412,2.203;3.21863,-0.957246,-0.102791;4.30426,-1.59526,-0.643779;5.12411,-0.893658,-1.51723;4.85003,0.417886,-1.8326;5.89923,1.29707,-2.9417;6.59479,0.266104,-3.80483;5.14827,2.28511,-3.74963;7.08009,2.05388,-2.03558;3.75958,1.07193,-1.29428;2.96242,0.354879,-0.431309),wD:1.0,wU:5.4,8.8|"_smiles;
  REQUIRE(queryMol);
  auto results = space.shapeSearch(*queryMol, params);
  std::cout << "Number of hits : " << results.getHitMolecules().size()
            << std::endl;
  for (const auto &mol : results.getHitMolecules()) {
    std::cout << MolToCXSmiles(*mol) << std::endl;
  }
  auto &bestHit = results.getBestHit();
  REQUIRE(bestHit);
  std::cout << "Best hit : " << MolToCXSmiles(*bestHit) << std::endl;
  std::cout << bestHit->getProp<std::string>(common_properties::_Name) << " : "
            << bestHit->getProp<double>("Similarity") << std::endl;
}

TEST_CASE("Two piece query") {
  {
    auto queryMol =
        "C1CCNC1.c1ccccc1 |(-2.73583,1.04666,0.696614;-3.34033,-0.283345,0.290893;-2.46516,-0.679843,-0.874401;-1.14837,-0.261434,-0.446583;-1.26073,0.836916,0.520219;3.47088,-1.00454,0.585539;2.18964,-1.35506,0.144714;1.37975,-0.430412,-0.483603;1.85597,0.849751,-0.670531;3.1267,1.22739,-0.242406;3.93803,0.297573,0.388032)|"_smiles;
    REQUIRE(queryMol);
    std::string fName(rdbase);
    std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
    std::string dbName = fullRoot + "amide_space_shapes.spc";
    SynthonSpace synthonspace;
    synthonspace.readDBFile(dbName);
    SynthonSpaceSearchParams params;
    params.similarityCutoff = 0.7;
    params.numThreads = 1;
    params.confRMSThreshold = 0.5;
    params.timeOut = 0;
    params.randomSeed = 1;
    params.bestHit = true;
    params.shapeOverlayOptions.simAlpha = 0.95;
    params.shapeOverlayOptions.simBeta = 0.05;
    auto results = synthonspace.shapeSearch(*queryMol, params);
    CHECK(results.getHitMolecules().size() == 3);
    std::vector<double> expScores{0.721, 0.714, 0.711};
    std::cout << "Number of results : " << results.getHitMolecules().size()
              << std::endl;
    for (unsigned int i = 0; i < results.getHitMolecules().size(); ++i) {
      auto &mol = results.getHitMolecules()[i];
      std::cout << "hit mol : " << MolToCXSmiles(*mol) << std::endl;
      std::cout << "hit name : "
                << mol->getProp<std::string>(common_properties::_Name) << " at "
                << mol->getProp<double>("Similarity") << std::endl;
      CHECK_THAT(mol->getProp<double>("Similarity"),
                 Catch::Matchers::WithinAbs(expScores[i], 0.001));
    }
  }

  {
    auto queryMol =
        "C1CCCC1.C[C@H]1CCNC1 |(-1.11549,-0.643316,-1.37025;-2.26431,-1.65079,-1.54776;-2.58926,-1.96617,-0.0975393;-2.00229,-0.836476,0.740437;-1.71735,0.224133,-0.333621;3.88187,-1.9608,1.02401;3.40947,-0.556473,0.685633;3.49763,-0.278493,-0.787477;2.18424,0.313597,-1.21035;1.39217,0.480256,0.0110759;1.90383,-0.551243,0.906815),wD:6.5|"_smiles;
    REQUIRE(queryMol);
    std::string fName(rdbase);
    std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
    std::string dbName = fullRoot + "triazole_space_shapes.spc";
    SynthonSpace synthonspace;
    synthonspace.readDBFile(dbName);
    SynthonSpaceSearchParams params;
    params.similarityCutoff = 0.7;
    params.fragSimilarityAdjuster = 0.2;
    params.numThreads = 1;
    params.confRMSThreshold = 0.5;
    params.timeOut = 0;
    params.randomSeed = 1;
    params.bestHit = true;
    params.shapeOverlayOptions.simAlpha = 0.95;
    params.shapeOverlayOptions.simBeta = 0.05;
    auto results = synthonspace.shapeSearch(*queryMol, params);
    std::cout << "Number of results : " << results.getHitMolecules().size()
              << std::endl;
    CHECK(results.getHitMolecules().size() == 6);
    std::vector<double> expScores{0.816, 0.790, 0.745, 0.735, 0.718, 0.709};
    for (unsigned int i = 0; i < results.getHitMolecules().size(); ++i) {
      auto &mol = results.getHitMolecules()[i];
      std::cout << "hit mol : " << MolToCXSmiles(*mol) << std::endl;
      std::cout << "hit name : "
                << mol->getProp<std::string>(common_properties::_Name) << " at "
                << mol->getProp<double>("Similarity") << std::endl;
      CHECK_THAT(mol->getProp<double>("Similarity"),
                 Catch::Matchers::WithinAbs(expScores[i], 0.001));
    }
  }
}

std::unique_ptr<RWMol> generateConfs(const std::string &smiles,
                                     unsigned int numConformers) {
  auto retMol = v2::SmilesParse::MolFromSmiles(smiles);
  MolOps::addHs(*retMol);
  auto dgParams = DGeomHelpers::ETKDGv3;
  dgParams.numThreads = 1;
  dgParams.pruneRmsThresh = -1;
  dgParams.randomSeed = 0xdac;
  auto cids =
      DGeomHelpers::EmbedMultipleConfs(*retMol, numConformers, dgParams);
  if (cids.empty() || !retMol->getNumConformers()) {
    retMol.reset();
  } else {
    MolOps::removeHs(*retMol);
  }
  return retMol;
}

TEST_CASE("User-defined conformer generator") {
  // This is the urea space from above.
  std::string spaceText(
      R"(SMILES	synton_id	synton#	reaction_id
O=C(NC1COC1)[1*]	277310376-742385846	0	urea-3
#Oc(cc(cc1)NC([1*])=O)c1Cl	287123986-010598048	0	urea-3
O=C(Nc1c(CN[1*])cc[s]1)[2*]	584456271-623025187	1	urea-3
C[C@H](CC(C1)N[1*])CN1C([2*])=O	182115391-684092275	1	urea-3
FC(Oc(cc1)cnc1[2*])F	441848376-976230122	2	urea-3
Fc1nccnc1[2*]	487354835-896308859	2	urea-3
N#CCc(cncc1)c1[2*]	689988332-107515102	2	urea-3)");
  std::istringstream iss(spaceText);
  ShapeBuildParams shapeBuildParamsWith;
  shapeBuildParamsWith.numThreads = 1;
  shapeBuildParamsWith.numConfs = 10;
  shapeBuildParamsWith.userConformerGenerator = generateConfs;
  shapeBuildParamsWith.shapeSimThreshold = -1.0;
  bool cancelled = false;
  SynthonSpace space;
  space.readStream(iss, cancelled);
  space.buildSynthonShapes(cancelled, shapeBuildParamsWith);
  // There should be 10 conformers of each synthon except 182115391-684092275
  // which has an unspecified stereocentre so gets 20.
  auto react = space.getReaction(("urea-3"));
  for (auto sst : react->getSynthons()) {
    for (auto &[sn, s] : sst) {
      if (sn == "182115391-684092275") {
        CHECK(s->getShapes()->getShapes().getNumShapes() == 20);
      } else {
        CHECK(s->getShapes()->getShapes().getNumShapes() == 10);
      }
    }
  }

  // This is the same query as in Shape Small tests.
  auto queryMol =
      "C[C@@H]1C[C@H](NC(=O)NC2COC2)CN(C(=O)c2nccnc2F)C1 |(-0.346914,-0.986206,-4.28744;-0.686863,-0.0357247,-3.13265;0.429505,-0.1946,-2.14134;0.21099,0.659676,-0.907145;1.06526,0.0812473,0.104663;2.29297,0.75201,0.373712;2.5837,1.80373,-0.246936;3.23325,0.27478,1.33777;4.47647,0.99197,1.57602;4.94347,1.01294,2.99117;5.59284,-0.21541,2.82613;5.71049,0.107583,1.47157;-1.19052,0.721766,-0.417623;-2.14086,-0.0964663,-1.12904;-3.25312,-0.745252,-0.540367;-3.95877,-1.43507,-1.38825;-3.7227,-0.763533,0.803759;-4.9481,-0.204581,1.05395;-5.52492,-0.18107,2.24654;-4.86554,-0.748644,3.30451;-3.63585,-1.32731,3.13734;-3.08234,-1.32608,1.89052;-1.84839,-1.9059,1.69441;-2.02702,-0.329978,-2.57998),wD:1.0,wU:3.3|"_smiles;
  REQUIRE(queryMol);
  SynthonSpaceSearchParams spaceSearchParams;
  spaceSearchParams.userConformerGenerator = generateConfs;
  spaceSearchParams.useProgressBar = 60;
  spaceSearchParams.timeOut = 0;
  spaceSearchParams.similarityCutoff = 0.8;
  spaceSearchParams.bestHit = true;
  auto results = space.shapeSearch(*queryMol, spaceSearchParams);
  CHECK(results.getHitMolecules().size() == 1);
  // This set of conformer parameters doesn't give a decent hit, but we
  // get something because of bestHit=true.
  for (const auto &mol : results.getHitMolecules()) {
    CHECK_THAT(mol->getProp<double>("Similarity"),
               Catch::Matchers::WithinAbs(0.728, 0.001));
    auto scores = GaussianShape::ScoreMolecule(*queryMol, *mol);
    CHECK_THAT(mol->getProp<double>("Similarity"),
               Catch::Matchers::WithinAbs(scores[0], 0.001));
  }
}

TEST_CASE("Another test") {
#if 1
  std::string spaceText(
      R"(SMILES	synton_id	synton#	reaction_id
O=C(NC1(CC2CC2)CCC1)[1*]	016642174-703335348	0	urea-3
CN(CC1(CC1)N[1*])C([2*])=O	009516132-856691534	1	urea-3
CSc1cc(C(O)=O)nc([2*])n1	821136904-635555130	2	urea-3
)");
  std::istringstream iss(spaceText);
  ShapeBuildParams shapeBuildParamsWith;
  shapeBuildParamsWith.numThreads = 1;
  shapeBuildParamsWith.numConfs = 10;
  shapeBuildParamsWith.userConformerGenerator = generateConfs;
  shapeBuildParamsWith.shapeSimThreshold = -1.0;
  bool cancelled = false;
  SynthonSpace space;
  space.readStream(iss, cancelled);
  space.buildSynthonShapes(cancelled, shapeBuildParamsWith);
#endif
  auto queryMol =
      "O=CNCc1ccsc1NC(=O)c1cnccn1 |(-3.13593,2.53671,-0.964221;-2.92978,1.52511,-1.69943;-2.53338,0.26403,-1.14888;-2.37819,0.196142,0.278959;-1.94716,-1.21639,0.708687;-2.82815,-2.05481,1.26358;-2.27563,-3.33898,1.61555;-0.605283,-3.29034,1.15883;-0.664978,-1.59895,0.520406;0.426145,-0.897222,-0.0427085;1.31101,-0.0672121,0.677347;1.15187,0.0882736,1.90787;2.40272,0.602377,-0.0218247;3.28333,1.42622,0.662823;4.30841,2.05393,0.0205924;4.47886,1.87718,-1.3136;3.60984,1.05805,-2.01765;2.62331,0.465647,-1.34378)|"_smiles;
  REQUIRE(queryMol);
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
  std::string dbName = fullRoot + "idorsia_toy_space_shapes.spc";
  SynthonSpaceSearchParams params;
  params.similarityCutoff = 0.6;
  params.numConformers = 100;
  params.numThreads = -1;
  params.confRMSThreshold = 1.0;
  params.timeOut = 0;
  params.randomSeed = 0xdac;
  params.bestHit = true;
  params.useProgressBar = 60;
  SynthonSpace synthonspace;
  synthonspace.readDBFile(dbName);
  auto results = synthonspace.shapeSearch(*queryMol, params);
  // auto results = space.shapeSearch(*queryMol, params);
  std::cout << "Num hits : " << results.getHitMolecules().size() << std::endl;
  for (const auto &mol : results.getHitMolecules()) {
    std::cout << mol->getProp<std::string>("_Name") << " : "
              << mol->getProp<double>("Similarity") << std::endl;
    auto scores = GaussianShape::ScoreMolecule(*queryMol, *mol);
    CHECK_THAT(mol->getProp<double>("Similarity"),
               Catch::Matchers::WithinAbs(scores[0], 0.001));
    std::cout << "check scores : " << scores[0] << ", " << scores[1] << ", "
              << scores[2] << std::endl;
  }
  SDWriter sdw2("another_test_hits.sdf");
  for (const auto &mol : results.getHitMolecules()) {
    sdw2.write(*mol);
  }
}

TEST_CASE("Unmatched synthon") {
  std::string spaceText(
      R"(SMILES	synton_id	synton#	reaction_id
c1ccccc1[1*]	1-1	0	test1
c1cccnc1[1*]	1-2	0	test1
[1*]C[2*]	2-1	1	test1
[1*]S[2*]	2-2	1	test1
[2*]C1CCCCC1 	3-1	2	test1
[2*]C1CCCC1	3-2	2	test1
)");
  std::istringstream iss(spaceText);
  ShapeBuildParams shapeBuildParams;
  shapeBuildParams.numThreads = 1;
  shapeBuildParams.numConfs = 10;
  shapeBuildParams.userConformerGenerator = generateConfs;
  shapeBuildParams.shapeSimThreshold = -1.0;
  shapeBuildParams.randomSeed = 0xdac;
  bool cancelled = false;
  SynthonSpace synthonSpace;
  synthonSpace.readStream(iss, cancelled);
  synthonSpace.buildSynthonShapes(cancelled, shapeBuildParams);
  auto queryMol =
      "c1ccc(C2CCCCC2)cc1 |(-2.862,-0.0009,0.0236;-2.1138,0.001,1.211;-0.7111,0,1.2027;0,0,;1.486,0,;2.0193,1.264,-0.6833;3.5453,1.2619,-0.7532;4.076,0,-1.4269;3.5453,-1.262,-0.7532;2.0193,-1.264,-0.6833;-0.7114,-0.0194,-1.1991;-2.1126,-0.0155,-1.1687)|"_smiles;
  SynthonSpaceSearchParams params;
  params.similarityCutoff = 0.7;
  params.fragSimilarityAdjuster = 0.3;
  params.approxSimilarityAdjuster = 0.3;
  params.numConformers = 100;
  params.numThreads = 1;
  params.confRMSThreshold = 1.0;
  params.timeOut = 0;
  params.randomSeed = 0xdac;
  params.bestHit = true;
  auto results = synthonSpace.shapeSearch(*queryMol, params);
  std::cout << "Num hits : " << results.getHitMolecules().size() << std::endl;
  CHECK(results.getHitMolecules().size() == 6);
  std::vector<std::string> expNames{
      "1-1;2-1;3-2;test1", "1-1;2-1;3-1;test1", "1-2;2-1;3-1;test1",
      "1-1;2-2;3-1;test1", "1-1;2-2;3-2;test1", "1-2;2-1;3-2;test1",
  };
  std::vector<double> expScores{0.835, 0.781, 0.752, 0.732, 0.730, 0.729};
  for (size_t i = 0; i < results.getHitMolecules().size(); ++i) {
    const auto &mol = results.getHitMolecules()[i];
    std::cout << mol->getProp<std::string>("_Name") << " : "
              << mol->getProp<double>("Similarity") << std::endl;
    std::cout << MolToCXSmiles(*mol) << std::endl;
  }
  for (size_t i = 0; i < expNames.size(); ++i) {
    const auto &mol = results.getHitMolecules()[i];
    CHECK(mol->getProp<std::string>("_Name") == expNames[i]);
    CHECK_THAT(mol->getProp<double>("Similarity"),
               Catch::Matchers::WithinAbs(expScores[i], 0.001));
  }
}

TEST_CASE("Trim sample molecules") {
  std::vector<std::tuple<std::string, unsigned int, std::string>> smiles{{R"(N#CC1=C(SCC2(O)CCN(C(=O)OCc3ccccc3)CC2)NC(N)=C(c2nc(-c3ccccc3)cs2)C12CCCCC2 |atomProp:0.molNum.0:0.idx.0:1.molNum.0:1.idx.1:2.molNum.0:2.idx.2:3.molNum.0:3.idx.3:4.molNum.0:4.idx.4:4.HoldsJoin.1:5.molNum.1:5.idx.15:5.HoldsJoin.1:6.molNum.1:6.idx.13:7.molNum.1:7.idx.14:8.molNum.1:8.idx.12:9.molNum.1:9.idx.11:10.molNum.1:10.idx.10:11.molNum.1:11.idx.1:12.molNum.1:12.idx.0:13.molNum.1:13.idx.2:14.molNum.1:14.idx.3:15.molNum.1:15.idx.4:16.molNum.1:16.idx.5:17.molNum.1:17.idx.6:18.molNum.1:18.idx.7:19.molNum.1:19.idx.8:20.molNum.1:20.idx.9:21.molNum.1:21.idx.18:22.molNum.1:22.idx.17:23.molNum.0:23.idx.6:24.molNum.0:24.idx.7:25.molNum.0:25.idx.8:26.molNum.0:26.idx.9:27.molNum.0:27.idx.10:28.molNum.0:28.idx.11:29.molNum.0:29.idx.12:30.molNum.0:30.idx.13:31.molNum.0:31.idx.14:32.molNum.0:32.idx.15:33.molNum.0:33.idx.16:34.molNum.0:34.idx.17:35.molNum.0:35.idx.18:36.molNum.0:36.idx.19:37.molNum.0:37.idx.20:38.molNum.0:38.idx.21:39.molNum.0:39.idx.22:40.molNum.0:40.idx.23:41.molNum.0:41.idx.24:42.molNum.0:42.idx.25:43.molNum.0:43.idx.26|)",
                                                                          1,
                                                                          R"(O=C(OCc1ccccc1)N1CCC(O)(CS)CC1)"},
                                                                         {R"(C=CCn1c(COc2ccccc2)nnc1S(=O)(=O)CC(=O)N(Cc1ccc(C)cc1)C1CC1 |atomProp:0.molNum.0:0.idx.0:1.molNum.0:1.idx.1:2.molNum.0:2.idx.2:3.molNum.0:3.idx.3:4.molNum.0:4.idx.4:5.molNum.0:5.idx.5:6.molNum.0:6.idx.6:7.molNum.0:7.idx.7:8.molNum.0:8.idx.8:9.molNum.0:9.idx.9:10.molNum.0:10.idx.10:11.molNum.0:11.idx.11:12.molNum.0:12.idx.12:13.molNum.0:13.idx.13:14.molNum.0:14.idx.14:15.molNum.0:15.idx.15:16.molNum.0:16.idx.16:16.HoldsJoin.1:17.molNum.0:17.idx.17:18.molNum.0:18.idx.18:19.molNum.1:19.idx.9:19.HoldsJoin.1:20.molNum.1:20.idx.7:21.molNum.1:21.idx.8:22.molNum.1:22.idx.6:23.molNum.1:23.idx.5:24.molNum.1:24.idx.4:25.molNum.1:25.idx.3:26.molNum.1:26.idx.2:27.molNum.1:27.idx.1:28.molNum.1:28.idx.0:29.molNum.1:29.idx.15:30.molNum.1:30.idx.14:31.molNum.1:31.idx.11:32.molNum.1:32.idx.12:33.molNum.1:33.idx.13|)", 1,
                                                                          "Cc1ccc(CN(C(=O)CS(=O)(=O)c2nnc[nH]2)C2CC2)cc1"},
                                                                         {R"(CC1CCN(c2nnc(CO)n2C2CCCC2)C1 |atomProp:0.molNum.2:0.idx.0:1.molNum.2:1.idx.1:2.molNum.2:2.idx.2:3.molNum.2:3.idx.3:4.molNum.2:4.idx.4:5.molNum.2:5.idx.6:5.HoldsJoin.2_3:6.molNum.0:6.idx.5:6.HoldsJoin.2:7.molNum.0:7.idx.4:8.molNum.0:8.idx.2:8.HoldsJoin.1:9.molNum.0:9.idx.1:10.molNum.0:10.idx.0:11.molNum.1:11.idx.5:11.HoldsJoin.3_1:12.molNum.1:12.idx.4:13.molNum.1:13.idx.0:14.molNum.1:14.idx.1:15.molNum.1:15.idx.2:16.molNum.1:16.idx.3:17.molNum.2:17.idx.5|)",
                                                                          2,
                                                                          R"(Cc1nnc(N2CCC(C)C2)n1C1CCCC1)"}};
  for (auto &[smi, molNum, expSmi] : smiles) {
    auto mol = v2::SmilesParse::MolFromSmiles(smi);
    REQUIRE(mol);
    auto newMol = SynthonSpaceSearch::details::trimSampleMol(*mol, molNum);
    CHECK(MolToSmiles(*newMol) == expSmi);
  }
}

TEST_CASE("Zero-length vector") {
  std::string dbName(
      "/home/dave/Projects/SynthonSpaceTests/REAL/random_real_0_conforge.spc");
  SynthonSpace synthonspace;
  synthonspace.readDBFile(dbName);

  auto queryMol =
      "CC(OC(=O)[C@H]1CC[C@@H](CN(C)C)O1)c1ccc(S(=O)(=O)F)cc1 |(2.90369,-2.16344,2.44104;2.07287,-1.7564,1.24069;1.06154,-0.810496,1.59795;-0.274374,-1.06478,1.34641;-0.643259,-2.13112,0.797722;-1.25553,-0.0302441,1.74923;-0.900332,1.25268,1.05924;-1.90239,1.25072,-0.0924938;-3.06709,0.550454,0.566467;-3.81933,-0.177119,-0.498406;-4.40849,0.666278,-1.51173;-5.08574,-0.25176,-2.4403;-5.41599,1.54263,-1.00978;-2.55682,-0.356408,1.46659;2.9964,-1.00911,0.319988;3.45629,-1.57333,-0.854797;4.29829,-0.834768,-1.65047;4.67403,0.452586,-1.27287;5.75692,1.33001,-2.35179;5.52648,2.79645,-2.25924;5.53034,0.850543,-3.74582;7.33965,0.982198,-1.87731;4.21515,1.02107,-0.096445;3.35774,0.270642,0.71504),wD:5.4,8.8| bA2k2xlIb7clbAeeDnjr3Q;quHGZvNovRYxJzoZJbdcSQ;m_1458cgb"_smiles;
  SynthonSpaceSearchParams params;
  params.similarityCutoff = 0.7;
  params.fragSimilarityAdjuster = 0.2;
  params.approxSimilarityAdjuster = 0.2;
  params.numConformers = 100;
  params.numThreads = 1;
  params.confRMSThreshold = 1.0;
  params.timeOut = 0;
  params.randomSeed = 0xdac;
  params.bestHit = true;

  auto hits = synthonspace.shapeSearch(*queryMol, params);
}

TEST_CASE("Bad Shape Mol") {
  std::string spaceText(
      R"(SMILES	synton_id	synton#	reaction_id	release
Fc1ccc(C2(C(N=[1*])=NO[2*])CC2)cc1	4RFh0fOYsX5pYFeRfm3hxw	1	m_265764cbb	3
Cn1nc(S(C)(=O)=O)cc1C(=[1*])[2*]	kxWX6dr7Fa59vLm4rVNo0g	2	m_265764cbb	3
)");
  std::istringstream iss(spaceText);
  ShapeBuildParams shapeBuildParams;
  shapeBuildParams.numThreads = 1;
  shapeBuildParams.numConfs = 10;
  shapeBuildParams.userConformerGenerator = generateConfs;
  shapeBuildParams.shapeSimThreshold = -1.0;
  shapeBuildParams.randomSeed = 0xdac;
  bool cancelled = false;
  SynthonSpace synthonSpace;
  synthonSpace.readStream(iss, cancelled);
  synthonSpace.buildSynthonShapes(cancelled, shapeBuildParams);
}