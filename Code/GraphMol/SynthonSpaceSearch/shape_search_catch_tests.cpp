//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <filesystem>

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
#if 0
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
      {0.962, 0.954, 0.916},
      {0.955, 0.949, 0.934, 0.903, 0.860, 0.854, 0.845, 0.807},
      {0.614},
      {0.962, 0.955, 0.934, 0.903, 0.855, 0.854, 0.851, 0.817},
      {0.643}};
  std::vector<std::vector<std::string>> expNames{
      {"1-1;2-1;amide-1", "1-2;2-1;amide-1", "1-3;2-1;amide-1"},
      {"1-1;2-1;3-1;triazole-1", "1-1;2-2;3-1;triazole-1",
       "1-1;2-1;3-2;triazole-1", "1-1;2-2;3-2;triazole-1",
       "1-2;2-1;3-1;triazole-1", "1-2;2-1;3-2;triazole-1",
       "1-2;2-2;3-1;triazole-1", "1-2;2-2;3-2;triazole-1"},
      {"277310376-742385846;182115391-684092275;487354835-896308859;urea-3"},
      {"1-1;2-2;3-1;triazole-1", "1-1;2-1;3-1;triazole-1",
       "1-1;2-1;3-2;triazole-1", "1-1;2-2;3-2;triazole-1",
       "1-2;2-1;3-1;triazole-1", "1-2;2-1;3-2;triazole-2",
       "1-2;2-2;3-1;triazole-1", "1-2;2-2;3-2;triazole-1"}};
  ShapeBuildParams shapeBuildOptions;
  shapeBuildOptions.numConfs = 100;
  shapeBuildOptions.rmsThreshold = 0.5;
  shapeBuildOptions.numThreads = 1;
  shapeBuildOptions.shapeSimThreshold = 0.95;
  shapeBuildOptions.randomSeed = 0xdac;

  for (size_t i = 0; i < dbNames.size(); i++) {
    SynthonSpace synthonspace;
#if 0
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
    for (const auto &mol : results.getHitMolecules()) {
      // Some of the CI tests get slightly different answers for the triazole
      // test.  It seems to be the conformer generator, even with the same
      // seed.
      if (i == 1) {
        CHECK((mol->getProp<std::string>(common_properties::_Name) ==
                   expNames[i][j] ||
               mol->getProp<std::string>(common_properties::_Name) ==
                   expNames[3][j]));
        CHECK_THAT(mol->getProp<double>("Similarity"),
                   Catch::Matchers::WithinAbs(expScores[1][j], 0.005) ||
                       Catch::Matchers::WithinAbs(expScores[3][j], 0.005));
      } else if (i == 2) {
        CHECK_THAT(mol->getProp<double>("Similarity"),
                   Catch::Matchers::WithinAbs(expScores[2][j], 0.005) ||
                       Catch::Matchers::WithinAbs(expScores[4][j], 0.005));
      } else {
        CHECK(mol->getProp<std::string>(common_properties::_Name) ==
              expNames[i][j]);
        CHECK_THAT(mol->getProp<double>("Similarity"),
                   Catch::Matchers::WithinAbs(expScores[i][j], 0.005));
      }
      auto scores = GaussianShape::ScoreMolecule(*queryMol, *mol);
      CHECK_THAT(mol->getProp<double>("Similarity"),
                 Catch::Matchers::WithinAbs(scores[0], 0.001));
      ++j;
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

TEST_CASE("Shape Callback Version") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
  std::string dbName = fullRoot + "amide_space_shapes.spc";
  SynthonSpace synthonspace;
  synthonspace.readDBFile(dbName);

  ShapeBuildParams shapeBuildOptions;
  shapeBuildOptions.numConfs = 100;
  shapeBuildOptions.rmsThreshold = 0.5;
  shapeBuildOptions.numThreads = 2;
  shapeBuildOptions.shapeSimThreshold = 0.95;
  shapeBuildOptions.randomSeed = 0xdac;

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
  auto queryMol =
      "O=C(c1ccccc1)N1CCCC1 |(0.0443291,-1.81486,-1.76886;0.0506321,-0.858174,-0.921491;1.37975,-0.430412,-0.483603;2.18964,-1.35506,0.144714;3.47088,-1.00454,0.585539;3.93803,0.297573,0.388032;3.1267,1.22739,-0.242406;1.85597,0.849751,-0.670531;-1.14837,-0.261434,-0.446583;-1.26073,0.836916,0.520219;-2.73583,1.04666,0.696614;-3.34033,-0.283345,0.290893;-2.46516,-0.679843,-0.874401)|"_smiles;

  auto results = synthonspace.shapeSearch(*queryMol, params);
  CHECK(results.getHitMolecules().size() == 3);

  // A somewhat pointless example, as a 2D SMILES isn't hugely useful
  // for a 3D search.
  std::set<std::string> resSmi;
  SearchResultCallback cb =
      [&resSmi](const std::vector<std::unique_ptr<ROMol>> &r) {
        for (auto &elem : r) {
          resSmi.insert(MolToSmiles(*elem));
        }
        return false;
      };
  synthonspace.shapeSearch(*queryMol, cb, params);
  CHECK(resSmi.size() == 3);
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

  SynthonSpace space;
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
  std::string dbName = fullRoot + "amide_space_shapes.spc";
  space.readDBFile(dbName);
  SynthonSpaceSearchParams params;
  params.similarityCutoff = 0.8;
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
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/amide_space_shapes.spc";

  SynthonSpace space;
  space.readDBFile(spaceName);
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  params.numThreads = -1;
  params.similarityCutoff = 1.0;

  auto mol =
      "O=C(c1ccccc1)N1CCCC1 |(0.0443291,-1.81486,-1.76886;0.0506321,-0.858174,-0.921491;1.37975,-0.430412,-0.483603;2.18964,-1.35506,0.144714;3.47088,-1.00454,0.585539;3.93803,0.297573,0.388032;3.1267,1.22739,-0.242406;1.85597,0.849751,-0.670531;-1.14837,-0.261434,-0.446583;-1.26073,0.836916,0.520219;-2.73583,1.04666,0.696614;-3.34033,-0.283345,0.290893;-2.46516,-0.679843,-0.874401)|"_smiles;
  REQUIRE(mol);
  auto results = space.shapeSearch(*mol, params);
  CHECK(results.getHitMolecules().empty());
  auto &bestHit = results.getBestHit();
  CHECK(bestHit);
  CHECK_NOTHROW(bestHit->getProp<std::string>(common_properties::_Name));
  CHECK(bestHit->getProp<double>("Similarity") < 1.0);
  CHECK(bestHit->getConformer().is3D());
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
    CHECK(results.getHitMolecules().size() == 2);
    std::vector<double> expScores{0.721, 0.715};
    for (unsigned int i = 0; i < results.getHitMolecules().size(); ++i) {
      auto &mol = results.getHitMolecules()[i];
      CHECK_THAT(mol->getProp<double>("Similarity"),
                 Catch::Matchers::WithinAbs(expScores[i], 0.005));
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
    CHECK(results.getHitMolecules().size() == 7);
    std::vector<double> expScores{0.801, 0.796, 0.778, 0.766,
                                  0.744, 0.717, 0.711};
    for (unsigned int i = 0; i < results.getHitMolecules().size(); ++i) {
      auto &mol = results.getHitMolecules()[i];
      CHECK_THAT(mol->getProp<double>("Similarity"),
                 Catch::Matchers::WithinAbs(expScores[i], 0.01));
    }
  }
}

// There is only 1 conformer generator available, so just do it with
// different parameters.
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
  spaceSearchParams.useProgressBar = 0;
  spaceSearchParams.timeOut = 0;
  spaceSearchParams.similarityCutoff = 0.7;
  spaceSearchParams.bestHit = true;
  auto results = space.shapeSearch(*queryMol, spaceSearchParams);
  CHECK(results.getHitMolecules().size() == 1);
  // This set of conformer parameters doesn't give a decent hit, but we
  // get something because of bestHit=true.
  for (const auto &mol : results.getHitMolecules()) {
    CHECK_THAT(mol->getProp<double>("Similarity"),
               Catch::Matchers::WithinAbs(0.728, 0.005));
    auto scores = GaussianShape::ScoreMolecule(*queryMol, *mol);
    CHECK_THAT(mol->getProp<double>("Similarity"),
               Catch::Matchers::WithinAbs(scores[0], 0.001));
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

unsigned int calcNumClashes(const ROMol &mol,
                            const GaussianShape::ShapeInput &excVol) {
  static constexpr double doubleCsq =
      4 * GaussianShape::CARBON_RAD * GaussianShape::CARBON_RAD;
  // std::cout << MolToSmiles(mol) << std::endl;
  const auto shpCoords = excVol.getCoords();
  boost::dynamic_bitset<> clashAtoms(mol.getNumAtoms());
  for (const auto atom : mol.atoms()) {
    auto aPos = mol.getConformer().getAtomPos(atom->getIdx());
    for (unsigned int i = 0; i < shpCoords.size(); i += 3) {
      const RDGeom::Point3D sPos{shpCoords[i], shpCoords[i + 1],
                                 shpCoords[i + 2]};
      if ((aPos - sPos).lengthSq() < doubleCsq) {
        // std::cout << atom->getIdx() << " -> " << i / 3 << " : "
        //           << (aPos - sPos).length() << " vs "
        //           << 2 * GaussianShape::CARBON_RAD << " : "
        //           << atom->getAtomicNum() << std::endl;
        clashAtoms[atom->getIdx()] = true;
        break;
      }
    }
  }
  return clashAtoms.count();
}

unsigned int countFileLines(const std::string &filename) {
  std::ifstream ifs(filename.c_str());
  ifs.unsetf(std::ios_base::skipws);
  unsigned int numLines = std::count(std::istream_iterator<char>(ifs),
                                     std::istream_iterator<char>(), '\n');
  ifs.close();
  return numLines;
}

TEST_CASE("Excluded volume") {
  // This is a piece of 4AJL, created by taking the 2 ligands from 4AJL and 4AJI
  // and dropping all atoms in 4AJL that were further than 5A from an atom
  // in either of the 2 ligands.  The RCSB PDB files of the 2 proteins are
  // aligned so the active sites are in the correct place.
  auto excVolMol =
      "C.C=O.CC.CC.CC(N)=O.CCC.CCC.CCC.CCC(=O)N[C@@H](C)C(=O)NCC=O.CCC(C)C.CCC(C)C.CCCC(=O)N[C@H](CN)CO.CCCCNC=O.CCCN[C@H](C=O)[C@@H](C)O.CCNC(=N)N.CNCC(N)=O.CNCCN.N.N.N=C(N)N.NC[C@H](CC(=O)O)NC=O.O.O.c1c[nH]cn1.CCC(C)C |(23.735,-5.551,4.831;10.558,-20.222,6.049;10.782,-19.593,7.082;25.756,-7.961,1.508;24.4,-8.595,1.992;7.488,-19.618,10.529;8.946,-19.661,10.098;16.988,-18.147,6.134;16.576,-18.196,7.564;15.432,-17.604,7.854;17.24,-18.815,8.408;15.229,-10.345,10.204;14.148,-9.909,9.204;13.16,-11.048,8.967;22.155,-3.456,-2.463;20.831,-2.784,-2.822;19.665,-3.765,-2.706;13.808,-17.366,17.137;12.337,-17.174,16.743;12.24,-16.466,15.411;18.52,-8.706,8.462;19.494,-9.564,7.605;19.097,-9.56,6.108;17.957,-9.849,5.775;20.037,-9.291,5.223;19.811,-9.168,3.789;20.952,-8.374,3.183;19.66,-10.496,3.038;20.161,-11.523,3.493;19.033,-10.416,1.856;18.853,-11.525,0.927;17.429,-11.999,0.744;16.5,-11.418,1.302;20.813,-11.741,-2.325;22.348,-11.427,-2.65;23.303,-11.616,-1.434;24.76,-11.817,-1.924;23.177,-10.46,-0.471;15.598,-14.535,14.266;14.867,-13.332,13.467;13.366,-13.102,13.804;13.147,-12.641,15.286;12.71,-12.136,12.79;18.012,-12.649,11.567;19.478,-12.385,11.172;20.187,-13.687,10.696;19.563,-14.252,9.404;18.599,-15.004,9.477;20.098,-13.86,8.229;19.648,-14.29,6.885;19.589,-15.803,6.809;18.634,-16.353,6.048;20.625,-13.781,5.83;20.759,-12.365,5.889;10.588,-11.836,3.288;11.37,-12.826,2.444;12.677,-13.19,3.171;13.48,-14.365,2.585;14.129,-13.976,1.342;15.412,-14.212,1.087;16.186,-14.727,1.907;7.421,-11.667,6.162;6.675,-12.99,6.581;6.654,-13.078,8.098;7.823,-13.191,8.731;7.885,-13.096,10.19;8.84,-11.926,10.495;9.77,-11.686,9.72;8.317,-14.401,10.877;7.872,-15.676,10.12;9.729,-14.357,11.106;13.973,-20.975,6.768;14.488,-21.058,8.204;13.543,-20.563,9.217;12.621,-21.314,9.805;11.859,-20.81,10.766;12.435,-22.565,9.423;14.659,-3.563,5.5;14.605,-4.881,5.293;13.35,-5.606,5.104;12.661,-5.673,6.457;11.536,-6.4,6.552;13.138,-5.069,7.429;21.154,-5.393,6.387;20.372,-4.647,5.612;19.099,-5.165,5.119;17.966,-4.53,5.886;17.138,-3.732,5.2;12.496,-8.081,8.665;12.367,-16.371,3.412;9.009,-17.463,14.046;9.374,-18.5,14.785;8.763,-18.749,15.936;10.373,-19.272,14.373;20.368,-1.798,-0.556;20.957,-1.177,0.463;20.595,-1.63,1.881;19.074,-1.555,2.174;18.319,-2.781,1.667;18.177,-2.926,0.433;17.852,-3.574,2.5;21.396,-0.863,2.83;21.647,-1.266,4.08;21.129,-2.256,4.588;26.022,-4.63,2.497;26.497,-10.537,-0.853;16.554,-20.899,12.518;16.264,-20.015,11.532;15.125,-19.344,11.93;14.762,-19.853,13.111;15.582,-20.814,13.498;27.715,-4.512,-2.12;27.605,-5.904,-2.016;28.301,-6.74,-2.887;29.098,-6.158,-3.87;28.228,-8.243,-2.741),wD:25.17,35.26,65.51,68.55,96.75,wU:40.30,49.38|"_smiles;
  // Query - ligands from 4AJI and 4AJ1, with proteins overlaid by RCSB.
  REQUIRE(excVolMol);
  GaussianShape::ShapeInputOptions excVolOptions;
  excVolOptions.useColors = false;
  auto excVolShape = std::make_unique<GaussianShape::ShapeInput>(*excVolMol, -1,
                                                                 excVolOptions);
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
  std::string dbName = fullRoot + "4ala_shapes.spc";
  SynthonSpace synthonSpace;
#if 0
  // This space is the ligand from 4al4, chopped up.  2-2 is the same as 2-1 but
  // with an extra phenyl group to ensure a big clash.  Keeping it here in
  // case it needs rebuilding.
  std::string spaceText(
      R"(SMILES	synton_id	synton#	reaction_id	release
Cc1nc2ccc(NC(=O)[1*])cc2s1	1-1	0	4al4
[1*]CCNC(=O)NCC[2*]	2-1	1	4al4
[1*]CCNC(=O)NC(c1ccccc1)C[2*]	2-2	1	4al4
[2*]Oc1ccc(cc1)CC(C(=O)O)C(=O)O	3-1	2	4al4
)");
  std::istringstream iss(spaceText);
  ShapeBuildParams shapeBuildParams;
  shapeBuildParams.numThreads = 1;
  shapeBuildParams.numConfs = 10;
  shapeBuildParams.shapeSimThreshold = -1.0;
  shapeBuildParams.randomSeed = 0xdac;
  bool cancelled = false;
  synthonSpace.readStream(iss, cancelled);
  synthonSpace.buildSynthonShapes(cancelled, shapeBuildParams);
  synthonSpace.writeDBFile(dbName);
#endif
  synthonSpace.readDBFile(dbName);

  auto comb_4aji_4aj1 =
      "CNc1nc2ccc(NC(C)=O)cc2s1.COc1ccc(CC(C(=O)O)C(=O)O)cc1OC |(23.956,-7.951,-4.055;24.081,-7.139,-2.841;23.068,-6.904,-1.967;23.237,-6.265,-0.838;22.03,-6.164,-0.135;21.811,-5.548,1.089;20.55,-5.531,1.66;19.465,-6.131,1.009;18.155,-6.121,1.583;17.318,-7.197,1.736;16.262,-7.032,2.809;17.43,-8.221,1.058;19.663,-6.743,-0.207;20.934,-6.758,-0.767;21.422,-7.467,-2.285;16.331,-12.235,6.625;15.378,-13.287,6.54;14.89,-13.803,7.711;15.669,-14.364,8.712;15.068,-14.88,9.87;13.688,-14.822,10.03;13.008,-15.402,11.226;12.32,-16.745,10.932;13.341,-17.807,10.61;14.399,-17.852,11.169;12.961,-18.68,9.696;11.481,-17.202,12.095;11.812,-18.168,12.768;10.397,-16.444,12.327;12.917,-14.275,9.014;13.494,-13.757,7.872;12.793,-13.19,6.852;11.466,-12.683,7.147)|"_smiles;

  SynthonSpaceSearchParams params;
  params.fragSimilarityAdjuster = 0.2;
  params.approxSimilarityAdjuster = 0.2;
  params.numConformers = 100;
  params.numThreads = 1;
  params.confRMSThreshold = 0.5;
  params.timeOut = 0;
  params.randomSeed = 0xdac;
  params.bestHit = true;
  params.maxHits = 10;
  params.similarityCutoff = 0.5;
  params.shapeOverlayOptions.simAlpha = 0.95;
  params.shapeOverlayOptions.simBeta = 0.05;
  params.excludedVolume = excVolShape.get();

  // This is one of those rare occasions where Mac and Linux give different
  // conformations even with the same parameters.  The results are similar
  // but ordered differently.
  std::vector<std::vector<double>> expVols{
      {74.0, 189.7}, {198.5, 74.2}, {74.3, 177.9}, {138.0, 64.6}};
  std::vector<std::vector<double>> expMeanVols{
      {3.0, 6.8}, {6.4, 3.0}, {5.7, 3.0}, {4.2, 2.7}};
  {
    params.possibleHitsFile = "poss_hits_1.txt";
    params.maxExcludedVolume = -1.0;
    auto results = synthonSpace.shapeSearch(*comb_4aji_4aj1, params);
    CHECK(results.getHitMolecules().size() == 2);
    unsigned int i = 0;
    for (const auto &mol : results.getHitMolecules()) {
      auto excVol = mol->getProp<double>("ExcludedVolume");
      CHECK_THAT(excVol, Catch::Matchers::WithinAbs(expVols[0][i], 0.1) ||
                             Catch::Matchers::WithinAbs(expVols[1][i], 0.1) ||
                             Catch::Matchers::WithinAbs(expVols[2][i], 0.1) ||
                             Catch::Matchers::WithinAbs(expVols[3][i], 0.1));
      CHECK_THAT(mol->getProp<double>("MeanExcludedVolume"),
                 Catch::Matchers::WithinAbs(expMeanVols[0][i], 0.1) ||
                     Catch::Matchers::WithinAbs(expMeanVols[1][i], 0.1) ||
                     Catch::Matchers::WithinAbs(expMeanVols[2][i], 0.1) ||
                     Catch::Matchers::WithinAbs(expMeanVols[3][i], 0.1));
      ++i;
    }
    CHECK(countFileLines("poss_hits_1.txt") == 2);
    std::remove("poss_hits_1.txt");
  }

  {
    params.possibleHitsFile = "poss_hits_2.txt";
    params.maxExcludedVolume = 100.0;
    auto results = synthonSpace.shapeSearch(*comb_4aji_4aj1, params);
    CHECK(results.getHitMolecules().size() == 1);
    CHECK(results.getHitMolecules()[0]->getProp<std::string>(
              common_properties::_Name) == "1-1;2-1;3-1;4al4");
    CHECK_THAT(results.getHitMolecules()[0]->getProp<double>("ExcludedVolume"),
               Catch::Matchers::WithinAbs(expVols[0][0], 0.1) ||
                   Catch::Matchers::WithinAbs(expVols[1][1], 0.1) ||
                   Catch::Matchers::WithinAbs(expVols[2][0], 0.1) ||
                   Catch::Matchers::WithinAbs(expVols[3][1], 0.1));
    CHECK_THAT(
        results.getHitMolecules()[0]->getProp<double>("MeanExcludedVolume"),
        Catch::Matchers::WithinAbs(expMeanVols[0][0], 0.1) ||
            Catch::Matchers::WithinAbs(expMeanVols[1][1], 0.1) ||
            Catch::Matchers::WithinAbs(expMeanVols[2][0], 0.1) ||
            Catch::Matchers::WithinAbs(expMeanVols[3][1], 0.1));
    CHECK(countFileLines("poss_hits_2.txt") == 2);
    std::remove("poss_hits_2.txt");
  }
}

TEST_CASE("Shape Write possible hits") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
  std::string dbName = fullRoot + "amide_space_shapes.spc";
  SynthonSpace synthonspace;
  synthonspace.readDBFile(dbName);

  SynthonSpaceSearchParams params;
  params.similarityCutoff = 0.8;
  params.fragSimilarityAdjuster = 0.2;
  params.approxSimilarityAdjuster = 0.2;
  params.numConformers = 100;
  params.numThreads = 2;
  params.confRMSThreshold = 0.5;
  params.timeOut = 0;
  params.randomSeed = 0xdac;
  params.bestHit = true;
  params.possibleHitsFile = "amide_space_shapes_poss_hits.txt";
  params.writePossibleHitsAndStop = false;
  auto queryMol =
      "O=C(c1ccccc1)N1CCCC1 |(0.0443291,-1.81486,-1.76886;0.0506321,-0.858174,-0.921491;1.37975,-0.430412,-0.483603;2.18964,-1.35506,0.144714;3.47088,-1.00454,0.585539;3.93803,0.297573,0.388032;3.1267,1.22739,-0.242406;1.85597,0.849751,-0.670531;-1.14837,-0.261434,-0.446583;-1.26073,0.836916,0.520219;-2.73583,1.04666,0.696614;-3.34033,-0.283345,0.290893;-2.46516,-0.679843,-0.874401)|"_smiles;

  auto results = synthonspace.shapeSearch(*queryMol, params);
  CHECK(results.getHitMolecules().size() == 3);

  std::remove(params.possibleHitsFile.c_str());
  params.writePossibleHitsAndStop = true;
  auto noresults = synthonspace.shapeSearch(*queryMol, params);
  CHECK(noresults.getHitMolecules().empty());
  CHECK(countFileLines("amide_space_shapes_poss_hits.txt") == 4);

  auto checkResults = synthonspace.shapeSearch(
      *queryMol, params, 0, std::numeric_limits<std::uint64_t>::max());
  CHECK(checkResults.getHitMolecules().size() == 3);

  auto shortResults = synthonspace.shapeSearch(*queryMol, params, 1, 3);
  CHECK(shortResults.getHitMolecules().size() == 2);
  std::remove(params.possibleHitsFile.c_str());

  params.maxPossibleHitsToWrite = 3;
  auto newResults = synthonspace.shapeSearch(*queryMol, params);
  CHECK(countFileLines("amide_space_shapes_poss_hits.txt") == 3);
  std::remove(params.possibleHitsFile.c_str());
}
