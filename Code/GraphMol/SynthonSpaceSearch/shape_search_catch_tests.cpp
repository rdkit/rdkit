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
      {1.0, 0.825, 0.813},
      {0.928, 0.873, 0.846, 0.844, 0.825, 0.820, 0.804, 0.804},
      {0.997}};
  ShapeBuildParams shapeBuildOptions;
  shapeBuildOptions.numConfs = 100;
  shapeBuildOptions.rmsThreshold = 0.5;
  shapeBuildOptions.numThreads = 1;

  for (size_t i = 0; i < dbNames.size(); i++) {
    if (i != 0) {
      continue;
    }
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
    params.numConformers = shapeBuildOptions.numConfs;
    params.numThreads = shapeBuildOptions.numThreads;
    params.confRMSThreshold = shapeBuildOptions.rmsThreshold;
    params.timeOut = 0;
    params.randomSeed = 1;
    // params.bestHit = true;
    auto queryMol = v2::SmilesParse::MolFromSmiles(querySmis[i]);
    auto results = synthonspace.shapeSearch(*queryMol, params);
    unsigned int j = 0;
    for (const auto &mol : results.getHitMolecules()) {
      std::cout << "Hit : " << MolToCXSmiles(*mol)
                << " sim = " << mol->getProp<double>("Similarity") << "  "
                << mol->getProp<std::string>("_Name") << std::endl;
      CHECK_THAT(mol->getProp<double>("Similarity"),
                 Catch::Matchers::WithinAbs(expScores[i][j++], 0.001));
      auto hitQuery = v2::SmilesParse::MolFromSmiles(
          mol->getProp<std::string>("Query_CXSmiles"));
      auto scores = GaussianShape::ScoreMolecule(*hitQuery, *mol);
      std::cout << "After score : " << scores[0] << ", " << scores[1] << ", "
                << scores[2] << std::endl;
      CHECK_THAT(mol->getProp<double>("Similarity"),
                 Catch::Matchers::WithinAbs(scores[0], 0.001));
      std::cout << "Query Conf : "
                << mol->getProp<std::string>("Query_CXSmiles") << std::endl;
    }
    CHECK(expNumHits[i] == results.getHitMolecules().size());
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
  synthonspace.buildSynthonShapes(cancelled);

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
      REQUIRE(irxn->getSynthons()[i][j].second->getShapes()->getNumShapes() ==
              orxn->getSynthons()[i][j].second->getShapes()->getNumShapes());
      for (size_t k = 0;
           k < irxn->getSynthons()[i][j].second->getShapes()->getNumShapes();
           ++k) {
        const auto ishape = irxn->getSynthons()[i][j].second->getShapes().get();
        const auto oshape = orxn->getSynthons()[i][j].second->getShapes().get();
        CHECK_THAT(fabs(ishape->getShapeVolume(k) - oshape->getShapeVolume(k)),
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
  SynthonSpace space;
  // This bit of FreedomSpace 2024-09 threw an exception originally
  // due to the Pd atom not being recognised by the shape builder.
  std::istringstream iss(R"(SMILES	synton_id	synton#	reaction_id
Cc1sc2ncnc(NN[U])c2c1C	100000182719	1	20a	2024-09
ClC(Cl)=C(Cl)c1cccc(N[U])c1	100016989384	1	20a	2024-09
Cl[Pd]c1cccc(-c2ccccc2)c1N[U]	114813040699	1	20a	2024-09
C=CCC1(S(=O)(=O)[U])CC1	100000101377	2	20a	2024-09
C=Cc1ccc(S(=O)(=O)[U])cc1	100000034458	2	20a	2024-09
)");
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
