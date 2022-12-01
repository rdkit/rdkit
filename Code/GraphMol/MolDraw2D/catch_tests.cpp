//
//  Copyright (C) 2019-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>

#include <RDGeneral/hash/hash.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/PNGParser.h>
#include <boost/algorithm/string/split.hpp>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <regex>

#ifdef RDK_BUILD_CAIRO_SUPPORT
#include <cairo.h>
#include "MolDraw2DCairo.h"
#endif

// a lot of the tests check <text> flags in the SVG.  That doesn't
// happen with the Freetype versions
static const bool NO_FREETYPE = true;

namespace {

// if the generated SVG hashes to the value we're expecting, delete
// the file.  That way, only the files that need inspection will be
// left at the end of the run.
// The hand-drawn pictures will fail this frequently due to the use
// of random numbers to draw the lines.  As well as all the testHandDrawn
// files, this includes testBrackets-5a.svg and testPositionVariation-1b.svg
static const bool DELETE_WITH_GOOD_HASH = true;
// The expected hash code for a file may be included in these maps, or
// provided in the call to check_file_hash().
// These values are for a build with FreeType, so expect them all to be
// wrong when building without.
static const std::map<std::string, std::hash_result_t> SVG_HASHES = {
    {"testAtomTags_1.svg", 146691388U},
    {"testAtomTags_2.svg", 3210393969U},
    {"testAtomTags_3.svg", 2131854465U},
    {"contourMol_1.svg", 3218870758U},
    {"contourMol_2.svg", 2353351393U},
    {"contourMol_3.svg", 3493070184U},
    {"contourMol_4.svg", 764999893U},
    {"testDativeBonds_1.svg", 221028862U},
    {"testDativeBonds_2.svg", 852819536U},
    {"testDativeBonds_3.svg", 2438158464U},
    {"testDativeBonds_2a.svg", 625232974U},
    {"testDativeBonds_2b.svg", 2879476699U},
    {"testDativeBonds_2c.svg", 388074377U},
    {"testDativeBonds_2d.svg", 1004854048U},
    {"testZeroOrderBonds_1.svg", 582365640U},
    {"testFoundations_1.svg", 767448647U},
    {"testFoundations_2.svg", 1248494165U},
    {"testTest_1.svg", 1248494165U},
    {"testKekulizationProblems_1.svg", 3747357148U},
    {"testAtomBondIndices_1.svg", 3247225756U},
    {"testAtomBondIndices_2.svg", 2812688864U},
    {"testAtomBondIndices_3.svg", 3653054459U},
    {"testAtomBondIndices_4.svg", 635195563U},
    {"testAtomBondIndices_5.svg", 23594241U},
    {"testAtomBondIndices_6.svg", 1127540948U},
    {"testGithub3226_1.svg", 4099803989U},
    {"testGithub3226_2.svg", 3134615187U},
    {"testGithub3226_3.svg", 1488097417U},
    {"testGithub3369_1.svg", 4165525787U},
    {"testIncludeRadicals_1a.svg", 2528551797U},
    {"testIncludeRadicals_1b.svg", 3075507489U},
    {"testLegendsAndDrawing-1.svg", 1693176512U},
    {"testGithub3577-1.svg", 3974438540U},
    {"testHandDrawn-1.svg", 799391905U},
    {"testHandDrawn-2.svg", 2605087576U},
    {"testHandDrawn-3.svg", 1015633173U},
    {"testHandDrawn-4.svg", 830784921U},
    {"testHandDrawn-5a.svg", 2845825621U},
    {"testHandDrawn-5b.svg", 476521352U},
    {"testBrackets-1a.svg", 3257646535U},
    {"testBrackets-1b.svg", 776088825U},
    {"testBrackets-1c.svg", 3257646535U},
    {"testBrackets-1d.svg", 776088825U},
    {"testBrackets-1e.svg", 1202405256U},
    {"testBrackets-2a.svg", 728321376U},
    {"testBrackets-2b.svg", 1408188695U},
    {"testBrackets-2c.svg", 728321376U},
    {"testBrackets-2d.svg", 1408188695U},
    {"testBrackets-3a.svg", 791450653U},
    {"testBrackets-4a.svg", 769125635U},
    {"testBrackets-4b.svg", 4066682338U},
    {"testBrackets-5a.svg", 1388227932U},
    {"testBrackets-5768.svg", 3070888879U},
    {"testSGroupData-1a.svg", 1463366807U},
    {"testSGroupData-1b.svg", 223883202U},
    {"testSGroupData-2a.svg", 3547547260U},
    {"testSGroupData-2b.svg", 2573013307U},
    {"testSGroupData-3a.svg", 2220120573U},
    {"testPositionVariation-1.svg", 4185441744U},
    {"testPositionVariation-1b.svg", 2588110577U},
    {"testPositionVariation-2.svg", 2026425280U},
    {"testPositionVariation-3.svg", 56671878U},
    {"testPositionVariation-4.svg", 886758688U},
    {"testNoAtomLabels-1.svg", 2648234379U},
    {"testNoAtomLabels-2.svg", 3213096674U},
    {"testQueryBonds-1a.svg", 3288272531U},
    {"testQueryBonds-1b.svg", 1706839957U},
    {"testQueryBonds-1c.svg", 333519907U},
    {"testQueryBonds-2.svg", 69341882U},
    {"testLinkNodes-2-0.svg", 2952965907U},
    {"testLinkNodes-2-30.svg", 4117540200U},
    {"testLinkNodes-2-60.svg", 520576199U},
    {"testLinkNodes-2-90.svg", 1403605120U},
    {"testLinkNodes-2-120.svg", 3607355853U},
    {"testLinkNodes-2-150.svg", 177350824U},
    {"testLinkNodes-2-180.svg", 3809030739U},
    {"testMolAnnotations-1.svg", 1091624544U},
    {"testMolAnnotations-2a.svg", 2203886283U},
    {"testMolAnnotations-2b.svg", 400443600U},
    {"testMolAnnotations-2c.svg", 3954034822U},
    {"testMolAnnotations-3a.svg", 1752047273U},
    {"testMolAnnotations-3b.svg", 2068377089U},
    {"testMolAnnotations-3c.svg", 4288169182U},
    {"testMolAnnotations-3d.svg", 1514775509U},
    {"testMolAnnotations-4a.svg", 569955128U},
    {"testLinkNodes-1-0.svg", 2929724949U},
    {"testLinkNodes-1-30.svg", 800625899U},
    {"testLinkNodes-1-60.svg", 609508172U},
    {"testLinkNodes-1-90.svg", 2665766032U},
    {"testLinkNodes-1-120.svg", 992989277U},
    {"testLinkNodes-1-150.svg", 1392691166U},
    {"testLinkNodes-1-180.svg", 130695597U},
    {"testGithub3744.svg", 2774492807U},
    {"testAtomLists-1.svg", 2751373083U},
    {"testAtomLists-2.svg", 385738799U},
    {"testIsoDummyIso.svg", 1696129196U},
    {"testNoIsoDummyIso.svg", 2004687512U},
    {"testIsoNoDummyIso.svg", 2734544682U},
    {"testNoIsoNoDummyIso.svg", 918094584U},
    {"testDeuteriumTritium.svg", 2634768249U},
    {"testHydrogenBonds1.svg", 4137715598U},
    {"testHydrogenBonds2.svg", 2044702263U},
    {"testGithub3912.1.svg", 3081580881U},
    {"testGithub3912.2.svg", 1662866562U},
    {"testGithub2976.svg", 703667023U},
    {"testReactionCoords.svg", 2325796920U},
    {"testAnnotationColors.svg", 445523422U},
    {"testGithub4323_1.svg", 1993234598U},
    {"testGithub4323_2.svg", 2933922429U},
    {"testGithub4323_3.svg", 1773544359U},
    {"testGithub4323_4.svg", 213795827U},
    {"testGithub4238_1.svg", 629357140U},
    {"testGithub4508_1.svg", 3784765069U},
    {"testGithub4508_1b.svg", 3433942203U},
    {"testGithub4508_2.svg", 326155865U},
    {"testGithub4508_2b.svg", 662225995U},
    {"testGithub4538.svg", 3198623323U},
    {"testDarkMode.1.svg", 1977391752U},
    {"testMonochrome.1.svg", 1776897420U},
    {"testMonochrome.2.svg", 399259780U},
    {"testAvalon.1.svg", 1614166818U},
    {"testCDK.1.svg", 3108685638U},
    {"testGithub4519_1.svg", 473230604U},
    {"testGithub4519_2.svg", 2515716875U},
    {"testGithub4519_3.svg", 1017109741U},
    {"testGithub4519_4.svg", 645908829U},
    {"testBaseFontSize.1a.svg", 3939288880U},
    {"testBaseFontSize.1b.svg", 2617787443U},
    {"testBaseFontSize.2a.svg", 1031690455U},
    {"testBaseFontSize.2b.svg", 3440038194U},
    {"testFlexiCanvas.1a.svg", 3145560884U},
    {"testFlexiCanvas.1b.svg", 1140847713U},
    {"testFlexiCanvas.1c.svg", 2832891200U},
    {"testFlexiCanvas.1d.svg", 4220526884U},
    {"testFlexiCanvas.2.svg", 1185770886U},
    {"testSemiFlexiCanvas.1a.svg", 414967968U},
    {"testSemiFlexiCanvas.1b.svg", 367831852U},
    {"testSemiFlexiCanvas.1c.svg", 316673185U},
    {"testFlexiCanvas.3.svg", 1164132085U},
    {"testFlexiCanvas.4a.svg", 438150211U},
    {"testFlexiCanvas.4b.svg", 2015277207U},
    {"testFlexiCanvas.4c.svg", 3138663789U},
    {"testFlexiCanvas.4d.svg", 1950746506U},
    {"testFlexiCanvas.5a.svg", 1204456580U},
    {"testFlexiCanvas.5b.svg", 4164471763U},
    {"testFlexiCanvas.5c.svg", 2381227232U},
    {"testFlexiCanvas.5d.svg", 2157866153U},
    {"testFlexiCanvas.6a.svg", 4104973953U},
    {"testFlexiCanvas.6b.svg", 2392263541U},
    {"testFlexiCanvas.6c.svg", 4104973953U},
    {"testFlexiCanvas.6d.svg", 4104973953U},
    {"testFlexiCanvas.7a.svg", 918094125U},
    {"testFlexiCanvas.7b.svg", 4094511140U},
    {"testFlexiCanvas.7c.svg", 918094125U},
    {"testFlexiCanvas.7d.svg", 918094125U},
    {"testGithub4764.sz1.svg", 493786705U},
    {"testGithub4764.sz2.svg", 2704253898U},
    {"testGithub4764.sz3.svg", 1328896014U},
    {"testDrawArc1.svg", 4039810147U},
    {"testMetalWedges.svg", 3278785383U},
    {"testVariableLegend_1.svg", 3914441319U},
    {"testVariableLegend_2.svg", 3458084009U},
    {"testVariableLegend_3.svg", 1996551457U},
    {"testGithub_5061.svg", 1947248304U},
    {"testGithub_5185.svg", 2944445711U},
    {"testGithub_5269_1.svg", 2368496794U},
    {"testGithub_5269_2.svg", 567813292U},
    {"test_classes_wavy_bonds.svg", 1271445012U},
    {"testGithub_5383_1.svg", 1391972140U},
    {"github5156_1.svg", 695855770U},
    {"github5156_2.svg", 2606649270U},
    {"github5156_3.svg", 3284451122U},
    {"test_molblock_wedges.svg", 1106580037U},
    {"github5383_1.svg", 2353351393U},
    {"acs1996_1.svg", 51426601U},
    {"acs1996_2.svg", 833573044U},
    {"acs1996_3.svg", 4007912653U},
    {"acs1996_4.svg", 3372558370U},
    {"acs1996_5.svg", 2883542240U},
    {"acs1996_6.svg", 1380727178U},
    {"acs1996_7.svg", 763391533U},
    {"acs1996_8.svg", 939325262U},
    {"acs1996_9.svg", 2607143500U},
    {"acs1996_10.svg", 199499735U},
    {"acs1996_11.svg", 2121789178U},
    {"acs1996_12.svg", 2233727631U},
    {"test_unspec_stereo.svg", 599119798U},
    {"light_blue_h_no_label_1.svg", 3735371135U},
    {"test_github_5534.svg", 574501211U},
    {"bond_highlights_1.svg", 1426179967U},
    {"bond_highlights_2.svg", 3654242474U},
    {"bond_highlights_3.svg", 2068128924U},
    {"bond_highlights_4.svg", 4115973245U},
    {"bond_highlights_5.svg", 4115973245U},
    {"bond_highlights_6.svg", 1566801788U},
    {"bond_highlights_7.svg", 2101261688U},
    {"bond_highlights_8.svg", 3826056528U},
    {"bond_highlights_9.svg", 2915809284U},
    {"testGithub5486_1.svg", 1149144091U},
    {"testGithub5511_1.svg", 940106456U},
    {"testGithub5511_2.svg", 1448975272U},
    {"test_github5767.svg", 3153964439U}};

// These PNG hashes aren't completely reliable due to floating point cruft,
// but they can still reduce the number of drawings that need visual
// inspection.  At present, the files
// testPNGMetadata_2.png
// give different results on my MBP and Ubuntu 20.04 VM.  The SVGs work
// better because the floats are all output to only 1 decimal place so there
// is a much smaller chance of different systems producing different files.
static const std::map<std::string, std::hash_result_t> PNG_HASHES = {
    {"testGithub3226_1.png", 284815097U},
    {"testGithub3226_2.png", 2460913971U},
    {"testGithub3226_3.png", 993799198U},
    {"testPNGMetadata_1.png", 2022143293U},
    {"testPNGMetadata_2.png", 3078435362U},
    {"testHandDrawn-1.png", 1551605661U},
    {"testHandDrawn-2.png", 2979412913U},
    {"testHandDrawn-3.png", 1765396301U},
    {"testHandDrawn-4.png", 2989933219U},
    {"testHandDrawn-5.png", 1526220279U},
    {"testGithub4323_1.png", 3711520691U},
    {"testGithub4323_3.png", 2300228708U},
    {"testFlexiCanvas.2a.png", 3618977786U},
    {"testFlexiCanvas.2b.png", 2780757414U},
    {"testGithub4764.sz1.png", 2320783268U},
    {"testGithub4764.sz2.png", 3297570843U},
    {"testGithub4764.sz3.png", 2178018272U},
    {"testGithub4238_1.png", 458925131U},
    {"github5383_1.png", 2963331215U},
    {"acs1996_1.png", 2674458798U},
    {"acs1996_2.png", 83755168U}};

std::hash_result_t hash_file(const std::string &filename) {
  std::ifstream ifs(filename, std::ios_base::binary);
  std::string file_contents(std::istreambuf_iterator<char>{ifs}, {});
  if (filename.substr(filename.length() - 4) == ".svg") {
    // deal with MSDOS newlines.
    file_contents.erase(
        remove(file_contents.begin(), file_contents.end(), '\r'),
        file_contents.end());
  }
  return gboost::hash_range(file_contents.begin(), file_contents.end());
}

void check_file_hash(const std::string &filename,
                     std::hash_result_t exp_hash = 0U) {
  //    std::cout << filename << " : " << hash_file(filename) << "U" <<
  //    std::endl;

  std::map<std::string, std::hash_result_t>::const_iterator it;
  if (filename.substr(filename.length() - 4) == ".svg") {
    it = SVG_HASHES.find(filename);
  } else {
    it = PNG_HASHES.find(filename);
  }
  std::hash_result_t file_hash = hash_file(filename);
  if (exp_hash == 0U) {
    exp_hash = it == SVG_HASHES.end() ? 0U : it->second;
  }
  if (it != SVG_HASHES.end() && file_hash == exp_hash) {
    if (DELETE_WITH_GOOD_HASH) {
      std::remove(filename.c_str());
    }
  } else {
    std::cout << "file " << filename << " gave hash " << file_hash
              << "U not the expected " << exp_hash << "U" << std::endl;
  }
}
}  // namespace

using namespace RDKit;

TEST_CASE("prepareAndDrawMolecule", "[drawing]") {
  SECTION("basics") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);

    // we will be able to recognize that the prep worked because there
    // will be an H in the output:
    MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    CHECK(text.find(">H</text>") != std::string::npos);
  }
  SECTION("kekulize") {
    auto m1 = "c1ccccc1"_smiles;
    REQUIRE(m1);

    {
      MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      CHECK(text.find("stroke-dasharray") == std::string::npos);
    }
    {
      MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1, "", nullptr, nullptr,
                                             nullptr, nullptr, nullptr, -1,
                                             false);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      CHECK(text.find("stroke-dasharray") != std::string::npos);
    }
  }
}

TEST_CASE("tag atoms in SVG", "[drawing][SVG]") {
  SECTION("basics") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);

    MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawMolecule(*m1);
    std::map<std::string, std::string> actions;
    actions["onclick"] = "alert";
    double radius = 0.2;
    drawer.tagAtoms(*m1, radius, actions);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testAtomTags_1.svg");
    outs << text;
    outs.close();
    check_file_hash("testAtomTags_1.svg");

    CHECK(text.find("<circle") != std::string::npos);
    CHECK(text.find("<circle") != std::string::npos);
    CHECK(text.find("atom-selector") != std::string::npos);
    CHECK(text.find("bond-selector") != std::string::npos);
  }
  SECTION("inject prop to class") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);

    for (auto atom : m1->atoms()) {
      auto prop = boost::format("__prop_class_atom_%d") % atom->getIdx();
      atom->setProp("_tagClass", prop.str());
    }
    for (auto bond : m1->bonds()) {
      auto prop = boost::format("__prop_class_bond_%d") % bond->getIdx();
      bond->setProp("_tagClass", prop.str());
    }

    MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawMolecule(*m1);
    drawer.tagAtoms(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testAtomTags_2.svg");
    outs << text;
    outs.close();
    check_file_hash("testAtomTags_2.svg");

    size_t i = 0;
    size_t c = 0;
    while (true) {
      auto i2 = text.find("__prop_class_atom_", i);
      if (i2 == std::string::npos) {
        break;
      }
      i = i2 + 1;
      c++;
    }
    CHECK(c == 6);

    i = 0;
    c = 0;
    while (true) {
      auto i2 = text.find("__prop_class_bond_", i);
      if (i2 == std::string::npos) {
        break;
      }
      i = i2 + 1;
      c++;
    }
    CHECK(c == 7);
  }
}

TEST_CASE("metadata in SVG", "[drawing][SVG]") {
  SECTION("inject prop to metada") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);

    for (auto atom : m1->atoms()) {
      auto prop = boost::format("__prop_metadata_atom_%d") % atom->getIdx();
      atom->setProp("_metaData-atom-inject-prop", prop.str());
    }
    for (auto bond : m1->bonds()) {
      auto prop = boost::format("__prop_metadata_bond_%d") % bond->getIdx();
      bond->setProp("_metaData-bond-inject-prop", prop.str());
    }

    MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawMolecule(*m1);
    drawer.addMoleculeMetadata(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testAtomTags_3.svg");
    outs << text;
    outs.close();
    check_file_hash("testAtomTags_3.svg");

    size_t i = 0;
    size_t c = 0;
    while (true) {
      auto i2 = text.find("atom-inject-prop=\"__prop_metadata_atom_", i);
      if (i2 == std::string::npos) {
        break;
      }
      i = i2 + 1;
      c++;
    }
    CHECK(c == 6);

    i = 0;
    c = 0;
    while (true) {
      auto i2 = text.find("bond-inject-prop=\"__prop_metadata_bond_", i);
      if (i2 == std::string::npos) {
        break;
      }
      i = i2 + 1;
      c++;
    }
    CHECK(c == 7);
  }
}

TEST_CASE("contour data", "[drawing][conrec]") {
  auto m1 = "C1N[C@@H]2OCC12"_smiles;
  REQUIRE(m1);
  SECTION("grid basics") {
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m1);

    const size_t gridSz = 100;
    auto *grid = new double[gridSz * gridSz];
    std::vector<double> xps(gridSz);
    std::vector<double> yps(gridSz);

    double minX = 1000, minY = 1000, maxX = -1000, maxY = -1000;
    const auto conf = m1->getConformer();
    for (size_t i = 0; i < conf.getNumAtoms(); ++i) {
      minX = std::min(minX, conf.getAtomPos(i).x);
      minY = std::min(minY, conf.getAtomPos(i).y);
      maxX = std::max(maxX, conf.getAtomPos(i).x);
      maxY = std::max(maxY, conf.getAtomPos(i).y);
    }
    double x1 = minX - 0.5, y1 = minY - 0.5, x2 = maxX + 0.5, y2 = maxY + 0.5;
    double dx = (x2 - x1) / gridSz, dy = (y2 - y1) / gridSz;
    double maxV = 0.0;
    for (size_t ix = 0; ix < gridSz; ++ix) {
      auto px = x1 + ix * dx;
      xps[ix] = px;
      for (size_t iy = 0; iy < gridSz; ++iy) {
        auto py = y1 + iy * dy;
        if (ix == 0) {
          yps[iy] = py;
        }
        RDGeom::Point2D loc(px, py);
        double val = 0.0;
        for (size_t ia = 0; ia < conf.getNumAtoms(); ++ia) {
          auto dv = loc - RDGeom::Point2D(conf.getAtomPos(ia).x,
                                          conf.getAtomPos(ia).y);
          auto r = dv.length();
          if (r > 0.1) {
            val += 1 / r;
          }
        }
        maxV = std::max(val, maxV);
        grid[ix * gridSz + iy] = val;
      }
    }

    std::vector<double> levels;
    drawer.clearDrawing();
    MolDraw2DUtils::contourAndDrawGrid(drawer, grid, xps, yps, 10, levels,
                                       MolDraw2DUtils::ContourParams(),
                                       m1.get());
    drawer.drawOptions().clearBackground = false;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("contourMol_1.svg");
    outs << text;
    outs.close();
    check_file_hash("contourMol_1.svg");
    delete[] grid;
  }
  SECTION("gaussian basics") {
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawOptions().padding = 0.1;

    const auto conf = m1->getConformer();
    std::vector<Point2D> cents(conf.getNumAtoms());
    std::vector<double> weights(conf.getNumAtoms());
    std::vector<double> widths(conf.getNumAtoms());
    for (size_t i = 0; i < conf.getNumAtoms(); ++i) {
      cents[i] = Point2D(conf.getAtomPos(i).x, conf.getAtomPos(i).y);
      weights[i] = 1;
      widths[i] = 0.4 * PeriodicTable::getTable()->getRcovalent(
                            m1->getAtomWithIdx(i)->getAtomicNum());
    }

    std::vector<double> levels;
    drawer.clearDrawing();
    MolDraw2DUtils::contourAndDrawGaussians(
        drawer, cents, weights, widths, 10, levels,
        MolDraw2DUtils::ContourParams(), m1.get());

    drawer.drawOptions().clearBackground = false;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("contourMol_2.svg");
    outs << text;
    outs.close();
    check_file_hash("contourMol_2.svg");
  }
  SECTION("gaussian fill") {
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawOptions().padding = 0.1;

    const auto conf = m1->getConformer();
    std::vector<Point2D> cents(conf.getNumAtoms());
    std::vector<double> weights(conf.getNumAtoms());
    std::vector<double> widths(conf.getNumAtoms());
    for (size_t i = 0; i < conf.getNumAtoms(); ++i) {
      cents[i] = Point2D(conf.getAtomPos(i).x, conf.getAtomPos(i).y);
      weights[i] = i % 2 ? -0.5 : 1;
      widths[i] = 0.4 * PeriodicTable::getTable()->getRcovalent(
                            m1->getAtomWithIdx(i)->getAtomicNum());
    }

    std::vector<double> levels;
    MolDraw2DUtils::ContourParams cps;
    cps.fillGrid = true;
    drawer.clearDrawing();
    MolDraw2DUtils::contourAndDrawGaussians(drawer, cents, weights, widths, 10,
                                            levels, cps, m1.get());

    drawer.drawOptions().clearBackground = false;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("contourMol_3.svg");
    outs << text;
    outs.close();
    check_file_hash("contourMol_3.svg");
  }

  SECTION("gaussian fill 2") {
    auto m2 = "C1N[C@@H]2OCC12C=CC"_smiles;
    REQUIRE(m2);

    MolDraw2DSVG drawer(450, 250, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m2);
    drawer.drawOptions().padding = 0.1;

    const auto conf = m2->getConformer();
    std::vector<Point2D> cents(conf.getNumAtoms());
    std::vector<double> weights(conf.getNumAtoms());
    std::vector<double> widths(conf.getNumAtoms());
    for (size_t i = 0; i < conf.getNumAtoms(); ++i) {
      cents[i] = Point2D(conf.getAtomPos(i).x, conf.getAtomPos(i).y);
      weights[i] = i % 2 ? -0.5 : 1;
      widths[i] = 0.3 * PeriodicTable::getTable()->getRcovalent(
                            m2->getAtomWithIdx(i)->getAtomicNum());
    }

    std::vector<double> levels;
    MolDraw2DUtils::ContourParams cps;
    cps.fillGrid = true;
    cps.gridResolution = 0.5;
    drawer.clearDrawing();
    MolDraw2DUtils::contourAndDrawGaussians(drawer, cents, weights, widths, 10,
                                            levels, cps, m2.get());

    drawer.drawOptions().clearBackground = false;
    drawer.drawMolecule(*m2);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("contourMol_4.svg");
    outs << text;
    outs.close();
    check_file_hash("contourMol_4.svg");
  }
}

TEST_CASE("dative bonds", "[drawing][organometallics]") {
  SECTION("basics") {
    auto m1 = "N->[Pt]"_smiles;
    REQUIRE(m1);
    MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testDativeBonds_1.svg");
    outs << text;
    outs.close();
    check_file_hash("testDativeBonds_1.svg");

    CHECK(text.find("d='M 122.5,88.4 L 85.6,88.4' "
                    "style='fill:none;fill-rule:evenodd;stroke:#0000FF") !=
          std::string::npos);
  }
  SECTION("more complex") {
    auto m1 = "N->1[C@@H]2CCCC[C@H]2N->[Pt]11OC(=O)C(=O)O1"_smiles;
    REQUIRE(m1);
    MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testDativeBonds_2.svg");
    outs << text;
    outs.close();
    check_file_hash("testDativeBonds_2.svg");

    CHECK(text.find("-8' d='M 101.1,79.8 L 95.8,87.1' "
                    "style='fill:none;fill-rule:evenodd;stroke:#0000FF;") !=
          std::string::npos);
  }
  SECTION("test colours") {
    // the dative bonds point the wrong way, but the point is to test
    // if the tip of the arrow is blue.
    auto m1 = "[Cu++]->1->2.N1CCN2"_smiles;
    REQUIRE(m1);
    MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testDativeBonds_3.svg");
    outs << text;
    outs.close();
    check_file_hash("testDativeBonds_3.svg");

    CHECK(
        text.find(
            "<path class='bond-2 atom-3 atom-4' d='M 50.4,140.6 L 77.9,149.5' "
            "style='fill:none;fill-rule:evenodd;stroke:#0000FF;") !=
        std::string::npos);
  }
  SECTION("dative series") {
    auto m1 = "N->1[C@@H]2CCCC[C@H]2N->[Pt]11OC(=O)C(=O)O1"_smiles;
    REQUIRE(m1);
    {
      MolDraw2DSVG drawer(150, 150, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareMolForDrawing(*m1);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testDativeBonds_2a.svg");
      outs << text;
      outs.close();
      check_file_hash("testDativeBonds_2a.svg");
    }
    {
      MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareMolForDrawing(*m1);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testDativeBonds_2b.svg");
      outs << text;
      outs.close();
      check_file_hash("testDativeBonds_2b.svg");
    }
    {
      MolDraw2DSVG drawer(350, 350, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareMolForDrawing(*m1);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testDativeBonds_2c.svg");
      outs << text;
      outs.close();
      check_file_hash("testDativeBonds_2c.svg");
    }
    {
      MolDraw2DSVG drawer(450, 450, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareMolForDrawing(*m1);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testDativeBonds_2d.svg");
      outs << text;
      outs.close();
      check_file_hash("testDativeBonds_2d.svg");
    }
  }
}

TEST_CASE("zero-order bonds", "[drawing][organometallics]") {
  SECTION("basics") {
    auto m1 = "N-[Pt]"_smiles;
    REQUIRE(m1);
    m1->getBondWithIdx(0)->setBondType(Bond::ZERO);
    MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testZeroOrderBonds_1.svg");
    outs << text;
    outs.close();
    check_file_hash("testZeroOrderBonds_1.svg");

    CHECK(text.find("stroke-dasharray:2,2") != std::string::npos);
  }
}

TEST_CASE("copying drawing options", "[drawing]") {
  auto m1 = "C1N[C@@H]2OCC12"_smiles;
  REQUIRE(m1);
  SECTION("foundations") {
    {
      MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testFoundations_1.svg");
      outs << text;
      outs.close();
      check_file_hash("testFoundations_1.svg");
      CHECK(text.find("fill:#0000FF' >N</text>") != std::string::npos);
    }
    {
      MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
      assignBWPalette(drawer.drawOptions().atomColourPalette);
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testFoundations_2.svg");
      outs << text;
      outs.close();
      check_file_hash("testFoundations_2.svg");
      CHECK(text.find("fill:#0000FF' >N</text>") == std::string::npos);
      CHECK(text.find("fill:#000000' >N</text>") != std::string::npos);
    }
  }
  SECTION("test") {
    {
      MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
      MolDrawOptions options = drawer.drawOptions();
      assignBWPalette(options.atomColourPalette);
      drawer.drawOptions() = options;
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testTest_1.svg");
      outs << text;
      outs.close();
      check_file_hash("testTest_1.svg");
      CHECK(text.find("fill:#0000FF' >N</text>") == std::string::npos);
      CHECK(text.find("fill:#000000' >N</text>") != std::string::npos);
    }
  }
}

TEST_CASE("bad DrawMolecules() when molecules are not kekulized",
          "[drawing][bug]") {
  auto m1 = "CCN(CC)CCn1nc2c3ccccc3sc3c(CNS(C)(=O)=O)ccc1c32"_smiles;
  REQUIRE(m1);
  SECTION("foundations") {
    MolDraw2DSVG drawer(500, 200, 250, 200, NO_FREETYPE);
    drawer.drawOptions().prepareMolsBeforeDrawing = false;
    RWMol dm1(*m1);
    RWMol dm2(*m1);
    bool kekulize = false;
    MolDraw2DUtils::prepareMolForDrawing(dm1, kekulize);
    kekulize = true;
    MolDraw2DUtils::prepareMolForDrawing(dm2, kekulize);
    MOL_PTR_VECT ms{&dm1, &dm2};
    drawer.drawMolecule(dm1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testKekulizationProblems_1.svg");
    outs << text;
    outs.close();
    check_file_hash("testKekulizationProblems_1.svg");

    // this is a very crude test - really we just need to look at the SVG - but
    // it's better than nothing.
    CHECK(text.find(
              "<path class='bond-18' d='M 169.076,79.056 L 191.285,69.2653' "
              "style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:"
              "2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;"
              "stroke-dasharray:6,6' />") == std::string::npos);
  }
}
TEST_CASE("draw atom/bond indices", "[drawing]") {
  auto m1 = "C[C@H](F)N"_smiles;
  auto m2 = "C[C@@H](F)N"_smiles;
  REQUIRE(m1);
  REQUIRE(m2);
  SECTION("foundations") {
    {
      MolDraw2DSVG drawer(250, 200, -1, -1, NO_FREETYPE);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_1.svg");
      outs << text;
      outs.close();
      check_file_hash("testAtomBondIndices_1.svg");
      CHECK(text.find(">1</text>") == std::string::npos);
      CHECK(text.find(">(</text>") == std::string::npos);
      CHECK(text.find(">S</text>") == std::string::npos);
      CHECK(text.find(">)</text>") == std::string::npos);
    }
    {
      MolDraw2DSVG drawer(250, 200, -1, -1, NO_FREETYPE);
      drawer.drawOptions().addAtomIndices = true;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_2.svg");
      outs << text;
      outs.close();
      check_file_hash("testAtomBondIndices_2.svg");
      CHECK(text.find(">1</text>") != std::string::npos);
      // it only appears once though:
      CHECK(text.find(">1</text>", text.find(">1</text>") + 1) ==
            std::string::npos);
      CHECK(text.find("1,(S)") == std::string::npos);
    }
    {
      MolDraw2DSVG drawer(250, 200, -1, -1, NO_FREETYPE);
      drawer.drawOptions().addBondIndices = true;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_3.svg");
      outs << text;
      outs.close();
      check_file_hash("testAtomBondIndices_3.svg");
      CHECK(text.find(">1</text>") != std::string::npos);
      // it only appears once though:
      CHECK(text.find(">1</text>", text.find(">1</text>") + 1) ==
            std::string::npos);
    }
    {
      MolDraw2DSVG drawer(250, 200, -1, -1, NO_FREETYPE);
      drawer.drawOptions().addAtomIndices = true;
      drawer.drawOptions().addBondIndices = true;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_4.svg");
      outs << text;
      outs.close();
      check_file_hash("testAtomBondIndices_4.svg");
      CHECK(text.find(">1</text>") != std::string::npos);
      // it appears twice:
      CHECK(text.find(">1</text>", text.find(">1</text>") + 1) !=
            std::string::npos);
    }
    {
      MolDraw2DSVG drawer(250, 200, -1, -1, NO_FREETYPE);
      m1->getAtomWithIdx(2)->setProp(common_properties::atomNote, "foo");
      drawer.drawOptions().addAtomIndices = true;
      drawer.drawOptions().addStereoAnnotation = true;
      drawer.drawMolecule(*m1);
      m1->getAtomWithIdx(2)->clearProp(common_properties::atomNote);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_5.svg");
      outs << text;
      outs.close();
      check_file_hash("testAtomBondIndices_5.svg");
      CHECK(text.find(">1</text>") != std::string::npos);
      CHECK(text.find(">,</text>") != std::string::npos);
      CHECK(text.find(">(</text>") != std::string::npos);
      CHECK(text.find(">S</text>") != std::string::npos);
      CHECK(text.find(")</text>") != std::string::npos);
      CHECK(text.find(">2</text>") != std::string::npos);
      CHECK(text.find(">f</text>") != std::string::npos);
      CHECK(text.find(">o</text>") != std::string::npos);
    }
    {
      // Make sure it works for solid wedges as well.
      MolDraw2DSVG drawer(250, 200, -1, -1, NO_FREETYPE);
      drawer.drawOptions().addAtomIndices = true;
      drawer.drawMolecule(*m2);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_6.svg");
      outs << text;
      outs.close();
      check_file_hash("testAtomBondIndices_6.svg");
      CHECK(text.find(">1</text>") != std::string::npos);
      // it only appears once though:
      CHECK(text.find(">1</text>", text.find(">1</text>") + 1) ==
            std::string::npos);
      CHECK(text.find("1,(S)") == std::string::npos);
    }
  }
}

TEST_CASE("Github #3226: Lines in wedge bonds being drawn too closely together",
          "[drawing]") {
  auto m1 =
      "C[C@H](C1=C(C=CC(=C1Cl)F)Cl)OC2=C(N=CC(=C2)C3=CN(N=C3)C4CCNCC4)N"_smiles;
  REQUIRE(m1);
  SECTION("larger SVG") {
    {
      MolDraw2DSVG drawer(450, 400);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub3226_1.svg");
      outs << text;
      outs.close();
      check_file_hash("testGithub3226_1.svg");
      std::vector<std::string> tkns;
      boost::algorithm::find_all(tkns, text, "bond-0");
      CHECK(tkns.size() == 10);
    }
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  SECTION("larger PNG") {
    {
      MolDraw2DCairo drawer(450, 400);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      drawer.writeDrawingText("testGithub3226_1.png");
      check_file_hash("testGithub3226_1.png");
    }
  }
#endif
  SECTION("smaller SVG") {
    {
      MolDraw2DSVG drawer(200, 150);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub3226_2.svg");
      outs << text;
      outs.close();
      check_file_hash("testGithub3226_2.svg");
      std::vector<std::string> tkns;
      boost::algorithm::find_all(tkns, text, "bond-0");
      CHECK(tkns.size() == 5);
    }
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  SECTION("smaller PNG") {
    {
      MolDraw2DCairo drawer(200, 150);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      drawer.writeDrawingText("testGithub3226_2.png");
      check_file_hash("testGithub3226_2.png");
    }
  }
#endif
  SECTION("middle SVG") {
    {
      MolDraw2DSVG drawer(300, 200);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub3226_3.svg");
      outs << text;
      outs.close();
      check_file_hash("testGithub3226_3.svg");
      std::vector<std::string> tkns;
      boost::algorithm::find_all(tkns, text, "bond-0");
      CHECK(tkns.size() == 7);
    }
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  SECTION("middle PNG") {
    {
      MolDraw2DCairo drawer(250, 200);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      drawer.writeDrawingText("testGithub3226_3.png");
      check_file_hash("testGithub3226_3.png");
    }
  }
#endif
}

TEST_CASE("github #3258: ", "[drawing][bug]") {
  auto m1 = "CCN"_smiles;
  REQUIRE(m1);
  SECTION("foundations") {
    MolDraw2DSVG drawer(500, 200, 250, 200, NO_FREETYPE);
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawOptions().addBondIndices = true;
    RWMol dm1(*m1);
    RWMol dm2(*m1);
    MOL_PTR_VECT ms{&dm1, &dm2};
    drawer.drawMolecules(ms);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    CHECK(text.find(">,</text>") == std::string::npos);
    CHECK(!dm1.hasProp("_atomIndicesAdded"));
    CHECK(!dm1.hasProp("_bondIndicesAdded"));
  }
}

#ifdef RDK_BUILD_CAIRO_SUPPORT
TEST_CASE("adding png metadata", "[drawing][png]") {
  SECTION("molecule") {
    auto m1 = R"CTAB(
  Mrv2014 08172015242D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 2.31 -1.3337 0 0
M  V30 2 C 3.6437 -2.1037 0 0
M  V30 3 O 4.9774 -1.3337 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m1);
    {
      MolDraw2DCairo drawer(250, 200);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      auto png = drawer.getDrawingText();
      drawer.writeDrawingText("testPNGMetadata_1.png");
      check_file_hash("testPNGMetadata_1.png");
      CHECK(png.find(PNGData::smilesTag) != std::string::npos);
      CHECK(png.find(PNGData::molTag) != std::string::npos);
      CHECK(png.find(PNGData::pklTag) != std::string::npos);
      std::unique_ptr<ROMol> newmol(PNGStringToMol(png));
      REQUIRE(newmol);
      CHECK(MolToCXSmiles(*m1) == MolToCXSmiles(*newmol));
    }
    {  // disable metadata output
      MolDraw2DCairo drawer(250, 200);
      drawer.drawOptions().includeMetadata = false;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      auto png = drawer.getDrawingText();
      CHECK(png.find(PNGData::smilesTag) == std::string::npos);
      CHECK(png.find(PNGData::molTag) == std::string::npos);
      CHECK(png.find(PNGData::pklTag) == std::string::npos);
    }
    {  // draw multiple molecules
      MolDraw2DCairo drawer(250, 200);
      drawer.drawMolecule(*m1);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      auto png = drawer.getDrawingText();
      CHECK(png.find(PNGData::smilesTag) != std::string::npos);
      CHECK(png.find(PNGData::molTag) != std::string::npos);
      CHECK(png.find(PNGData::pklTag) != std::string::npos);
      CHECK(png.find(PNGData::smilesTag + "1") != std::string::npos);
      CHECK(png.find(PNGData::molTag + "1") != std::string::npos);
      CHECK(png.find(PNGData::pklTag + "1") != std::string::npos);
    }
  }
  SECTION("reaction") {
    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
        "[N:1][C:2][C:3](=[O:4])[O:5].[N:6][C:7][C:8](=[O:9])[O:10]>>[N:1]1[C:"
        "2][C:3](=[O:4])[N:6][C:7][C:8]1=[O:9].[O:5][O:10]"));
    REQUIRE(rxn);
    {
      MolDraw2DCairo drawer(600, 200);
      drawer.drawReaction(*rxn);
      drawer.finishDrawing();
      auto png = drawer.getDrawingText();
      drawer.writeDrawingText("testPNGMetadata_2.png");
      check_file_hash("testPNGMetadata_2.png");
      CHECK(png.find(PNGData::smilesTag) == std::string::npos);
      CHECK(png.find(PNGData::molTag) == std::string::npos);
      CHECK(png.find(PNGData::pklTag) == std::string::npos);
      CHECK(png.find(PNGData::rxnPklTag) != std::string::npos);
      CHECK(png.find(PNGData::rxnSmartsTag) != std::string::npos);
      std::unique_ptr<ChemicalReaction> rxn2(PNGStringToChemicalReaction(png));
      REQUIRE(rxn2);
      CHECK(ChemicalReactionToRxnSmarts(*rxn) ==
            ChemicalReactionToRxnSmarts(*rxn2));
    }
    {  // disable metadata
      MolDraw2DCairo drawer(600, 200);
      drawer.drawOptions().includeMetadata = false;
      drawer.drawReaction(*rxn);
      drawer.finishDrawing();
      auto png = drawer.getDrawingText();
      CHECK(png.find(PNGData::smilesTag) == std::string::npos);
      CHECK(png.find(PNGData::molTag) == std::string::npos);
      CHECK(png.find(PNGData::pklTag) == std::string::npos);
      CHECK(png.find(PNGData::rxnPklTag) == std::string::npos);
      CHECK(png.find(PNGData::rxnSmartsTag) == std::string::npos);
    }
  }
}

#endif

TEST_CASE(
    "github #3392: prepareMolForDrawing() incorrectly adds chiral Hs if no "
    "ring info is present",
    "[bug]") {
  SECTION("foundations") {
    SmilesParserParams ps;
    ps.sanitize = false;
    ps.removeHs = false;
    std::unique_ptr<RWMol> m1(SmilesToMol("C[C@H](F)Cl", ps));
    REQUIRE(m1);
    m1->updatePropertyCache();
    CHECK(m1->getNumAtoms() == 4);
    const bool kekulize = false;
    const bool addChiralHs = true;
    MolDraw2DUtils::prepareMolForDrawing(*m1, kekulize, addChiralHs);
    CHECK(m1->getNumAtoms() == 4);
  }
}

TEST_CASE(
    "github #3369: support new CIP code and StereoGroups in "
    "addStereoAnnotation()",
    "[chirality]") {
  auto m1 =
      "C[C@@H]1N[C@H](C)[C@@H]([C@H](C)[C@@H]1C)C1[C@@H](C)O[C@@H](C)[C@@H](C)[C@H]1C/C=C/C |a:5,o1:1,8,o2:14,16,&1:11,18,&2:3,6,r|"_smiles;
  REQUIRE(m1);
  SECTION("defaults") {
    ROMol m2(*m1);
    MolDraw2D_detail::addStereoAnnotation(m2);

    std::string txt;
    CHECK(m2.getAtomWithIdx(5)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "abs (S)");
    CHECK(m2.getAtomWithIdx(3)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "and4");
  }
  SECTION("including CIP with relative stereo") {
    ROMol m2(*m1);
    bool includeRelativeCIP = true;
    MolDraw2D_detail::addStereoAnnotation(m2, includeRelativeCIP);

    std::string txt;
    CHECK(m2.getAtomWithIdx(5)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "abs (S)");
    CHECK(m2.getAtomWithIdx(3)->getPropIfPresent(common_properties::atomNote,
                                                 txt));
    CHECK(txt == "and4 (R)");
  }
  SECTION("new CIP labels") {
    ROMol m2(*m1);
    REQUIRE(m2.getBondBetweenAtoms(20, 21));
    m2.getBondBetweenAtoms(20, 21)->setStereo(Bond::BondStereo::STEREOTRANS);
    // initially no label is assigned since we have TRANS
    MolDraw2D_detail::addStereoAnnotation(m2);
    CHECK(
        !m2.getBondBetweenAtoms(20, 21)->hasProp(common_properties::bondNote));

    CIPLabeler::assignCIPLabels(m2);
    std::string txt;
    CHECK(m2.getBondBetweenAtoms(20, 21)->getPropIfPresent(
        common_properties::_CIPCode, txt));
    CHECK(txt == "E");
    MolDraw2D_detail::addStereoAnnotation(m2);
    CHECK(m2.getBondBetweenAtoms(20, 21)->getPropIfPresent(
        common_properties::bondNote, txt));
    CHECK(txt == "(E)");
  }
  SECTION("works with the drawing code") {
    MolDraw2DSVG drawer(300, 250);
    RWMol dm1(*m1);
    bool includeRelativeCIP = true;
    MolDraw2D_detail::addStereoAnnotation(dm1, includeRelativeCIP);
    drawer.drawMolecule(dm1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub3369_1.svg");
    outs << text;
    outs.close();
    check_file_hash("testGithub3369_1.svg");
  }
}

TEST_CASE("includeRadicals", "[options]") {
  SECTION("basics") {
    auto m = "[O][C]"_smiles;
    REQUIRE(m);
    int panelHeight = -1;
    int panelWidth = -1;
    bool noFreeType = true;
    {
      MolDraw2DSVG drawer(250, 200, panelWidth, panelHeight, noFreeType);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testIncludeRadicals_1a.svg");
      outs << text;
      outs.close();
      check_file_hash("testIncludeRadicals_1a.svg");
      CHECK(text.find("<path class='atom-0' d='M") != std::string::npos);
    }
    {
      MolDraw2DSVG drawer(250, 200, panelWidth, panelHeight, noFreeType);
      drawer.drawOptions().includeRadicals = false;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testIncludeRadicals_1b.svg");
      outs << text;
      outs.close();
      check_file_hash("testIncludeRadicals_1b.svg");
      CHECK(text.find("<path class='atom-0' d='M") == std::string::npos);
    }
  }
}

TEST_CASE("including legend in drawing results in offset drawing later",
          "[bug]") {
  SECTION("basics") {
    auto m = "c1ccccc1"_smiles;
    REQUIRE(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    auto &conf = m->getConformer();
    std::vector<Point2D> polyg;
    for (const auto &pt : conf.getPositions()) {
      polyg.emplace_back(pt);
    }
    MolDraw2DSVG drawer(350, 300);
    drawer.drawMolecule(*m, "molecule legend");
    drawer.setFillPolys(true);
    drawer.setColour(DrawColour(1.0, 0.3, 1.0));
    drawer.drawPolygon(polyg);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testLegendsAndDrawing-1.svg");
    outs << text;
    outs.close();
    outs.close();
    check_file_hash("testLegendsAndDrawing-1.svg");

    // make sure the polygon starts at a bond
    CHECK(text.find("<path class='bond-0 atom-0 atom-1' d='M 316.7,135.0") !=
          std::string::npos);
    CHECK(text.find("<path d='M 316.7,135.0") != std::string::npos);
  }
}

TEST_CASE("Github #3577", "[bug]") {
  SECTION("basics") {
    auto m = "CCC"_smiles;
    REQUIRE(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    m->getAtomWithIdx(1)->setProp("atomNote", "CCC");
    m->getAtomWithIdx(2)->setProp("atomNote", "ccc");
    m->getBondWithIdx(0)->setProp("bondNote", "CCC");

    MolDraw2DSVG drawer(350, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testGithub3577-1.svg");
    outs << text;
    outs.close();
    check_file_hash("testGithub3577-1.svg");
  }
}
TEST_CASE("hand drawn", "[play]") {
  SECTION("basics") {
    auto m =
        "CC[CH](C)[CH]1NC(=O)[CH](Cc2ccc(O)cc2)NC(=O)[CH](N)CSSC[CH](C(=O)N2CCC[CH]2C(=O)N[CH](CC(C)C)C(=O)NCC(N)=O)NC(=O)[CH](CC(N)=O)NC(=O)[CH](CCC(N)=O)NC1=O"_smiles;
    REQUIRE(m);
    RDDepict::preferCoordGen = true;
    MolDraw2DUtils::prepareMolForDrawing(*m);

    std::string fName = getenv("RDBASE");
    fName += "/Data/Fonts/ComicNeue-Regular.ttf";

    {
      MolDraw2DSVG drawer(450, 400);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m, "Oxytocin (flat)");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testHandDrawn-1.svg");
      outs << text;
      outs.close();
      check_file_hash("testHandDrawn-1.svg");
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(450, 400);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m, "Oxytocin (flat)");
      drawer.finishDrawing();
      drawer.writeDrawingText("testHandDrawn-1.png");
      check_file_hash("testHandDrawn-1.png");
    }
#endif
  }
  SECTION("with chirality") {
    auto m =
        "CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H](N)CSSC[C@@H](C(=O)N2CCC[C@H]2C(=O)N[C@@H](CC(C)C)C(=O)NCC(N)=O)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCC(N)=O)NC1=O"_smiles;
    REQUIRE(m);
    RDDepict::preferCoordGen = true;
    MolDraw2DUtils::prepareMolForDrawing(*m);

    std::string fName = getenv("RDBASE");
    fName += "/Data/Fonts/ComicNeue-Regular.ttf";

    {
      MolDraw2DSVG drawer(450, 400);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m, "Oxytocin");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testHandDrawn-2.svg");
      outs << text;
      outs.close();
      check_file_hash("testHandDrawn-2.svg");
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(450, 400);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m, "Oxytocin");
      drawer.finishDrawing();
      drawer.writeDrawingText("testHandDrawn-2.png");
      check_file_hash("testHandDrawn-2.png");
    }
#endif
  }
  SECTION("smaller") {
    auto m = "N=c1nc([C@H]2NCCCC2)cc(N)n1O"_smiles;
    REQUIRE(m);
    RDDepict::preferCoordGen = true;
    MolDraw2DUtils::prepareMolForDrawing(*m);

    std::string fName = getenv("RDBASE");
    fName += "/Data/Fonts/ComicNeue-Regular.ttf";

    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testHandDrawn-3.svg");
      outs << text;
      outs.close();
      check_file_hash("testHandDrawn-3.svg");
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(350, 300);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText("testHandDrawn-3.png");
      check_file_hash("testHandDrawn-3.png");
    }
#endif
  }
  SECTION("another one") {
    auto m =
        "CCCc1nn(C)c2c(=O)nc(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)[nH]c12"_smiles;
    REQUIRE(m);
    RDDepict::preferCoordGen = true;
    MolDraw2DUtils::prepareMolForDrawing(*m);

    std::string fName = getenv("RDBASE");
    fName += "/Data/Fonts/ComicNeue-Regular.ttf";

    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testHandDrawn-4.svg");
      outs << text;
      outs.close();
      check_file_hash("testHandDrawn-4.svg");
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(350, 300);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText("testHandDrawn-4.png");
      check_file_hash("testHandDrawn-4.png");
    }
#endif
  }
  SECTION("large") {
    auto m =
        "CC[C@H](C)[C@@H](C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc1ccc(cc1)O)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)NCC(=O)N[C@@H](CCC(=O)N)C(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCCCN)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](CC(C)C)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H]3CCCN3C(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CO)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@H](CO)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCSC)NC(=O)[C@H](Cc4ccccc4)NC(=O)CNC(=O)CNC(=O)[C@H](Cc5ccc(cc5)O)N"_smiles;
    REQUIRE(m);
    RDDepict::preferCoordGen = true;
    MolDraw2DUtils::prepareMolForDrawing(*m);

    std::string fName = getenv("RDBASE");
    fName += "/Data/Fonts/ComicNeue-Regular.ttf";

    {
      MolDraw2DSVG drawer(900, 450);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testHandDrawn-5a.svg");
      outs << text;
      outs.close();
      check_file_hash("testHandDrawn-5a.svg");
    }
    {
      MolDraw2DSVG drawer(900, 450);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testHandDrawn-5b.svg");
      outs << text;
      outs.close();
      check_file_hash("testHandDrawn-5b.svg");
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(900, 450);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText("testHandDrawn-5.png");
      check_file_hash("testHandDrawn-5.png");
    }
#endif
  }
}

TEST_CASE("drawMoleculeBrackets", "[extras]") {
  SECTION("basics") {
    auto m = R"CTAB(
  ACCLDraw11042015112D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 7 -6.7813 0 0 
M  V30 2 C 8.0229 -6.1907 0 0 CFG=3 
M  V30 3 C 8.0229 -5.0092 0 0 
M  V30 4 C 9.046 -6.7814 0 0 
M  V30 5 C 10.0692 -6.1907 0 0 
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 
M  V30 2 1 2 3 
M  V30 3 1 2 4 
M  V30 4 1 4 5 
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 1 ATOMS=(3 3 2 4) XBONDS=(2 1 4) BRKXYZ=(9 7.51 -7.08 0 7.51 -
M  V30 -5.9 0 0 0 0) BRKXYZ=(9 9.56 -5.9 0 9.56 -7.08 0 0 0 0) -
M  V30 CONNECT=HT LABEL=n 
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-1a.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-1a.svg");
    }
    {  // rotation
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().rotate = 90;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-1b.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-1b.svg");
    }
    {  // centering
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().centreMoleculesBeforeDrawing = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-1c.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-1c.svg");
    }
    {  // rotation + centering
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().centreMoleculesBeforeDrawing = true;
      drawer.drawOptions().rotate = 90;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-1d.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-1d.svg");
    }
    {  // rotation
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().rotate = 180;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-1e.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-1e.svg");
    }
  }
  SECTION("three brackets") {
    auto m = R"CTAB(three brackets
  Mrv2014 11052006542D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 1 0 0
M  V30 BEGIN ATOM
M  V30 1 * -1.375 3.1667 0 0
M  V30 2 C -0.0413 3.9367 0 0
M  V30 3 C 1.2924 3.1667 0 0
M  V30 4 * 2.626 3.9367 0 0
M  V30 5 C 0.0003 5.6017 0 0
M  V30 6 * 1.334 6.3717 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 2 5
M  V30 5 1 5 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(3 2 3 5) XBONDS=(3 1 3 5) BRKXYZ=(9 0.0875 6.7189 0 -
M  V30 1.0115 5.1185 0 0 0 0) BRKXYZ=(9 1.3795 4.2839 0 2.3035 2.6835 0 0 0 -
M  V30 0) BRKXYZ=(9 -0.1285 2.8194 0 -1.0525 4.4198 0 0 0 0) CONNECT=HT -
M  V30 LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-2a.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-2a.svg");
    }
    {  // rotation
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().rotate = 90;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-2b.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-2b.svg");
    }
    {  // centering
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().centreMoleculesBeforeDrawing = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-2c.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-2c.svg");
    }
    {  // rotation + centering
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().centreMoleculesBeforeDrawing = true;
      drawer.drawOptions().rotate = 90;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-2d.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-2d.svg");
    }
  }
  SECTION("ChEBI 59342") {
    // thanks to John Mayfield for pointing out the example
    auto m = R"CTAB(ChEBI59342 
Marvin  05041012302D          

 29 30  0  0  1  0            999 V2000
   10.1615   -7.7974    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.7305   -6.9763    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.7309   -7.8004    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    9.4464   -8.2109    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    8.0153   -8.2225    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.4464   -9.0437    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    8.0138   -9.0500    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    8.7293   -9.4606    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
   10.1669   -9.4529    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.3058   -9.4590    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.7368  -10.2801    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    8.0263  -10.6992    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.0339  -11.5241    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    7.3081  -10.2933    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.7305   -5.3264    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.0159   -5.7369    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    8.0159   -6.5618    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    7.2936   -5.3263    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.2936   -6.9762    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    6.5751   -5.7368    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    6.5751   -6.5618    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    7.2973   -7.8049    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.8681   -5.3263    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.8680   -6.9762    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    5.1510   -6.5684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4392   -6.9856    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.1455   -5.7435    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.4142   -5.3560    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
   11.5590   -7.8297    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
  3  2  1  6  0  0  0
  3  4  1  0  0  0  0
  3  5  1  0  0  0  0
  4  6  1  0  0  0  0
  4  1  1  1  0  0  0
  5  7  1  0  0  0  0
  6  8  1  0  0  0  0
  6  9  1  1  0  0  0
  7 10  1  1  0  0  0
  8 11  1  6  0  0  0
  7  8  1  0  0  0  0
 13 12  1  0  0  0  0
 14 12  2  0  0  0  0
 11 12  1  0  0  0  0
 16 15  1  6  0  0  0
 16 17  1  0  0  0  0
 16 18  1  0  0  0  0
 17 19  1  0  0  0  0
 17  2  1  1  0  0  0
 18 20  1  0  0  0  0
 19 21  1  0  0  0  0
 19 22  1  1  0  0  0
 20 23  1  1  0  0  0
 21 24  1  6  0  0  0
 20 21  1  0  0  0  0
 26 25  1  0  0  0  0
 27 25  2  0  0  0  0
 24 25  1  0  0  0  0
 15 28  1  0  0  0  0
  1 29  1  0  0  0  0
M  STY  1   1 SRU
M  SCN  1   1 HT 
M  SAL   1 15   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
M  SAL   1 12  16  17  18  19  20  21  22  23  24  25  26  27
M  SDI   1  4    9.4310   -4.9261    9.4165   -5.7510
M  SDI   1  4   10.7464   -7.3983   10.7274   -8.2231
M  SBL   1  2  30  29
M  SMT   1 n
M  END)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-3a.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-3a.svg");
    }
  }
  SECTION("pathological bracket orientation") {
    {  // including the bonds
      auto m = R"CTAB(bogus
  Mrv2014 11202009512D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 1 0 1
M  V30 BEGIN ATOM
M  V30 1 C 23.5462 -14.464 0 0
M  V30 2 C 20.8231 -13.0254 0 0
M  V30 3 C 20.8776 -14.5628 0 0
M  V30 4 C 22.2391 -15.2819 0 0
M  V30 5 C 16.2969 -9.9426 0 0
M  V30 6 C 14.963 -10.7089 0 0
M  V30 7 C 19.463 -12.2987 0 0
M  V30 8 * 19.4398 -9.9979 0 0
M  V30 9 * 26.1554 -14.4332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 3 4
M  V30 2 1 6 7
M  V30 3 1 5 8
M  V30 4 1 1 9
M  V30 5 1 7 2
M  V30 6 1 6 5
M  V30 7 1 4 1
M  V30 8 1 3 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(7 4 3 7 6 5 2 1) XBONDS=(2 3 4) BRKXYZ=(9 17.6045 -
M  V30 -9.1954 0 17.5775 -10.7352 0 0 0 0) BRKXYZ=(9 24.6113 -13.6813 0 -
M  V30 24.6296 -15.2213 0 0 0 0) CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-4a.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-4a.svg");
    }

    {  // no bonds in the sgroup, the bracket should point the other way
       // (towards the majority of the atoms in the sgroup)
      auto m = R"CTAB(bogus
  Mrv2014 11202009512D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 1 0 1
M  V30 BEGIN ATOM
M  V30 1 C 23.5462 -14.464 0 0
M  V30 2 C 20.8231 -13.0254 0 0
M  V30 3 C 20.8776 -14.5628 0 0
M  V30 4 C 22.2391 -15.2819 0 0
M  V30 5 C 16.2969 -9.9426 0 0
M  V30 6 C 14.963 -10.7089 0 0
M  V30 7 C 19.463 -12.2987 0 0
M  V30 8 * 19.4398 -9.9979 0 0
M  V30 9 * 26.1554 -14.4332 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 3 4
M  V30 2 1 6 7
M  V30 3 1 5 8
M  V30 4 1 1 9
M  V30 5 1 7 2
M  V30 6 1 6 5
M  V30 7 1 4 1
M  V30 8 1 3 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(7 4 3 7 6 5 2 1) BRKXYZ=(9 17.6045 -
M  V30 -9.1954 0 17.5775 -10.7352 0 0 0 0) BRKXYZ=(9 24.6113 -13.6813 0 -
M  V30 24.6296 -15.2213 0 0 0 0) CONNECT=HT LABEL=n
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-4b.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-4b.svg");
    }
  }
  SECTION("comic brackets (no font though)") {
    auto m = R"CTAB(
  ACCLDraw11042015112D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 7 -6.7813 0 0 
M  V30 2 C 8.0229 -6.1907 0 0 CFG=3 
M  V30 3 C 8.0229 -5.0092 0 0 
M  V30 4 C 9.046 -6.7814 0 0 
M  V30 5 C 10.0692 -6.1907 0 0 
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 
M  V30 2 1 2 3 
M  V30 3 1 2 4 
M  V30 4 1 4 5 
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 1 ATOMS=(3 3 2 4) XBONDS=(2 1 4) BRKXYZ=(9 7.51 -7.08 0 7.51 -
M  V30 -5.9 0 0 0 0) BRKXYZ=(9 9.56 -5.9 0 9.56 -7.08 0 0 0 0) -
M  V30 CONNECT=HT LABEL=n 
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-5a.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-5a.svg");
    }
  }
  SECTION("Github5768 - rightmost bracket wrong way round.)") {
    auto m = R"CTAB(
  Marvin  10140911012D

 19 18  0  0  0  0            999 V2000
   -2.0296    1.6372    0.0000 Si  0  0  0  0  0  0  0  0  0  0  0  0
   -2.0296    2.4622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0296    0.8122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8546    1.6372    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2046    1.6372    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3796    1.6372    0.0000 Si  0  0  0  0  0  0  0  0  0  0  0  0
    0.4454    1.6372    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3796    0.8122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3796    2.4622    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3349    0.3997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0494    0.8122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7638    0.3997    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2704    1.6372    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0
    2.4783    0.8122    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1928    0.3996    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9072    0.8121    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.6217    0.3996    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3362    0.8120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0506    0.3995    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  6  8  1  0  0  0  0
  6  9  1  0  0  0  0
  8 10  1  0  0  0  0
 10 11  1  0  0  0  0
  7 13  1  0  0  0  0
 11 12  1  0  0  0  0
 12 14  1  0  0  0  0
 14 15  1  0  0  0  0
 15 16  1  0  0  0  0
 16 17  1  0  0  0  0
 17 18  1  0  0  0  0
 18 19  1  0  0  0  0
M  STY  2   1 SRU   2 SRU
M  SCN  1   1 HT
M  SAL   1  4   1   2   3   5
M  SDI   1  4   -0.8649    2.0497   -0.8649    1.2247
M  SDI   1  4   -2.3693    1.2247   -2.3693    2.0497
M  SBL   1  2   3   5
M  SMT   1 n
M  SCN  1   2 HT
M  SAL   2 13   6   7   8   9  10  11  12  14  15  16  17  18  19
M  SDI   2  4    0.7851    2.0497    0.7851    1.2247
M  SDI   2  4   -0.7193    1.2247   -0.7193    2.0497
M  SBL   2  2   5  11
M  SMT   2 m
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().addAtomIndices = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-5768.svg");
      outs << text;
      outs.close();
      check_file_hash("testBrackets-5768.svg");
    }
  }
}

#ifdef RDK_BUILD_CAIRO_SUPPORT
TEST_CASE("github #3543: Error adding PNG metadata when kekulize=False",
          "[bug][metadata][png]") {
  SECTION("basics") {
    auto m = "n1cccc1"_smarts;
    m->updatePropertyCache(false);
    MolDraw2DCairo drawer(350, 300);
    bool kekulize = false;
    MolDraw2DUtils::prepareMolForDrawing(*m, kekulize);
    drawer.drawOptions().prepareMolsBeforeDrawing = false;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    auto png = drawer.getDrawingText();
  }
  SECTION("as reported") {
    auto m = "n1cnc2c(n)ncnc12"_smarts;
    m->updatePropertyCache(false);
    MolDraw2DCairo drawer(350, 300);
    bool kekulize = false;
    MolDraw2DUtils::prepareMolForDrawing(*m, kekulize);
    drawer.drawOptions().prepareMolsBeforeDrawing = false;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    auto png = drawer.getDrawingText();
  }
}
#endif

TEST_CASE("SGroup Data") {
  SECTION("ABS") {
    auto m = R"CTAB(
  Mrv2014 12072015352D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6.5833 4.3317 0 0
M  V30 2 C -7.917 3.5617 0 0
M  V30 3 C -7.917 2.0216 0 0
M  V30 4 C -6.5833 1.2516 0 0
M  V30 5 C -5.2497 2.0216 0 0
M  V30 6 C -5.2497 3.5617 0 0
M  V30 7 C -3.916 4.3317 0 0
M  V30 8 O -3.916 5.8717 0 0
M  V30 9 O -2.5823 3.5617 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 6 7
M  V30 8 2 7 8
M  V30 9 1 7 9
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 9) FIELDNAME=pKa -
M  V30 FIELDDISP="   -2.2073    2.3950    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=4.2
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m, "abs");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testSGroupData-1a.svg");
      outs << text;
      outs.close();
      check_file_hash("testSGroupData-1a.svg");
    }
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().centreMoleculesBeforeDrawing = true;
      drawer.drawOptions().rotate = 90;
      drawer.drawMolecule(*m, "centered, rotated");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testSGroupData-1b.svg");
      outs << text;
      outs.close();
      check_file_hash("testSGroupData-1b.svg");
    }
  }
  SECTION("REL") {
    auto m = R"CTAB(
  Mrv2014 12072015352D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -6.5833 4.3317 0 0
M  V30 2 C -7.917 3.5617 0 0
M  V30 3 C -7.917 2.0216 0 0
M  V30 4 C -6.5833 1.2516 0 0
M  V30 5 C -5.2497 2.0216 0 0
M  V30 6 C -5.2497 3.5617 0 0
M  V30 7 C -3.916 4.3317 0 0
M  V30 8 O -3.916 5.8717 0 0
M  V30 9 O -2.5823 3.5617 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 6 7
M  V30 8 2 7 8
M  V30 9 1 7 9
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 9) FIELDNAME=pKa -
M  V30 FIELDDISP="    0.2000    0.2000    DRU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=4.2
M  V30 END SGROUP
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m, "rel");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testSGroupData-2a.svg");
      outs << text;
      outs.close();
      check_file_hash("testSGroupData-2a.svg");
    }
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().centreMoleculesBeforeDrawing = true;
      drawer.drawOptions().rotate = 90;
      drawer.drawMolecule(*m, "rel, centered, rotated");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testSGroupData-2b.svg");
      outs << text;
      outs.close();
      check_file_hash("testSGroupData-2b.svg");
    }
  }
  {
    auto m = R"CTAB(random example found on internet
   JSDraw204221719232D

 20 21  0  0  0  0              0 V2000
   10.1710   -5.6553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.9428   -4.2996    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    8.6110   -5.6647    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.9591   -7.0015    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   12.5190   -6.9921    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   13.3072   -8.3384    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   13.2909   -5.6364    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.5028   -4.2902    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   13.2746   -2.9345    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   14.8508   -5.6270    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   15.6226   -4.2713    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   20.3026   -4.2431    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   19.5307   -5.5987    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   21.8625   -4.2336    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   19.5144   -2.8968    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   17.9544   -2.9062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.1663   -1.5600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   17.1826   -4.2619    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.9708   -5.6082    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   17.1989   -6.9638    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  2  0  0  0  0
  1  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  5  7  1  0  0  0  0
  7  8  1  0  0  0  0
  8  9  2  0  0  0  0
  8  2  1  0  0  0  0
  7 10  1  0  0  0  0
 10 11  2  0  0  0  0
 12 13  1  0  0  0  0
 12 14  2  0  0  0  0
 12 15  1  0  0  0  0
 15 16  1  0  0  0  0
 16 17  2  0  0  0  0
 16 18  1  0  0  0  0
 18 19  1  0  0  0  0
 19 20  1  0  0  0  0
 19 13  2  0  0  0  0
 11 18  1  0  0  0  0
M  STY  1   1 DAT
M  SDT   1 UNKNOWN                        F
M  SDD   1    16.0856   -8.1573    DA    ALL  1       5
M  SED   1 Ni-complex
M  END)CTAB"_ctab;
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testSGroupData-3a.svg");
      outs << text;
      outs.close();
      check_file_hash("testSGroupData-3a.svg");
    }
  }
}

TEST_CASE("position variation bonds", "[extras]") {
  SECTION("simple") {
    auto m = R"CTAB(
  Mrv2014 12092006072D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.7083 4.915 0 0
M  V30 2 C -6.042 4.145 0 0
M  V30 3 C -6.042 2.605 0 0
M  V30 4 C -4.7083 1.835 0 0
M  V30 5 C -3.3747 2.605 0 0
M  V30 6 C -3.3747 4.145 0 0
M  V30 7 * -3.8192 3.8883 0 0
M  V30 8 O -3.8192 6.1983 0 0
M  V30 9 C -2.4855 6.9683 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 8 ENDPTS=(3 1 6 5) ATTACH=ANY
M  V30 8 1 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m, "variations");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testPositionVariation-1.svg");
      outs << text;
      outs.close();
      check_file_hash("testPositionVariation-1.svg");
    }
    {  // make sure comic mode doesn't screw this up
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m, "comic variations");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testPositionVariation-1b.svg");
      outs << text;
      outs.close();
      check_file_hash("testPositionVariation-1b.svg");
    }
  }
  SECTION("multiple") {
    auto m = R"CTAB(
  Mrv2014 12092006082D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 15 14 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.7083 4.915 0 0
M  V30 2 C -6.042 4.145 0 0
M  V30 3 C -6.042 2.605 0 0
M  V30 4 C -4.7083 1.835 0 0
M  V30 5 C -3.3747 2.605 0 0
M  V30 6 C -3.3747 4.145 0 0
M  V30 7 * -3.8192 3.8883 0 0
M  V30 8 O -3.8192 6.1983 0 0
M  V30 9 C -2.4855 6.9683 0 0
M  V30 10 C -7.3757 4.915 0 0
M  V30 11 C -8.7093 4.145 0 0
M  V30 12 C -8.7093 2.605 0 0
M  V30 13 C -7.3757 1.835 0 0
M  V30 14 * -8.7093 3.375 0 0
M  V30 15 O -10.2922 3.375 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 8 ENDPTS=(3 1 6 5) ATTACH=ANY
M  V30 8 1 8 9
M  V30 9 1 10 11
M  V30 10 2 11 12
M  V30 11 1 12 13
M  V30 12 2 10 2
M  V30 13 2 13 3
M  V30 14 1 14 15 ENDPTS=(2 11 12) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m, "multiple variations");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testPositionVariation-2.svg");
      outs << text;
      outs.close();
      check_file_hash("testPositionVariation-2.svg");
    }
  }
  SECTION("non-contiguous") {
    auto m = R"CTAB(
  Mrv2014 12092006102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.875 8.7484 0 0
M  V30 2 C -2.2087 7.9784 0 0
M  V30 3 C -2.2087 6.4383 0 0
M  V30 4 C -0.875 5.6683 0 0
M  V30 5 C 0.4587 6.4383 0 0
M  V30 6 C 0.4587 7.9784 0 0
M  V30 7 * -0.4304 6.9517 0 0
M  V30 8 O -0.4304 4.6417 0 0
M  V30 9 C -1.7641 3.8717 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 8 ENDPTS=(3 1 5 4) ATTACH=ANY
M  V30 8 1 8 9
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m, "non-contiguous atoms");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testPositionVariation-3.svg");
      outs << text;
      outs.close();
      check_file_hash("testPositionVariation-3.svg");
    }
  }
  SECTION("larger mol") {
    auto m = R"CTAB(
  Mrv2014 12092009152D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 23 24 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.875 8.7484 0 0
M  V30 2 C -2.2087 7.9784 0 0
M  V30 3 C -2.2087 6.4383 0 0
M  V30 4 C -0.875 5.6683 0 0
M  V30 5 N 0.4587 6.4383 0 0
M  V30 6 C 0.4587 7.9784 0 0
M  V30 7 * -0.4304 6.9517 0 0
M  V30 8 O -0.4304 4.6417 0 0
M  V30 9 C -1.7641 3.8717 0 0
M  V30 10 C -3.5423 8.7484 0 0
M  V30 11 C -4.876 7.9784 0 0
M  V30 12 C -4.876 6.4383 0 0
M  V30 13 C -3.5423 5.6683 0 0
M  V30 14 C -4.876 11.0584 0 0
M  V30 15 C -6.2097 10.2884 0 0
M  V30 16 C -6.2097 8.7484 0 0
M  V30 17 C -3.5423 10.2884 0 0
M  V30 18 C -6.2097 13.3685 0 0
M  V30 19 C -7.5433 12.5985 0 0
M  V30 20 C -7.5433 11.0584 0 0
M  V30 21 C -4.876 12.5985 0 0
M  V30 22 * -5.5428 9.1334 0 0
M  V30 23 C -7.3712 7.7304 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 8 ENDPTS=(3 1 4 5) ATTACH=ANY
M  V30 8 1 8 9
M  V30 9 2 10 11
M  V30 10 1 11 12
M  V30 11 2 12 13
M  V30 12 1 10 2
M  V30 13 1 13 3
M  V30 14 1 14 15
M  V30 15 2 15 16
M  V30 16 2 14 17
M  V30 17 1 10 17
M  V30 18 1 16 11
M  V30 19 1 18 19
M  V30 20 2 19 20
M  V30 21 2 18 21
M  V30 22 1 14 21
M  V30 23 1 20 15
M  V30 24 1 22 23 ENDPTS=(2 15 11) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(250, 200);
      drawer.drawMolecule(*m, "smaller");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testPositionVariation-4.svg");
      outs << text;
      outs.close();
      check_file_hash("testPositionVariation-4.svg");
    }
  }
}

TEST_CASE("disable atom labels", "[feature]") {
  SECTION("basics") {
    {
      auto m = "NCC(=O)O"_smiles;
      MolDraw2DSVG drawer(350, 300);
      MolDraw2DUtils::prepareMolForDrawing(*m);
      drawer.drawOptions().noAtomLabels = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testNoAtomLabels-1.svg");
      outs << text;
      outs.close();
      check_file_hash("testNoAtomLabels-1.svg");
      CHECK(text.find("class='atom-0") == std::string::npos);
      CHECK(text.find("class='atom-3") == std::string::npos);
    }
    {
      auto m = "F[C@H](O)C[C@@H](Cl)I"_smiles;
      MolDraw2DSVG drawer(350, 300);
      MolDraw2DUtils::prepareMolForDrawing(*m);
      drawer.drawOptions().noAtomLabels = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testNoAtomLabels-2.svg");
      outs << text;
      outs.close();
      check_file_hash("testNoAtomLabels-2.svg");
    }
  }
}

TEST_CASE("drawing query bonds", "[queries]") {
  SECTION("basics") {
    auto m = R"CTAB(
  Mrv2014 12072005332D          
  
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 14 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 3.7917 -2.96 0 0
M  V30 2 C 2.458 -3.73 0 0
M  V30 3 C 2.458 -5.27 0 0
M  V30 4 C 3.7917 -6.04 0 0
M  V30 5 C 5.1253 -5.27 0 0
M  V30 6 C 5.1253 -3.73 0 0
M  V30 7 C 6.459 -2.96 0 0
M  V30 8 C 3.7917 -7.58 0 0
M  V30 9 C 4.8806 -8.669 0 0
M  V30 10 C 4.482 -10.1565 0 0
M  V30 11 C 6.459 -6.04 0 0
M  V30 12 C 7.7927 -5.27 0 0
M  V30 13 C 9.1263 -6.0399 0 0
M  V30 14 C 9.1263 -7.5799 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 4 5
M  V30 3 1 1 6
M  V30 4 5 1 2
M  V30 5 6 5 6
M  V30 6 7 3 4
M  V30 7 8 6 7
M  V30 8 1 4 8
M  V30 9 1 8 9 TOPO=1
M  V30 10 1 9 10 TOPO=2
M  V30 11 1 5 11
M  V30 12 1 12 13
M  V30 13 2 11 12 TOPO=1
M  V30 14 2 13 14 TOPO=2
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testQueryBonds-1a.svg");
      outs << text;
      outs.close();
      check_file_hash("testQueryBonds-1a.svg");
    }
    {
      MolDraw2DSVG drawer(350, 300);
      m->getBondWithIdx(3)->setProp("bondNote", "S/D");
      m->getBondWithIdx(4)->setProp("bondNote", "S/A");
      m->getBondWithIdx(5)->setProp("bondNote", "D/A");
      m->getBondWithIdx(6)->setProp("bondNote", "Any");
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testQueryBonds-1b.svg");
      outs << text;
      outs.close();
      check_file_hash("testQueryBonds-1b.svg");
    }
    {
      MolDraw2DSVG drawer(350, 300);
      std::vector<int> highlightAtoms = {0, 1, 2, 3, 4, 5, 7, 8, 9};
      std::vector<int> highlightBonds = {0, 3, 2, 4, 1, 5, 8, 9};

      drawer.drawMolecule(*m, "", &highlightAtoms, &highlightBonds);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testQueryBonds-1c.svg");
      outs << text;
      outs.close();
      check_file_hash("testQueryBonds-1c.svg");
    }
  }
  SECTION("smaller drawing") {
    auto m = R"CTAB(
  Mrv2014 12012004302D          
  
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 26 29 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 3.7917 -2.96 0 0
M  V30 2 C 2.458 -3.73 0 0
M  V30 3 C 2.458 -5.27 0 0
M  V30 4 N 3.7917 -6.04 0 0
M  V30 5 N 5.1253 -5.27 0 0
M  V30 6 C 5.1253 -3.73 0 0
M  V30 7 C 6.459 -2.96 0 0
M  V30 8 C 3.7917 -7.58 0 0
M  V30 9 C 4.8806 -8.669 0 0
M  V30 10 C 4.482 -10.1565 0 0
M  V30 11 C 1.1243 -2.9599 0 0
M  V30 12 C -0.2093 -3.73 0 0
M  V30 13 C -0.2093 -5.27 0 0
M  V30 14 C 1.1243 -6.04 0 0
M  V30 15 C -0.2093 -0.6499 0 0
M  V30 16 C -1.543 -1.4199 0 0
M  V30 17 C -1.543 -2.9599 0 0
M  V30 18 C 1.1243 -1.4199 0 0
M  V30 19 C -2.8767 -0.6499 0 0
M  V30 20 C -4.2103 -1.4199 0 0
M  V30 21 C -4.2103 -2.9599 0 0
M  V30 22 C -2.8767 -3.73 0 0
M  V30 23 C -5.544 -3.7299 0 0
M  V30 24 C -6.8777 -2.9599 0 0
M  V30 25 C -8.2114 -3.7299 0 0
M  V30 26 C -9.5451 -2.9599 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 3
M  V30 2 1 4 5
M  V30 3 1 1 6
M  V30 4 5 1 2
M  V30 5 6 5 6
M  V30 6 7 3 4
M  V30 7 8 6 7
M  V30 8 1 4 8
M  V30 9 1 8 9 TOPO=1
M  V30 10 1 9 10 TOPO=2
M  V30 11 1 12 13
M  V30 12 1 13 14
M  V30 13 1 14 3
M  V30 14 1 11 2
M  V30 15 1 15 16
M  V30 16 1 16 17
M  V30 17 2 15 18
M  V30 18 1 11 18
M  V30 19 1 17 12
M  V30 20 2 12 11
M  V30 21 1 19 20
M  V30 22 2 20 21
M  V30 23 1 21 22
M  V30 24 2 19 16
M  V30 25 2 22 17
M  V30 26 1 21 23
M  V30 27 1 23 24
M  V30 28 1 24 25
M  V30 29 1 25 26
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(250, 200);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testQueryBonds-2.svg");
      outs << text;
      outs.close();
      check_file_hash("testQueryBonds-2.svg");
    }
  }
  SECTION("two linknodes") {
    auto m = R"CTAB(two linknodes
  Mrv2014 07072016412D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.25 12.1847 0 0
M  V30 2 C 6.9164 12.9547 0 0
M  V30 3 C 7.2366 14.4611 0 0
M  V30 4 C 8.7681 14.622 0 0
M  V30 5 C 9.3945 13.2151 0 0
M  V30 6 O 8.25 10.6447 0 0
M  V30 7 F 9.5382 15.9557 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 5
M  V30 4 1 1 5
M  V30 5 1 3 4
M  V30 6 1 1 6
M  V30 7 1 4 7
M  V30 END BOND
M  V30 LINKNODE 1 3 2 1 2 1 5
M  V30 LINKNODE 1 4 2 4 3 4 5
M  V30 END CTAB
M  END)CTAB"_ctab;
    std::vector<int> rotns = {0, 30, 60, 90, 120, 150, 180};
    for (auto rotn : rotns) {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().rotate = (double)rotn;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::string filename(
          (boost::format("testLinkNodes-2-%d.svg") % rotn).str());
      std::ofstream outs(filename);
      outs << text;
      outs.close();
      check_file_hash(filename);
    }
  }
}

TEST_CASE("molecule annotations", "[extra]") {
  int panelHeight = -1;
  int panelWidth = -1;
  bool noFreeType = false;
  SECTION("basics") {
    auto m = "NCC(=O)O"_smiles;
    MolDraw2DSVG drawer(350, 300, panelHeight, panelWidth, noFreeType);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    m->setProp(common_properties::molNote, "molecule note");
    drawer.drawMolecule(*m, "with note");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testMolAnnotations-1.svg");
    outs << text;
    outs.close();
    check_file_hash("testMolAnnotations-1.svg");
    CHECK(text.find("class='note'") != std::string::npos);
  }
  SECTION("chiral flag") {
    auto m = R"CTAB(
  Mrv2014 12152012512D          
 
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -0.6317 0.6787 0 0 CFG=2
M  V30 2 C -1.7207 1.7677 0 0
M  V30 3 C 0.4571 1.7677 0 0
M  V30 4 C -0.6317 2.8566 0 0 CFG=1
M  V30 5 C 0.1729 4.1698 0 0
M  V30 6 N -0.5619 5.5231 0 0
M  V30 7 C -1.4364 4.1698 0 0
M  V30 8 C -0.6316 -0.8613 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=3
M  V30 2 1 1 3
M  V30 3 1 4 3
M  V30 4 1 4 2
M  V30 5 1 4 5
M  V30 6 1 5 6
M  V30 7 1 4 7 CFG=1
M  V30 8 1 1 8
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    {
      MolDraw2DSVG drawer(350, 300, panelHeight, panelWidth, noFreeType);
      drawer.drawMolecule(*m, "chiral flag set, option disabled");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testMolAnnotations-2a.svg");
      outs << text;
      outs.close();
      check_file_hash("testMolAnnotations-2a.svg");
      CHECK(text.find("class='note'") == std::string::npos);
    }
    {
      MolDraw2DSVG drawer(350, 300, panelHeight, panelWidth, noFreeType);
      drawer.drawOptions().includeChiralFlagLabel = true;
      drawer.drawMolecule(*m, "chiral flag set, option enabled");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testMolAnnotations-2b.svg");
      outs << text;
      outs.close();
      check_file_hash("testMolAnnotations-2b.svg");
      CHECK(text.find("class='note'") != std::string::npos);
    }
    {
      MolDraw2DSVG drawer(350, 300, panelHeight, panelWidth, noFreeType);
      drawer.drawOptions().includeChiralFlagLabel = true;
      m->clearProp(common_properties::_MolFileChiralFlag);
      drawer.drawMolecule(*m, "chiral flag not set, option enabled");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testMolAnnotations-2c.svg");
      outs << text;
      outs.close();
      check_file_hash("testMolAnnotations-2c.svg");
      CHECK(text.find("class='note'") == std::string::npos);
    }
  }
  SECTION("simplified stereo 1") {
    {
      auto m = "C[C@H](F)[C@@H](F)[C@@H](C)Cl |o1:3,5,1|"_smiles;
      MolDraw2DSVG drawer(350, 300, panelHeight, panelWidth, noFreeType);
      MolDraw2DUtils::prepareMolForDrawing(*m);
      drawer.drawOptions().addStereoAnnotation = true;
      drawer.drawMolecule(*m, "enhanced no flag");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testMolAnnotations-3a.svg");
      outs << text;
      outs.close();
      check_file_hash("testMolAnnotations-3a.svg");
    }
    {
      auto m = "C[C@H](F)[C@@H](F)[C@@H](C)Cl |o1:3,5,1|"_smiles;
      MolDraw2DSVG drawer(350, 300, panelHeight, panelWidth, noFreeType);
      MolDraw2DUtils::prepareMolForDrawing(*m);
      drawer.drawOptions().addStereoAnnotation = true;
      drawer.drawOptions().simplifiedStereoGroupLabel = true;
      drawer.drawMolecule(*m, "enhanced with flag");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testMolAnnotations-3b.svg");
      outs << text;
      outs.close();
      check_file_hash("testMolAnnotations-3b.svg");
    }
    {
      auto m = "C[C@H](F)[C@@H](F)[C@@H](C)Cl |&1:3,5,1|"_smiles;
      MolDraw2DSVG drawer(350, 300, panelHeight, panelWidth, noFreeType);
      MolDraw2DUtils::prepareMolForDrawing(*m);
      drawer.drawOptions().addStereoAnnotation = true;
      drawer.drawOptions().simplifiedStereoGroupLabel = true;
      drawer.drawMolecule(*m, "enhanced & with flag");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testMolAnnotations-3c.svg");
      outs << text;
      outs.close();
      check_file_hash("testMolAnnotations-3c.svg");
    }
  }
  SECTION("simplified stereo 2") {
    auto m = "C[C@H](F)[C@@H](F)[C@@H](C)Cl |o1:3,5,o2:1|"_smiles;
    MolDraw2DSVG drawer(350, 300, panelHeight, panelWidth, noFreeType);
    drawer.drawOptions().addStereoAnnotation = true;
    drawer.drawOptions().simplifiedStereoGroupLabel = true;
    MolDraw2DUtils::prepareMolForDrawing(*m);
    drawer.drawMolecule(*m, "multi-groups");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testMolAnnotations-3d.svg");
    outs << text;
    outs.close();
    check_file_hash("testMolAnnotations-3d.svg");
  }
  SECTION("label placement") {
    auto m = R"CTAB(
  Mrv2014 12162004412D          
 
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 16 15 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -9.2917 3.5833 0 0
M  V30 2 C -7.958 4.3533 0 0 CFG=2
M  V30 3 C -6.6243 3.5833 0 0 CFG=1
M  V30 4 C -5.2906 4.3533 0 0 CFG=2
M  V30 5 Cl -7.958 5.8933 0 0
M  V30 6 F -6.6243 2.0433 0 0
M  V30 7 F -3.957 3.5833 0 0
M  V30 8 C -5.2906 5.8933 0 0
M  V30 9 C -3.957 6.6633 0 0
M  V30 10 C -3.957 8.2033 0 0
M  V30 11 C -2.6233 8.9733 0 0
M  V30 12 C -2.6233 5.8933 0 0
M  V30 13 C -5.2906 8.9733 0 0
M  V30 14 C -2.6233 10.5133 0 0
M  V30 15 C -1.2896 8.2033 0 0
M  V30 16 C -1.2896 6.6633 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 2 5 CFG=1
M  V30 5 1 3 6 CFG=1
M  V30 6 1 4 7 CFG=1
M  V30 7 1 4 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 1 9 12
M  V30 12 1 10 13
M  V30 13 1 11 14
M  V30 14 1 11 15
M  V30 15 1 12 16
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEREL1 ATOMS=(3 2 3 4)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)CTAB"_ctab;
    MolDraw2DSVG drawer(350, 300, panelHeight, panelWidth, noFreeType);
    drawer.drawOptions().addStereoAnnotation = true;
    drawer.drawOptions().simplifiedStereoGroupLabel = true;
    drawer.drawMolecule(*m, "label crowding");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testMolAnnotations-4a.svg");
    outs << text;
    outs.close();
    check_file_hash("testMolAnnotations-4a.svg");
  }
}

TEST_CASE("draw link nodes", "[extras]") {
  SECTION("one linknode") {
    auto m = R"CTAB(one linknode
  Mrv2007 06222005102D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.25 12.1847 0 0
M  V30 2 C 6.9164 12.9547 0 0
M  V30 3 C 6.9164 14.4947 0 0
M  V30 4 C 9.5836 14.4947 0 0
M  V30 5 C 9.5836 12.9547 0 0
M  V30 6 O 8.25 10.6447 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 4 5
M  V30 4 1 1 5
M  V30 5 1 3 4
M  V30 6 1 1 6
M  V30 END BOND
M  V30 LINKNODE 1 4 2 1 2 1 5
M  V30 END CTAB
M  END)CTAB"_ctab;
    std::vector<int> rotns = {0, 30, 60, 90, 120, 150, 180};
    for (auto rotn : rotns) {
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().rotate = (double)rotn;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs(
          (boost::format("testLinkNodes-1-%d.svg") % rotn).str());
      outs << text;
      outs.close();
      check_file_hash((boost::format("testLinkNodes-1-%d.svg") % rotn).str());
    }
  }
}

TEST_CASE("Github #3744: Double bonds incorrectly drawn outside the ring",
          "[drawing]") {
  SECTION("SVG") {
    ROMOL_SPTR m1(MolBlockToMol(R"CTAB(
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0684   -1.2135    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4949   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4949    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0684    1.2135    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8133    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3133   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  1  0
  5  1  1  0
M  END)CTAB"));
    REQUIRE(m1);
    MolDraw2DSVG drawer(400, 300);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub3744.svg");
    outs << text;
    outs.close();
    check_file_hash("testGithub3744.svg");
    std::vector<std::string> bond0;
    std::vector<std::string> bond2;
    std::istringstream ss(text);
    std::string line;
    while (std::getline(ss, line)) {
      if (line.find("bond-0") != std::string::npos) {
        bond0.push_back(line);
      } else if (line.find("bond-2") != std::string::npos) {
        bond2.push_back(line);
      }
    }
    CHECK(bond0.size() == 2);
    CHECK(bond2.size() == 2);
    std::regex regex(
        "^.*d='M\\s+(\\d+\\.\\d+),(\\d+\\.\\d+)\\s+L\\s+(\\d+\\.\\d+),(\\d+\\."
        "\\d+)'.*$");
    std::smatch bond0OuterMatch;
    REQUIRE(std::regex_match(bond0[0], bond0OuterMatch, regex));
    REQUIRE(bond0OuterMatch.size() == 5);
    std::smatch bond0InnerMatch;
    REQUIRE(std::regex_match(bond0[1], bond0InnerMatch, regex));
    REQUIRE(bond0InnerMatch.size() == 5);
    std::smatch bond2OuterMatch;
    REQUIRE(std::regex_match(bond2[0], bond2OuterMatch, regex));
    REQUIRE(bond2OuterMatch.size() == 5);
    std::smatch bond2InnerMatch;
    REQUIRE(std::regex_match(bond2[1], bond2InnerMatch, regex));
    REQUIRE(bond2InnerMatch.size() == 5);
    RDGeom::Point2D bond0InnerCtd(
        RDGeom::Point2D(std::stof(bond0InnerMatch[1]),
                        std::stof(bond0InnerMatch[2])) +
        RDGeom::Point2D(std::stof(bond0InnerMatch[3]),
                        std::stof(bond0InnerMatch[4])) /
            2.0);
    RDGeom::Point2D bond0OuterCtd(
        RDGeom::Point2D(std::stof(bond0OuterMatch[1]),
                        std::stof(bond0OuterMatch[2])) +
        RDGeom::Point2D(std::stof(bond0OuterMatch[3]),
                        std::stof(bond0OuterMatch[4])) /
            2.0);
    RDGeom::Point2D bond2InnerCtd(
        RDGeom::Point2D(std::stof(bond2InnerMatch[1]),
                        std::stof(bond2InnerMatch[2])) +
        RDGeom::Point2D(std::stof(bond2InnerMatch[3]),
                        std::stof(bond2InnerMatch[4])) /
            2.0);
    RDGeom::Point2D bond2OuterCtd(
        RDGeom::Point2D(std::stof(bond2OuterMatch[1]),
                        std::stof(bond2OuterMatch[2])) +
        RDGeom::Point2D(std::stof(bond2OuterMatch[3]),
                        std::stof(bond2OuterMatch[4])) /
            2.0);
    // we look at the two double bonds of pyrrole
    // we check that the ratio between the distance of the centroids of the
    // outer bonds and the distance of the centroids of the inner bonds is at
    // least 1.275, otherwise the inner bonds are not actually inside the ring.
    float outerBondsDistance = (bond0OuterCtd - bond2OuterCtd).length();
    float innerBondsDistance = (bond0InnerCtd - bond2InnerCtd).length();
    CHECK(outerBondsDistance / innerBondsDistance > 1.275f);
  }
}

TEST_CASE("draw atom list queries", "[extras]") {
  SECTION("atom list") {
    auto m = R"CTAB(
  Mrv2102 02112115002D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 [N,O,S] 9.2083 12.8058 0 0
M  V30 2 C 8.4383 11.4721 0 0
M  V30 3 C 9.9783 11.4721 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 1
M  V30 3 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    MolDraw2DSVG drawer(350, 300);
    drawer.drawMolecule(*m, "atom list");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testAtomLists-1.svg");
    outs << text;
    outs.close();
    check_file_hash("testAtomLists-1.svg");
  }

  SECTION("NOT atom list") {
    auto m = R"CTAB(
  Mrv2102 02112115032D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 "NOT [N,O,S]" 9.2083 12.8058 0 0
M  V30 2 C 8.4383 11.4721 0 0
M  V30 3 C 9.9783 11.4721 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 1
M  V30 3 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);
    MolDraw2DSVG drawer(350, 300);
    drawer.drawMolecule(*m, "NOT atom list");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testAtomLists-2.svg");
    outs << text;
    outs.close();
    check_file_hash("testAtomLists-2.svg");
  }
}

TEST_CASE("test the options that toggle isotope labels", "[drawing]") {
  SECTION("test all permutations") {
    auto m = "[1*]c1cc([2*])c([3*])c[14c]1"_smiles;
    REQUIRE(m);
    std::regex regex(R"regex(<text\s+.*>\d</text>)regex");
    std::smatch match;
    std::string line;
    {
      MolDraw2DSVG drawer(300, 300, -1, -1, true);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string textIsoDummyIso = drawer.getDrawingText();
      std::ofstream outs("testIsoDummyIso.svg");
      outs << textIsoDummyIso;
      outs.close();
      check_file_hash("testIsoDummyIso.svg");
      size_t nIsoDummyIso = std::distance(
          std::sregex_token_iterator(textIsoDummyIso.begin(),
                                     textIsoDummyIso.end(), regex),
          std::sregex_token_iterator());
      CHECK(nIsoDummyIso == 5);
    }
    {
      MolDraw2DSVG drawer(300, 300, -1, -1, true);
      drawer.drawOptions().isotopeLabels = false;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string textNoIsoDummyIso = drawer.getDrawingText();
      std::ofstream outs("testNoIsoDummyIso.svg");
      outs << textNoIsoDummyIso;
      outs.close();
      check_file_hash("testNoIsoDummyIso.svg");
      size_t nNoIsoDummyIso = std::distance(
          std::sregex_token_iterator(textNoIsoDummyIso.begin(),
                                     textNoIsoDummyIso.end(), regex, 1),
          std::sregex_token_iterator());
      CHECK(nNoIsoDummyIso == 3);
    }
    {
      MolDraw2DSVG drawer(300, 300, -1, -1, true);
      drawer.drawOptions().dummyIsotopeLabels = false;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string textIsoNoDummyIso = drawer.getDrawingText();
      std::ofstream outs("testIsoNoDummyIso.svg");
      outs << textIsoNoDummyIso;
      outs.close();
      check_file_hash("testIsoNoDummyIso.svg");
      size_t nIsoNoDummyIso = std::distance(
          std::sregex_token_iterator(textIsoNoDummyIso.begin(),
                                     textIsoNoDummyIso.end(), regex, 1),
          std::sregex_token_iterator());
      CHECK(nIsoNoDummyIso == 2);
    }
    {
      MolDraw2DSVG drawer(300, 300, -1, -1, true);
      drawer.drawOptions().isotopeLabels = false;
      drawer.drawOptions().dummyIsotopeLabels = false;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string textNoIsoNoDummyIso = drawer.getDrawingText();
      std::ofstream outs("testNoIsoNoDummyIso.svg");
      outs << textNoIsoNoDummyIso;
      outs.close();
      check_file_hash("testNoIsoNoDummyIso.svg");
      size_t nNoIsoNoDummyIso = std::distance(
          std::sregex_token_iterator(textNoIsoNoDummyIso.begin(),
                                     textNoIsoNoDummyIso.end(), regex, 1),
          std::sregex_token_iterator());
      CHECK(nNoIsoNoDummyIso == 0);
    }
  }
  SECTION("test that D/T show up even if isotope labels are hidden") {
    auto m = "C([1H])([2H])([3H])[H]"_smiles;
    std::regex regex(R"regex(<text\s+.*>[DT]</text>)regex");
    std::smatch match;
    REQUIRE(m);
    std::string line;
    MolDraw2DSVG drawer(300, 300, -1, -1, true);
    drawer.drawOptions().isotopeLabels = false;
    drawer.drawOptions().dummyIsotopeLabels = false;
    drawer.drawOptions().atomLabelDeuteriumTritium = true;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string textDeuteriumTritium = drawer.getDrawingText();
    std::ofstream outs("testDeuteriumTritium.svg");
    outs << textDeuteriumTritium;
    outs.close();
    check_file_hash("testDeuteriumTritium.svg");
    size_t nDeuteriumTritium = std::distance(
        std::sregex_token_iterator(textDeuteriumTritium.begin(),
                                   textDeuteriumTritium.end(), regex, 1),
        std::sregex_token_iterator());
    CHECK(nDeuteriumTritium == 2);
  }
}

TEST_CASE("draw hydrogen bonds", "[drawing]") {
  SECTION("basics") {
    auto m = R"CTAB(
  Mrv2014 03022114422D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -5.4583 -0.125 0 0
M  V30 2 C -4.1247 0.645 0 0
M  V30 3 C -2.791 -0.125 0 0
M  V30 4 C -1.4573 0.645 0 0
M  V30 5 O -2.791 -1.665 0 0
M  V30 6 C -6.792 0.645 0 0
M  V30 7 O -5.4583 -1.665 0 0
M  V30 8 H -4.1247 -2.435 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 2 3 5
M  V30 5 1 1 6
M  V30 6 1 1 7
M  V30 7 1 7 8
M  V30 8 10 5 8
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m);

    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::ofstream outs("testHydrogenBonds1.svg");
    outs << drawer.getDrawingText();
    outs.close();
    check_file_hash("testHydrogenBonds1.svg");
  }
  SECTION("from CXSMILES") {
    auto m = "CC1O[H]O=C(C)C1 |H:4.3|"_smiles;
    REQUIRE(m);

    MolDraw2DSVG drawer(300, 300);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::ofstream outs("testHydrogenBonds2.svg");
    outs << drawer.getDrawingText();
    outs.close();
    check_file_hash("testHydrogenBonds2.svg");
  }
}

TEST_CASE("github #3912: cannot draw atom lists from SMARTS", "[query][bug]") {
  SECTION("original") {
    auto m = "C-[N,O]"_smarts;
    REQUIRE(m);
    int panelWidth = -1;
    int panelHeight = -1;
    bool noFreeType = true;
    MolDraw2DSVG drawer(300, 300, panelWidth, panelHeight, noFreeType);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::ofstream outs("testGithub3912.1.svg");
    auto txt = drawer.getDrawingText();
    outs << txt;
    outs.close();
    check_file_hash("testGithub3912.1.svg");
    CHECK(txt.find(">N<") != std::string::npos);
    CHECK(txt.find(">O<") != std::string::npos);
    CHECK(txt.find(">!<") == std::string::npos);
  }
  SECTION("negated") {
    auto m = "C-[N,O]"_smarts;
    REQUIRE(m);
    REQUIRE(m->getAtomWithIdx(1)->hasQuery());
    m->getAtomWithIdx(1)->getQuery()->setNegation(true);
    int panelWidth = -1;
    int panelHeight = -1;
    bool noFreeType = true;
    MolDraw2DSVG drawer(300, 300, panelWidth, panelHeight, noFreeType);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::ofstream outs("testGithub3912.2.svg");
    auto txt = drawer.getDrawingText();
    outs << txt;
    outs.close();
    check_file_hash("testGithub3912.2.svg");
    CHECK(txt.find(">N<") != std::string::npos);
    CHECK(txt.find(">O<") != std::string::npos);
    CHECK(txt.find(">!<") != std::string::npos);
  }
}

TEST_CASE("github #2976: kekulizing reactions when drawing", "[reactions]") {
  SECTION("basics") {
    bool asSmiles = true;
    std::unique_ptr<ChemicalReaction> rxn{
        RxnSmartsToChemicalReaction("c1ccccc1>>c1ncccc1", nullptr, asSmiles)};
    MolDraw2DSVG drawer(450, 200);
    drawer.drawReaction(*rxn);
    drawer.finishDrawing();
    std::ofstream outs("testGithub2976.svg");
    auto txt = drawer.getDrawingText();
    outs << txt;
    outs.close();
    check_file_hash("testGithub2976.svg");
  }
}

TEST_CASE("preserve Reaction coordinates", "[reactions]") {
  SECTION("basics") {
    std::string data = R"RXN($RXN

  Mrv16822    031301211645

  2  2  1
$MOL

  Mrv1682203132116452D          

  3  2  0  0  0  0            999 V2000
   -4.3304    2.5893    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3304    1.7643    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5054    1.7643    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$MOL

  Mrv1682203132116452D          

  2  1  0  0  0  0            999 V2000
   -2.1652    2.6339    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1652    1.8089    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
$MOL

  Mrv1682203132116452D          

  3  2  0  0  0  0            999 V2000
    3.6109    1.9512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7859    1.9512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.7859    2.7762    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  3  2  1  0  0  0  0
M  END
$MOL

  Mrv1682203132116452D          

  2  1  0  0  0  0            999 V2000
    4.9511    1.9959    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9511    2.8209    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
M  END
$MOL

  Mrv1682203132116452D          

  2  1  0  0  0  0            999 V2000
   -0.3571    2.7232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4003    3.5471    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
)RXN";
    std::unique_ptr<ChemicalReaction> rxn{RxnBlockToChemicalReaction(data)};
    MolDraw2DSVG drawer(450, 200);
    drawer.drawReaction(*rxn);
    drawer.finishDrawing();
    std::ofstream outs("testReactionCoords.svg");
    auto txt = drawer.getDrawingText();
    outs << txt;
    outs.close();
    check_file_hash("testReactionCoords.svg");

    // the reaction is drawn with some bonds vertical, make sure they remain
    // vertical
    {
      std::regex regex("class='bond-0.*? d='M (\\d+\\.\\d+).* L (\\d+\\.\\d+)");
      std::smatch bondMatch;
      REQUIRE(std::regex_search(txt, bondMatch, regex));
      REQUIRE(bondMatch.size() == 3);  // match both halves of the bond
      CHECK(bondMatch[1].str() == bondMatch[2].str());
    }
    {
      std::regex regex("class='bond-2.*? d='M (\\d+\\.\\d+).* L (\\d+\\.\\d+)");
      std::smatch bondMatch;
      REQUIRE(std::regex_search(txt, bondMatch, regex));
      REQUIRE(bondMatch.size() == 3);  // match both halves of the bond
      CHECK(bondMatch[1].str() == bondMatch[2].str());
    }
    {
      std::regex regex("class='bond-4.*? d='M (\\d+\\.\\d+).* L (\\d+\\.\\d+)");
      std::smatch bondMatch;
      REQUIRE(std::regex_search(txt, bondMatch, regex));
      REQUIRE(bondMatch.size() == 3);  // match both halves of the bond
      CHECK(bondMatch[1].str() == bondMatch[2].str());
    }
  }
}
TEST_CASE("support annotation colors", "[drawing]") {
  SECTION("basics") {
    auto m = "CCCO"_smiles;
    REQUIRE(m);
    int panelWidth = -1;
    int panelHeight = -1;
    bool noFreeType = true;
    MolDraw2DSVG drawer(300, 300, panelWidth, panelHeight, noFreeType);
    drawer.drawOptions().annotationColour = DrawColour{0, 0, 1, 1};
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawMolecule(*m, "blue annotations");
    drawer.finishDrawing();
    std::ofstream outs("testAnnotationColors.svg");
    auto txt = drawer.getDrawingText();
    outs << txt;
    outs.close();
    check_file_hash("testAnnotationColors.svg");
    CHECK(txt.find("fill:#0000FF' >2<") != std::string::npos);
  }
}

TEST_CASE("Github #4238: prepareMolForDrawing and wavy bonds") {
  {
    auto mol = "CC=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(1)->setStereoAtoms(0, 3);
    mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOANY);
    bool kekulize = true;
    bool addChiralHs = true;
    bool wedgeBonds = true;
    bool forceCoords = true;
    bool wavyBonds = false;
    MolDraw2DUtils::prepareMolForDrawing(*mol, kekulize, addChiralHs,
                                         wedgeBonds, forceCoords, wavyBonds);
    CHECK(mol->getBondWithIdx(0)->getBondDir() == Bond::BondDir::NONE);
    CHECK(mol->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREOANY);

    RWMol mol2(*mol);
    wavyBonds = true;
    MolDraw2DUtils::prepareMolForDrawing(mol2, kekulize, addChiralHs,
                                         wedgeBonds, forceCoords, wavyBonds);
    CHECK(mol2.getBondWithIdx(0)->getBondDir() == Bond::BondDir::UNKNOWN);
    CHECK(mol2.getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREONONE);

    MOL_PTR_VECT ms{mol.get(), &mol2};
    {
      MolDraw2DSVG drawer(500, 200, 250, 200);
      // drawer.drawOptions().prepareMolsBeforeDrawing = false;
      std::vector<std::string> legends = {"before", "after"};
      drawer.drawMolecules(ms, &legends);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub4238_1.svg");
      outs << text;
      outs.flush();
      check_file_hash("testGithub4238_1.svg");
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(500, 200, 250, 200);
      std::vector<std::string> legends = {"before", "after"};
      drawer.drawMolecules(ms, &legends);
      drawer.finishDrawing();
      drawer.writeDrawingText("testGithub4238_1.png");
      check_file_hash("testGithub4238_1.png");
    }
#endif
  }
}

TEST_CASE("Github #4323: support providing RGBA colors") {
  auto mol = "CCCO"_smiles;
  REQUIRE(mol);
#ifdef RDK_BUILD_FREETYPE_SUPPORT
  SECTION("with alpha") {
    MolDraw2DSVG drawer(200, 150);
    drawer.drawOptions().legendColour = DrawColour(1, 0, 1, 0.3);
    drawer.drawOptions().backgroundColour = DrawColour(0.5, 0.5, 0.5, 0.3);
    drawer.drawMolecule(*mol, "partially transparent legend/background");
    drawer.finishDrawing();

    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub4323_1.svg");
    outs << text;
    outs.flush();
    // background
    CHECK(text.find("fill:#7F7F7F4C;") != std::string::npos);
    CHECK(text.find("fill:#7F7F7F;") == std::string::npos);
    // legend
    CHECK(text.find("fill='#FF00FF4C'") != std::string::npos);
    CHECK(text.find("fill='#FF00FF'") == std::string::npos);
    check_file_hash("testGithub4323_1.svg");
  }
  SECTION("without alpha") {
    MolDraw2DSVG drawer(200, 150);
    drawer.drawOptions().legendColour = DrawColour(1, 0, 1);
    drawer.drawOptions().backgroundColour = DrawColour(0.5, 0.5, 0.5);
    drawer.drawMolecule(*mol, "no transparency");
    drawer.finishDrawing();

    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub4323_2.svg");
    outs << text;
    outs.flush();
    // background
    CHECK(text.find("fill:#7F7F7F4C;") == std::string::npos);
    CHECK(text.find("fill:#7F7F7F;") != std::string::npos);
    // legend
    CHECK(text.find("fill='#FF00FF4C'") == std::string::npos);
    CHECK(text.find("fill='#FF00FF'") != std::string::npos);
    check_file_hash("testGithub4323_2.svg");
  }
#endif
  SECTION("no FT with alpha") {
    MolDraw2DSVG drawer(200, 150, -1, -1, NO_FREETYPE);
    drawer.drawOptions().legendColour = DrawColour(1, 0, 1, 0.3);
    drawer.drawOptions().backgroundColour = DrawColour(0.5, 0.5, 0.5, 0.3);
    drawer.drawMolecule(*mol, "partially transparent legend/background");
    drawer.finishDrawing();

    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub4323_3.svg");
    outs << text;
    outs.flush();
    // background
    CHECK(text.find("fill:#7F7F7F4C;") != std::string::npos);
    CHECK(text.find("fill:#7F7F7F;") == std::string::npos);
    // legend
    CHECK(text.find("fill:#FF00FF4C'") != std::string::npos);
    CHECK(text.find("fill:#FF00FF'") == std::string::npos);
    check_file_hash("testGithub4323_3.svg");
  }
  SECTION("no FT without alpha") {
    MolDraw2DSVG drawer(200, 150, -1, -1, NO_FREETYPE);
    drawer.drawOptions().legendColour = DrawColour(1, 0, 1);
    drawer.drawOptions().backgroundColour = DrawColour(0.5, 0.5, 0.5);
    drawer.drawMolecule(*mol, "no transparency");
    drawer.finishDrawing();

    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub4323_4.svg");
    outs << text;
    outs.flush();
    // background
    CHECK(text.find("fill:#7F7F7F4C;") == std::string::npos);
    CHECK(text.find("fill:#7F7F7F;") != std::string::npos);
    // legend
    CHECK(text.find("fill:#FF00FF4C'") == std::string::npos);
    CHECK(text.find("fill:#FF00FF'") != std::string::npos);
    check_file_hash("testGithub4323_4.svg");
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
#ifdef RDK_BUILD_FREETYPE_SUPPORT
  SECTION("Cairo with alpha") {
    MolDraw2DCairo drawer(200, 150);
    drawer.drawOptions().legendColour = DrawColour(1, 0, 1, 0.3);
    drawer.drawOptions().backgroundColour = DrawColour(0.5, 0.5, 0.5, 0.3);
    drawer.drawMolecule(*mol, "partially transparent legend/background");
    drawer.finishDrawing();
    drawer.writeDrawingText("testGithub4323_1.png");
    check_file_hash("testGithub4323_1.png");
  }
#endif
  SECTION("No FT Cairo with alpha") {
    MolDraw2DCairo drawer(200, 150, -1, -1, NO_FREETYPE);
    drawer.drawOptions().legendColour = DrawColour(1, 0, 1, 0.3);
    drawer.drawOptions().backgroundColour = DrawColour(0.5, 0.5, 0.5, 0.3);
    drawer.drawMolecule(*mol, "partially transparent legend/background");
    drawer.finishDrawing();
    drawer.writeDrawingText("testGithub4323_3.png");
    check_file_hash("testGithub4323_3.png");
  }
#endif
}

TEST_CASE(
    "Github #4508: SubstanceGroup labels sometimes overlap with atoms in image "
    "generation") {
  SECTION("Basics") {
    auto mol = R"CTAB(
  Mrv2114 09132120172D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 1 0 1
M  V30 BEGIN ATOM
M  V30 1 C -0.5878 0.8085 0 0
M  V30 2 C -1.9434 0.078 0 0
M  V30 3 C -1.9884 -1.4614 0 0
M  V30 4 C -0.6778 -2.2702 0 0
M  V30 5 C 0.6778 -1.5394 0 0
M  V30 6 C 0.7228 -0.0001 0 0
M  V30 7 N -0.5428 2.3478 0 0
M  V30 8 O 1.9884 -2.3479 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 1 7
M  V30 8 1 5 8
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 7) FIELDNAME=UV FIELDINFO=nm -
M  V30 FIELDDISP="    0.0000    0.0000    DRU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=340
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol);

    {
      MolDraw2DSVG drawer(300, 250);
      drawer.drawMolecule(*mol, "data label with DRU");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub4508_1.svg");
      outs << text;
      outs.flush();
      check_file_hash("testGithub4508_1.svg");
    }

    // remove the sgroup-atom atom... the SGroup will not be drawn
    auto &sgs = getSubstanceGroups(*mol);
    REQUIRE(sgs.size() == 1);
    sgs[0].setAtoms(std::vector<unsigned int>());
    {
      MolDraw2DSVG drawer(300, 250);
      drawer.drawMolecule(*mol, "no data label drawn");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub4508_1b.svg");
      outs << text;
      outs.flush();
      check_file_hash("testGithub4508_1b.svg");
    }
  }
  SECTION("Absolute") {
    auto mol = R"CTAB(
  Mrv2114 09132120172D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 1 0 1
M  V30 BEGIN ATOM
M  V30 1 C -0.5878 0.8085 0 0
M  V30 2 C -1.9434 0.078 0 0
M  V30 3 C -1.9884 -1.4614 0 0
M  V30 4 C -0.6778 -2.2702 0 0
M  V30 5 C 0.6778 -1.5394 0 0
M  V30 6 C 0.7228 -0.0001 0 0
M  V30 7 N -0.5428 2.3478 0 0
M  V30 8 O 1.9884 -2.3479 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 1 7
M  V30 8 1 5 8
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 7) FIELDNAME=UV FIELDINFO=nm -
M  V30 FIELDDISP="    0.0000    0.0000    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=340
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol);

    {
      MolDraw2DSVG drawer(300, 250);
      drawer.drawMolecule(*mol, "data label with DAU\n(expect odd placement)");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub4508_2.svg");
      outs << text;
      outs.flush();
      check_file_hash("testGithub4508_2.svg");
    }

    // remove the sgroup-atom atom... the SGroup will still be drawn
    auto &sgs = getSubstanceGroups(*mol);
    REQUIRE(sgs.size() == 1);
    sgs[0].setAtoms(std::vector<unsigned int>());
    {
      MolDraw2DSVG drawer(300, 250);
      drawer.drawMolecule(*mol,
                          "DAU, no associated atom\n(expect odd placement)");
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub4508_2b.svg");
      outs << text;
      outs.flush();
      check_file_hash("testGithub4508_2b.svg");
    }
  }
}

TEST_CASE("Github #4538 drawMolecules crash") {
  auto m = "CCc1ccccc1"_smiles;
  REQUIRE(m);
  RDDepict::compute2DCoords(*m);
  ROMol m1(*m);
  ROMol m2(*m);
  std::vector<ROMol *> mols{&m1, &m2};
  SECTION("basics") {
    MolDraw2DSVG drawer(500, 200, 250, 200);
    drawer.drawOptions().prepareMolsBeforeDrawing = false;
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testGithub4538.svg");
    outs << text;
    outs.flush();
    check_file_hash("testGithub4538.svg");
  }
}

TEST_CASE("dark mode mol drawing") {
  SECTION("Basics") {
    auto m =
        "CS(=O)(=O)COC(=N)c1cc(Cl)cnc1[NH3+] |SgD:7:note:some extra text:=:::|"_smiles;
    REQUIRE(m);
    MolDraw2DSVG drawer(350, 300);
    setDarkMode(drawer);
    drawer.drawMolecule(*m, "dark mode!");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testDarkMode.1.svg");
    outs << text;
    outs.flush();
    check_file_hash("testDarkMode.1.svg");
  }
}
TEST_CASE("monochrome mol drawing") {
  SECTION("Basics") {
    auto m =
        "CS(=O)(=O)COC(=N)c1cc(Cl)cnc1[NH3+] |SgD:7:note:some extra text:=:::|"_smiles;
    REQUIRE(m);
    MolDraw2DSVG drawer(350, 300);
    setMonochromeMode(drawer, DrawColour{0.1, 0.1, 0.6},
                      DrawColour{0.75, 0.75, 0.75});
    drawer.drawMolecule(*m, "monochrome");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testMonochrome.1.svg");
    outs << text;
    outs.flush();
    check_file_hash("testMonochrome.1.svg");
  }
  SECTION("Basics inverted") {
    auto m =
        "CS(=O)(=O)COC(=N)c1cc(Cl)cnc1[NH3+] |SgD:7:note:some extra text:=:::|"_smiles;
    REQUIRE(m);
    MolDraw2DSVG drawer(350, 300);
    setMonochromeMode(drawer, DrawColour{0.75, 0.75, 0.75},
                      DrawColour{0.1, 0.1, 0.6});
    drawer.drawMolecule(*m, "monochrome");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testMonochrome.2.svg");
    outs << text;
    outs.flush();
    check_file_hash("testMonochrome.2.svg");
  }
}
TEST_CASE("other palettes") {
  auto m =
      "CS(=O)(=O)COC(=N)c1c(I)c(Cl)c(Br)nc1[NH2+]CP(=O) |SgD:7:note:some extra text:=:::|"_smiles;
  REQUIRE(m);
  SECTION("Avalon") {
    MolDraw2DSVG drawer(350, 300);
    assignAvalonPalette(drawer.drawOptions().atomColourPalette);
    drawer.drawMolecule(*m, "Avalon");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testAvalon.1.svg");
    outs << text;
    outs.flush();
    check_file_hash("testAvalon.1.svg");
  }
  SECTION("CDK") {
    MolDraw2DSVG drawer(350, 300);
    assignCDKPalette(drawer.drawOptions().atomColourPalette);
    drawer.drawMolecule(*m, "CDK");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testCDK.1.svg");
    outs << text;
    outs.flush();
    check_file_hash("testCDK.1.svg");
  }
}

TEST_CASE("SDD record parsing") {
  auto mol = R"CTAB(
  Mrv2008 11122110292D

  6  6  0  0  0  0            999 V2000
    9.3527    2.5661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6382    2.1536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6382    1.3286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.3527    0.9161    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0671    1.3286    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0671    2.1536    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  1  6  2  0  0  0  0
M  STY  1   1 DAT
M  SLB  1   1   1
M  SAL   1  1   1
M  SDT   1 NAME
M  SDD   1 -2345.1234-2345.1234    DR    ALL  1       0
M  SED   1 Hello World
M  END
)CTAB"_ctab;
  // SDD record has format
  // M  SDD sss xxxxx.xxxxyyyyy.yyyy eeefgh i jjjkkk ll m noo
  MolDraw2DSVG drawer(350, 300, -1, -1, 1);
  drawer.drawMolecule(*mol);
  drawer.finishDrawing();
  auto text = drawer.getDrawingText();
  std::string name("Hello World");
  for (auto &c : name) {
    std::stringstream ss;
    ss << " >" << c << "</text>";
    auto pos = text.find(ss.str());
    CHECK(pos != std::string::npos);
  }
}

TEST_CASE("Github #4519 bad placement of datafield labels") {
  auto mol1 = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 C 2.598076 -0.000000 0.000000 0
M  V30 4 C 1.299038 2.250000 0.000000 0
M  V30 5 C 2.598076 3.000000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 2 4
M  V30 4 2 4 5
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(5 2 4 5 3 1) FIELDNAME="Lambda Max" FIELDINFO=nm -
M  V30 FIELDDATA="2222"
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
  REQUIRE(mol1);

  auto mol2 = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 8 1 0 0
M  V30 BEGIN ATOM
M  V30 1 N 3.000000 0.000000 0.000000 0
M  V30 2 C 1.500000 0.000000 0.000000 0
M  V30 3 C 0.750000 -1.299038 0.000000 0
M  V30 4 C -0.750000 -1.299038 0.000000 0
M  V30 5 C -1.500000 0.000000 0.000000 0
M  V30 6 C -0.750000 1.299038 0.000000 0
M  V30 7 O -1.500000 2.598076 0.000000 0
M  V30 8 C 0.750000 1.299038 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 2 6 8
M  V30 8 1 8 2
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 1) FIELDNAME=UV FIELDINFO=nm -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  0       0" -
M  V30 FIELDDATA="340"
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
  REQUIRE(mol2);

  auto mol3 = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.750000 -1.299038 0.000000 0
M  V30 2 C 0.000000 0.000000 0.000000 0
M  V30 3 C 1.500000 0.000000 0.000000 0
M  V30 4 C 2.250000 1.299038 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 3) FIELDNAME=Stereo -
M  V30 FIELDDATA="Cis"
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
  REQUIRE(mol3);

  std::vector<std::string> legends = {
      "datafield label bad placement1", "datafield label bad placement2",
      "datafield label bad placement3"};  //  std::vector<std::string> legends =
                                          //  {"datafield label bad
                                          //  placement2"};
  {
    MolDraw2DSVG drawer(300, 250);
    drawer.drawMolecule(*mol1, legends[0]);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub4519_1.svg");
    outs << text;
    outs.flush();
    check_file_hash("testGithub4519_1.svg");
  }
  {
    MolDraw2DSVG drawer(300, 250);
    drawer.drawMolecule(*mol2, legends[1]);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub4519_2.svg");
    outs << text;
    outs.flush();
    check_file_hash("testGithub4519_2.svg");
  }
  {
    MolDraw2DSVG drawer(300, 250);
    drawer.drawMolecule(*mol3, legends[2]);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub4519_3.svg");
    outs << text;
    outs.flush();
    check_file_hash("testGithub4519_3.svg");
  }

  {
    std::vector<ROMol *> mols;
    mols.push_back(mol1.get());
    mols.push_back(mol2.get());
    mols.push_back(mol3.get());
    MolDraw2DSVG drawer(900, 250, 300, 250);
    drawer.drawMolecules(mols, &legends);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testGithub4519_4.svg");
    outs << text;
    outs.flush();
    outs.close();
    check_file_hash("testGithub4519_4.svg");
  }
}

TEST_CASE("changing baseFontSize") {
  RDDepict::preferCoordGen = false;
  auto mol1 =
      "CC(C)C[C@H](NC(=O)[C@H](CCCCN)NC(=O)[C@H](CS)NC(=O)CNC(=O)[C@H](C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CC(C)C)NC(=O)CNC(=O)[C@H](C)NC(=O)[C@H](CS)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@@H](NC(=O)[C@H](CS)NC(=O)CNC(=O)[C@H](C)NC(=O)[C@H](CCCCN)NC(=O)CNC(=O)[C@H](C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)N)[C@@H](C)O)C(=O)O"_smiles;
  REQUIRE(mol1);
  MolDraw2DUtils::prepareMolForDrawing(*mol1);
  auto mol2 = "C[C@H](N)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](C)C(=O)NCC(=O)O"_smiles;
  REQUIRE(mol2);
  MolDraw2DUtils::prepareMolForDrawing(*mol2);
  SECTION("basics-large") {
    MolDraw2DSVG drawer(350, 300, -1, -1, 1);
    drawer.drawMolecule(*mol1);
    drawer.finishDrawing();
    CHECK(drawer.fontSize() == Approx(6.0).margin(0.1));
    auto text = drawer.getDrawingText();
    std::ofstream outs("testBaseFontSize.1a.svg");
    outs << text;
    outs.flush();
    check_file_hash("testBaseFontSize.1a.svg");
  }
  SECTION("increase size - large") {
    // here we change the base font size, but it doesn't matter since the
    // structure is big enough we end up stuck with the minimum font size.
    MolDraw2DSVG drawer(350, 300, -1, -1, 1);
    drawer.drawOptions().baseFontSize = 0.9;
    drawer.drawMolecule(*mol1);
    drawer.finishDrawing();
    CHECK(drawer.fontSize() == Approx(5.5).margin(.1));
    auto text = drawer.getDrawingText();
    std::ofstream outs("testBaseFontSize.1b.svg");
    outs << text;
    outs.flush();
    check_file_hash("testBaseFontSize.1b.svg");
  }
  SECTION("basics-small") {
    MolDraw2DSVG drawer(350, 300, -1, -1, 1);
    drawer.drawMolecule(*mol2);
    drawer.finishDrawing();
    CHECK(drawer.fontSize() == Approx(14.0).margin(0.1));
    auto text = drawer.getDrawingText();
    std::ofstream outs("testBaseFontSize.2a.svg");
    outs << text;
    outs.flush();
    check_file_hash("testBaseFontSize.2a.svg");
  }
  SECTION("increase size - smaller") {
    MolDraw2DSVG drawer(350, 300, -1, -1, 1);
    drawer.drawOptions().baseFontSize = 0.9;
    drawer.drawMolecule(*mol2);
    drawer.finishDrawing();
    CHECK(drawer.fontSize() == Approx(20.4).margin(0.1));
    auto text = drawer.getDrawingText();
    std::ofstream outs("testBaseFontSize.2b.svg");
    outs << text;
    outs.flush();
    check_file_hash("testBaseFontSize.2b.svg");
  }
}

TEST_CASE("flexicanvas: set canvas size automatically") {
  // note that these examples use Freetype if it's available.
  auto mol1 = "CCN(CC)CCn1nc2c3ccccc3sc3c(CNS(C)(=O)=O)ccc1c32"_smiles;
  REQUIRE(mol1);
  MolDraw2DUtils::prepareMolForDrawing(*mol1);

  auto mol2 = R"CTAB(
  Mrv2108 11192104292D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -5.2 -1.4 0 0
M  V30 2 O -5.2 -2.8 0 0
M  V30 3 C -3.7 -1.4 0 0
M  V30 4 C -3.7 -2.8 0 0 CFG=1
M  V30 5 N -2.5994 -3.9839 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 2 4
M  V30 4 1 3 4
M  V30 5 1 4 5 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
  REQUIRE(mol2);
  MolDraw2DUtils::prepareMolForDrawing(*mol2);
  SECTION("fixed canvas") {
    MolDraw2DSVG drawer(308, 223, -1, -1);
    drawer.drawMolecule(*mol1);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testFlexiCanvas.1a.svg");
    outs << text;
    outs.flush();
    check_file_hash("testFlexiCanvas.1a.svg");
  }
  SECTION("flexicanvas1") {
    MolDraw2DSVG drawer(-1, -1, -1, -1);
    drawer.drawMolecule(*mol1);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testFlexiCanvas.1b.svg");
    outs << text;
    outs.flush();
    check_file_hash("testFlexiCanvas.1b.svg");
  }
  SECTION("flexicanvas1") {
    MolDraw2DSVG drawer(-1, -1, -1, -1);
    drawer.drawOptions().scalingFactor = 30;
    drawer.drawOptions().baseFontSize = 0.6;
    drawer.drawMolecule(*mol1);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testFlexiCanvas.1c.svg");
    outs << text;
    outs.flush();
    check_file_hash("testFlexiCanvas.1c.svg");
  }
  SECTION("flexicanvas1") {
    MolDraw2DSVG drawer(-1, -1, -1, -1);
    drawer.drawOptions().scalingFactor = 30;
    drawer.drawOptions().fixedFontSize = 32;
    drawer.drawMolecule(*mol1);
    drawer.finishDrawing();
    CHECK(drawer.fontSize() == Approx(32).margin(0.1));
    auto text = drawer.getDrawingText();
    std::ofstream outs("testFlexiCanvas.1d.svg");
    outs << text;
    outs.flush();
    check_file_hash("testFlexiCanvas.1d.svg");
  }
  SECTION("square") {
    MolDraw2DSVG drawer(-1, -1, -1, -1);
    drawer.drawOptions().baseFontSize = 0.8;
    drawer.drawMolecule(*mol2);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testFlexiCanvas.2.svg");
    outs << text;
    outs.flush();
    check_file_hash("testFlexiCanvas.2.svg");
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  SECTION("square PNG no freetype") {
    MolDraw2DCairo drawer(-1, -1, -1, -1, true);
    drawer.drawOptions().baseFontSize = 0.8;
    drawer.drawMolecule(*mol2);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testFlexiCanvas.2a.png");
    outs << text;
    outs.flush();
    check_file_hash("testFlexiCanvas.2a.png");
  }
  SECTION("square PNG with freetype") {
    MolDraw2DCairo drawer(-1, -1, -1, -1, false);
    drawer.drawOptions().baseFontSize = 0.8;
    drawer.drawMolecule(*mol2);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testFlexiCanvas.2b.png");
    outs << text;
    outs.flush();
    check_file_hash("testFlexiCanvas.2b.png");
  }
#endif
  // semiflexicanvas - with freetype
  SECTION("semiflexicanvas1") {
    MolDraw2DSVG drawer(308, -1, -1, -1, false);
    drawer.drawOptions().scalingFactor = 30;
    drawer.drawOptions().baseFontSize = 0.6;
    drawer.drawMolecule(*mol1);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testSemiFlexiCanvas.1a.svg");
    outs << text;
    outs.flush();
    check_file_hash("testSemiFlexiCanvas.1a.svg");
  }
  SECTION("semiflexicanvas2") {
    MolDraw2DSVG drawer(-1, 223, -1, -1, false);
    drawer.drawOptions().scalingFactor = 30;
    drawer.drawOptions().baseFontSize = 0.6;
    drawer.drawMolecule(*mol1);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testSemiFlexiCanvas.1b.svg");
    outs << text;
    outs.flush();
    check_file_hash("testSemiFlexiCanvas.1b.svg");
  }
  SECTION("semiflexicanvas3") {
    auto mol3 = "ON"_smiles;
    REQUIRE(mol3);
    MolDraw2DSVG drawer(-1, 150, -1, -1, false);
    drawer.drawOptions().scalingFactor = 30;
    drawer.drawOptions().baseFontSize = 0.6;
    drawer.drawMolecule(*mol3);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testSemiFlexiCanvas.1c.svg");
    outs << text;
    outs.flush();
    check_file_hash("testSemiFlexiCanvas.1c.svg");
  }
  SECTION("reaction") {
    std::unique_ptr<ChemicalReaction> rxn(RxnSmartsToChemicalReaction(
        "[N:1]-[C:2]-[C:3](=[O:4])-[O:5].[N:6]-[C:7]-[C:8](=[O:9])-[O:10]>>[N:"
        "1]1-[C:2]-[C:3](=[O:4])-[N:6]-[C:7]-[C:8]-1=[O:9].[O:5]=[O:10]"));
    MolDraw2DSVG drawer(-1, -1, -1, -1, true);
    drawer.drawReaction(*rxn);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testFlexiCanvas.3.svg");
    outs << text;
    outs.flush();
    check_file_hash("testFlexiCanvas.3.svg");
  }
  SECTION("data labels") {
    auto mol1 = R"CTAB(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0
M  V30 2 C 1.299038 0.750000 0.000000 0
M  V30 3 C 2.598076 -0.000000 0.000000 0
M  V30 4 C 1.299038 2.250000 0.000000 0
M  V30 5 C 2.598076 3.000000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 2 4
M  V30 4 2 4 5
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(5 2 4 5 3 1) FIELDNAME="Lambda Max" FIELDINFO=nm -
M  V30 FIELDDATA="2222"
M  V30 END SGROUP
M  V30 END CTAB
M  END)CTAB"_ctab;
    REQUIRE(mol1);
    {
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawMolecule(*mol1);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.4a.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.4a.svg");
    }
    {
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawMolecule(*mol1, "legendary");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.4b.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.4b.svg");
    }
    {
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawMolecule(*mol1, "doubly\nlegendary");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.4c.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.4c.svg");
    }
    {
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawOptions().legendFraction = 0.25;
      drawer.drawOptions().legendFontSize = 32;
      drawer.drawMolecule(*mol1, "Hugely\nLegendary");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.4d.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.4d.svg");
    }
  }
  SECTION("including legends") {
    // add an atomNote so that we can compare font sizes
    mol1->getAtomWithIdx(0)->setProp(common_properties::atomNote, "n1");
    {
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawMolecule(*mol1, "legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.5a.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.5a.svg");
    }
    {
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawMolecule(*mol1, "legend\nwith two lines");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.5b.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.5b.svg");
    }
    {
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawOptions().scalingFactor = 45;
      drawer.drawMolecule(*mol1, "legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.5c.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.5c.svg");
    }
    {
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawOptions().scalingFactor = 10;
      drawer.drawMolecule(*mol1, "legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.5d.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.5d.svg");
    }
  }

  SECTION("partially flexicanvas (height) + legends") {
    // add an atomNote so that we can compare font sizes
    mol1->getAtomWithIdx(0)->setProp(common_properties::atomNote, "n1");
    {
      MolDraw2DSVG drawer(-1, 200);
      drawer.drawMolecule(*mol1, "legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.6a.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.6a.svg");
    }
    {
      MolDraw2DSVG drawer(-1, 200);
      drawer.drawMolecule(*mol1, "legend\nwith two lines");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.6b.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.6b.svg");
    }
    {
      MolDraw2DSVG drawer(-1, 200);
      drawer.drawOptions().scalingFactor = 45;
      drawer.drawMolecule(*mol1, "legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.6c.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.6c.svg");
    }
    {
      MolDraw2DSVG drawer(-1, 200);
      drawer.drawOptions().scalingFactor = 10;
      drawer.drawMolecule(*mol1, "legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.6d.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.6d.svg");
    }
  }

  SECTION("partially flexicanvas (width) + legends") {
    // add an atomNote so that we can compare font sizes
    mol1->getAtomWithIdx(0)->setProp(common_properties::atomNote, "n1");
    {
      MolDraw2DSVG drawer(300, -1);
      drawer.drawMolecule(*mol1, "legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.7a.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.7a.svg");
    }
    {
      MolDraw2DSVG drawer(300, -1);
      drawer.drawMolecule(*mol1, "legend\nwith two lines");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.7b.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.7b.svg");
    }
    {
      MolDraw2DSVG drawer(300, -1);
      drawer.drawOptions().scalingFactor = 45;
      drawer.drawMolecule(*mol1, "legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.7c.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.7c.svg");
    }
    {
      MolDraw2DSVG drawer(300, -1);
      drawer.drawOptions().scalingFactor = 10;
      drawer.drawMolecule(*mol1, "legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testFlexiCanvas.7d.svg");
      outs << text;
      outs.flush();
      check_file_hash("testFlexiCanvas.7d.svg");
    }
  }
}

TEST_CASE("Github #4764") {
  SECTION("basics") {
    auto mol = "c1ccccc1-C1CCCCC1"_smiles;
    REQUIRE(mol);
    std::vector<int> highlights{6, 7, 8, 9, 10, 11};
    {
      MolDraw2DSVG drawer(200, 150);
      drawer.drawMolecule(*mol, "highlight", &highlights);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testGithub4764.sz1.svg");
      outs << text;
      outs.flush();
      check_file_hash("testGithub4764.sz1.svg");
    }
    {
      MolDraw2DSVG drawer(400, 350);
      drawer.drawMolecule(*mol, "highlight", &highlights);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testGithub4764.sz2.svg");
      outs << text;
      outs.flush();
      check_file_hash("testGithub4764.sz2.svg");
    }
    {
      MolDraw2DSVG drawer(800, 700);
      drawer.drawMolecule(*mol, "highlight", &highlights);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testGithub4764.sz3.svg");
      outs << text;
      outs.flush();
      check_file_hash("testGithub4764.sz3.svg");
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(200, 150);
      drawer.drawMolecule(*mol, "highlight", &highlights);
      drawer.finishDrawing();
      drawer.writeDrawingText("testGithub4764.sz1.png");
      check_file_hash("testGithub4764.sz1.png");
    }
    {
      MolDraw2DCairo drawer(400, 350);
      drawer.drawMolecule(*mol, "highlight", &highlights);
      drawer.finishDrawing();
      drawer.writeDrawingText("testGithub4764.sz2.png");
      check_file_hash("testGithub4764.sz2.png");
    }
    {
      MolDraw2DCairo drawer(800, 700);
      drawer.drawMolecule(*mol, "highlight", &highlights);
      drawer.finishDrawing();
      drawer.writeDrawingText("testGithub4764.sz3.png");
      check_file_hash("testGithub4764.sz3.png");
    }
#endif
    // check_file_hash("testGithub4538.svg");
  }
}

TEST_CASE("drawArc starting from wrong angle") {
  SECTION("basics") {
    auto mol = R"CTAB(
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -1.2135   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5844    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.2135   -0.7027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    0.7238    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    0.7238    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6317    1.9374    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6401   -1.1663    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    2.6401   -1.1663    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.6317    1.9374    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  5  1  2  0
  1  7  1  0
  3  8  1  0
  4  9  1  0
M  END)CTAB"_ctab;
    REQUIRE(mol);
    {
      MolDraw2DSVG drawer(400, 350);
      drawer.drawOptions().noAtomLabels = true;
      drawer.drawMolecule(*mol, "drawArc");
      drawer.setFillPolys(false);
      drawer.setColour({1, 0, 0});
      drawer.drawArc(mol->getConformer().getAtomPos(3), 0.3, -72, 54);
      drawer.drawArc(mol->getConformer().getAtomPos(0), 0.3, -162, -36);
      drawer.drawArc(mol->getConformer().getAtomPos(4), 0.3, 126, 252);
      drawer.drawArc(mol->getConformer().getAtomPos(2), 0.3, -18, 108);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testDrawArc1.svg");
      outs << text;
      outs.flush();
      check_file_hash("testDrawArc1.svg");
    }
  }
}

TEST_CASE("wedged bonds to metals drawn in the wrong direction") {
  SECTION("basics") {
    auto m = R"CTAB(
  Mrv2108 01092205442D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 F 10.6667 -0.75 0 0
M  V30 2 Pt 10.6667 -2.29 0 0 CFG=1
M  V30 3 Cl 12.2067 -2.29 0 0
M  V30 4 C 10.6667 -3.83 0 0
M  V30 5 O 9.1267 -2.29 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 2 4 CFG=1
M  V30 4 1 2 5 CFG=3
M  V30 END BOND
M  V30 END CTAB
M  END)CTAB"_ctab;
    m->getBondWithIdx(2)->setBondDir(Bond::BondDir::BEGINWEDGE);
    m->getBondWithIdx(3)->setBondDir(Bond::BondDir::BEGINDASH);
    MolDraw2DSVG drawer(250, 200);
    assignBWPalette(drawer.drawOptions().atomColourPalette);
    drawer.drawMolecule(*m, "check wedges");
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testMetalWedges.svg");
    outs << text;
    outs.flush();
    check_file_hash("testMetalWedges.svg");
  }
}

TEST_CASE("vary proportion of panel for legend", "[drawing]") {
  SECTION("basics") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);
    // These look a bit pants with NO_FREETYPE=true, but much better with
    // Freetype.
    {
      // default legend
      MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1, "default legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testVariableLegend_1.svg");
      outs << text;
      outs.flush();
      CHECK(text.find("<text x='34.5' y='195.0' class='legend' "
                      "style='font-size:16px;") != std::string::npos);
      check_file_hash("testVariableLegend_1.svg");
    }
    {
      // 1/4 of panel
      MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
      drawer.drawOptions().legendFraction = 0.25;
      drawer.drawOptions().legendFontSize = 32;
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1, "massive legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testVariableLegend_2.svg");
      outs << text;
      outs.flush();
      CHECK(text.find("<text x='1.6' y='195.0' class='legend' "
                      "style='font-size:31px;") != std::string::npos);
      check_file_hash("testVariableLegend_2.svg");
    }
    {
      // tiny
      MolDraw2DSVG drawer(200, 200, -1, -1, NO_FREETYPE);
      drawer.drawOptions().legendFraction = 0.05;
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1, "small legend");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testVariableLegend_3.svg");
      outs << text;
      outs.flush();
      CHECK(text.find("<text x='84.7' y='195.0' class='legend' "
                      "style='font-size:6px;") != std::string::npos);
      check_file_hash("testVariableLegend_3.svg");
    }
  }
}

TEST_CASE(
    "Github 5061 - draw reaction with no reagents and scaleBondWidth true") {
  SECTION("basics") {
    std::string data = R"RXN($RXN

  Mrv16425    091201171606

  0  1
$MOL

  Mrv1642509121716062D

  2  1  0  0  0  0            999 V2000
    3.5357    0.0000    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    2.7107    0.0000    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  RGP  2   1   1   2   2


M  END)RXN";
    {
      std::unique_ptr<ChemicalReaction> rxn{RxnBlockToChemicalReaction(data)};
      MolDraw2DSVG drawer(450, 200);
      drawer.drawOptions().scaleBondWidth = true;
      drawer.drawReaction(*rxn);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testGithub_5061.svg");
      outs << text;
      outs.flush();
      check_file_hash("testGithub_5061.svg");
    }
  }
}

TEST_CASE("Github 5185 - don't draw atom indices between double bond") {
  SECTION("basics") {
    auto m1 = "OC(=O)CCCC(=O)O"_smiles;
    REQUIRE(m1);
    {
      // default legend
      MolDraw2DSVG drawer(400, 200, -1, -1);
      drawer.drawOptions().addAtomIndices = true;
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testGithub_5185.svg");
      outs << text;
      outs.flush();
#ifdef RDK_BUILD_FREETYPE_SUPPORT
      // the 2nd note
      CHECK(text.find("<path class='note' d='M 93.4 129.9") !=
            std::string::npos);
      check_file_hash("testGithub_5185.svg");
#else
      CHECK(text.find("<text x='90.4' y='130.3' class='note' ") !=
            std::string::npos);
#endif
    }
  }
}

TEST_CASE(
    "Github 5259 - drawReaction should not fail when prepareMolsBeforeDrawing "
    "is false") {
  SECTION("basics") {
    auto rxn = "[CH3:1][OH:2]>>[CH2:1]=[OH0:2]"_rxnsmarts;
    REQUIRE(rxn);
    MolDraw2DSVG drawer(400, 200, -1, -1);
    drawer.drawOptions().prepareMolsBeforeDrawing = false;
    REQUIRE_NOTHROW(drawer.drawReaction(*rxn));
  }
}

TEST_CASE("Github 5269 - bad index positions with highlights") {
  SECTION("basics") {
    auto m1 = "CC(=O)Oc1c(C(=O)O)cccc1"_smiles;
    auto q1 = "CC(=O)Oc1c(C(=O)O)cccc1"_smarts;
    REQUIRE(m1);
    REQUIRE(q1);
    {
      std::vector<int> hit_atoms;
      std::vector<MatchVectType> hits_vect;
      SubstructMatch(*m1, *q1, hits_vect);
      for (size_t i = 0; i < hits_vect.size(); ++i) {
        for (size_t j = 0; j < hits_vect[i].size(); ++j) {
          hit_atoms.push_back(hits_vect[i][j].second);
        }
      }
      std::vector<int> hit_bonds;
      for (int i : hit_atoms) {
        for (int j : hit_atoms) {
          if (i > j) {
            Bond *bnd = m1->getBondBetweenAtoms(i, j);
            if (bnd) {
              hit_bonds.push_back(bnd->getIdx());
            }
          }
        }
      }
      {
        MolDraw2DSVG drawer(400, 400, -1, -1);
        drawer.drawOptions().addAtomIndices = true;
        drawer.drawMolecule(*m1, &hit_atoms, &hit_bonds);
        drawer.finishDrawing();
        auto text = drawer.getDrawingText();
        std::ofstream outs("testGithub_5269_1.svg");
        outs << text;
        outs.flush();
#ifdef RDK_BUILD_FREETYPE_SUPPORT
        check_file_hash("testGithub_5269_1.svg");
#endif
      }
    }
  }
  {
    auto m2 = "CN(C)C(C)C=O"_smiles;
    REQUIRE(m2);
    std::vector<int> hit_atoms{0, 1, 2};
    auto atom = m2->getAtomWithIdx(0);
    atom->setProp(common_properties::atomNote, "0.91");
    atom = m2->getAtomWithIdx(1);
    atom->setProp(common_properties::atomNote, "1.03");
    atom = m2->getAtomWithIdx(2);
    atom->setProp(common_properties::atomNote, "0.74");
    MolDraw2DSVG drawer(400, 400, -1, -1);
    drawer.drawMolecule(*m2, &hit_atoms);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testGithub_5269_2.svg");
    outs << text;
    outs.flush();
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    check_file_hash("testGithub_5269_2.svg");
#endif
  }
}

#ifdef RDK_BUILD_CAIRO_SUPPORT
TEST_CASE("drawing doesn't destroy reaction properties", "[drawing]") {
  auto rxn = "[CH3:1][OH:2]>>[CH2:1]=[OH0:2]"_rxnsmarts;
  REQUIRE(rxn);
  MolDraw2DCairo drawer(400, 200);
  bool highlightByReactant = true;
  drawer.drawReaction(*rxn, highlightByReactant);
  drawer.finishDrawing();
  auto png = drawer.getDrawingText();
  std::unique_ptr<ChemicalReaction> rxn2{PNGStringToChemicalReaction(png)};
  REQUIRE(rxn2);
  CHECK(rxn->getReactants()[0]->getAtomWithIdx(0)->getAtomMapNum() == 1);
  CHECK(rxn->getReactants()[0]->getAtomWithIdx(1)->getAtomMapNum() == 2);
  CHECK(rxn2->getReactants()[0]->getAtomWithIdx(0)->getAtomMapNum() == 1);
  CHECK(rxn2->getReactants()[0]->getAtomWithIdx(1)->getAtomMapNum() == 2);
}
#endif

TEST_CASE("Class values in SVG for wavy bonds.") {
  SECTION("basics") {
    auto m1 = R"CTAB(mol1
  ChemDraw05162216032D

 11 11  0  0  0  0  0  0  0  0999 V2000
    1.1514    0.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1514    0.0788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9360   -0.1762    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.4209    0.4913    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9360    1.1587    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4369   -0.3337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2775    0.0788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9920   -0.3337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7065    0.0788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4209   -0.3337    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4369   -1.1587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  1  1  0
  2  6  1  0
  6  7  1  0
  7  8  2  3
  8  9  1  0
  9 10  3  0
  6 11  1  4
M  END)CTAB"_ctab;
    REQUIRE(m1);
    auto b10 = m1->getBondWithIdx(10);
    b10->setBondDir(Bond::UNKNOWN);
    MolDraw2DSVG drawer(400, 400, -1, -1);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    CHECK(text.find("<path class='bond-10 atom-5 atom-10'") !=
          std::string::npos);
    std::ofstream outs("test_classes_wavy_bonds.svg");
    outs << text;
    outs.flush();
#ifdef RDK_BUILD_FREETYPE_SUPPORT
    check_file_hash("test_classes_wavy_bonds.svg");
#endif
  }
}

TEST_CASE("GitHub #5383: cairo error when using similarity maps", "") {
  auto m1 = "C1N[C@@H]2OCC12"_smiles;
  REQUIRE(m1);
  MolDraw2DUtils::prepareMolForDrawing(*m1);
  const auto conf = m1->getConformer();
  std::vector<Point2D> cents(conf.getNumAtoms());
  std::vector<double> weights(conf.getNumAtoms());
  std::vector<double> widths(conf.getNumAtoms());
  for (size_t i = 0; i < conf.getNumAtoms(); ++i) {
    cents[i] = Point2D(conf.getAtomPos(i).x, conf.getAtomPos(i).y);
    weights[i] = 1;
    widths[i] = 0.4 * PeriodicTable::getTable()->getRcovalent(
                          m1->getAtomWithIdx(i)->getAtomicNum());
  }

  SECTION("svg basics") {
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    drawer.drawOptions().padding = 0.1;

    drawer.clearDrawing();
    std::vector<double> levels;
    MolDraw2DUtils::contourAndDrawGaussians(
        drawer, cents, weights, widths, 10, levels,
        MolDraw2DUtils::ContourParams(), m1.get());

    drawer.drawOptions().clearBackground = false;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    CHECK(text.find("width='250px' height='250px' viewBox='0 0 250 250'>") !=
          std::string::npos);
    std::ofstream outs("github5383_1.svg");
    outs << text;
    outs.flush();
    check_file_hash("github5383_1.svg");
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  SECTION("cairo basics") {
    MolDraw2DCairo drawer(250, 250, -1, -1, NO_FREETYPE);
    drawer.drawOptions().padding = 0.1;

    drawer.clearDrawing();
    std::vector<double> levels;
    MolDraw2DUtils::contourAndDrawGaussians(
        drawer, cents, weights, widths, 10, levels,
        MolDraw2DUtils::ContourParams(), m1.get());

    drawer.drawOptions().clearBackground = false;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    drawer.writeDrawingText("github5383_1.png");
    check_file_hash("github5383_1.png");
  }
#endif
}

TEST_CASE("github #5156") {
  SECTION("basics") {
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> m{SmilesToMol("c1ccnc1", ps)};
    REQUIRE(m);
    unsigned int failed;
    MolOps::sanitizeMol(*m, failed,
                        MolOps::SANITIZE_ALL ^ MolOps::SANITIZE_KEKULIZE);
    MolDraw2DSVG d2d(200, 200);
    d2d.drawOptions().prepareMolsBeforeDrawing = false;
    d2d.drawMolecule(*m);
    d2d.finishDrawing();
    auto text = d2d.getDrawingText();
    // CHECK(text.find("width='250px' height='250px' viewBox='0 0 250 250'>") !=
    //       std::string::npos);
    std::ofstream outs("github5156_1.svg");
    outs << text;
    outs.flush();
    check_file_hash("github5156_1.svg");
  }
  SECTION("as reported") {
    auto m =
        "[#6](:,-[#6]-,:[#7]-,:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1):,-[#6]:,-[#7]:,-[#6]"_smarts;
    REQUIRE(m);
    MolDraw2DSVG d2d(200, 200);
    d2d.drawOptions().prepareMolsBeforeDrawing = false;
    d2d.drawMolecule(*m);
    d2d.finishDrawing();
    auto text = d2d.getDrawingText();
    // CHECK(text.find("width='250px' height='250px' viewBox='0 0 250 250'>") !=
    //       std::string::npos);
    std::ofstream outs("github5156_2.svg");
    outs << text;
    outs.flush();
    check_file_hash("github5156_2.svg");
  }
  SECTION("check no wedging") {
    // if we aren't preparing molecules, we won't end up with wedging in this
    // case
    auto m = "C[C@H](F)Cl"_smiles;
    REQUIRE(m);
    MolDraw2DSVG d2d(200, 200);
    d2d.drawOptions().prepareMolsBeforeDrawing = false;
    d2d.drawMolecule(*m);
    d2d.finishDrawing();
    auto text = d2d.getDrawingText();
    CHECK(text.find(" Z' style='fill=#000000") == std::string::npos);
    std::ofstream outs("github5156_3.svg");
    outs << text;
    outs.flush();
    check_file_hash("github5156_3.svg");
  }
}

TEST_CASE("ACS 1996 mode") {
  SECTION("basics") {
    std::string nameBase = "acs1996_";
#if 1

    {
      auto m = R"CTAB(mol1
  ChemDraw05162216032D

 11 11  0  0  0  0  0  0  0  0999 V2000
    1.1514    0.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1514    0.0788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9360   -0.1762    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.4209    0.4913    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9360    1.1587    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4369   -0.3337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2775    0.0788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9920   -0.3337    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7065    0.0788    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4209   -0.3337    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.4369   -1.1587    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  1  1  0
  2  6  1  0
  6  7  1  0
  7  8  2  3
  8  9  1  0
  9 10  3  0
  6 11  1  4
M  END)CTAB"_ctab;
      REQUIRE(m);
      {
        MolDraw2DSVG drawer(-1, -1);
        MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 1", nullptr, nullptr);
        drawer.finishDrawing();
        std::string text = drawer.getDrawingText();
        std::ofstream outs(nameBase + "1.svg");
        outs << text;
        outs.flush();
        outs.close();
        check_file_hash(nameBase + "1.svg");
      }
#ifdef RDK_BUILD_CAIRO_SUPPORT
      {
        MolDraw2DCairo drawer(-1, -1);
        MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 1", nullptr, nullptr);
        drawer.finishDrawing();
        drawer.writeDrawingText(nameBase + "1.png");
        check_file_hash(nameBase + "1.png");
      }
#endif
    }
#endif
#if 1
    {
      auto m = R"CTAB(mol2
  ChemDraw06062216302D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 11 11 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 1.151400 0.903800 0.000000 0
M  V30 2 C 1.151400 0.078800 0.000000 0
M  V30 3 N 1.935999 -0.176199 0.000000 0
M  V30 4 C 2.420899 0.491300 0.000000 0
M  V30 5 N 1.935999 1.158700 0.000000 0
M  V30 6 C 0.436900 -0.333700 0.000000 0
M  V30 7 C -0.277500 0.078800 0.000000 0
M  V30 8 C -0.992000 -0.333700 0.000000 0
M  V30 9 C -1.706500 0.078800 0.000000 0
M  V30 10 N -2.420899 -0.333700 0.000000 0
M  V30 11 C 0.436900 -1.158700 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 1 5 1
M  V30 6 1 2 6
M  V30 7 1 6 7
M  V30 8 2 7 8 CFG=2
M  V30 9 1 8 9
M  V30 10 3 9 10
M  V30 11 1 6 11 CFG=3
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 6)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      MolDraw2DSVG drawer(-1, -1);
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 2", nullptr, nullptr);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "2.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "2.svg");
#ifdef RDK_BUILD_CAIRO_SUPPORT
      {
        MolDraw2DCairo drawer(-1, -1);
        MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 2", nullptr, nullptr);
        drawer.finishDrawing();
        drawer.writeDrawingText(nameBase + "2.png");
        check_file_hash(nameBase + "2.png");
      }
#endif
    }
#endif
#if 1
    {
      auto m = "C[C@H](I)CC(Cl)C[C@@H](F)C"_smiles;
      m->setProp<std::string>("_Name", "mol3");
      REQUIRE(m);
      MolDraw2DUtils::prepareMolForDrawing(*m);
      MolDraw2DSVG drawer(-1, -1);
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 3", nullptr, nullptr);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "3.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "3.svg");
    }
#endif
#if 1
    {
      auto m = "CC(I)CC(Cl)CC(F)C"_smiles;
      m->setProp<std::string>("_Name", "mol4");
      REQUIRE(m);
      MolDraw2DUtils::prepareMolForDrawing(*m);
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawOptions().unspecifiedStereoIsUnknown = true;
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 4", nullptr, nullptr);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "4.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "4.svg");
    }
#endif
#if 1
    {
      auto m = R"CTAB(mol5
  ChemDraw06112209342D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 30 33 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -2.240810 -1.031250 0.000000 0
M  V30 2 C -2.240810 -0.206250 0.000000 0
M  V30 3 C -2.955281 0.206250 0.000000 0
M  V30 4 C -2.955281 1.031250 0.000000 0
M  V30 5 C -3.669752 1.443750 0.000000 0
M  V30 6 C -4.384224 1.031250 0.000000 0
M  V30 7 C -4.384224 0.206250 0.000000 0
M  V30 8 C -3.669752 -0.206250 0.000000 0
M  V30 9 Cl -3.669752 -1.031250 0.000000 0
M  V30 10 F -5.098694 -0.206250 0.000000 0
M  V30 11 Cl -2.240810 1.443750 0.000000 0
M  V30 12 O -1.526340 0.206250 0.000000 0
M  V30 13 C -0.811869 -0.206250 0.000000 0
M  V30 14 C -0.811869 -1.031250 0.000000 0
M  V30 15 N -0.097397 -1.443750 0.000000 0
M  V30 16 C 0.617074 -1.031250 0.000000 0
M  V30 17 C 0.617074 -0.206250 0.000000 0
M  V30 18 C -0.097397 0.206250 0.000000 0
M  V30 19 C 1.331544 0.206250 0.000000 0
M  V30 20 C 2.085220 -0.129308 0.000000 0
M  V30 21 N 2.637252 0.483787 0.000000 0
M  V30 22 N 2.224752 1.198258 0.000000 0
M  V30 23 C 1.417781 1.026730 0.000000 0
M  V30 24 C 3.457733 0.397551 0.000000 0
M  V30 25 C 3.942655 1.064990 0.000000 0
M  V30 26 C 4.763136 0.978754 0.000000 0
M  V30 27 N 5.098694 0.225079 0.000000 0
M  V30 28 C 4.613771 -0.442361 0.000000 0
M  V30 29 C 3.793290 -0.356124 0.000000 0
M  V30 30 N -1.526340 -1.443750 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 8
M  V30 8 1 3 8
M  V30 9 1 8 9
M  V30 10 1 7 10
M  V30 11 1 4 11
M  V30 12 1 2 12
M  V30 13 1 12 13
M  V30 14 2 13 14
M  V30 15 1 14 15
M  V30 16 2 15 16
M  V30 17 1 16 17
M  V30 18 2 17 18
M  V30 19 1 13 18
M  V30 20 1 17 19
M  V30 21 2 19 20
M  V30 22 1 20 21
M  V30 23 1 21 22
M  V30 24 2 22 23
M  V30 25 1 19 23
M  V30 26 1 21 24
M  V30 27 1 24 25
M  V30 28 1 25 26
M  V30 29 1 26 27
M  V30 30 1 27 28
M  V30 31 1 28 29
M  V30 32 1 24 29
M  V30 33 1 14 30
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 2)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)CTAB"_ctab;
      MolDraw2DSVG drawer(-1, -1);
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 5", nullptr, nullptr);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "5.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "5.svg");
    }
#endif
#if 1
    {
      auto m = R"CTAB(mol6
  ChemDraw06132212082D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 16 15 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.427702 0.413216 0.000000 0
M  V30 2 C -0.715712 0.824284 0.000000 0
M  V30 3 C -2.139692 0.824284 0.000000 0
M  V30 4 C -0.003721 0.413216 0.000000 0
M  V30 5 C 0.708270 0.824284 0.000000 0
M  V30 6 C 1.420260 0.413216 0.000000 0
M  V30 7 C 0.708270 1.646420 0.000000 0
M  V30 8 C -1.427702 -0.408920 0.000000 0
M  V30 9 C -0.715712 -0.819988 0.000000 0
M  V30 10 C -2.139692 -0.819988 0.000000 0
M  V30 11 C -0.003721 -0.408920 0.000000 0
M  V30 12 C 2.134731 0.825716 0.000000 0
M  V30 13 C 0.710751 -0.821420 0.000000 0
M  V30 14 C 1.425221 -0.408920 0.000000 0
M  V30 15 C 2.139692 -0.821420 0.000000 0
M  V30 16 C 2.139692 -1.646420 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 1 3
M  V30 3 1 2 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 1 5 7
M  V30 7 1 1 8
M  V30 8 2 8 9 CFG=2
M  V30 9 1 8 10
M  V30 10 1 9 11
M  V30 11 2 6 12
M  V30 12 2 11 13
M  V30 13 1 13 14
M  V30 14 2 14 15
M  V30 15 1 15 16
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      MolDraw2DSVG drawer(-1, -1);
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 6", nullptr, nullptr);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "6.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "6.svg");
    }
#endif
#if 1
    {
      auto m = R"CTAB(mol7
  ChemDraw06192209312D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 18 17 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.147146 0.827156 0.000000 0
M  V30 2 C -1.432676 1.239655 0.000000 0
M  V30 3 O -2.861616 1.239655 0.000000 0
M  V30 4 C -0.718205 0.827156 0.000000 0
M  V30 5 C -0.003735 1.239655 0.000000 0
M  V30 6 C 0.710736 0.827156 0.000000 0
M  V30 7 C -0.003735 2.064654 0.000000 0
M  V30 8 C -2.147146 0.002156 0.000000 0
M  V30 9 C -1.432676 -0.410344 0.000000 0
M  V30 10 C -2.861616 -0.410344 0.000000 0
M  V30 11 C -0.718205 0.002156 0.000000 0
M  V30 12 S 1.427695 1.241093 0.000000 0
M  V30 13 C -0.001244 -0.411781 0.000000 0
M  V30 14 C 0.715714 0.002156 0.000000 0
M  V30 15 C 1.432674 -0.411781 0.000000 0
M  V30 16 C 1.432674 -1.239654 0.000000 0
M  V30 17 C 2.147145 -1.652154 0.000000 0
M  V30 18 C 2.861616 -2.064654 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 1 3
M  V30 3 1 2 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 1 5 7
M  V30 7 1 1 8
M  V30 8 2 8 9 CFG=2
M  V30 9 1 8 10
M  V30 10 1 9 11
M  V30 11 2 6 12
M  V30 12 2 11 13
M  V30 13 1 13 14
M  V30 14 2 14 15
M  V30 15 1 15 16
M  V30 16 2 16 17
M  V30 17 2 17 18
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      MolDraw2DSVG drawer(-1, -1);
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 7", nullptr, nullptr);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "7.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "7.svg");
    }
#endif
#if 1
    {
      auto m = R"CTAB(mol8
  ChemDraw07042207302D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 18 17 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -2.500648 -0.206250 0.000000 0
M  V30 2 C -1.786177 0.206250 0.000000 0
M  V30 3 C -1.071707 -0.206250 0.000000 0
M  V30 4 C -0.357237 0.206250 0.000000 0
M  V30 5 C 0.357236 -0.206250 0.000000 0
M  V30 6 C 1.071707 0.206250 0.000000 0
M  V30 7 C 1.786177 -0.206250 0.000000 0
M  V30 8 C 2.500648 0.206250 0.000000 0
M  V30 9 C -1.786177 1.031251 0.000000 0
M  V30 10 C -2.500648 1.443750 0.000000 0
M  V30 11 C -0.357237 1.031251 0.000000 0
M  V30 12 C -1.071707 1.443750 0.000000 0
M  V30 13 C 0.357236 1.443750 0.000000 0
M  V30 14 C 1.786177 -1.031250 0.000000 0
M  V30 15 C 2.500648 -1.443750 0.000000 0
M  V30 16 Cl 0.357236 -1.031250 0.000000 0
M  V30 17 C -1.071707 -1.031250 0.000000 0
M  V30 18 C 1.071707 1.031250 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 2 9 CFG=3
M  V30 9 1 9 10
M  V30 10 1 4 11 CFG=3
M  V30 11 1 11 12
M  V30 12 1 11 13
M  V30 13 1 7 14 CFG=3
M  V30 14 2 14 15
M  V30 15 1 5 16 CFG=3
M  V30 16 1 3 17
M  V30 17 1 6 18 CFG=3
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(5 2 4 5 6 7)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      MolDraw2DSVG drawer(-1, -1);
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 8", nullptr, nullptr);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "8.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "8.svg");
    }
#endif
#if 1
    {
      auto m = R"CTAB(mol9
  ChemDraw06302215142D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 22 21 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C -2.857884 -0.412500 0.000000 0
M  V30 2 C -2.143413 0.000000 0.000000 0
M  V30 3 C -1.428942 -0.412500 0.000000 0
M  V30 4 C -0.714471 0.000000 0.000000 0
M  V30 5 C 0.000000 -0.412500 0.000000 0
M  V30 6 C 0.714471 0.000000 0.000000 0
M  V30 7 C 1.428941 -0.412500 0.000000 0
M  V30 8 C 2.143413 0.000000 0.000000 0
M  V30 9 C -2.143413 0.825000 0.000000 0
M  V30 10 C -2.857884 1.237500 0.000000 0
M  V30 11 C -0.714471 0.825000 0.000000 0
M  V30 12 C -1.428942 1.237500 0.000000 0
M  V30 13 C 0.000000 1.237500 0.000000 0
M  V30 14 C 1.428941 -1.237500 0.000000 0
M  V30 15 C 2.143413 -1.650000 0.000000 0
M  V30 16 Cl 0.000000 -1.237500 0.000000 0
M  V30 17 C -1.428942 -1.237500 0.000000 0
M  V30 18 C 2.857884 -0.412500 0.000000 0
M  V30 19 C 2.143413 0.825000 0.000000 0
M  V30 20 C 1.428941 1.237500 0.000000 0
M  V30 21 C 2.857884 1.237500 0.000000 0
M  V30 22 C 2.143413 1.650000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 2 9 CFG=1
M  V30 9 1 9 10
M  V30 10 1 4 11 CFG=1
M  V30 11 1 11 12
M  V30 12 1 11 13
M  V30 13 1 7 14 CFG=1
M  V30 14 2 14 15
M  V30 15 1 5 16 CFG=1
M  V30 16 1 3 17
M  V30 17 1 8 18
M  V30 18 1 8 19 CFG=1
M  V30 19 1 19 20
M  V30 20 1 19 21
M  V30 21 1 19 22
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(5 2 4 5 7 8)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawOptions().useMolBlockWedging = true;
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 9", nullptr, nullptr);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "9.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "9.svg");
    }
#endif
#if 1
    {
      auto m = R"CTAB(mol10
  ChemDraw07062213362D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 8 7 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.357235 -1.031250 0.000000 0
M  V30 2 C -0.357235 -0.206250 0.000000 0
M  V30 3 N -1.071706 0.206250 0.000000 0
M  V30 4 O -1.786177 -0.206250 0.000000 0
M  V30 5 C 0.357236 0.206250 0.000000 0
M  V30 6 N 1.071706 -0.206250 0.000000 0
M  V30 7 O 1.786177 0.206250 0.000000 0
M  V30 8 C 0.357236 1.031250 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3 CFG=2
M  V30 3 1 3 4
M  V30 4 1 2 5
M  V30 5 2 5 6 CFG=2
M  V30 6 1 6 7
M  V30 7 1 5 8
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
      REQUIRE(m);
      MolDraw2DSVG drawer(-1, -1);
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 10", nullptr, nullptr);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "10.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "10.svg");
    }
#endif
#if 1
    {
      auto m = "CCOC(=O)Nc1ccc(SCC2COC(Cn3ccnc3)(c3ccc(Cl)cc3Cl)O2)cc1"_smiles;
      REQUIRE(m);
      MolDraw2DUtils::prepareMolForDrawing(*m);
      std::vector<int> highlight_atoms{17, 18, 19, 20, 21, 6, 7, 8, 9, 31, 32};
      std::map<int, DrawColour> atom_highlight_colors;
      atom_highlight_colors[8] = DrawColour(1.0, 1.0, 0.0);
      atom_highlight_colors[31] = DrawColour(0.0, 1.0, 1.0);
      std::vector<int> highlight_bonds{0, 1, 2, 11, 15, 19};
      std::map<int, DrawColour> bond_highlight_colors;
      bond_highlight_colors[0] = DrawColour(1.0, 1.0, 0.0);
      bond_highlight_colors[11] = DrawColour(0.0, 1.0, 1.0);
      MolDrawOptions options;
      options.circleAtoms = true;
      options.highlightColour = DrawColour(1, .5, .5);
      options.continuousHighlight = true;
      MolDraw2DSVG drawer(-1, -1);
      drawer.drawOptions() = options;
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 11", &highlight_atoms,
                                     &highlight_bonds, &atom_highlight_colors,
                                     &bond_highlight_colors);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "11.svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + "11.svg");
    }
#endif
#if 1
    auto drawnBondLength = [&](const std::string &r1,
                               const std::string &t) -> double {
      std::regex regex1(r1);
      auto match_begin = std::sregex_iterator(t.begin(), t.end(), regex1);
      auto match_end = std::sregex_iterator();
      std::vector<Point2D> ends;
      for (std::sregex_iterator i = match_begin; i != match_end; ++i) {
        std::smatch match = *i;
        ends.push_back(Point2D(std::stod(match[1]), std::stod(match[2])));
        ends.push_back(Point2D(std::stod(match[3]), std::stod(match[4])));
      }
      return (ends[0] - ends[1]).length();
    };

    {
      // make sure it also works with an arbitrarily sized drawer.
      auto m = "c1ccccc1"_smiles;
      REQUIRE(m);
      MolDraw2DUtils::prepareMolForDrawing(*m);
      MolDraw2DSVG drawer(500, 500);
      MolDraw2DUtils::drawMolACS1996(drawer, *m, "Mol 12", nullptr, nullptr);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + "12.svg");
      outs << text;
      outs.flush();
      outs.close();

      std::string regex =
          R"(class='bond-0 atom-0 atom-1' d='M ([\d.]*),([\d.]*) L ([\d.]*),([\d.]*)')";
      double dbl = drawnBondLength(regex, text);
      // the bonds should all be 14.4 long, but the SVG is only written
      // to 1 decimal place, so rounding errors are largish.
      REQUIRE(dbl == Approx(14.4253));
      regex =
          R"(class='bond-1 atom-1 atom-2' d='M ([\d.]*),([\d.]*) L ([\d.]*),([\d.]*)')";
      dbl = drawnBondLength(regex, text);
      REQUIRE(dbl == Approx(14.4));
      check_file_hash(nameBase + "12.svg");
    }
#endif
  }
}

TEST_CASE("Unspecified stereochemistry means unknown.", "") {
  auto m1 = "ClC(I)(F)C=CC"_smiles;
  REQUIRE(m1);
  MolDraw2DUtils::prepareMolForDrawing(*m1);
  MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
  drawer.drawOptions().unspecifiedStereoIsUnknown = true;
  drawer.drawMolecule(*m1);
  drawer.finishDrawing();
  auto text = drawer.getDrawingText();
  std::ofstream outs("test_unspec_stereo.svg");
  outs << text;
  outs.flush();

  std::regex regex1("class='bond-2 atom-1 atom-3' .*M.*C");
  std::smatch wavyMatch;
  REQUIRE(std::regex_search(text, wavyMatch, regex1));
  REQUIRE(wavyMatch.size() == 1);

  std::regex regex2(
      "class='bond-4 atom-4 atom-5' d='M 118\\..*,135\\..* L "
      "68\\..*,107\\..*'");
  std::smatch cross1Match;
  REQUIRE(std::regex_search(text, cross1Match, regex2));
  REQUIRE(cross1Match.size() == 1);

  std::regex regex3(
      "class='bond-4 atom-4 atom-5' d='M 117\\..*,145\\..* L "
      "74\\..*,100\\..*'");
  std::smatch cross2Match;
  REQUIRE(std::regex_search(text, cross2Match, regex3));
  REQUIRE(cross1Match.size() == 1);

  check_file_hash("test_unspec_stereo.svg");
}

TEST_CASE("Colour H light blue with no atom labels", "") {
  auto m1 = "C[C@]12CCCC[C@H]1OCCC2"_smiles;
  MolDraw2DUtils::prepareMolForDrawing(*m1);
  MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
  drawer.drawOptions().noAtomLabels = true;
  drawer.drawMolecule(*m1);
  drawer.finishDrawing();
  auto text = drawer.getDrawingText();
  std::ofstream outs("light_blue_h_no_label_1.svg");
  outs << text;
  outs.flush();
  std::regex regex1(R"(class='bond-12 atom-6 atom-11'.*fill:#ADD8E5)");
  std::smatch regex1Match;
  REQUIRE(std::regex_search(text, regex1Match, regex1));
  REQUIRE(regex1Match.size() == 1);
  check_file_hash("light_blue_h_no_label_1.svg");
}

TEST_CASE("Bond Highlights", "") {
  auto m1 = "c1c(OCC)cncc1CCCC=O"_smiles;
  REQUIRE(m1);
  MolDraw2DUtils::prepareMolForDrawing(*m1);
#if 1
  {
    // only bonds highlighted, continuous highlighting, highlights
    // joining neatly.
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    drawer.drawOptions().addBondIndices = true;
    std::vector<int> highAts{};
    std::vector<int> highBnds{0, 1, 4, 5, 6, 7, 13};
    //    std::vector<int> highBnds{0, 1, 4};
    drawer.drawMolecule(*m1, &highAts, &highBnds);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("bond_highlights_1.svg");
    outs << text;
    outs.flush();
    check_file_hash("bond_highlights_1.svg");
  }
#endif
#if 1
  {
    // same as 1, but with highlighting as coloured bonds.  The O for
    // atom 2 is red because it is not highlighted, though bond 1 from
    // the pyridyl is.
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    std::vector<int> highAts{};
    std::vector<int> highBnds{0, 1, 4, 5, 6, 7, 13};
    drawer.drawOptions().continuousHighlight = false;
    drawer.drawOptions().circleAtoms = false;
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawMolecule(*m1, &highAts, &highBnds);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("bond_highlights_2.svg");
    outs << text;
    outs.flush();
    check_file_hash("bond_highlights_2.svg");
  }
#endif
#if 1
  {
    // same bonds highlighted, but some atoms highlighted with
    // different colours.  Where an atom and a bond off it are
    // highlighted in different colours, the bond colour takes
    // precedence and the atom highlight is lost unless it has
    // an atom symbol drawn or there's a non-highlighted bond
    // off it.  Thus half of bond 8 should be green, as is
    // the N of atom 6 and the two half bonds off atom 10.
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    std::vector<int> highAts{0, 1, 5, 6, 7, 8, 10};
    std::vector<int> highBnds{0, 1, 4, 5, 6, 7, 13};
    std::map<int, DrawColour> atom_highlight_colors;
    atom_highlight_colors[0] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[5] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[6] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[7] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[8] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[10] = DrawColour(0.0, 1.0, 0.0);
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawOptions().addBondIndices = true;
    drawer.drawOptions().continuousHighlight = false;
    drawer.drawOptions().circleAtoms = false;
    drawer.drawMolecule(*m1, &highAts, &highBnds, &atom_highlight_colors);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("bond_highlights_3.svg");
    outs << text;
    outs.flush();
    check_file_hash("bond_highlights_3.svg");
  }
#endif
#if 1
  {
    // same as 3, except that the N on atom 6 isn't highlighted,
    // but both bonds off it are, so it gets the highlight colour.
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    std::vector<int> highAts{0, 1, 5, 7, 8, 10};
    std::vector<int> highBnds{0, 1, 4, 5, 6, 7, 13};
    std::map<int, DrawColour> atom_highlight_colors;
    atom_highlight_colors[0] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[5] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[7] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[8] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[10] = DrawColour(0.0, 1.0, 0.0);
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawOptions().addBondIndices = true;
    drawer.drawOptions().continuousHighlight = false;
    drawer.drawOptions().circleAtoms = false;
    drawer.drawMolecule(*m1, &highAts, &highBnds, &atom_highlight_colors);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("bond_highlights_4.svg");
    outs << text;
    outs.flush();
    check_file_hash("bond_highlights_4.svg");
  }
#endif
#if 1
  {
    // same as 4, except that atom 6 has a highlight colour assigned
    // in the map, but isn't highlighted.  It just happens that it's
    // the default highlight colour.
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    std::vector<int> highAts{0, 1, 5, 7, 8, 10};
    std::vector<int> highBnds{0, 1, 4, 5, 6, 7, 13};
    std::map<int, DrawColour> atom_highlight_colors;
    atom_highlight_colors[0] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[5] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[6] = drawer.drawOptions().highlightColour;
    atom_highlight_colors[7] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[8] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[10] = DrawColour(0.0, 1.0, 0.0);
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawOptions().addBondIndices = true;
    drawer.drawOptions().continuousHighlight = false;
    drawer.drawOptions().circleAtoms = false;
    drawer.drawMolecule(*m1, &highAts, &highBnds, &atom_highlight_colors);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("bond_highlights_5.svg");
    outs << text;
    outs.flush();
    check_file_hash("bond_highlights_5.svg");
  }
#endif
#if 1
  {
    // same as 3, but showing that atom circles can be used
    // to rescue the missing atom highlights.
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    std::vector<int> highAts{0, 1, 5, 6, 7, 8, 10};
    std::vector<int> highBnds{0, 1, 4, 5, 6, 7, 13};
    std::map<int, DrawColour> atom_highlight_colors;
    atom_highlight_colors[0] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[5] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[6] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[7] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[8] = DrawColour(0.0, 1.0, 0.0);
    atom_highlight_colors[10] = DrawColour(0.0, 1.0, 0.0);
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawOptions().addBondIndices = true;
    drawer.drawOptions().continuousHighlight = false;
    drawer.drawOptions().circleAtoms = true;
    drawer.drawMolecule(*m1, &highAts, &highBnds, &atom_highlight_colors);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("bond_highlights_6.svg");
    outs << text;
    outs.flush();
    check_file_hash("bond_highlights_6.svg");
  }
#endif
#if 1
  {
    // same as 1, but in ACS1996 mode.
    MolDraw2DSVG drawer(-1, -1, -1, -1, NO_FREETYPE);
    std::vector<int> highAts{};
    std::vector<int> highBnds{0, 1, 4, 5, 6, 7, 13};
    drawer.drawOptions().continuousHighlight = false;
    drawer.drawOptions().circleAtoms = false;
    MolDraw2DUtils::drawMolACS1996(drawer, *m1, "", &highAts, &highBnds);
    drawer.drawMolecule(*m1, &highAts, &highBnds);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("bond_highlights_7.svg");
    outs << text;
    outs.flush();
    check_file_hash("bond_highlights_7.svg");
  }
#endif
#if 1
  {
    // check 3- and 4-way intersections of continuous highlights are ok
    auto m = "c1c(C(C)(C)C)cccc1"_smiles;
    REQUIRE(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    drawer.drawOptions().continuousHighlight = true;
    drawer.drawOptions().addBondIndices = true;
    std::vector<int> highBnds{0, 1, 2, 3, 4, 5, 6};
    std::map<int, DrawColour> bond_highlight_colors;
    bond_highlight_colors[0] = DrawColour(1.0, 0.0, 0.0);
    bond_highlight_colors[1] = DrawColour(0.0, 1.0, 0.0);
    bond_highlight_colors[3] = DrawColour(1.0, 0.0, 0.0);
    bond_highlight_colors[4] = DrawColour(0.0, 1.0, 0.0);
    bond_highlight_colors[5] = DrawColour(0.0, 0.0, 1.0);
    drawer.drawMolecule(*m, nullptr, &highBnds, nullptr,
                        &bond_highlight_colors);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("bond_highlights_8.svg");
    outs << text;
    outs.flush();
    check_file_hash("bond_highlights_8.svg");
  }
#endif
#if 1
  {
    // cyclopropane (Github5592)
    auto m = "CC1CC1"_smiles;
    REQUIRE(m);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    std::vector<int> highBnds{1, 2, 3};
    std::map<int, DrawColour> bond_highlight_colors;
    bond_highlight_colors[1] = DrawColour(1.0, 0.0, 0.0);
    bond_highlight_colors[2] = DrawColour(0.0, 1.0, 0.0);
    bond_highlight_colors[3] = DrawColour(0.0, 0.0, 1.0);
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    drawer.drawOptions().continuousHighlight = true;
    drawer.drawMolecule(*m, nullptr, &highBnds, nullptr,
                        &bond_highlight_colors);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("bond_highlights_9.svg");
    outs << text;
    std::regex regex1(
        R"(class='bond-[\d*] atom-[\d] atom-[\d]' d='M ([\d.]*),([\d.]*) L ([\d.]*),([\d.]*) L ([\d.]*),([\d.]*) L ([\d.]*),([\d.]*) Z')");
    auto match_begin = std::sregex_iterator(text.begin(), text.end(), regex1);
    auto match_end = std::sregex_iterator();
    std::vector<Point2D> ends;
    for (std::sregex_iterator i = match_begin; i != match_end; ++i) {
      // with this result, match[0] is the whole string that matched,
      // match[1] is the 1st float (the x coord of the M), match[2]
      // is the 2nd float, etc.
      std::smatch match = *i;
      ends.push_back(Point2D(std::stod(match[1]), std::stod(match[2])));
      ends.push_back(Point2D(std::stod(match[3]), std::stod(match[4])));
      ends.push_back(Point2D(std::stod(match[5]), std::stod(match[6])));
      ends.push_back(Point2D(std::stod(match[7]), std::stod(match[8])));
    }
    REQUIRE(ends.size() == 12);
    // When this had a bug in it, it drew a butterfly-type motif, because
    // ends 2 and 3 were the wrong way round.
    for (int i = 0; i < 12; i += 4) {
      REQUIRE(!MolDraw2D_detail::doLinesIntersect(
          ends[i], ends[i + 3], ends[i + 1], ends[i + 2], nullptr));
    }
    outs.flush();
    check_file_hash("bond_highlights_9.svg");
  }
#endif
}

TEST_CASE("drawMolecules should not crash on null molecules",
          "[drawing][bug]") {
  auto m1 = "c1ccccc1"_smiles;
  auto m2 = "c1ccncc1"_smiles;
  REQUIRE(m1);
  REQUIRE(m2);
  MolDraw2DSVG drawer(1000, 200, 100, 100, NO_FREETYPE);
  RWMol dm1(*m1);
  RWMol dm2(*m2);
  MOL_PTR_VECT ms{&dm1,    nullptr, nullptr, nullptr, nullptr, nullptr,
                  nullptr, nullptr, nullptr, nullptr, &dm2};
  drawer.drawMolecules(ms);
  drawer.finishDrawing();
  auto text = drawer.getDrawingText();
  std::regex regex1("<path d=");
  auto nMatches =
      std::distance(std::sregex_iterator(text.begin(), text.end(), regex1),
                    std::sregex_iterator());
  REQUIRE(nMatches == 11);
}

TEST_CASE("Crossed bonds in transdecene") {
  SECTION("basics") {
    std::string nameBase = "testGithub5486_";
    {
      auto m = R"CTAB(
  ChemDraw08042214332D

 10 10  0  0  0  0  0  0  0  0999 V2000
   -1.4289    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4289   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7145    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145   -0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4289   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4289    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  5  6  2  0
  6  1  1  0
  4  7  1  0
  7  8  1  0
  8  9  1  0
  9 10  1  0
 10  5  1  0
M  END)CTAB"_ctab;
      REQUIRE(m);
      {
        MolDraw2DSVG drawer(300, 300);
        drawer.drawMolecule(*m);
        drawer.finishDrawing();
        std::string text = drawer.getDrawingText();
        std::ofstream outs(nameBase + "1.svg");
        outs << text;
        outs.flush();
        outs.close();
        std::regex regex1(
            R"(class='bond-3 atom-4 atom-5' d='M ([\d.]*),([\d.]*) L ([\d.]*),([\d.]*)')");
        auto match_begin =
            std::sregex_iterator(text.begin(), text.end(), regex1);
        auto match_end = std::sregex_iterator();
        std::vector<Point2D> ends;
        for (std::sregex_iterator i = match_begin; i != match_end; ++i) {
          std::smatch match = *i;
          ends.push_back(Point2D(std::stod(match[1]), std::stod(match[2])));
          ends.push_back(Point2D(std::stod(match[3]), std::stod(match[4])));
        }
        REQUIRE(ends.size() == 4);
        REQUIRE(!MolDraw2D_detail::doLinesIntersect(ends[0], ends[1], ends[2],
                                                    ends[3], nullptr));
        check_file_hash(nameBase + "1.svg");
      }
    }
  }
}

TEST_CASE("Bad O position in aldehydes", "") {
  std::string nameBase("testGithub5511_");
  {
    auto m1 = "O=Cc1cccc(C=O)c1"_smiles;
    REQUIRE(m1);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::string filename = nameBase + "1.svg";
    std::ofstream outs(filename);
    outs << text;
    outs.flush();
    check_file_hash(filename);
  }
  {
    auto m = R"CTAB(
     RDKit          2D

 11 11  0  0  0  0  0  0  0  0999 V2000
   -4.2885    0.5445    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9895    1.2945    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6818    0.5394    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3738    1.2945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9342    0.5394    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2422    1.2945    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2422    2.7945    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9342   -0.9708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3738   -1.7262    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6818   -0.9708    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9895   -1.7262    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  1  0
  6  7  2  0
  5  8  2  0
  8  9  1  0
  9 10  2  0
 10 11  1  0
 10  3  1  0
M  END
)CTAB"_ctab;
    REQUIRE(m);
    MolDraw2DSVG drawer(250, 250, -1, -1);
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::string filename = nameBase + "2.svg";
    std::ofstream outs(filename);
    outs << text;
    outs.flush();
    check_file_hash(filename);
  }
}

TEST_CASE("Github5534") {
  std::string nameBase = "test_github_5534";
  {
    auto m = R"CTAB(
  INFOCHEM          2D 1   1.00000     0.00000     0

 45 45  0  0  0  0  0  0  0  0999 V2000
   -1.3388    1.6804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0303    1.6804    0.0000 Si  0  0  0  0  0  0  0  0  0  0  0  0
   -1.0303    1.9890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0303    1.3609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0496   -0.6997    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0496   -0.3581    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0496   -0.0386    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0496    0.2810    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0496    0.6116    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0826    0.9532    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0496   -1.0193    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0496   -1.3499    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0496   -1.6694    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.6116   -0.6997    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.6116   -0.3581    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.6116   -0.0386    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.6116    0.2810    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.6116    0.6116    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.6116    0.9532    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.6116   -1.0193    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.6116   -1.3499    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.6116   -1.6694    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.0303    1.6804    0.0000 Si  0  0  0  0  0  0  0  0  0  0  0  0
    1.0303    2.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0303    1.3609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3388    1.6804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2810   -0.6997    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2810   -0.3581    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2810   -0.0386    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2810    0.2810    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2810    0.6116    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2810    0.9532    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3361    1.2948    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1708    1.4601    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3581    1.6804    0.0000 Si  0  0  0  0  0  0  0  0  0  0  0  0
    0.0275    1.6804    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3581    1.9890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2810   -1.0193    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2810   -1.3499    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2810   -1.6694    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2810   -2.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6997    1.6804    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3691    1.6804    0.0000 Si  0  0  0  0  0  0  0  0  0  0  0  0
   -0.3691    1.9890    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3691    1.3609    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1 23  1  0  0  0  0
  1 35  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2 42  1  0  0  0  0
  5 27  1  0  0  0  0
  6 28  1  0  0  0  0
  7 29  1  0  0  0  0
  8 30  1  0  0  0  0
  9 31  1  0  0  0  0
 10 32  1  0  0  0  0
 11 38  1  0  0  0  0
 12 39  1  0  0  0  0
 13 40  1  0  0  0  0
 14 27  1  0  0  0  0
 15 28  1  0  0  0  0
 16 29  1  0  0  0  0
 17 30  1  0  0  0  0
 18 31  1  0  0  0  0
 19 32  1  0  0  0  0
 20 38  1  0  0  0  0
 21 39  1  0  0  0  0
 22 40  1  0  0  0  0
 23 24  1  0  0  0  0
 23 25  1  0  0  0  0
 23 26  1  0  0  0  0
 27 28  1  0  0  0  0
 27 38  1  0  0  0  0
 28 29  1  0  0  0  0
 29 30  1  0  0  0  0
 30 31  1  0  0  0  0
 31 32  1  0  0  0  0
 32 33  1  0  0  0  0
 33 34  1  0  0  0  0
 34 35  1  0  0  0  0
 35 36  1  0  0  0  0
 35 37  1  0  0  0  0
 36 43  1  0  0  0  0
 38 39  1  0  0  0  0
 39 40  1  0  0  0  0
 40 41  1  0  0  0  0
 42 43  1  0  0  0  0
 43 44  1  0  0  0  0
 43 45  1  0  0  0  0
M  STY  2   1 GEN   2 GEN
M  SAL   1 15  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41
M  SBL   1 15   7  16   8  17   9  18  10  19  11  20  12  21   3  39  13
M  SBL   1  5  22  14  23  15  24
M  SDI   1  4   -0.0386    2.0661   -0.0386    1.1846
M  SDI   1  4    0.4683    2.0661    0.4683    1.1846
M  SAL   2  4  42  43  44  45
M  SBL   2  2   6  39
M  SDI   2  4   -0.7658    2.0661   -0.7658    1.1846
M  SDI   2  4   -0.2149    2.0661   -0.2149    1.1846
M  END
)CTAB"_ctab;
    REQUIRE(m);
    {
      MolDraw2DSVG drawer(300, 300, 300, 300, NO_FREETYPE);
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs(nameBase + ".svg");
      outs << text;
      outs.flush();
      outs.close();
      check_file_hash(nameBase + ".svg");
    }
  }
}

TEST_CASE(
    "Github5704: set bond highlight color when atom highlight color changes") {
  std::string nameBase = "test_github5704";

  auto m =
      "CCCO |(-1.97961,-0.1365,;-0.599379,0.450827,;0.599379,-0.450827,;1.97961,0.1365,)|"_smiles;
  REQUIRE(m);
  std::regex redline(R"RE(<path .*fill:#FF7F7F)RE");
  std::regex blueline(R"RE(<path .*fill:#4C4CFF)RE");

  SECTION("no atom colors specified, default behavior") {
    MolDraw2DSVG drawer(300, 300, 300, 300, NO_FREETYPE);
    std::vector<int> aids{0, 1, 2};
    drawer.drawMolecule(*m, "red bond highlight", &aids);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();

    std::smatch rematch;
    CHECK(std::regex_search(text, rematch, redline));
    CHECK(!std::regex_search(text, rematch, blueline));

    std::ofstream outs(nameBase + "_1.svg");
    outs << text;
    outs.flush();
    outs.close();
  }

  SECTION("both ends specified") {
    MolDraw2DSVG drawer(300, 300, 300, 300, NO_FREETYPE);
    std::vector<int> aids{0, 1, 2};
    std::map<int, DrawColour> acolors{
        {0, {.3, .3, 1}}, {1, {.3, .3, 1}}, {2, {.3, .3, 1}}};
    drawer.drawMolecule(*m, "blue bond highlight", &aids, &acolors);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();

    std::smatch rematch;
    CHECK(!std::regex_search(text, rematch, redline));
    CHECK(std::regex_search(text, rematch, blueline));

    std::ofstream outs(nameBase + "_2.svg");
    outs << text;
    outs.flush();
    outs.close();
  }

  SECTION("color just on begin") {
    MolDraw2DSVG drawer(300, 300, 300, 300, NO_FREETYPE);
    std::vector<int> aids{0, 1};
    std::map<int, DrawColour> acolors{{0, {.3, .3, 1}}};
    drawer.drawMolecule(*m, "blue bond highlight", &aids, &acolors);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();

    std::smatch rematch;
    CHECK(!std::regex_search(text, rematch, redline));
    CHECK(std::regex_search(text, rematch, blueline));

    std::ofstream outs(nameBase + "_3.svg");
    outs << text;
    outs.flush();
    outs.close();
  }
  SECTION("color just on end") {
    MolDraw2DSVG drawer(300, 300, 300, 300, NO_FREETYPE);
    std::vector<int> aids{0, 1};
    std::map<int, DrawColour> acolors{{1, {.3, .3, 1}}};
    drawer.drawMolecule(*m, "blue bond highlight", &aids, &acolors);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();

    std::smatch rematch;

    CHECK(!std::regex_search(text, rematch, redline));
    CHECK(std::regex_search(text, rematch, blueline));

    std::ofstream outs(nameBase + "_4.svg");
    outs << text;
    outs.flush();
    outs.close();
  }
}

TEST_CASE("Github5767: monomer label missing for MON SGroups ") {
  std::string nameBase = "test_github5767";
  auto m = R"CTAB(
  Marvin  06091012252D

 13 11  0  0  0  0            999 V2000
   -3.5063    2.1509    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7918    2.5634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0773    2.1509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3628    2.5634    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6484    2.1509    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0661    2.5634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7806    2.1509    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9984   -0.4714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2839   -0.0589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5695   -0.4714    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8550   -0.0589    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -2.9984    0.3536    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.6777   -0.1775    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  7  1  0  0  0  0
  8  9  1  0  0  0  0
  9 10  1  0  0  0  0
 10 11  1  0  0  0  0
  9 12  1  0  0  0  0
 12  8  1  0  0  0  0
M  STY  2   1 MON   2 MON
M  SAL   1  7   1   2   3   4   5   6   7
M  SDI   1  4   -3.9263    1.7309   -3.9263    2.9834
M  SDI   1  4    1.2006    2.9834    1.2006    1.7309
M  SAL   2  5   8   9  10  11  12
M  SDI   2  4   -3.4184   -0.8914   -3.4184    0.7736
M  SDI   2  4   -0.4350    0.7736   -0.4350   -0.8914
M  END
)CTAB"_ctab;
  REQUIRE(m);
  {
    MolDraw2DSVG drawer(300, 300, 300, 300, true);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs(nameBase + ".svg");
    outs << text;
    outs.flush();
    outs.close();
    // there should be 2 each of " >[mon]</text>"
    std::vector<std::string> needed{" >m</text>", " >o</text>", " >n</text>"};
    for (const auto &n : needed) {
      std::regex rn(n);
      std::ptrdiff_t const match_count(
          std::distance(std::sregex_iterator(text.begin(), text.end(), rn),
                        std::sregex_iterator()));
      REQUIRE(match_count == 2);
    }
    check_file_hash(nameBase + ".svg");
  }
}
