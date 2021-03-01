//
//  Copyright (C) 2019-2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>

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
    outs.flush();

    CHECK(text.find("<circle") != std::string::npos);
    CHECK(text.find("<circle") != std::string::npos);
    CHECK(text.find("atom-selector") != std::string::npos);
    CHECK(text.find("bond-selector") != std::string::npos);
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
    outs.flush();
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
    outs.flush();
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
    outs.flush();
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
    outs.flush();
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
    outs.flush();

    CHECK(text.find("<path class='bond-0 atom-0 atom-1' d='M 126.052,100 L 85.9675,100'"
                    " style='fill:none;fill-rule:evenodd;"
                    "stroke:#0000FF;") != std::string::npos);
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
    outs.flush();

    CHECK(text.find("<path class='bond-7 atom-7 atom-8' d='M 101.307,79.424 "
                    "L 95.669,87.1848' style='fill:none;"
                    "fill-rule:evenodd;stroke:#0000FF;") != std::string::npos);
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
    outs.flush();

    CHECK(text.find("<path class='bond-2 atom-3 atom-4' d='M 53.289,140.668"
                    " L 81.0244,149.68' style='fill:none;"
                    "fill-rule:evenodd;stroke:#0000FF;") != std::string::npos);
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
      outs.flush();
    }
    {
      MolDraw2DSVG drawer(250, 250, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareMolForDrawing(*m1);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testDativeBonds_2b.svg");
      outs << text;
      outs.flush();
    }
    {
      MolDraw2DSVG drawer(350, 350, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareMolForDrawing(*m1);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testDativeBonds_2c.svg");
      outs << text;
      outs.flush();
    }
    {
      MolDraw2DSVG drawer(450, 450, -1, -1, NO_FREETYPE);
      MolDraw2DUtils::prepareMolForDrawing(*m1);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testDativeBonds_2d.svg");
      outs << text;
      outs.flush();
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
    outs.flush();

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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
    outs.flush();

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
  REQUIRE(m1);
  SECTION("foundations") {
    {
      MolDraw2DSVG drawer(250, 200, -1, -1, NO_FREETYPE);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_1.svg");
      outs << text;
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
      CHECK(text.find(">1</text>") != std::string::npos);
      CHECK(text.find(">,</text>") != std::string::npos);
      CHECK(text.find(">(</text>") != std::string::npos);
      CHECK(text.find(">S</text>") != std::string::npos);
      CHECK(text.find(")</text>") != std::string::npos);
      CHECK(text.find(">2</text>") != std::string::npos);
      CHECK(text.find(">f</text>") != std::string::npos);
      CHECK(text.find(">o</text>") != std::string::npos);
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
      outs.flush();
      std::vector<std::string> tkns;
      boost::algorithm::find_all(tkns, text, "bond-0");
      CHECK(tkns.size() == 6);
    }
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  SECTION("larger PNG") {
    {
      MolDraw2DCairo drawer(450, 400);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      drawer.writeDrawingText("testGithub3226_1.png");
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
      outs.flush();
      std::vector<std::string> tkns;
      boost::algorithm::find_all(tkns, text, "bond-0");
      CHECK(tkns.size() == 4);
    }
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  SECTION("smaller PNG") {
    {
      MolDraw2DCairo drawer(200, 150);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      drawer.writeDrawingText("testGithub3226_2.png");
    }
  }
#endif
  SECTION("middle SVG") {
    {
      MolDraw2DSVG drawer(250, 200);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub3226_3.svg");
      outs << text;
      outs.flush();
      std::vector<std::string> tkns;
      boost::algorithm::find_all(tkns, text, "bond-0");
      CHECK(tkns.size() == 4);
    }
  }
#ifdef RDK_BUILD_CAIRO_SUPPORT
  SECTION("middle PNG") {
    {
      MolDraw2DCairo drawer(250, 200);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      drawer.writeDrawingText("testGithub3226_3.png");
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
    outs.flush();
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
      outs.flush();
      CHECK(text.find("<path d='M") != std::string::npos);
    }
    {
      MolDraw2DSVG drawer(250, 200, panelWidth, panelHeight, noFreeType);
      drawer.drawOptions().includeRadicals = false;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testIncludeRadicals_1b.svg");
      outs << text;
      outs.flush();
      CHECK(text.find("<path d='M") == std::string::npos);
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
    outs.flush();

    // make sure the polygon starts at a bond
    CHECK(text.find("<path class='bond-0 atom-0 atom-1' d='M 321.962,140") !=
          std::string::npos);
    CHECK(text.find("<path d='M 321.962,140") != std::string::npos);
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
    outs.flush();
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
      outs.flush();
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(450, 400);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m, "Oxytocin (flat)");
      drawer.finishDrawing();
      drawer.writeDrawingText("testHandDrawn-1.png");
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
      outs.flush();
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(450, 400);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m, "Oxytocin");
      drawer.finishDrawing();
      drawer.writeDrawingText("testHandDrawn-2.png");
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
      outs.flush();
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(350, 300);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText("testHandDrawn-3.png");
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
      outs.flush();
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(350, 300);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText("testHandDrawn-4.png");
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
      outs.flush();
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
      outs.flush();
    }
#ifdef RDK_BUILD_CAIRO_SUPPORT
    {
      MolDraw2DCairo drawer(900, 450);
      drawer.drawOptions().fontFile = fName;
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      drawer.writeDrawingText("testHandDrawn-5.png");
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
      outs.flush();
    }
    {  // rotation
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().rotate = 90;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-1b.svg");
      outs << text;
      outs.flush();
    }
    {  // centering
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().centreMoleculesBeforeDrawing = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-1c.svg");
      outs << text;
      outs.flush();
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
      outs.flush();
    }
    {  // rotation
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().rotate = 180;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-1e.svg");
      outs << text;
      outs.flush();
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
      outs.flush();
    }
    {  // rotation
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().rotate = 90;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-2b.svg");
      outs << text;
      outs.flush();
    }
    {  // centering
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().centreMoleculesBeforeDrawing = true;
      drawer.drawMolecule(*m);
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testBrackets-2c.svg");
      outs << text;
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
    }
    {  // make sure comic mode doesn't screw this up
      MolDraw2DSVG drawer(350, 300);
      drawer.drawOptions().comicMode = true;
      drawer.drawMolecule(*m, "comic variations");
      drawer.finishDrawing();
      auto text = drawer.getDrawingText();
      std::ofstream outs("testPositionVariation-1b.svg");
      outs << text;
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
    }
  }
}

TEST_CASE("disable atom labels", "[feature]") {
  SECTION("basics") {
    auto m = "NCC(=O)O"_smiles;
    MolDraw2DSVG drawer(350, 300);
    MolDraw2DUtils::prepareMolForDrawing(*m);
    drawer.drawOptions().noAtomLabels = true;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs("testNoAtomLabels-1.svg");
    outs << text;
    outs.flush();
    CHECK(text.find("class='atom-0") == std::string::npos);
    CHECK(text.find("class='atom-3") == std::string::npos);
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
    std::vector<int> rotns={0,30,60,90,120,150,180};
    for(auto rotn : rotns){
    MolDraw2DSVG drawer(350, 300);
    drawer.drawOptions().rotate = (double)rotn;
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    auto text = drawer.getDrawingText();
    std::ofstream outs((boost::format("testLinkNodes-2-%d.svg")%rotn).str());
    outs << text;
    outs.flush();
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
    outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
    outs.flush();
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
    outs.flush();
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
      outs.flush();
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
    outs.flush();
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
    // least 1.3, otherwise the inner bonds are not actually inside the ring.
    float outerBondsDistance = (bond0OuterCtd - bond2OuterCtd).length();
    float innerBondsDistance = (bond0InnerCtd - bond2InnerCtd).length();
    CHECK(outerBondsDistance / innerBondsDistance > 1.3f);
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
    outs.flush();
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
    outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
      outs.flush();
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
    outs.flush();
    size_t nDeuteriumTritium = std::distance(
        std::sregex_token_iterator(textDeuteriumTritium.begin(),
                                   textDeuteriumTritium.end(), regex, 1),
        std::sregex_token_iterator());
    CHECK(nDeuteriumTritium == 2);
  }
}
