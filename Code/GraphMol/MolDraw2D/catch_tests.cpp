//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>

using namespace RDKit;

TEST_CASE("prepareAndDrawMolecule", "[drawing]") {
  SECTION("basics") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);

    // we will be able to recognize that the prep worked because there
    // will be an H in the output:
    MolDraw2DSVG drawer(200, 200);
    MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    CHECK(text.find("<tspan>H</tspan>") != std::string::npos);
  }
}

TEST_CASE("tag atoms in SVG", "[drawing, SVG]") {
  SECTION("basics") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);

    MolDraw2DSVG drawer(200, 200);
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
TEST_CASE("contour data", "[drawing, conrec]") {
  auto m1 = "C1N[C@@H]2OCC12"_smiles;
  REQUIRE(m1);
  SECTION("grid basics") {
    MolDraw2DSVG drawer(250, 250);
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
    MolDraw2DUtils::contourAndDrawGrid(drawer, grid, xps, yps, 10, levels);
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
    MolDraw2DSVG drawer(250, 250);
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
    MolDraw2DUtils::contourAndDrawGaussians(drawer, cents, weights, widths, 10,
                                            levels);

    drawer.drawOptions().clearBackground = false;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("contourMol_2.svg");
    outs << text;
    outs.flush();
  }
  SECTION("gaussian fill") {
    MolDraw2DSVG drawer(250, 250);
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
                                            levels, cps);

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

    MolDraw2DSVG drawer(450, 250);
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
                                            levels, cps);

    drawer.drawOptions().clearBackground = false;
    drawer.drawMolecule(*m2);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("contourMol_4.svg");
    outs << text;
    outs.flush();
  }
}

TEST_CASE("dative bonds", "[drawing, organometallics]") {
  SECTION("basics") {
    auto m1 = "N->[Pt]"_smiles;
    REQUIRE(m1);
    MolDraw2DSVG drawer(200, 200);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testDativeBonds_1.svg");
    outs << text;
    outs.flush();

    CHECK(text.find("<path class='bond-0' d='M 55.2063,102.675"
                    " L 48.7489,100 L 55.2063,97.3253")
          != std::string::npos);
  }
  SECTION("more complex") {
    auto m1 = "N->1[C@@H]2CCCC[C@H]2N->[Pt]11OC(=O)C(=O)O1"_smiles;
    REQUIRE(m1);
    MolDraw2DSVG drawer(200, 200);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testDativeBonds_2.svg");
    outs << text;
    outs.flush();

    CHECK(text.find("<path class='bond-7' d='M 93.3102,93.7849"
                    " L 91.7317,94.7522 L 92.1639,92.952") !=
          std::string::npos);
  }
}

TEST_CASE("zero-order bonds", "[drawing, organometallics]") {
  SECTION("basics") {
    auto m1 = "N-[Pt]"_smiles;
    REQUIRE(m1);
    m1->getBondWithIdx(0)->setBondType(Bond::ZERO);
    MolDraw2DSVG drawer(200, 200);
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
      MolDraw2DSVG drawer(200, 200);
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      CHECK(text.find("fill:#0000FF' ><tspan>HN") != std::string::npos);
    }
    {
      MolDraw2DSVG drawer(200, 200);
      assignBWPalette(drawer.drawOptions().atomColourPalette);
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      CHECK(text.find("fill:#0000FF' ><tspan>HN") == std::string::npos);
      CHECK(text.find("fill:#000000' ><tspan>HN") != std::string::npos);
    }
  }
  SECTION("test") {
    {
      MolDraw2DSVG drawer(200, 200);
      MolDrawOptions options = drawer.drawOptions();
      assignBWPalette(options.atomColourPalette);
      drawer.drawOptions() = options;
      MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      CHECK(text.find("fill:#0000FF' ><tspan>HN") == std::string::npos);
      CHECK(text.find("fill:#000000' ><tspan>HN") != std::string::npos);
    }
  }
}
