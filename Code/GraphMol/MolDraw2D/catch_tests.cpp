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
#include <GraphMol/FileParsers/FileParsers.h>
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

    CHECK(text.find("<path class='bond-0' d='M 55.1495,101.204"
                    " L 52.2436,100 L 55.1495,98.7964") != std::string::npos);
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

    CHECK(text.find("<path class='bond-7' d='M 92.9861,93.568"
                    " L 92.2758,94.0033 L 92.4703,93.1932") !=
          std::string::npos);
  }
  SECTION("test colours") {
    // the dative bonds point the wrong way, but the point is to test
    // if the tip of the arrow is blue.
    auto m1 = "[Cu++]->1->2.N1CCN2"_smiles;
    REQUIRE(m1);
    MolDraw2DSVG drawer(200, 200);
    MolDraw2DUtils::prepareMolForDrawing(*m1);
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs("testDativeBonds_3.svg");
    outs << text;
    outs.flush();

    CHECK(text.find("<path class='bond-2' d='M 57.7741,143.825"
                    " L 85.7826,152.925' style='fill:none;"
                    "fill-rule:evenodd;stroke:#0000FF") != std::string::npos);
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

TEST_CASE("bad DrawMolecules() when molecules are not kekulized",
          "[drawing,bug]") {
  auto m1 = "CCN(CC)CCn1nc2c3ccccc3sc3c(CNS(C)(=O)=O)ccc1c32"_smiles;
  REQUIRE(m1);
  SECTION("foundations") {
    MolDraw2DSVG drawer(500, 200, 250, 200);
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
      MolDraw2DSVG drawer(250, 200);
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_1.svg");
      outs << text;
      outs.flush();
      CHECK(text.find("<tspan>1</tspan>") == std::string::npos);
      CHECK(text.find("<tspan>1,(S)</tspan>") == std::string::npos);
    }
    {
      MolDraw2DSVG drawer(250, 200);
      drawer.drawOptions().addAtomIndices = true;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_2.svg");
      outs << text;
      outs.flush();
      CHECK(text.find("<tspan>1</tspan>") != std::string::npos);
      // it only appears once though:
      CHECK(text.find("<tspan>1</tspan>", text.find("<tspan>1</tspan>") + 1) ==
            std::string::npos);
      CHECK(text.find("1,(S)") == std::string::npos);
    }
    {
      MolDraw2DSVG drawer(250, 200);
      drawer.drawOptions().addBondIndices = true;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_3.svg");
      outs << text;
      outs.flush();
      CHECK(text.find("<tspan>1</tspan>") != std::string::npos);
      // it only appears once though:
      CHECK(text.find("<tspan>1</tspan>", text.find("<tspan>1</tspan>") + 1) ==
            std::string::npos);
    }
    {
      MolDraw2DSVG drawer(250, 200);
      drawer.drawOptions().addAtomIndices = true;
      drawer.drawOptions().addBondIndices = true;
      drawer.drawMolecule(*m1);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testAtomBondIndices_4.svg");
      outs << text;
      outs.flush();
      CHECK(text.find("<tspan>1</tspan>") != std::string::npos);
      // it appears twice:
      CHECK(text.find("<tspan>1</tspan>", text.find("<tspan>1</tspan>") + 1) !=
            std::string::npos);
    }
    {
      MolDraw2DSVG drawer(250, 200);
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
      CHECK(text.find("<tspan>1,(S)</tspan>") != std::string::npos);
      CHECK(text.find("<tspan>2,foo</tspan>") != std::string::npos);
      CHECK(text.find("<tspan>1</tspan>") == std::string::npos);
    }
  }
}

TEST_CASE("Github #3126: DrawMolecules does not center molecules",
          "[drawing][bug]") {
  SECTION("basics") {
    auto m1 = R"CTAB(
  Mrv2007 05012010452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O -4.2917 5.3517 0 0
M  V30 2 C -5.5375 4.4464 0 0
M  V30 3 C -5.0617 2.9819 0 0
M  V30 4 C -3.5217 2.9819 0 0
M  V30 5 C -3.0458 4.4464 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
      )CTAB"_ctab;
    REQUIRE(m1);
    auto m2 = R"CTAB(
  Mrv2007 05012010452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 S -3.2917 6.3517 0 0
M  V30 2 C -4.5375 5.4464 0 0
M  V30 3 C -4.0617 3.9819 0 0
M  V30 4 C -2.5217 3.9819 0 0
M  V30 5 C -2.0458 5.4464 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
      )CTAB"_ctab;
    REQUIRE(m2);
    {
      MolDraw2DSVG drawer(500, 200, 250, 200);
      drawer.drawOptions().prepareMolsBeforeDrawing = false;
      drawer.drawOptions().centreMoleculesBeforeDrawing = true;
      MOL_PTR_VECT ms{m1.get(), m2.get()};
      drawer.drawMolecules(ms);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub3126_1.svg");
      outs << text;
      outs.flush();

      // should have two atom labels drawn at the same height:
      std::string tgt = "y='31.6464' style='font-size";
      auto idx1 = text.find(tgt);
      CHECK(idx1 != std::string::npos);
      auto idx2 = text.find(tgt, idx1 + 1);
      CHECK(idx2 != std::string::npos);
    }
    {
      // if we don't center we get different results:
      MolDraw2DSVG drawer(500, 200, 250, 200);
      drawer.drawOptions().prepareMolsBeforeDrawing = false;
      drawer.drawOptions().centreMoleculesBeforeDrawing = false;
      MOL_PTR_VECT ms{m1.get(), m2.get()};
      drawer.drawMolecules(ms);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub3126_2.svg");
      outs << text;
      outs.flush();

      // should NOT have two atom labels drawn at the same height:
      std::string tgt = "y='75.0166' style='font-size";
      auto idx1 = text.find(tgt);
      CHECK(idx1 != std::string::npos);
      auto idx2 = text.find(tgt, idx1 + 1);
      CHECK(idx2 == std::string::npos);
    }
  }

  SECTION("more complex") {
    auto m1 = R"CTAB(
  Mrv2007 05012010452D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O -4.2917 5.3517 0 0
M  V30 2 C -5.5375 4.4464 0 0
M  V30 3 C -5.0617 2.9819 0 0
M  V30 4 C -3.5217 2.9819 0 0
M  V30 5 C -3.0458 4.4464 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
      )CTAB"_ctab;
    REQUIRE(m1);
    auto m2 = R"CTAB(
  Mrv2007 05012011202D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 16 16 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 0.875 5.3517 0 0
M  V30 2 C -0.3709 4.4464 0 0
M  V30 3 C 0.105 2.9819 0 0
M  V30 4 C 1.645 2.9819 0 0
M  V30 5 C 2.1209 4.4464 0 0
M  V30 6 C -0.8002 1.736 0 0
M  V30 7 C -0.1738 0.3291 0 0
M  V30 8 C -1.079 -0.9167 0 0
M  V30 9 C -0.4526 -2.3236 0 0
M  V30 10 C -1.3578 -3.5695 0 0
M  V30 11 C -0.7314 -4.9763 0 0
M  V30 12 C -1.6366 -6.2222 0 0
M  V30 13 C -1.0103 -7.6291 0 0
M  V30 14 C -1.9154 -8.875 0 0
M  V30 15 C -1.2891 -10.2818 0 0
M  V30 16 C -2.1943 -11.5277 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 1 5
M  V30 6 1 3 6
M  V30 7 1 6 7
M  V30 8 1 7 8
M  V30 9 1 8 9
M  V30 10 1 9 10
M  V30 11 1 10 11
M  V30 12 1 11 12
M  V30 13 1 12 13
M  V30 14 1 13 14
M  V30 15 1 14 15
M  V30 16 1 15 16
M  V30 END BOND
M  V30 END CTAB
M  END
)CTAB"_ctab;
    REQUIRE(m2);
    {
      MolDraw2DSVG drawer(500, 200, 250, 200);
      drawer.drawOptions().prepareMolsBeforeDrawing = false;
      drawer.drawOptions().centreMoleculesBeforeDrawing = true;
      MOL_PTR_VECT ms{m1.get(), m2.get()};
      drawer.drawMolecules(ms);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub3126_3.svg");
      outs << text;
      outs.flush();

      auto idx1 = text.find("y='78.0281' style='font-size");
      CHECK(idx1 != std::string::npos);
      auto idx2 = text.find("y='12.5406' style='font-size", idx1 + 1);
      CHECK(idx2 != std::string::npos);
    }
    {
      MolDraw2DSVG drawer(500, 200, 250, 200);
      drawer.drawOptions().prepareMolsBeforeDrawing = false;
      drawer.drawOptions().centreMoleculesBeforeDrawing = false;
      MOL_PTR_VECT ms{m1.get(), m2.get()};
      drawer.drawMolecules(ms);
      drawer.finishDrawing();
      std::string text = drawer.getDrawingText();
      std::ofstream outs("testGithub3126_4.svg");
      outs << text;
      outs.flush();

      auto idx1 = text.find("y='12.5406' style='font-size");
      CHECK(idx1 != std::string::npos);
      auto idx2 = text.find("y='12.5406' style='font-size", idx1 + 1);
      CHECK(idx2 != std::string::npos);
    }
  }
}
