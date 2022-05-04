//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/Qt/MolDraw2DQt.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Depictor/RDDepictor.h>

#include <QImage>
#include <QPainter>
#include <QGuiApplication>

using namespace RDKit;

TEST_CASE("basic generate PNGs", "[drawing][Qt]") {
  SECTION("with freetype") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);
    QImage qimg(250, 200, QImage::Format_RGB32);
    QPainter qpt(&qimg);
    MolDraw2DQt drawer(qimg.width(), qimg.height(), &qpt);
    MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
    qimg.save("qttest-1a.png");
  }
  SECTION("no freetype") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);
    QImage qimg(250, 200, QImage::Format_RGB32);
    QPainter qpt(&qimg);
    bool no_freetype = true;
    MolDraw2DQt drawer(qimg.width(), qimg.height(), &qpt, -1, -1, no_freetype);
    MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
    qimg.save("qttest-1b.png");
  }
}

TEST_CASE("Github #4764") {
  SECTION("basics") {
    auto mol = "c1ccccc1-C1CCCCC1"_smiles;
    REQUIRE(mol);
    std::vector<int> highlights{6, 7, 8, 9, 10, 11};
    {
      QImage qimg(200, 150, QImage::Format_RGB32);
      QPainter qpt(&qimg);
      MolDraw2DQt drawer(qimg.width(), qimg.height(), &qpt);
      drawer.drawMolecule(*mol, "highlight", &highlights);
      qimg.save("testGithub4764.qt.sz1.png");
    }
    {
      QImage qimg(400, 350, QImage::Format_RGB32);
      QPainter qpt(&qimg);
      MolDraw2DQt drawer(qimg.width(), qimg.height(), &qpt);
      drawer.drawMolecule(*mol, "highlight", &highlights);
      qimg.save("testGithub4764.qt.sz2.png");
    }
    {
      QImage qimg(800, 700, QImage::Format_RGB32);
      QPainter qpt(&qimg);
      MolDraw2DQt drawer(qimg.width(), qimg.height(), &qpt);
      drawer.drawMolecule(*mol, "highlight", &highlights);
      qimg.save("testGithub4764.qt.sz3.png");
    }
  }
}

TEST_CASE("Github #5122: bad highlighting when bondLineWidth is increased") {
  SECTION("basics") {
    auto mol = "CC1=CC=C(C=C1)C1=CC=CC=C1"_smiles;
    REQUIRE(mol);
    std::vector<int> highlightAtoms{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    std::vector<int> highlightBonds{1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13};

    {
      QImage qimg(400, 400, QImage::Format_ARGB32);
      QPainter qpt(&qimg);
      MolDraw2DQt drawer(qimg.width(), qimg.height(), &qpt);
      drawer.drawOptions().bondLineWidth = 3;
      drawer.drawOptions().scaleHighlightBondWidth = true;
      drawer.drawMolecule(*mol, "highlight both", &highlightAtoms,
                          &highlightBonds);
      qimg.save("testGithub5122.qt.1.png");
    }
    {
      QImage qimg(400, 400, QImage::Format_ARGB32);
      QPainter qpt(&qimg);
      MolDraw2DQt drawer(qimg.width(), qimg.height(), &qpt);
      drawer.drawOptions().bondLineWidth = 3;
      drawer.drawOptions().scaleHighlightBondWidth = true;
      drawer.drawMolecule(*mol, "highlight bonds", nullptr, &highlightBonds);
      qimg.save("testGithub5122.qt.2.png");
    }
  }
}

int main(int argc, char* argv[]) {
  QGuiApplication app(argc, argv);

  int result = Catch::Session().run(argc, argv);

  return result;
}