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
#include <catch2/catch_all.hpp>

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

TEST_CASE("wavy bonds and Qt") {
  {
    auto mol = "CC=CC"_smiles;
    REQUIRE(mol);
    mol->getBondWithIdx(1)->setStereoAtoms(0, 3);
    mol->getBondWithIdx(1)->setStereo(Bond::BondStereo::STEREOANY);
    bool kekulize = true;
    bool addChiralHs = true;
    bool wedgeBonds = true;
    bool forceCoords = true;
    bool wavyBonds = true;
    MolDraw2DUtils::prepareMolForDrawing(*mol, kekulize, addChiralHs,
                                         wedgeBonds, forceCoords, wavyBonds);
    CHECK(mol->getBondWithIdx(0)->getBondDir() == Bond::BondDir::UNKNOWN);
    CHECK(mol->getBondWithIdx(1)->getStereo() == Bond::BondStereo::STEREONONE);

    {
      QImage qimg(400, 400, QImage::Format_ARGB32);
      QPainter qpt(&qimg);
      MolDraw2DQt drawer(qimg.width(), qimg.height(), &qpt);
      drawer.drawMolecule(*mol, "wavy bonds");
      qimg.save("testWavyBonds.qt.1.png");
    }
  }
}

TEST_CASE("drawMoleculeBrackets") {
  SECTION("basics") {
    auto mol = R"CTAB(
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
    REQUIRE(mol);
    {
      QImage qimg(400, 400, QImage::Format_ARGB32);
      QPainter qpt(&qimg);
      MolDraw2DQt drawer(qimg.width(), qimg.height(), &qpt);
      drawer.drawMolecule(*mol, "brackets");
      qimg.save("testBrackets-qt-1a.png");
    }
  }
}
int main(int argc, char* argv[]) {
  QGuiApplication app(argc, argv);

  int result = Catch::Session().run(argc, argv);

  return result;
}