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
#include <GraphMol/MolDraw2D/MolDraw2DQt.h>
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
    MolDraw2DQt drawer(qimg.width(), qimg.height(), qpt);
    MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
    qimg.save("qttest-1a.png");
  }
  SECTION("no freetype") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);
    QImage qimg(250, 200, QImage::Format_RGB32);
    QPainter qpt(&qimg);
    bool no_freetype = true;
    MolDraw2DQt drawer(qimg.width(), qimg.height(), qpt, -1, -1, no_freetype);
    MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
    qimg.save("qttest-1b.png");
  }
}

int main(int argc, char* argv[]) {
  QGuiApplication app(argc, argv);

  int result = Catch::Session().run(argc, argv);

  return result;
}