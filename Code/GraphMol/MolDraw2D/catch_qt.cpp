//
//  Copyright (C) 2019 Greg Landrum
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
#include <GraphMol/MolDraw2D/MolDraw2DQt.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolDraw2D/MolDraw2DDetails.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Depictor/RDDepictor.h>

#include <QImage>
#include <QPainter>
#include <QtCore/QCoreApplication>

using namespace RDKit;

TEST_CASE("basic generate PNGs", "[drawing][Qt]") {
  SECTION("basics") {
    auto m1 = "C1N[C@@H]2OCC12"_smiles;
    REQUIRE(m1);
    QCoreApplication app();
    // QPixmap qpm(250, 200);
    // QPainter qpt(&qpm);
    QImage qimg(250, 200, QImage::Format_RGB32);
    QPainter qpt(&qimg);
    MolDraw2DQt drawer(qimg.width(), qimg.height(), qpt);
    MolDraw2DUtils::prepareAndDrawMolecule(drawer, *m1);
    qimg.save("qttest-1a.png");
  }
}
