//
//  Copyright (C) 2015 Greg Landrum 
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>

#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace RDKit;

void test1(){
  {
    std::string smiles="CO[C@@H](O)C1=C(O[C@H](F)Cl)C=C1ONNC[NH3+]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m,&(m->getConformer()));
    std::ofstream outs("test1_1.svg");
    MolDraw2DSVG drawer(300,300,outs);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.flush();
  }
  {
    std::string smiles="Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m,&(m->getConformer()));
    std::ofstream outs("test1_2.svg");
    MolDraw2DSVG drawer(300,300,outs);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.flush();
  }
}

#ifdef RDK_CAIRO_BUILD
#include <cairo.h>
#include "MolDraw2DCairo.h"
void test2(){
  {
    std::string smiles="CO[C@@H](O)C1=C(O[C@H](F)Cl)C=C1ONNC[NH3+]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m,&(m->getConformer()));

    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);

    MolDraw2DCairo drawer(300,300,cr);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "test2_1.png");
    cairo_surface_destroy (surface);
  }
  {
    std::string smiles="Cc1c(C(=O)NCCO)[n+](=O)c2ccccc2n1[O-]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    WedgeMolBonds(*m,&(m->getConformer()));

    cairo_surface_t *surface =
      cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 300, 300);
    cairo_t *cr = cairo_create (surface);

    MolDraw2DCairo drawer(300,300,cr);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();

    cairo_destroy (cr);
    cairo_surface_write_to_png (surface, "test2_2.png");
    cairo_surface_destroy (surface);
  }
}
#else // RDK_CAIRO_BUILD
void test2(){
}
#endif

int main(){
  RDLog::InitLogs();
  test1();
  test2();
}
