//
// file test1.cc
// Greg Landrum
// 26 Dec 2014
//
//

#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>


#include "MolDraw2D.H"
#include "MolDraw2DSVG.H"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace RDKit;

void test1(){
  {
    std::string smiles="C[C@H](O)C1=C(O[C@H](F)Cl)C=C1ONNC[NH3+]";
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
}

#ifdef RDK_CAIRO_BUILD
#include <cairo.h>
#include "MolDraw2DCairo.H"
void test2(){
  {
    std::string smiles="C[C@H](O)C1=C(O[C@H](F)Cl)C=C1ONNC[NH3+]";
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
