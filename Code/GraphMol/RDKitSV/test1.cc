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


#include "MolDraw2D.H"
#include "MolDraw2DSVG.H"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace RDKit;

void test1(){
  {
    std::string smiles="C1CCC1ONNC[NH3+]";
    ROMol *m = SmilesToMol(smiles);
    TEST_ASSERT(m);
    RDDepict::compute2DCoords(*m);
    std::ofstream outs("test1_1.svg");
    MolDraw2DSVG drawer(300,300,outs);
    drawer.drawMolecule(*m);
    drawer.finishDrawing();
    outs.flush();
  }
  
}

int main(){
  RDLog::InitLogs();
  test1();
}
