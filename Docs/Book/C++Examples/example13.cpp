//
// Generating depictions - example12.cpp

#include <fstream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>

int main( int argc , char **argv ) {

  RDKit::SDMolSupplier mol_supplier( "data/cdk2.sdf" , true );
  RDKit::ROMol *mol1 = mol_supplier.next();
  RDDepict::compute2DCoords( *mol1 );
  std::ofstream outs("cdk_mol1.svg");
  RDKit::MolDraw2DSVG svg_drawer(300, 300, outs);
  svg_drawer.drawMolecule( *mol1 );
  svg_drawer.finishDrawing();
  outs.close();

  RDKit::MolDraw2DCairo cairo_drawer(300, 300);
  cairo_drawer.drawMolecule(*mol1);
  cairo_drawer.finishDrawing();
  cairo_drawer.writeDrawingText("cdk_mol1.png");
  
}
