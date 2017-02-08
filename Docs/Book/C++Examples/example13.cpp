//
// Generating depictions - example13.cpp

#include <fstream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>

int main( int argc , char **argv ) {

  std::string file_root = getenv( "RDBASE" );
  file_root += "/Docs/Book";

  std::string sdf_file = file_root + "/data/cdk2.sdf";
  RDKit::SDMolSupplier mol_supplier( sdf_file , true );
  RDKit::ROMol *mol1 = mol_supplier.next();
  RDDepict::compute2DCoords( *mol1 );
  std::string svg_file = file_root + "/data/cdk_mol1.svg";
  std::ofstream outs( svg_file.c_str() );
  RDKit::MolDraw2DSVG svg_drawer(300, 300, outs);
  svg_drawer.drawMolecule( *mol1 );
  svg_drawer.finishDrawing();
  outs.close();

  RDKit::MolDraw2DCairo cairo_drawer(300, 300);
  cairo_drawer.drawMolecule(*mol1);
  cairo_drawer.finishDrawing();
  std::string png_file = file_root + "/data/cdk_mol1.png";
  cairo_drawer.writeDrawingText( png_file );
  
}
