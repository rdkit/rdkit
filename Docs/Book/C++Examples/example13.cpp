//
// Generating depictions - example13.cpp

#include <fstream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#ifdef RDK_BUILD_CAIRO_SUPPORT
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#endif

using namespace RDKit;

int main( int argc , char **argv ) {

  std::string file_root = getenv( "RDBASE" );
  file_root += "/Docs/Book";

  std::string sdf_file = file_root + "/data/cdk2.sdf";
  std::shared_ptr<ROMol> mol1( MolFileToMol( sdf_file ) );
  RDDepict::compute2DCoords( *mol1 );
  std::string svg_file = file_root + "/data/cdk_mol1.svg";
  std::ofstream outs( svg_file.c_str() );
  MolDraw2DSVG svg_drawer(300, 300, outs);
  svg_drawer.drawMolecule( *mol1 );
  svg_drawer.finishDrawing();
  outs.close();

#ifdef RDK_BUILD_CAIRO_SUPPORT
  MolDraw2DCairo cairo_drawer(300, 300);
  cairo_drawer.drawMolecule(*mol1);
  cairo_drawer.finishDrawing();
  std::string png_file = file_root + "/data/cdk_mol1.png";
  cairo_drawer.writeDrawingText( png_file );
 #endif

  // drawing multiple molecules to a grid
  {
    std::string base_smi("c1ccccc1");
    std::vector<std::string> extras = {"F", "Cl", "Br", "OC", "C(=O)O"};
    std::vector<ROMol *> mols;
    for(auto extra: extras) {
      mols.push_back(SmilesToMol(base_smi + extra));
    }
    MolDraw2DSVG drawer(750, 400, 250, 200);
    drawer.drawMolecules(mols);
    drawer.finishDrawing();
    std::ofstream grids(file_root + "/data/example_13_grid.svg");
    grids << drawer.getDrawingText();
    grids.flush();
    grids.close();
    
    for(auto mol: mols) {
      delete mol;
    }
  }

  // adding annotation to the molecule
  {
    auto m1 = "Cl[C@H](F)NC\\C=C\\C"_smiles;
    MolDraw2DSVG drawer(250, 200);
    m1->getAtomWithIdx(2)->setProp(common_properties::atomNote, "foo");
    m1->getBondWithIdx(0)->setProp(common_properties::bondNote, "bar");
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawOptions().addStereoAnnotation = true;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs(file_root + "/images/example_13_note.svg");
    outs << text;
    outs.flush();
  }

}
