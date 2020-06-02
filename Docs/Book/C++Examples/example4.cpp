//
// Writing molecules - example4.cpp

#include <fstream>
#include <iostream>
#include <string>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>

int main( int argc , char **argv ) {

  std::shared_ptr<RDKit::ROMol> mol1( RDKit::SmilesToMol( "C1CCC1" ) );
  std::cout << RDKit::MolToMolBlock( *mol1 ) << std::endl;

  mol1->setProp( "_Name" , "cyclobutane" );
  std::cout << RDKit::MolToMolBlock( *mol1 ) << std::endl;

  RDDepict::compute2DCoords( *mol1 );
  std::cout << RDKit::MolToMolBlock( *mol1 ) << std::endl;

  std::shared_ptr<RDKit::ROMol> mol2( RDKit::SmilesToMol( "C1CCC1" ) );
  mol2->setProp( "_Name" , "cyclobutane3D" );
  RDKit::DGeomHelpers::EmbedMolecule( *mol2 );
  RDKit::MMFF::MMFFOptimizeMolecule( *mol2 , 1000 , "MMFF94s" );
  std::cout << RDKit::MolToMolBlock( *mol2 ) << std::endl;

  std::shared_ptr<RDKit::ROMol> mol3( RDKit::MolOps::addHs( *mol2 ) );
  RDKit::MMFF::MMFFOptimizeMolecule( *mol3 , 1000 , "MMFF94s" );
  std::cout << RDKit::MolToMolBlock( *mol3 ) << std::endl;

  std::shared_ptr<RDKit::RWMol> mol4( new RDKit::RWMol( *mol3 ) );
  RDKit::MolOps::addHs( *mol4 );

  std::shared_ptr<RDKit::ROMol> mol3sp( RDKit::MolOps::addHs( *mol2 ) );
  mol3sp->setProp( "_Name" , "cyclobutaneSP" );
  RDKit::MMFF::MMFFOptimizeMolecule( *mol3sp , 1000 , "MMFF94s" );
  std::cout << RDKit::MolToMolBlock( *mol3sp ) << std::endl;


  std::shared_ptr<RDKit::ROMol> mol5( RDKit::MolOps::removeHs( *mol3 ) );
  RDKit::MolOps::removeHs( *mol4 );

  std::string file_root = getenv( "RDBASE" );
  file_root += "/Docs/Book";

  std::string mol_file = file_root + "/data/foo.mol";
  std::ofstream ofs( mol_file.c_str() );
  ofs << RDKit::MolToMolBlock( *mol5 );

}
