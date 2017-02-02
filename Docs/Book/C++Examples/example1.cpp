//
// Reading molecules - example1.cpp

#include <iostream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>

int main( int argc , char **argv ) {

  RDKit::ROMol *mol1 = RDKit::SmilesToMol( "Cc1ccccc1" );
  std::cout << mol1 << std::endl;

  std::string file_root = getenv( "RDBASE" );
  file_root += "/Docs/Book";

  std::string mol_file = file_root + "/data/input.mol";
  RDKit::SDMolSupplier mol_supplier( mol_file , true );
  RDKit::ROMol *mol2 = mol_supplier.next();
  std::cout << mol2 << std::endl;

  try {
    RDKit::ROMol *mol3 = RDKit::SmilesToMol( "CO(C)C" );
    std::cout << mol3 << std::endl;
  } catch( std::exception &e ) {
    // std::cout << e.what() << std::endl;
  }
  try {
    RDKit::ROMol *mol4 = RDKit::SmilesToMol( "c1cc1" );
    std::cout << mol4 << std::endl;
  } catch( std::exception &e ) {
    // std::cout << e.what() << std::endl;
  }
}

