//
// Reading molecules - example1.cpp

#include <iostream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

int main( int argc , char **argv ) {

  RDKit::ROMol *mol1 = RDKit::SmilesToMol( "Cc1ccccc1" );

  std::string file_root = getenv( "RDBASE" );
  file_root += "/Docs/Book";

  std::string mol_file = file_root + "/data/input.mol";
  RDKit::ROMOL_SPTR mol2( RDKit::MolFileToMol( mol_file ) );

  try {
    RDKit::ROMOL_SPTR mol3( RDKit::SmilesToMol( "CO(C)C" ) );
  } catch( RDKit::MolSanitizeException &e ) {
    // std::cout << e.what() << std::endl;
  }
  try {
    RDKit::ROMOL_SPTR mol4( RDKit::SmilesToMol( "c1cc1" ) );
  } catch( RDKit::MolSanitizeException &e ) {
    // std::cout << e.what() << std::endl;
  }
}

