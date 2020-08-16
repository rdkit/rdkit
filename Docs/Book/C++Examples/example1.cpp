//
// Reading molecules - example1.cpp

#include <iostream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

int main( int argc , char **argv ) {

  RDKit::ROMol *mol1 = RDKit::SmilesToMol( "Cc1ccccc1" );
  std::cout << "Number of atoms " << mol1->getNumAtoms() << std::endl;

  using namespace RDKit;
  auto mol = "C[C@H](F)c1ccc(C#N)cc1"_smiles;
  std::cout << "Number of atoms : " << mol->getNumAtoms() << std::endl;
  
  std::string file_root = getenv( "RDBASE" );
  file_root += "/Docs/Book";

  std::string mol_file = file_root + "/data/input.mol";
  std::shared_ptr<ROMol> mol2( RDKit::MolFileToMol( mol_file ) );

  try {
    std::shared_ptr<ROMol> mol3( RDKit::SmilesToMol( "CO(C)C" ) );
  } catch( RDKit::MolSanitizeException &e ) {
    std::cout << e.what() << std::endl;
  }
  try {
    std::unique_ptr<ROMol> mol4( RDKit::SmilesToMol( "c1cc1" ) );
  } catch( RDKit::MolSanitizeException &e ) {
    std::cout << e.what() << std::endl;
  }
}

