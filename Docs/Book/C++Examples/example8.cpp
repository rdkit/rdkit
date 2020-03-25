//
// Modifying molecules example8.cpp

#include <iostream>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>

int main( int argc , char **argv ) {

  std::shared_ptr<RDKit::ROMol> mol1( RDKit::SmilesToMol( "CCO" ) );
  std::cout << "Number of atoms : " << mol1->getNumAtoms() << std::endl;
  std::shared_ptr<RDKit::ROMol> mol2( RDKit::MolOps::addHs( *mol1 ) );
  std::cout << "Number of atoms : " << mol2->getNumAtoms() << std::endl;

  std::shared_ptr<RDKit::RWMol> mol3( new RDKit::RWMol( *mol2 ) );
  RDKit::MolOps::removeHs( *mol3 );
  std::cout << "Number of atoms : " << mol3->getNumAtoms() << std::endl;
  
}
