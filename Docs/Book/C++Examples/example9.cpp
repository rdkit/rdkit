//
// Modifying molecules example9.cpp

#include <iostream>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>

int main( int argc , char **argv ) {

  std::shared_ptr<RDKit::RWMol> mol( new RDKit::RWMol( *RDKit::SmilesToMol( "c1ccccc1" ) ) );
  std::cout << "Order : " << mol->getBondWithIdx( 0 )->getBondType() << std::endl;
  std::cout << "Aromatic : " << mol->getBondWithIdx( 0 )->getIsAromatic() << std::endl;

  RDKit::MolOps::Kekulize( *mol );
  std::cout << "After default Kekulize : Order : " << mol->getBondWithIdx( 0 )->getBondType() << std::endl;
  std::cout << "After default Kekulize : Aromatic : " << mol->getBondWithIdx( 0 )->getIsAromatic() << std::endl;
 
  std::shared_ptr<RDKit::RWMol> mol1( new RDKit::RWMol( *RDKit::SmilesToMol( "c1ccccc1" ) ) );
  RDKit::MolOps::Kekulize( *mol1 , false );
  std::cout << "After Kekulize, markAtomsBonds false : Aromatic : " << mol1->getBondWithIdx( 0 )->getIsAromatic() << std::endl;

  std::shared_ptr<RDKit::RWMol> mol2( new RDKit::RWMol( *RDKit::SmilesToMol( "c1ccccc1" ) ) );
  RDKit::MolOps::Kekulize( *mol2 , true );
  std::cout << "After Kekulize, markAtomsBonds true : Aromatic : " << mol2->getBondWithIdx( 0 )->getIsAromatic() << std::endl;

  RDKit::MolOps::sanitizeMol( *mol );
  std::cout << "Order : " << mol->getBondWithIdx( 0 )->getBondType() << std::endl;
  std::cout << "Aromatic : " << mol->getBondWithIdx( 0 )->getIsAromatic() << std::endl;
  
}
