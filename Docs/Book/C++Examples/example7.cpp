//
// Working with molecules example7.cpp

#include <iostream>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>

int main( int argc , char **argv ) {

  RDKit::ROMOL_SPTR mol( RDKit::SmilesToMol( "OC1C2C1CC2" ) );

  if( !mol->getRingInfo()->isInitialized() ) {
    RDKit::MolOps::findSSSR( *mol );
  }
  for( unsigned int i = 0 , is = mol->getNumAtoms() ; i < is ; ++i ) {
    const RDKit::Atom *atom = mol->getAtomWithIdx( i );
    std::cout << mol->getRingInfo()->numAtomRings( atom->getIdx() ) << " ";
  }
  std::cout << std::endl;

  for( unsigned int i = 0 , is = mol->getNumBonds() ; i < is ; ++i ) {
    const RDKit::Bond *bond = mol->getBondWithIdx( i );
    std::cout << mol->getRingInfo()->numBondRings( bond->getIdx() ) << " ";
  }
  std::cout << std::endl;

  const RDKit::Bond *bond = mol->getBondWithIdx( 1 );
  if( mol->getRingInfo()->numBondRings( bond->getIdx() )) {
    std::cout <<  "Bond " << bond->getIdx() << " is in a ring" << std::endl;;
  }

  std::cout << "Atom 2 is in ring of size 3 : "
	    << mol->getRingInfo()->isAtomInRingOfSize( 2 , 3 ) << std::endl;
  std::cout << "Atom 2 is in ring of size 4 : "
	    << mol->getRingInfo()->isAtomInRingOfSize( 2 , 4 ) << std::endl;
  std::cout << "Atom 2 is in ring of size 5 : "
	    << mol->getRingInfo()->isAtomInRingOfSize( 2 , 5 ) << std::endl;
  std::cout << "Bond 1 is in ring of size 3 : "
	    << mol->getRingInfo()->isBondInRingOfSize( 1 , 3 ) << std::endl;
  std::cout << "Atom 1 is in ring of size 5 : "
	    << mol->getRingInfo()->isAtomInRingOfSize( 1 , 5 ) << std::endl;

  RDKit::VECT_INT_VECT rings;
  RDKit::MolOps::symmetrizeSSSR( *mol , rings );
  std::cout << "Number of symmetric SSSR rings : " << rings.size() << std::endl;
  for( RDKit::VECT_INT_VECT::iterator it1 = rings.begin() , it1_end = rings.end() ; it1 != it1_end ; ++it1 ) {
    for( RDKit::INT_VECT::iterator it2 = it1->begin() , it2_end = it1->end() ; it2 != it2_end ; ++it2 ) {
      std::cout << *it2 << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Number of SSSR rings : " << RDKit::MolOps::findSSSR( *mol ) << std::endl;
  
}
