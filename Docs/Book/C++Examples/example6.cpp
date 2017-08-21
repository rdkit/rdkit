//
// Working with molecules example6.cpp

#include <iostream>

#include <GraphMol/GraphMol.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

int main( int argc , char **argv ) {

  RDKit::ROMOL_SPTR mol( RDKit::SmilesToMol( "C1OC1" ) );

  RDKit::ROMol::VERTEX_ITER it , end;
  boost::tie( it , end ) = mol->getVertices();
  while( it != end ) {
    const RDKit::Atom *atom = (*mol)[*it].get();
    std::cout << atom->getAtomicNum() << " ";
    ++it;
  }
  std::cout << std::endl;   

  for( RDKit::ROMol::AtomIterator ai = mol->beginAtoms() ; ai != mol->endAtoms() ; ++ai) {
    std::cout << (*ai)->getAtomicNum() << " ";
  }
  std::cout << std::endl;   
    
  for( unsigned int i = 0 , is = mol->getNumAtoms() ; i < is ; ++i ) {
    const RDKit::Atom *atom = mol->getAtomWithIdx( i ); 
    std::cout << atom->getAtomicNum() << " ";
  }
  std::cout << std::endl;   

  RDKit::ROMol::EDGE_ITER bond_it , bond_end;
  boost::tie( bond_it , bond_end ) = mol->getEdges();
  while( bond_it != bond_end ) {
    const RDKit::Bond *bond = (*mol)[*bond_it].get();
    std::cout << bond->getBondType() << " ";
    ++bond_it;
  }
  std::cout << std::endl;   

  for( unsigned int i = 0 , is = mol->getNumBonds() ; i < is ; ++i ) {
    const RDKit::Bond *bond = mol->getBondWithIdx( i ); 
    std::cout << bond->getIsAromatic() << " ";   
  }
  std::cout << std::endl;   

  RDKit::ROMOL_SPTR mol2( RDKit::SmilesToMol( "C1OC1Cl" ) );
  const RDKit::Bond *bond = mol2->getBondBetweenAtoms( 0 , 1 );
  std::cout << bond->getBeginAtomIdx() << " to "
	    << bond->getBeginAtomIdx() << " is "
	    << bond->getBondType() << std::endl;
  if( !mol2->getBondBetweenAtoms( 0 , 3 ) ) {
    std::cout << "No bond between 0 and 3" << std::endl;
  }

  const RDKit::Atom *atom = mol2->getAtomWithIdx( 2 );
  RDKit::ROMol::ADJ_ITER nbr , end_nbr;
  boost::tie( nbr , end_nbr ) = mol2->getAtomNeighbors( atom );
  while( nbr != end_nbr ) {
    const RDKit::Atom *nbr_atom = (*mol2)[*nbr].get();
    std::cout << nbr_atom->getIdx() << " : " << nbr_atom->getAtomicNum() << std::endl;
    ++nbr;
  }
  
}

