//
// Working with molecules example6.cpp

#include <iostream>

#include <GraphMol/GraphMol.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

int main( int argc , char **argv ) {

  std::shared_ptr<RDKit::ROMol> mol( RDKit::SmilesToMol( "C1OC1" ) );

  for(auto atom: mol->atoms()) {
    std::cout << atom->getAtomicNum() << " ";
  }
  std::cout << std::endl;   

  for( unsigned int i = 0; i < mol->getNumAtoms() ; ++i ) {
    const RDKit::Atom *atom = mol->getAtomWithIdx( i ); 
    std::cout << atom->getAtomicNum() << " ";
  }
  std::cout << std::endl;   

  for(auto bond: mol->bonds()) {
    std::cout << bond->getBondType() << " ";
  }
  std::cout << std::endl;   

  for( unsigned int i = 0 , is = mol->getNumBonds() ; i < is ; ++i ) {
    const RDKit::Bond *bond = mol->getBondWithIdx( i ); 
    std::cout << bond->getIsAromatic() << " ";   
  }
  std::cout << std::endl;   

  std::shared_ptr<RDKit::ROMol> mol2( RDKit::SmilesToMol( "C1OC1Cl" ) );
  const RDKit::Bond *bond = mol2->getBondBetweenAtoms( 0 , 1 );
  std::cout << bond->getBeginAtomIdx() << " to "
	    << bond->getBeginAtomIdx() << " is "
	    << bond->getBondType() << std::endl;
  if( !mol2->getBondBetweenAtoms( 0 , 3 ) ) {
    std::cout << "No bond between 0 and 3" << std::endl;
  }

  const RDKit::Atom *atom = mol2->getAtomWithIdx( 2 );
  for(const auto &nbri: make_iterator_range(mol->getAtomBonds(atom))) {
    const RDKit::Bond *bond = (*mol)[nbri];
    unsigned int nbr_idx = bond->getOtherAtomIdx(atom->getIdx());
    int nbr_atnum = bond->getOtherAtom(atom)->getAtomicNum();
    std::cout << nbr_idx << " : " << nbr_atnum << std::endl;
  }
  
}

