//
// Substructure searching with stereochemistry - example15.cpp

#include <iostream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

int main( int argc , char **argv ) {

  std::shared_ptr<RDKit::ROMol> mol1( RDKit::SmilesToMol( "CC[C@H](F)Cl" ) );
  std::shared_ptr<RDKit::RWMol> patt1( RDKit::SmartsToMol( "C[C@H](F)Cl" ) );
  RDKit::MatchVectType res;
  if( RDKit::SubstructMatch( *mol1 , *patt1 , res ) ) {
    std::cout << "SMARTS 1 match" << std::endl;
  } else {
    std::cout << "Not SMARTS 1 match" << std::endl;
  }
  std::shared_ptr<RDKit::RWMol> patt2( RDKit::SmartsToMol( "C[C@@H](F)Cl" ) );
  if( RDKit::SubstructMatch( *mol1 , *patt2 , res ) ) {
    std::cout << "SMARTS 2 match" << std::endl;
  } else {
    std::cout << "Not SMARTS 2 match" << std::endl;
  }
  std::shared_ptr<RDKit::RWMol> patt3( RDKit::SmartsToMol( "CC(F)Cl" ) );
  if( RDKit::SubstructMatch( *mol1 , *patt3 , res ) ) {
    std::cout << "SMARTS 3 match" << std::endl;
  } else {
    std::cout << "Not SMARTS 3 match" << std::endl;
  }

  if( RDKit::SubstructMatch( *mol1 , *patt1 , res , true , true ) ) {
    std::cout << "SMARTS 1 chiral match" << std::endl;
  } else {
    std::cout << "Not SMARTS 1 chiral match" << std::endl;
  }
  if( RDKit::SubstructMatch( *mol1 , *patt2 , res , true , true ) ) {
    std::cout << "SMARTS 2 chiral match" << std::endl;
  } else {
    std::cout << "Not SMARTS 2 chiral match" << std::endl;
  }
  if( RDKit::SubstructMatch( *mol1 , *patt3 , res , true , true ) ) {
    std::cout << "SMARTS 3 chiral match" << std::endl;
  } else {
    std::cout << "Not SMARTS 3 chiral match" << std::endl;
  }

  std::shared_ptr<RDKit::RWMol> mol2( RDKit::SmilesToMol( "CC(F)Cl" ) );
  if( RDKit::SubstructMatch( *mol1 , *mol2 , res , true , true ) ) {
    std::cout << "Chiral mol, non-chiral query : match" << std::endl;
  } else {
    std::cout << "Chiral mol, non-chiral query : NO match" << std::endl;
  }
  
  std::shared_ptr<RDKit::RWMol> patt5( RDKit::SmilesToMol( "C[C@H](F)Cl" ) );
  if( RDKit::SubstructMatch( *mol2 , *patt5 , res , true , true ) ) {
    std::cout << "Non-chiral mol, chiral query : match" << std::endl;
  } else {
    std::cout << "Non-chiral mol, chiral query : NO match" << std::endl;
  }

}

