//
// Substructure searching - example14.cpp

#include <iostream>
#include <vector>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

int main( int argc , char **argv ) {

  RDKit::ROMOL_SPTR mol1( RDKit::SmilesToMol( "c1ccccc1O" ) );
  RDKit::RWMOL_SPTR patt( RDKit::SmartsToMol( "ccO" ) );
  RDKit::MatchVectType res;
  if( RDKit::SubstructMatch( *mol1 , *patt , res ) ) {
    std::cout << "Pattern matched molecule : " << std::endl;
    for( size_t i = 0 ; i < res.size() ; ++i ) {
      std::cout << "(" << res[i].first << "," << res[i].second << ")";
    }
    std::cout << std::endl;
  }


  std::vector<RDKit::MatchVectType> hits_vect;
  if( RDKit::SubstructMatch( *mol1 , *patt , hits_vect ) ) {
    for( size_t i = 0 ; i < hits_vect.size() ; ++i ) {
      std::cout << "Match " << i + 1 << " : ";
      for( size_t j = 0 ; j < hits_vect[i].size() ; ++j ) {
  	std::cout << "(" << hits_vect[i][j].first << ","
  		  << hits_vect[i][j].second << ")";
      }
      std::cout << std::endl;
    }
  }

  std::string file_root = getenv( "RDBASE" );
  file_root += "/Docs/Book";

  std::string sdf_file = file_root + "/data/actives_5ht3.sdf";
  RDKit::SDMolSupplier mol_supplier( sdf_file , true );
  RDKit::RWMOL_SPTR patt1( RDKit::SmartsToMol( "c[NH1]" ) );
  std::vector<RDKit::ROMOL_SPTR> matches;
  while( !mol_supplier.atEnd() ) {
    RDKit::ROMOL_SPTR mol3( mol_supplier.next() );
    if( mol3 && RDKit::SubstructMatch( *mol3 , *patt1 , res ) ) {
      matches.push_back( mol3 );
    }
  }
  std::cout << "There were " << matches.size() << " hits in the file." << std::endl;

  RDKit::ROMOL_SPTR mol4( RDKit::SmilesToMol( "C1=CC=CC=C1OC" ) );
  RDKit::RWMOL_SPTR smi_mol1( RDKit::SmilesToMol( "CO" ) );
  if( RDKit::SubstructMatch( *mol4 , *smi_mol1 , res ) ) {
    std::cout << "SMILES match" << std::endl;
  } else {
    std::cout << "Not SMILES match" << std::endl;
  }
  RDKit::RWMOL_SPTR smt_mol1( RDKit::SmartsToMol( "CO" ) );
  if( RDKit::SubstructMatch( *mol4 , *smt_mol1 , res ) ) {
    std::cout << "SMARTS match" << std::endl;
  } else {
    std::cout << "Not SMARTS match" << std::endl;
  }

  RDKit::RWMOL_SPTR smi_mol2( RDKit::SmilesToMol( "COC" ) );
  if( RDKit::SubstructMatch( *mol4 , *smi_mol2 , res ) ) {
    std::cout << "SMILES match" << std::endl;
  } else {
    std::cout << "Not SMILES match" << std::endl;
  }
  RDKit::RWMOL_SPTR smt_mol2( RDKit::SmartsToMol( "COC" ) );
  if( RDKit::SubstructMatch( *mol4 , *smt_mol2 , res ) ) {
    std::cout << "SMARTS match" << std::endl;
  } else {
    std::cout << "Not SMARTS match" << std::endl;
  }
  // Needs aromatic C
  RDKit::RWMOL_SPTR smt_mol3( RDKit::SmartsToMol( "COc" ) );
  if( RDKit::SubstructMatch( *mol4 , *smt_mol3 , res ) ) {
    std::cout << "SMARTS match" << std::endl;
  } else {
    std::cout << "Not SMARTS match" << std::endl;
  }
}
