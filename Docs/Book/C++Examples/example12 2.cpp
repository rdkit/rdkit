//
// Preserving molecules - example12.cpp

#include <fstream>
#include <iostream>
#include <string>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

int main( int argc , char **argv ) {

  RDKit::ROMOL_SPTR mol1( RDKit::SmilesToMol( "c1ccncc1" ) );
  std::string pickle;
  RDKit::MolPickler::pickleMol( *mol1 , pickle );
  RDKit::ROMol mol2;
  RDKit::MolPickler::molFromPickle( pickle , mol2 );
  std::cout << RDKit::MolToSmiles( mol2 ) << std::endl;

  // writing to pickle file
  std::string smi_file = getenv("RDBASE");
  smi_file += "/Code/GraphMol/test_data/canonSmiles.long.smi";
  std::string pkl_name = "canonSmiles.long.bin";
  
  // tab-delimited file, SMILES in column 0, name in 1, no title line
  RDKit::SmilesMolSupplier suppl( smi_file , "\t" , 0 , 1 , false );
  std::ofstream pickle_ostream( pkl_name.c_str() , std::ios_base::binary );
  int write_cnt = 0;
  while( !suppl.atEnd() ) {
    RDKit::ROMOL_SPTR mol( suppl.next() );
    RDKit::MolPickler::pickleMol( *mol , pickle_ostream );
    ++write_cnt;
  }
  pickle_ostream.close();  
  std::cout << "Wrote " << write_cnt << " molecules" << std::endl;
  
  // reading from pickle file
  std::ifstream pickle_istream( pkl_name.c_str() , std::ios_base::binary );
  int read_cnt = 0;
  while( !pickle_istream.eof() ) {
    RDKit::ROMol mol3;
    try {
      RDKit::MolPickler::molFromPickle( pickle_istream , mol3 );
    } catch( RDKit::MolPicklerException &e ) {
      break;
    }
    ++read_cnt;
  }
  pickle_istream.close();
  std::cout << "Read " << read_cnt << " molecules." << std::endl;
  
}


