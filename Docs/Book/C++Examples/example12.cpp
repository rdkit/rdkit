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

  std::shared_ptr<RDKit::ROMol> mol1( RDKit::SmilesToMol( "c1ccncc1" ) );
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
    std::shared_ptr<RDKit::ROMol> mol( suppl.next() );
    // write all props (including _Name) to the pickle.  Default is not to:
    // RDKit::MolPickler::pickleMol( *mol , pickle_ostream);
    // RDKit::MolPickler::pickleMol( *mol , pickle_ostream, RDKit::PicklerOps::NoProps);
    RDKit::MolPickler::pickleMol( *mol , pickle_ostream, RDKit::PicklerOps::AllProps );
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
    if(!read_cnt) {
      try {
	std::cout << "Read mol of name : "
		  << mol3.getProp<std::string>("_Name") << std::endl;
      } catch(KeyErrorException &e) {
	std::cerr << "Pickle read error : " << e.what() << std::endl;
      }
    }
    ++read_cnt;
  }
  pickle_istream.close();
  std::cout << "Read " << read_cnt << " molecules." << std::endl;
  
}


