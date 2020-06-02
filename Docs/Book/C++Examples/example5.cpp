//
// Writing molecules example5.cpp

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>

int main( int argc , char **argv ) {

  std::string file_root = getenv( "RDBASE" );
  file_root += "/Docs/Book";

  std::string sdf_file = file_root + "/data/5ht3ligs.sdf";
  bool takeOwnership = true;
  RDKit::SDMolSupplier mol_supplier( sdf_file , takeOwnership );
  std::vector<std::shared_ptr<RDKit::ROMol>> mols;
  while( !mol_supplier.atEnd() ) {
    std::shared_ptr<RDKit::ROMol> mol( mol_supplier.next() );
    if( mol ) {
      mols.push_back( mol );
    }
  }

  std::string pdb_file = file_root + "/data/5ht3ligs.pdb";
  RDKit::PDBWriter pdb_writer( pdb_file );
  for( std::size_t i = 0 , is = mols.size() ; i < is ; ++i ) {
    pdb_writer.write( *mols[i] );
  }

  std::ostringstream oss;
  // takeOwnership must be false for this, as we don't want the SDWriter trying
  // to delete the std::ostringstream.
  takeOwnership = false;
  boost::shared_ptr<RDKit::SDWriter> sdf_writer( new RDKit::SDWriter( &oss , takeOwnership ) );
  for( auto mol: mols ) {
    sdf_writer->write( *mol );
  }

  std::cout << oss.str() << std::endl;
  
}
