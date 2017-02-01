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

  RDKit::ROMol *mol;
  RDKit::SDMolSupplier mol_supplier( "data/5ht3ligs.sdf" , true );
  std::vector<RDKit::ROMol *> mols;
  while( !mol_supplier.atEnd() ) {
    mol = mol_supplier.next();
    if( mol ) {
      mols.push_back( mol );
    }
  }

  RDKit::PDBWriter pdb_writer( "data/5ht3ligs.pdb" );
  for( std::size_t i = 0 , is = mols.size() ; i < is ; ++i ) {
    pdb_writer.write( *mols[i] );
  }

  std::ostringstream oss;
  RDKit::SDWriter *sdf_writer = new RDKit::SDWriter( &oss , false );
  // Note that this requires a C++11 compliant compiler
  for( std::vector<RDKit::ROMol *>::iterator it = mols.begin() ; it != mols.end() ; ++it ) {
    sdf_writer->write( *(*it) );
  }

  std::cout << oss.str() << std::endl;
  
}
