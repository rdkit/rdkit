//
// Reading molecules - example2.cpp

#include <fstream>
#include <iostream>
#include <string>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>

int main( int argc , char **argv ) {

  RDKit::ROMol *mol;
  RDKit::SDMolSupplier mol_supplier( "data/5ht3ligs.sdf" , true );
  while( !mol_supplier.atEnd() ) {
    mol = mol_supplier.next();
    std::cout << mol->getProp<std::string>( "_Name" ) << " has " << mol->getNumAtoms() << " atoms." << std::endl;
    delete mol;
  }

  for( int i = int( mol_supplier.length() ) - 1 ; i >= 0  ; --i ) {
    RDKit::ROMol *mol = mol_supplier[i];
    if( mol ) {
      std::cout << mol->getProp<std::string>( "_Name" ) << " has " << mol->getNumAtoms() << " atoms." << std::endl;
      delete mol;
    }
  }

  boost::iostreams::filtering_istream ins;
  ins.push( boost::iostreams::gzip_decompressor() );
  ins.push( boost::iostreams::file_source( "data/actives_5ht3.sdf.gz" ) );

  RDKit::ForwardSDMolSupplier forward_supplier( &ins , true );
  while( !forward_supplier.atEnd() ) {
    mol = forward_supplier.next();
    std::cout << mol->getProp<std::string>( "_Name" ) << " has " << mol->getNumAtoms() << " atoms." << std::endl;
    delete mol;
  }

  // This is not allowed, and will give a compiler error:
  // mol = forward_supplier[1];

}

