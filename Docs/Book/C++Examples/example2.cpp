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

  std::string file_root = getenv( "RDBASE" );
  file_root += "/Docs/Book";

  std::unique_ptr<RDKit::ROMol> mol;
  std::string sdf_file = file_root + "/data/5ht3ligs.sdf";
  bool takeOwnership = true;
  RDKit::SDMolSupplier mol_supplier( sdf_file , takeOwnership );
  while( !mol_supplier.atEnd() ) {
    mol.reset( mol_supplier.next() );
    std::cout << mol->getProp<std::string>( "_Name" ) << " has " << mol->getNumAtoms() << " atoms." << std::endl;
  }

  for( int i = int( mol_supplier.length() ) - 1 ; i >= 0  ; --i ) {
    std::unique_ptr<RDKit::ROMol> mol( mol_supplier[i] );
    if( mol ) {
      std::cout << mol->getProp<std::string>( "_Name" ) << " has " << mol->getNumAtoms() << " atoms." << std::endl;
    }
  }

  boost::iostreams::filtering_istream ins;
  ins.push( boost::iostreams::gzip_decompressor() );
  std::string comp_sdf_file = file_root + "/data/actives_5ht3.sdf.gz";
  ins.push( boost::iostreams::file_source( comp_sdf_file ) );
  // takeOwnership must be false for this, as we don't want the SDWriter trying
  // to delete the boost::iostream
  takeOwnership = false;
  RDKit::ForwardSDMolSupplier forward_supplier( &ins , takeOwnership );
  while( !forward_supplier.atEnd() ) {
    mol.reset( forward_supplier.next() );
    if( mol ) {
      std::cout << mol->getProp<std::string>( "_Name" ) << " has " << mol->getNumAtoms() << " atoms." << std::endl;
    }
  }

  // This is not allowed, and will give a compiler error:
  // mol.reset(forward_supplier[1]);

}
