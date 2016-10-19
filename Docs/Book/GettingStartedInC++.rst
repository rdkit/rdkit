Getting Started with the RDKit in C++
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



What is this?
*************

This document is intended to provied an overview of how one can use
the RDKit functionality from C++. Like the 'Getting Started with the
RDKit in Python' it is not comprehensive and it's not a manual.  It is
modelled very closely on the Python version, and most of the text will
be similar if not identical.

Building and Running C++ RDKit Programs
***************************************

Unlike Python scripts, which, once the environment is setup up
correctly, can be run directly from a script or Jupyter notebook
session, the C++ programs must be compiled and linked before they can
be run. This creates an additional step, which varies with operating
system.  Using the program cmake makes it easier, but you may need to
experiment and change settings to get this working.  You will
need a reasonably modern C++ compiler. On linux systems this will most
likely be the GNU compiler gcc, although it could be Clang.  On Macs,
the opposite is true, and on Windows machines it will probably be
Visual C++.  At present, RDKit uses a relatively old-fashioned dialect
of C++, so gcc version 3.4 or higher is adequate, but this is due to
change within the next year.

The functions of the RDKit system are declared in a large number of
different header files spread across several directories in the
system, and defined across a number of different libraries.  The right
headers will need to be included in the source code, and libraries
linked to during linking.  Whilst it's possible to include all headers
and libraries in all executables this will result in bloated
executables which will take up a lot more disk and memory, especially
if you are doing static linking.  When linking to the static (.a)
libraries rather than the shared-object (.so) ones, the order the
libraries appear in linking list can be important.  See the
CMakeLists.txt file in C++Examples directory for a good order.  In
this case, the same library list is used for all examples, so some
will be unnecessary for some of the programs. The first 3 programs
don't need the Depictor and SubstructMatch libraries, for instance
although on my Ubuntu 16.04 system, the RDGeometryLib appears to need
to be included twice. Working out which libraries need to be linked to
and in what order can involve a tedious amount of trial and error!

Should You Use C++ or Python?
*****************************

There is no doubt that is is much easier to start learning with
Python.  If you follow the installation instructions, you will be able
to start programming and using scripts straightaway.  If all you are
going to do is use scripts to do relatively simple things,
essentially stitching RDKit function calls together, there should be
little or no speed issues with using the Python interpreted language,
as all the RDKit functions are compiled C++ and well optimised.
However, if you are going to do more complicated things, using a lot
of your own programming logic and only using the RDKit for peripheral
things like I/O, SMARTS matching, preparing 2D images and the like,
then it is likely that you will have good performance gains if you
write in C++ and compile to a native executable.

The Molecule Objects
********************

Unlike in the Python libraries, in C++ there are two different
molecule objects, `RDKit::ROMol` and `RDKit::RWMol`.  They are both
declared in GraphMol.h. ROMol (the Read-Only molecule) is used in
most instances. It can't be edited and therefore can be used as a
`const` parameter in function calls, which allows the compiler more
freedom to optimise code.  On those occasions where you will need to
edit the molecule, you'll need to use the RWMol (Read-Write).

Reading and Writing Molecules
*****************************

The majority of basic molecular functionality is found in the RDKit
namespace, and only a small number of header files will need to be included
to cover most use cases: ::

    #include <GraphMol/GraphMol.h>
    #include <GraphMol/FileParsers/MolSupplier.h>
    #include <GraphMol/FileParsers/MolWriters.h>

Individual molecules can be constructed using a variety of approaches: ::

    RDKit::ROMol *mol1 = RDKit::SmilesToMol( "Cc1ccccc1" );

    RDKit::MolSupplier mol_supplier( "data/input.mol" , true );
    RDKit::ROMol *mol2 = mol_supplier.next();

    RDKit::ROMol *mol3 = RDKit::SmilesToMol( "Cc1cccc" );

All these return a pointer to an ROMol on success, or 0 (aka Null) on
failure. Obviously the object must be deleted when finished with to
prevent memory leaks.

If the SMILES parsing fails, an `RDKit::MolSanitizeException` (derived
from `std::exception`) is thrown, and an attempt is made to provide
sensible error messages: :: 

   try {
      RDKit::ROMol *mol = RDKit::SmilesToMol( "CO(C)C" )
   } catch( std::exception &e ) {
      // empty catch
   }

displays something like ``[15:58:22] Explicit valence of atom # 1 O, 3,
is greater than permitted`` and ::
  
   try {
      RDKit::ROMol *mol = RDKit::SmilesToMol( "c1cc1" )
   } catch( std::exception &e ) {
      // empty catch
   }
 
displays something like: ``[12:20:41] Can't kekulize mol``.

Reading sets of molecules
=========================

Groups of molecule are read using a Supplier (for example, an
`RDKit::SDMolSupplier` or an `RDKit::SmilesMolSupplier` ::
  
   RDKit::SDMolSupplier mol_supplier( "data/5ht3ligs.sdf" , true );

   while( !mol_supplier.atEnd() ) {
	  mol = mol_supplier.next();
	  std::cout << mol->getNumAtoms() << std::endl;
	  delete mol;
   }
   20
   24
   24
   26

The supplier can be treated as a random-access object: ::

   RDKit::SDMolSupplier mol_supplier( "data/5ht3ligs.sdf" , true );

   for( int i = int( mol_supplier.length() ) - 1 ; i >= 0  ; --i ) {
	RDKit::ROMol *mol = mol_supplier[i];
	std::cout << mol->getProp<std::string>( "_Name" ) << " has "
	          << mol->getNumAtoms() << " atoms." << std::endl;
	delete mol;
   }
   
   mol-732 has 26 atoms.
   mol-15 has 24 atoms.
   mol-54 has 24 atoms.
   mol-295 has 20 atoms.

A good practice is to test each molecule to see if it was correctly
read before working with it: ::

   RDKit::SDMolSupplier *mol_supplier = new RDKit::SDMolSupplier( "data/5ht3ligs.sdf" , true );

   for( int i = int( mol_supplier->length() ) - 1 ; i >= 0 ; --i ) {
	RDKit::ROMol *mol = (*mol_supplier)[i];
	if( !mol ) {
		continue;
	}
	std::cout << mol->getProp<std::string>( "_Name" ) << " has "
	          << mol->getNumAtoms() << " atoms." << std::endl;
	delete mol;
   }

An alternative type of Supplier, the `RDKit::ForwardMolSupplier`
can be used to read from file-like objects.  This allows the reading
of compressed files, using, for example, the `boost::iostreams`
objects: ::

   boost::iostreams::filtering_istream ins;
   ins.push( boost::iostreams::gzip_decompressor() );
   ins.push( boost::iostreams::file_source( "data/actives_5ht3.sdf.gz" ) );

   RDKit::ForwardSDMolSupplier forward_supplier( &ins , true );
   while( !forward_supplier.atEnd() ) {
     mol = forward_supplier.next();
     std::cout << mol->getProp<std::string>( "_Name" ) << " has " << mol->getNumAtoms() << " atoms." << std::endl;
     delete mol;
   }

Note that the forward suppliers cannot be used in random-access mode,
and a compile-time error will result if you attempt to: ::

   error: no match for ‘operator[]’ (operand types are
   ‘RDKit::ForwardSDMolSupplier’ and ‘int’) 
      mol = forward_supplier[1];


Writing molecules
=================

Single molecules can be converted to text using several functions
present in the `RDKit` namespace.

For example, for SMILES: ::

   #include <GraphMol/SmilesParse/SmilesWrite.h>
   .
   .
   RDKit::ROMol *mol = RDKit::MolFromMolFile( "data/chiral.mol" );
   std::cout << RDKit::MolToSmiles( *mol ) << std::endl;

   CC(O)c1ccccc1
   
   std::cout << RDKit::MolToSmiles( *mol , true ) << std::endl;

   C[C@H](O)c1ccccc1

Note that the SMILES provided is canonical, so the output should be
the same no matter how a particular molecule is input.  For example ::

   RDKit::ROMol *mol1 = RDKit::SmilesToMol( "C1=CC=CN=C1" );
   std::cout << RDKit::MolToSmiles( *mol1 ) << std::endl;

   RDKit::ROMol *mol2 = RDKit::SmilesToMol( "c1cccnc1" );
   std::cout << RDKit::MolToSmiles( *mol2 ) << std::endl;

   RDKit::ROMol *mol3 = RDKit::SmilesToMol( "n1ccccc1" );
   std::cout << RDKit::MolToSmiles( *mol3 ) << std::endl;

all produce c1ccncc1 as output.

If you'd like to have the Kekule form of the SMILES, you need to Kekulise
a RWMol copy the molecule, using the Kekulize function declared in
MolOps.h: ::

   #include <GraphMol/MolOps.h>
   .
   .
   RDKit::RWMol *mol4 = new RDKit::RWMol( *mol );
   RDKit::MolOps::Kekulize( *mol4 );
   std::cout << RDKit::MolToSmiles( *mol4 ) << std::endl;

   CC(O)C1=CC=CC=C1
   
Note: as of this writing (Aug 2008), the smiles provided when one
requests kekuleSmiles are not canonical. The limitation is not in the
SMILES generation, but in the kekulization itself.

MDL Mol blocks are also available: ::
  
    RDKit::ROMol *mol1 = RDKit::SmilesToMol( "C1=CC=CN=C1" );
    std::cout << RDKit::MolToMolBlock( *mol1 ) << std::endl;

    <BLANKLINE>    
        RDKit          
    <BLANKLINE>    
    
      6  6  0  0  0  0  0  0  0  0999 V2000
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  2  0
      2  3  1  0
      3  4  2  0
      4  5  1  0
      5  6  2  0
      6  1  1  0
    M  END
    <BLANKLINE>
 
To include names in the mol blocks, set the molecule's “_Name”
property: ::
  
   mol1 = RDKit::SmilesToMol( "C1CCC1" );
   mol1->setProp( "_Name" , "cyclobutane" );
   std::cout << RDKit::MolToMolBlock( *mol1 ) << std::endl;

    cyclobutane
         RDKit          
    <BLANKLINE>
      4  4  0  0  0  0  0  0  0  0999 V2000
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0
      2  3  1  0
      3  4  1  0
      4  1  1  0
    M  END
    <BLANKLINE>
    
Note that setProp, which is a general function, can be called on an
ROMol as well as an RWMol which came as a surprise to me.

In order for atom or bond stereochemistry to be recognised correctly by most
software, it's essential that the Mol block have atomic coordinates.
It's also convenient for many reasons, such as drawing the molecules.

You can either include 2D coordinates (i.e. a depiction), using the
function in the RDDepict namespace and declared in RDDepictor.h ::

    #include <GraphMol/Depictor/RDDepictor.h>
    .
    .
    RDKit::ROMol *mol1 = RDKit::SmilesToMol( "C1CCC1" );
    RDDepict::compute2DCoords( *mol1 );
    std::cout << RDKit::MolToMolBlock( *mol1 ) << std::endl;


    <BLANKLINE>
         RDKit          2D
    <BLANKLINE>
      4  4  0  0  0  0  0  0  0  0999 V2000
        1.0607    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -0.0000   -1.0607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
       -1.0607    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    1.0607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0
      2  3  1  0
      3  4  1  0
      4  1  1  0
    M  END

Or you can add 3D coordinates by embedding the molecule: ::
  
    #include <GraphMol/DistGeomHelpers/Embedder.h>
    #include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
    .
    .
    RDKit::ROMol *mol2 = RDKit::SmilesToMol( "C1CCC1" );
    mol2->setProp( "_Name" , "cyclobutane3D" );
    RDKit::DGeomHelpers::EmbedMolecule( *mol2 );
    RDKit::MMFF::MMFFOptimizeMolecule( *mol2 , 1000 , "MMFF94s" );
    std::cout << RDKit::MolToMolBlock( *mol2 ) << std::endl;
    
    cyclobutane3D
         RDKit          3D
    <BLANKLINE>    
      4  4  0  0  0  0  0  0  0  0999 V2000
       -0.8321    0.5405   -0.1981 C   0  0  0  0  0  0  0  0  0  0  0  0
       -0.3456   -0.8799   -0.2639 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.7190   -0.5613    0.7314 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.4587    0.9006    0.5008 C   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0
      2  3  1  0
      3  4  1  0
      4  1  1  0
    M  END
    
The optimization step isn't necessary, but it substantially improves
the quality of the conformation.

To get good 3D conformations, it's almost always a good idea to add
hydrogens to the molecule first: ::

    RDKit::ROMol *mol3 = RDKit::MolOps::addHs( *mol2 );
    RDKit::MMFF::MMFFOptimizeMolecule( *mol3 , 1000 , "MMFF94s" );

    RDKit::RWMol *mol4 = new RDKit::RWMol( *mol3 );
    RDKit::MolOps::addHs( *mol4 );

Note that there are 2 overloaded versions of addHs. The first takes an
ROMol and, because that can't be edited, returns a pointer to a new
ROMol with the result.  If you use this version be careful not to leak
memory by not deleting mol2 when you are finished with it. The second
takes an RWMol which it is able to modify in place.

Once the optimisation is complete, the hydrogens can be removed
again: ::

    RDKit::ROMol *mol5 = RDKit::MolOps::removeHs( *mol3 );
    RDKit::MolOps::removeHs( *mol4 );

 Again, there are two versions, one of which has an opportunity for a
 memory leak.

 If you'd like write the molecules to file, use the normal C++
 streams: ::

    #include <fstream>
    .
    .
    std::ofstream ofs( "data/foo.mol" );
    ofs << RDKit::MolToMolBlock( *mol5 );

Writing sets of molecules
=========================

Multiple molecules can be written to a file using an object of a
concrete subclass of the `MolWriter` class: ::

    #include <GraphMol/FileParsers/MolWriters.h>
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

A MolWriter can also be initialised to a file-like object, so
compressed files can be written or it can be written to a string in
memory: ::

    #include <sstream>
    .
    .
    std::ostringstream oss;
    RDKit::SDWriter *sdf_writer = new RDKit::SDWriter( &oss , false );
    // Note that this requires a C++11 compliant compiler
    for( auto it = mols.begin() ; it != mols.end() ; ++it ) {
      sdf_writer->write( *(*it) );
    }

    std::cout << oss.str() << std::endl;

Other available writers include SmilesWriter and TDTWriter (for those
of you with an interest in historical Cheminformatics!)
