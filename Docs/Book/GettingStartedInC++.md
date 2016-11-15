# Getting Started with the RDKit in C++

## What is this?

This document is intended to provied an overview of how one can use
the RDKit functionality from C++. Like the 'Getting Started with the
RDKit in Python' it is not comprehensive and it's not a manual.  It is
modelled very closely on the Python version, and most of the text will
be similar if not identical.

## Building and Running C++ RDKit Programs

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

C++ programs generally have quite a lot of extra verbiage in them
compared to, for example, python scripts, that is similar or identical
in all programs.  The would produce extraneous clutter in the example
code in this document.  Because of this, only the minimum code to
exemplify the point being made is shown in the text. Full programs for
all the examples (generally with multiple examples in the same
program) are given in the `$RDBASE/Docs/Book/C++Examples` directory in
the distribution. The particular program will be mentioned in the text.

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
don't need the Depictor and SubstructMatch libraries, for instance,
although on my Ubuntu 16.04 system, the RDGeometryLib appears to need
to be included twice. Working out which libraries need to be linked to
and in what order can involve a tedious amount of trial and error!

## Should You Use C++ or Python?

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

## Memory Management

Memory leaks (where objects are created using `new` but never
destroyed using `delete`) are particularly insidious bugs that can be
difficult to track down
([Valgrind is your friend!](http://valgrind.org)) and may be a
surprise to people used to managed-memory languages like python and
java.  In RDKit, many of the functions return pointers to molecules
and accept such pointers as arguments.  The
[Boost libraries](http://boost.org) provide tools such as shared
(`boost::shared_ptr`) and scoped (`boost::scoped_ptr`) pointers that
can be very helpful in preventing memory leaks.  If an object is
created and a pointer to it is stored in a shared pointer, for
example, then when the shared pointer goes out of scope, the object
is automatically deleted.  Otherwise, the shared pointer is used
exactly as one would use an ordinary pointer.  To save some typing,
RDKit has a number of typedefs for shared pointers to its objects,
defined in the relevant header files.  Two particularly useful ones
are `RDKit::ROMOL_SPTR` and `RDKit::RWMOL_SPTR`, for `RDKit::ROMol`
and `RDKit::RWMol` objects respectively.

## The Molecule Objects

Unlike in the Python libraries, in C++ there are two different
molecule objects, `RDKit::ROMol` and `RDKit::RWMol`.  They are both
declared in GraphMol.h. ROMol (the Read-Only molecule) is used in
most instances. It can't be edited and therefore can be used as a
`const` parameter in function calls, which allows the compiler more
freedom to optimise code.  On those occasions where you will need to
edit the molecule, you'll need to use the RWMol (Read-Write).

## Reading and Writing Molecules

The majority of basic molecular functionality is found in the RDKit
namespace, and only a small number of header files will need to be included
to cover most use cases:

```c++
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
```

### Reading Single Molecules
Individual molecules can be constructed using a variety of approaches [(example1)](./C++Examples/example1.cpp):

```c++
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "Cc1ccccc1" );

RDKit::MolSupplier mol_supplier( "data/input.mol" , true );
RDKit::ROMol *mol2 = mol_supplier.next();

RDKit::ROMol *mol3 = RDKit::SmilesToMol( "Cc1cccc" );
```

All these return a pointer to an ROMol on success, or 0 (aka Null) on
failure. Obviously the object must be deleted when finished with to
prevent memory leaks.

If the SMILES parsing fails, an `RDKit::MolSanitizeException` (derived
from `std::exception`) is thrown, and an attempt is made to provide
sensible error messages [(example1)](./C++Examples/example1.cpp):

```c++
try {
   RDKit::ROMol *mol = RDKit::SmilesToMol( "CO(C)C" )
} catch( std::exception &e ) {
   // empty catch
}
```

displays something like `[15:58:22] Explicit valence of atom # 1 O, 3,
is greater than permitted` and [(example1)](./C++Examples/example1.cpp)

```c++
try {
   RDKit::ROMol *mol = RDKit::SmilesToMol( "c1cc1" )
} catch( std::exception &e ) {
   // empty catch
}
```

displays something like: `[12:20:41] Can't kekulize mol`.

### Reading sets of molecules

Groups of molecule are read using a Supplier (for example, an
`RDKit::SDMolSupplier` or an `RDKit::SmilesMolSupplier`)
[(example2)](./C++Examples/example2.cpp):

```c++
RDKit::SDMolSupplier mol_supplier( "data/5ht3ligs.sdf" , true );

while( !mol_supplier.atEnd() ) {
	  mol = mol_supplier.next();
	  std::cout << mol->getNumAtoms() << std::endl;
	  delete mol;
}
```
gives
```
20
24
24
26
```

The supplier can be treated as a random-access object [(example2)](./C++Examples/example2.cpp):

```c++
RDKit::SDMolSupplier mol_supplier( "data/5ht3ligs.sdf" , true );

for( int i = int( mol_supplier.length() ) - 1 ; i >= 0  ; --i ) {
	RDKit::ROMol *mol = mol_supplier[i];
	std::cout << mol->getProp<std::string>( "_Name" ) << " has "
	          << mol->getNumAtoms() << " atoms." << std::endl;
	delete mol;
}
```
gives
```
mol-732 has 26 atoms.
mol-15 has 24 atoms.
mol-54 has 24 atoms.
mol-295 has 20 atoms.
```

A good practice is to test each molecule to see if it was correctly
read before working with it [(example2)](./C++Examples/example2.cpp):

```c++
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
```

An alternative type of Supplier, the `RDKit::ForwardMolSupplier`
can be used to read from file-like objects.  This allows the reading
of compressed files, using, for example, the `boost::iostreams`
objects [(example2)](./C++Examples/example2.cpp):

```c++
boost::iostreams::filtering_istream ins;
ins.push( boost::iostreams::gzip_decompressor() );
ins.push( boost::iostreams::file_source( "data/actives_5ht3.sdf.gz" ) );

RDKit::ForwardSDMolSupplier forward_supplier( &ins , true );
while( !forward_supplier.atEnd() ) {
  mol = forward_supplier.next();
  std::cout << mol->getProp<std::string>( "_Name" ) << " has " << mol->getNumAtoms() << " atoms." << std::endl;
  delete mol;
}
```

Note that the forward suppliers cannot be used in random-access mode,
and a compile-time error will result if you attempt to [(example2)](./C++Examples/example2.cpp):

```
   error: no match for ‘operator[]’ (operand types are
   ‘RDKit::ForwardSDMolSupplier’ and ‘int’) 
      mol = forward_supplier[1];
```

### Writing molecules

Single molecules can be converted to text using several functions
present in the `RDKit` namespace.

For example, for SMILES [(example3)](./C++Examples/example3.cpp):

```c++
#include <GraphMol/SmilesParse/SmilesWrite.h>
.
.
RDKit::ROMol *mol = RDKit::MolFromMolFile( "data/chiral.mol" );
std::cout << RDKit::MolToSmiles( *mol ) << std::endl;
```
gives
```
CC(O)c1ccccc1
```
and [(example3)](./C++Examples/example3.cpp)
```c++
std::cout << RDKit::MolToSmiles( *mol , true ) << std::endl;
```
produces
```
C[C@H](O)c1ccccc1
```
where the `true` in the second function call specifies that isomeric
SMILES should be produced.
Note that the SMILES produced is canonical, so the output should be
the same no matter how a particular molecule is input.  For example
[(example3)](./C++Examples/example3.cpp)

```c++
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "C1=CC=CN=C1" );
std::cout << RDKit::MolToSmiles( *mol1 ) << std::endl;

RDKit::ROMol *mol2 = RDKit::SmilesToMol( "c1cccnc1" );
std::cout << RDKit::MolToSmiles( *mol2 ) << std::endl;

RDKit::ROMol *mol3 = RDKit::SmilesToMol( "n1ccccc1" );
std::cout << RDKit::MolToSmiles( *mol3 ) << std::endl;
```

all produce `c1ccncc1` as output.

If you'd like to have the Kekule form of the SMILES, you need to Kekulize
a RWMol copy of the molecule, using the Kekulize function declared in
MolOps.h [(example3)](./C++Examples/example3.cpp):

```c++
#include <GraphMol/MolOps.h>
.
.
RDKit::RWMol *mol4 = new RDKit::RWMol( *mol );
RDKit::MolOps::Kekulize( *mol4 );
std::cout << RDKit::MolToSmiles( *mol4 ) << std::endl;
```
gives
```
CC(O)C1=CC=CC=C1
```

Note: as of this writing (Aug 2008), the SMILES provided when one
requests kekuleSmiles are not canonical. The limitation is not in the
SMILES generation, but in the kekulization itself.

MDL Mol blocks are also available [(example3)](./C++Examples/example3.cpp):

```c++
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "C1=CC=CN=C1" );
std::cout << RDKit::MolToMolBlock( *mol1 ) << std::endl;
```
gives
```

    RDKit
	
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
```

To include names in the mol blocks, set the molecule's “\_Name”
property [(example3)](./C++Examples/example3.cpp):
  
```c++
mol1 = RDKit::SmilesToMol( "C1CCC1" );
mol1->setProp( "_Name" , "cyclobutane" );
std::cout << RDKit::MolToMolBlock( *mol1 ) << std::endl;
```
gives
```
cyclobutane
     RDKit
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
```

Note that setProp, which is a general function, can be called on an
ROMol as well as an RWMol, which came as a surprise to me.

In order for atom or bond stereochemistry to be recognised correctly by most
software, it's essential that the Mol block have atomic coordinates.
It's also convenient for many reasons, such as drawing the molecules.

You can either include 2D coordinates (i.e. a depiction), using the
function in the RDDepict namespace and declared in RDDepictor.h
[(example4)](./C++Examples/example4.cpp):

```c++
#include <GraphMol/Depictor/RDDepictor.h>
.
.
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "C1CCC1" );
RDDepict::compute2DCoords( *mol1 );
std::cout << RDKit::MolToMolBlock( *mol1 ) << std::endl;
```
gives
```

     RDKit          2D

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
```

Or you can add 3D coordinates by embedding the molecule [(example4)](./C++Examples/example4.cpp):

```c++
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
.
.
RDKit::ROMol *mol2 = RDKit::SmilesToMol( "C1CCC1" );
mol2->setProp( "_Name" , "cyclobutane3D" );
RDKit::DGeomHelpers::EmbedMolecule( *mol2 );
RDKit::MMFF::MMFFOptimizeMolecule( *mol2 , 1000 , "MMFF94s" );
std::cout << RDKit::MolToMolBlock( *mol2 ) << std::endl;
```
gives
```
cyclobutane3D
     RDKit          3D

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
```

The optimization step isn't necessary, but it substantially improves
the quality of the conformation.

To get good 3D conformations, it's almost always a good idea to add
hydrogens to the molecule first [(example4)](./C++Examples/example4.cpp):

```c++
RDKit::ROMol *mol3 = RDKit::MolOps::addHs( *mol2 );
RDKit::MMFF::MMFFOptimizeMolecule( *mol3 , 1000 , "MMFF94s" );

RDKit::RWMol *mol4 = new RDKit::RWMol( *mol3 );
RDKit::MolOps::addHs( *mol4 );
```

<a name="twoAddHs"></a>Note that there are 2 overloaded versions of
addHs. The first takes an
ROMol and, because that can't be edited, returns a pointer to a new
ROMol with the result.  If you use this version be careful not to leak
memory by not deleting mol2 when you are finished with it. The second
takes an RWMol which it is able to modify in place.  With shared
pointers, memory leaks can be avoided [(example4)](./C++Examples/example4.cpp):
```c++
RDKit::ROMOL_SPTR mol3sp( RDKit::MolOps::addHs( *mol2 ) );
mol3sp->setProp( "_Name" , "cyclobutaneSP" );
RDKit::MMFF::MMFFOptimizeMolecule( *mol3sp , 1000 , "MMFF94s" );
```

Once the optimisation is complete, the hydrogens can be removed
again [(example4)](./C++Examples/example4.cpp):

```c++
RDKit::ROMol *mol5 = RDKit::MolOps::removeHs( *mol3 );
RDKit::MolOps::removeHs( *mol4 );
```

Again, there are two versions, one of which has an opportunity for a
memory leak.

If you'd like write the molecules to file, use the normal C++
streams [(example4)](./C++Examples/example4.cpp):

```c++
    #include <fstream>
    .
    .
    std::ofstream ofs( "data/foo.mol" );
    ofs << RDKit::MolToMolBlock( *mol5 );
```

### Writing sets of molecules

Multiple molecules can be written to a file using an object of a
concrete subclass of the `MolWriter` class [(example5)](./C++Examples/example5.cpp):

```c++
#include <GraphMol/FileParsers/MolWriters.h>
.
.
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
```

A MolWriter can also be initialised to a file-like object, so
compressed files can be written or molecules can be written to a
string in memory [(example5)](./C++Examples/example5.cpp):

```c++
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
```

Other available writers include SmilesWriter and TDTWriter (for those
of you with an interest in historical Cheminformatics!)

## Working with Molecules

### Looping over Atoms and Bonds

Once you have a molecule, it's relatively easy to loop over its atoms
and bonds so long as you remember that an atom is a Vertex and a bond
is an Edge and accept the odd syntax [(example6)](./C++Examples/example6.cpp):

```c++
RDKit::ROMOL_SPTR mol( RDKit::SmilesToMol( "C1OC1" ) );
RDKit::ROMol::VERTEX_ITER it , end;
boost::tie( it , end ) = mol->getVertices();
while( it != end ) {
    const RDKit::Atom *atom = (*mol)[*it].get();
    std::cout << atom->getAtomicNum() << std::endl;
    ++it;
}
```
gives
```
6 8 6
```
An alternative method uses the fact that atoms and bonds can be
selected by index number [(example6)](./C++Examples/example6.cpp):
```c++
for( unsigned int i = 0 , is = mol->getNumAtoms() ; i < is ; ++i ) {
    const RDKit::Atom *atom = mol->getAtomWithIdx( i ); 
    std::cout << atom->getAtomicNum() << std::endl;
}
```
Likewise with bonds [(example6)](./C++Examples/example6.cpp):
```c++
RDKit::ROMol::EDGE_ITER bond_it , bond_end;
boost::tie( bond_it , bond_end ) = mol->getEdges();
while( bond_it != bond_end ) {
	const RDKit::Bond *bond = (*mol)[*bond_it].get();
	std::cout << bond->getBondType() << std::endl;
    ++bond_it;
}

for( unsigned int i = 0 , is = mol->getNumBonds() ; i < is ; ++i ) {
	const RDKit::Bond *bond = mol->getBondWithIdx( i ); 
    std::cout << bond->getIsAromatic() << std::endl;   
}
```
gives
```
1 1 1
0 0 0
```

A bond can be specified by the atoms at its ends, with
a zero pointer being returned if there isn't one [(example6)](./C++Examples/example6.cpp):
```c++
RDKit::ROMOL_SPTR mol2( RDKit::SmilesToMol( "C1OC1Cl" ) );
const RDKit::Bond *bond = mol2->getBondBetweenAtoms( 0 , 1 );
std::cout << bond->getBeginAtomIdx() << " to "
          << bond->getBeginAtomIdx() << " is "
		  << bond->getBondType() << std::endl;
if( !mol2->getBondBetweenAtoms( 0 , 3 ) ) {
    std::cout << "No bond between 0 and 3" << std::endl;
}
```

The neighbours of an atom can also be extracted, but note that you
need an ADJ\_ITER rather than a VERTEX\_ITER [(example6)](./C++Examples/example6.cpp):
```c++
const RDKit::Atom *atom = mol2->getAtomWithIdx( 2 );
RDKit::ROMol::ADJ_ITER nbr , end_nbr;
boost::tie( nbr , end_nbr ) = mol2->getAtomNeighbors( atom );
while( nbr != end_nbr ) {
    const RDKit::Atom *nbr_atom = (*mol2)[*nbr].get();
    std::cout << nbr_atom->getIdx() << " : " << nbr_atom->getAtomicNum() << std::endl;
    ++nbr;
}
```
gives
```
1 : 8
3 : 17
0 : 6
```

### Ring Information

It is relatively easy to obtain ring information for atoms and bonds
[(example7)](./C++Examples/example7.cpp):
```c++
#include <GraphMol/MolOps.h>
.
.
RDKit::ROMOL_SPTR mol( RDKit::SmilesToMol( "OC1C2C1CC2" ) );

if( !mol->getRingInfo()->isInitialized() ) {
    RDKit::MolOps::findSSSR( *mol );
}
for( unsigned int i = 0 , is = mol->getNumAtoms() ; i < is ; ++i ) {
	const RDKit::Atom *atom = mol->getAtomWithIdx( i );
    std::cout << mol->getRingInfo()->numAtomRings( atom->getIdx() ) << " ";
}
std::cout << std::endl;

for( unsigned int i = 0 , is = mol->getNumBonds() ; i < is ; ++i ) {
    const RDKit::Bond *bond = mol->getBondWithIdx( i );
	std::cout << mol->getRingInfo()->numBondRings( bond->getIdx() ) << " ";
}
std::cout << std::endl;
```
gives
```
0 1 2 1 1 1
0 1 2 1 1 1 1
```
Obviously, findSSSR only needs to be called once for the molecule. If
you only need to know whether the atom or bond is in a ring, just test
whether or not the return value is zero [(example7)](./C++Examples/example7.cpp):
```c++
const RDKit::Bond *bond = mol->getBondWithIdx( 1 );
if( mol->getRingInfo()->numBondRings( bond->getIdx() )) {
    std::cout <<  "Bond " << bond->getIdx() << " is in a ring" << std::endl;;
}

```
gives
```
Bond 1 is in a ring
```
Other information about presence in smallest rings can also be
obtained from the RingInfo object of the molecule [(example7)](./C++Examples/example7.cpp):
```c++
std::cout << "Atom 2 is in ring of size 3 : "
          << mol->getRingInfo()->isAtomInRingOfSize( 2 , 3 ) << std::endl;
std::cout << "Atom 2 is in ring of size 4 : "
	      << mol->getRingInfo()->isAtomInRingOfSize( 2 , 4 ) << std::endl;
std::cout << "Atom 2 is in ring of size 5 : "
	      << mol->getRingInfo()->isAtomInRingOfSize( 2 , 5 ) << std::endl;
std::cout << "Bond 1 is in ring of size 3 : "
	      << mol->getRingInfo()->isBondInRingOfSize( 1 , 3 ) << std::endl;

```
gives
```
Atom 2 is in ring of size 3 : 1
Atom 2 is in ring of size 4 : 1
Atom 2 is in ring of size 5 : 0
Bond 1 is in ring of size 3 : 1
```
More detail about the smallest set of smallest rings (SSSR) is
available [(example7)](./C++Examples/example7.cpp):
```c++
RDKit::VECT_INT_VECT rings;
RDKit::MolOps::symmetrizeSSSR( *mol , rings );
std::cout << "Number of symmetric SSSR rings : " << rings.size() << std::endl;
for( auto it1 = rings.begin() , it1_end = rings.end() ; it1 != it1_end ; ++it1 ) {
	for( auto it2 = it1->begin() , it2_end = it1->end() ; it2 != it2_end ; ++it2 ) {
      std::cout << *it2 << " ";
    }
    std::cout << std::endl;
}
```
gives
```
Number of symmetric SSSR rings : 2
1 2 3 
4 5 2 3
```
As the name suggests, this is a symmetrized SSSR; if you are
interested in the number of "true" SSSR, use the `findSSSR` function
[(example7)](./C++Examples/example7.cpp):
```c++
std::cout << "Number of SSSR rings : " << RDKit::MolOps::findSSSR( *mol ) << std::endl;

```
gives
```
2
```
The distinction between symmetrized and non-symmetrized SSSR is
discussed in more detail below in the section
[The SSSR Problem](#TheSSSRProblem).

## Modifying molecules

Normally molecules are stored in the RDKit with the hydrogen atoms
implicit (i.e. not explicitly present in the molecular graph).  When
it is useful to have the hydrogens explicitly persent, for example
wehn generating or optimizing the 3D geometry, the
`RDKit::MolOps::addHs` function can be used
[(example8)](./C++Examples/example8.cpp). 
```c++
RDKit::ROMOL_SPTR mol1( RDKit::SmilesToMol( "CCO" ) );
std::cout << "Number of atoms : " << mol1->getNumAtoms() << std::endl;
RDKit::ROMOL_SPTR mol2( RDKit::MolOps::addHs( *mol1 ) );
std::cout << "Number of atoms : " << mol2->getNumAtoms() << std::endl;
```
gives
```
Number of atoms : 3
Number of atoms : 9
```
Recall that there are two versions of `RDKit::MolOps::addHs`, as
described [above](#twoAddHs).
The Hs can be removed again using the `RDKit::MolOps::RemoveHs`
function, which again has two forms
[(example8)](./C++Examples/example8.cpp):
```c++
RDKit::RWMOL_SPTR mol3( new RDKit::RWMol( *mol2 ) );
RDKit::MolOps::removeHs( *mol3 );
std::cout << "Number of atoms : " << mol3->getNumAtoms() << std::endl;
```
which returns the atom count to 3.

RDKit molecules are usually stored with the bonds in aromatic rings
having aromatic bond types. This can be changed with the
`RDKit::MolOps::Kekulize` function, which must be called with an RWMol
[(example9)](./C++Examples/example9.cpp):
```c++
  RDKit::RWMOL_SPTR mol( new RDKit::RWMol( *RDKit::SmilesToMol( "c1ccccc1" ) ) );
  std::cout << "Order : " << mol->getBondWithIdx( 0 )->getBondType() << std::endl;
  std::cout << "Aromatic : " << mol->getBondWithIdx( 0 )->getIsAromatic() << std::endl;

  RDKit::MolOps::Kekulize( *mol );
  std::cout << "After default Kekulize : Order : " << mol->getBondWithIdx( 0 )->getBondType() << std::endl;
  std::cout << "After default Kekulize : Aromatic : " << mol->getBondWithIdx( 0 )->getIsAromatic() << std::endl;
```
gives
```
Order : 12
Aromatic : 1
After default Kekulize : Order : 2
After default Kekulize : Aromatic : 0
```
The bond orders are defined as the enum BondType in
[Bond.h](../../Code/GraphMol/Bond.h), and an aromatic bond
currenly has the value 12.
Note that by default, the Kekulize function clears the aromatic flags
on the atoms and bonds. **This is in contrast to the Python version of
Kekulize, which preserves the flags by default.**  The behaviour can be
forced explicitly [(example9.cpp)](./C++Examples/example9.cpp): 
```c++
RDKit::RWMOL_SPTR mol1( new RDKit::RWMol( *RDKit::SmilesToMol( "c1ccccc1" ) ) );
RDKit::MolOps::Kekulize( *mol1 , false );
std::cout << "After Kekulize, markAtomsBonds false : Aromatic : " << mol1->getBondWithIdx( 0 )->getIsAromatic() << std::endl;

RDKit::RWMOL_SPTR mol2( new RDKit::RWMol( *RDKit::SmilesToMol( "c1ccccc1" ) ) );
RDKit::MolOps::Kekulize( *mol2 , true );
std::cout << "After Kekulize, markAtomsBonds true : Aromatic : " << mol2->getBondWithIdx( 0 )->getIsAromatic() << std::endl;
```
gives
```
After Kekulize, markAtomsBonds false : Aromatic : 1
After Kekulize, markAtomsBonds true : Aromatic : 0
```

Bonds can be restored to the aromatic bond type using the
`RDKit::MolOps::sanitizeMol` function:
```c++
RDKit::MolOps::sanitizeMol( *mol );
std::cout << "Order : " << mol->getBondWithIdx( 0 )->getBondType() << std::endl;
std::cout << "Aromatic : " << mol->getBondWithIdx( 0 )->getIsAromatic() << std::endl;
```
gives
```
Order : 12
Aromatic : 1
```
once more.

## Working with 2D molecules: Generating Depictions

The RDKit has a library for generating depictions (sets of 2D)
coordinates for molecules.  This library, which is part of the
RDDepict namespace, is accessed via the `RDDepict::Compute2DCoords`
function [(example10.cpp)](./C++Examples/example10.cpp):
```c++
#include <GraphMol/Depictor/RDDepictor.h>
.
.
RDKit::RWMOL_SPTR mol( new RDKit::RWMol( *RDKit::SmilesToMol( "c1nccc2n1ccc2" ) ) );
RDDepict::compute2DCoords( *mol );
```

The 2D conformation is constructed to minimize intramolecular clashes,
i.e. to maximize the clarity of the drawing.  Unlike the Python
equivalent, the depiction is not placed in a canonical orientation by
default. This can be forced by passing `true` as the third parameter
[(example10.cpp)](./C++Examples/example10.cpp):
```c++
#include <Geometry/point.h>
.
.
RDDepict::compute2DCoords( *mol , static_cast<RDGeom::INT_POINT2D_MAP *>( 0 ) ,
	                       true );
```

The `point.h` must be included for the typedef that defines
`INT_POINT2D_MAP`, of which more later..
By default, all existing conformations are removed when the 2D
coordinates are created.  This can be changed by passing false as a
4th parameter.  The 2D coordinates are added as another conformation
of the molecule so it's a bit tricky combining them both in the same
molecule, and probably best avoided.

The Python API has a convenience function
`GenerateDepictionMatching2DStructure` which forces the 2D coordinate
generation to orientate molecules according to a template structure.
There is currently no such function in C++, but it is relatively
simple to reproduce the effects
[(example10.cpp)](./C++Examples/example10.cpp):
```c++
#include <GraphMol/Substruct/SubstructMatch.h>
.
.
RDKit::ROMol *templ = RDKit::SmilesToMol( "c1nccc2n1ccc2" );
RDDepict::compute2DCoords( *templ );
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "c1cccc2ncn3cccc3c21" );

RDKit::MatchVectType matchVect;
if( RDKit::SubstructMatch( *mol1 , *templ , matchVect ) ) {
  RDKit::Conformer &conf = templ->getConformer();
  RDGeom::INT_POINT2D_MAP coordMap;
  for( RDKit::MatchVectType::const_iterator mv = matchVect.begin() ;
	 mv != matchVect.end() ; ++mv ) {
    RDGeom::Point3D pt3 = conf.getAtomPos( mv->first );
    RDGeom::Point2D pt2( pt3.x , pt3.y );
    coordMap[mv->second] = pt2;
  }
  RDDepict::compute2DCoords( *mol1 , &coordMap );
}
```
Here, `coordMap` maps the coordinates of atoms in the target
molecule templ onto corresponding atoms in the reference molecule.

It is also possible to produce a 2D picture that attempts to mimic as
closely as possible a 3D conformation.  Again, an equivalent of the
Python function
`rdkit.Chem.AllChem.GenerateDepictionMatching3DStructure` doesn't
exist in C++, but if you examine
[AllChem.py](../../rdkit/Chem/AllChem.py) you will get an idea of how
to reproduce the effect.

## Working with 3D Molecules

The RDKit can generate conformations for molecules using two different
methods.  The original method used distance geometry. [#blaney]_
The algorithm followed is:

1. The molecule's distance bounds matrix is calculated based on the
   connection table and a set of rules. 

2. The bounds matrix is smoothed using a triangle-bounds smoothing
   algorithm. 

3. A random distance matrix that satisfies the bounds matrix is
   generated. 

4. This distance matrix is embedded in 3D dimensions (producing
   coordinates for each atom). 

5. The resulting coordinates are cleaned up somewhat using a crude
   force field and the bounds matrix. 

Note that the conformations that result from this procedure tend to be
fairly ugly. They should be cleaned up using a force field.
This can be done within the RDKit using its implementation of the
Universal Force Field (UFF). [#rappe]_ 

More recently, there is an implementation of the method of Riniker and
Landrum [#riniker2]_ which uses torsion angle preferences from the
Cambridge Structural Database (CSD) to correct the conformers after
distance geometry has been used to generate them.  With this method,
there should be no need to use a minimisation step to clean up the
structures.

The full process of embedding and optimizing a molecule is easier than
all the above verbiage makes it sound
[(example11.cpp)](./C++Examples/example11.cpp):
```c++
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
.
.
RDKit::ROMOL_SPTR mol( RDKit::SmilesToMol( "C1CCC1OC" ) );
RDKit::ROMOL_SPTR mol1( RDKit::MolOps::addHs( *mol ) );
// Original distance geometry embedding
RDKit::DGeomHelpers::EmbedMolecule( *mol1 , 0 , 1234 );
RDKit::UFF::UFFOptimizeMolecule( *mol1 );

// new Riniker and Landrum CSD-based method
RDKit::ROMOL_SPTR mol2( RDKit::MolOps::addHs( *mol ) );
RDKit::DGeomHelpers::EmbedMolecule( *mol2 , 0 , 1234 , true , false ,
                                    2.0 , true , 1 ,
                                    static_cast<const std::map<int,RDGeom::Point3D> *> ( 0 ) ,
                                    1e-3 , false , true , true , true );
```
In the second call to `EmbedMolecule`, it is the last two `true`
parameters that enforce the second embedding method.  Apart from the
random number seed (1234) the other parameters are set to default
values.  Setting the random number seed to other than the default -1
ensures that the same conformations are produced each time the code is
run.  This is convenient when testing to ensure reproducibility of
results but disguises the inherently non-deterministic nature of the
algorithm.

The RDKit also has an implementation of the MMFF94 force field
available. [#mmff1]_, [#mmff2]_, [#mmff3]_, [#mmff4]_, [#mmffs]
[(example11.cpp)](./C++Examples/example11.cpp):
```c++
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
.
.
RDKit::MMFF::MMFFOptimizeMolecule( *mol2 , 1000 , "MMFF94s" );
```
Please note that the MMFF atom typing code uses its own aromaticity model,
so the aromaticity flags of the molecule will be modified after calling
MMFF-related methods.
Note the calls to `RDKit::MolOps::addHs()` in the examples above. By
default RDKit molecules do not have H atoms explicitly present in the
graph, but they are important for getting realistic geometries, so
they generally should be added.  They can always be removed afterwards
if necessary with a call to `RDKit::MolOps::removeHs()`

With the RDKit, multiple conformers can also be generated using the two
different embedding methods. In both cases this is simply a matter of
running the distance geometry calculation multiple times from
different random start points. The 2nd parameter to
`EmbedMultipleConfs` allows the user to
set the number of conformers that should be generated.  Otherwise the
procedures are similar to before
[(example11.cpp)](./C++Examples/example11.cpp):
```c++
RDKit::INT_VECT mol1_cids = RDKit::DGeomHelpers::EmbedMultipleConfs( *mol1 , 10 );
std::cout << "Number of conformations : " << mol1_cids.size() << std::endl;

RDKit::INT_VECT mol2_cids;
RDKit::DGeomHelpers::EmbedMultipleConfs( *mol2 , mol2_cids , 20 , 1 , 30 , 1234 ,
                                         true , false , 2.0 , true , 1 , -1.0 , 
                                         static_cast<const std::map<int,RDGeom::Point3D> *> ( 0 ) ,
                                         1e-3 , false , true , true , true );
std::cout << "Number of conformations : " << mol2_cids.size()
          << std::endl;
```
The conformer ids are returned in `mol1_cids` and `mol2_cids` and
there are two overloaded functions with different ways of supplying
the information.  As before, the CSD-based method is invoked by that
last two `true` parameters, and in the example above the the default
number of conformations to be produced has been changed from 10 to 20.
The conformers so generated can be aligned
to each other and the RMS values calculated
[(example11.cpp)](./C++Examples/example11.cpp):
```c++
#include <GraphMol/MolAlign/AlignMolecules.h>
.
.
std::vector<double> rms_list;
std::vector<unsigned int> m2cids( mol2_cids.begin() , mol2_cids.end() );
RDKit::MolAlign::alignMolConformers( *mol2 ,
                                     static_cast<const std::vector<unsigned int> *>( 0 ) ,
                                     &m2cids ,
                                     static_cast<const RDNumeric::DoubleVector *>( 0 ) ,
                                     false , 50 , &rms_list );
```
The RMS values for the overlays will be fed into rms_list on return.
Note the somewhat inconvenient issue that `EmbedMultipleConfs` returns
a vector of `ints` for the conformer ids, but `alignMolConformers`
requires a vector of `unsigned ints`. The first vector of `unsigned
ints` in the `alignMolConformers` declaration is atom ids, and allowd
the alignment to be performed on just a subset of atoms which can be
convenient for overlaying a core and seeing how the other bits of the
molecule varied in the different conformations.

There is no C++ equivalent to the Python function
`AllChem.GetConformerRMS()` to compute the RMS between two specific
conformers (e.g. 1 and 9).

*Disclaimer/Warning*: Conformation generation is a difficult and
subtle task. The original, default, 2D->3D conversion provided with
the RDKit is not intended to be a replacement for a “real”
conformational analysis tool; it merely provides quick 3D structures
for cases when they are required. On the other hand, the second
method, when a sufficiently large number of conformers are generated,
should be adequate for most purposes.

## <a name="TheSSSRProblem"></a>The SSSR Problem

As others have ranted about with more energy and eloquence than I
intend to, the definition of a molecule's smallest set of smallest
rings is not unique.  In some high symmetry molecules, a “true” SSSR
will give results that are unappealing.  For example, the SSSR for
cubane only contains 5 rings, even though there are
“obviously” 6. This problem can be fixed by implementing a *small*
(instead of *smallest*) set of smallest rings algorithm that returns
symmetric results.  This is the approach that we took with the RDKit.

Because it is sometimes useful to be able to count how many SSSR rings
are present in the molecule, there is a
`rdkit.Chem.rdmolops.GetSSSR` function, but this only returns the
SSSR count, not the potentially non-unique set of rings.


