# Getting Started with the RDKit in C++

## What is this?

This document is intended to provide an overview of how one can use
the RDKit functionality from C++. Like the 'Getting Started with the
RDKit in Python' it is not comprehensive and it's not a manual.  It is
modelled very closely on the Python version, and most of the text will
be similar if not identical.  It is a work-in-progress, and will be
added to over the coming months.

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
Visual C++.  As of the March 2018 release, RDKit uses modern C++,
which currently means up to C++14.  You will probably need to tell
your compiler that you want to use the new features.  On Clang and
gcc, use the flag `-std=c++14`.

C++ programs generally have quite a lot of extra verbiage in them
compared to, for example, python scripts, that is similar or identical
in all programs.  This would produce extraneous clutter in the example
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
and libraries in all executables this will result in slower compile
times, especially
if you are doing static linking.  When linking to the static (.a)
libraries rather than the shared-object (.so) ones, the order the
libraries appear in the linking list can be important.  See the
CMakeLists.txt file in C++Examples directory for a good order.  In
this case, the same library list is used for all examples, so some
will be unnecessary for some of the programs. The first 3 programs
don't need the Depictor and SubstructMatch libraries, for instance,
although on my Ubuntu 16.04 system, the RDGeometryLib appears to need
to be included twice. Working out which libraries need to be linked to
and in what order can involve a tedious amount of trial and error.

## Should You Use C++ or Python?

There is no doubt that it is much easier get started with
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
write in C++ and compile to a native executable.  One reason for
faster executables from C++ is that because the code is only compiled
once, it can be worth the compiler spending more time optimising the
code for speed.  In Python, where the compilation is done each time at
run-time, this overhead is less acceptable.  Writing inefficient code
is relatively easy in any language, but the C++ compiler can save you
from yourself-at higher optimisation levels, it will re-arrange loops,
factorise expressions etc., so that the final executable may be
difficult to align with the original source code.  For example
(Huw!), the gcc will change sqrt(a) * sqrt(b) to sqrt(a*b) removing an
expensive square root operation.

Another consideration is the completeness of the API.  A lot of the
higher level functionality in RDKit is developed in Python, and
back-porting to C++ occurs on a demand-driven basis.  There are
therefore examples of quite useful functionality, such as computing
the RMS differences between conformers (GetConformerRMSMatrix) that
are not available in C++. 
Of course, if this affects you you can always implement the C++ version
and submit a Pull Request. Indeed, the RMS calculation is on its way.

## Memory Management

Memory leaks (where objects are created using `new` but never
destroyed using `delete`) are particularly insidious bugs that can be
difficult to track down
([Valgrind is your friend!](http://valgrind.org)) and may be a
surprise to people used to managed-memory languages like python and
java.  In RDKit, many of the functions return pointers to molecules
and accept such pointers as arguments.  Modern C++ provides tools such
as shared 
(`std::shared_ptr`) and single use (`std::unique_ptr`) (collectively
'smart') pointers that can be very helpful in preventing memory leaks.
If an object is created and a pointer to it is stored in either of
these, then when the smart pointer goes out of scope, the object 
is automatically deleted.  Otherwise, the shared pointer is used
exactly as one would use an ordinary pointer.  The unique_ptr is a bit
different in that it cannot be copied (you will get a compiler message
such as
```
/Users/david/Tmp/rdkit-my-fork/Docs/Book/C++Examples/example1.cpp:16:8: error: call to
      implicitly-deleted copy constructor of 'std::__1::unique_ptr<RDKit::RWMol,
      std::__1::default_delete<RDKit::RWMol> >'
  auto m2 = m1;
       ^    ~~
/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../include/c++/v1/memory:2494:3: note: 
      copy constructor is implicitly deleted because 'unique_ptr<RDKit::RWMol,
      std::__1::default_delete<RDKit::RWMol> >' has a user-declared move constructor
  unique_ptr(unique_ptr&& __u) noexcept
```
if you try to) although they can be `move`d.  There is
some anecdotal evidence that VisualStudio doesn't cope with vectors of
vectors of `unique_ptr` (`std::vector<std::vector<unique_ptr<RWMol>>>`,
for example).


## The Molecule Classes

Unlike in the Python libraries, in C++ there are two different
molecule classes, `RDKit::ROMol` and `RDKit::RWMol`.  They are both
declared in GraphMol.h. ROMol (the Read-Only molecule) is used in
most instances. It can't be edited.  On those occasions where you will
need to edit the molecule, you'll need to use the RWMol (Read-Write).

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

std::string file_root = getenv( "RDBASE" );
file_root += "/Docs/Book";

std::string mol_file = file_root + "/data/input.mol";
std::shared_ptr<RDKit::ROMol> mol2( RDKit::MolFileToMol( mol_file ) );
std::cout << *mol2 << std::endl;

std::shared_ptr<RDKit::ROMol> mol3( RDKit::SmilesToMol( "Cc1cccc" ) );
```

All these return a pointer to an ROMol on success, or NULL on
failure. Obviously, the object must be deleted when finished with to
prevent memory leaks. In the example above, and henceforth in this
document, the molecules, apart from mol1, are wrapped in shared
pointers so that the objects are deleted as soon as the shared pointer
goes out of scope.

There is also an overload of `operator""` as a convenient shorthand
method of parsing a SMILES string into a molecule
[(example1)](./C++Examples/example1.cpp):
```
using namespace RDKit;
auto mol = "C[C@H](F)c1ccc(C#N)cc1"_smiles;"
``` 
Normally the keyword `auto` means you don't need to know the return
type of the function - the compiler works it out and away you go.
However, in this case it is useful to know that you get a
`std::unique_ptr<RWMol>`, which tells you it can't be copied, the
molecule can be edited, the memory won't leak and you need to use the
`->` operator to use functions on the RWMol object.

If the molecule can't be sanitized after SMILES parsing, an
`RDKit::MolSanitizeException` (derived
from `std::exception`) is thrown, and an attempt is made to provide
sensible error messages [(example1)](./C++Examples/example1.cpp):

```c++
try {
   std::shared_ptr<RDKit::ROMol> mol( RDKit::SmilesToMol( "CO(C)C" ) );
} catch( RDKit::MolSanitizeException &e ) {
   // empty catch
}
```

displays something like `[15:58:22] Explicit valence of atom # 1 O, 3,
is greater than permitted` and [(example1)](./C++Examples/example1.cpp)

```c++
try {
   std::shared_ptr<RDKit::ROMol> mol( RDKit::SmilesToMol( "c1cc1" ) );
} catch( RDKit::MolSanitizeException &e ) {
   // empty catch
}
```

displays something like: `[12:20:41] Can't kekulize mol`.

### Reading sets of molecules

Groups of molecules are read using a Supplier (for example, an
`RDKit::SDMolSupplier` or an `RDKit::SmilesMolSupplier`)
[(example2)](./C++Examples/example2.cpp):

```c++
std::unqiue_ptr<RDKit::ROMol> mol;
std::string file_root = getenv( "RDBASE" );
file_root += "/Docs/Book";
std::string sdf_file = file_root + "/data/5ht3ligs.sdf";
bool takeOwnership = true;
RDKit::SDMolSupplier mol_supplier( sdf_file , takeOwnership );
while( !mol_supplier.atEnd() ) {
  mol.reset( mol_supplier.next() );
  std::cout << mol->getProp<std::string>( "_Name" ) << " has " << mol->getNumAtoms() << " atoms." << std::endl;
}
```
gives
```
mol-295 has 20 atoms.
mol-54 has 24 atoms.
mol-15 has 24 atoms.
mol-732 has 26 atoms.
```

The supplier can be treated as a random-access object [(example2)](./C++Examples/example2.cpp):

```c++
RDKit::SDMolSupplier mol_supplier( "data/5ht3ligs.sdf" , takeOwnership );
for( int i = int( mol_supplier.length() ) - 1 ; i >= 0  ; --i ) {
  unique_ptr<RDKit::ROMol> mol( mol_supplier[i] );
  if( mol ) {
    std::cout << mol->getProp<std::string>( "_Name" ) << " has "
	          << mol->getNumAtoms() << " atoms." << std::endl;
  }
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
bool takeOwnership = true;
RDKit::SDMolSupplier *mol_supplier = new RDKit::SDMolSupplier( "data/5ht3ligs.sdf" , takeOwnership );

for( int i = int( mol_supplier->length() ) - 1 ; i >= 0 ; --i ) {
  std::unique_ptr<RDKit::ROMol> mol( (*mol_supplier)[i] );
  if( !mol ) {
    continue;
  }
  std::cout << mol->getProp<std::string>( "_Name" ) << " has "
            << mol->getNumAtoms() << " atoms." << std::endl;
}
```

An alternative type of Supplier, the `RDKit::ForwardMolSupplier`
can be used to read from file-like objects.  This allows the reading
of compressed files, using, for example, the `boost::iostreams`
objects [(example2)](./C++Examples/example2.cpp):

```c++
std::string file_root = getenv( "RDBASE" );
file_root += "/Docs/Book";
boost::iostreams::filtering_istream ins;
ins.push( boost::iostreams::gzip_decompressor() );
std::string comp_sdf_file = file_root + "/data/actives_5ht3.sdf.gz";
ins.push( boost::iostreams::file_source( comp_sdf_file ) );
// takeOwnership must be false for this, as we don't want the SDWriter trying
// to delete the boost::iostream
bool takeOwnership = false;
RDKit::ForwardSDMolSupplier forward_supplier( &ins , takeOwnership );
while( !forward_supplier.atEnd() ) {
  mol.reset( forward_supplier.next() );
  if( mol ) {
    std::cout << mol->getProp<std::string>( "_Name" ) << " has " << mol->getNumAtoms() << " atoms." << std::endl;
  }
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
std::shared_ptr<RDKit::ROMol> mol( RDKit::MolFromMolFile( "data/chiral.mol" ) );
std::cout << RDKit::MolToSmiles( *mol ) << std::endl;
```
gives
```
C[C@H](O)c1ccccc1
```
and [(example3)](./C++Examples/example3.cpp)
```c++
bool isomeric = false;
std::cout << RDKit::MolToSmiles( *mol , isomeric ) << std::endl;
```
produces
```
CC(O)c1ccccc1
```
where the `isomeric` in the second function call specifies that isomeric
SMILES should not be produced.
Note that the SMILES produced is canonical, so the output should be
the same no matter how a particular molecule is input.  For example
[(example3)](./C++Examples/example3.cpp)

```c++
std::shared_ptr<RDKit::ROMol> mol1( RDKit::SmilesToMol( "C1=CC=CN=C1" ) );
std::cout << RDKit::MolToSmiles( *mol1 ) << std::endl;

std::shared_ptr<RDKit::ROMol> mol2( RDKit::SmilesToMol( "c1cccnc1" ) );
std::cout << RDKit::MolToSmiles( *mol2 ) << std::endl;

std::shared_ptr<RDKit::ROMol> mol3( RDKit::SmilesToMol( "n1ccccc1" ) );
std::cout << RDKit::MolToSmiles( *mol3 ) << std::endl;
```

all produce `c1ccncc1` as output.

If you'd like to have the Kekule form of the SMILES, you need to Kekulize
an RWMol copy of the molecule, using the Kekulize function declared in
MolOps.h [(example3)](./C++Examples/example3.cpp):

```c++
#include <GraphMol/MolOps.h>
.
.
std::shared_ptr<RDKit::RWMol> mol4( new RDKit::RWMol( *mol ) );
RDKit::MolOps::Kekulize( *mol4 );
std::cout << RDKit::MolToSmiles( *mol4, true, true ) << std::endl;
```
gives
```
C[C@H](O)C1=CC=CC=C1
```

Note: as of March 2017, the SMILES provided when one
requests kekuleSmiles are not canonical. The limitation is not in the
SMILES generation, but in the kekulization itself.

MDL Mol blocks are also available [(example3)](./C++Examples/example3.cpp):

```c++
std::shared_ptr<RDKit::ROMol> mol1( RDKit::SmilesToMol( "C1=CC=CN=C1" ) );
std::cout << RDKit::MolToMolBlock( *mol1 ) << std::endl;
```
gives
```

     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
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

Note that setProp, which is a general function, can be called on an
ROMol as well as an RWMol, which came as a surprise to me as I had
assumed a read-only molecule would be less changeable than that.

In order for atom or bond stereochemistry to be recognised correctly by most
software, it's essential that the mol block have atomic coordinates.
It's also convenient for many reasons, such as drawing the molecules.
Generating a mol block for a molecule that does not have coordinates will, by
default, automatically cause coordinates to be generated. These are not,
however, added to the molecule, but re-generated each time.

You can either include 2D coordinates (i.e. a depiction), using the
function in the RDDepict namespace and declared in RDDepictor.h
[(example4)](./C++Examples/example4.cpp):

```c++
#include <GraphMol/Depictor/RDDepictor.h>
.
.
std::shared_ptr<RDKit::ROMol> mol1( RDKit::SmilesToMol( "C1CCC1" ) );
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
std::shared_ptr<RDKit::ROMol> mol2( RDKit::SmilesToMol( "C1CCC1" ) );
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
std::shared_ptr<RDKit::ROMol> mol3( RDKit::MolOps::addHs( *mol2 ) );
RDKit::MMFF::MMFFOptimizeMolecule( *mol3 , 1000 , "MMFF94s" );

std::shared_ptr<RDKit::RWMol> mol4( new RDKit::RWMol( *mol3 ) );
RDKit::MolOps::addHs( *mol4 );
```

Note that there are 2 overloaded versions of
addHs. The first takes an
ROMol and, because that can't be edited, returns a pointer to a new
ROMol with the result.  If you use this version be careful not to leak
memory by not deleting mol2 when you are finished with it. The second
takes an RWMol which it is able to modify in place.  With shared
pointers, memory leaks can be avoided [(example4)](./C++Examples/example4.cpp):
```c++
std::shared_ptr<RDKit::ROMol> mol3sp( RDKit::MolOps::addHs( *mol2 ) );
mol3sp->setProp( "_Name" , "cyclobutaneSP" );
RDKit::MMFF::MMFFOptimizeMolecule( *mol3sp , 1000 , "MMFF94s" );
```

Once the optimisation is complete, the hydrogens can be removed
again [(example4)](./C++Examples/example4.cpp):

```c++
std::shared_ptr<RDKit::ROMol> mol5( RDKit::MolOps::removeHs( *mol3 );
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
// takeOwnership must be false for this, as we don't want the SDWriter trying
// to delete the std::ostringstream.
takeOwnership = false;
boost::shared_ptr<RDKit::SDWriter> sdf_writer( new RDKit::SDWriter( &oss , takeOwnership ) );
for( auto mol: mols ) {
  sdf_writer->write( *(*it) );
}
std::cout << oss.str() << std::endl;
```

Other available writers include SmilesWriter and TDTWriter (for those
of you with an interest in historical Cheminformatics).

## Working with Molecules

### Looping over Atoms and Bonds

Since the adoption of c++11, looping over atoms and bonds has become
extremely simple:

```c++
std::shared_ptr<RDKit::ROMol> mol( RDKit::SmilesToMol( "C1OC1" ) );
for(auto atom: mol->atoms()) {
  std::cout << atom->getAtomicNum() << " ";
}
std::cout << std::endl;
```
gives
```
6 8 6
```
As well as  this, there are a method that uses the fact that atoms and bonds
can be selected by index number
[(example6)](./C++Examples/example6.cpp):
```c++
for( unsigned int i = 0; i < mol->getNumAtoms() ; ++i ) {
  const RDKit::Atom *atom = mol->getAtomWithIdx( i ); 
  std::cout << atom->getAtomicNum() << " ";
}
```
Likewise with bonds [(example6)](./C++Examples/example6.cpp):
```c++
for(auto bond: mol->bonds()) {
  std::cout << bond->getBondType() << " ";
}
std::cout << std::endl;
for( unsigned int i = 0 , is = mol->getNumBonds() ; i < is ; ++i ) {
  const RDKit::Bond *bond = mol->getBondWithIdx( i ); 
  std::cout << bond->getIsAromatic() << " ";   
}
std::cout << std::endl;   

```
gives
```
1 1 1
0 0 0
```

A bond can be specified by the atoms at its ends, with
a NULL pointer being returned if there isn't one [(example6)](./C++Examples/example6.cpp):
```c++
std::shared_ptr<RDKit::ROMol> mol2( RDKit::SmilesToMol( "C1OC1Cl" ) );
const RDKit::Bond *bond = mol2->getBondBetweenAtoms( 0 , 1 );
std::cout << bond->getBeginAtomIdx() << " to "
          << bond->getBeginAtomIdx() << " is "
          << bond->getBondType() << std::endl;
if( !mol2->getBondBetweenAtoms( 0 , 3 ) ) {
  std::cout << "No bond between 0 and 3" << std::endl;
}
```

The neighbours of an atom can also be extracted [(example6)](./C++Examples/example6.cpp):
```c++
  const RDKit::Atom *atom = mol2->getAtomWithIdx( 2 );
  for(const auto &nbri: make_iterator_range(mol->getAtomBonds(atom))) {
    const RDKit::Bond *bond = (*mol)[nbri];
    unsigned int nbr_idx = bond->getOtherAtomIdx(atom->getIdx());
    int nbr_atnum = bond->getOtherAtom(atom)->getAtomicNum();
    std::cout << nbr_idx << " : " << nbr_atnum << std::endl;
  }
```
gives
```
1 : 8
3 : 17
0 : 6
```
The `auto' keyword here is hiding the fact that `getBondBetweenAtoms`
returns the index of the bond, an `unsigned int`.

### Ring Information

It is relatively easy to obtain ring information for atoms and bonds
[(example7)](./C++Examples/example7.cpp):
```c++
#include <GraphMol/MolOps.h>
.
.
std::shared_ptr<RDKit::ROMol> mol( RDKit::SmilesToMol( "OC1C2C1CC2" ) );

if( !mol->getRingInfo()->isInitialized() ) {
  RDKit::MolOps::findSSSR( *mol );
}
for( unsigned int i = 0; i < mol->getNumAtoms() ; ++i ) {
  const RDKit::Atom *atom = mol->getAtomWithIdx( i );
  std::cout << mol->getRingInfo()->numAtomRings( atom->getIdx() ) << " ";
}
std::cout << std::endl;

for( unsigned int i = 0; i < mol->getNumBonds() ; ++i ) {
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
for( auto ring: rings ) {
  for( auto ringat: ring ) {
    std::cout << ringat << " ";
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

### Modifying molecules

Normally molecules are stored in the RDKit with the hydrogen atoms
implicit (i.e. not explicitly present in the molecular graph).  When
it is useful to have the hydrogens explicitly present, for example
when generating or optimizing the 3D geometry, the
`RDKit::MolOps::addHs` function can be used
[(example8)](./C++Examples/example8.cpp).
```c++
std::shared_ptr<RDKit::ROMol> mol1( RDKit::SmilesToMol( "CCO" ) );
std::cout << "Number of atoms : " << mol1->getNumAtoms() << std::endl;
std::shared_ptr<RDKit::ROMol> mol2( RDKit::MolOps::addHs( *mol1 ) );
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
std::shared_ptr<RDKit::RWMol> mol3( new RDKit::RWMol( *mol2 ) );
RDKit::MolOps::removeHs( *mol3 );
std::cout << "Number of atoms : " << mol3->getNumAtoms() << std::endl;
```
which returns the atom count to 3.

RDKit molecules are usually stored with the bonds in aromatic rings
having aromatic bond types. This can be changed with the
`RDKit::MolOps::Kekulize` function, which must be called with an RWMol
[(example9)](./C++Examples/example9.cpp):
```c++
std::shared_ptr<RDKit::RWMol> mol( new RDKit::RWMol( *RDKit::SmilesToMol( "c1ccccc1" ) ) );
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
currently has the value 12.
Note that by default the Kekulize function clears the aromatic flags
on the atoms and bonds. **This is in contrast to the Python version of
Kekulize, which preserves the flags by default.**  The behaviour can be
forced explicitly [(example9.cpp)](./C++Examples/example9.cpp):
```c++
std::shared_ptr<RDKit::RWMol> mol1( new RDKit::RWMol( *RDKit::SmilesToMol( "c1ccccc1" ) ) );
RDKit::MolOps::Kekulize( *mol1 , false );
std::cout << "After Kekulize, markAtomsBonds false : Aromatic : " << mol1->getBondWithIdx( 0 )->getIsAromatic() << std::endl;

std::shared_ptr<RDKit::RWMol> mol2( new RDKit::RWMol( *RDKit::SmilesToMol( "c1ccccc1" ) ) );
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

### Working with 2D molecules: Generating Depictions

The RDKit has 2 different methods for generating depictions (sets of 2D
coordinates) for molecules.  They are part of the
RDDepict namespace, is accessed via the `RDDepict::Compute2DCoords`
function [(example10.cpp)](./C++Examples/example10.cpp):
```c++
#include <GraphMol/Depictor/RDDepictor.h>
.
.
std::shared_ptr<RDKit::RWMol>xs mol( new RDKit::RWMol( *RDKit::SmilesToMol( "c1nccc2n1ccc2" ) ) );
RDDepict::compute2DCoords( *mol );
```
The default method is the original RDKit algorithm.  The alternative
is CoordGen, kindly contributed by Schrodinger.  The latter is
reckoned to be slower but better at congested structures and
macrocycles.  CoordGen isn't always built into the RDKit libraries, so
it's a reasonable idea to use it conditionally
[(example10.cpp)](./C++Examples/example10.cpp):  
```c++
#ifdef RDK_BUILD_COORDGEN_SUPPORT 
  RDDepict::preferCoordGen = true;
#else
  std::cout << "CoordGen support not available" << std::endl;
#endif
```
The [CMakeLists.txt](./C++Examples/CMakeLists.txt) with the examples
defines `RDK_BUILD_COORDGEN_SUPPORT`.

The 2D conformation is constructed to minimize intramolecular clashes,
i.e. to maximize the clarity of the drawing.  Unlike the Python
equivalent, the depiction is not placed in a canonical orientation by
default. This can be forced by passing `true` as the third parameter
[(example10.cpp)](./C++Examples/example10.cpp):
```c++
RDDepict::compute2DCoords( *mol , nullptr , true );
```

By default, all existing conformations are removed when the 2D
coordinates are created.  This can be changed by passing false as a
4th parameter.  The 2D coordinates are added as another conformation
of the molecule so it's a bit tricky combining them both in the same
molecule, and probably best avoided.

The Python API has a convenience function
`GenerateDepictionMatching2DStructure` which forces the 2D coordinate
generation to orientate molecules according to a template structure.
A C++ version of the function, generateDepictionMatching2DStructure
was included in late December 2016.  If that is later than the version
of RDKit you are using, then the effect can be achieved thus:
[(example10.cpp)](./C++Examples/example10.cpp):
```c++
#include <Geometry/point.h>
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
  for(auto mv: matchVect) {
    RDGeom::Point3D pt3 = conf.getAtomPos( mv.first );
    RDGeom::Point2D pt2( pt3.x , pt3.y );
    coordMap[mv.second] = pt2;
  }
  RDDepict::compute2DCoords( *mol1 , &coordMap );
}
```
Here, `coordMap` maps the coordinates of atoms in the target
molecule templ onto corresponding atoms in the reference molecule.

It is also possible to produce a 2D picture that attempts to mimic as
closely as possible a 3D conformation.  Again, an equivalent of the
Python function
`rdkit.Chem.AllChem.GenerateDepictionMatching3DStructure` was
incorporated in December 2016.

### Working with 3D Molecules

The RDKit can generate conformations for molecules using two different
methods.  The original method uses distance geometry [[1]](#blaney).
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
Universal Force Field UFF[[2]](#rappe).

More recently, there is an implementation of the method of Riniker and
Landrum [[3]](#riniker2) which uses torsion angle preferences from the
Cambridge Structural Database (CSD) to correct the conformers after
distance geometry has been used to generate them.  With this method,
there should be no need to use a minimisation step to clean up the
structures; indeed, it is often undesirable as it may move the
torsions away from the CSD-based distributions, somewhat negating the
point.

The full process of embedding and optimizing a molecule is easier than
all the above verbiage makes it sound
[(example11.cpp)](./C++Examples/example11.cpp):
```c++
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
.
.
std::shared_ptr<RDKit::ROMol> mol( RDKit::SmilesToMol( "C1CCC1OC" ) );
std::shared_ptr<RDKit::ROMol>xs mol1( RDKit::MolOps::addHs( *mol ) );
// Original distance geometry embedding
RDKit::DGeomHelpers::EmbedMolecule( *mol1 , 0 , 1234 );
RDKit::UFF::UFFOptimizeMolecule( *mol1 );

// new Riniker and Landrum CSD-based method
// using the parameters class
RDKit::DGeomHelpers::EmbedParameters params( RDKit::DGeomHelpers::ETKDG );
params.randomSeed = 1234;
RDKit::DGeomHelpers::EmbedMolecule( *mol2 , params );
```
The Riniker and Landrum method has a number of parameters that may be
altered for various reasons beyond the scope of this document. One
that you may want to alter is the random number seed; setting the
random number seed to other than the default -1 ensures that the same
conformations are produced each time the code is run.  This is
convenient when testing to ensure reproducibility of results.  To make
it easier to vary the parameters, there is the EmbedParameters class
which is initialised to the default values on construction, and whose
individual values can be varied as desired.

The RDKit also has an implementation of the MMFF94 force field
available [[4]](#mmff1), [[5]](#mmff2), [[6]](#mmff3),
[[7]](#mmff4), [[8]](#mmffs).
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
int numConfs = 20;
RDKit::DGeomHelpers::EmbedMultipleConfs( *mol2 , mol2_cids , numConfs , params );
std::cout << "Number of conformations : " << mol2_cids.size()
          << std::endl;
```
The conformer ids are returned in `mol1_cids` and `mol2_cids` and
there are two overloaded functions with different ways of supplying
the information.  As before, the CSD-based method is invoked by
EmbedParameters object, and in the example above the default
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
requires a vector of `unsigned ints`. The reason for this is that
`EmbedMultipleConfs` uses -1 to denote a failed embedding.  The first
vector of `unsigned
ints` in the `alignMolConformers` declaration is atom ids, and allows
the alignment to be performed on just a subset of atoms which can be
convenient for overlaying a core and seeing how the other bits of the
molecule varied in the different conformations.

There is no C++ equivalent to the Python function
`AllChem.GetConformerRMS()` to compute the RMS between two specific
conformers (e.g. 1 and 9) although it is coming.

It is important to remember that unless you specify a random number
seed, you will not necessarily get the same conformations each time
you run the embedding on the same molecule, especially if you only
generate a small number of conformations relative to the number of
torsions in the structure.  If it's important in your use-case that
you have a good sampling of the conformations including all the
low-energy ones, you should be sure to
specify a large maximum number of conformations.

*Disclaimer/Warning*: Conformation generation is a difficult and
subtle task. The original, default, 2D->3D conversion provided with
the RDKit is not intended to be a replacement for a “real”
conformational analysis tool; it merely provides quick 3D structures
for cases when they are required. On the other hand, the second
method, when a sufficiently large number of conformers are generated,
should be adequate for most purposes. It is probably better to ignore
the first, historical, method entirely. It is only left as the default
method to avoid breaking existing code.

### Preserving Molecules

Molecules can be preserved, or serialised, or pickled, using the class
MolPickler in the namespace of the same name:
[(example12.cpp)](./C++Examples/example12.cpp):
```c++
#include <GraphMol/MolPickler.h>
.
.
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "c1ccncc1" );
std::string pickle;
RDKit::MolPickler::pickleMol( *mol1 , pickle );
RDKit::ROMol mol2;
RDKit::MolPickler::molFromPickle( pickle , mol2 );
std::cout << RDKit::MolToSmiles( mol2 ) << std::endl;
```
Note that the string is in binary format and will appear
as gibberish if printed to a screen. The RDKit pickle format is fairly
compact and it is much, much faster to build a molecule from a pickle
than from a Mol file or SMILES string, so storing pickles of molecules
you will be working with repeatedly can be a good idea:
[(example12.cpp)](./C++Examples/example12.cpp):
```c++
// writing to pickle file
std::string smi_file = getenv("RDBASE");
smi_file += "/Code/GraphMol/test_data/canonSmiles.long.smi";
std::string pkl_name = "canonSmiles.long.pkl";

// tab-delimited file, SMILES in column 0, name in 1, no title line
RDKit::SmilesMolSupplier suppl( smi_file , "\t" , 0 , 1 , false );
std::ofstream pickle_ostream( pkl_name , std::ios_base::binary );
int write_cnt = 0;
while( !suppl.atEnd() ) {
  RDKit::ROMol *mol = suppl.next();
  RDKit::MolPickler::pickleMol( *mol , pickle_ostream);
  delete mol;
  ++write_cnt;
}
pickle_ostream.close();  
std::cout << "Wrote " << write_cnt << " molecules" << std::endl;

// reading from pickle file
std::ifstream pickle_istream( pkl_name , std::ios_base::binary );
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
```
By default, the pickling process does not preserve any
properties attached to the molecule, which includes the molecule name
(property "_Name").  This can be forced:
```c++
RDKit::MolPickler::pickleMol( *mol , pickle_ostream, RDKit::PicklerOps::AllProps );
```

### Drawing Molecules
The RDKit has some built-in functionality for drawing molecules, found
in the RDKit namespace, with header files in
`$RDBASE/Code/GraphMol/MolDraw2D`.  There is an abstract base class
MolDraw2D which defines the interface and does the drawing, with
concrete classes for drawing to SVG or PNG files and Qt and wx
widgets.  Only the SVG output is built by default, Cairo support
requires the argument `-DRDK_BUILD_CAIRO_SUPPORT=ON` to cmake, and Qt
support `-DRDK_BUILD_QT_SUPPORT=ON`. Note that the Qt code hasn't been
looked at in years, is not tested as part of the release cycle and so
is not known to work.  To create an SVG file:
[(example13.cpp)](./C++Examples/example13.cpp):
```c++
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
.
.
RDKit::SDMolSupplier mol_supplier( "data/cdk2.sdf" , true );
RDKit::ROMol *mol1 = mol_supplier.next();
RDDepict::compute2DCoords( *mol1 );
std::ofstream outs("cdk_mol1.svg");
RDKit::MolDraw2DSVG svg_drawer(300, 300, outs);
svg_drawer.drawMolecule( *mol1 );
svg_drawer.finishDrawing();
outs.close();
```
The procedure for a PNG is slightly different:
[(example13.cpp)](./C++Examples/example13.cpp):
```c++
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
.
.
RDKit::MolDraw2DCairo cairo_drawer(300, 300);
cairo_drawer.drawMolecule(*mol1);
cairo_drawer.finishDrawing();
cairo_drawer.writeDrawingText("cdk_mol1.png");
```
As in Python, you can draw a set of molecules to a grid
[(example13.cpp)](./C++Examples/example13.cpp):
 ```c++
  std::string base_smi("c1ccccc1");
  std::vector<std::string> extras = {"F", "Cl", "Br", "OC", "C(=O)O"};
  std::vector<ROMol *> mols;
  for(auto extra: extras) {
    mols.push_back(SmilesToMol(base_smi + extra));
  }
  MolDraw2DSVG drawer(750, 400, 250, 200);
  drawer.drawMolecules(mols);
  drawer.finishDrawing();
  std::ofstream grids(file_root + "/data/example_13_grid.svg");
  grids << drawer.getDrawingText();
  grids.flush();
  grids.close();
```

Atoms in a molecule can be highlighted by drawing a coloured solid or
open circle around them, and bonds likewise can have a coloured
outline applied.  An obvious use is to show atoms and bonds that have
matched a substructure query
[(example17.cpp)](./C++Examples/example17.cpp):
```c++
drawer.drawMolecule(*mol, &hit1_atoms, &hit1_bonds);
```
In the example above, `hit1_atoms` and `hit1_bonds` are vectors of
ints giving the atoms and bonds to be highlighted.  They will all be
coloured the same default colour.  It is possible to specify the
colours for individual atoms and bonds:
```c++
std::map<int, DrawColour> atom_cols;
atom_cols.insert(std::make_pair(13, DrawColour(1.0, 0.0, 0.0)));
atom_cols.insert(std::make_pair(14, DrawColour(0.0, 1.0, 0.0)));
atom_cols.insert(std::make_pair(15, DrawColour(0.0, 0.0, 1.0)));
atom_cols.insert(std::make_pair(16, DrawColour(1.0, 0.55, 0.0)));
std::map<int, DrawColour> bond_cols;
bond_cols.insert(std::make_pair(13, DrawColour(0.8, 0.8, 0.0)));
bond_cols.insert(std::make_pair(14, DrawColour(0.0, 1.0, 1.0)));
bond_cols.insert(std::make_pair(15, DrawColour(1.0, 0.0, 1.0)));
drawer.drawMolecule(*mol, &hit1_atoms, &hit1_bonds, &atom_cols, &bond_cols);
```

Atoms and bonds can also be highlighted with multiple colours if they
fall into multiple sets, for example if they are matched by more
than 1 substructure pattern.  This is too complicated to show in this
simple introduction, but there is an example at the bottom of
[example17.cpp](./C++Examples/example17.cpp).

As of version 2020.03, it is possible to add arbitrary small strings
to annotate atoms and bonds in the drawing.  The strings are added as
properties `common_properties::atomNote` and
`common_properties::bondNote` and they will be placed automatically
close to the atom or bond in question in a manner intended to minimise
their clash with the rest of the drawing.  There are 3 flags in
`MolDraw2DOptions` that will add stereo information (R/S to atoms, E/Z
to bonds) and atom and bond sequence numbers
[(example13.cpp)](./C++Examples/example13.cpp):
 ```c++
    auto m1 = "Cl[C@H](F)NC\\C=C\\C"_smiles;
    MolDraw2DSVG drawer(250, 200);
    m1->getAtomWithIdx(2)->setProp(common_properties::atomNote, "foo");
    m1->getBondWithIdx(0)->setProp(common_properties::bondNote, "bar");	
    drawer.drawOptions().addAtomIndices = true;
    drawer.drawOptions().addStereoAnnotation = true;
    drawer.drawMolecule(*m1);
    drawer.finishDrawing();
    std::string text = drawer.getDrawingText();
    std::ofstream outs(file_root + "/data/example_13_note.svg");
    outs << text;
    outs.flush();
```

### Substructure Searching

Substructure matching can be done using query molecules built from
SMARTS.
[(example14.cpp)](./C++Examples/example14.cpp):
```c++
#include <GraphMol/Substruct/SubstructMatch.h>
.
.
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "c1ccccc1O" );
RDKit::RWMol *patt = RDKit::SmartsToMol( "ccO" );
RDKit::MatchVectType res;
if( RDKit::SubstructMatch( *mol1 , *patt , res ) ) {
  std::cout << "Pattern matched molecule" << std::endl;
}
for( size_t i = 0 ; i < res.size() ; ++i ) {
  std::cout << "(" << res[i].first << "," << res[i].second << ") ";
}
std::cout << std::endl;
```
`SubstructMatch` returns a bool to flag whether there was a match, and
replaces the contents of `res` with a mapping of the atom indices in
the pattern and a set of atoms that match in the molecule. In the
above example, the output is:
```
Pattern matched molecule
(0,0)(1,5)(2,6)
```
showing that atoms 0, 5 and 6 in the phenol matched the query. If the
pattern matches multiple times (as in this case, where 4, 5, 6 is also
a match), a single arbitrary set is returned.

All possible matches can also be returned:
[(example14.cpp)](./C++Examples/example14.cpp):
```c++
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
```
This gives
```
Match 1 : (0,0)(1,5)(2,6)
Match 2 : (0,4)(1,5)(2,6)

```
It is easy to filter lists of molecules:
[(example14.cpp)](./C++Examples/example14.cpp):
```c++
RDKit::SDMolSupplier mol_supplier( "data/actives_5ht3.sdf" , true );
RDKit::RWMol *patt1 = RDKit::SmartsToMol( "c[NH1]" );
std::vector<RDKit::ROMol *> matches;
while( !mol_supplier.atEnd() ) {
  RDKit::ROMol *mol3 = mol_supplier.next();
  if( mol3 && RDKit::SubstructMatch( *mol3 , *patt1 , res ) ) {
    matches.push_back( mol3 );
  } else {
    delete mol3;
  }
}
std::cout << "There were " << matches.size() << " hits in the file." << std::endl;
```
There should be 22 matches in the file.

Substructure matching can also be done using molecules built from
SMILES instead of SMARTS:
[(example14.cpp)](./C++Examples/example14.cpp):
```c++
RDKit::ROMol *mol4 = RDKit::SmilesToMol( "C1=CC=CC=C1OC" );
RDKit::RWMol *smi_mol1 = RDKit::SmilesToMol( "CO" );
if( RDKit::SubstructMatch( *mol4 , *smi_mol1 , res ) ) {
  std::cout << "SMILES match" << std::endl;
} else {
  std::cout << "Not SMILES match" << std::endl;
}
RDKit::RWMol *smt_mol1 = RDKit::SmartsToMol( "CO" );
if( RDKit::SubstructMatch( *mol4 , *smt_mol1 , res ) ) {
  std::cout << "SMARTS match" << std::endl;
} else {
  std::cout << "Not SMARTS match" << std::endl;
}
```
But don't forget that the semantics of the two languages are not
exactly equivalent:
[(example14.cpp)](./C++Examples/example14.cpp):
```c++
RDKit::ROMol *mol4 = RDKit::SmilesToMol( "C1=CC=CC=C1OC" );
RDKit::RWMol *smi_mol2 = RDKit::SmilesToMol( "COC" );
if( RDKit::SubstructMatch( *mol4 , *smi_mol2 , res ) ) {
  std::cout << "SMILES match" << std::endl;
} else {
  std::cout << "Not SMILES match" << std::endl;
}
RDKit::RWMol *smt_mol2 = RDKit::SmartsToMol( "COC" );
if( RDKit::SubstructMatch( *mol4 , *smt_mol2 , res ) ) {
  std::cout << "SMARTS match" << std::endl;
} else {
  std::cout << "Not SMARTS match" << std::endl;
}
// Needs aromatic C
RDKit::RWMol *smt_mol3 = RDKit::SmartsToMol( "COc" );
if( RDKit::SubstructMatch( *mol4 , *smt_mol3 , res ) ) {
  std::cout << "SMARTS match" << std::endl;
} else {
  std::cout << "Not SMARTS match" << std::endl;
}
```
gives
```
SMILES match
Not SMARTS match
SMARTS match
```

### Stereochemistry in substructure matches
By default, information about stereochemistry is not used in
substructure searches:
[(example15.cpp)](./C++Examples/example15.cpp):
```c++
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "CC[C@H](F)Cl" );
RDKit::RWMol *patt1 = RDKit::SmartsToMol( "C[C@H](F)Cl" );
RDKit::MatchVectType res;
if( RDKit::SubstructMatch( *mol1 , *patt1 , res ) ) {
  std::cout << "SMARTS 1 match" << std::endl;
} else {
  std::cout << "Not SMARTS 1 match" << std::endl;
}
RDKit::RWMol *patt2 = RDKit::SmartsToMol( "C[C@@H](F)Cl" );
if( RDKit::SubstructMatch( *mol1 , *patt2 , res ) ) {
  std::cout << "SMARTS 2 match" << std::endl;
} else {
  std::cout << "Not SMARTS 2 match" << std::endl;
}
RDKit::RWMol *patt3 = RDKit::SmartsToMol( "CC(F)Cl" );
if( RDKit::SubstructMatch( *mol1 , *patt3 , res ) ) {
  std::cout << "SMARTS 3 match" << std::endl;
} else {
  std::cout << "Not SMARTS 3 match" << std::endl;
}
```
All 3 SMARTS patterns match the molecule. To use the chirality
information, you need to pass `true` as the optional fourth parameter,
corresponding to useChirality:
[(example15.cpp)](./C++Examples/example15.cpp):
```c++
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "CC[C@H](F)Cl" );
RDKit::RWMol *patt1 = RDKit::SmartsToMol( "C[C@H](F)Cl" );
RDKit::MatchVectType res;
if( RDKit::SubstructMatch( *mol1 , *patt1 , res , true , true ) ) {
  std::cout << "SMARTS 1 chiral match" << std::endl;
} else {
  std::cout << "Not SMARTS 1 chiral match" << std::endl;
}
RDKit::RWMol *patt2 = RDKit::SmartsToMol( "C[C@@H](F)Cl" );
if( RDKit::SubstructMatch( *mol1 , *patt2 , res , true , true ) ) {
  std::cout << "SMARTS 2 chiral match" << std::endl;
} else {
  std::cout << "Not SMARTS 2 chiral match" << std::endl;
}
RDKit::RWMol *patt3 = RDKit::SmartsToMol( "CC(F)Cl" );
if( RDKit::SubstructMatch( *mol1 , *patt3 , res , true , true ) ) {
  std::cout << "SMARTS 3 chiral match" << std::endl;
} else {
  std::cout << "Not SMARTS 3 chiral match" << std::endl;
}
```
gives
```
SMARTS 1 chiral match
Not SMARTS 2 chiral match
SMARTS 3 chiral match
```
Notice that when useChirality is true, a non-chiral query __does__
match a chiral molecule. The same is not true for a chiral query and a
non-chiral molecule:
[(example15.cpp)](./C++Examples/example15.cpp):
```c++
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "CC[C@H](F)Cl" );
RDKit::RWMol *mol2 = RDKit::SmilesToMol( "CC(F)Cl" );
if( RDKit::SubstructMatch( *mol1 , *mol2 , res , true , true ) ) {
  std::cout << "Chiral mol, non-chiral query : match" << std::endl;
} else {
  std::cout << "Chiral mol, non-chiral query : NO match" << std::endl;
}

RDKit::RWMol *patt5 = RDKit::SmilesToMol( "C[C@H](F)Cl" );
if( RDKit::SubstructMatch( *mol2 , *patt5 , res , true , true ) ) {
  std::cout << "Non-chiral mol, chiral query : match" << std::endl;
} else {
  std::cout << "Non-chiral mol, chiral query : NO match" << std::endl;
}
```
gives
```
Chiral mol, non-chiral query : match
Non-chiral mol, chiral query : NO match
```

### Atom Map Indices in SMARTS

It is possible to attach indices to the atoms in the SMARTS
pattern. This is most often done in reaction SMARTS (see Chemical
Reactions), but is more general than that. For example, in the SMARTS
patterns for torsion angle analysis published by Guba et
al. [[9]](#guba) indices are used to define the four atoms of the
torsion of interest. This allows additional atoms to be used to define
the environment of the four torsion atoms, as in
`[cH0:1][c:2]([cH0])!@[CX3!r:3]=[NX2!r:4]` for an aromatic C=N
torsion. We might wonder in passing why they didn’t use recursive
SMARTS for this, which would have made life easier, but it is what it
is. The atom lists from GetSubstructureMatches are guaranteed to be in
order of the SMARTS, but in this case we’ll get five atoms so we need
a way of picking out, in the correct order, the four of interest. When
the SMARTS is parsed, the relevant atoms are assigned an atom map
number property that we can easily extract:
[(example16.cpp)](./C++Examples/example16.cpp):
```c++
RDKit::RWMol *patt1 = RDKit::SmartsToMol( "[cH0:1][c:2]([cH0])!@[CX3!r:3]=[NX2!r:4]" );
std::map<int,unsigned int> ind_map;
for(auto atom: patt1->atoms()) {
  int map_num = atom->getAtomMapNum();
  if( map_num ) {
    ind_map[map_num-1] = atom->getIdx();
  }
}
std::vector<unsigned int> map_list;
for(auto im: ind_map) {
  map_list.push_back(im.second);
}
for( size_t i = 0 , is = map_list.size() ; i < is ; ++i ) {
  std::cout << map_list[i] << " ";
}
std::cout << std::endl;
```
gives
```
0 1 3 4
```
Then, when using the query on a molecule, you can get the indices of
the four matching atoms like this:
[(example16.cpp)](./C++Examples/example16.cpp):
```c++
RDKit::ROMol *mol1 = RDKit::SmilesToMol( "Cc1cccc(C)c1C(C)=NC" );
std::vector<RDKit::MatchVectType> hits_vect;
if( RDKit::SubstructMatch( *mol1 , *patt1 , hits_vect ) ) {
  for( size_t i = 0 ; i < hits_vect.size() ; ++i ) {
    std::cout << "Match " << i + 1 << " : ";
    for( size_t j = 0 ; j < map_list.size() ; ++j ) {
	std::cout << hits_vect[i][map_list[j]].second << " ";
    }
    std::cout << std::endl;
  }
}
```
gives
```
Match 1 : 1 7 8 10
```

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

## Footnotes
1. <a name="blaney"></a>Blaney, J. M.; Dixon, J. S. "Distance Geometry
in Molecular Modeling".  *Reviews in Computational Chemistry*; VCH:
New York, 1994.
2. <a name="rappe"></a>Rappé, A. K.; Casewit, C. J.; Colwell, K. S.;
Goddard III, W. A.; Skiff, W. M. "UFF, a full periodic table force
field for molecular mechanics and molecular dynamics
simulations". *J. Am. Chem. Soc.* **114**:10024-35 (1992) .
3. <a name="riniker2"></a>Riniker, S.; Landrum, G. A. "Better Informed
Distance Geometry: Using What We Know To Improve Conformation
Generation" *J. Chem. Inf. Comp. Sci.* **55**:2562-74 (2015)
4. <a name="mmff1"></a>Halgren, T. A. "Merck molecular force
field. I. Basis, form, scope, parameterization, and performance of
MMFF94." *J. Comp. Chem.* **17**:490–19 (1996).
5. <a name="mmff2"></a>Halgren, T. A. "Merck molecular force
field. II. MMFF94 van der Waals and electrostatic parameters for
intermolecular interactions." *J. Comp. Chem.* **17**:520–52 (1996).
6. <a name="mmff3"></a>Halgren, T. A. "Merck molecular force
field. III. Molecular geometries and vibrational frequencies for
MMFF94." *J. Comp. Chem.* **17**:553–86 (1996).
7. <a name="mmff4"></a>Halgren, T. A. & Nachbar, R. B. "Merck
molecular force field. IV. conformational energies and geometries
for MMFF94." *J. Comp. Chem.* **17**:587-615 (1996).
8. <a name="mmffs"></a>Halgren, T. A. "MMFF VI. MMFF94s option for
energy minimization studies." *J. Comp. Chem.* **20**:720–9
(1999).
9. <a name="guba"></a>Guba, W.; Meyder, A.; Rarey, M.; Hert,
J. "Torsion Library Reloaded: A New Version of Expert-Derived
SMARTS Rules for Assessing Conformations of Small
Molecules". *J. Chem. Inf. Model.* ** 56**:1-5 (2016)
