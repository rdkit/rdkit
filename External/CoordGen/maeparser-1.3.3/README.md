# maeparser
[![Azure_Build_Status](https://dev.azure.com/patlorton/maeparser/_apis/build/status/schrodinger.maeparser?branchName=master)](https://dev.azure.com/patlorton/maeparser/_build/latest?definitionId=1&branchName=master)
[![Build_Status](https://ci.appveyor.com/api/projects/status/github/schrodinger/maeparser?branch=master&svg=true)](https://ci.appveyor.com/project/torcolvin/maeparser)

maeparser is a parser for [Schrodinger](https://www.schrodinger.com/) Maestro
files.

Structure files (.mae,.maegz,.mae.gz) can contain multiple structures
delimited by "f_m_ct".  See [MaeConstants.hpp](MaeConstants.hpp) 
for standard block and property names.

To read a structure,

```C++
#include "Reader.hpp"

...

FILE* f = fopen("test.mae", "r");
schrodinger::mae::Reader r(f);

std::shared_ptr<schrodinger::mae::Block> b;
while ((b = r.next(schrodinger::mae::CT_BLOCK)) != nullptr) {
  // Parse structure
}
fclose(f);
```

See also test/UsageDemo.cpp, which reads an example structure and
stores it in a dummy Molecule class.

Background
==========

Why do we have a recursive descent parser for mae files instead of using a
parser generator like ANTLR or lex/yacc? The main reasons are that the mae
format language is 1) pretty simple, 2) unlikely to change significantly,
and 3) not a context free grammar. In addition, speed of parsing these
files is important.

In what way is the current version of the language not a CFG? Special tokens
like block opener `{` and key/value separator `:::` can also be string
values because the quotes on string values are not required. This results in
complication and pain in attempts to define a grammar.

Why this format?
=================

There are many molecule formats out there, and the significant strength of this
one is that it exactly fits the use case of Schrödinger's physics-based
modeling.  As the primary (and only lossless) Schrödinger output format, any
package that wishes to implement lossless data extraction of Maestro output
needs to interact directly with this format.  This is not an intentional
limitation, but is due to the nature of chemical storage formats: it's
extremely hard to get a format that's both 1) Flexible enough to hold any type
of data and 2) Not so flexible that each user has to implement their own rules.
The Maestro format avoids this paradox by having the exact flexibility
that Schrödinger's physics based backends require, without additional
flexibility that other use cases might demand.

In supporting Schrödinger's backend suite, maeparser is able
to handle output from:
* Molecular Dynamics applications, such as Desmond and FEP+
* Ligand-Protein Docking applications, such as Glide
* Homology Modeling and folding applications, such as Prime
* Ligand-based search applications, such as Phase and Phase Shape
* Quantum Mechanics applications, such as Jaguar
* Protein-Protein Docking applications
* Many other backends used in both Life and Material Sciences

Installation
============

Command line installation on a
Unix-like operating system follows a typical configure, build, and test procedure.
Configuration is via [CMake](https://cmake.org/).  Here is an
example command sequence:

```bash
git clone git@github.com:schrodinger/maeparser.git maeparser

# Set up the build configuration
cd maeparser
mkdir build
cd build
export CC=gcc
cmake --verbose ..

# Build it
make

# Run the custom testing
ctest
```

Defining CC ensures that the specified
compiler is used in the build
(in the example, it will be the first instance of `gcc` in one's
PATH), 
and the `--verbose` argumentenables viewing the gory details of compiling and
linking that will be necessary for debugging or reporting issues if the build fails.
