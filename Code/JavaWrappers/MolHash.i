/*
 *  Copyright (c) 2019 Greg Landrum
 *  All rights reserved.
 *
 *  This file is part of the RDKit.
 *  The contents are covered by the terms of the BSD license
 *  which is included in the file license.txt, found at the root
 *  of the RDKit source tree.
 */

%include "std_string.i"

%{
#include <GraphMol/MolHash/MolHash.h>
%}

%ignore RDKit::MolHash::generateMoleculeHashCode;
%ignore RDKit::MolHash::CodeFlags;
%ignore RDKit::MolHash::fillAtomBondCodes;
%ignore RDKit::MolHash::HashSet;
%ignore RDKit::MolHash::generateMoleculeHashSet;
%ignore RDKit::MolHash::encode;

%include<GraphMol/MolHash/nmmolhash.h>

