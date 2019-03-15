/*
*
*  Copyright (c) 2019, Greg Landrum
*  All rights reserved.
*
*  This file is part of the RDKit.
*  The contents are covered by the terms of the BSD license
*  which is included in the file license.txt, found at the root
*  of the RDKit source tree.
*
*/
%{
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
%}

%template(SparseIntVect64) RDKit::SparseIntVect<boost::int64_t>;

%template(StringMolMap) std::map<std::string, boost::shared_ptr<RDKit::ROMol>>;
%template(ROMol_Vect) std::vector<boost::shared_ptr<RDKit::ROMol>>;
%template(StringMolMap_Vect) std::vector<std::map<std::string, boost::shared_ptr<RDKit::ROMol>>>;
%template(StringROMol_VectMap) std::map<std::string,std::vector<boost::shared_ptr<RDKit::ROMol>>>;

%include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
