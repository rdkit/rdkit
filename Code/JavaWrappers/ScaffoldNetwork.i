/*
*
*  Copyright (c) 2019, Greg Landrum and T5 Informatics GmbH
*  All rights reserved.
*
*  This file is part of the RDKit.
*  The contents are covered by the terms of the BSD license
*  which is included in the file license.txt, found at the root
*  of the RDKit source tree.
*
*/
%{
#include <GraphMol/ScaffoldNetwork/ScaffoldNetwork.h>
typedef std::vector<std::string> STR_VECT;
typedef std::vector<unsigned> UINT_VECT;
%}


%template(ROMol_Vect) std::vector<boost::shared_ptr<RDKit::ROMol>>;
%template(NetworkEdge_Vect) std::vector<RDKit::ScaffoldNetwork::NetworkEdge>;
%include <GraphMol/ScaffoldNetwork/ScaffoldNetwork.h>
%template(createScaffoldNetwork) RDKit::ScaffoldNetwork::createScaffoldNetwork<std::vector<boost::shared_ptr<RDKit::ROMol>>>;

