//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include <vector>
#include "FMCS.h"
#include "Graph.h"
#include "MatchTable.h"

namespace RDKit {
namespace FMCS {
typedef std::vector<std::pair<FMCS::Graph::vertex_descriptor,
                              FMCS::Graph::vertex_descriptor> > match_V_t;
const unsigned int NotSet = (unsigned int)-1;

bool SubstructMatchCustomTable(
    const FMCS::Graph& target, const ROMol& target_mol,
    const FMCS::Graph& query,
    const ROMol& querySrc  // seed and full source query molecules
    ,
    const MatchTable& atomMatchTable, const MatchTable& bondMatchTable,
    const MCSParameters* parameters = 0  // for final checker (CHIRALITY)
    ,
    match_V_t* match = 0);

bool SubstructMatchCustom(
    const FMCS::Graph& target, const ROMol& mol, const FMCS::Graph& query,
    const ROMol& querySrc  // seed and full source query molecules
    ,
    MCSAtomCompareFunction atomCompare, MCSBondCompareFunction bondCompare,
    MCSFinalMatchCheckFunction finalCompare,
    const MCSAtomCompareParameters& acp, const MCSBondCompareParameters& bcp,
    void* user_data, match_V_t* match = 0);
}
}
