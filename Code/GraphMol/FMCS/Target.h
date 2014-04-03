#pragma once
#include <vector>
#include <stdexcept>
#include "../RDKitBase.h"
#include "Graph.h"
#include "MatchTable.h"
#include "DebugTrace.h" // algorithm filter definitions

namespace RDKit
{
 namespace FMCS
 {
#ifdef FAST_INCREMENTAL_MATCH
        struct AtomAdjacency
        {
            const Bond* Bond;
            unsigned    BondIdx;
            unsigned    ConnectedAtomIdx;   // another end on the bond
        };
        typedef std::vector<AtomAdjacency>  AtomAdjacencyList;
#endif // FAST_INCREMENTAL_MATCH

    struct Target
    {
        const ROMol*    Molecule;
        Graph           Topology;
#ifdef PRECOMPUTED_TABLES_MATCH
        MatchTable      AtomMatchTable;
        MatchTable      BondMatchTable;
#endif
#ifdef FAST_INCREMENTAL_MATCH
        std::vector<AtomAdjacencyList> AtomAdjacency;
#endif
    };

 }}