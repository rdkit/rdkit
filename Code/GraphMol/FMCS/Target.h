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
#include <stdexcept>
#include "../RDKitBase.h"
#include "Graph.h"
#include "MatchTable.h"
#include "DebugTrace.h" // algorithm filter definitions

namespace RDKit {
    namespace FMCS {
        struct Target {
            const ROMol*    Molecule;
            Graph           Topology;
            MatchTable      AtomMatchTable;
            MatchTable      BondMatchTable;
        };

    }
}
