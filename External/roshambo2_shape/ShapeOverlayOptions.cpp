//
//  Copyright (C) 2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#include <RDGeneral/export.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "ShapeOverlayOptions.h"

namespace RDKit {
namespace ShapeAlign {

ShapeOverlayOptions::ShapeOverlayOptions() {
  buildPh4Patterns();
  buildInteractionMatrices();
}

// Definitions for feature points adapted from:
// Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)
const std::vector<std::vector<std::string>> smartsPatterns = {
    {"[$([N;!H0;v3,v4&+1]),\
$([O,S;H1;+0]),\
n&H1&+0]"},                                // Donor
    {"[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
$([O,S;H0;v2]),\
$([O,S;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),\
n&H0&+0,\
$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]"},  // Acceptor
    {
        "[r]1[r][r]1",
        "[r]1[r][r][r]1",
        "[r]1[r][r][r][r]1",
        "[r]1[r][r][r][r][r]1",
        "[r]1[r][r][r][r][r][r]1",
    },  // rings
        //    "[a]",                                                  //
        //    Aromatic
        //    "[F,Cl,Br,I]",                                          // Halogen
    {"[#7;+,\
$([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]"},  // Basic
    {"[$([C,S](=[O,S,P])-[O;H1,-1])]"}                       // Acidic
};
void ShapeOverlayOptions::buildPh4Patterns() {
  for (const auto &smartsV : smartsPatterns) {
    std::vector<std::shared_ptr<ROMol>> v;
    for (const auto &smarts : smartsV) {
      v.emplace_back(v2::SmilesParse::MolFromSmarts(smarts));
    }
    d_ph4Patterns.emplace_back(v);
  }
  d_nTypes = d_ph4Patterns.size();
}

// Make the r and p matrices for the features in the shapes.  Returned
// as linearised square matrices, of size nFeats * nFeats.  Currently both
// parameters are set to 1.0 for all types that are used.  Types are numbered
// from 1 - type 0 is ignored in the color volume calculation.
void ShapeOverlayOptions::buildInteractionMatrices() {
  d_rMat.resize((d_nTypes + 1) * (d_nTypes + 1), 1.0);
  d_pMat.resize((d_nTypes + 1) * (d_nTypes + 1), 1.0);
}

}  // namespace ShapeAlign
}  // namespace RDKit
