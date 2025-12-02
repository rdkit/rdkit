//
//  Copyright (c) 2024, Glysade Inc
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// #include "node.h"
#include "utils.h"
#include "bracket.h"

namespace RDKit {
namespace ChemDraw {
// This is currently unimplemented waiting on full bracket support in the rdkit
//  or support for expansion inside the RDChemDrawLib
bool parseBracket(CDXBracketedGroup &bracket, PageData & /*pagedata*/) {
  // Get the contained atoms/bonds in the bracket
  for (auto &attachment : bracket.ContainedObjects()) {
    auto childid = (CDXDatumID)attachment.second->GetTag();
    if (childid == kCDXObj_BracketAttachment) {
      auto &bracketattachment = (CDXBracketAttachment &)(*attachment.second);
      for (auto &bracketdata : bracketattachment.ContainedObjects()) {
        auto bracketid = (CDXDatumID)bracketdata.second->GetTag();
        if (bracketid == kCDXObj_CrossingBond) {
          // CDXCrossingBond &crossingbond =
          //     (CDXCrossingBond &)(*attachment.second);
          //  XX unimplmented crossingbond.m_bondID;       // bond that crosses
          //  brackets XX unimplmented crossingbond.m_innerAtomID;  // atom
          //  within brackets
        }
      }
    }
  }

  // SubstanceGroup sgroup;
  switch (bracket.m_usage) {
    case kCDXBracketUsage_Unspecified:
      break;
    case kCDXBracketUsage_Anypolymer:
      break;
    case kCDXBracketUsage_Component:
      break;
    case kCDXBracketUsage_Copolymer:
      break;
    case kCDXBracketUsage_CopolymerAlternating:
      break;
    case kCDXBracketUsage_CopolymerBlock:
      break;
    case kCDXBracketUsage_CopolymerRandom:
      break;
    case kCDXBracketUsage_Crosslink:
      break;
    case kCDXBracketUsage_Generic:
      break;
    case kCDXBracketUsage_Graft:
      break;
    case kCDXBracketUsage_Mer:
    case kCDXBracketUsage_MixtureOrdered:
      break;
    case kCDXBracketUsage_MixtureUnordered:
      break;
    case kCDXBracketUsage_Modification:
      break;
    case kCDXBracketUsage_Monomer:  // repeat head-to-tail, head-to-head (check
                                    // flip)
      break;
    case kCDXBracketUsage_MultipleGroup:
      break;
    case kCDXBracketUsage_MultipleGroupOverride:
      break;
    case kCDXBracketUsage_SRU:  // Structural repeating unit, repeat pattern
                                // head-to-tail (default) head-to-head (check
                                // flip?)
      break;
    case kCDXBracketUsage_Unused1:
      break;
    case kCDXBracketUsage_Unused2:
      break;
  }
  return true;
}
}  // namespace ChemDraw
}  // namespace RDKit
