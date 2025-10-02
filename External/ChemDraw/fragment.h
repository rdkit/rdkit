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
//
#ifndef CHEMDRAW_FRAGMENT_H
#define CHEMDRAW_FRAGMENT_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>

#include "ChemDrawStartInclude.h"
#include "chemdraw/CDXStdObjects.h"
#include "ChemDrawEndInclude.h"

#include "reaction.h"
#include "utils.h"

namespace RDKit {
namespace ChemDraw {
struct PageData {
  PageData()
      : atomIds(),
        bondIds(),
        mols(),
        fragmentLookup(),
        groupedFragments(),
        schemes() {}

  PageData(const PageData &) = delete;

  std::map<unsigned int, Atom *> atomIds;
  std::map<unsigned int, Bond *> bondIds;
  std::vector<std::unique_ptr<RWMol>> mols;  // All molecules found in the doc
  std::map<unsigned int, size_t> fragmentLookup;  // fragment.id->molecule index
  std::map<unsigned int, std::vector<int>>
      groupedFragments;               // grouped.id -> [fragment.id]
  std::vector<ReactionInfo> schemes;  // reaction schemes found

  void clearCDXProps() {
    for (auto &mol : mols) {
      for (auto atom : mol->atoms()) {
        atom->clearProp(CDX_ATOM_ID);
        atom->clearProp(CDX_BOND_ORDERING);
        atom->clearProp(CDX_CIP);
      }
      for (auto bond : mol->bonds()) {
        bond->clearProp(CDX_BOND_ID);
      }
    }
  }
};
//! Parse a CDX fragment record
//! params
//! RWMol mol : molecule to parse the fragment into
//! CDXFragment fragment : fragment to read
//! std::map<unsigned int, Atom*> ids: atom lookup, used for bonding and fusing
//! fragments int missing_frag_id: if the fragment id is missing, this is what
//! to use.  n.b. may be obsolete, everything needs an id to be valid int
//! external_attachment:: if this fragment has a external node, this it it's id,
//! otherwise -1
//!                   external node's are normally NickNames or  new Fragments
bool parseFragment(RWMol &mol, CDXFragment &fragment, PageData &pagedata,
                   int &missingFragId, int externalAttachment = -1);
}  // namespace ChemDraw
}  // namespace RDKit

#endif
