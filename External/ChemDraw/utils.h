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
#ifndef CHEMDRAW_UTILS_H
#define CHEMDRAW_UTILS_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>

#include "ChemDrawStartInclude.h"
#include "chemdraw/CDXStdObjects.h"
#include "ChemDrawEndInclude.h"

namespace RDKit {
constexpr double RDKIT_DEPICT_BONDLENGTH = 1.5;
const std::string NEEDS_FUSE("CDX_NEEDS_FUSE");
const std::string CDX_FRAG_ID("CDX_FRAG_ID");
const std::string CDX_GROUP_ID("CDX_GROUP_ID");
const std::string FUSE_LABEL("CDX_NODE_ID");
const std::string CDX_SCHEME_ID("CDX_SCHEME_ID");
const std::string CDX_STEP_ID("CDX_STEP_ID");
const std::string CDX_REAGENT_ID("CDX_REAGENT_ID");
const std::string CDX_PRODUCT_ID("CDX_PRODUCT_ID");
const std::string CDX_AGENT_ID("CDX_AGENT_ID");
const std::string CDX_ATOM_POS("CDX_ATOM_POS");
const std::string CDX_ATOM_ID("_CDX_ATOM_ID");
const std::string CDX_BOND_ID("_CDX_BOND_ID");
const std::string CDX_BOND_ORDERING("CDX_BOND_ORDERING");
const std::string CDX_CIP("CDX_CIP");
const std::string CDX_IMPLICIT_HYDROGEN_STEREO("CDX_ATOM_STEREO");

// Convert a ChemDrawNode to a string
std::string NodeType(CDXNodeType nodetype);

// Scale the bonds to the targetBondLength.  If bondLength is zero
//  use the average bond length in the molecule
void scaleBonds(const ROMol &mol, Conformer &conf, double targetBondLength,
                double bondLength);

// Indicate which atoms should be fused together from various
//  fragments in the ChemDraw file

unsigned int get_fuse_label(Atom *atm);
void set_fuse_label(Atom *atm, unsigned int idx);

// Replace fragments that are not possible with molzip
bool replaceFragments(RWMol &mol);

// Add a Query to a molecule
template <typename Q>
Atom *addquery(Q *qry, std::string symbol, RWMol &mol, unsigned int idx) {
  PRECONDITION(qry, "bad query");
  auto *atm = mol.getAtomWithIdx(idx);
  auto qa = std::make_unique<QueryAtom>(*atm);
  qa->setQuery(qry);
  qa->setNoImplicit(true);
  mol.replaceAtom(idx, qa.get());
  Atom *res = mol.getAtomWithIdx(idx);
  if (symbol != "") {
    res->setProp(common_properties::atomLabel, symbol);
  }
  return res;
}

// Simple Structure for keeping track of Stereo Groups
struct StereoGroupInfo {
  int sgroup = -1;
  bool conflictingSgroupTypes = false;
  StereoGroupType grouptype;
  std::vector<Atom *> atoms;
};

// check to see if we have a tetrahedral flag and ChemDraw CIP set but no
//  stereo assigned, if so check the bond ordering for CW and CCW
void checkChemDrawTetrahedralGeometries(RWMol &mol);
}  // namespace RDKit

#endif
