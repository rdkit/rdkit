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
#include "fragment.h"
#include "bond.h"
#include "node.h"

namespace RDKit {
namespace {
const char * sequenceTypeToName(CDXSeqType seqtype) {
  switch(seqtype) {
      case kCDXSeqType_Unknown : return "Unknown";
      case kCDXSeqType_Peptide: return "Peptide (Helm)";      // HELM peptides
      case kCDXSeqType_Peptide1 : return "Peptide1 (Single Letter Amino Acid)";      // Single letter amino acids (Legacy biopolymer support)
      case kCDXSeqType_Peptide3 : return "Peptide3 (Three letter amino acid)";      // Three letter amino acids (Legacy biopolymer support)
      case kCDXSeqType_DNA : return "DNA";
      case kCDXSeqType_RNA : return "RNA";
      case kCDXSeqType_Biopolymer : return "Biopolymer";
      default:
        return "";
    }
}
}
bool parse_fragment(RWMol &mol, CDXFragment &fragment,
                    PageData &pagedata, int &missing_frag_id,
                    int external_attachment) {
  int frag_id = fragment.GetObjectID();
  if(fragment.m_sequenceType != kCDXSeqType_Unknown) {
    BOOST_LOG(rdWarningLog)
    << "Unhandled chemdraw sequence type " << sequenceTypeToName(fragment.m_sequenceType)
    << std::endl;
    return false;
  }
  if (frag_id == -1) {
    // ChemDraw simply assigns a new one
    BOOST_LOG(rdWarningLog)
        << "Invalid or missing fragment id from CDXML fragment, assigning new one..."
        << std::endl;
    frag_id = missing_frag_id;
    missing_frag_id--;
  }
  mol.setProp(CDXML_FRAG_ID, frag_id);

  // for atom in frag
  std::map<std::pair<int, StereoGroupType>, StereoGroupInfo> sgroups;

  // nodetypes =
  // https://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/properties/Node_Type.htm
  bool skip_fragment =
      false;  // is there an irrecoverable error for this fragment

  for (auto child : fragment.ContainedObjects()) {
    CDXDatumID id = (CDXDatumID)child.second->GetTag();
#ifdef DEBUG
    std::cerr << "Data Type: " << id << std::endl;
#endif
    switch (id) {
      case kCDXObj_Node: {
        CDXNode &node = (CDXNode &)(*child.second);
        if (!parse_node(mol, frag_id, node, pagedata, sgroups, missing_frag_id,
                        external_attachment)) {
          skip_fragment = true;
        }
        break;
      }
      case kCDXObj_Bond: {
        CDXBond &bond = (CDXBond &)(*child.second);
        if (!parse_bond(mol, frag_id, bond, pagedata)) {
          skip_fragment = true;
          break;
        }
      }
    }
  }

  // Add the stereo groups
  if (!sgroups.empty()) {
    std::vector<StereoGroup> stereo_groups;
    for (auto &sgroup : sgroups) {
      unsigned gId = 0;
      if (sgroup.second.grouptype != StereoGroupType::STEREO_ABSOLUTE &&
          sgroup.second.sgroup > 0) {
        gId = sgroup.second.sgroup;
      }
      std::vector<Bond *> newBonds;
      stereo_groups.emplace_back(sgroup.second.grouptype, sgroup.second.atoms,
                                 newBonds, gId);
    }
    mol.setStereoGroups(std::move(stereo_groups));
  }

  return !skip_fragment;
}

}  // namespace RDKit
