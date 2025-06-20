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
//#include "node.h"
#include "fragment.h"
#include "utils.h"

namespace RDKit {
namespace ChemDraw {
bool parse_node(
    RWMol &mol, unsigned int fragment_id, CDXNode &node, PageData &pagedata,
    std::map<std::pair<int, StereoGroupType>, StereoGroupInfo> &sgroups,
    int &missing_frag_id, int external_attachment) {
  int atom_id = node.GetObjectID();
  int elemno = node.m_elementNum;  // default to carbon
  // UINT16 max is not addigned?
  int num_hydrogens =
      node.m_numHydrogens == kNumHydrogenUnspecified ? 0 : node.m_numHydrogens;
  bool explicitHs = node.m_numHydrogens != kNumHydrogenUnspecified;
  int charge = 0;
  if ((node.m_charge & 0x00FFFFFF) == 0)
    charge = node.m_charge >> 24;
  else
    charge = node.m_charge;
  int atommap = 0;
  int rgroup_num = -1;
  int isotope = node.m_isotope;

  bool checkForRGroup = false;
  ;
  std::string query_label;
  std::vector<int16_t> elementlist;

  // position node.m_2dPosition;
#ifdef DEBUG
  std::cerr << NodeType(node.m_nodeType) << std::endl;
#endif
  switch (node.m_nodeType) {
    case kCDXNodeType_Element: {
      break;
    }
    case kCDXNodeType_ElementList: {
      if (node.m_elementList) {
        elementlist = *node.m_elementList;
        query_label = "ElementList";
      }
      break;
    }
    case kCDXNodeType_Nickname: {
      elemno = 0;
      atommap = atom_id;
      break;
    }
    case kCDXNodeType_Fragment: {
      elemno = 0;
      atommap = atom_id;
      break;
    }
    case kCDXNodeType_ExternalConnectionPoint: {
      if (external_attachment <= 0) {
        // sometimes this is a dummy atom, but I don't know when.
        if (node.m_externalConnectionType == kCDXExternalConnection_Diamond) {
          elemno = 0;
        }
        atommap = atom_id;
      } else {
        elemno = 0;
        atommap = external_attachment;
      }
      break;
    }
    case kCDXNodeType_GenericNickname: {
      if (node.m_genericNickname.size()) {
        switch (node.m_genericNickname[0]) {
          case 'R': {
            checkForRGroup = true;
            elemno = 0;
            query_label = node.m_genericNickname;
            break;
          }
          case 'A':
          case 'Q':
          case 'X':
          case 'M': {
            elemno = 0;
            query_label = node.m_genericNickname;
          } break;
          default:
            std::cerr << "Unhandled generic nickname: "
                      << node.m_genericNickname << std::endl;
        }
      }
      break;
    }
    case kCDXNodeType_Unspecified:
      break;
    case kCDXNodeType_ElementListNickname:
      break;
    case kCDXNodeType_Formula:
      break;
    case kCDXNodeType_AnonymousAlternativeGroup:
      break;
    case kCDXNodeType_NamedAlternativeGroup:
      break;
    case kCDXNodeType_MultiAttachment:
      break;
    case kCDXNodeType_VariableAttachment:
      break;
    case kCDXNodeType_LinkNode:
      break;
    case kCDXNodeType_Monomer:
      break;
  }

  for (auto &child : node.ContainedObjects()) {
    if (child.second->GetTag() == kCDXObj_Text) {
      const std::string &text = ((CDXText *)child.second)->GetText().str();
      if (text.size() > 0 && text[0] == 'R') {
        try {
          if (checkForRGroup)
            rgroup_num = text.size() > 1 ? stoi(text.substr(1)) : 0;
          else
            isotope = text.size() > 1 ? stoi(text.substr(1)) : 0;
        } catch (const std::invalid_argument &e) {
          if (rgroup_num)
            BOOST_LOG(rdWarningLog)
                << "RGroupError: Invalid argument - Cannot convert '" << text
                << "' to an integer." << std::endl;
        } catch (const std::out_of_range &e) {
          if (rgroup_num)
            BOOST_LOG(rdWarningLog)
                << "RGroupError: Out of range - The number '" << text
                << "' is too large or too small." << std::endl;
        }
      }
    }
  }

  StereoGroupType grouptype = StereoGroupType::STEREO_ABSOLUTE;
  switch (node.m_enhancedStereoType) {
    case kCDXEnhancedStereo_Absolute:
      grouptype = StereoGroupType::STEREO_ABSOLUTE;
      break;
    case kCDXEnhancedStereo_And:
      grouptype = StereoGroupType::STEREO_AND;
      break;
    case kCDXEnhancedStereo_Or:
      grouptype = StereoGroupType::STEREO_OR;
      break;
    default:
      break;
  }

  CHECK_INVARIANT(atom_id != -1, "Uninitialized atom id in cdxml.");
  Atom *rd_atom = new Atom(elemno);
  rd_atom->setFormalCharge(charge);
  rd_atom->setNumExplicitHs(num_hydrogens);
  rd_atom->setNoImplicit(explicitHs);

  rd_atom->setIsotope(isotope);
  if (rgroup_num >= 0) {
    rd_atom->setAtomMapNum(rgroup_num);
  }
  set_fuse_label(rd_atom, atommap);
  switch (node.m_hStereo) {
    case kCDXProp_Atom_HDot:  // this atom has an implicit hydrogen with a
                              // wedged bond
      rd_atom->setProp<char>(CDX_IMPLICIT_HYDROGEN_STEREO, 'w');
      break;
    case kCDXProp_Atom_HDash:  // this atom has an implicit hydrogen with a
                               // hashed bond
      rd_atom->setProp<char>(CDX_IMPLICIT_HYDROGEN_STEREO, 'h');
      break;
  }
  
  if (node.m_bondOrdering) {
    // This node may be completely replaced by the fragment
    // i.e. [*:1]C[*:1].C[*:1]C => CCC
    rd_atom->setProp<std::vector<int>>(CDX_BOND_ORDERING, *node.m_bondOrdering);
  }
  if (node.m_geometry == kCDXAtomGeometry_Tetrahedral) {
    // std::cerr << "tetrahedral" << std::endl;
    //  if we have a cip type we can interpret, set it, otherwise don't

    switch (node.m_CIP) {
      case kCDXCIPAtom_R:
      case kCDXCIPAtom_r:
      case kCDXCIPAtom_S:
      case kCDXCIPAtom_s:
        rd_atom->setProp<CDXAtomCIPType>(CDX_CIP, node.m_CIP);
        break;
      default:
        rd_atom->setProp<CDXAtomCIPType>(CDX_CIP, kCDXCIPAtom_Undetermined);
        break;
    }
  }

  std::vector<double> atom_coords;
  if (node.KnownPosition3D()) {
    atom_coords.reserve(3);
    atom_coords.push_back(node.m_3dPosition.x);
    atom_coords.push_back(node.m_3dPosition.y);
    atom_coords.push_back(node.m_3dPosition.z);
  } else {
    atom_coords.reserve(2);
    atom_coords.push_back(node.m_2dPosition.x);
    atom_coords.push_back(node.m_2dPosition.y);
  }
  rd_atom->setProp<std::vector<double>>(CDX_ATOM_POS, atom_coords);
  rd_atom->setProp<unsigned int>(CDX_ATOM_ID, atom_id);

  const bool updateLabels = true;
  const bool takeOwnership = true;
  auto idx = mol.addAtom(rd_atom, updateLabels, takeOwnership);
  if (query_label.size()) {
    if (query_label[0] == 'R') {
      rd_atom = addquery(makeAtomNullQuery(), query_label, mol, idx);
    } else if (query_label == "A") {
      rd_atom = addquery(makeAAtomQuery(), query_label, mol, idx);
    } else if (query_label == "Q") {
      rd_atom = addquery(makeQAtomQuery(), query_label, mol, idx);
    } else if (query_label == "M") {
      rd_atom = addquery(makeMAtomQuery(), query_label, mol, idx);
    } else if (query_label == "MH") {
      rd_atom = addquery(makeMHAtomQuery(), query_label, mol, idx);
    } else if (query_label == "X") {
      rd_atom = addquery(makeXAtomQuery(), query_label, mol, idx);
    } else if (query_label == "ElementList") {
      if (!elementlist.size()) {
        BOOST_LOG(rdWarningLog)
            << "ElementList is empty, ignoring..." << std::endl;
      } else {
        auto *q = new ATOM_OR_QUERY;
        q->setDescription("AtomOr");
        for (auto atNum : elementlist) {
          q->addChild(
              QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(makeAtomNumQuery(atNum)));
        }
        rd_atom = addquery(q, query_label, mol, idx);
        rd_atom->setAtomicNum(elementlist.front());
      }
    } else if (query_label.size()) {
      std::cerr << "Unhandled generic nickname: " << query_label << std::endl;
    } else {
      rd_atom->setProp(common_properties::atomLabel, query_label);
    }
  }

  switch (node.m_radical) {
    case kCDXRadical_None:
      break;
    case kCDXRadical_Singlet:
      rd_atom->setNumRadicalElectrons(2);
      break;
    case kCDXRadical_Doublet: {
      rd_atom->setNumRadicalElectrons(1);
      break;
    }
    case kCDXRadical_Triplet: {
      rd_atom->setNumRadicalElectrons(2);
      break;
    }
  }

  if (node.m_enhancedStereoGroupNum > 0) {
    auto key = std::make_pair(node.m_enhancedStereoGroupNum, grouptype);
    auto &stereo = sgroups[key];
    stereo.sgroup = node.m_enhancedStereoGroupNum;
    stereo.grouptype = grouptype;
    stereo.atoms.push_back(rd_atom);
  }

  pagedata.atom_ids[atom_id] =
      rd_atom;  // The mol has ownership so this can't leak
  if (node.m_nodeType == kCDXNodeType_Nickname ||
      node.m_nodeType == kCDXNodeType_Fragment) {
    // This fragment needs to be expanded and joined to the current one
    //  the external_id is the node's atom_id
    for (auto fragment : node.ContainedObjects()) {
      if (fragment.second->GetTag() == kCDXObj_Fragment) {
        if (!parse_fragment(mol, (CDXFragment &)(*fragment.second), pagedata,
                            missing_frag_id, atom_id)) {
          return false;
        }
        mol.setProp<bool>(NEEDS_FUSE, true);
        // might need to reset to OUR frag_id since parse_fragment will
        // set
        //  it to the fragments
        mol.setProp(CDX_FRAG_ID, fragment_id);
      }
    }
  }
  return true;
}
}
}  // namespace RDKit
