//
//  Copyright (c) 2024 Glysade Inc and other RDkit contributors
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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <CDXMLParser.h>
#include <CDXStdObjects.h>

#include "chemdraw.h"
#include <RDGeneral/BadFileException.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/Atropisomers.h>

// #define DEBUG 1
namespace {
using namespace RDKit;

const std::string NEEDS_FUSE("CDXML_NEEDS_FUSE");
const std::string CDXML_FRAG_ID("CDXML_FRAG_ID");
const std::string CDXML_GROUP_ID("CDXML_GROUP_ID");
const std::string FUSE_LABEL("CDXML_NODE_ID");
const std::string CDX_SCHEME_ID("CDX_SCHEME_ID");
const std::string CDX_STEP_ID("CDX_STEP_ID");
const std::string CDX_REAGENT_ID("CDX_REAGENT_ID");
const std::string CDX_PRODUCT_ID("CDX_PRODUCT_ID");
const std::string CDX_AGENT_ID("CDX_AGENT_ID");
const std::string CDX_ATOM_POS("CDX_ATOM_POS");
const std::string CDX_ATOM_ID("_CDX_ATOM_ID");
const std::string CDX_BOND_ID("_CDX_BOND_ID");
const std::string CDX_BOND_ORDERING("CDXML_BOND_ORDERING");

constexpr double RDKIT_DEPICT_BONDLENGTH = 1.5;
struct BondInfo {
  int bond_id = -1;
  int start = -1;
  int end = -1;
  Bond::BondType order;
  std::string display;
  Bond::BondType getBondType() { return order; }
  bool validate(const std::map<unsigned int, Atom *> &ids,
                unsigned int num_atoms) const {
    auto s = ids.find(start);
    auto e = ids.find(end);
    if (s == ids.end() || e == ids.end()) {
      return false;
    }
    auto st = s->second->getIdx();
    auto et = e->second->getIdx();

    return (st < num_atoms && et < num_atoms && st != et);
  }
};

std::string NodeType(CDXNodeType nodetype) {
  switch (nodetype) {
    case kCDXNodeType_Unspecified: return "Unspecified";
    case kCDXNodeType_Element: return "Element";
    case kCDXNodeType_ElementList: return "ElementList";
    case kCDXNodeType_ElementListNickname: return "ElementListNickname";
    case kCDXNodeType_Nickname: return "Nickname";
    case kCDXNodeType_Fragment: return "Fragment";
    case kCDXNodeType_Formula: return "Forumla";
    case kCDXNodeType_GenericNickname: return "GenericNickname";
    case kCDXNodeType_AnonymousAlternativeGroup: return "Anonymous Alternative Group";
    case kCDXNodeType_NamedAlternativeGroup: return "Named Alternative Group";
    case kCDXNodeType_MultiAttachment: return "MultiAttachment";
    case kCDXNodeType_VariableAttachment: return "Variable Attachment";
    case kCDXNodeType_ExternalConnectionPoint: return "ExternalConnectionPoint";
    case kCDXNodeType_LinkNode: return "LinkNode";
    case kCDXNodeType_Monomer: return "Monomer";
      default:
      return "?";
      }
  }
  
struct StereoGroupInfo {
  int sgroup = -1;
  bool conflictingSgroupTypes = false;
  StereoGroupType grouptype;
  std::vector<Atom *> atoms;
};

struct SchemeInfo {
  int scheme_id;
  int step_id;
  std::vector<int> ReactionStepProducts;
  std::vector<int> ReactionStepReactants;
  std::vector<int> ReactionStepObjectsAboveArrow;
  std::vector<int> ReactionStepObjectsBelowArrow;
  std::vector<std::pair<int,int>> ReactionStepAtomMap;
};

unsigned int get_fuse_label(Atom *atm) {
  // return atm->getAtomMapNum(); easier debugging
  unsigned int label = 0;  // default is no label
  atm->getPropIfPresent<unsigned int>(FUSE_LABEL, label);
  return label;
}

void set_fuse_label(Atom *atm, unsigned int idx) {
  // atm->setAtomMapNum(idx); //for debugging
  if (idx) {
    atm->setProp<unsigned int>(FUSE_LABEL, idx);
  } else {
    atm->clearProp(FUSE_LABEL);
  }
}

void scaleBonds(const ROMol &mol, Conformer &conf, double targetBondLength,
                double bondLength) {
  double avg_bond_length = 0.0;
  if (bondLength < 0) {
    // If we don't have a bond length for any reason, just scale the avgerage
    // bond length
    for (auto &bond : mol.bonds()) {
      avg_bond_length += (conf.getAtomPos(bond->getBeginAtomIdx()) -
                          conf.getAtomPos(bond->getEndAtomIdx()))
                             .length();
    }
    avg_bond_length /= mol.getNumBonds();
  } else {
    avg_bond_length = bondLength;
  }

  if (avg_bond_length > 0) {
    double scale = targetBondLength / avg_bond_length;
    for (auto &pos : conf.getPositions()) {
      pos *= scale;
    }
  }
}

template <class T>
std::vector<T> to_vec(const std::string &s) {
  std::vector<T> n;
  std::stringstream ss(s);
  std::copy(std::istream_iterator<T>(ss), std::istream_iterator<T>(),
            std::back_inserter(n));
  return n;
}

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

bool parse_fragment(RWMol &mol, CDXFragment &fragment,
                    std::map<unsigned int, Atom *> &ids, int &missing_frag_id,
                    int external_attachment = -1) {
  int frag_id = fragment.GetObjectID();
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
  int atom_id = -1;
  std::vector<BondInfo> bonds;
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
        CDXNode *node = (CDXNode *)(child.second);
        int atom_id = node->GetObjectID();
        int elemno = node->m_elementNum;  // default to carbon
        // UINT16 max is not addigned?
        int num_hydrogens = node->m_numHydrogens == kNumHydrogenUnspecified
                                ? 0
                                : node->m_numHydrogens;
        bool explicitHs = node->m_numHydrogens != kNumHydrogenUnspecified;
        int charge = node->m_charge;
        int atommap = 0;
        int mergeparent = -1;
        int rgroup_num = -1;
        int isotope = node->m_isotope;
        std::string query_label;
        std::vector<int16_t> elementlist;

        // position node->m_2dPosition;
#ifdef DEBUG
        std::cerr << NodeType(node->m_nodeType) << std::endl;
#endif
        switch (node->m_nodeType) {
          case kCDXNodeType_Element: {
            break;
          }
          case kCDXNodeType_ElementList: {
            if (node->m_elementList) {
              elementlist = *node->m_elementList;
              query_label = "ElementList";
            }
            break;
          }
          case kCDXNodeType_Nickname:
          case kCDXNodeType_Fragment: {
            elemno = 0;
            atommap = atom_id;
            break;
          }
          case kCDXNodeType_ExternalConnectionPoint: {
            if (external_attachment <= 0) {
                           BOOST_LOG(rdErrorLog)
                 << "External Connection Point is not set skipping fragment";
             skip_fragment = true;
             break;
           }
           elemno = 0;
           atommap = external_attachment;
           mergeparent = external_attachment;
            break;
          }
          case kCDXNodeType_GenericNickname: {
            if (node->m_genericNickname.size()) {
              switch(node->m_genericNickname[0]) {
                  case 'R': {
                    for(auto &child : node->ContainedObjects()) {
                      if(child.second->GetTag() == kCDXObj_Text) {
                        const std::string & text = ((CDXText*)child.second)->GetText().str();
                        
                        //std::string &legacyText = (CDXText*)(child.second)->m_legacyText;
                        rgroup_num = text.size() > 1 ? stoi(text.substr(1)) : 0;
                      }
                   }
                  }
                  case 'A':
                  case 'Q': {
                    elemno = 0;
                    query_label = node->m_genericNickname;
                  }
                  break;
                  default:
                     std::cerr << "Unhandled generic nickname: " << node->m_genericNickname << std::endl;
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

        StereoGroupType grouptype = StereoGroupType::STEREO_ABSOLUTE;
        switch (node->m_enhancedStereoType) {
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
        if (mergeparent > 0) {
          rd_atom->setProp<int>("MergeParent", mergeparent);
        }

        std::vector<double> atom_coords;
        atom_coords.reserve(2);
        atom_coords.push_back(node->m_2dPosition.x);
        atom_coords.push_back(node->m_2dPosition.y);

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
          } else if (query_label == "ElementList") {
            if (!elementlist.size()) {
              BOOST_LOG(rdWarningLog)
                  << "ElementList is empty, ignoring..." << std::endl;
            } else {
              auto *q = new ATOM_OR_QUERY;
              q->setDescription("AtomOr");
              for (auto atNum : elementlist) {
                q->addChild(QueryAtom::QUERYATOM_QUERY::CHILD_TYPE(
                    makeAtomNumQuery(atNum)));
              }
              rd_atom = addquery(q, query_label, mol, idx);
              rd_atom->setAtomicNum(elementlist.front());
            }
          } else {
            rd_atom->setProp(common_properties::atomLabel, query_label);
          }
        }

        if (node->m_enhancedStereoGroupNum > 0) {
          auto key = std::make_pair(node->m_enhancedStereoGroupNum, grouptype);
          auto &stereo = sgroups[key];
          stereo.sgroup = node->m_enhancedStereoGroupNum;
          stereo.grouptype = grouptype;
          stereo.atoms.push_back(rd_atom);
        }

        ids[atom_id] = rd_atom;  // The mol has ownership so this can't leak
        if (node->m_nodeType == kCDXNodeType_Nickname ||
            node->m_nodeType == kCDXNodeType_Fragment) {
          for (auto fragment : node->ContainedObjects()) {
            if (fragment.second->GetTag() == kCDXObj_Fragment) {
              if (!parse_fragment(mol, (CDXFragment &)(*fragment.second), ids,
                                  missing_frag_id, atom_id)) {
                skip_fragment = true;
                break;
              }
              mol.setProp<bool>(NEEDS_FUSE, true);
              // might need to reset to OUR frag_id since parse_fragment will
              // set
              //  it to the fragments
              mol.setProp(CDXML_FRAG_ID, frag_id);
            }
          }
        }
        break;
      }

      case kCDXObj_Bond: {
        CDXBond &bond = *(CDXBond *)(child.second);
        int bond_id = bond.GetObjectID();
        Atom *start_atom = ids[bond.m_beginNodeID];
        Atom *end_atom = ids[bond.m_endNodeID];
        if ( (start_atom && !end_atom) ||
             (end_atom && !start_atom)) {
          BOOST_LOG(rdErrorLog) << "Bad bond in CDXML skipping fragment "
                                << frag_id << "..." << std::endl;
          skip_fragment = true;
          // ERROR
          break;
        }
        Bond::BondType order = Bond::UNSPECIFIED;
        std::unique_ptr<QueryBond> qb;
        switch (bond.m_bondOrder) {
          case kCDXBondOrder_Single:
            order = Bond::BondType::SINGLE;
            break;
          case kCDXBondOrder_Double:
            order = Bond::BondType::DOUBLE;
            break;
          case kCDXBondOrder_Triple:
            order = Bond::BondType::TRIPLE;
            break;
          case kCDXBondOrder_Quadruple:
            order = Bond::BondType::QUADRUPLE;
            break;
          case kCDXBondOrder_Quintuple:
            order = Bond::BondType::QUINTUPLE;
            break;
          case kCDXBondOrder_Sextuple:
            order = Bond::BondType::HEXTUPLE;
            break;
          case kCDXBondOrder_OneHalf:
            order = Bond::BondType::AROMATIC;
            break;
          case kCDXBondOrder_TwoHalf:
            order = Bond::BondType::TWOANDAHALF;
            break;
          case kCDXBondOrder_ThreeHalf:
            order = Bond::BondType::THREEANDAHALF;
            break;
          case kCDXBondOrder_FourHalf:
            order = Bond::BondType::FOURANDAHALF;
            break;
          case kCDXBondOrder_FiveHalf:
            order = Bond::BondType::FIVEANDAHALF;
            break;
          case kCDXBondOrder_Dative:
            order = Bond::BondType::DATIVE;
            break;
          case kCDXBondOrder_Ionic:
            order = Bond::BondType::IONIC;
            break;
          case kCDXBondOrder_SingleOrDouble: {
            order = Bond::BondType::SINGLE;
            qb = std::make_unique<QueryBond>();
            qb->setQuery(makeSingleOrDoubleBondQuery());
            break;
          }
          case kCDXBondOrder_SingleOrAromatic: {
            order = Bond::BondType::SINGLE;
            qb = std::make_unique<QueryBond>();
            qb->setQuery(makeSingleOrAromaticBondQuery());
            break;
          }
          case kCDXBondOrder_DoubleOrAromatic: {
            order = Bond::BondType::DOUBLE;
            qb = std::make_unique<QueryBond>();
            qb->setQuery(makeDoubleOrAromaticBondQuery());
            break;
          }
          case kCDXBondOrder_Any: {
            qb = std::make_unique<QueryBond>();
            qb->setQuery(makeBondNullQuery());
            break;
          }
          case kCDXBondOrder_Hydrogen:
            std::cerr << "Unhandled bond order Hydrogen" << std::endl;
          case kCDXBondOrder_ThreeCenter:
            std::cerr << "Unhandled bond order ThreeCenter" << std::endl;
          case kCDXBondOrder_Half:
            std::cerr << "Unhandled bond order Half" << std::endl;
            std::cerr << "Bad bond, skipping fragment" << std::endl;
            skip_fragment = true;
            break;
        };

        // The RDKit only supports one direction for wedges so
        //  normalize it
        bool swap = false;
        switch (bond.m_display) {
          case kCDXBondDisplay_Solid:
            break;
          case kCDXBondDisplay_Dash:
            break;
          case kCDXBondDisplay_Hash:
            break;
          case kCDXBondDisplay_WedgedHashBegin:
            break;
          case kCDXBondDisplay_WedgedHashEnd:
            swap = true;
            break;
          case kCDXBondDisplay_Bold:
            break;
          case kCDXBondDisplay_WedgeBegin:
            break;
          case kCDXBondDisplay_WedgeEnd:
            swap = true;
            break;
          case kCDXBondDisplay_Wavy:
            break;
          case kCDXBondDisplay_HollowWedgeBegin:
            break;
          case kCDXBondDisplay_HollowWedgeEnd:
            break;
          case kCDXBondDisplay_WavyWedgeBegin:
            break;
          case kCDXBondDisplay_WavyWedgeEnd:
            break;
          case kCDXBondDisplay_Dot:
            break;
          case kCDXBondDisplay_DashDot:
            break;
          case kCDXBondDisplay_DottedHydrogen:
            break;
        }

        unsigned int bondIdx = 0;
        auto startIdx = start_atom->getIdx();
        auto endIdx = end_atom->getIdx();
        if (swap) std::swap(startIdx, endIdx);

        if (qb) {
          qb->setBeginAtomIdx(startIdx);
          qb->setEndAtomIdx(endIdx);
          bondIdx = mol.addBond(qb.release(), true) - 1;
        } else {
          bondIdx = mol.addBond(startIdx, endIdx, order) - 1;
        }

        Bond *bnd = mol.getBondWithIdx(bondIdx);
        if (order == Bond::BondType::AROMATIC) {
          bnd->setIsAromatic(true);
          bnd->getBeginAtom()->setIsAromatic(true);
          bnd->getEndAtom()->setIsAromatic(true);
        }
        bnd->setProp("CDX_BOND_ID", bondIdx);

        switch (bond.m_display) {
          case kCDXBondDisplay_WedgedHashBegin:
          case kCDXBondDisplay_WedgedHashEnd: {
            bnd->setBondDir(Bond::BondDir::BEGINDASH);
            bnd->setProp(common_properties::_MolFileBondCfg, 3);
          } break;
          case kCDXBondDisplay_WedgeBegin:
          case kCDXBondDisplay_WedgeEnd: {
            bnd->setBondDir(Bond::BondDir::BEGINWEDGE);
            bnd->setProp(common_properties::_MolFileBondCfg, 1);
          } break;
          case kCDXBondDisplay_Wavy: {
            switch (order) {
              case Bond::BondType::SINGLE:
                bnd->setBondDir(Bond::BondDir::UNKNOWN);
                bnd->setProp(common_properties::_MolFileBondCfg, 2);
                break;
              case Bond::BondType::DOUBLE:
                bnd->setBondDir(Bond::BondDir::EITHERDOUBLE);
                bnd->setStereo(Bond::STEREOANY);
                break;
              default:
                BOOST_LOG(rdWarningLog)
                    << "ignoring Wavy bond set on a non double bond id: "
                    << bond_id << std::endl;
            }
            break;
          
          default:
            break;
          }
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
      stereo_groups.emplace_back(sgroup.second.grouptype,
                                 sgroup.second.atoms, newBonds, gId);
    }
    mol.setStereoGroups(std::move(stereo_groups));
  }
  
  return !skip_fragment;
}

void set_reaction_data(
    std::string type, std::string prop, SchemeInfo &scheme,
    const std::vector<int> &frag_ids,
    const std::map<unsigned int, size_t> &fragments,
    std::map<unsigned int, std::vector<int>> &grouped_fragments,
    const std::vector<std::unique_ptr<RWMol>> &mols) {
  unsigned int reagent_idx = 0;
  for (auto idx : frag_ids) {
    auto iter = grouped_fragments.find(idx);
    if (iter == grouped_fragments.end()) {
      BOOST_LOG(rdWarningLog)
          << "CDXMLParser: Schema " << scheme.scheme_id << " step "
          << scheme.step_id << " " << type << " reaction fragment " << idx
          << " not found in document." << std::endl;
      continue;
    }
    for (auto reaction_fragment_id : iter->second) {
      auto fragment = fragments.find(reaction_fragment_id);
      if (fragment == fragments.end()) {
        BOOST_LOG(rdWarningLog)
            << "CDXMLParser: Schema " << scheme.scheme_id << " step "
            << scheme.step_id << " " << type << " fragment " << idx
            << " not found in document." << std::endl;
        continue;
      }
      auto &mol = mols[fragment->second];
      mol->setProp(CDX_SCHEME_ID, scheme.scheme_id);
      mol->setProp(CDX_STEP_ID, scheme.step_id);
      mol->setProp(prop, reagent_idx);
    }
    reagent_idx += 1;
  }
}


// The parsing of fragments needed to be moved to a recursive function since
// they may be
//  embedded further in the document, i.e. a group may hold multiple
//  fragments
//
// Additionally, a grouped_fragments map is included to group fragments together
// for the purposes of reactions.
//
// Ungrouped fragments will end up as vectors of size 1 in the grouped_fragement
// list. The reaction schemes in the CDXML docs appear to use the fragment id
// for ungrouped fragments and the grouped id for grouped fragments, so the
// grouped_fragments holds both for ease of bookkeeping.
template <class T>
void visit_children(
    T &node, std::map<unsigned int, Atom *> &ids,
    std::vector<std::unique_ptr<RWMol>>
        &mols,  // All molecules found in the doc
    std::map<unsigned int, size_t>
        &fragment_lookup,  // fragment.id->molecule index
    std::map<unsigned int, std::vector<int>>
        &grouped_fragments,            // grouped.id -> [fragment.id]
    std::vector<SchemeInfo> &schemes,  // reaction schemes found
    int &missing_frag_id,  // if we don't have a fragment id, start at -1 and
                           // decrement
    double bondLength,  // bond length of the document for assigning coordinates
    const ChemDrawParserParams &params,  // parser parameters
    int group_id = -1) {  // current group id for this set of subnodes
  MolzipParams molzip_params;
  molzip_params.label = MolzipLabel::AtomProperty;
  molzip_params.atomProperty = FUSE_LABEL;
  molzip_params.enforceValenceRules = false;
  
  for (auto frag : node.ContainedObjects()) {
    CDXDatumID id = (CDXDatumID)frag.second->GetTag();
    if (id == kCDXObj_Fragment ) {
      std::unique_ptr<RWMol> mol = std::make_unique<RWMol>();
      if (!parse_fragment(*mol, (CDXFragment&)(*frag.second), ids, missing_frag_id)) {
        continue;
      }
      unsigned int frag_id = mol->getProp<int>(CDXML_FRAG_ID);
      fragment_lookup[frag_id] = mols.size();
      if (group_id != -1) {
        grouped_fragments[group_id].push_back(frag_id);
      } else {
        grouped_fragments[frag_id].push_back(frag_id);
      }
      if (mol->hasProp(NEEDS_FUSE)) {
        mol->clearProp(NEEDS_FUSE);
        std::unique_ptr<ROMol> fused;
        try {
          fused = molzip(*mol, molzip_params);
        } catch (Invar::Invariant &) {
          BOOST_LOG(rdWarningLog) << "Failed fusion of fragment skipping... "
                                  << frag_id << std::endl;
          // perhaps have an option to extract all fragments?
          // mols.push_back(std::move(mol));
          continue;
        }
        fused->setProp<int>(CDXML_FRAG_ID, static_cast<int>(frag_id));
        mols.emplace_back(dynamic_cast<RWMol *>(fused.release()));
      } else {
        mols.push_back(std::move(mol));
      }
      RWMol *res = mols.back().get();
      auto conf = std::make_unique<Conformer>(res->getNumAtoms());
      conf->set3D(false);

      bool hasConf = false;
      for (auto &atm : res->atoms()) {
        RDGeom::Point3D p{0.0, 0.0, 0.0};

        if (atm->hasProp(CDX_ATOM_POS)) {
          hasConf = true;
          const std::vector<double> coord =
              atm->getProp<std::vector<double>>(CDX_ATOM_POS);

          if (coord.size() == 2) {
            p.x = coord[0];
            p.y = -1 * coord[1];  // CDXML uses an inverted coordinate
            // system, so we need to reverse that
            p.z = 0.0;
          }
        }
        conf->setAtomPos(atm->getIdx(), p);
        atm->clearProp(CDX_ATOM_POS);
      }

      if (hasConf) {
        scaleBonds(*res, *conf, RDKIT_DEPICT_BONDLENGTH, bondLength);
        auto confidx = res->addConformer(conf.release());
        MolOps::assignChiralTypesFromBondDirs(*res, confidx, true);
        Atropisomers::detectAtropisomerChirality(*res,
                                                 &res->getConformer(confidx));
      } else {  // no Conformer
        Atropisomers::detectAtropisomerChirality(*res, nullptr);
      }

      // now that atom stereochem has been perceived, the wedging
      // information is no longer needed, so we clear
      // single bond dir flags:
      MolOps::clearSingleBondDirFlags(*res);

      if (params.sanitize) {
        try {
          if (params.removeHs) {
            // Bond stereo detection must happen before H removal, or
            // else we might be removing stereogenic H atoms in double
            // bonds (e.g. imines). But before we run stereo detection,
            // we need to run mol cleanup so don't have trouble with
            // e.g. nitro groups. Sadly, this a;; means we will find
            // run both cleanup and ring finding twice (a fast find
            // rings in bond stereo detection, and another in
            // sanitization's SSSR symmetrization).
            unsigned int failedOp = 0;
            MolOps::sanitizeMol(*res, failedOp, MolOps::SANITIZE_CLEANUP);
            MolOps::detectBondStereochemistry(*res);
            MolOps::removeHs(*res, false, false);
          } else {
            MolOps::sanitizeMol(*res);
            MolOps::detectBondStereochemistry(*res);
          }
        } catch (...) {
          BOOST_LOG(rdWarningLog)
              << "CDXMLParser: failed sanitizing skipping fragment " << frag_id
              << std::endl;
          mols.pop_back();
          continue;
        }
        MolOps::assignStereochemistry(*res, true, true, true);
      } else {
        MolOps::detectBondStereochemistry(*res);
      }
    } else if (id == kCDXObj_ReactionScheme) {  // get the reaction info
      CDXReactionScheme &scheme = (CDXReactionScheme&)(*frag.second);
      int scheme_id = scheme.GetObjectID();   //frag.second.template get<int>("<xmlattr>.id", -1);
      for (auto &rxnNode : scheme.ContainedObjects()) {
        CDXDatumID type_id = (CDXDatumID)rxnNode.second->GetTag();
        if (type_id == kCDXObj_ReactionStep) {
          CDXReactionStep &step = (CDXReactionStep&)(*rxnNode.second);
          auto step_id = step.GetObjectID();
          SchemeInfo scheme;
          scheme.scheme_id = scheme_id;
          scheme.step_id = step_id;
          scheme.ReactionStepProducts = step.m_products;
          scheme.ReactionStepReactants = step.m_reactants;
          scheme.ReactionStepObjectsBelowArrow = step.m_objectsBelowArrow;
          scheme.ReactionStepAtomMap = step.m_aamap;
          schemes.push_back(scheme);
        }
      }
    } else if (id == kCDXObj_Group) {
      CDXGroup &group = (CDXGroup&)(*frag.second);
      group_id = frag.second->GetObjectID();
      visit_children(group, ids, mols, fragment_lookup, grouped_fragments,
                     schemes, missing_frag_id, bondLength, params, group_id);
    }
  }
}

std::vector<std::unique_ptr<RWMol>> MolsFromCDXMLDataStream(
    std::istream &inStream, const ChemDrawParserParams &params) {
  CDXMLParser parser;
  // populate tree structure pt
  std::string data = std::string(std::istreambuf_iterator<char>(inStream),
                                  std::istreambuf_iterator<char>());
  const bool HaveAllXml = true;
  if(XML_STATUS_OK != parser.XML_Parse(data.c_str(), static_cast<int>(data.size()), HaveAllXml)) {
    auto error = XML_GetErrorCode(parser);
    BOOST_LOG(rdErrorLog) << "Failed parsing XML with error code " << error;
    throw FileParseException("Bad Input File");
  }
  
  std::unique_ptr<CDXDocument> document = parser.ReleaseDocument();
  if (!document) {
    // error
    return std::vector<std::unique_ptr<RWMol>>();
  }
  std::map<unsigned int, Atom *> ids;  // id->Atom* in fragment (used for linkages)
  std::vector<std::unique_ptr<RWMol>> mols;  // All molecules found in the doc
  std::map<unsigned int, size_t> fragment_lookup;  // fragment.id->molecule index
  std::map<unsigned int, std::vector<int>> grouped_fragments; // grouped.id -> [fragment.id]
  std::vector<SchemeInfo> schemes;  // reaction schemes found
  auto bondLength = document->m_bondLength;
  
  int missing_frag_id = -1;
  for (auto node : document->ContainedObjects()) {
    CDXDatumID id = (CDXDatumID)node.second->GetTag();
    switch (id) {
      case kCDXObj_Page:
        visit_children(*node.second, ids, mols, fragment_lookup, grouped_fragments,
                       schemes, missing_frag_id, bondLength, params);
        break;
      default:
        break;
    }
  }
  
    // Apply schemes
    if (schemes.size()) {
      std::map<unsigned int, size_t> fragments;
      std::map<unsigned int, size_t> agents;
      std::map<unsigned int, size_t> products;
      std::map<unsigned int, Atom *> atoms;
      size_t mol_idx = 0;
      for (auto &mol : mols) {
        auto idx = mol->getProp<unsigned int>(CDXML_FRAG_ID);
        fragments[idx] = mol_idx++;
        for (auto &atom : mol->atoms()) {
          unsigned int idx = atom->getProp<unsigned int>(CDX_ATOM_ID);
          atoms[idx] = atom;
        }
      }

      for (auto &scheme : schemes) {
        // Set the molecule properties
        set_reaction_data("ReactionStepReactants", CDX_REAGENT_ID, scheme,
                          scheme.ReactionStepReactants, fragments,
                          grouped_fragments, mols);
        set_reaction_data("ReactionStepProducts", CDX_PRODUCT_ID, scheme,
                          scheme.ReactionStepProducts, fragments,
                          grouped_fragments, mols);
        auto agents = scheme.ReactionStepObjectsAboveArrow;
        agents.insert(agents.end(),
                      scheme.ReactionStepObjectsBelowArrow.begin(),
                      scheme.ReactionStepObjectsBelowArrow.end());
        set_reaction_data("ReactionStepAgents", CDX_AGENT_ID, scheme, agents,
                          fragments, grouped_fragments, mols);
        // Set the Atom Maps
        int atommap = 0;
        for(auto mapping: scheme.ReactionStepAtomMap) {
          ++atommap;
          unsigned int idx1 = mapping.first;
          unsigned int idx2 = mapping.second;
          if (atoms.find(idx1) != atoms.end()) {
            atoms[idx1]->setAtomMapNum(atommap);
          } else {
            BOOST_LOG(rdWarningLog)
                << "CDXMLParser: Schema " << scheme.scheme_id << " step "
                << scheme.step_id
                << " ReactionStepAtomMap cannot find atom with node id " << idx1
                << "skipping schema..." << std::endl;
          }
          if (atoms.find(idx2) != atoms.end()) {
            atoms[idx2]->setAtomMapNum(atommap);
          } else {
            // XXX log error
            BOOST_LOG(rdWarningLog)
                << "CDXMLParser: Schema " << scheme.scheme_id << " step "
                << scheme.step_id
                << " ReactionStepAtomMap cannot find atom with node id " << idx2
                << " skipping schema..." << std::endl;
          }
        }
      }
    }

  // what do we do with the reaction schemes here???
  return mols;
}
}  // namespace

namespace RDKit {
std::vector<std::unique_ptr<ROMol>> ChemDrawToMols(std::istream &inStream, const ChemDrawParserParams &params) {
  auto chemdrawmols = MolsFromCDXMLDataStream(inStream, params);
  std::vector<std::unique_ptr<ROMol>> mols;
  mols.reserve(chemdrawmols.size());
  for(auto &mol : chemdrawmols) {
    ROMol *m = (ROMol*)mol.release();
    mols.push_back(std::unique_ptr<ROMol>(m));
  }
  return mols;
}

std::vector<std::unique_ptr<ROMol>> ChemDrawToMols(const std::string &filename, const ChemDrawParserParams &params) {
  CDXMLParser parser;
  std::vector<std::unique_ptr<ROMol>> mols;
  
  std::fstream chemdrawfile(filename);
  if (!chemdrawfile) {
    throw BadFileException(filename + " does not exist");
    return mols;
  }
  auto chemdrawmols = MolsFromCDXMLDataStream(chemdrawfile, params);
  
  
  mols.reserve(chemdrawmols.size());
  for(auto &mol : chemdrawmols) {
    ROMol *m = (ROMol*)mol.release();
    mols.push_back(std::unique_ptr<ROMol>(m));
  }
  return mols;
}
}  // namespace RDKit
