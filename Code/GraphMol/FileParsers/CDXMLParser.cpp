//
//  Copyright (c) 2022 Brian P Kelley
//  All rights reserved.
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "CDXMLParser.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/BadFileException.h>
#include <fstream>
#include <sstream>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/Atropisomers.h>

using boost::property_tree::ptree;
namespace RDKit {

namespace {
const std::string NEEDS_FUSE("CDXML_NEEDS_FUSE");
const std::string CDXML_FRAG_ID("CDXML_FRAG_ID");
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

struct StereoGroupInfo {
  int sgroup = -1;
  bool conflictingSgroupTypes = false;
  StereoGroupType grouptype;
  std::vector<Atom *> atoms;
};

struct SchemeInfo {
  int scheme_id;
  int step_id;
  std::vector<unsigned int> ReactionStepProducts;
  std::vector<unsigned int> ReactionStepReactants;
  std::vector<unsigned int> ReactionStepObjectsAboveArrow;
  std::vector<unsigned int> ReactionStepObjectsBelowArrow;
  std::vector<unsigned int> ReactionStepAtomMap;
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

template <class T>
std::vector<T> to_vec(const std::string &s) {
  std::vector<T> n;
  std::stringstream ss(s);
  std::copy(std::istream_iterator<T>(ss), std::istream_iterator<T>(),
            std::back_inserter(n));
  return n;
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
bool parse_fragment(RWMol &mol, ptree &frag,
                    std::map<unsigned int, Atom *> &ids, int &missing_frag_id,
                    int external_attachment = -1) {
  // XXX Need to put the fragid on the molecule so we can properly make
  // reactions

  int frag_id = frag.get<int>("<xmlattr>.id", -1);
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
  for (auto &node : frag) {
    if (node.first == "n") {  // atom node
      int elemno = 6;         // default to carbon
      int num_hydrogens = 0;
      int charge = 0;
      int atommap = 0;
      int mergeparent = -1;
      int rgroup_num = -1;
      int isotope = 0;
      int sgroup = -1;
      bool explicitHs = false;
      StereoGroupType grouptype = StereoGroupType::STEREO_ABSOLUTE;
      std::string query_label;
      std::vector<int> bond_ordering;
      std::vector<int> elementlist;
      std::vector<double> atom_coords;
      std::string nodetype = "";
      for (auto &attr : node.second.get_child("<xmlattr>")) {
        try {
          if (attr.first == "id") {
            atom_id = stoi(attr.second.data());
            if (ids.find(atom_id) != ids.end()) {
              BOOST_LOG(rdErrorLog) << "Warning, duplicated atom id " << atom_id
                                    << " skipping fragment" << std::endl;
              skip_fragment = true;
              break;
            }
          } else if (attr.first == "Element") {
            elemno = stoi(attr.second.data());
          } else if (attr.first == "NumHydrogens") {
            num_hydrogens = stoi(attr.second.data());
            explicitHs = true;
          } else if (attr.first == "Charge") {
            charge = stoi(attr.second.data());
          } else if (attr.first == "Isotope") {
            isotope = stoi(attr.second.data());
          } else if (attr.first == "NodeType") {
            nodetype = attr.second.data();
            if (nodetype == "Nickname" || nodetype == "Fragment") {
              elemno = 0;
              atommap = atom_id;
            } else if (nodetype == "ExternalConnectionPoint") {
              if (external_attachment <= 0) {
                BOOST_LOG(rdErrorLog)
                    << "External Connection Point is not set skipping fragment";
                skip_fragment = true;
                break;
              }
              elemno = 0;
              atommap = external_attachment;
              mergeparent = external_attachment;
            } else if (nodetype == "GenericNickname") {
              // RGroups for example?
              for (auto &tnode : node.second) {
                if (tnode.first == "t") {
                  for (auto &snode : tnode.second) {
                    if (snode.first == "s") {
                      auto s = snode.second.data();
                      if (s.size()) {
                        if (s[0] == 'R') {
                          if (s.size() > 1) {
                            rgroup_num = stoi(s.substr(1));
                          }
                          elemno = 0;
                          query_label = s;
                        } else if (s == "A") {
                          query_label = s;
                          elemno = 0;
                        } else if (s == "Q") {
                          query_label = s;
                          elemno = 0;
                        }
                      }
                      break;
                    }
                  }
                  break;
                }
              }
            } else if (nodetype == "ElementList") {
              query_label = "ElementList";
            }
          } else if (attr.first == "ElementList") {
            elementlist = to_vec<int>(attr.second.data());

          } else if (attr.first == "p") {
            atom_coords = to_vec<double>(attr.second.data());
          } else if (attr.first == "EnhancedStereoGroupNum") {
            sgroup = stoi(attr.second.data());
          } else if (attr.first == "EnhancedStereoType") {
            auto stereo_type = attr.second.data();
            if (stereo_type == "And") {
              grouptype = StereoGroupType::STEREO_AND;
            } else if (stereo_type == "Or") {
              grouptype = StereoGroupType::STEREO_OR;
            } else if (stereo_type == "Absolute") {
              grouptype = StereoGroupType::STEREO_ABSOLUTE;
            } else {
              BOOST_LOG(rdWarningLog)
                  << "Unhandled enhanced stereo type " << stereo_type
                  << " ignoring" << std::endl;
            }
          }
        } catch (...) {
          BOOST_LOG(rdErrorLog)
              << "Failed to parse XML fragment " << frag_id
              << " node: " << node.first << " attribute: " << attr.first << ": "
              << attr.second.data() << std::endl;
          return false;
        }
      }
      // add the atom
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
      if (sgroup != -1) {
        auto key = std::make_pair(sgroup, grouptype);
        auto &stereo = sgroups[key];
        stereo.sgroup = sgroup;
        stereo.grouptype = grouptype;
        stereo.atoms.push_back(rd_atom);
      }
      ids[atom_id] = rd_atom;  // The mol has ownership so this can't leak
      if (nodetype == "Nickname" || nodetype == "Fragment") {
        for (auto &fragment : node.second) {
          if (fragment.first == "fragment") {
            if (!parse_fragment(mol, fragment.second, ids, missing_frag_id,
                                atom_id)) {
              skip_fragment = true;
              break;
            }
            mol.setProp<bool>(NEEDS_FUSE, true);
            // might need to reset to OUR frag_id since parse_fragment will set
            //  it to the fragments
            mol.setProp(CDXML_FRAG_ID, frag_id);
          }
        }
      }
    } else if (node.first == "b") {  // bond
      int bond_id = -1;
      int start_atom = -1;
      int end_atom = -1;
      Bond::BondType order = Bond::SINGLE;
      std::string display;
      for (auto &attr : node.second.get_child("<xmlattr>")) {
        try {
          if (attr.first == "id") {
            bond_id = stoi(attr.second.data());
          } else if (attr.first == "B") {
            start_atom = stoi(attr.second.data());
          } else if (attr.first == "E") {
            end_atom = stoi(attr.second.data());
          } else if (attr.first == "Order") {
            if (attr.second.data() == "1.5") {
              order = Bond::BondType::AROMATIC;
            } else if (attr.second.data() == "any") {
              order = Bond::BondType::UNSPECIFIED;
            } else {
              int bond_order = stoi(attr.second.data());

              switch (bond_order) {
                case 1:
                  order = Bond::BondType::SINGLE;
                  break;
                case 2:
                  order = Bond::BondType::DOUBLE;
                  break;
                case 3:
                  order = Bond::BondType::TRIPLE;
                  break;
                case 4:
                  order = Bond::BondType::QUADRUPLE;
                  break;
                default:
                  throw std::invalid_argument("Unhandled bond order");
              }
            }
          } else if (attr.first ==
                     "Display") {  // gets wedge/hash stuff and probably more
            display = attr.second.data();
          }
        } catch (...) {
          BOOST_LOG(rdErrorLog)
              << "Failed to parse XML fragment " << frag_id
              << " node: " << node.first << " attribute: " << attr.first << ": "
              << attr.second.data() << std::endl;
          return false;
        }
      }
      // CHECK_INVARIANT(start_atom>=0 && end_atom>=0 && start_atom != end_atom,
      // "Bad bond in CDXML");
      BondInfo bond{bond_id, start_atom, end_atom, order, display};
      if (!bond.validate(ids, mol.getNumAtoms())) {
        BOOST_LOG(rdErrorLog) << "Bad bond in CDXML skipping fragment "
                              << frag_id << "..." << std::endl;
        ;
        skip_fragment =
            true;  // ChemDraw doesn't skip, it just ignores bad bonds...
        break;
      } else {
        bonds.push_back(bond);
      }
      // end if atom or bond
    }
  }  // for node

  // add bonds
  if (!skip_fragment) {
    for (auto &bond : bonds) {
      bool swap = false;
      if (bond.display == "WedgeEnd") {
        swap = true;
        bond.display = "WedgeBegin";
      }
      if (bond.display == "WedgedHashEnd") {
        swap = true;
        bond.display = "WedgedHashBegin";
      }

      auto startIdx = ids[bond.start]->getIdx();
      auto endIdx = ids[bond.end]->getIdx();
      if (swap) {
        std::swap(startIdx, endIdx);
      }
      unsigned bondIdx = 0;
      if (bond.order == Bond::BondType::UNSPECIFIED) {
        auto qb = new QueryBond();
        qb->setQuery(makeBondNullQuery());
        qb->setBeginAtomIdx(startIdx);
        qb->setEndAtomIdx(endIdx);
        bondIdx = mol.addBond(qb, true) - 1;
      } else {
        bondIdx = mol.addBond(startIdx, endIdx, bond.getBondType()) - 1;
      }
      Bond *bnd = mol.getBondWithIdx(bondIdx);
      if (bond.order == Bond::BondType::AROMATIC) {
        bnd->setIsAromatic(true);
        ids[bond.end]->setIsAromatic(true);
        ids[bond.start]->setIsAromatic(true);
      }
      bnd->setProp("CDX_BOND_ID", bond.bond_id);
      if (bond.display == "WedgeBegin") {
        bnd->setBondDir(Bond::BondDir::BEGINWEDGE);
        bnd->setProp(common_properties::_MolFileBondCfg, 1);
      } else if (bond.display == "WedgedHashBegin") {
        bnd->setBondDir(Bond::BondDir::BEGINDASH);
        bnd->setProp(common_properties::_MolFileBondCfg, 3);
      } else if (bond.display == "Wavy") {
        switch (bond.getBondType()) {
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
                << bond.bond_id << std::endl;
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

void set_reaction_data(std::string type, std::string prop, SchemeInfo &scheme,
                       const std::vector<unsigned int> &frag_ids,
                       const std::map<unsigned int, size_t> &fragments,
                       const std::vector<std::unique_ptr<RWMol>> &mols) {
  unsigned int reagent_idx = 0;
  for (auto idx : frag_ids) {
    auto iter = fragments.find(idx);
    if (iter == fragments.end()) {
      BOOST_LOG(rdWarningLog)
          << "CDXMLParser: Schema " << scheme.scheme_id << " step "
          << scheme.step_id << " " << type << " fragment " << idx
          << " not found in document." << std::endl;
      continue;
    }
    if (iter->second >= mols.size()) {
      // shouldn't get here
      continue;
    }
    auto &mol = mols[iter->second];
    mol->setProp(CDX_SCHEME_ID, scheme.scheme_id);
    mol->setProp(CDX_STEP_ID, scheme.step_id);
    mol->setProp(prop, reagent_idx++);
  }
}
}  // namespace

namespace v2 {
namespace CDXMLParser {

std::vector<std::unique_ptr<RWMol>> MolsFromCDXMLDataStream(
    std::istream &inStream, const CDXMLParserParams &params) {
  // populate tree structure pt
  using boost::property_tree::ptree;
  ptree pt;
  try {
    read_xml(inStream, pt);
  } catch (boost::property_tree::ptree_error &e) {
    auto xml = dynamic_cast<boost::property_tree::file_parser_error *>(&e);
    if (xml != nullptr) {
      auto msg = std::string(xml->message()) +
                 " at line: " + boost::lexical_cast<std::string>(xml->line());
      throw FileParseException(msg);
    }

    throw FileParseException(e.what());
  }
  std::map<unsigned int, Atom *> ids;
  std::vector<std::unique_ptr<RWMol>> mols;
  std::map<unsigned int, size_t> fragment_lookup;
  std::vector<SchemeInfo> schemes;

  MolzipParams molzip_params;
  molzip_params.label = MolzipLabel::AtomProperty;
  molzip_params.atomProperty = FUSE_LABEL;
  molzip_params.enforceValenceRules = false;
  int missing_frag_id = -1;
  for (auto &cdxml : pt) {
    if (cdxml.first == "CDXML") {
      double bondLength = cdxml.second.get<double>("<xmlattr>.BondLength");
      for (auto &node : cdxml.second) {
        if (node.first == "page") {
          for (auto &frag : node.second) {
            if (frag.first == "fragment") {  // chemical matter
              std::unique_ptr<RWMol> mol = std::make_unique<RWMol>();
              if (!parse_fragment(*mol, frag.second, ids, missing_frag_id)) {
                continue;
              }
              unsigned int frag_id = mol->getProp<int>(CDXML_FRAG_ID);
              fragment_lookup[frag_id] = mols.size();
              if (mol->hasProp(NEEDS_FUSE)) {
                mol->clearProp(NEEDS_FUSE);
                std::unique_ptr<ROMol> fused;
                try {
                  fused = molzip(*mol, molzip_params);
                } catch (Invar::Invariant &) {
                  BOOST_LOG(rdWarningLog)
                      << "Failed fusion of fragment skipping... " << frag_id
                      << std::endl;
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
                DetectAtomStereoChemistry(*res, &res->getConformer(confidx));

                Atropisomers::detectAtropisomerChirality(
                    *res, &res->getConformer(confidx));
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
                    MolOps::sanitizeMol(*res, failedOp,
                                        MolOps::SANITIZE_CLEANUP);
                    MolOps::detectBondStereochemistry(*res);
                    MolOps::removeHs(*res, false, false);
                  } else {
                    MolOps::sanitizeMol(*res);
                    MolOps::detectBondStereochemistry(*res);
                  }
                } catch (...) {
                  BOOST_LOG(rdWarningLog)
                      << "CDXMLParser: failed sanitizing skipping fragment "
                      << frag_id << std::endl;
                  mols.pop_back();
                  continue;
                }
                MolOps::assignStereochemistry(*res, true, true, true);
              } else {
                MolOps::detectBondStereochemistry(*res);
              }
            } else if (frag.first == "scheme") {  // get the reaction info
              int scheme_id = frag.second.get<int>("<xmlattr>.id", -1);
              for (auto &node : frag.second) {
                if (node.first == "step") {
                  auto step_id = node.second.get<int>("<xmlattr>.id", -1);
                  SchemeInfo scheme;
                  scheme.scheme_id = scheme_id;
                  scheme.step_id = step_id;
                  for (auto &attrib : node.second.get_child("<xmlattr>")) {
                    if (attrib.first == "ReactionStepProducts") {
                      scheme.ReactionStepProducts =
                          to_vec<unsigned int>(attrib.second.data());
                    } else if (attrib.first == "ReactionStepReactants") {
                      scheme.ReactionStepReactants =
                          to_vec<unsigned int>(attrib.second.data());
                    } else if (attrib.first ==
                               "ReactionStepObjectsAboveArrow") {
                      scheme.ReactionStepObjectsAboveArrow =
                          to_vec<unsigned int>(attrib.second.data());
                    } else if (attrib.first ==
                               "ReactionStepObjectsBelowArrow") {
                      scheme.ReactionStepObjectsBelowArrow =
                          to_vec<unsigned int>(attrib.second.data());
                    } else if (attrib.first == "ReactionStepAtomMap") {
                      scheme.ReactionStepAtomMap =
                          to_vec<unsigned int>(attrib.second.data());
                    }
                  }
                  schemes.push_back(std::move(scheme));
                }
              }
            }
          }
        }
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
                          scheme.ReactionStepReactants, fragments, mols);
        set_reaction_data("ReactionStepProducts", CDX_PRODUCT_ID, scheme,
                          scheme.ReactionStepProducts, fragments, mols);
        auto agents = scheme.ReactionStepObjectsAboveArrow;
        agents.insert(agents.end(),
                      scheme.ReactionStepObjectsBelowArrow.begin(),
                      scheme.ReactionStepObjectsBelowArrow.end());
        set_reaction_data("ReactionStepAgents", CDX_AGENT_ID, scheme, agents,
                          fragments, mols);
        // Set the Atom Maps
        int sz = scheme.ReactionStepAtomMap.size();
        if (sz % 2) {
          BOOST_LOG(rdWarningLog)
              << "CDXMLParser: Schema " << scheme.scheme_id << " step "
              << scheme.step_id
              << " ReactionStepAtomMap has odd number of entries, skipping schema..."
              << std::endl;
          continue;
        }
        CHECK_INVARIANT(sz % 2 == 0, "bad size");
        for (int i = 0; i < sz / 2; ++i) {
          unsigned int idx1 = scheme.ReactionStepAtomMap[i * 2];
          unsigned int idx2 = scheme.ReactionStepAtomMap[i * 2 + 1];
          if (atoms.find(idx1) != atoms.end()) {
            atoms[idx1]->setAtomMapNum(i + 1);
          } else {
            BOOST_LOG(rdWarningLog)
                << "CDXMLParser: Schema " << scheme.scheme_id << " step "
                << scheme.step_id
                << " ReactionStepAtomMap cannot find atom with node id " << idx1
                << "skipping schema..." << std::endl;
          }
          if (atoms.find(idx2) != atoms.end()) {
            atoms[idx2]->setAtomMapNum(i + 1);
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
  }

  // what do we do with the reaction schemes here???
  return mols;
}

std::vector<std::unique_ptr<RWMol>> MolsFromCDXMLFile(
    const std::string &fileName, const CDXMLParserParams &params) {
  std::ifstream ifs(fileName);
  if (!ifs || ifs.bad()) {
    std::ostringstream errout;
    errout << "Bad input file " << fileName;
    throw BadFileException(errout.str());
  }
  return MolsFromCDXMLDataStream(ifs, params);
}

std::vector<std::unique_ptr<RWMol>> MolsFromCDXML(
    const std::string &cdxml, const CDXMLParserParams &params) {
  std::stringstream iss(cdxml);
  return MolsFromCDXMLDataStream(iss, params);
}
}  // namespace CDXMLParser
}  // namespace v2
}  // namespace RDKit
