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
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/BadFileException.h>
#include <fstream>
#include <sstream>
#include "MolFileStereochem.h"

// TODO
//  I'm currently using atom map numbers to join structures
//   this is not very forward thinking as chemdraw can also have
//   atom map numbers in addition to the ids in the xml file.
//  [] enhance molzip to be able to zip on on an atom tag
//  [] what happens to stereo chemistry when the fuse_frags below is called
//    get cases in the "wild" and see what the expectation is.
//  [] add coordinates for the atoms.  Might not work for nicknames like Boc?

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
const std::string CDX_ATOM_ID("CDX_ATOM_ID");
const std::string CDX_BOND_ID("CDX_BOND_ID");
const std::string CDX_BOND_ORDERING("DXML_BOND_ORDERING");

struct BondInfo {
    int bond_id;
    int start;
    int end;
    int order;
    std::string display;
    Bond::BondType getBondType() {
        if(order && order <= 4)
            return static_cast<Bond::BondType>(order);
        // handle aromatic bond type
        return Bond::BondType::UNSPECIFIED;
    }
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

unsigned int get_fuse_label(Atom*atm) {
    // return atm->getAtomMapNum(); easier debugging
    unsigned int label = 0; // default is no label
    atm->getPropIfPresent<unsigned int>(FUSE_LABEL, label);
    return label;
}

void set_fuse_label(Atom*atm, unsigned int idx) {
    //atm->setAtomMapNum(idx); //for debugging
    if(idx) {
        atm->setProp<unsigned int>(FUSE_LABEL, idx);
    } else {
        atm->clearProp(FUSE_LABEL);
    }
}

template<class T>
std::vector<T> to_vec(const std::string &s) {
    std::vector<T> n;
    std::stringstream ss(s);
    std::copy(std::istream_iterator<T>(ss), std::istream_iterator<T>(), std::back_inserter(n));
    return n;
}

void parse_fragment(RWMol &mol,
                    ptree &frag,
		            std::map<unsigned int, Atom *> &ids,
                    int external_attachment = -1) {
  // XXX Need to put the fragid on the molecule so we can properly make reactions
    
  int frag_id = frag.get<int>("<xmlattr>.id", -1);
  if(frag_id == -1) {
      CHECK_INVARIANT(frag_id >= 0, "Invalid or missing id from CDXML fragment");
  }
  mol.setProp<int>(CDXML_FRAG_ID, frag_id);
  int elemno = -1;
  int fragmap = 10000;
  bool fuse_fragments = false;
  // for atom in frag
  int atom_id;
  std::vector<BondInfo> bonds;
  std::map<unsigned int, Bond*> bond_ids;
  // nodetypes = https://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/properties/Node_Type.htm
  for(auto &node: frag) {
    if(node.first == "n") { // atom node
      int elemno = 6; // default to carbon
      int num_hydrogens = 0;
      int charge = 0;
      int atommap = 0;
      int mergeparent = -1;
      int rgroup_num = -1;
      int isotope = 0;
      bool has_atom_stereo = false;
      std::vector<int> bond_ordering;
      std::vector<double> atom_coords;
      std::string nodetype = "";
      for(auto &attr: node.second.get_child("<xmlattr>")) {
          if(attr.first == "id") {
              atom_id = stoi(attr.second.data());
          } else if (attr.first == "Element") {
              elemno = stoi(attr.second.data());
          } else if(attr.first == "NumHydrogens") {
              num_hydrogens = stoi(attr.second.data());
          } else if (attr.first == "Charge") {
              charge = stoi(attr.second.data());
          } else if (attr.first == "Isotope") {
              isotope = stoi(attr.second.data());
          } else if(attr.first == "NodeType") {
              nodetype = attr.second.data();
              if (nodetype == "Nickname" || nodetype == "Fragment") {
                  elemno = 0;
                  atommap = atom_id;
              } else if (nodetype == "ExternalConnectionPoint") {
                  PRECONDITION(external_attachment > 0, "External Connection Point is not set");
                  elemno = 0;
                  atommap = external_attachment;
                  mergeparent = external_attachment;
              } else if (nodetype == "GenericNickname") {
                  // RGroups for example?
                  for(auto &tnode : node.second) {
                      if(tnode.first == "t") {
                          for(auto &snode : tnode.second) {
                              if(snode.first == "s") {
                                  auto s = snode.second.data();
                                  if(s.size()  && s[0] == 'R') {
                                      rgroup_num = stoi(s.substr(1));
                                      elemno = 0;
                                  }
                                  break;
                              }
                          }
                          break;
                      }
                  }
              } else if (nodetype == "ElementList") {
                  // query atom?
              }
          } else if (attr.first == "Geometry") {
              if (attr.second.data() == "Tetrahedral") {
                  has_atom_stereo = true;
              }
          } else if (attr.first == "BondOrdering") {
              std::vector<int> order = to_vec<int>(attr.second.data());
              // remove any idx that is 0, they can appear anywhere in the BondOrdering
              //  I really don't know why...
              for(auto idx : order) {
                  if(idx) {
                      bond_ordering.push_back(idx);
                  }
              }
          } else if (attr.first == "p") {
              atom_coords = to_vec<double>(attr.second.data());
          }
      }
      // add the atom
      Atom *rd_atom = new Atom(elemno);
      rd_atom->setFormalCharge(charge);
      rd_atom->setNumExplicitHs(num_hydrogens);
      rd_atom->setIsotope(isotope);
      if(rgroup_num >= 0) {
          rd_atom->setAtomMapNum(rgroup_num);
      }
      set_fuse_label(rd_atom, atommap);
      if(mergeparent > 0) {
          rd_atom->setProp<int>("MergeParent", mergeparent);
      }
      if(has_atom_stereo) {
          rd_atom->setProp<std::vector<int>>(CDX_BOND_ORDERING, bond_ordering);
      }
      
      if(ids.find(atom_id) != ids.end()) {
          // error fail processing
      }
      rd_atom->setProp<std::vector<double>>(CDX_ATOM_POS, atom_coords);
      rd_atom->setProp<unsigned int>(CDX_ATOM_ID, atom_id);
      ids[atom_id] = rd_atom;
      const bool updateLabels=true;
      const bool takeOwnership=true;
      auto idx = mol.addAtom(rd_atom, updateLabels, takeOwnership);
        assert (idx == rd_atom->getIdx());
      
      if (nodetype == "Nickname" || nodetype == "Fragment") {
          for(auto &fragment: node.second) {
              if(fragment.first == "fragment") {
                  parse_fragment(mol, fragment.second, ids, atom_id);
                  mol.setProp<bool>(NEEDS_FUSE, true);
                  mol.setProp<unsigned int>(CDXML_FRAG_ID, frag_id);
              }
          }
      }
    } else if (node.first == "b") { // bond
          int bond_id = -1;
          int start_atom = -1;
          int end_atom = -1;
          int order = 1;
          std::string display;
          for(auto &attr: node.second.get_child("<xmlattr>")) {
              if(attr.first == "id") {
                  bond_id = stoi(attr.second.data());
              } else if(attr.first == "B") {
                  start_atom = stoi(attr.second.data());
              } else if (attr.first == "E") {
                  end_atom = stoi(attr.second.data());
              } else if (attr.first == "Order") {
                  order = stoi(attr.second.data());
              } else if (attr.first == "Display") { // gets wedge/hash stuff and probably more
                  display = attr.second.data();
              }
          }
        CHECK_INVARIANT(start_atom>=0 && end_atom>=0 && start_atom != end_atom, "Bad bond in CDXML");
        bonds.push_back({bond_id, start_atom, end_atom, order, display});
          
    } // end if atom or bond
  } // for node
    
  // add bonds
  for(auto &bond: bonds) {
      unsigned int bond_idx;
      if(bond.display == "WedgeEnd" || bond.display == "WedgedHashEnd") {
          // here The "END" of the bond is really our Beginning.
          // swap atom direction
          bond_idx = mol.addBond(ids[bond.end]->getIdx(), ids[bond.start]->getIdx(),
                       bond.getBondType()) - 1;
      } else {
          bond_idx = mol.addBond(ids[bond.start]->getIdx(), ids[bond.end]->getIdx(),
                       bond.getBondType()) - 1;
      }
      Bond * bnd = mol.getBondWithIdx(bond_idx);
      bnd->setProp("CDX_BOND_ID", bond.bond_id);
      if(bond.display == "WedgeEnd" || bond.display == "WedgeBegin") {
          bnd->setBondDir(Bond::BondDir::BEGINDASH);
      } else if (bond.display == "WedgedHashBegin" || bond.display == "WedgedHashEnd") {
          bnd->setBondDir(Bond::BondDir::BEGINWEDGE);
      }
      bond_ids[bond_idx] = bnd;
  }
}

void set_reaction_data(std::string type,
                       std::string prop,
                       SchemeInfo &scheme,
                       const std::vector<unsigned int> &frag_ids,
                       const std::map<unsigned int, size_t> &fragments,
                       const std::vector<std::unique_ptr<RWMol>> &mols) {
    unsigned int reagent_idx = 0;
    for(auto idx : frag_ids) {
        auto iter = fragments.find(idx);
        if (iter == fragments.end()) {
            BOOST_LOG(rdWarningLog) << "CDXMLParser: Schema " << scheme.scheme_id << " step " << scheme.step_id <<
            " " << type << " fragment " << idx << " not found in document." << std::endl;
            continue;
        }
        if(iter->second >= mols.size() ) {
            // shouldn't get here
            continue;
        }
        auto &mol = mols[iter->second];
        mol->setProp(CDX_SCHEME_ID, scheme.scheme_id);
        mol->setProp(CDX_STEP_ID, scheme.step_id);
        mol->setProp(prop, reagent_idx++);
    }
}
} // namepspace

std::vector<std::unique_ptr<RWMol>> CDXMLDataStreamToMols(
			       std::istream &inStream,
			       bool sanitize,
			       bool removeHs) {
  
  // populate tree structure pt
  using boost::property_tree::ptree;
  ptree pt;
  read_xml(inStream, pt);
  std::map<unsigned int, Atom*> ids;
  std::vector<std::unique_ptr<RWMol>> mols;
  std::map<unsigned int, size_t> fragment_lookup;
  std::vector<SchemeInfo> schemes;
    
  MolzipParams molzip_params;
  molzip_params.label = MolzipLabel::AtomProperty;
  molzip_params.atomProperty = FUSE_LABEL;
  molzip_params.enforceValenceRules = false;
  
  for( auto &cdxml : pt ) {
      if(cdxml.first == "CDXML") {
          for( auto &node : cdxml.second) {
              if (node.first == "page") {
                  for( auto &frag: node.second ) {
                      if (frag.first == "fragment") { // chemical matter
                          std::unique_ptr<RWMol> mol = std::make_unique<RWMol>();
                          parse_fragment(*mol, frag.second, ids);
                          unsigned int frag_id = mol->getProp<unsigned int>(CDXML_FRAG_ID);
                          fragment_lookup[frag_id] = mols.size();
                          if(mol->hasProp(NEEDS_FUSE)) {
                              mol->clearProp(NEEDS_FUSE);
                              auto fused = molzip(*mol, molzip_params);
                              fused->setProp<unsigned int>(CDXML_FRAG_ID, frag_id);
                              mols.emplace_back(dynamic_cast<RWMol*>(fused.release()));
                          } else {
                              mols.push_back(std::move(mol));
                          }
                          RWMol *res = mols.back().get();
                          Conformer *conf = new Conformer(res->getNumAtoms());
                            conf->set3D(false);
                          for(auto &atm: res->atoms()) {
                              const std::vector<double> coord =  atm->getProp<std::vector<double>>(CDX_ATOM_POS);
                          
                                RDGeom::Point3D p;
                                if(coord.size() == 2) {
                                    p.x = coord[0];
                                    p.y = coord[1];
                                    p.z = 0.0;
                                }
                                conf->setAtomPos(atm->getIdx(), p);
                          }
                          res->addConformer(conf);
                          DetectAtomStereoChemistry(*res, conf);
                          
                          if (sanitize) {
                            if (removeHs) {
                              MolOps::removeHs(*res, false, false);
                            } else {
                              MolOps::sanitizeMol(*res);
                            }
                            // now that atom stereochem has been perceived, the wedging
                            // information is no longer needed, so we clear
                            // single bond dir flags:

                            ClearSingleBondDirFlags(*res);
                            MolOps::detectBondStereochemistry(*res);
                            MolOps::assignStereochemistry(*res, true, true, true);
                        }
                      } else if (frag.first == "scheme") { // get the reaction info
                          int scheme_id = frag.second.get<int>("<xmlattr>.id", -1);
                          for(auto &node : frag.second) {
                              if(node.first == "step") {
                                  auto step_id = node.second.get<int>("<xmlattr>.id", -1);
                                  SchemeInfo scheme;
                                  scheme.scheme_id = scheme_id;
                                  scheme.step_id = step_id;
                                  for(auto &attrib : node.second.get_child("<xmlattr>")) {
                                      if(attrib.first == "ReactionStepProducts") {
                                          scheme.ReactionStepProducts = to_vec<unsigned int>(attrib.second.data());
                                      } else if (attrib.first == "ReactionStepReactants") {
                                          scheme.ReactionStepReactants = to_vec<unsigned int>(attrib.second.data());
                                      } else if (attrib.first == "ReactionStepObjectsAboveArrow") {
                                          scheme.ReactionStepObjectsAboveArrow = to_vec<unsigned int>(attrib.second.data());
                                      } else if (attrib.first == "ReactionStepObjectsBelowArrow") {
                                          scheme.ReactionStepObjectsBelowArrow = to_vec<unsigned int>(attrib.second.data());
                                      }  else if (attrib.first == "ReactionStepAtomMap") {
                                          scheme.ReactionStepAtomMap = to_vec<unsigned int>(attrib.second.data());
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
      if(schemes.size()) {
          std::map<unsigned int, size_t> fragments;
          std::map<unsigned int, size_t> agents;
          std::map<unsigned int, size_t> products;
          std::map<unsigned int, Atom*> atoms;
          size_t mol_idx=0;
          for(auto &mol:mols) {
              auto idx = mol->getProp<unsigned int>(CDXML_FRAG_ID);
              fragments[idx] = mol_idx++;
              for(auto &atom:mol->atoms()) {
                  unsigned int idx = atom->getProp<unsigned int>(CDX_ATOM_ID);
                  atoms[idx] = atom;
              }
          }
    
          for(auto &scheme:schemes) {
              // Set the molecule properties
              set_reaction_data("ReactionStepReactants",
                                CDX_REAGENT_ID,
                                scheme,
                                scheme.ReactionStepReactants,
                                fragments,
                                mols);
              set_reaction_data("ReactionStepProducts",
                                CDX_PRODUCT_ID,
                                scheme,
                                scheme.ReactionStepProducts,
                                fragments,
                                mols);
              auto agents = scheme.ReactionStepObjectsAboveArrow;
              agents.insert(agents.end(),
                         scheme.ReactionStepObjectsBelowArrow.begin(),
                         scheme.ReactionStepObjectsBelowArrow.end());
              set_reaction_data("ReactionStepAgents",
                                CDX_AGENT_ID,
                                scheme,
                                agents,
                                fragments,
                                mols);
              // Set the Atom Maps
              int sz = scheme.ReactionStepAtomMap.size();
              if(sz%2 == 1) {
                  BOOST_LOG(rdWarningLog) << "CDXMLParser: Schema " << scheme.scheme_id << " step " << scheme.step_id << " ReactionStepAtomMap has odd number of entries, skipping." << std::endl;
                  continue;
              }
              assert (sz%2 == 0);
              for(int i=0; i<sz/2; ++i) {
                  unsigned int idx1 = scheme.ReactionStepAtomMap[i*2];
                  unsigned int idx2 = scheme.ReactionStepAtomMap[i*2 + 1];
                  if( atoms.find(idx1) != atoms.end() ) {
                      atoms[idx1]->setAtomMapNum(i+1);
                  } else {
                      BOOST_LOG(rdWarningLog) << "CDXMLParser: Schema " << scheme.scheme_id << " step " << scheme.step_id << " ReactionStepAtomMap cannot find atom with node id " << idx1 << std::endl;
                  }
                  if( atoms.find(idx2) != atoms.end() ) {
                      atoms[idx2]->setAtomMapNum(i+1);
                  } else {
                      // XXX log error
                      BOOST_LOG(rdWarningLog) << "CDXMLParser: Schema " << scheme.scheme_id << " step " << scheme.step_id << " ReactionStepAtomMap cannot find atom with node id " << idx2 << std::endl;
                  }
              }
          }
      }
  } 

  // what do we do with the reaction schemes here???
  return mols;
}

std::vector<std::unique_ptr<RWMol>> CDXMLFileToMols(
			       const std::string &fileName,
			       bool sanitize,
			       bool removeHs) {
  std::ifstream ifs(fileName);
  if (!ifs || ifs.bad()) {
    std::ostringstream errout;
    errout << "Bad input file " << fileName;
    throw BadFileException(errout.str());
  }
  return CDXMLDataStreamToMols(ifs, sanitize, removeHs);
}
  
std::vector<std::unique_ptr<RWMol>> CDXMLToMols(
			       const std::string &cdxml,
			       bool sanitize,
			       bool removeHs) {
  std::stringstream iss(cdxml);
  return CDXMLDataStreamToMols(iss, sanitize, removeHs);
}

}
