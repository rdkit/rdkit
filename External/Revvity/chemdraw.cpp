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
#include "fragment.h"
#include "reaction.h"
#include "utils.h"

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

// The parsing of fragments needed to be moved to a recursive function since
// they may be embedded further in the document, i.e. a group may hold multiple
//  fragments
//
// Additionally, a grouped_fragments map is included to group fragments together
// for the purposes of reactions.
//
// Ungrouped fragments will end up as vectors of size 1 in the grouped_fragement
// list. The reaction schemes in the CDXML docs appear to use the fragment id
// for ungrouped fragments and the grouped id for grouped fragments, so the
// grouped_fragments holds both for ease of bookkeeping.
void visit_children(
    CDXObject &node, std::map<unsigned int, Atom *> &ids,
    std::vector<std::unique_ptr<RWMol>>
        &mols,  // All molecules found in the doc
    std::map<unsigned int, size_t>
        &fragment_lookup,  // fragment.id->molecule index
    std::map<unsigned int, std::vector<int>>
        &grouped_fragments,            // grouped.id -> [fragment.id]
    std::vector<ReactionInfo> &schemes,  // reaction schemes found
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
      schemes.emplace_back(scheme);
      /*
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
      */
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
  std::vector<ReactionInfo> schemes;  // reaction schemes found
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
  for(auto &scheme: schemes) {
    scheme.set_reaction_steps(grouped_fragments, mols);
  }
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
