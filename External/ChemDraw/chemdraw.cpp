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

#include "ChemDrawStartInclude.h"
#include "chemdraw/CDXMLParser.h"
#include "chemdraw/CDXStdObjects.h"
#include "ChemDrawEndInclude.h"

#include "bracket.h"
#include "chemdraw.h"
#include "chemdraw_doc.h"
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
#include <boost/algorithm/string.hpp>
#include <filesystem>

// #define DEBUG 1
namespace {
using namespace RDKit;
using namespace RDKit::v2;
using namespace RDKit::ChemDraw;
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
    CDXObject &node, PageData &pagedata,
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
    if (id == kCDXObj_Fragment) {
      std::unique_ptr<RWMol> mol = std::make_unique<RWMol>();
      if (!parseFragment(*mol, (CDXFragment &)(*frag.second), pagedata,
                          missing_frag_id)) {
        continue;
      }
      unsigned int frag_id = mol->getProp<int>(CDX_FRAG_ID);
      pagedata.fragmentLookup[frag_id] = pagedata.mols.size();
      if (group_id != -1) {
        pagedata.groupedFragments[group_id].push_back(frag_id);
      } else {
        pagedata.groupedFragments[frag_id].push_back(frag_id);
      }

      if (mol->hasProp(NEEDS_FUSE)) {
        mol->clearProp(NEEDS_FUSE);
        std::unique_ptr<ROMol> fused;
        try {
          replaceFragments(*mol);
          fused = molzip(*mol, molzip_params);
        } catch (Invar::Invariant &) {
          BOOST_LOG(rdWarningLog) << "Failed fusion of fragment skipping... "
                                  << frag_id << std::endl;
          // perhaps have an option to extract all fragments?
          // mols.push_back(std::move(mol));
          continue;
        }
        fused->setProp<int>(CDX_FRAG_ID, static_cast<int>(frag_id));
        pagedata.mols.emplace_back(dynamic_cast<RWMol *>(fused.release()));
      } else {
        pagedata.mols.push_back(std::move(mol));
      }
      RWMol *res = pagedata.mols.back().get();
      auto conf = std::make_unique<Conformer>(res->getNumAtoms());
      conf->set3D(false);

      bool hasConf = false;
      bool is3D = false;
      for (auto &atm : res->atoms()) {
        RDGeom::Point3D p{0.0, 0.0, 0.0};

        if (atm->hasProp(CDX_ATOM_POS)) {
          hasConf = true;
          const std::vector<double> coord =
              atm->getProp<std::vector<double>>(CDX_ATOM_POS);

          p.x = coord[0];
          p.y = -1 * coord[1];  // CDXML uses an inverted coordinate
          // system, so we need to reverse that
          if (coord.size() == 2) {
            p.z = 0.0;
          } else {
            p.z = coord[2];
            is3D = true;
          }
        }
        conf->setAtomPos(atm->getIdx(), p);
        atm->clearProp(CDX_ATOM_POS);
      }

      if (hasConf) {
        if (!is3D) {
          scaleBonds(*res, *conf, RDKIT_DEPICT_BONDLENGTH, bondLength);
        }
        conf->set3D(is3D);

        auto confidx = res->addConformer(conf.release());

        if (is3D) {
          res->updatePropertyCache(false);
          MolOps::assignChiralTypesFrom3D(*res, confidx, true);
        } else {
          MolOps::assignChiralTypesFromBondDirs(*res, confidx, true);
        }
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
            MolOps::removeHs(*res);
          } else {
            MolOps::sanitizeMol(*res);
            MolOps::detectBondStereochemistry(*res);
          }

        } catch (...) {
          BOOST_LOG(rdWarningLog)
              << "CDXMLParser: failed sanitizing skipping fragment " << frag_id
              << std::endl;
          pagedata.mols.pop_back();
          continue;
        }
        MolOps::assignStereochemistry(*res, true, true, true);
        // Sometimes ChemDraw just marks with R and S, so let's assign
        //  these as long as they were not already determined
        checkChemDrawTetrahedralGeometries(*res);
      } else {
        MolOps::detectBondStereochemistry(*res);
      }
    } else if (id == kCDXObj_ReactionScheme) {  // get the reaction info
      CDXReactionScheme &scheme = (CDXReactionScheme &)(*frag.second);
      pagedata.schemes.emplace_back(scheme);
      /*
      int scheme_id = scheme.GetObjectID();   //frag.second.template
      get<int>("<xmlattr>.id", -1); for (auto &rxnNode :
      scheme.ContainedObjects()) { CDXDatumID type_id =
      (CDXDatumID)rxnNode.second->GetTag(); if (type_id == kCDXObj_ReactionStep)
      { CDXReactionStep &step = (CDXReactionStep&)(*rxnNode.second); auto
      step_id = step.GetObjectID(); SchemeInfo scheme; scheme.scheme_id =
      scheme_id; scheme.step_id = step_id; scheme.ReactionStepProducts =
      step.m_products; scheme.ReactionStepReactants = step.m_reactants;
          scheme.ReactionStepObjectsBelowArrow = step.m_objectsBelowArrow;
          scheme.ReactionStepAtomMap = step.m_aamap;
          schemes.push_back(scheme);
        }
      }
      */
    } else if (id == kCDXObj_Group) {
      CDXGroup &group = (CDXGroup &)(*frag.second);
      group_id = frag.second->GetObjectID();
      visit_children(group, pagedata, missing_frag_id, bondLength, params,
                     group_id);
    } else if (id == kCDXObj_BracketedGroup) {
      CDXBracketedGroup &bracketgroup = (CDXBracketedGroup &)(*frag.second);
      parseBracket(bracketgroup, pagedata);
    }
  }
}

std::unique_ptr<CDXDocument> streamToCDXDocument(std::istream &inStream,
                                                 CDXFormat format) {
  if (format == CDXFormat::CDXML) {
    CDXMLParser parser;
    // populate tree structure pt
    std::string data = std::string(std::istreambuf_iterator<char>(inStream),
                                   std::istreambuf_iterator<char>());
    const bool HaveAllXml = true;
    if (XML_STATUS_OK != parser.XML_Parse(data.c_str(),
                                          static_cast<int>(data.size()),
                                          HaveAllXml)) {
      auto error = XML_GetErrorCode(parser);
      std::stringstream msg;
      msg << "Failed parsing XML with error code " << error;
      BOOST_LOG(rdErrorLog) << msg.str() << std::endl;
      throw FileParseException(msg.str());
    }

    return parser.ReleaseDocument();
  } else {
    CDXistream input(inStream);
    const bool doThrow = true;
    std::unique_ptr<CDXDocument> doc(CDXReadDocFromStorage(input, doThrow));
    return doc;
  }
}

// may raise FileParseException
std::vector<std::unique_ptr<RWMol>> molsFromCDXMLDataStream(
    std::istream &inStream, const ChemDrawParserParams &params) {
  std::unique_ptr<CDXDocument> document =
    streamToCDXDocument(inStream, params.format);
  if (!document) {
    // error
    return std::vector<std::unique_ptr<RWMol>>();
  }
  PageData pagedata;
  auto bondLength = document->m_bondLength;

  int missing_frag_id = -1;
  for (auto node : document->ContainedObjects()) {
    CDXDatumID id = (CDXDatumID)node.second->GetTag();
    switch (id) {
      case kCDXObj_Page:
        visit_children(*node.second, pagedata, missing_frag_id, bondLength,
                       params);
        break;
      default:
        break;
    }
  }
  for (auto &scheme : pagedata.schemes) {
    scheme.set_reaction_steps(pagedata.groupedFragments, pagedata.mols);
  }
  pagedata.clearCDXProps();

  return std::move(pagedata.mols);
}
}  // namespace

namespace RDKit {
namespace ChemDraw {
std::unique_ptr<CDXDocument> ChemDrawToDocument(std::istream &inStream,
                                                CDXFormat format) {
  return streamToCDXDocument(inStream, format);
}

std::unique_ptr<CDXDocument> ChemDrawToDocument(const std::string &filename) {
  std::fstream chemdrawfile(filename);
  std::string ext = std::filesystem::path(filename).extension().string();
  boost::algorithm::to_lower(ext);
  if (ext == ".cdxml")
    return streamToCDXDocument(chemdrawfile, CDXFormat::CDXML);
  else if (ext == ".cdx") {
    return streamToCDXDocument(chemdrawfile, CDXFormat::CDX);
  }
  std::string msg =
      std::string("Unknoen filetype ") +
      (std::string)std::filesystem::path(filename).extension().string();
  throw FileParseException(msg.c_str());
}
}

namespace v2 {
std::vector<std::unique_ptr<RWMol>> MolsFromChemDrawDataStream(
    std::istream &inStream, const ChemDrawParserParams &params) {
  auto chemdrawmols = molsFromCDXMLDataStream(inStream, params);
  std::vector<std::unique_ptr<RWMol>> mols;
  mols.reserve(chemdrawmols.size());
  for (auto &mol : chemdrawmols) {
    RWMol *m = (RWMol *)mol.release();
    mols.push_back(std::unique_ptr<RWMol>(m));
  }
  return mols;
}

std::vector<std::unique_ptr<RWMol>> MolsFromChemDrawBlock(
    const std::string &block, const ChemDrawParserParams &params) {
  std::stringstream ss;
  ss << block;
  return MolsFromChemDrawDataStream(ss, params);
}

std::vector<std::unique_ptr<RWMol>> MolsFromChemDrawFile(
    const std::string &filename, const ChemDrawParserParams &params) {
  CDXMLParser parser;
  std::vector<std::unique_ptr<RWMol>> mols;

  std::fstream chemdrawfile(filename);  // FIX ME CHECK CDX versus CDXML
  if (!chemdrawfile) {
    throw BadFileException(filename + " does not exist");
    return mols;
  }
  auto chemdrawmols = molsFromCDXMLDataStream(chemdrawfile, params);

  mols.reserve(chemdrawmols.size());
  for (auto &mol : chemdrawmols) {
    RWMol *m = (RWMol *)mol.release();
    mols.push_back(std::unique_ptr<RWMol>(m));
  }
  return mols;
}
}
}  // namespace RDKit
