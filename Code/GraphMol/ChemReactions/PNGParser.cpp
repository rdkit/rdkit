//
//  Copyright (c) 2020 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/FileParsers/PNGParser.h>
#include <boost/algorithm/string.hpp>

namespace RDKit {

namespace PNGData {
const std::string rxnSmilesTag = "ReactionSmiles";
const std::string rxnSmartsTag = "ReactionSmarts";
const std::string rxnRxnTag = "ReactionRxn";
const std::string rxnPklTag = "rdkitReactionPKL";
}  // namespace PNGData

std::string addChemicalReactionToPNGStream(const ChemicalReaction &rxn,
                                           std::istream &iStream,
                                           bool includePkl, bool includeSmiles,
                                           bool includeSmarts,
                                           bool includeRxn) {
  std::vector<std::pair<std::string, std::string>> metadata;
  if (includePkl) {
    std::string pkl;
    ReactionPickler::pickleReaction(rxn, pkl);
    metadata.push_back(std::make_pair(augmentTagName(PNGData::rxnPklTag), pkl));
  }
  if (includeSmiles) {
    std::string smi = ChemicalReactionToRxnSmiles(rxn);
    metadata.push_back(
        std::make_pair(augmentTagName(PNGData::rxnSmilesTag), smi));
  }
  if (includeSmarts) {
    std::string smi = ChemicalReactionToRxnSmarts(rxn);
    metadata.push_back(
        std::make_pair(augmentTagName(PNGData::rxnSmartsTag), smi));
  }
  if (includeRxn) {
    std::string mb = ChemicalReactionToRxnBlock(rxn);
    metadata.push_back(std::make_pair(augmentTagName(PNGData::rxnRxnTag), mb));
  }
  return addMetadataToPNGStream(iStream, metadata);
};

ChemicalReaction *PNGStreamToChemicalReaction(std::istream &inStream) {
  ChemicalReaction *res = nullptr;
  auto metadata = PNGStreamToMetadata(inStream);
  bool formatFound = false;
  for (const auto &pr : metadata) {
    if (boost::starts_with(pr.first, PNGData::rxnPklTag)) {
      res = new ChemicalReaction(pr.second);
      formatFound = true;
    } else if (boost::starts_with(pr.first, PNGData::rxnSmilesTag)) {
      std::map<std::string, std::string> *replacements = nullptr;
      bool useSmiles = true;
      res = RxnSmartsToChemicalReaction(pr.second, replacements, useSmiles);
      formatFound = true;
    } else if (boost::starts_with(pr.first, PNGData::rxnSmartsTag)) {
      res = RxnSmartsToChemicalReaction(pr.second);
      formatFound = true;
    } else if (boost::starts_with(pr.first, PNGData::rxnRxnTag)) {
      res = RxnBlockToChemicalReaction(pr.second);
      formatFound = true;
    }
    if (formatFound) {
      break;
    }
  }
  if (!formatFound) {
    throw FileParseException("No suitable metadata found.");
  }
  return res;
}

}  // namespace RDKit
