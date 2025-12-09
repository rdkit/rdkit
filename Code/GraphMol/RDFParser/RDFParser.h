//
//  Copyright (C) 2025 RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_RDFPARSER_H
#define RD_RDFPARSER_H

#include <RDGeneral/export.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>

#include <istream>
#include <memory>
#include <string>
#include <vector>

namespace RDKit {
namespace RDF {

// Basic container for a parsed RDF reaction section
struct RDKIT_RDFPARSER_EXPORT RdfReactionEntry {
  // canonical reaction SMILES: reactants>reagents>products
  std::string reactionSmiles;

  // Reactant/product/reagent molecules as mol blocks
  std::vector<std::string> reactantMolBlocks;
  std::vector<std::string> productMolBlocks;
  std::vector<std::string> reagentMolBlocks;

  // Optional metadata (minimal scaffold; can be expanded)
  std::string reactionName;
  std::string comment;
  std::string registryNumber;
};

// Parameters controlling parsing behavior
struct RDKIT_RDFPARSER_EXPORT RdfParserParams {
  bool exceptOnInvalidMolecule = true;
  bool exceptOnInvalidReaction = true;
  bool parseConditions = true;
};

// Quick checks
RDKIT_RDFPARSER_EXPORT bool RdfFileIsReaction(const std::string &fname);
RDKIT_RDFPARSER_EXPORT bool RdfDataStreamIsReaction(std::istream &inStream);
RDKIT_RDFPARSER_EXPORT bool RdfBlockIsReaction(const std::string &rdfText);

// Parse helpers returning entries
RDKIT_RDFPARSER_EXPORT std::vector<RdfReactionEntry> EntriesFromRdfDataStream(
    std::istream &inStream, const RdfParserParams &params = RdfParserParams());
RDKIT_RDFPARSER_EXPORT std::vector<RdfReactionEntry> EntriesFromRdfBlock(
    const std::string &rdfText,
    const RdfParserParams &params = RdfParserParams());
RDKIT_RDFPARSER_EXPORT std::vector<RdfReactionEntry> EntriesFromRdfFile(
    const std::string &fName, const RdfParserParams &params = RdfParserParams());

// Convenience builders (similar to MarvinParser)
RDKIT_RDFPARSER_EXPORT std::unique_ptr<ChemicalReaction>
ReactionFromRdfEntry(const RdfReactionEntry &entry,
                     const RdfParserParams &params = RdfParserParams());

RDKIT_RDFPARSER_EXPORT std::vector<std::unique_ptr<ChemicalReaction>>
ReactionsFromRdfBlock(const std::string &rdfText,
                      const RdfParserParams &params = RdfParserParams());

RDKIT_RDFPARSER_EXPORT std::vector<std::unique_ptr<ChemicalReaction>>
ReactionsFromRdfFile(const std::string &fName,
                     const RdfParserParams &params = RdfParserParams());

// Utility to split reaction SMILES into components
RDKIT_RDFPARSER_EXPORT void SplitReactionSmiles(
    const std::string &rxSmiles, std::vector<std::string> &reactants,
    std::vector<std::string> &reagents, std::vector<std::string> &products);

}  // namespace RDF
}  // namespace RDKit

#endif  // RD_RDFPARSER_H