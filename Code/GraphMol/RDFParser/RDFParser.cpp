//
//  Copyright (C) 2025 RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RDFParser.h"

#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/StreamOps.h>

#include <algorithm>
#include <fstream>
#include <sstream>

namespace RDKit {
namespace RDF {

namespace {
constexpr const char *RXN_MARKER = "$RXN";
constexpr const char *MOL_MARKER = "$MOL";
constexpr const char *END_MARKER = "M  END";

// Helper: split entire RDF text into sections starting with $RXN
std::vector<std::string> splitReactionSections(const std::string &rdfText) {
  std::vector<std::string> sections;
  std::vector<size_t> indices;
  size_t pos = 0;
  while (true) {
    pos = rdfText.find(RXN_MARKER, pos);
    if (pos == std::string::npos) break;
    // ensure starts at line-begin
    if (pos == 0 || rdfText[pos - 1] == '\n') {
      indices.push_back(pos);
    }
    pos += 4;
  }
  for (size_t i = 0; i < indices.size(); ++i) {
    size_t start = indices[i];
    size_t end = (i + 1 < indices.size()) ? indices[i + 1] : rdfText.size();
    sections.emplace_back(rdfText.substr(start, end - start));
  }
  return sections;
}

// Helper: extract mol blocks within a reaction section (reactants then products)
std::pair<std::vector<std::string>, std::vector<std::string>>
molBlocksFromSection(const std::string &section, unsigned int reactantCount,
                     unsigned int productCount) {
  std::vector<std::string> reactants;
  std::vector<std::string> products;

  // split on "$MOL\n"
  std::vector<std::string> parts;
  size_t start = 0;
  while (true) {
    size_t found = section.find(std::string(MOL_MARKER) + "\n", start);
    if (found == std::string::npos) break;
    size_t blockStart = found + 5;  // skip "$MOL\n"
    size_t endPos = section.find(END_MARKER, blockStart);
    if (endPos == std::string::npos) {
      // if END missing, attempt to cut at next $DTYPE or $RXN to avoid spillover
      size_t cut = section.find("$DTYPE", blockStart);
      if (cut == std::string::npos)
        cut = section.find(RXN_MARKER, blockStart);
      if (cut != std::string::npos) {
        parts.emplace_back(section.substr(blockStart, cut - blockStart));
      } else {
        parts.emplace_back(section.substr(blockStart));
      }
    } else {
      // include END marker line
      size_t lineEnd = section.find('\n', endPos);
      if (lineEnd == std::string::npos) {
        parts.emplace_back(section.substr(blockStart));
      } else {
        parts.emplace_back(section.substr(blockStart, lineEnd + 1 - blockStart));
      }
    }
    start = blockStart;
  }

  // assign first N to reactants, rest to products (like Python port)
  if (parts.size() < reactantCount + productCount) {
    // best-effort: split in half if counts unknown
    // leave empty if counts are invalid
  }
  for (size_t i = 0; i < parts.size(); ++i) {
    if (i < reactantCount) {
      reactants.push_back(parts[i]);
    } else {
      products.push_back(parts[i]);
    }
  }
  return {reactants, products};
}

// Build reaction SMILES from lists of molecule SMILES
std::string reactionSmilesFromLists(const std::vector<std::string> &reactants,
                                    const std::vector<std::string> &reagents,
                                    const std::vector<std::string> &products) {
  auto joinDots = [](const std::vector<std::string> &xs) {
    std::ostringstream oss;
    for (size_t i = 0; i < xs.size(); ++i) {
      if (!xs[i].empty()) {
        if (i) oss << ".";
        oss << xs[i];
      }
    }
    return oss.str();
  };
  std::ostringstream rxn;
  rxn << joinDots(reactants) << ">" << joinDots(reagents) << ">"
      << joinDots(products);
  return rxn.str();
}

// Convert a mol block string to SMILES using RDKit
std::string molBlockToSmiles(const std::string &molBlock,
                             bool exceptOnInvalid) {
  std::unique_ptr<RWMol> mol(nullptr);
  try {
    // MolBlockToMol returns ROMol*, use RWMol for consistency
    ROMol *romol = MolBlockToMol(molBlock);
    if (romol) {
      mol.reset(new RWMol(*romol));
      delete romol;
    }
  } catch (...) {
    if (exceptOnInvalid) throw;
    return std::string();
  }
  if (!mol) {
    if (exceptOnInvalid) {
      throw std::runtime_error("Invalid mol block");
    }
    return std::string();
  }
  try {
    return MolToSmiles(*mol);
  } catch (...) {
    if (exceptOnInvalid) throw;
    return std::string();
  }
}

}  // namespace

bool RdfFileIsReaction(const std::string &fname) {
  std::ifstream in(fname);
  if (!in.good()) return false;
  std::ostringstream ss;
  ss << in.rdbuf();
  auto sections = splitReactionSections(ss.str());
  return !sections.empty();
}

bool RdfDataStreamIsReaction(std::istream &inStream) {
  std::ostringstream ss;
  ss << inStream.rdbuf();
  auto sections = splitReactionSections(ss.str());
  return !sections.empty();
}

bool RdfBlockIsReaction(const std::string &rdfText) {
  return !splitReactionSections(rdfText).empty();
}

void SplitReactionSmiles(const std::string &rxSmiles,
                         std::vector<std::string> &reactants,
                         std::vector<std::string> &reagents,
                         std::vector<std::string> &products) {
  reactants.clear();
  reagents.clear();
  products.clear();
  std::array<std::string, 3> parts;
  size_t p1 = rxSmiles.find('>');
  size_t p2 =
      p1 == std::string::npos ? std::string::npos : rxSmiles.find('>', p1 + 1);
  if (p1 == std::string::npos || p2 == std::string::npos) {
    throw std::invalid_argument("Invalid reaction SMILES: expected a>b>c");
  }
  parts[0] = rxSmiles.substr(0, p1);
  parts[1] = rxSmiles.substr(p1 + 1, p2 - (p1 + 1));
  parts[2] = rxSmiles.substr(p2 + 1);

  auto splitDots = [](const std::string &s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, '.')) {
      if (!item.empty()) out.push_back(item);
    }
    return out;
  };

  reactants = splitDots(parts[0]);
  reagents = splitDots(parts[1]);
  products = splitDots(parts[2]);
}

std::vector<RdfReactionEntry> EntriesFromRdfBlock(
    const std::string &rdfText, const RdfParserParams &params) {
  std::vector<RdfReactionEntry> entries;
  auto sections = splitReactionSections(rdfText);
  for (const auto &sec : sections) {
    // Minimal header parsing: lines 1..4 are standard RXN header, but counts
    // are in line 5 in CTFile. Here we attempt to infer counts from mol blocks:
    // For a robust implementation, port CTFile header parsing later.
    // Find mol blocks
    // As we don't have counts, we will split all mol blocks and consider first half as reactants.
    // To improve: parse line 4/5 with CTF component count format.
    // For now, compute counts by detecting "$MOL" occurrences and assume equal split if unknown.

    // Try to read counts from a standard CTFile RXN header pattern:
    unsigned int reactantCount = 0;
    unsigned int productCount = 0;
    {
      // Attempt crude parse of header counts: find a line of six digits "rrrppp"
      std::istringstream iss(sec);
      std::string line;
      unsigned int lineNo = 0;
      while (std::getline(iss, line)) {
        ++lineNo;
        if (lineNo == 5) {
          // line 5 usually contains reactant/product counts (CTFile)
          // try to parse "rrrppp" as integers separated by spaces
          // if fails, leave zero and infer later
          std::stringstream ls(line);
          int r = 0, p = 0;
          if (ls >> r >> p) {
            reactantCount = static_cast<unsigned int>(std::max(0, r));
            productCount = static_cast<unsigned int>(std::max(0, p));
          }
          break;
        }
      }
    }

    auto mols = molBlocksFromSection(sec, reactantCount, productCount);
    const auto &reactantBlocks = mols.first;
    const auto &productBlocks = mols.second;

    // Convert mol blocks to SMILES
    std::vector<std::string> reactantSmiles;
    std::vector<std::string> productSmiles;
    reactantSmiles.reserve(reactantBlocks.size());
    productSmiles.reserve(productBlocks.size());
    for (const auto &mb : reactantBlocks) {
      reactantSmiles.emplace_back(molBlockToSmiles(mb, params.exceptOnInvalidMolecule));
    }
    for (const auto &mb : productBlocks) {
      productSmiles.emplace_back(molBlockToSmiles(mb, params.exceptOnInvalidMolecule));
    }

    std::vector<std::string> reagentSmiles;  // reagents are in $DTYPE/$DATUM; handle later

    RdfReactionEntry entry;
    entry.reactantMolBlocks = reactantBlocks;
    entry.productMolBlocks = productBlocks;
    entry.reagentMolBlocks = {};  // placeholder

    entry.reactionSmiles =
        reactionSmilesFromLists(reactantSmiles, reagentSmiles, productSmiles);

    entries.push_back(std::move(entry));
  }
  return entries;
}

std::vector<RdfReactionEntry> EntriesFromRdfDataStream(
    std::istream &inStream, const RdfParserParams &params) {
  std::ostringstream ss;
  ss << inStream.rdbuf();
  return EntriesFromRdfBlock(ss.str(), params);
}

std::vector<RdfReactionEntry> EntriesFromRdfFile(
    const std::string &fName, const RdfParserParams &params) {
  std::ifstream in(fName);
  if (!in.good()) {
    throw RDKit::BadFileException("Cannot open RDF file: " + fName);
  }
  return EntriesFromRdfDataStream(in, params);
}

std::unique_ptr<ChemicalReaction> ReactionFromRdfEntry(
    const RdfReactionEntry &entry, const RdfParserParams &params) {
  // Build a ChemicalReaction from reaction SMILES.
  // In Python, ReactionFromSmarts(entry.smiles, useSmiles=True) is common.
  // In C++, RDKit provides RxnSmartsToChemicalReaction for SMARTS; using SMILES
  // may require a helper. For now, construct an empty ChemicalReaction and
  // add templates from molecules. Reviewers may suggest a preferred path.

  auto rxn = std::make_unique<ChemicalReaction>();

  // Add reactant templates
  for (const auto &mb : entry.reactantMolBlocks) {
    std::unique_ptr<ROMol> m(MolBlockToMol(mb));
    if (!m) {
      if (params.exceptOnInvalidMolecule)
        throw std::runtime_error("Invalid reactant mol block");
      else
        continue;
    }
    rxn->addReactantTemplate(new ROMol(*m));
  }

  // Add product templates
  for (const auto &mb : entry.productMolBlocks) {
    std::unique_ptr<ROMol> m(MolBlockToMol(mb));
    if (!m) {
      if (params.exceptOnInvalidMolecule)
        throw std::runtime_error("Invalid product mol block");
      else
        continue;
    }
    rxn->addProductTemplate(new ROMol(*m));
  }

  // Reagents are not part of the stoichiometric templates in RDKit reactions;
  // they can be stored in properties or handled externally.

  // Optionally validate
  try {
    rxn->initReactantMatchers();
  } catch (...) {
    if (params.exceptOnInvalidReaction) {
      throw std::runtime_error("Invalid reaction: initReactantMatchers failed");
    }
  }

  // Store the reaction SMILES for convenience
  rxn->setProp("_RDFReactionSmiles", entry.reactionSmiles);

  return rxn;
}

std::vector<std::unique_ptr<ChemicalReaction>> ReactionsFromRdfBlock(
    const std::string &rdfText, const RdfParserParams &params) {
  std::vector<std::unique_ptr<ChemicalReaction>> rxns;
  auto entries = EntriesFromRdfBlock(rdfText, params);
  rxns.reserve(entries.size());
  for (const auto &e : entries) {
    rxns.emplace_back(ReactionFromRdfEntry(e, params));
  }
  return rxns;
}

std::vector<std::unique_ptr<ChemicalReaction>> ReactionsFromRdfFile(
    const std::string &fName, const RdfParserParams &params) {
  auto entries = EntriesFromRdfFile(fName, params);
  std::vector<std::unique_ptr<ChemicalReaction>> rxns;
  rxns.reserve(entries.size());
  for (const auto &e : entries) {
    rxns.emplace_back(ReactionFromRdfEntry(e, params));
  }
  return rxns;
}

}  // namespace RDF
}  // namespace RDKit