// Minimal tests for RDFParse module using an inline RDF example.
//
// Add this test to the GraphMol test target in CMake.
//
// NOTE: This assumes the new RDFParse API exists:
//   - RDKit::RDF::RdfBlockIsReaction
//   - RDKit::RDF::EntriesFromRdfBlock
//   - RDKit::RDF::ReactionFromRdfEntry
//
// If these functions are in a different namespace or path, update includes accordingly.

#include <GraphMol/ChemReactions/Reaction.h>
#include <catch2/catch_all.hpp>

#include "../RDFParse/RDFParser.h"  // adjust include path as needed

using RDKit::RDF::RdfBlockIsReaction;
using RDKit::RDF::EntriesFromRdfBlock;
using RDKit::RDF::ReactionFromRdfEntry;
using RDKit::RDF::RdfParserParams;

static const char *kSampleRdf = R"RDF($RDFILE 1
$DATM 12/9/2025, 1:46:53 PM
$RFMT
$RXN

      RDKit

  2  1
$MOL

     RDKit          2D

 10 10  0  0  0  0  0  0  0  0999 V2000
   -0.2579   -0.8086    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6063   -0.3056    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6027    0.6944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4741   -0.8024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3383   -0.2994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2063   -0.7964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0703   -0.2932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0669    0.7068    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1989    1.2038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3349    0.7006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
 10  5  1  0
M  END
$MOL

     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
    6.6267    1.2014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6287    0.2014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7637   -0.3004    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8967    0.1980    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7657   -1.3004    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  3  5  1  0
M  END
$MOL

     RDKit          2D

 15 15  0  0  0  0  0  0  0  0999 V2000
    7.2537   -1.4305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.1201   -0.9313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.1209    0.0687    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9857   -1.4319    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.8519   -0.9327    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7177   -1.4333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.5841   -0.9341    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.5849    0.0661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.4513    0.5653    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   13.3169    0.0645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.4521    1.5653    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.5865    2.0659    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   13.3185    2.0645    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.7193    0.5665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.8527    0.0673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  1  0
  9 11  1  0
 11 12  2  0
 11 13  1  0
  8 14  1  0
 14 15  2  0
 15  5  1  0
M  END

$DTYPE Name
$DATUM Published reaction
$DTYPE Reference
$DATUM [10.1002/ANIE.200903055]: https://doi.org/10.1002/ANIE.200903055
$DTYPE Reaction Conditions
$DATUM Not available
$DTYPE SMILES
$DATUM CC(C)Cc1ccccc1.CCC(=O)O>>CC(C)Cc1ccc(C(C)C(=O)O)cc1
$DTYPE Protections
$DATUM Not available
$DTYPE Inventory
$DATUM Not available

END)RDF";

TEST_CASE("RDFParse detects and parses simple RDF", "[RDFParse]") {
  REQUIRE(RdfBlockIsReaction(kSampleRdf) == true);

  RdfParserParams params;
  params.exceptOnInvalidMolecule = true;
  params.exceptOnInvalidReaction = true;

  auto entries = EntriesFromRdfBlock(kSampleRdf, params);
  REQUIRE(entries.size() == 1);

  const auto &entry = entries.front();

  // Basic sanity: we expect three mol blocks in the section (2 reactants, 1 product)
  REQUIRE(entry.reactantMolBlocks.size() == 2);
  REQUIRE(entry.productMolBlocks.size() == 1);

  // Reaction SMILES should have two '>' separators
  REQUIRE(entry.reactionSmiles.find('>') != std::string::npos);
  size_t gt1 = entry.reactionSmiles.find('>');
  REQUIRE(gt1 != std::string::npos);
  size_t gt2 = entry.reactionSmiles.find('>', gt1 + 1);
  REQUIRE(gt2 != std::string::npos);

  // Build a ChemicalReaction object from the entry
  auto rxn = ReactionFromRdfEntry(entry, params);
  REQUIRE(rxn != nullptr);

  // Check template counts align with reactants/products
  REQUIRE(rxn->getNumReactantTemplates() == entry.reactantMolBlocks.size());
  REQUIRE(rxn->getNumProductTemplates() == entry.productMolBlocks.size());

  // Optional: the property storing reaction SMILES should exist
  std::string storedSmiles;
  REQUIRE(rxn->getPropIfPresent("_RDFReactionSmiles", storedSmiles));
  REQUIRE_FALSE(storedSmiles.empty());
}