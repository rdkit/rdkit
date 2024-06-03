//
//  Copyright (c) 2007-2024, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
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

#include <RDGeneral/export.h>
#ifndef RD_REACTIONPARSER_H_21Aug2006
#define RD_REACTIONPARSER_H_21Aug2006

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <utility>
#include <boost/format.hpp>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {
class ROMol;
class ChemicalReaction;

//! used to indicate an error in parsing reaction data
class RDKIT_CHEMREACTIONS_EXPORT ChemicalReactionParserException
    : public std::exception {
 public:
  //! construct with an error message
  explicit ChemicalReactionParserException(const char *msg) : _msg(msg) {}
  //! construct with an error message
  explicit ChemicalReactionParserException(std::string msg)
      : _msg(std::move(msg)) {}
  //! get the error message
  const char *what() const noexcept override { return _msg.c_str(); }
  ~ChemicalReactionParserException() noexcept override = default;

 private:
  std::string _msg;
};

namespace v2 {
namespace ReactionParser {
struct RDKIT_CHEMREACTIONS_EXPORT ReactionSmartsParserParams {
  bool sanitize = false; /**< sanitize the molecules after building them */
  std::map<std::string, std::string>
      replacements;          /**< allows SMILES "macros" */
  bool allowCXSMILES = true; /**< recognize and parse CXSMILES*/
  bool strictCXSMILES =
      true; /**< throw an exception if the CXSMILES parsing fails */
};
RDKIT_CHEMREACTIONS_EXPORT std::unique_ptr<ChemicalReaction> ReactionFromSmarts(
    const std::string &smarts,
    const ReactionSmartsParserParams &params = ReactionSmartsParserParams());

RDKIT_CHEMREACTIONS_EXPORT std::unique_ptr<ChemicalReaction> ReactionFromSmiles(
    const std::string &smarts,
    const ReactionSmartsParserParams &params = ReactionSmartsParserParams());
}  // namespace ReactionParser
}  // namespace v2

inline namespace v1 {
//---------------------------------------------------------------------------
//! \name Reaction SMARTS/SMILES Support
//! @{

//! Parse a string containing "Reaction SMARTS" into a ChemicalReaction
/*!
   Our definition of Reaction SMARTS is something that looks a lot like reaction
   SMILES, except that SMARTS queries are allowed on the reactant side and that
   atom-map numbers are required (at least for now)

   \param text          the SMARTS to convert

   \param replacements  a string->string map of replacement strings. \see
   SmilesToMol for more information about replacements

   \param useSmiles     if set, the SMILES parser will be used instead of the
   SMARTS parserfor the individual components

   \param allowCXSMILES     if set, any CXSMILES extensions present will be
   parsed, otherwise it will be ignored
 */
inline ChemicalReaction *RxnSmartsToChemicalReaction(
    const std::string &text,
    std::map<std::string, std::string> *replacements = nullptr,
    bool useSmiles = false, bool allowCXSMILES = true) {
  v2::ReactionParser::ReactionSmartsParserParams params;
  if (replacements) {
    params.replacements = *replacements;
  }
  params.allowCXSMILES = allowCXSMILES;
  if (useSmiles) {
    return v2::ReactionParser::ReactionFromSmiles(text, params).release();
  } else {
    return v2::ReactionParser::ReactionFromSmarts(text, params).release();
  }
}
}  // namespace v1
//! returns the reaction SMARTS for a reaction
RDKIT_CHEMREACTIONS_EXPORT std::string ChemicalReactionToRxnSmarts(
    const ChemicalReaction &rxn, const SmilesWriteParams &params);
//! \overload
inline std::string ChemicalReactionToRxnSmarts(const ChemicalReaction &rxn) {
  SmilesWriteParams params;
  params.canonical = false;
  return ChemicalReactionToRxnSmarts(rxn, params);
}

//! returns the reaction SMILES for a reaction
RDKIT_CHEMREACTIONS_EXPORT std::string ChemicalReactionToRxnSmiles(
    const ChemicalReaction &rxn,
    const SmilesWriteParams &params = SmilesWriteParams());
//! \overload
inline std::string ChemicalReactionToRxnSmiles(const ChemicalReaction &rxn,
                                               bool canonical) {
  SmilesWriteParams params;
  params.canonical = canonical;
  return ChemicalReactionToRxnSmiles(rxn, params);
}
//! @}

//---------------------------------------------------------------------------
//! \name Reaction Mol Support
//! @{

//! Parse a ROMol into a ChemicalReaction, RXN role must be set before
/*!
   Alternative to build a reaction from a molecule (fragments) which have RXN
   roles set as atom properties: common_properties::molRxnRole (1=reactant,
   2=product, 3=agent)

   \param mol           ROMol with RXN roles set
 */
RDKIT_CHEMREACTIONS_EXPORT ChemicalReaction *RxnMolToChemicalReaction(
    const ROMol &mol);

//! returns a ROMol with RXN roles used to describe the reaction
RDKIT_CHEMREACTIONS_EXPORT ROMol *ChemicalReactionToRxnMol(
    const ChemicalReaction &rxn);
//! @}

//---------------------------------------------------------------------------
//! \name MDL rxn Support
//! @{
namespace v2 {
namespace ReactionParser {

RDKIT_CHEMREACTIONS_EXPORT std::unique_ptr<ChemicalReaction>
ReactionFromRxnBlock(const std::string &rxnBlock,
                     const FileParsers::MolFileParserParams &params =
                         FileParsers::MolFileParserParams());
RDKIT_CHEMREACTIONS_EXPORT std::unique_ptr<ChemicalReaction>
ReactionFromRxnFile(const std::string &fileName,
                    const FileParsers::MolFileParserParams &params =
                        FileParsers::MolFileParserParams());
RDKIT_CHEMREACTIONS_EXPORT std::unique_ptr<ChemicalReaction>
ReactionFromRxnDataStream(std::istream &rxnStream, unsigned int &line,
                          const FileParsers::MolFileParserParams &params =
                              FileParsers::MolFileParserParams());

}  // namespace ReactionParser
}  // namespace v2
inline namespace v1 {
//! Parse a text block in MDL rxn format into a ChemicalReaction
inline ChemicalReaction *RxnBlockToChemicalReaction(const std::string &rxnBlock,
                                                    bool sanitize = false,
                                                    bool removeHs = false,
                                                    bool strictParsing = true) {
  v2::FileParsers::MolFileParserParams params;
  params.sanitize = sanitize;
  params.removeHs = removeHs;
  params.strictParsing = strictParsing;
  return v2::ReactionParser::ReactionFromRxnBlock(rxnBlock, params).release();
}
//! Parse a file in MDL rxn format into a ChemicalReaction
inline ChemicalReaction *RxnFileToChemicalReaction(const std::string &fileName,
                                                   bool sanitize = false,
                                                   bool removeHs = false,
                                                   bool strictParsing = true) {
  v2::FileParsers::MolFileParserParams params;
  params.sanitize = sanitize;
  params.removeHs = removeHs;
  params.strictParsing = strictParsing;
  return v2::ReactionParser::ReactionFromRxnFile(fileName, params).release();
}
//! Parse a text stream in MDL rxn format into a ChemicalReaction
inline ChemicalReaction *RxnDataStreamToChemicalReaction(
    std::istream &rxnStream, unsigned int &line, bool sanitize = false,
    bool removeHs = false, bool strictParsing = true) {
  v2::FileParsers::MolFileParserParams params;
  params.sanitize = sanitize;
  params.removeHs = removeHs;
  params.strictParsing = strictParsing;
  return v2::ReactionParser::ReactionFromRxnDataStream(rxnStream, line, params)
      .release();
}
}  // namespace v1
//! returns an rxn block for a reaction
/*!
   \param rxn            chemical reaction

   \param separateAgents flag to decide if agents are put in a separate block,
                         otherwise they are included in the reactants block
                         (default)

   \param forceV3000     flag to cause the V3000 format to be used instead of
                         V2000
 */
RDKIT_CHEMREACTIONS_EXPORT std::string ChemicalReactionToRxnBlock(
    const ChemicalReaction &rxn, bool separateAgents = false,
    bool forceV3000 = false);
//! returns an V3000 rxn block for a reaction
/*!
   \param rxn            chemical reaction

   \param separateAgents flag to decide if agents are put in a separate block,
                         otherwise they are included in the reactants block
                         (default)
*/
RDKIT_CHEMREACTIONS_EXPORT std::string ChemicalReactionToV3KRxnBlock(
    const ChemicalReaction &rxn, bool separateAgents = false);

//! @}

//---------------------------------------------------------------------------
//! \name PNG Support
//! @{

//! Tags used for PNG metadata
namespace PNGData {
RDKIT_CHEMREACTIONS_EXPORT extern const std::string rxnSmilesTag;
RDKIT_CHEMREACTIONS_EXPORT extern const std::string rxnSmartsTag;
RDKIT_CHEMREACTIONS_EXPORT extern const std::string rxnRxnTag;
RDKIT_CHEMREACTIONS_EXPORT extern const std::string rxnPklTag;
}  // namespace PNGData

namespace v2 {
namespace ReactionParser {

//! \brief constructs a ChemicalReaction from the metadata in a PNG stream
/*!

Looks through the metadata in the PNG to find the first tag that matches one
of the tags in \c RDKit::PNGData. A reaction is constructed from this chunk.

Throws a \c FileParseException if no suitable tag is found.

 */
RDKIT_CHEMREACTIONS_EXPORT std::unique_ptr<ChemicalReaction>
ReactionFromPNGStream(std::istream &pngStream);
//! \brief constructs a ChemicalReaction from the metadata in a PNG string
//! See \c PNGStreamToChemicalReaction() for more details
inline std::unique_ptr<ChemicalReaction> ReactionFromPNGString(
    const std::string &data) {
  std::stringstream inStream(data);
  return ReactionFromPNGStream(inStream);
};
//! \brief constructs a ChemicalReaction from the metadata in a PNG file
//! See \c PNGStreamToChemicalReaction() for more details
inline std::unique_ptr<ChemicalReaction> ReactionFromPNGFile(
    const std::string &fname) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  if (!inStream || (inStream.bad())) {
    throw BadFileException((boost::format("Bad input file %s") % fname).str());
  }
  return ReactionFromPNGStream(inStream);
};
}  // namespace ReactionParser
}  // namespace v2

inline namespace v1 {
//! \brief constructs a ChemicalReaction from the metadata in a PNG stream
/*!

Looks through the metadata in the PNG to find the first tag that matches one
of the tags in \c RDKit::PNGData. A reaction is constructed from this chunk.

Throws a \c FileParseException if no suitable tag is found.

The caller is responsible for the returned pointer.

 */
inline ChemicalReaction *PNGStreamToChemicalReaction(std::istream &pngStream) {
  return v2::ReactionParser::ReactionFromPNGStream(pngStream).release();
}
//! \brief constructs a ChemicalReaction from the metadata in a PNG string
//! See \c PNGStreamToChemicalReaction() for more details
inline ChemicalReaction *PNGStringToChemicalReaction(const std::string &data) {
  return v2::ReactionParser::ReactionFromPNGString(data).release();
}
//! \brief constructs a ChemicalReaction from the metadata in a PNG file
//! See \c PNGStreamToChemicalReaction() for more details
inline ChemicalReaction *PNGFileToChemicalReaction(const std::string &fname) {
  return v2::ReactionParser::ReactionFromPNGFile(fname).release();
}
}  // namespace v1

//! \brief adds metadata for a ChemicalReaction to the data from a PNG stream.
//! The modified PNG data is returned.
/*!

  \param rxn            the reaction to add
  \param iStream        the stream to read from
  \param includePkl     include a reaction pickle
  \param includeSmiles  include reaction SMILES for the reaction
  \param includeSmarts  include reaction SMARTS for the reaction
  \param includeRxn     include an RXN block for the reaction

*/
RDKIT_CHEMREACTIONS_EXPORT std::string addChemicalReactionToPNGStream(
    const ChemicalReaction &rxn, std::istream &iStream, bool includePkl = true,
    bool includeSmiles = true, bool includeSmarts = false,
    bool includeRxn = false);
//! \brief adds metadata for a ChemicalReaction to the data from a PNG string.
//! See addChemicalReactionToPNGStream() for more details.
inline std::string addChemicalReactionToPNGString(const ChemicalReaction &rxn,
                                                  const std::string &pngString,
                                                  bool includePkl = true,
                                                  bool includeSmiles = true,
                                                  bool includeSmarts = false,
                                                  bool includeRxn = false) {
  std::stringstream inStream(pngString);
  return addChemicalReactionToPNGStream(
      rxn, inStream, includePkl, includeSmiles, includeSmarts, includeRxn);
}
//! \brief adds metadata for a ChemicalReaction to the data from a PNG string.
//! See addChemicalReactionToPNGStream() for more details.
inline std::string addChemicalReactionToPNGFile(const ChemicalReaction &rxn,
                                                const std::string &fname,
                                                bool includePkl = true,
                                                bool includeSmiles = true,
                                                bool includeSmarts = false,
                                                bool includeRxn = false) {
  std::ifstream inStream(fname.c_str(), std::ios::binary);
  return addChemicalReactionToPNGStream(
      rxn, inStream, includePkl, includeSmiles, includeSmarts, includeRxn);
}
//! @}

inline std::unique_ptr<ChemicalReaction> operator"" _rxnsmarts(const char *text,
                                                               size_t len) {
  std::string sma(text, len);
  std::unique_ptr<ChemicalReaction> ptr;
  try {
    ptr = v2::ReactionParser::ReactionFromSmarts(sma);
  } catch (...) {
    ptr = nullptr;
  }
  return ptr;
}
inline std::unique_ptr<ChemicalReaction> operator"" _rxnsmiles(const char *text,
                                                               size_t len) {
  std::string sma(text, len);
  std::unique_ptr<ChemicalReaction> ptr;
  try {
    ptr = v2::ReactionParser::ReactionFromSmiles(sma);
  } catch (...) {
    ptr = nullptr;
  }
  return ptr;
}

//---------------------------------------------------------------------------
//! \name CDXML rxn Support
///@{

//! Parse text in CDXML rxn format into a vector of ChemicalReactions
RDKIT_CHEMREACTIONS_EXPORT std::vector<std::unique_ptr<ChemicalReaction>>
CDXMLToChemicalReactions(const std::string &rxnBlock, bool sanitize = false,
                         bool removeHs = false);
//! Parse a file in CDXML rxn format into a vector of ChemicalReactions
RDKIT_CHEMREACTIONS_EXPORT std::vector<std::unique_ptr<ChemicalReaction>>
CDXMLFileToChemicalReactions(const std::string &fileName, bool sanitize = false,
                             bool removeHs = false);
//! Parse a text stream in CDXML rxn format into a vector of ChemicalReactions
RDKIT_CHEMREACTIONS_EXPORT std::vector<std::unique_ptr<ChemicalReaction>>
CDXMLDataStreamToChemicalReactions(std::istream &rxnStream,
                                   bool sanitize = false,
                                   bool removeHs = false);

}  // namespace RDKit
#endif
