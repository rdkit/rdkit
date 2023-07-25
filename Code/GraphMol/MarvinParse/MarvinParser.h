//
//  Copyright (C) 2002-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MARVINPARSER_H
#define RD_MARVINPARSER_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>

#include <string>
#include <iostream>

namespace RDKit {

RDKIT_FILEPARSERS_EXPORT void *MrvFileParser(const std::string &fname,
                                             bool &isReaction,
                                             bool sanitize = false,
                                             bool removeHs = false);
RDKIT_FILEPARSERS_EXPORT void *MrvDataStreamParser(std::istream *inStream,
                                                   bool &isReaction,
                                                   bool sanitize,
                                                   bool removeHs);

RDKIT_FILEPARSERS_EXPORT void *MrvDataStreamParser(std::istream &inStream,
                                                   bool &isReaction,
                                                   bool sanitize = false,
                                                   bool removeHs = false);
RDKIT_FILEPARSERS_EXPORT void *MrvBlockParser(const std::string &molmrvText,
                                              bool &isReaction,
                                              bool sanitize = false,
                                              bool removeHs = false);

RDKIT_FILEPARSERS_EXPORT RWMol *MrvMolDataStreamParser(std::istream *inStream,
                                                       bool sanitize = false,
                                                       bool removeHs = false);
RDKIT_FILEPARSERS_EXPORT RWMol *MrvMolDataStreamParser(std::istream &inStream,
                                                       bool sanitize = false,
                                                       bool removeHs = false);
RDKIT_FILEPARSERS_EXPORT RWMol *MrvMolStringParser(
    const std::string &molmrvText, bool sanitize = false,
    bool removeHs = false);
RDKIT_FILEPARSERS_EXPORT RWMol *MrvMolFileParser(const std::string &fName,
                                                 bool sanitize = false,
                                                 bool removeHs = false);

RDKIT_FILEPARSERS_EXPORT ChemicalReaction *MrvRxnDataStreamParser(
    std::istream *inStream, bool sanitize = false, bool removeHs = false);
RDKIT_FILEPARSERS_EXPORT ChemicalReaction *MrvRxnDataStreamParser(
    std::istream &inStream, bool sanitize = false, bool removeHs = false);
RDKIT_FILEPARSERS_EXPORT ChemicalReaction *MrvRxnStringParser(
    const std::string &molmrvText, bool sanitize = false,
    bool removeHs = false);
RDKIT_FILEPARSERS_EXPORT ChemicalReaction *MrvRxnFileParser(
    const std::string &fName, bool sanitize = false, bool removeHs = false);

RDKIT_FILEPARSERS_EXPORT void MolToMrvFile(const ROMol &mol,
                                           const std::string &fName,
                                           bool includeStereo, int confId,
                                           bool kekulize,
                                           bool prettyPrint = false);
RDKIT_FILEPARSERS_EXPORT std::string MolToMrvBlock(const ROMol &mol,
                                                   bool includeStereo,
                                                   int confId, bool kekulize,
                                                   bool prettyPrint = false);

RDKIT_FILEPARSERS_EXPORT std::string ChemicalReactionToMrvBlock(
    const ChemicalReaction &rxn, bool prettyPrint = false);

RDKIT_FILEPARSERS_EXPORT void ChemicalReactionToMrvFile(
    const ChemicalReaction &rxn, const std::string &fName,
    bool prettyPrint = false);
}  // namespace RDKit

#endif
