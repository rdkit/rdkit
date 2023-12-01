//
//  Copyright (C) 2022-2023 Tad Hurst, Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>

#ifndef RD_MARVINPARSER_H
#define RD_MARVINPARSER_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemReactions/Reaction.h>

#include <string>
#include <iostream>

namespace RDKit {

RDKIT_MARVINPARSER_EXPORT bool MrvFileIsReaction(const std::string &fname);
RDKIT_MARVINPARSER_EXPORT bool MrvDataStreamIsReaction(std::istream *inStream);

RDKIT_MARVINPARSER_EXPORT bool MrvDataStreamIsReaction(std::istream &inStream);
RDKIT_MARVINPARSER_EXPORT bool MrvBlockIsReaction(
    const std::string &molmrvText);

RDKIT_MARVINPARSER_EXPORT RWMol *MrvDataStreamToMol(std::istream *inStream,
                                                    bool sanitize = false,
                                                    bool removeHs = false);
RDKIT_MARVINPARSER_EXPORT RWMol *MrvDataStreamToMol(std::istream &inStream,
                                                    bool sanitize = false,
                                                    bool removeHs = false);
RDKIT_MARVINPARSER_EXPORT RWMol *MrvBlockToMol(const std::string &molmrvText,
                                               bool sanitize = false,
                                               bool removeHs = false);
RDKIT_MARVINPARSER_EXPORT RWMol *MrvFileToMol(const std::string &fName,
                                              bool sanitize = false,
                                              bool removeHs = false);

RDKIT_MARVINPARSER_EXPORT ChemicalReaction *MrvDataStreamToChemicalReaction(
    std::istream *inStream, bool sanitize = false, bool removeHs = false);
RDKIT_MARVINPARSER_EXPORT ChemicalReaction *MrvDataStreamToChemicalReaction(
    std::istream &inStream, bool sanitize = false, bool removeHs = false);
RDKIT_MARVINPARSER_EXPORT ChemicalReaction *MrvBlockToChemicalReaction(
    const std::string &molmrvText, bool sanitize = false,
    bool removeHs = false);
RDKIT_MARVINPARSER_EXPORT ChemicalReaction *MrvFileToChemicalReaction(
    const std::string &fName, bool sanitize = false, bool removeHs = false);

RDKIT_MARVINPARSER_EXPORT void MolToMrvFile(const ROMol &mol,
                                            const std::string &fName,
                                            bool includeStereo = true,
                                            int confId = -1,
                                            bool kekulize = true,
                                            bool prettyPrint = false);
RDKIT_MARVINPARSER_EXPORT std::string MolToMrvBlock(const ROMol &mol,
                                                    bool includeStereo = true,
                                                    int confId = -1,
                                                    bool kekulize = true,
                                                    bool prettyPrint = false);

RDKIT_MARVINPARSER_EXPORT std::string ChemicalReactionToMrvBlock(
    const ChemicalReaction &rxn, bool prettyPrint = false);

RDKIT_MARVINPARSER_EXPORT void ChemicalReactionToMrvFile(
    const ChemicalReaction &rxn, const std::string &fName,
    bool prettyPrint = false);
}  // namespace RDKit

#endif
