//
//  Copyright (C) 2002-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_FILEWRITERS_H
#define RD_FILEWRITERS_H

#include <RDGeneral/types.h>
#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>

namespace RDKit {

struct RDKIT_FILEPARSERS_EXPORT MolWriterParams {
  bool includeStereo = true;  /**< toggles inclusion of stereochemistry
                                   information*/
  bool kekulize = true;       /**< triggers kekulization of the molecule before
                                   it is written*/
  bool forceV3000 = false;    /**< force generation a V3000 mol block (happens
                                   automatically with more than 999 atoms or
                                   bonds)*/
  unsigned int precision = 6; /**< precision of coordinates (only available in
                                   V3000)*/
};

// \brief generates an MDL mol block for a molecule
/*!
 *   \param mol             - the molecule in question
 *   \param MolWriterParams - parmeter struct with write options
 *   \param confId          - selects the conformer to be used
 *                            (default=-1 - find first in mol)
 */
RDKIT_FILEPARSERS_EXPORT std::string MolToMolBlock(
    const ROMol &mol, const MolWriterParams &params, int confId = -1);

// \brief generates an MDL mol block for a molecule
/*!
 *   \param mol           - the molecule in question
 *   \param includeStereo - toggles inclusion of stereochemistry information
 *                          (default=true)
 *   \param confId        - selects the conformer to be used
 *                          (default=-1 - find first in mol)
 *   \param kekulize      - triggers kekulization
 *                          of the molecule before it is written (default=true)
 *   \param forceV3000    - force generation a V3000 mol block (happens
 *                          automatically with more than 999 atoms or
 *                          bonds)(default=false)
 */
inline std::string MolToMolBlock(const ROMol &mol, bool includeStereo = true,
                                 int confId = -1, bool kekulize = true,
                                 bool forceV3000 = false) {
  MolWriterParams params{includeStereo, kekulize, forceV3000};
  return MolToMolBlock(mol, params, confId);
}

// \brief generates an MDL v3000 mol block for a molecule
/*!
 *   \param mol             - the molecule in question
 *   \param MolWriterParams - parameter struct with write options
 *   \param confId          - selects the conformer to be used
 *                            (default=-1 - find first in mol)
 */
inline std::string MolToV3KMolBlock(const ROMol &mol,
                                    const MolWriterParams &params,
                                    int confId = -1) {
  // have to set forceV300, prefer copy over mutable params argument
  MolWriterParams v3KParams{params};
  v3KParams.forceV3000 = true;
  return MolToMolBlock(mol, v3KParams, confId);
}

// \brief generates an MDL v3000 mol block for a molecule
/*!
 *   \param mol           - the molecule in question
 *   \param includeStereo - toggles inclusion of stereochemistry information
 *   \param confId        - selects the conformer to be used
 *                          (default=-1 - find first in mol)
 *   \param kekulize      - triggers kekulization of the molecule before it is
 *                        - written
 */
inline std::string MolToV3KMolBlock(const ROMol &mol, bool includeStereo = true,
                                    int confId = -1, bool kekulize = true) {
  MolWriterParams params{includeStereo, kekulize, true};
  return MolToMolBlock(mol, params, confId);
}

// \brief generates an MDL v2000 mol block for a molecule
/*!
 *   \param mol             - the molecule in question
 *   \param MolWriterParams - parameter struct with write options
 *   \param confId          - selects the conformer to be used
 *                            (default=-1 - find first in mol)
 *
 *  \note This function will throw a ValueError exception if the molecule has
 *   more than 999 atoms, bonds, or SGroups.
 */
RDKIT_FILEPARSERS_EXPORT std::string MolToV2KMolBlock(
    const ROMol &mol, const MolWriterParams &params = MolWriterParams(),
    int confId = -1);

// \brief Writes a molecule to an MDL mol file
/*!
 *   \param mol             - the molecule in question
 *   \param fName           - the name of the file to use
 *   \param MolWriterParams - parameter struct with write options
 *   \param confId          - selects the conformer to be used
 */
RDKIT_FILEPARSERS_EXPORT void MolToMolFile(const ROMol &mol,
                                           const std::string &fName,
                                           const MolWriterParams &params,
                                           int confId = -1);

// \brief Writes a molecule to an MDL mol file
/*!
 *   \param mol           - the molecule in question
 *   \param fName         - the name of the file to use
 *   \param includeStereo - toggles inclusion of stereochemistry information
 *   \param confId        - selects the conformer to be used
 *   \param kekulize      - triggers kekulization of the molecule before it is
 * written
 *   \param forceV3000    - force generation a V3000 mol block (happens
 * automatically with
 *                          more than 999 atoms or bonds)
 */
inline void MolToMolFile(const ROMol &mol, const std::string &fName,
                         bool includeStereo = true, int confId = -1,
                         bool kekulize = true, bool forceV3000 = false) {
  MolWriterParams params{includeStereo, kekulize, forceV3000};
  MolToMolFile(mol, fName, params, confId);
}

// \brief Writes a molecule to an MDL V3000 mol file
/*!
 *   \param mol             - the molecule in question
 *   \param fName           - the name of the file to use
 *   \param MolWriterParams - parameter struct with write options
 *   \param confId          - selects the conformer to be used
 */
inline void MolToV3KMolFile(const ROMol &mol, const std::string &fName,
                            const MolWriterParams &params, int confId = -1) {
  // have to set forceV300, prefer copy over mutable params argument
  MolWriterParams v3KParams{params};
  v3KParams.forceV3000 = true;
  MolToMolFile(mol, fName, v3KParams, confId);
}

// \brief Writes a molecule to an MDL V3000 mol file
/*!
 *   \param mol           - the molecule in question
 *   \param fName         - the name of the file to use
 *   \param includeStereo - toggles inclusion of stereochemistry information
 *   \param confId        - selects the conformer to be used
 *   \param kekulize      - triggers kekulization of the molecule before it is
 * written
 */
inline void MolToV3KMolFile(const ROMol &mol, const std::string &fName,
                            bool includeStereo = true, int confId = -1,
                            bool kekulize = true) {
  MolWriterParams params{includeStereo, kekulize, true};
  MolToMolFile(mol, fName, params, confId);
}

RDKIT_FILEPARSERS_EXPORT std::string MolToCMLBlock(const ROMol &mol,
                                                   int confId = -1,
                                                   bool kekulize = true);

RDKIT_FILEPARSERS_EXPORT void MolToCMLFile(const ROMol &mol,
                                           const std::string &fName,
                                           int confId = -1,
                                           bool kekulize = true);

// \brief Writes a molecule to an XYZ block
/*!
 *   \param mol       - the molecule in question
 *   \param confId    - selects which conformation to output
 *   \param precision - precision of the coordinates
 */
RDKIT_FILEPARSERS_EXPORT std::string MolToXYZBlock(const ROMol &mol,
                                                   int confId = -1,
                                                   unsigned int precision = 6);

// \brief Writes a molecule to an XYZ block
/*!
 *   \param mol       - the molecule in question
 *   \param fName     - the file to write to
 *   \param confId    - selects which conformation to output
 *   \param precision - precision of the coordinates
 */
RDKIT_FILEPARSERS_EXPORT void MolToXYZFile(const ROMol &mol,
                                           const std::string &fName,
                                           int confId = -1,
                                           unsigned int precision = 6);

RDKIT_FILEPARSERS_EXPORT std::string MolToTPLText(
    const ROMol &mol, const std::string &partialChargeProp = "_GasteigerCharge",
    bool writeFirstConfTwice = false);
RDKIT_FILEPARSERS_EXPORT void MolToTPLFile(
    const ROMol &mol, const std::string &fName,
    const std::string &partialChargeProp = "_GasteigerCharge",
    bool writeFirstConfTwice = false);

// \brief generates an PDB block for a molecule
/*!
 *   \param mol           - the molecule in question
 *   \param confId        - selects the conformer to be used
 *   \param flavor        - controls what gets written:
 *         flavor & 1 : Write MODEL/ENDMDL lines around each record
 *         flavor & 2 : Don't write single CONECT records
 *         flavor & 4 : Write CONECT records in both directions
 *         flavor & 8 : Don't use multiple CONECTs to encode bond order
 *         flavor & 16 : Write MASTER record
 *         flavor & 32 : Write TER record
 */
RDKIT_FILEPARSERS_EXPORT std::string MolToPDBBlock(const ROMol &mol,
                                                   int confId = -1,
                                                   unsigned int flavor = 0);
// \brief Writes a molecule to an MDL mol file
/*!
 *   \param mol           - the molecule in question
 *   \param fName         - the name of the file to use
 *   \param confId        - selects the conformer to be used
 *   \param flavor        - controls what gets written:
 *         flavor & 1 : Write MODEL/ENDMDL lines around each record
 *         flavor & 2 : Don't write single CONECT records
 *         flavor & 4 : Write CONECT records in both directions
 *         flavor & 8 : Don't use multiple CONECTs to encode bond order
 *         flavor & 16 : Write MASTER record
 *         flavor & 32 : Write TER record
 */
RDKIT_FILEPARSERS_EXPORT void MolToPDBFile(const ROMol &mol,
                                           const std::string &fname,
                                           int confId = -1,
                                           unsigned int flavor = 0);

}  // namespace RDKit

#endif
