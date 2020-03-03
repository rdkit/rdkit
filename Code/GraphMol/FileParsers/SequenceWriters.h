//
//  Copyright (C) 2015 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_SEQUENCEWRITE_H_
#define _RD_SEQUENCEWRITE_H_
#include <string>

namespace RDKit {
class ROMol;

// \brief construct a sequence string from a molecule (currently only supports
// peptides)
/*!
 *   \param mol - the molecule to work with
 *
 *   \note \c mol should contain monomer information in \c AtomMonomerInfo
 *structures
 */
RDKIT_FILEPARSERS_EXPORT std::string MolToSequence(const ROMol &mol);
// \brief construct a FASTA string from a molecule (currently only supports
// peptides)
/*!
 *   \param mol - the molecule to work with
 *
 *   \note \c mol should contain monomer information in \c AtomMonomerInfo
 *structures
 */
RDKIT_FILEPARSERS_EXPORT std::string MolToFASTA(const ROMol &mol);
// \brief construct a HELM string from a molecule (currently only supports
// peptides)
/*!
 *   \param mol - the molecule to work with
 *
 *   \note \c mol should contain monomer information in \c AtomMonomerInfo
 *structures
 */
RDKIT_FILEPARSERS_EXPORT std::string MolToHELM(const ROMol &mol);
}  // namespace RDKit

#endif
