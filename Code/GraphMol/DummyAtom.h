//
// Copyright (C) 2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_DUMMY_ATOM_H
#define RD_DUMMY_ATOM_H

#include <memory>
#include <optional>
#include <string_view>
#include <GraphMol/Atom.h>

namespace RDKit {

//! The atomic number used for dummy atoms (wildcards, R-groups, attachment
//! points).
constexpr int dummyAtomicNum = 0;

//! The atomLabel property prefix used to identify attachment point dummy atoms.
inline constexpr std::string_view attachmentPointLabelPrefix = "_AP";

//! The atomLabel property prefix used for R-group dummy atoms.
inline constexpr std::string_view rGroupLabelPrefix = "_R";

//! Creates a dummy atom properly configured for SMARTS and MDL MOL file output.
//!
//! The atom has atomic number 0, a null query (so it appears as * rather than
//! [#0] in SMARTS output), and no implicit valence (so no VAL line is written
//! in MDL MOL files).
[[nodiscard]] RDKIT_GRAPHMOL_EXPORT std::unique_ptr<Atom> createDummyAtom();

//! Returns whether the atom is an attachment point dummy.
//!
//! An attachment point dummy has atomic number 0, exactly one neighbor, and an
//! atomLabel property whose value starts with "_AP".
RDKIT_GRAPHMOL_EXPORT bool isAttachmentPointDummy(const Atom &atom);

//! Returns the R-group number of the atom, or an empty optional if the atom is
//! not an R-group atom.
//!
//! \throws std::invalid_argument if the atom's R-group number is 0.
[[nodiscard]] RDKIT_GRAPHMOL_EXPORT std::optional<unsigned int> getRGroupNumber(
    const Atom *atom);

//! Creates a new R-group atom with the specified R-group number.
//!
//! The atom has atomic number 0 and its atomLabel, dummyLabel, _MolFileRLabel,
//! and isotope properties set appropriately for the given R-group number.
//!
//! \throws std::invalid_argument if rGroupNum is 0.
[[nodiscard]] RDKIT_GRAPHMOL_EXPORT std::unique_ptr<Atom> makeNewRGroup(
    unsigned int rGroupNum);

}  // namespace RDKit

#endif
