//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MACROMOL_H
#define RD_MACROMOL_H

#include <string>

#include "MacroAtomInfo.h"
#include "RWMol.h"

namespace RDKit {

//! Classes of monomer that a macro atom can represent.
/*!
  Supported monomer classes for macro atoms.

  A regular enum is used here because MonomerClass appears in public MacroMol
  method signatures. RDKit's BETTER_ENUM macro can expand to different C++
  types depending on per-source-file preprocessor settings, which can cause
  link errors in shared-library and unity builds. String conversion is kept in
  MacroMol.cpp so the public header stays stable.
*/
enum MonomerClass : int { AA, NA, CHEM, OTHER };

class RDKIT_GRAPHMOL_EXPORT MacroMol : public RWMol {
 public:
  using RWMol::addBond;

  //! Adds a new macro atom to the molecule.
  /*!
    \param monomerClass the class of monomer the macro atom represents
    \param symbol       the symbol (dummy label) used to identify the monomer

    \return the index of the newly added atom
  */
  unsigned int addMacroAtom(MonomerClass monomerClass, std::string symbol);

  //! Adds a bond between two macro atoms.
  /*!
    At least one of the two atoms must be a macro atom.

    \param beginAtomIdx  index of the atom where the bond begins
    \param endAtomIdx    index of the atom where the bond ends
    \param beginAttachPt the attachment point on the begin (macro) atom
    \param endAttachPt   the attachment point on the end (macro) atom
    \param bondType      the type of the bond

    \return the index of the newly added bond
  */
  unsigned int addMacroBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                            int beginAttachPt, int endAttachPt,
                            Bond::BondType bondType = Bond::BondType::SINGLE);

  //! Adds a bond from a regular atom to a macro atom.
  /*!
    \param beginAtomIdx     index of the regular atom where the bond begins
    \param endMacroAtomIdx  index of the macro atom where the bond ends
    \param endAttachPt      the attachment point on the end (macro) atom
    \param bondType         the type of the bond

    \return the index of the newly added bond
  */
  unsigned int addAtomToMacroAtomBond(
      unsigned int beginAtomIdx, unsigned int endMacroAtomIdx, int endAttachPt,
      Bond::BondType bondType = Bond::BondType::SINGLE);

  //! Adds a bond from a macro atom to a regular atom.
  /*!
    \param beginMacroAtomIdx index of the macro atom where the bond begins
    \param endAtomIdx        index of the regular atom where the bond ends
    \param beginAttachPt     the attachment point on the begin (macro) atom
    \param bondType          the type of the bond

    \return the index of the newly added bond
  */
  unsigned int addMacroAtomToAtomBond(
      unsigned int beginMacroAtomIdx, unsigned int endAtomIdx,
      int beginAttachPt, Bond::BondType bondType = Bond::BondType::SINGLE);

  //! Adds a bond between two regular atoms.
  /*!
    \b Note: neither atom may be a macro atom; use addMacroBond(),
    addAtomToMacroAtomBond(), or addMacroAtomToAtomBond() for bonds involving
    macro atoms.

    \param beginAtomIdx index of the atom where the bond begins
    \param endAtomIdx   index of the atom where the bond ends
    \param bondType     the type of the bond

    \return the index of the newly added bond
  */
  unsigned int addBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                       Bond::BondType bondType = Bond::BondType::SINGLE);
};
}  // namespace RDKit

#endif
