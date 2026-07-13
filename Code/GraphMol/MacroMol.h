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

class RDKIT_GRAPHMOL_EXPORT MacroMol : public RWMol {
 public:
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
    The graph bond's bond type is unspecified; use MacroBondInfo for the
    macro bond type.

    \param beginAtomIdx  index of the atom where the bond begins
    \param endAtomIdx    index of the atom where the bond ends
    \param beginAttachPt the attachment point on the begin (macro) atom
    \param endAttachPt   the attachment point on the end (macro) atom
    \param bondType      the type of the bond

    \return the new number of bonds
  */
  unsigned int addMacroBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                            int beginAttachPt, int endAttachPt,
                            Bond::BondType bondType = Bond::UNSPECIFIED);

  //! Adds a bond from a regular atom to a macro atom.
  /*!
    The graph bond's bond type is unspecified; use MacroBondInfo for the
    macro bond type.

    \param beginAtomIdx     index of the regular atom where the bond begins
    \param endMacroAtomIdx  index of the macro atom where the bond ends
    \param endAttachPt      the attachment point on the end (macro) atom
    \param bondType         the type of the bond

    \return the new number of bonds
  */
  unsigned int addAtomToMacroAtomBond(
      unsigned int beginAtomIdx, unsigned int endMacroAtomIdx, int endAttachPt,
      Bond::BondType bondType = Bond::BondType::UNSPECIFIED);

  //! Adds a bond from a macro atom to a regular atom.
  /*!
    The graph bond's bond type is unspecified; use MacroBondInfo for the
    macro bond type.

    \param beginMacroAtomIdx index of the macro atom where the bond begins
    \param endAtomIdx        index of the regular atom where the bond ends
    \param beginAttachPt     the attachment point on the begin (macro) atom
    \param bondType          the type of the bond

    \return the new number of bonds
  */
  unsigned int addMacroAtomToAtomBond(
      unsigned int beginMacroAtomIdx, unsigned int endAtomIdx,
      int beginAttachPt, Bond::BondType bondType = Bond::BondType::UNSPECIFIED);

  //! Adds a bond between two regular atoms.
  /*!
    \b Note: neither atom may be a macro atom; use addMacroBond(),
    addAtomToMacroAtomBond(), or addMacroAtomToAtomBond() for bonds involving
    macro atoms.

    \param beginAtomIdx index of the atom where the bond begins
    \param endAtomIdx   index of the atom where the bond ends
    \param bondType     the type of the bond

    \return the new number of bonds
  */
  unsigned int addBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                       Bond::BondType bondType = Bond::UNSPECIFIED);
  //! \overload
  unsigned int addBond(Atom *beginAtom, Atom *endAtom,
                       Bond::BondType bondType = Bond::UNSPECIFIED);
  //! \overload
  unsigned int addBond(Bond *bond, bool takeOwnership = false);

 private:
  unsigned int addMacroBondHelper(unsigned int beginAtomIdx,
                                  unsigned int endAtomIdx, int beginAttachPt,
                                  int endAttachPt,
                                  Bond::BondType bondType = Bond::UNSPECIFIED);
};
}  // namespace RDKit

#endif
