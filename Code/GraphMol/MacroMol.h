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

#include <memory>
#include <string>

#include "MacroAtomInfo.h"
#include "MacroMolTemplate.h"
#include "RWMol.h"

namespace RDKit {

class RDKIT_GRAPHMOL_EXPORT MacroMol : public RWMol {
 public:
  //! Constructs a MacroMol with an empty local template library.
  MacroMol();

  //! Constructs a MacroMol that takes ownership of a local template library.
  /*!
    \param localTemplateLibrary the non-null library to take ownership of
  */
  explicit MacroMol(
      std::unique_ptr<MacroMolTemplateLibrary> localTemplateLibrary);

  MacroMol(const MacroMol &other);
  MacroMol &operator=(const MacroMol &other);
  MacroMol(MacroMol &&other) noexcept = default;
  MacroMol &operator=(MacroMol &&other) noexcept = default;

  //! Replaces the local template library and takes ownership of it.
  /*!
    \param newTemplateLibrary the non-null library to take ownership of
  */
  void setLocalTemplateLibrary(
      std::unique_ptr<MacroMolTemplateLibrary> newTemplateLibrary);

  //! Returns this molecule's local template library.
  MacroMolTemplateLibrary &getLocalTemplateLibrary();
  //! Returns this molecule's local template library.
  const MacroMolTemplateLibrary &getLocalTemplateLibrary() const;

  //! Checks that every macro atom resolves in the local template library.
  /*!
    Macro atoms are looked up by monomer class and symbol, then by monomer
    class and template name. Ordinary atoms are ignored.
  */
  bool checkLocalTemplateReferences() const;

  //! Adds a new macro atom to the molecule.
  /*!
    \param symbol       the symbol (dummy label) used to identify the monomer
    \param monomerClass the class of monomer the macro atom represents

    \return the index of the newly added atom
  */
  unsigned int addMacroAtom(std::string symbol, MonomerClass monomerClass);

  //! Adds a bond between two macro atoms.
  /*!
    Both atoms must be macro atoms.
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
                            Bond::BondType bondType = Bond::SINGLE);

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
      Bond::BondType bondType = Bond::BondType::SINGLE);

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
      int beginAttachPt, Bond::BondType bondType = Bond::BondType::SINGLE);

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
                       Bond::BondType bondType = Bond::SINGLE);
  //! \overload
  unsigned int addBond(Atom *beginAtom, Atom *endAtom,
                       Bond::BondType bondType = Bond::SINGLE);
  //! \overload
  unsigned int addBond(Bond *bond, bool takeOwnership = false);

 private:
  unsigned int addMacroBondHelper(unsigned int beginAtomIdx,
                                  unsigned int endAtomIdx, int beginAttachPt,
                                  int endAttachPt,
                                  Bond::BondType bondType = Bond::SINGLE);

  std::unique_ptr<MacroMolTemplateLibrary> dp_localTemplateLibrary;
};
}  // namespace RDKit

#endif
