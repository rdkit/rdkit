//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MOLFROMMACROMOL_H
#define RD_MOLFROMMACROMOL_H

#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <RDGeneral/export.h>

#include <memory>

namespace RDKit {

//! Expands a MacroMol into an ordinary molecule.
/*!
  Each macro atom is replaced by a copy of its template, looked up in
  \c templates by monomer class and symbol, then by monomer class and name.
  The leaving groups of the attachment points used by macro bonds are removed,
  and each macro atom's main group is recorded as a SUP SGroup. Ordinary atoms
  and bonds are copied unchanged. The result is not sanitized.

  \param macroMol  the molecule to expand
  \param templates the library used to look up the macro atoms' templates

  \throws ValueErrorException if
    - a macro atom has an empty symbol or no template,
    - an attachment point is used by more than one macro bond,
    - a macro bond uses an attachment point its template does not have, or
    - a macro bond would change the valence of an attachment atom: its bond
      type must match that of the bond to the leaving atom it replaces or, if
      the attachment point has no leaving atom, be HYDROGEN or ZERO
*/
RDKIT_GRAPHMOL_EXPORT std::unique_ptr<RWMol> MolFromMacroMol(
    const MacroMol &macroMol, const MacroMolTemplateLibrary &templates);

}  // namespace RDKit

#endif
