//
//  Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//
#ifndef RD_MACROMOLTEMPLATELIBRARY_H
#define RD_MACROMOLTEMPLATELIBRARY_H

#include <GraphMol/MacroMolTemplate.h>
#include <RDGeneral/export.h>

namespace RDKit {

// The built-in MacroMol templates are read from the JSON file named by the
// RDKIT_MONOMER_LIBRARY environment variable if it is set, and otherwise from
// default definitions compiled into the library. An error is thrown if the
// file cannot be opened or contains invalid definitions.

//! Returns the process-wide library of built-in MacroMol templates.
/*!
  The library is read-only; use addBuiltinMacroMolTemplates() to populate a
  caller-owned library that can be extended.
*/
RDKIT_MACROMOLTEMPLATELIBRARY_EXPORT const MacroMolTemplateLibrary &
getGlobalMacroMolTemplateLibrary();

//! Adds copies of all built-in MacroMol templates to a caller-owned library.
RDKIT_MACROMOLTEMPLATELIBRARY_EXPORT void addBuiltinMacroMolTemplates(
    MacroMolTemplateLibrary &templates);

}  // namespace RDKit

#endif
