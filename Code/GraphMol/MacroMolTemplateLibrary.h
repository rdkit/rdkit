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

//! Returns the process-wide library of built-in MacroMol templates.
RDKIT_GRAPHMOL_EXPORT MacroMolTemplateLibrary &
getGlobalMacroMolTemplateLibrary();

//! Adds copies of all built-in MacroMol templates to a caller-owned library.
RDKIT_GRAPHMOL_EXPORT void addBuiltinMacroMolTemplates(
    MacroMolTemplateLibrary &templates);

}  // namespace RDKit

#endif
