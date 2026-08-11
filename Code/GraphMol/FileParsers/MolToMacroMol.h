//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
#ifndef RD_MOLTOMACROMOL_H
#define RD_MOLTOMACROMOL_H

#include <GraphMol/MacroMol.h>
#include <GraphMol/MacroMolTemplate.h>
#include <GraphMol/RWMol.h>
#include <RDGeneral/export.h>

#include <memory>

namespace RDKit {

//! Annotates template matches on a complete atomistic molecule.
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<RWMol> MolToMacroMol(
    const ROMol &mol, const MacroMolTemplateLibrary &templates);

//! Collapses MacroMol SubstanceGroups into macro atoms and macro bonds.
RDKIT_FILEPARSERS_EXPORT std::unique_ptr<MacroMol> CollapseMacroMol(
    const ROMol &groupedMol, const MacroMolTemplateLibrary &templates);

}  // namespace RDKit

#endif
