/* -------------------------------------------------------------------------
 * Declares RDKit atomistic ROMol <-> MonomerMol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#pragma once

#include <RDGeneral/export.h>

#include <memory>

namespace RDKit
{
class ROMol;
class RWMol;
class MonomerMol;

/**
 * Converts an atomistic ROMol into a MonomerMol.
 *
 * @param atomistic_mol Atomistic molecule to convert to MonomerMol
 * @return MonomerMol
 */
RDKIT_MONOMERMOL_EXPORT std::unique_ptr<MonomerMol>
toMonomeric(const ROMol& atomistic_mol);

/**
 * Build an atomistic molecule from a MonomerMol using monomers
 * in MonomerLibrary.cpp
 *
 * @param monomer_mol Monomeric molecule to convert to atomistic
 * @return Atomistic molecule
 */
RDKIT_MONOMERMOL_EXPORT std::unique_ptr<RWMol>
toAtomistic(const ROMol& monomer_mol);

} // namespace RDKit
