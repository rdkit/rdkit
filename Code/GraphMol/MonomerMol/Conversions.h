/* -------------------------------------------------------------------------
 * Declares RDKit atomistic ROMol <-> MonomerMol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#pragma once

#include <RDGeneral/export.h>

#include <boost/shared_ptr.hpp>

namespace RDKit
{
class ROMol;
class RWMol;

/**
 * Converts an atomistic ROMol into a Monomer-based ROMol that can
 * be written to HELM format using the HELM writer.
 *
 * @param atomistic_mol Atomistic molecule to convert to MonomerMol
 * @param use_residue_info Whether to use PDBAtomResidueInfo to determine
 * monomer boundaries
 * @return MonomerMol
 */
RDKIT_MONOMERMOL_EXPORT boost::shared_ptr<RDKit::RWMol>
atomisticToMonomerMol(const RDKit::ROMol& atomistic_mol, bool use_residue_info = true);

/**
 * Build an atomistic molecule from a MonomerMol
 *
 * @param monomer_mol Monomeristic molecule to convert to atomistic
 * @return Atomistic molecule
 */
RDKIT_MONOMERMOL_EXPORT boost::shared_ptr<RDKit::RWMol>
monomerMolToAtomsitic(const RDKit::ROMol& monomer_mol);

} // namespace RDKit