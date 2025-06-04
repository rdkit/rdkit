/* -------------------------------------------------------------------------
 * Declares tools for Monomeristic ROMols
 *
 * A "MonomerMol" is an ROMol that uses RDKit atoms to represent monomers. Chains
 * are represented via the PDBAtomResidueInfo structs on atoms, and linkages are
 * stord as a LINKAGE property on bonds in the form of RX-RY, where X is the attachment
 * point used on the begin monomer and Y is the attachment point used on the monomer.
 *
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <string>
#include <string_view>

#include <RDGeneral/export.h>

namespace RDKit
{

// Forward declarations
class ROMol;
class RWMol;
class Atom;
class Bond;
class SubstanceGroup;

const std::string ANNOTATION{"ANNOTATION"};
const std::string LINKAGE{"attachmentPoints"};
const std::string ATOM_LABEL{"atomLabel"};
const std::string SUPPLEMENTARY_INFORMATION{"SUPPLEMENTARY_INFORMATION"};
const std::string BRANCH_LINKAGE{"R3-R1"};
const std::string BACKBONE_LINKAGE{"R2-R1"};
const std::string HELM_MODEL{"HELM_MODEL"};
const std::string MONOMER_LIST{"MONOMER_LIST"};
const std::string UNKNOWN_MONOMER{"UNKNOWN_MONOMER"};
const std::string REPETITION_DUMMY_ID{"REPETITION_DUMMY_ID"};

// NOTE: These are to allow replacement of the python api
const std::string BRANCH_MONOMER{"isBranchMonomer"};
const std::string SMILES_MONOMER{"isSmilesMonomer"};
const std::string CUSTOM_BOND{"customBond"};
enum class ChainType { PEPTIDE, RNA, DNA, CHEM };
enum class ConnectionType { FORWARD, SIDECHAIN };
enum class MonomerType { REGULAR, /* LIST, WILDCARD, */ SMILES };

RDKIT_MONOMERMOL_EXPORT ChainType toChainType(std::string_view chain_type);

RDKIT_MONOMERMOL_EXPORT std::string toString(ChainType chain_type);

/*
 * Add a monomer to the molecule
 *
 * @param monomer_mol The monomeristic molecule to add the monomer to
 * @param name The name of the monomer
 * @param residue_number The residue number of the monomer
 * @param chain_id The chain ID of the monomer
 * @param monomer_type The type of monomer to add
 *
 * @return The index of the added monomer
 */
RDKIT_MONOMERMOL_EXPORT size_t addMonomer(
    RDKit::RWMol& monomer_mol, std::string_view name, int residue_number,
    std::string_view chain_id, MonomerType monomer_type = MonomerType::REGULAR);

/*
 * Add a monomer to the molecule. Overload that uses the last monomer
 * added to the molecule to determine the chain ID and residue number.
 *
 * @param monomer_mol The monomeristic molecule to add the monomer to
 * @param name The name of the monomer
 * @param monomer_type The type of monomer to add
 *
 * @return The index of the added monomer
 */
RDKIT_MONOMERMOL_EXPORT size_t
addMonomer(RDKit::RWMol& monomer_mol, std::string_view name,
           MonomerType monomer_type = MonomerType::REGULAR);

/*
 * Add a connection between two monomers in the molecule. The connection has
 * directionality that starts at monomer1 and ends at monomer2.
 *
 * @param mol The molecule to add the connection to
 * @param monomer1 The index of the first monomer
 * @param monomer2 The index of the second monomer
 * @param connection_type The type of connection to add
 */
RDKIT_MONOMERMOL_EXPORT void
addConnection(RDKit::RWMol& mol, size_t monomer1, size_t monomer2,
              ConnectionType connection_type = ConnectionType::FORWARD);

// overload for helm writer
RDKIT_MONOMERMOL_EXPORT void addConnection(RDKit::RWMol& mol, size_t monomer1,
                                        size_t monomer2,
                                        const std::string& linkage,
                                        const bool is_custom_bond = false);

RDKIT_MONOMERMOL_EXPORT [[nodiscard]] std::string
getPolymerId(const RDKit::Atom* atom);

RDKIT_MONOMERMOL_EXPORT [[nodiscard]] std::vector<std::string>
getPolymerIds(const RDKit::ROMol& monomer_mol);

RDKIT_MONOMERMOL_EXPORT [[nodiscard]] unsigned int
getResidueNumber(RDKit::Atom* atom);

// Discards existing chains and reassigns monomers to sequential chains.
// (in HELM world, "chains" are called "polymers")
RDKIT_MONOMERMOL_EXPORT void assignChains(RDKit::RWMol& mol);

} // namespace RDKit
