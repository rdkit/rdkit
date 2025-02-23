/* -------------------------------------------------------------------------
 * Declares tools for "monomeristic" ROMols
 *
 * A "monomeristic"  ROMol uses RDKit atoms to represent monomers. Chains
 * are represented by COP Substance Groups on the ROMol.
 *
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <RDGeneral/export.h>

#include <string>
#include <string_view>
#include <vector>

// These are properties set on the atom, bond and molecule level to store HELM
// and monomeristic information.
// NOTE: Most of these won't be fully used until we add a HELM parser
const std::string ANNOTATION{"ANNOTATION"};
const std::string LINKAGE{"attachmentPoints"};
const std::string ATOM_LABEL{"atomLabel"};
const std::string BRANCH_LINKAGE{"R3-R1"};
const std::string BACKBONE_LINKAGE{"R2-R1"};
const std::string BRANCH_MONOMER{"isBranchMonomer"};
const std::string SMILES_MONOMER{"isSmilesMonomer"};


// Forward declarations:
namespace RDKit
{
class ROMol;
class RWMol;
class Atom;
class Bond;
class SubstanceGroup;
} // namespace RDKit

namespace RDKit {

using Monomer = RDKit::Atom;

// acts as a view to a polymer chain, and is only valid as long as the owning
// mol is valid.
struct Chain {
    std::vector<unsigned int> atoms;
    std::vector<unsigned int> bonds;
    std::string annotation;
    // std::string polymer_id;
};

enum class ChainType { PEPTIDE, RNA, DNA, CHEM };
enum class ConnectionType { FORWARD, SIDECHAIN }; /*Forward: R2-R1, Sidechain: R3-R1*/
enum class MonomerType { REGULAR, /* LIST, WILDCARD, */ SMILES };

/*
 * Add a monomer to the molecule
 *
 * @param monomer_mol The CG to add the monomer to
 * @param name The name of the monomer
 * @param residue_number The residue number of the monomer
 * @param chain_id The chain ID of the monomer
 * @param monomer_type The type of monomer to add
 *
 * @return The index of the added monomer
 */
RDKIT_MONOMERMOL_EXPORT size_t add_monomer(
    RDKit::RWMol& monomer_mol, std::string_view name, int residue_number,
    std::string_view chain_id, MonomerType monomer_type = MonomerType::REGULAR);

/*
 * Add a monomer to the molecule. Overload that uses the last monomer
 * added to the molecule to determine the chain ID and residue number.
 *
 * @param monomer_mol The CG to add the monomer to
 * @param name The name of the monomer
 * @param monomer_type The type of monomer to add
 *
 * @return The index of the added monomer
 */
RDKIT_MONOMERMOL_EXPORT size_t
add_monomer(RDKit::RWMol& monomer_mol, std::string_view name,
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
add_connection(RDKit::RWMol& monomer_mol, size_t monomer1, size_t monomer2,
               ConnectionType connection_type = ConnectionType::FORWARD);

// overload to use string to specify the linkage (R3-R1, R2-R1, etc)
RDKIT_MONOMERMOL_EXPORT void add_connection(RDKit::RWMol& monomer_mol, size_t monomer1,
                                         size_t monomer2,
                                         const std::string& linkage);

RDKIT_MONOMERMOL_EXPORT std::vector<std::string> get_polymer_ids(const RDKit::ROMol& monomer_mol);
RDKIT_MONOMERMOL_EXPORT std::string get_polymer_id(const Atom* atom);
RDKIT_MONOMERMOL_EXPORT int get_residue_number(const Atom* atom);
RDKIT_MONOMERMOL_EXPORT Chain get_polymer(const RDKit::ROMol& monomer_mol,
                                       std::string_view polymer_id);

RDKIT_MONOMERMOL_EXPORT ChainType to_chain_type(std::string_view chain_type);
RDKIT_MONOMERMOL_EXPORT std::string to_string(ChainType chain_type);
} // namespace RDKit