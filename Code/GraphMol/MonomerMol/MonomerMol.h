/* -------------------------------------------------------------------------
 * Declares tools for Monomeric ROMols
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include <RDGeneral/BetterEnums.h>
#include <RDGeneral/export.h>
#include <GraphMol/RWMol.h>

namespace RDKit
{

// Forward declarations
class Atom;
class Bond;
class SubstanceGroup;

const std::string ANNOTATION{"ANNOTATION"};
const std::string LINKAGE{"attachmentPoints"};
const std::string ATOM_LABEL{"atomLabel"};
const std::string SUPPLEMENTARY_INFORMATION{"SUPPLEMENTARY_INFORMATION"};
const std::string BRANCH_LINKAGE{"R3-R1"};
const std::string BACKBONE_LINKAGE{"R2-R1"};
const std::string CROSS_LINKAGE{"R3-R3"};
const std::string HELM_MODEL{"HELM_MODEL"};
const std::string MONOMER_LIST{"MONOMER_LIST"};
const std::string UNKNOWN_MONOMER{"UNKNOWN_MONOMER"};
const std::string REPETITION_DUMMY_ID{"REPETITION_DUMMY_ID"};
const std::string BRANCH_MONOMER{"isBranchMonomer"};
const std::string SMILES_MONOMER{"isSmilesMonomer"};
const std::string CUSTOM_BOND{"customBond"}; // TODO: remove?

// acts as a view to a polymer chain
struct Chain {
    std::vector<unsigned int> atoms;
    std::vector<unsigned int> bonds;
    std::string annotation;
};

enum class ChainType { PEPTIDE, RNA, DNA, CHEM };
enum class ConnectionType { FORWARD, SIDECHAIN, CROSSLINK };
enum class MonomerType { REGULAR, SMILES };

RDKIT_MONOMERMOL_EXPORT ChainType toChainType(std::string_view chain_type);

RDKIT_MONOMERMOL_EXPORT std::string toString(ChainType chain_type);

//! MonomerMol is an RWMol that represents molecules using both monomers and atoms.
/*!
 * Monomers are represented as atoms with specific properties set, and linkages are
 * stored as a LINKAGE property on bonds in the form of RX-RY, where X is the attachment
 * point used on the begin monomer and Y is the attachment point used on the
 * end monomer.
 */
class RDKIT_MONOMERMOL_EXPORT MonomerMol : public RWMol {
 public:
  //! Default constructor
  MonomerMol() : RWMol() {}

  //! Copy constructor
  MonomerMol(const MonomerMol& other) : RWMol(other) {}

  //! Move constructor
  MonomerMol(MonomerMol&& other) noexcept : RWMol(std::move(other)) {}

  //! Construct from an ROMol
  explicit MonomerMol(const ROMol& other) : RWMol(other) {}

  //! Copy assignment operator
  MonomerMol& operator=(const MonomerMol& other);

  //! Move assignment operator
  MonomerMol& operator=(MonomerMol&& other) noexcept;

  /*!
   * Add a monomer to the molecule
   *
   * @param name The name of the monomer
   * @param residue_number The residue number of the monomer
   * @param chain_id The chain ID of the monomer
   * @param monomer_type The type of monomer to add
   *
   * @return The index of the added monomer
   */
  size_t addMonomer(std::string_view name, int residue_number,
                    std::string_view chain_id,
                    MonomerType monomer_type = MonomerType::REGULAR);

  /*!
   * Add a monomer to the molecule. Overload that uses the last monomer
   * added to the molecule to determine the chain ID and residue number.
   *
   * @param name The name of the monomer
   * @param monomer_type The type of monomer to add
   *
   * @return The index of the added monomer
   */
  size_t addMonomer(std::string_view name,
                    MonomerType monomer_type = MonomerType::REGULAR);

  /*!
   * Add a connection between two monomers in the molecule. The connection has
   * directionality that starts at monomer1 and ends at monomer2.
   *
   * @param monomer1 The index of the first monomer
   * @param monomer2 The index of the second monomer
   * @param connection_type The type of connection to add
   */
  void addConnection(size_t monomer1, size_t monomer2,
                     ConnectionType connection_type = ConnectionType::FORWARD);

  /*!
   * Add a connection between two monomers with a specific linkage string.
   *
   * @param monomer1 The index of the first monomer
   * @param monomer2 The index of the second monomer
   * @param linkage The linkage string in RX-RY format
   * @param is_custom_bond Whether this is a custom bond
   */
  void addConnection(size_t monomer1, size_t monomer2,
                     const std::string& linkage, bool is_custom_bond = false);

  /*!
   * Check if an atom at the given index is a monomer.
   *
   * @param atomIdx The index of the atom to check
   * @return true if the atom is a monomer, false if it's a regular atom
   */
  bool isMonomer(unsigned int atomIdx) const;

  /*!
   * Add a bond between two atoms/monomers with automatic type detection.
   * - If both are regular atoms: creates a normal bond
   * - If both are monomers: throws (use addConnection instead)
   * - If one is a monomer and one is an atom: throws (use overload with
   *   attachment point)
   *
   * @param beginAtomIdx The index of the first atom
   * @param endAtomIdx The index of the second atom
   * @param order The bond type
   * @return The new number of bonds
   */
  unsigned int addBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                       Bond::BondType order = Bond::UNSPECIFIED);

  /*!
   * Add a bond between a monomer and an atom, specifying the attachment point
   * on the monomer.
   *
   * @param monomerIdx The index of the monomer
   * @param atomIdx The index of the atom
   * @param attachmentPoint The attachment point on the monomer (1 for R1, 2 for
   *   R2, etc.)
   * @param order The bond type
   * @return The new number of bonds
   */
  unsigned int addBond(unsigned int monomerIdx, unsigned int atomIdx,
                       unsigned int attachmentPoint,
                       Bond::BondType order = Bond::SINGLE);
};

// Free functions for querying - these work on any molecule with monomer info

[[nodiscard]] RDKIT_MONOMERMOL_EXPORT Chain
getPolymer(const RDKit::ROMol& monomer_mol, std::string_view polymer_id);

[[nodiscard]] RDKIT_MONOMERMOL_EXPORT std::string
getPolymerId(const RDKit::Atom* atom);

[[nodiscard]] RDKIT_MONOMERMOL_EXPORT std::vector<std::string>
getPolymerIds(const RDKit::ROMol& monomer_mol);

[[nodiscard]] RDKIT_MONOMERMOL_EXPORT unsigned int
getResidueNumber(const RDKit::Atom* atom);

//! Discards existing chains and reassigns monomers to sequential chains.
//! (in HELM world, "chains" are called "polymers")
RDKIT_MONOMERMOL_EXPORT void assignChains(MonomerMol& mol);

} // namespace RDKit
