/* -------------------------------------------------------------------------
 * Declares tools for Monomeric ROMols
 *
 * A "MonomerMol" is an ROMol that uses RDKit atoms to represent monomers. Chains
 * are represented via the AtomMonomerInfo structs on atoms, and linkages are
 * stord as a LINKAGE property on bonds in the form of RX-RY, where X is the attachment
 * point used on the begin monomer and Y is the attachment point used on the end monomer.
 *
 * Linkages between an atom and a monomer are represented similarly, but with just RX where
 * X is the attachment point on the monomer.
 *
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <RDGeneral/BetterEnums.h>
#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>

namespace RDKit
{

// Forward declarations
class Atom;
class Bond;
class SubstanceGroup;
class MonomerLibrary;

const std::string LINKAGE{"attachmentPoints"};
const std::string EXTRA_LINKAGE{"extraAttachmentPoints"};
const std::string ATOM_LABEL{"atomLabel"};

// Some default linkage options
const std::string BRANCH_LINKAGE{"R3-R1"};
const std::string BACKBONE_LINKAGE{"R2-R1"};
const std::string CROSS_LINKAGE{"R3-R3"};
const std::string HYDROGEN_LINKAGE{"pair-pair"};

// Monomer properties stored on atoms
const std::string BRANCH_MONOMER{"isBranchMonomer"};
const std::string SMILES_MONOMER{"isSmilesMonomer"};

// Substance group property to indicate an annotation on a chain
const std::string ANNOTATION{"annotation"};

// acts as a view to a polymer chain
struct Chain {
    std::vector<unsigned int> atoms;
    std::vector<unsigned int> bonds;
    std::string annotation;
};

// Types of polymer chains
enum class ChainType { PEPTIDE, RNA, DNA, CHEM, OTHER };
enum class MonomerType { REGULAR, SMILES };

// Free utility functions for working with monomer atoms

// Returns true if the atom represents a monomer in a MonomerMol, false otherwise
RDKIT_MONOMERMOL_EXPORT bool isMonomer(const Atom* atom);

// Get the polymer/chain ID for an atom
RDKIT_MONOMERMOL_EXPORT std::string getPolymerId(const Atom* atom);

// Get the residue number for an atom
RDKIT_MONOMERMOL_EXPORT unsigned int getResidueNumber(const Atom* atom);

//! MonomerMol is a molecule class for monomeric representations
/*!
    MonomerMol inherits from ROMol and provides methods for working with
    monomeric molecules where atoms represent monomers and bonds represent
    linkages between monomers.
*/
class RDKIT_MONOMERMOL_EXPORT MonomerMol : public ROMol {
 public:
  // Default constructor
  MonomerMol() : ROMol() {}

  //! Constructor with custom monomer library
  /*!
    \param library  shared pointer to a MonomerLibrary instance
  */
  explicit MonomerMol(std::shared_ptr<MonomerLibrary> library)
      : ROMol(), d_library(std::move(library)) {}

  // Copy constructor from ROMol.
  /*!
    \param other     the molecule to be copied
    \param quickCopy (optional) if this is true, the resulting MonomerMol will not
         copy any of the properties or bookmarks and conformers from \c other.
    \param confId if this is >=0, the resulting MonomerMol will contain only
         the specified conformer from \c other.
  */
  MonomerMol(const ROMol &other, bool quickCopy = false, int confId = -1)
      : ROMol(other, quickCopy, confId) {}
  MonomerMol(const MonomerMol &other);
  MonomerMol &operator=(const MonomerMol &other);
  MonomerMol(MonomerMol &&other) noexcept;
  MonomerMol &operator=(MonomerMol &&other) noexcept;
  MonomerMol(const std::string &binStr); // Pickle constructor (for deserialization)

  //! adds an Atom (not a monomer!) to our collection, effectively making ROMol::addAtom public
  /*!
    \param atom          pointer to the Atom to add
    \param updateLabel   (optional) if this is true, the new Atom will be
                         our \c activeAtom
    \param takeOwnership (optional) if this is true, we take ownership of \c
    atom
                         instead of copying it.

    \return the index of the added atom
  */
  unsigned int addAtom(Atom *atom, bool updateLabel = true,
                       bool takeOwnership = false) {
    return ROMol::addAtom(atom, updateLabel, takeOwnership);
  }

  //! adds a Bond to our collection
  /*!
    \param bond          pointer to the Bond to add
    \param takeOwnership (optional) if this is true, we take ownership of \c
    bond
                         instead of copying it.

    \return the new number of bonds
  */
  unsigned int addBond(Bond *bond, bool takeOwnership = false) {
    return ROMol::addBond(bond, takeOwnership);
  }

  //! \name Bonds
  //! @{

  //! adds a Bond between the indicated Atoms
  /*!
     \return the number of Bonds
  */
  unsigned int addBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                       Bond::BondType order = Bond::UNSPECIFIED);
  //! \overload
  unsigned int addBond(Atom *beginAtom, Atom *endAtom,
                       Bond::BondType order = Bond::UNSPECIFIED);

  // ---- Monomer Operations ----

  /*
   * Add a monomer to the molecule
   *
   * @param name The name of the monomer
   * @param residue_number The residue number of the monomer
   * @param chain_id The chain ID of the monomer, defaults to "A"
   * @param monomer_type The type of monomer to add
   *
   * @return The index of the added monomer
   */
  size_t addMonomer(std::string_view name, int residue_number,
                    std::string_view monomer_class,
                    std::string_view chain_id = "A",
                    MonomerType monomer_type = MonomerType::REGULAR);

  /*
   * Add a monomer to the molecule. Overload that uses the last monomer
   * added to the molecule to determine the chain ID, residue number, and monomer class.
   *
   * @param name The name of the monomer
   * @param monomer_type The type of monomer to add
   *
   * @return The index of the added monomer
   */
  size_t addMonomer(std::string_view name,
                    MonomerType monomer_type = MonomerType::REGULAR);

  // ---- Connection Operations ----

  /*
   * Add a connection between two monomers in the molecule. The connection has
   * directionality that starts at monomer1 and ends at monomer2.
   *
   * @param monomer1 The index of the first monomer
   * @param monomer2 The index of the second monomer
   * @param connection_type The type of connection to add
   */
  void addConnection(size_t monomer1, size_t monomer2,
                     const std::string &linkage);

  /*
   * Add a connection between an atom and a monomer in the molecule.
   * @param atom_idx The index of the atom
   * @param monomer_idx The index of the monomer
   * @param bond_type The type of bond to add
   * @param linkage The linkage information in the form RX where X is the
   * attachment point on the monomer.
  */
  void addAtomMonomerConnection(size_t atom_idx, size_t monomer_idx,
                             const std::string &linkage, Bond::BondType bond_type =
                                 Bond::BondType::SINGLE);

  // ---- Query Operations ----

  [[nodiscard]] Chain getPolymer(std::string_view polymer_id) const;

  [[nodiscard]] std::vector<std::string> getPolymerIds() const;

  // ---- Chain Assignment ----

  // Discards existing chains and reassigns monomers to sequential chains where monomers
  // are reordered based on connectivity.
  void assignChains();

  // ---- Monomer Library Access ----

  //! Get the monomer library for this molecule
  /*!
    Returns the instance library if one has been set, otherwise returns
    the global library.
    \return reference to the MonomerLibrary
  */
  MonomerLibrary& getMonomerLibrary();
  const MonomerLibrary& getMonomerLibrary() const;

  //! Set a custom monomer library for this molecule
  /*!
    \param library  shared pointer to a MonomerLibrary instance,
                    or nullptr to use the global library
  */
  void setMonomerLibrary(std::shared_ptr<MonomerLibrary> library);

  //! Check if this molecule has a custom library set
  [[nodiscard]] bool hasCustomLibrary() const { return d_library != nullptr; }

 private:
  std::shared_ptr<MonomerLibrary> d_library;  // nullptr means use global
};

} // namespace RDKit
