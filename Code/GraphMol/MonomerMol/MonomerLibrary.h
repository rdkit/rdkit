//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#pragma once

#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>

namespace RDKit {

//! Represents a single monomer definition
struct RDKIT_MONOMERMOL_EXPORT MonomerDef {
    std::string helmSymbol;       // Single letter code (e.g., "A")
    std::string pdbCode;          // 3-letter PDB code (e.g., "ALA") -- optional
    std::string monomerType;      // "PEPTIDE", "RNA", "DNA", "CHEM", free string
    std::string sourceFormat;     // Original format: "SMILES", "SDF", "MOL", etc.
    std::string sourceData;       // Original text data (SMILES string, SDF block, etc.)
    std::shared_ptr<ROMol> mol;   // Parsed RDKit molecule with attachment points
};

//! Library of monomer definitions for converting between representations
/*!
 * MonomerLibrary stores monomer definitions that map between:
 * - HELM single-letter symbols (e.g., "A" for Alanine)
 * - PDB 3-letter codes (e.g., "ALA") -- optional
 * - SMILES with attachment point annotations
 * - RDKit molecules with attachment point annotations
 *
 * Three usage patterns are supported:
 * 1. Per-molecule: Create a library during SCSR parsing, attach to MonomerMol
 * 2. Shared: Create once, pass to multiple MolFromSCSR calls
 * 3. Global: Loaded from file/environment, used as fallback
 */
class RDKIT_MONOMERMOL_EXPORT MonomerLibrary {
  public:
    using helm_info_t =
        std::optional<std::tuple<std::string, std::string, std::string>>;

    //! Construct an empty library
    MonomerLibrary();

    //! Move constructor
    MonomerLibrary(MonomerLibrary&& other) noexcept;

    //! Move assignment
    MonomerLibrary& operator=(MonomerLibrary&& other) noexcept;

    //! Deleted copy constructor and assignment
    MonomerLibrary(const MonomerLibrary&) = delete;
    MonomerLibrary& operator=(const MonomerLibrary&) = delete;

    ~MonomerLibrary();

    //! Add a monomer from SMILES
    /*!
     * @param helmSymbol The HELM symbol (e.g., "A")
     * @param monomerType The monomer type (e.g., "PEPTIDE")
     * @param smiles SMILES string with attachment points
     * @param pdbCode Optional PDB 3-letter code
     * @throws std::runtime_error if SMILES cannot be parsed or if a conflicting
     *         definition already exists
     */
    void addMonomerFromSmiles(std::string_view helmSymbol,
                              std::string_view monomerType,
                              std::string_view smiles,
                              std::string_view pdbCode = "");

    //! Add a monomer from an RDKit molecule
    /*!
     * @param helmSymbol The HELM symbol (e.g., "A")
     * @param monomerType The monomer type (e.g., "PEPTIDE")
     * @param mol The RDKit molecule with attachment points marked
     * @param pdbCode Optional PDB 3-letter code
     * @throws std::runtime_error if a conflicting definition already exists
     */
    void addMonomerFromMol(std::string_view helmSymbol,
                           std::string_view monomerType,
                           std::shared_ptr<ROMol> mol,
                           std::string_view pdbCode = "");

    //! Add a PDB code alias for an existing monomer
    /*!
     * Maps a PDB code to an existing HELM symbol. Used for protonation
     * variants like HID, HIE, HIP -> H.
     *
     * @param pdbCode The PDB code to add (e.g., "HID")
     * @param helmSymbol The HELM symbol it maps to (e.g., "H")
     * @param monomerType The monomer type (e.g., "PEPTIDE")
     * @throws std::runtime_error if the helmSymbol doesn't exist
     */
    void addPdbAlias(std::string_view pdbCode, std::string_view helmSymbol,
                     std::string_view monomerType);

    //! Check if a monomer is defined
    [[nodiscard]] bool hasMonomer(std::string_view helmSymbol,
                                  std::string_view monomerType) const;

    //! Get the RDKit molecule for a monomer
    /*!
     * Returns a copy of the monomer's molecule template.
     * @param helmSymbol The HELM symbol
     * @param monomerType The monomer type
     * @return A new RWMol copy of the monomer, or nullptr if not found
     */
    [[nodiscard]] std::unique_ptr<RWMol>
    getMonomerMol(std::string_view helmSymbol,
                  std::string_view monomerType) const;

    //! Get the original source data for a monomer (SMILES, SDF, etc.)
    [[nodiscard]] std::optional<std::string>
    getMonomerSourceData(std::string_view helmSymbol,
                         std::string_view monomerType) const;

    //! Get SMILES for a monomer by HELM symbol (for backwards compatibility)
    /*!
     * If the monomer was added from SMILES, returns the original SMILES.
     * Otherwise, generates SMILES from the stored molecule.
     */
    [[nodiscard]] std::optional<std::string>
    getMonomerSmiles(std::string_view monomerId,
                     std::string_view monomerType) const;

    //! Get HELM info (symbol, SMILES, type) from PDB code
    [[nodiscard]] helm_info_t getHelmInfo(std::string_view pdbCode) const;

    //! Get PDB code from HELM symbol
    [[nodiscard]] std::optional<std::string>
    getPdbCode(std::string_view helmSymbol, std::string_view monomerType) const;

    //! Get the full monomer definition
    [[nodiscard]] const MonomerDef*
    getMonomerDef(std::string_view helmSymbol,
                  std::string_view monomerType) const;

    //! Get the number of monomers in the library
    [[nodiscard]] size_t size() const;

    //! Check if library is empty
    [[nodiscard]] bool empty() const;

    //! Clear all monomer definitions
    void clear();

    // Static factory methods

    //! Create a library with default amino acid definitions
    [[nodiscard]] static std::shared_ptr<MonomerLibrary> createWithDefaults();

    // Global library management

    //! Get the global library (lazy-initialized with defaults)
    [[nodiscard]] static std::shared_ptr<MonomerLibrary> getGlobalLibrary();

    //! Set the global library (pass nullptr to reset to defaults)
    static void setGlobalLibrary(std::shared_ptr<MonomerLibrary> library);

  private:
    std::vector<MonomerDef> d_monomers;

    // Index: (helm_symbol + ":" + monomer_type) -> index in d_monomers
    std::unordered_map<std::string, size_t> d_helmIndex;

    // Index: pdb_code -> index in d_monomers
    std::unordered_map<std::string, size_t> d_pdbIndex;

    // Helper to create lookup key
    static std::string makeHelmKey(std::string_view symbol,
                                   std::string_view monomerType);

    // Helper to add a monomer definition (used by addMonomerFromSmiles/Mol)
    void addMonomerDef(MonomerDef def);

    // Static global library instance
    static std::shared_ptr<MonomerLibrary> s_globalLibrary;
    static std::mutex s_globalMutex;
};

//! Get default amino acid monomer definitions
RDKIT_MONOMERMOL_EXPORT std::vector<MonomerDef> getDefaultAminoAcids();

}  // namespace RDKit
