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

namespace RDKit {

//! Represents a single monomer definition
struct RDKIT_MONOMERMOL_EXPORT MonomerDef {
    std::string helmSymbol;   // Single letter code (e.g., "A")
    std::string pdbCode;      // 3-letter PDB code (e.g., "ALA") -- optional
    std::string smiles;       // SMILES with attachment points OR SDF with attachment points
    std::string monomerType;  // "PEPTIDE", "RNA", "DNA", "CHEM", free string
};

//! Library of monomer definitions for converting between representations
/*!
 * MonomerLibrary stores monomer definitions that map between:
 * - HELM single-letter symbols (e.g., "A" for Alanine)
 * - PDB 3-letter codes (e.g., "ALA") -- optional
 * - SMILES with attachment point annotations
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

    //! Add a single monomer definition
    /*!
     * @param def The monomer definition to add
     * @throws std::runtime_error if a monomer with the same HELM symbol and
     *         monomer type already exists with a different SMILES
     *
     * If an identical definition already exists, this is a no-op.
     */
    void addMonomer(const MonomerDef& def);

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

    //! Get SMILES for a monomer by HELM symbol
    [[nodiscard]] std::optional<std::string>
    getMonomerSmiles(std::string_view monomerId,
                     std::string_view monomerType) const;

    //! Get HELM info (symbol, SMILES, type) from PDB code
    [[nodiscard]] helm_info_t getHelmInfo(std::string_view pdbCode) const;

    //! Get PDB code from HELM symbol
    [[nodiscard]] std::optional<std::string>
    getPdbCode(std::string_view helmSymbol, std::string_view monomerType) const;

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

    // Static global library instance
    static std::shared_ptr<MonomerLibrary> s_globalLibrary;
    static std::mutex s_globalMutex;
};

//! Get default amino acid monomer definitions
RDKIT_MONOMERMOL_EXPORT std::vector<MonomerDef> getDefaultAminoAcids();

}  // namespace RDKit
