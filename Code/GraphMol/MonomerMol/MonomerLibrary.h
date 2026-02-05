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

#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>

#include <RDGeneral/export.h>

namespace RDKit {

class ROMol;

//! Represents a monomer definition in the library
struct RDKIT_MONOMERMOL_EXPORT MonomerEntry {
    std::string symbol;           // e.g., "A" for Alanine
    std::string original_data;    // Original definition (SMILES, SDF, etc.)
    std::shared_ptr<ROMol> mol;   // Parsed molecule
    std::string monomer_class;    // "PEPTIDE", "RNA", etc.
    std::string pdb_code;         // e.g., "ALA" (optional)
};

//! Library of monomer definitions for use with MonomerMol
/*!
    MonomerLibrary provides access to monomer definitions. It supports two modes:
    1. Global mode (default): A singleton with global monomer definitions
    2. Instance mode: Per-MonomerMol instance with local definitions

    Access is typically via MonomerMol::getMonomerLibrary().
*/
class RDKIT_MONOMERMOL_EXPORT MonomerLibrary {
 public:
    //! Return type for getMonomerInfo: {symbol, original_data, monomer_class}
    using monomer_info_t = std::optional<std::tuple<std::string, std::string, std::string>>;

    //! Default constructor - loads built-in definitions
    MonomerLibrary();

    //! Constructor with path to database/json file (future use)
    explicit MonomerLibrary(std::string_view database_path);

    ~MonomerLibrary();

    // Non-copyable but movable
    MonomerLibrary(const MonomerLibrary&) = delete;
    MonomerLibrary& operator=(const MonomerLibrary&) = delete;
    MonomerLibrary(MonomerLibrary&&) = default;
    MonomerLibrary& operator=(MonomerLibrary&&) = default;

    // --- Query operations ---

    //! Get original data (SMILES/SDF/etc) for a monomer by its symbol and class
    [[nodiscard]] std::optional<std::string> getMonomerData(
        const std::string& monomer_id,
        const std::string& monomer_class) const;

    //! Get full monomer info from a three-letter PDB code
    /*!
        \param pdb_code The three-letter PDB residue code (e.g., "ALA")
        \return tuple of {symbol, original_data, monomer_class} or nullopt if not found
    */
    [[nodiscard]] monomer_info_t
    getMonomerInfo(const std::string& pdb_code) const;

    //! Get PDB code for a monomer symbol
    [[nodiscard]] std::optional<std::string>
    getPdbCode(const std::string& symbol, const std::string& monomer_class) const;

    //! Get parsed molecule for a monomer
    [[nodiscard]] std::shared_ptr<ROMol> getMonomerMol(
        const std::string& monomer_id,
        const std::string& monomer_class) const;

    // --- Mutation operations (for instance libraries) ---

    //! Add a monomer from a SMILES string
    void addMonomerFromSmiles(const std::string& smiles,
                              const std::string& symbol,
                              const std::string& monomer_class,
                              const std::string& pdb_code = "");

    //! Add a monomer from an SDF/molblock string
    void addMonomerFromSDF(const std::string& sdf_data,
                           const std::string& symbol,
                           const std::string& monomer_class,
                           const std::string& pdb_code = "");

    //! Add a monomer with a pre-parsed molecule
    void addMonomer(std::shared_ptr<ROMol> mol,
                    const std::string& symbol,
                    const std::string& monomer_class,
                    const std::string& pdb_code = "");

    //! Check if a monomer exists
    [[nodiscard]] bool hasMonomer(const std::string& symbol,
                                  const std::string& monomer_class) const;

    // --- Global library configuration ---

    //! Enable/disable global library mode (default: true)
    static void useGlobalLibrary(bool use_global);

    //! Check if global library mode is enabled
    static bool isUsingGlobalLibrary();

    //! Get the global singleton instance
    static MonomerLibrary& getGlobalLibrary();

 private:
    //! Key is (symbol, monomer_class)
    using MonomerKey = std::pair<std::string, std::string>;

    //! Hash function for MonomerKey
    struct MonomerKeyHash {
        std::size_t operator()(const MonomerKey& key) const {
            std::size_t h1 = std::hash<std::string>{}(key.first);
            std::size_t h2 = std::hash<std::string>{}(key.second);
            return h1 ^ (h2 << 1);
        }
    };

    std::unordered_map<MonomerKey, MonomerEntry, MonomerKeyHash> d_monomers;

    //! Index from PDB code to (symbol, monomer_class)
    std::unordered_map<std::string, MonomerKey> d_pdbCodeIndex;

    //! Load built-in monomer definitions
    void loadBuiltinDefinitions();

    //! Static state for global mode
    static bool s_useGlobalLibrary;
    static std::unique_ptr<MonomerLibrary> s_globalLibrary;
    static std::once_flag s_globalLibraryOnce;
};

} // namespace RDKit
