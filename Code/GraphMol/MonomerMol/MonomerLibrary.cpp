//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "MonomerLibrary.h"

#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

namespace RDKit {

// Static member definitions
bool MonomerLibrary::s_useGlobalLibrary = true;
std::unique_ptr<MonomerLibrary> MonomerLibrary::s_globalLibrary;
std::once_flag MonomerLibrary::s_globalLibraryOnce;

// Built-in monomer definitions (symbol -> SMILES)
namespace {
const std::unordered_map<std::string, std::string> builtin_monomer_data({
    {"A", "C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CB :1.pdbName. CA "
          ":2.pdbName. N  :3.pdbName. H  :4.pdbName. C  :5.pdbName. O  "
          ":6.pdbName. OXT|"},
    {"C", "O=C([C@H](CS[H:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. SG "
          ":5.pdbName. HG :6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|"},
    {"D",
     "O=C([C@H](CC(=O)[OH:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "OD1:6.pdbName. OD2:7.pdbName. N  :8.pdbName. H  :9.pdbName. OXT|"},
    {"E", "O=C([C@H](CCC(=O)[OH:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD :6.pdbName. OE1:7.pdbName. OE2:8.pdbName. N  "
          ":9.pdbName. H  :10.pdbName. OXT|"},
    {"F",
     "O=C([C@H](Cc1ccccc1)N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C "
     " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. CE2:9.pdbName. "
     "CD2:10.pdbName. N  :11.pdbName. H  :12.pdbName. OXT|"},
    {"G", "O=C(CN[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
          ":2.pdbName. CA :3.pdbName. N  :4.pdbName. H  :5.pdbName. OXT|"},
    {"H", "O=C([C@H](Cc1cnc[nH]1)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD2:6.pdbName. NE2:7.pdbName. CE1:8.pdbName. "
          "ND1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|"},
    {"I",
     "CC[C@H](C)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CD1:1.pdbName. "
     "CG1:2.pdbName. CB :3.pdbName. CG2:4.pdbName. CA :5.pdbName. N  "
     ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"K", "O=C([C@H](CCCCN[H:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD :6.pdbName. CE :7.pdbName. NZ :8.pdbName. "
          "HZ1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|"},
    {"L", "CC(C)C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CD1:1.pdbName. "
          "CG :2.pdbName. CD2:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"M", "CSCC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CE :1.pdbName. "
          "SD :2.pdbName. CG :3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"N",
     "NC(=O)C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. ND2:1.pdbName. CG "
     ":2.pdbName. OD1:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  :6.pdbName. "
     "H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|"},
    {"O", "C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N[H:1])C(=O)[OH:2] "
          "|atomProp:0.pdbName. CB2:1.pdbName. CG2:2.pdbName. CD2:3.pdbName. "
          "CE2:4.pdbName. N2 :5.pdbName. CA2:6.pdbName. C2 :7.pdbName. O2 "
          ":8.pdbName. NZ :9.pdbName. CE :10.pdbName. CD :11.pdbName. CG "
          ":12.pdbName. CB :13.pdbName. CA :14.pdbName. N  :15.pdbName. H  "
          ":16.pdbName. C  :17.pdbName. O  :18.pdbName. OXT|"},
    {"P", "O=C([C@@H]1CCCN1[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C "
          " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD "
          ":6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|"},
    {"Q",
     "NC(=O)CC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. NE2:1.pdbName. CD "
     ":2.pdbName. OE1:3.pdbName. CG :4.pdbName. CB :5.pdbName. CA :6.pdbName. "
     "N  :7.pdbName. H  :8.pdbName. C  :9.pdbName. O  :10.pdbName. OXT|"},
    {"R", "N=C(N)NCCC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. "
          "NH2:1.pdbName. CZ :2.pdbName. NH1:3.pdbName. NE :4.pdbName. CD "
          ":5.pdbName. CG :6.pdbName. CB :7.pdbName. CA :8.pdbName. N  "
          ":9.pdbName. H  :10.pdbName. C  :11.pdbName. O  :12.pdbName. OXT|"},
    {"S", "O=C([C@H](CO)N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
          ":2.pdbName. CA :3.pdbName. CB :4.pdbName. OG :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. OXT|"},
    {"T", "C[C@@H](O)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. "
          "CG2:1.pdbName. CB :2.pdbName. OG1:3.pdbName. CA :4.pdbName. N  "
          ":5.pdbName. H  :6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|"},
    {"U", "O=C([C@H](C[SeH])N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. "
          "C  :2.pdbName. CA :3.pdbName. CB :4.pdbName.SE  :5.pdbName. N  "
          ":6.pdbName. H  :7.pdbName. OXT|"},
    {"V", "CC(C)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CG1:1.pdbName. "
          "CB :2.pdbName. CG2:3.pdbName. CA :4.pdbName. N  :5.pdbName. H  "
          ":6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|"},
    {"W", "O=C([C@H](Cc1c[nH]c2ccccc12)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
          ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
          ":5.pdbName. CD1:6.pdbName. NE1:7.pdbName. CE2:8.pdbName. "
          "CZ2:9.pdbName. CH2:10.pdbName. CZ3:11.pdbName. CE3:12.pdbName. "
          "CD2:13.pdbName. N  :14.pdbName. H  :15.pdbName. OXT|"},
    {"Y",
     "O=C([C@H](Cc1ccc(O)cc1)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
     ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
     "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. OH :9.pdbName. "
     "CE2:10.pdbName. CD2:11.pdbName. N  :12.pdbName. H  :13.pdbName. OXT|"}
});

// PDB three-letter code to single-letter symbol mapping
const std::unordered_map<std::string, std::string> pdb_to_symbol({
    {"ALA", "A"}, // Alanine
    {"ARG", "R"}, // Arginine
    {"ASH", "D"}, // Protonated Aspartic
    {"ASN", "N"}, // Asparagine
    {"ASP", "D"}, // Aspartic acid
    {"CYS", "C"}, // Cysteine
    {"GLN", "Q"}, // Glutamine
    {"GLU", "E"}, // Glutamic
    {"GLY", "G"}, // Glycine
    {"HIS", "H"}, // Histidine
    {"ILE", "I"}, // Isoleucine
    {"LEU", "L"}, // Leucine
    {"LYS", "K"}, // Lysine
    {"MET", "M"}, // Methionine
    {"PHE", "F"}, // Phenylalanine
    {"PRO", "P"}, // Proline
    {"SER", "S"}, // Serine
    {"THR", "T"}, // Threonine
    {"TRP", "W"}, // Tryptophan
    {"TYR", "Y"}, // Tyrosine
    {"VAL", "V"}, // Valine
    {"ARN", "R"}, // Neutral-Arginine
    {"GLH", "E"}, // Protonated Glutamic
    {"HID", "H"}, // Histidine (protonated at delta N)
    {"HIE", "H"}, // Histidine (protonated at epsilon N)
    {"HIP", "H"}, // Histidine (protonated at both N)
    {"HSD", "H"}, // Histidine (protonated at delta N, CHARMM name)
    {"HSE", "H"}, // Histidine (protonated at epsilon N, CHARMM name)
    {"HSP", "H"}, // Histidine (protonated at both N, CHARMM name)
    {"LYN", "K"}, // Protonated Lysine
    {"SRO", "S"}, // Ionized Serine
    {"THO", "T"}, // Ionized Threonine
    {"TYO", "Y"}, // Ionized Tyrosine
    {"SEC", "U"}, // Selenocysteine
    {"PYL", "O"}, // Pyrrolysine
    {"XXX", "X"} // Unknown
});

// Single-letter symbol to canonical PDB code mapping
const std::unordered_map<std::string, std::string> symbol_to_pdb({
    {"A", "ALA"}, // Alanine
    {"R", "ARG"}, // Arginine
    {"N", "ASN"}, // Asparagine
    {"D", "ASP"}, // Aspartic acid
    {"C", "CYS"}, // Cysteine
    {"Q", "GLN"}, // Glutamine
    {"E", "GLU"}, // Glutamic acid
    {"G", "GLY"}, // Glycine
    {"H", "HIS"}, // Histidine
    {"I", "ILE"}, // Isoleucine
    {"L", "LEU"}, // Leucine
    {"K", "LYS"}, // Lysine
    {"M", "MET"}, // Methionine
    {"F", "PHE"}, // Phenylalanine
    {"P", "PRO"}, // Proline
    {"S", "SER"}, // Serine
    {"T", "THR"}, // Threonine
    {"W", "TRP"}, // Tryptophan
    {"Y", "TYR"}, // Tyrosine
    {"V", "VAL"}, // Valine
    {"U", "SEC"}, // Selenocysteine
    {"O", "PYL"}  // Pyrrolysine
});

} // anonymous namespace

// MonomerLibrary implementation
MonomerLibrary::MonomerLibrary(bool loadBuiltins) {
    if (loadBuiltins) {
        loadBuiltinDefinitions();
    }
}

MonomerLibrary::MonomerLibrary([[maybe_unused]] std::string_view database_path) {
    // TODO: Load from database/JSON file
    loadBuiltinDefinitions();
}

MonomerLibrary::~MonomerLibrary() = default;

void MonomerLibrary::loadBuiltinDefinitions() {
    // Load all built-in peptide monomers. default PEPTIDE to begin
    for (const auto& [symbol, data] : builtin_monomer_data) {
        MonomerKey key{symbol, "PEPTIDE"};
        MonomerEntry entry{
            symbol,
            data,
            std::shared_ptr<ROMol>(SmilesToMol(data, 0, false)),
            "PEPTIDE",
            symbol_to_pdb.count(symbol) ? symbol_to_pdb.at(symbol) : ""
        };
        d_monomers[key] = std::move(entry);

        // Add to PDB code index
        if (symbol_to_pdb.count(symbol)) {
            d_pdbCodeIndex[symbol_to_pdb.at(symbol)] = key;
        }
    }

    // Add reverse mappings from PDB codes (for variants like HID, HIE, etc.)
    for (const auto& [pdb_code, symbol] : pdb_to_symbol) {
        if (d_pdbCodeIndex.find(pdb_code) == d_pdbCodeIndex.end()) {
            // Map variant PDB codes to the canonical symbol
            d_pdbCodeIndex[pdb_code] = {symbol, "PEPTIDE"};
        }
    }
}

std::optional<std::string> MonomerLibrary::getMonomerData(
    const std::string& monomer_id,
    const std::string& monomer_class) const {
    MonomerKey key{monomer_id, monomer_class};
    auto it = d_monomers.find(key);
    if (it != d_monomers.end()) {
        return it->second.original_data;
    }
    return std::nullopt;
}

MonomerLibrary::monomer_info_t
MonomerLibrary::getMonomerInfo(const std::string& pdb_code) const {
    auto it = d_pdbCodeIndex.find(pdb_code);
    if (it != d_pdbCodeIndex.end()) {
        const auto& [symbol, monomer_class] = it->second;
        auto monomer_it = d_monomers.find(it->second);
        if (monomer_it != d_monomers.end()) {
            return std::make_tuple(
                symbol,
                monomer_it->second.original_data,
                monomer_class
            );
        }
    }
    return std::nullopt;
}

std::optional<std::string>
MonomerLibrary::getPdbCode(const std::string& symbol,
                           [[maybe_unused]] const std::string& monomer_class) const {
    MonomerKey key{symbol, monomer_class};
    auto it = d_monomers.find(key);
    if (it != d_monomers.end() && !it->second.pdb_code.empty()) {
        return it->second.pdb_code;
    }
    return std::nullopt;
}

std::shared_ptr<ROMol> MonomerLibrary::getMonomer(
    const std::string& monomer_id,
    const std::string& monomer_class) const {
    MonomerKey key{monomer_id, monomer_class};
    auto it = d_monomers.find(key);
    if (it != d_monomers.end()) {
        return it->second.mol;
    }
    return nullptr;
}

void MonomerLibrary::addMonomerFromSmiles(const std::string& smiles,
                                          const std::string& symbol,
                                          const std::string& monomer_class,
                                          const std::string& pdb_code) {
    MonomerKey key{symbol, monomer_class};
    MonomerEntry entry{
        symbol,
        smiles,
        std::shared_ptr<ROMol>(SmilesToMol(smiles, 0, false)),
        monomer_class,
        pdb_code
    };
    d_monomers[key] = std::move(entry);

    if (!pdb_code.empty()) {
        d_pdbCodeIndex[pdb_code] = key;
    }
}

void MonomerLibrary::addMonomerFromSDF(const std::string& sdf_data,
                                       const std::string& symbol,
                                       const std::string& monomer_class,
                                       const std::string& pdb_code) {
    MonomerKey key{symbol, monomer_class};
    MonomerEntry entry{
        symbol,
        sdf_data,
        std::shared_ptr<ROMol>(MolBlockToMol(sdf_data, false, false)),
        monomer_class,
        pdb_code
    };
    d_monomers[key] = std::move(entry);

    if (!pdb_code.empty()) {
        d_pdbCodeIndex[pdb_code] = key;
    }
}

void MonomerLibrary::addMonomer(std::shared_ptr<ROMol> mol,
                                const std::string& symbol,
                                const std::string& monomer_class,
                                const std::string& pdb_code) {
    MonomerKey key{symbol, monomer_class};
    MonomerEntry entry{
        symbol,
        "",  // no original data
        std::move(mol),
        monomer_class,
        pdb_code
    };
    d_monomers[key] = std::move(entry);

    if (!pdb_code.empty()) {
        d_pdbCodeIndex[pdb_code] = key;
    }
}

bool MonomerLibrary::hasMonomer(const std::string& symbol,
                                const std::string& monomer_class) const {
    return d_monomers.find({symbol, monomer_class}) != d_monomers.end();
}

// Static methods for global library management
void MonomerLibrary::useGlobalLibrary(bool use_global) {
    s_useGlobalLibrary = use_global;
}

bool MonomerLibrary::isUsingGlobalLibrary() {
    return s_useGlobalLibrary;
}

MonomerLibrary& MonomerLibrary::getGlobalLibrary() {
    std::call_once(s_globalLibraryOnce, []() {
        s_globalLibrary = std::make_unique<MonomerLibrary>(true);  // global always has built-ins
    });
    return *s_globalLibrary;
}

} // namespace RDKit
