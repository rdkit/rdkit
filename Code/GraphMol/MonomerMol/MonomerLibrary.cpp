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

#include <stdexcept>

#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {

// Static member initialization
std::shared_ptr<MonomerLibrary> MonomerLibrary::s_globalLibrary = nullptr;
std::mutex MonomerLibrary::s_globalMutex;

std::string MonomerLibrary::makeHelmKey(std::string_view symbol,
                                         std::string_view monomerType)
{
    return std::string(symbol) + ":" + std::string(monomerType);
}

MonomerLibrary::MonomerLibrary() = default;

MonomerLibrary::MonomerLibrary(MonomerLibrary&& other) noexcept
    : d_monomers(std::move(other.d_monomers)),
      d_helmIndex(std::move(other.d_helmIndex)),
      d_pdbIndex(std::move(other.d_pdbIndex))
{
}

MonomerLibrary& MonomerLibrary::operator=(MonomerLibrary&& other) noexcept
{
    if (this != &other) {
        d_monomers = std::move(other.d_monomers);
        d_helmIndex = std::move(other.d_helmIndex);
        d_pdbIndex = std::move(other.d_pdbIndex);
    }
    return *this;
}

MonomerLibrary::~MonomerLibrary() = default;

void MonomerLibrary::addMonomerDef(MonomerDef def)
{
    auto key = makeHelmKey(def.helmSymbol, def.monomerType);
    auto it = d_helmIndex.find(key);

    if (it != d_helmIndex.end()) {
        // Check if the existing definition matches by comparing molecules
        const auto& existing = d_monomers[it->second];

        // Compare by canonical SMILES if both have molecules
        if (existing.mol && def.mol) {
            auto existingSmiles = MolToSmiles(*existing.mol);
            auto newSmiles = MolToSmiles(*def.mol);
            if (existingSmiles != newSmiles) {
                throw std::runtime_error(
                    "Conflicting monomer definition for '" + def.helmSymbol +
                    "': existing structure differs from new structure");
            }
        }
        // Identical definition, skip silently
        return;
    }

    // Add new monomer
    size_t index = d_monomers.size();
    d_monomers.push_back(std::move(def));
    d_helmIndex[key] = index;

    if (!d_monomers.back().pdbCode.empty()) {
        d_pdbIndex[d_monomers.back().pdbCode] = index;
    }
}

void MonomerLibrary::addMonomerFromSmiles(std::string_view helmSymbol,
                                           std::string_view monomerType,
                                           std::string_view smiles,
                                           std::string_view pdbCode)
{
    // Parse the SMILES without sanitization to preserve attachment points
    bool sanitize = false;
    std::unique_ptr<RWMol> mol(SmilesToMol(std::string(smiles), 0, sanitize));

    if (!mol) {
        throw std::runtime_error("Failed to parse SMILES for monomer '" +
                                 std::string(helmSymbol) + "': " +
                                 std::string(smiles));
    }

    MonomerDef def;
    def.helmSymbol = std::string(helmSymbol);
    def.pdbCode = std::string(pdbCode);
    def.monomerType = std::string(monomerType);
    def.sourceFormat = "SMILES";
    def.sourceData = std::string(smiles);
    def.mol = std::move(mol);

    addMonomerDef(std::move(def));
}

void MonomerLibrary::addMonomerFromMol(std::string_view helmSymbol,
                                        std::string_view monomerType,
                                        std::shared_ptr<ROMol> mol,
                                        std::string_view pdbCode)
{
    if (!mol) {
        throw std::runtime_error("Cannot add null molecule for monomer '" +
                                 std::string(helmSymbol) + "'");
    }

    MonomerDef def;
    def.helmSymbol = std::string(helmSymbol);
    def.pdbCode = std::string(pdbCode);
    def.monomerType = std::string(monomerType);
    def.sourceFormat = "MOL";
    def.sourceData = "";  // No text source for molecules added directly
    def.mol = std::move(mol);

    addMonomerDef(std::move(def));
}

void MonomerLibrary::addPdbAlias(std::string_view pdbCode,
                                  std::string_view helmSymbol,
                                  std::string_view monomerType)
{
    auto key = makeHelmKey(helmSymbol, monomerType);
    auto it = d_helmIndex.find(key);

    if (it == d_helmIndex.end()) {
        throw std::runtime_error("Cannot add PDB alias '" +
                                 std::string(pdbCode) +
                                 "': HELM symbol '" + std::string(helmSymbol) +
                                 "' does not exist");
    }

    d_pdbIndex[std::string(pdbCode)] = it->second;
}

bool MonomerLibrary::hasMonomer(std::string_view helmSymbol,
                                 std::string_view monomerType) const
{
    auto key = makeHelmKey(helmSymbol, monomerType);
    return d_helmIndex.find(key) != d_helmIndex.end();
}

std::unique_ptr<RWMol>
MonomerLibrary::getMonomerMol(std::string_view helmSymbol,
                               std::string_view monomerType) const
{
    auto key = makeHelmKey(helmSymbol, monomerType);
    auto it = d_helmIndex.find(key);

    if (it != d_helmIndex.end() && d_monomers[it->second].mol) {
        // Return a copy of the molecule
        return std::make_unique<RWMol>(*d_monomers[it->second].mol);
    }
    return nullptr;
}

std::optional<std::string>
MonomerLibrary::getMonomerSourceData(std::string_view helmSymbol,
                                      std::string_view monomerType) const
{
    auto key = makeHelmKey(helmSymbol, monomerType);
    auto it = d_helmIndex.find(key);

    if (it != d_helmIndex.end()) {
        return d_monomers[it->second].sourceData;
    }
    return std::nullopt;
}

std::optional<std::string>
MonomerLibrary::getMonomerSmiles(std::string_view monomerId,
                                  std::string_view monomerType) const
{
    auto key = makeHelmKey(monomerId, monomerType);
    auto it = d_helmIndex.find(key);

    if (it != d_helmIndex.end()) {
        const auto& def = d_monomers[it->second];
        // If originally from SMILES, return the original
        if (def.sourceFormat == "SMILES" && !def.sourceData.empty()) {
            return def.sourceData;
        }
        // Otherwise generate SMILES from the molecule
        if (def.mol) {
            return MolToSmiles(*def.mol);
        }
    }
    return std::nullopt;
}

MonomerLibrary::helm_info_t
MonomerLibrary::getHelmInfo(std::string_view pdbCode) const
{
    auto it = d_pdbIndex.find(std::string(pdbCode));

    if (it != d_pdbIndex.end()) {
        const auto& def = d_monomers[it->second];
        std::string smiles;
        if (def.sourceFormat == "SMILES" && !def.sourceData.empty()) {
            smiles = def.sourceData;
        } else if (def.mol) {
            smiles = MolToSmiles(*def.mol);
        }
        return std::make_tuple(def.helmSymbol, smiles, def.monomerType);
    }
    return std::nullopt;
}

std::optional<std::string>
MonomerLibrary::getPdbCode(std::string_view helmSymbol,
                            std::string_view monomerType) const
{
    auto key = makeHelmKey(helmSymbol, monomerType);
    auto it = d_helmIndex.find(key);

    if (it != d_helmIndex.end()) {
        const auto& pdbCode = d_monomers[it->second].pdbCode;
        if (!pdbCode.empty()) {
            return pdbCode;
        }
    }
    return std::nullopt;
}

const MonomerDef*
MonomerLibrary::getMonomerDef(std::string_view helmSymbol,
                               std::string_view monomerType) const
{
    auto key = makeHelmKey(helmSymbol, monomerType);
    auto it = d_helmIndex.find(key);

    if (it != d_helmIndex.end()) {
        return &d_monomers[it->second];
    }
    return nullptr;
}

size_t MonomerLibrary::size() const
{
    return d_monomers.size();
}

bool MonomerLibrary::empty() const
{
    return d_monomers.empty();
}

void MonomerLibrary::clear()
{
    d_monomers.clear();
    d_helmIndex.clear();
    d_pdbIndex.clear();
}

std::shared_ptr<MonomerLibrary> MonomerLibrary::createWithDefaults()
{
    auto lib = std::make_shared<MonomerLibrary>();

    for (const auto& def : getDefaultAminoAcids()) {
        lib->addMonomerFromSmiles(def.helmSymbol, def.monomerType,
                                   def.sourceData, def.pdbCode);
    }

    // Add PDB code aliases for protonation variants
    lib->addPdbAlias("ARN", "R", "PEPTIDE");  // Neutral-Arginine
    lib->addPdbAlias("ASH", "D", "PEPTIDE");  // Protonated Aspartic
    lib->addPdbAlias("GLH", "E", "PEPTIDE");  // Protonated Glutamic
    lib->addPdbAlias("HID", "H", "PEPTIDE");  // His (delta protonated)
    lib->addPdbAlias("HIE", "H", "PEPTIDE");  // His (epsilon protonated)
    lib->addPdbAlias("HIP", "H", "PEPTIDE");  // His (both protonated)
    lib->addPdbAlias("HSD", "H", "PEPTIDE");  // CHARMM delta
    lib->addPdbAlias("HSE", "H", "PEPTIDE");  // CHARMM epsilon
    lib->addPdbAlias("HSP", "H", "PEPTIDE");  // CHARMM both
    lib->addPdbAlias("LYN", "K", "PEPTIDE");  // Protonated Lysine
    lib->addPdbAlias("SRO", "S", "PEPTIDE");  // Ionized Serine
    lib->addPdbAlias("THO", "T", "PEPTIDE");  // Ionized Threonine
    lib->addPdbAlias("TYO", "Y", "PEPTIDE");  // Ionized Tyrosine
    lib->addPdbAlias("XXX", "X", "PEPTIDE");  // Unknown

    return lib;
}

std::shared_ptr<MonomerLibrary> MonomerLibrary::getGlobalLibrary()
{
    std::lock_guard<std::mutex> lock(s_globalMutex);

    if (!s_globalLibrary) {
        s_globalLibrary = createWithDefaults();
    }
    return s_globalLibrary;
}

void MonomerLibrary::setGlobalLibrary(std::shared_ptr<MonomerLibrary> library)
{
    std::lock_guard<std::mutex> lock(s_globalMutex);

    if (library) {
        s_globalLibrary = std::move(library);
    } else {
        // Reset to defaults
        s_globalLibrary = createWithDefaults();
    }
}

// Default amino acid definitions
// Note: These return MonomerDef with sourceData populated for SMILES parsing
std::vector<MonomerDef> getDefaultAminoAcids()
{
    return {
        {"A", "ALA", "PEPTIDE", "SMILES",
         "C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CB :1.pdbName. CA "
         ":2.pdbName. N  :3.pdbName. H  :4.pdbName. C  :5.pdbName. O  "
         ":6.pdbName. OXT|",
         nullptr},
        {"C", "CYS", "PEPTIDE", "SMILES",
         "O=C([C@H](CS[H:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. SG "
         ":5.pdbName. HG :6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|",
         nullptr},
        {"D", "ASP", "PEPTIDE", "SMILES",
         "O=C([C@H](CC(=O)[OH:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
         "OD1:6.pdbName. OD2:7.pdbName. N  :8.pdbName. H  :9.pdbName. OXT|",
         nullptr},
        {"E", "GLU", "PEPTIDE", "SMILES",
         "O=C([C@H](CCC(=O)[OH:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
         ":5.pdbName. CD :6.pdbName. OE1:7.pdbName. OE2:8.pdbName. N  "
         ":9.pdbName. H  :10.pdbName. OXT|",
         nullptr},
        {"F", "PHE", "PEPTIDE", "SMILES",
         "O=C([C@H](Cc1ccccc1)N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C "
         " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
         "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. CE2:9.pdbName. "
         "CD2:10.pdbName. N  :11.pdbName. H  :12.pdbName. OXT|",
         nullptr},
        {"G", "GLY", "PEPTIDE", "SMILES",
         "O=C(CN[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
         ":2.pdbName. CA :3.pdbName. N  :4.pdbName. H  :5.pdbName. OXT|",
         nullptr},
        {"H", "HIS", "PEPTIDE", "SMILES",
         "O=C([C@H](Cc1cnc[nH]1)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
         ":5.pdbName. CD2:6.pdbName. NE2:7.pdbName. CE1:8.pdbName. "
         "ND1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|",
         nullptr},
        {"I", "ILE", "PEPTIDE", "SMILES",
         "CC[C@H](C)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CD1:1.pdbName. "
         "CG1:2.pdbName. CB :3.pdbName. CG2:4.pdbName. CA :5.pdbName. N  "
         ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
         nullptr},
        {"K", "LYS", "PEPTIDE", "SMILES",
         "O=C([C@H](CCCCN[H:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
         ":5.pdbName. CD :6.pdbName. CE :7.pdbName. NZ :8.pdbName. "
         "HZ1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|",
         nullptr},
        {"L", "LEU", "PEPTIDE", "SMILES",
         "CC(C)C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CD1:1.pdbName. "
         "CG :2.pdbName. CD2:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
         ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
         nullptr},
        {"M", "MET", "PEPTIDE", "SMILES",
         "CSCC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CE :1.pdbName. "
         "SD :2.pdbName. CG :3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
         ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
         nullptr},
        {"N", "ASN", "PEPTIDE", "SMILES",
         "NC(=O)C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. ND2:1.pdbName. CG "
         ":2.pdbName. OD1:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  :6.pdbName. "
         "H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
         nullptr},
        {"O", "PYL", "PEPTIDE", "SMILES",
         "C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N[H:1])C(=O)[OH:2] "
         "|atomProp:0.pdbName. CB2:1.pdbName. CG2:2.pdbName. CD2:3.pdbName. "
         "CE2:4.pdbName. N2 :5.pdbName. CA2:6.pdbName. C2 :7.pdbName. O2 "
         ":8.pdbName. NZ :9.pdbName. CE :10.pdbName. CD :11.pdbName. CG "
         ":12.pdbName. CB :13.pdbName. CA :14.pdbName. N  :15.pdbName. H  "
         ":16.pdbName. C  :17.pdbName. O  :18.pdbName. OXT|",
         nullptr},
        {"P", "PRO", "PEPTIDE", "SMILES",
         "O=C([C@@H]1CCCN1[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C "
         " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD "
         ":6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|",
         nullptr},
        {"Q", "GLN", "PEPTIDE", "SMILES",
         "NC(=O)CC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. NE2:1.pdbName. CD "
         ":2.pdbName. OE1:3.pdbName. CG :4.pdbName. CB :5.pdbName. CA :6.pdbName. "
         "N  :7.pdbName. H  :8.pdbName. C  :9.pdbName. O  :10.pdbName. OXT|",
         nullptr},
        {"R", "ARG", "PEPTIDE", "SMILES",
         "N=C(N)NCCC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. "
         "NH2:1.pdbName. CZ :2.pdbName. NH1:3.pdbName. NE :4.pdbName. CD "
         ":5.pdbName. CG :6.pdbName. CB :7.pdbName. CA :8.pdbName. N  "
         ":9.pdbName. H  :10.pdbName. C  :11.pdbName. O  :12.pdbName. OXT|",
         nullptr},
        {"S", "SER", "PEPTIDE", "SMILES",
         "O=C([C@H](CO)N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
         ":2.pdbName. CA :3.pdbName. CB :4.pdbName. OG :5.pdbName. N  "
         ":6.pdbName. H  :7.pdbName. OXT|",
         nullptr},
        {"T", "THR", "PEPTIDE", "SMILES",
         "C[C@@H](O)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. "
         "CG2:1.pdbName. CB :2.pdbName. OG1:3.pdbName. CA :4.pdbName. N  "
         ":5.pdbName. H  :6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|",
         nullptr},
        {"U", "SEC", "PEPTIDE", "SMILES",
         "O=C([C@H](C[SeH])N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. "
         "C  :2.pdbName. CA :3.pdbName. CB :4.pdbName.SE  :5.pdbName. N  "
         ":6.pdbName. H  :7.pdbName. OXT|",
         nullptr},
        {"V", "VAL", "PEPTIDE", "SMILES",
         "CC(C)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CG1:1.pdbName. "
         "CB :2.pdbName. CG2:3.pdbName. CA :4.pdbName. N  :5.pdbName. H  "
         ":6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|",
         nullptr},
        {"W", "TRP", "PEPTIDE", "SMILES",
         "O=C([C@H](Cc1c[nH]c2ccccc12)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
         ":5.pdbName. CD1:6.pdbName. NE1:7.pdbName. CE2:8.pdbName. "
         "CZ2:9.pdbName. CH2:10.pdbName. CZ3:11.pdbName. CE3:12.pdbName. "
         "CD2:13.pdbName. N  :14.pdbName. H  :15.pdbName. OXT|",
         nullptr},
        {"Y", "TYR", "PEPTIDE", "SMILES",
         "O=C([C@H](Cc1ccc(O)cc1)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
         "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. OH :9.pdbName. "
         "CE2:10.pdbName. CD2:11.pdbName. N  :12.pdbName. H  :13.pdbName. OXT|",
         nullptr},
        // X for unknown amino acid - mapped to glycine structure
        {"X", "UNK", "PEPTIDE", "SMILES",
         "O=C(CN[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
         ":2.pdbName. CA :3.pdbName. N  :4.pdbName. H  :5.pdbName. OXT|",
         nullptr},
    };
}

}  // namespace RDKit
