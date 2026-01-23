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

void MonomerLibrary::addMonomer(const MonomerDef& def)
{
    auto key = makeHelmKey(def.helmSymbol, def.monomerType);
    auto it = d_helmIndex.find(key);

    if (it != d_helmIndex.end()) {
        // Check if the existing definition matches
        const auto& existing = d_monomers[it->second];
        if (existing.smiles != def.smiles) {
            throw std::runtime_error(
                "Conflicting monomer definition for '" + def.helmSymbol +
                "': existing SMILES '" + existing.smiles +
                "' differs from new SMILES '" + def.smiles + "'");
        }
        // Identical definition, skip silently
        return;
    }

    // Add new monomer
    size_t index = d_monomers.size();
    d_monomers.push_back(def);
    d_helmIndex[key] = index;

    if (!def.pdbCode.empty()) {
        d_pdbIndex[def.pdbCode] = index;
    }
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

std::optional<std::string>
MonomerLibrary::getMonomerSmiles(std::string_view monomerId,
                                  std::string_view monomerType) const
{
    auto key = makeHelmKey(monomerId, monomerType);
    auto it = d_helmIndex.find(key);

    if (it != d_helmIndex.end()) {
        return d_monomers[it->second].smiles;
    }
    return std::nullopt;
}

MonomerLibrary::helm_info_t
MonomerLibrary::getHelmInfo(std::string_view pdbCode) const
{
    auto it = d_pdbIndex.find(std::string(pdbCode));

    if (it != d_pdbIndex.end()) {
        const auto& def = d_monomers[it->second];
        return std::make_tuple(def.helmSymbol, def.smiles, def.monomerType);
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
        lib->addMonomer(def);
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
std::vector<MonomerDef> getDefaultAminoAcids()
{
    return {
        {"A", "ALA",
         "C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CB :1.pdbName. CA "
         ":2.pdbName. N  :3.pdbName. H  :4.pdbName. C  :5.pdbName. O  "
         ":6.pdbName. OXT|",
         "PEPTIDE"},
        {"C", "CYS",
         "O=C([C@H](CS[H:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. SG "
         ":5.pdbName. HG :6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|",
         "PEPTIDE"},
        {"D", "ASP",
         "O=C([C@H](CC(=O)[OH:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
         "OD1:6.pdbName. OD2:7.pdbName. N  :8.pdbName. H  :9.pdbName. OXT|",
         "PEPTIDE"},
        {"E", "GLU",
         "O=C([C@H](CCC(=O)[OH:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
         ":5.pdbName. CD :6.pdbName. OE1:7.pdbName. OE2:8.pdbName. N  "
         ":9.pdbName. H  :10.pdbName. OXT|",
         "PEPTIDE"},
        {"F", "PHE",
         "O=C([C@H](Cc1ccccc1)N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C "
         " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
         "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. CE2:9.pdbName. "
         "CD2:10.pdbName. N  :11.pdbName. H  :12.pdbName. OXT|",
         "PEPTIDE"},
        {"G", "GLY",
         "O=C(CN[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
         ":2.pdbName. CA :3.pdbName. N  :4.pdbName. H  :5.pdbName. OXT|",
         "PEPTIDE"},
        {"H", "HIS",
         "O=C([C@H](Cc1cnc[nH]1)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
         ":5.pdbName. CD2:6.pdbName. NE2:7.pdbName. CE1:8.pdbName. "
         "ND1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|",
         "PEPTIDE"},
        {"I", "ILE",
         "CC[C@H](C)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CD1:1.pdbName. "
         "CG1:2.pdbName. CB :3.pdbName. CG2:4.pdbName. CA :5.pdbName. N  "
         ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
         "PEPTIDE"},
        {"K", "LYS",
         "O=C([C@H](CCCCN[H:3])N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
         ":5.pdbName. CD :6.pdbName. CE :7.pdbName. NZ :8.pdbName. "
         "HZ1:9.pdbName. N  :10.pdbName. H  :11.pdbName. OXT|",
         "PEPTIDE"},
        {"L", "LEU",
         "CC(C)C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CD1:1.pdbName. "
         "CG :2.pdbName. CD2:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
         ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
         "PEPTIDE"},
        {"M", "MET",
         "CSCC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CE :1.pdbName. "
         "SD :2.pdbName. CG :3.pdbName. CB :4.pdbName. CA :5.pdbName. N  "
         ":6.pdbName. H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
         "PEPTIDE"},
        {"N", "ASN",
         "NC(=O)C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. ND2:1.pdbName. CG "
         ":2.pdbName. OD1:3.pdbName. CB :4.pdbName. CA :5.pdbName. N  :6.pdbName. "
         "H  :7.pdbName. C  :8.pdbName. O  :9.pdbName. OXT|",
         "PEPTIDE"},
        {"O", "PYL",
         "C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N[H:1])C(=O)[OH:2] "
         "|atomProp:0.pdbName. CB2:1.pdbName. CG2:2.pdbName. CD2:3.pdbName. "
         "CE2:4.pdbName. N2 :5.pdbName. CA2:6.pdbName. C2 :7.pdbName. O2 "
         ":8.pdbName. NZ :9.pdbName. CE :10.pdbName. CD :11.pdbName. CG "
         ":12.pdbName. CB :13.pdbName. CA :14.pdbName. N  :15.pdbName. H  "
         ":16.pdbName. C  :17.pdbName. O  :18.pdbName. OXT|",
         "PEPTIDE"},
        {"P", "PRO",
         "O=C([C@@H]1CCCN1[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C "
         " :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. CD "
         ":6.pdbName. N  :7.pdbName. H  :8.pdbName. OXT|",
         "PEPTIDE"},
        {"Q", "GLN",
         "NC(=O)CC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. NE2:1.pdbName. CD "
         ":2.pdbName. OE1:3.pdbName. CG :4.pdbName. CB :5.pdbName. CA :6.pdbName. "
         "N  :7.pdbName. H  :8.pdbName. C  :9.pdbName. O  :10.pdbName. OXT|",
         "PEPTIDE"},
        {"R", "ARG",
         "N=C(N)NCCC[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. "
         "NH2:1.pdbName. CZ :2.pdbName. NH1:3.pdbName. NE :4.pdbName. CD "
         ":5.pdbName. CG :6.pdbName. CB :7.pdbName. CA :8.pdbName. N  "
         ":9.pdbName. H  :10.pdbName. C  :11.pdbName. O  :12.pdbName. OXT|",
         "PEPTIDE"},
        {"S", "SER",
         "O=C([C@H](CO)N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
         ":2.pdbName. CA :3.pdbName. CB :4.pdbName. OG :5.pdbName. N  "
         ":6.pdbName. H  :7.pdbName. OXT|",
         "PEPTIDE"},
        {"T", "THR",
         "C[C@@H](O)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. "
         "CG2:1.pdbName. CB :2.pdbName. OG1:3.pdbName. CA :4.pdbName. N  "
         ":5.pdbName. H  :6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|",
         "PEPTIDE"},
        {"U", "SEC",
         "O=C([C@H](C[SeH])N[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. "
         "C  :2.pdbName. CA :3.pdbName. CB :4.pdbName.SE  :5.pdbName. N  "
         ":6.pdbName. H  :7.pdbName. OXT|",
         "PEPTIDE"},
        {"V", "VAL",
         "CC(C)[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CG1:1.pdbName. "
         "CB :2.pdbName. CG2:3.pdbName. CA :4.pdbName. N  :5.pdbName. H  "
         ":6.pdbName. C  :7.pdbName. O  :8.pdbName. OXT|",
         "PEPTIDE"},
        {"W", "TRP",
         "O=C([C@H](Cc1c[nH]c2ccccc12)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG "
         ":5.pdbName. CD1:6.pdbName. NE1:7.pdbName. CE2:8.pdbName. "
         "CZ2:9.pdbName. CH2:10.pdbName. CZ3:11.pdbName. CE3:12.pdbName. "
         "CD2:13.pdbName. N  :14.pdbName. H  :15.pdbName. OXT|",
         "PEPTIDE"},
        {"Y", "TYR",
         "O=C([C@H](Cc1ccc(O)cc1)N[H:1])[OH:2] |atomProp:0.pdbName. O  "
         ":1.pdbName. C  :2.pdbName. CA :3.pdbName. CB :4.pdbName. CG :5.pdbName. "
         "CD1:6.pdbName. CE1:7.pdbName. CZ :8.pdbName. OH :9.pdbName. "
         "CE2:10.pdbName. CD2:11.pdbName. N  :12.pdbName. H  :13.pdbName. OXT|",
         "PEPTIDE"},
        // X for unknown amino acid - mapped to glycine structure
        {"X", "UNK",
         "O=C(CN[H:1])[OH:2] |atomProp:0.pdbName. O  :1.pdbName. C  "
         ":2.pdbName. CA :3.pdbName. N  :4.pdbName. H  :5.pdbName. OXT|",
         "PEPTIDE"},
    };
}

}  // namespace RDKit
