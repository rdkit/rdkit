#include "MonomerDatabase.h"

#include <array>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>

#include "boost/algorithm/string/trim.hpp"
#include "boost/noncopyable.hpp"

#include "MonomerMol.h" // ChainType

namespace RDKit
{

const std::unordered_map<std::string, std::string> monomer_to_smiles(
    {{"A", "C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CB :1.pdbName. CA "
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

const std::unordered_map<std::string, std::string> three_character_codes({
    {"ALA", "A"}, // Alanine
    {"ARG", "R"}, // Arginine
    {"ASH", "D"}, // Protonated Aspartic
    {"ASN", "N"}, // Asparagine
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
    {"ASH", "D"}, // Protonated Aspartic
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
    {"XXX", "X"} // Unknown
});

MonomerDatabase::MonomerDatabase([[maybe_unused]] std::string_view database_path)
{
    // TODO: integration with database
    return;
}

MonomerDatabase::~MonomerDatabase()
{
    // TODO: close connection to database
    return;
}

[[nodiscard]] MonomerDatabase::monomer_smiles_t
MonomerDatabase::getMonomerSmiles(std::string monomer_id,
                                        [[maybe_unused]] ChainType monomer_type)
{
    // currently everything is a peptide
    if (monomer_to_smiles.find(monomer_id) != monomer_to_smiles.end()) {
        return monomer_to_smiles.at(monomer_id);
    }
    return std::nullopt;
}

[[nodiscard]] MonomerDatabase::helm_info_t
MonomerDatabase::getHelmInfo(const std::string& pdb_code)
{
    // currently everything is a peptide
    if (three_character_codes.find(pdb_code) != three_character_codes.end()) {
        return std::make_pair(three_character_codes.at(pdb_code), ChainType::PEPTIDE);
    }
    return std::nullopt;
}

} // namespace RDKit
