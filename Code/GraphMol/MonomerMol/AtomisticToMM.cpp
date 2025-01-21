#include "Conversions.h"

#include <queue>
#include <unordered_map>
#include <memory>

#include <GraphMol/RWMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "MonomerDatabase.h"
#include "MonomerMol.h"

namespace RDKit
{

namespace
{

const std::string ATTACH_FROM{"attachFrom"};
const std::string MONOMER_IDX1{"monomerIndex1"};
const std::string MONOMER_IDX2{"monomerIndex2"};
const std::string REFERENCE_IDX{"referenceIndex"};

constexpr int SIDECHAIN_IDX = 2;
constexpr int MIN_ATTCHPTS = 2;
constexpr auto NO_ATTACHMENT = std::numeric_limits<unsigned int>::max();

// 3-letter to 1-letter amino acid code mapping
// From mmpdb_get_three_to_one_letter_residue_map
// These are not currently in the monomer database, but some have
// symbols that are already in the monomer database. We may need to
// figure out how to have multiple 3-letter codes for a single symbol
// polymer_type pair (Histidine is the best example)
const std::unordered_map<std::string, char> backup_res_table({
    {"ARN", 'R'}, // Neutral-Arginine
    {"ASH", 'D'}, // Protonated Aspartic
    {"GLH", 'E'}, // Protonated Glutamic
    {"HID", 'H'}, // Histidine (protonated at delta N)
    {"HIE", 'H'}, // Histidine (protonated at epsilon N)
    {"HIP", 'H'}, // Histidine (protonated at both N)
    {"HSD", 'H'}, // Histidine (protonated at delta N, CHARMM name)
    {"HSE", 'H'}, // Histidine (protonated at epsilon N, CHARMM name)
    {"HSP", 'H'}, // Histidine (protonated at both N, CHARMM name)
    {"LYN", 'K'}, // Protonated Lysine
    {"SRO", 'S'}, // Ionized Serine
    {"THO", 'T'}, // Ionized Threonine
    {"TYO", 'Y'}, // Ionized Tyrosine
    {"XXX", 'X'}, // Unknown
});

struct Linkage {
    unsigned int monomer_idx1;
    unsigned int monomer_idx2;
    unsigned int attach_from;
    unsigned int attach_to;
    std::string to_string() const
    {
        return "R" + std::to_string(attach_from) + "-R" + std::to_string(attach_to);
    }
};

void remove_waters(RDKit::RWMol& mol)
{
    mol.beginBatchEdit();
    auto is_water = [](const RDKit::Atom* atom) {
        const auto res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        if (res_info && res_info->getResidueName() == "HOH ") {
            return true;
        }
        return false;
    };
    for (auto atom : mol.atoms()) {
        if (is_water(atom)) {
            mol.removeAtom(atom);
        }
    }
    mol.commitBatchEdit();
}

void find_chains_and_residues(
    const RDKit::ROMol& mol,
    std::map<std::string, std::map<int, std::vector<unsigned int>>>&
        chains_and_residues)
{
    // Find all chains and residues in the molecule
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        const auto atom = mol.getAtomWithIdx(i);
        const auto res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        const auto chain_id = res_info->getChainId();
        const auto res_num = res_info->getResidueNumber();
        chains_and_residues[chain_id][res_num].push_back(i);
    }
}

boost::shared_ptr<RDKit::RWMol>
annotated_atomistic_to_monomeristic(const RDKit::ROMol& input_mol)
{
    // Make RWMol and remove waters
    RDKit::RWMol mol(input_mol);
    remove_waters(mol);

    // Set reference index for SMILES fragments
    for (auto at : mol.atoms()) {
        at->setProp(REFERENCE_IDX, at->getIdx());
    }

    // Map chain_id -> {residue mols}
    std::map<std::string, std::map<int, std::vector<unsigned int>>>
        chains_and_residues;
    find_chains_and_residues(mol, chains_and_residues);

    // Monomer database connection to verify monomers and get HELM info
    // TODO: connect to actual database
    MonomerDatabase db("");
    std::map<ChainType, unsigned int> chain_counts = {{ChainType::PEPTIDE, 0},
                                                      {ChainType::RNA, 0},
                                                      {ChainType::DNA, 0},
                                                      {ChainType::CHEM, 0}};
    boost::shared_ptr<RDKit::RWMol> monomer_mol = boost::make_shared<RDKit::RWMol>();
    for (const auto& [chain_id, residues] : chains_and_residues) {
        // Use first residue to determine chain type. We assume that PDB data
        // is correct and there aren't multiple chain types in a single chain.
        // TODO: Actually check for this
        // What if the first residue is unknown?
        // note: residues are 1-indexed
        const auto first_res_info =
            dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                mol.getAtomWithIdx(residues.begin()->second[0])
                    ->getMonomerInfo());
        auto res_name = first_res_info->getResidueName();
        res_name.erase(
            std::remove_if(res_name.begin(), res_name.end(), ::isspace),
            res_name.end());

        // Default chain type is PEPTIDE if not specified.
        auto helm_info = db.get_helm_info(res_name);
        auto chain_type = helm_info ? helm_info->second : ChainType::PEPTIDE;
        std::string helm_chain_id = to_string(chain_type) + std::to_string(++chain_counts[chain_type]);
        int last_monomer = -1;
        size_t count = 0; // not always the same as res_num, sometimes res_num
                          // doesn't start at 1
        for (const auto& [res_num, atom_idxs] : residues) {
            ++count;
            const auto res_info =
                dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                    mol.getAtomWithIdx(atom_idxs[0])->getMonomerInfo());
            res_name = res_info->getResidueName();
            res_name.erase(
                std::remove_if(res_name.begin(), res_name.end(), ::isspace),
                res_name.end());

            // Determine whether every residue with the number has the same PDB
            // code
            bool same_code = true;
            for (const auto& atom_idx : atom_idxs) {
                const auto atom_res_info =
                    dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                        mol.getAtomWithIdx(atom_idx)->getMonomerInfo());
                if (atom_res_info->getResidueName() != res_name) {
                    same_code = false;
                    break;
                }
            }

            helm_info = db.get_helm_info(res_name);
            size_t this_monomer;
            if (helm_info) {
                // Standard residue in monomer DB
                // TODO: verify SMILES is correct with database
                this_monomer = add_monomer(*monomer_mol, helm_info->first, res_num,
                                           helm_chain_id);
            } else if (!same_code && backup_res_table.find(res_name) !=
                                         backup_res_table.end()) {
                // Standard residue not in monomer DB. 1-letter code is stored
                // through map
                this_monomer = add_monomer(
                    *monomer_mol, std::string{backup_res_table.at(res_name)},
                    res_num, helm_chain_id);
            } else {
                // TODO: add support for smiles monomer (we have internal support, need to move more code)
                throw std::runtime_error("All residues must have a PDB code (no SMILES monomer for now)");
            }
            if (last_monomer != -1) {
                // Add linkage. For now we assume all linkages are backbone
                // linkages and there are no cycles
                add_connection(*monomer_mol, last_monomer, this_monomer,
                               BACKBONE_LINKAGE);
            }
            last_monomer = this_monomer;
        }
    }

    return monomer_mol;
}

bool has_residue_info(const RDKit::ROMol& mol)
{
    return std::any_of(mol.atoms().begin(), mol.atoms().end(),
                       [](const auto& atom) {
                           return atom->getMonomerInfo() &&
                                  atom->getMonomerInfo()->getMonomerType() ==
                                      RDKit::AtomMonomerInfo::PDBRESIDUE;
                       });
}

} // unnamed namespace

boost::shared_ptr<RDKit::RWMol> atomisticToMonomerMol(const RDKit::ROMol& mol)
{

    if (!has_residue_info(mol)) {
        throw std::runtime_error(
            "No residue information found in molecule");
    }
    auto monomer_mol = annotated_atomistic_to_monomeristic(mol);
    assign_chains(*monomer_mol);
    return monomer_mol;
}

} // namespace RDKit