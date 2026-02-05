//
//  Copyright (C) 2025 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Conversions.h"

#include <queue>
#include <span>
#include <unordered_map>
#include <memory>

#include <GraphMol/RWMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "MonomerLibrary.h"
#include "MonomerMol.h"

namespace RDKit
{

namespace
{

static constexpr const char* MONOMER_IDX{"monomerIndex"};
static constexpr const char* MONOMER_MAP_NUM{"monomerMapNumber"};
static constexpr const char* REFERENCE_IDX{"referenceIndex"};

static constexpr int MIN_ATTCHPTS = 2;

// std::map to allow sequential/ordered iteration
using ChainsAndResidues =
    std::map<std::string,
             std::map<std::pair<int, std::string>, std::vector<unsigned int>>>;

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

void neutralizeAtoms(RDKit::ROMol& mol)
{
    // Algorithm for neutralizing molecules from
    // https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules by Noel
    // O’Boyle Will neutralize the molecule by adding or removing hydrogens as
    // needed. This will ensure SMILES can be used to match atomistic structures
    // to the correct monomer.
    static const std::unique_ptr<RDKit::RWMol> neutralize_query(
        RDKit::SmartsToMol(
            "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"));
    for (const auto& match : RDKit::SubstructMatch(mol, *neutralize_query)) {
        auto atom = mol.getAtomWithIdx(match[0].second);
        auto chg = atom->getFormalCharge();
        auto hcount = atom->getTotalNumHs();
        atom->setFormalCharge(0);
        atom->setNumExplicitHs(hcount - chg);
        atom->updatePropertyCache();
    }
}

void removeWaters(RDKit::RWMol& mol)
{
    mol.beginBatchEdit();
    auto is_water = [](const RDKit::Atom* atom) {
        const auto res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());

        // TODO: This seems like it shouldn't be the job of this function; there
        // should be some sort of separate preprocessing step that removes
        // waters and other unwanted residues
        if (!res_info) {
            return false;
        }
        // strip whitespace from residue name
        auto res_name = res_info->getResidueName();
        res_name.erase(std::remove(res_name.begin(), res_name.end(), ' '),
                       res_name.end());
        if (res_info && res_name == "HOH") {
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

void findChainsAndResidues(const RDKit::ROMol& mol,
                           ChainsAndResidues& chains_and_residues)
{
    // Find all chains and residues in the molecule
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        const auto atom = mol.getAtomWithIdx(i);
        const auto res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());

        if (!res_info) {
            if (atom->getAtomicNum() == 1) {
                // This is a hydrogen that is part of a chain, so add it to the
                // chain
                continue;
            }
            throw std::runtime_error(
                "Atom " + std::to_string(atom->getIdx()) + " does not have residue info");
        }
        const auto chain_id = res_info->getChainId();
        const auto res_num = res_info->getResidueNumber();
        const auto ins_code = res_info->getInsertionCode();
        chains_and_residues[chain_id][std::make_pair(res_num, ins_code)]
            .push_back(i);
    }
}

static void copy_selected_atoms_and_bonds(::RDKit::RWMol& extracted_mol,
                                          const RDKit::ROMol& reference_mol,
                                          const std::vector<bool>& selected_atoms,
                                          const std::vector<bool>& selected_bonds)
{
    std::unordered_map<unsigned int, unsigned int> atom_mapping;
    std::unordered_map<unsigned int, unsigned int> bond_mapping;
    for (const auto& ref_atom : reference_mol.atoms()) {
        if (!selected_atoms[ref_atom->getIdx()]) {
            continue;
        }

        std::unique_ptr<::RDKit::Atom> extracted_atom{ref_atom->copy()};

        static constexpr bool updateLabel = true;
        static constexpr bool takeOwnership = true;
        atom_mapping[ref_atom->getIdx()] = extracted_mol.addAtom(
            extracted_atom.release(), updateLabel, takeOwnership);
    }

    for (const auto& ref_bond : reference_mol.bonds()) {
        if (!selected_bonds[ref_bond->getIdx()]) {
            continue;
        }

        std::unique_ptr<::RDKit::Bond> extracted_bond{ref_bond->copy()};
        extracted_bond->setBeginAtomIdx(
            atom_mapping[ref_bond->getBeginAtomIdx()]);
        extracted_bond->setEndAtomIdx(atom_mapping[ref_bond->getEndAtomIdx()]);

        static constexpr bool takeOwnership = true;
        auto num_bonds =
            extracted_mol.addBond(extracted_bond.release(), takeOwnership);
        bond_mapping[ref_bond->getIdx()] = num_bonds - 1;
    }
}

std::unique_ptr<RDKit::RWMol>
extractMolFragment(const RDKit::ROMol& mol,
                   const std::vector<unsigned int>& atom_ids)
{

    const auto num_atoms = mol.getNumAtoms();
    std::vector<bool> selected_atoms(num_atoms);
    std::vector<bool> selected_bonds(mol.getNumBonds());
    for (const auto& atom_idx : atom_ids) {
        if (atom_idx < num_atoms) {
            selected_atoms[atom_idx] = true;
        }
    }
    for (const auto& bond : mol.bonds()) {
        if (selected_atoms[bond->getBeginAtomIdx()] &&
            selected_atoms[bond->getEndAtomIdx()]) {
            selected_bonds[bond->getIdx()] = true;
        }
    }

    auto extracted_mol = std::make_unique<::RDKit::RWMol>();
    copy_selected_atoms_and_bonds(*extracted_mol, mol, selected_atoms, selected_bonds);

    // NOTE: Bookmarks are currently not copied
    return extracted_mol;
}

std::string getMonomerSmiles(RDKit::ROMol& mol,
                             const std::vector<unsigned int>& atom_idxs,
                             const std::string& chain_id,
                             const std::pair<int, std::string>& current_key,
                             int res_num, bool end_of_chain)
{
    // Determine the atoms in current_res that connect to adjacent residues
    std::vector<std::pair<int, int>> attch_idxs; // adjacent res num, atom_idx
    for (auto idx : atom_idxs) {
        const auto at = mol.getAtomWithIdx(idx);
        for (const auto neigh : mol.atomNeighbors(at)) {
            if (neigh->getAtomicNum() == 1) {
                // skip Hs
                continue;
            }
            const auto res_info =
                dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                    neigh->getMonomerInfo());
            const auto key = std::make_pair(res_info->getResidueNumber(),
                                            res_info->getInsertionCode());
            if (key != current_key || res_info->getChainId() != chain_id) {
                // neighbor is in different residue, this will be an attachment
                // point
                attch_idxs.push_back(
                    {res_info->getResidueNumber(), at->getIdx()});
            }
        }
    }

    // This part is difficult; we need to figure out the order of the attachment
    // points on this monomer, Ideally we'd find the backbone and ensure to
    // label them directionally correctly using R2-R1
    std::ranges::sort(
        attch_idxs.begin(), attch_idxs.end(),
        [&res_num](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            // we want the most 'adjacent' residue first in terms of residue
            // ordering; this gets more complicated if we have different chains
            if (std::abs(a.first - res_num) == std::abs(b.first - res_num)) {
                return a.first < b.first;
            }
            return std::abs(a.first - res_num) < std::abs(b.first - res_num);
        });
    auto mol_fragment = extractMolFragment(mol, atom_idxs);

    // Add dummy atoms with an attachment point #s
    // For now, all monomers have at least attachment points 1 and 2 so
    // that backbone bonds can be formed (except the beginning monomer)
    int current_attchpt = 1;

    // If this is the beginning of the chain or a branch, start with attachment
    // point 2
    if (!end_of_chain && res_num == 1) {
        current_attchpt = 2;
    }

    static constexpr bool update_label = true;
    static constexpr bool take_ownership = true;
    for (const auto& [_, ref_idx] : attch_idxs) {
        for (auto at : mol_fragment->atoms()) {
            int orig_idx;
            if (at->getPropIfPresent(REFERENCE_IDX, orig_idx) &&
                orig_idx == ref_idx) {
                auto new_at = new RDKit::Atom(0);
                new_at->setProp(RDKit::common_properties::molAtomMapNumber,
                                current_attchpt);
                auto new_at_idx =
                    mol_fragment->addAtom(new_at, update_label, take_ownership);
                mol_fragment->addBond(new_at_idx, at->getIdx(),
                                      RDKit::Bond::SINGLE);
                mol.getAtomWithIdx(orig_idx)->setProp<unsigned int>(
                    MONOMER_MAP_NUM, current_attchpt);
                ++current_attchpt;
            }
        }
    }
    // There should always be enough attachment points so that
    // backbone connections can be made (R1 and R2)
    // TODO: Should this indicate a new chain?
    while (current_attchpt <= MIN_ATTCHPTS) {
        if (end_of_chain && current_attchpt > 1) {
            break;
        }
        auto new_at = new RDKit::Atom(0);
        new_at->setProp(RDKit::common_properties::molAtomMapNumber,
                        current_attchpt);
        mol_fragment->addAtom(new_at, update_label, take_ownership);
        ++current_attchpt;
    }

    // removing hydrogens to keep HELM string readable
    MolOps::RemoveHsParameters params;
    MolOps::removeHs(*mol_fragment, params, false);
    return RDKit::MolToSmiles(*mol_fragment);
}

bool sameMonomer(RDKit::RWMol& atomistic_mol,
                 const std::vector<unsigned int>& atom_idxs,
                 const std::string& db_smiles)
{
    // Occasionally SMILES that cannot be kekulized are extracted from atomistic
    // mol, so skip sanitization
    constexpr int debug = 0;
    constexpr bool sanitize = false;

    // Monomer from original atomistic mol and monomer defined by database
    auto monomer_frag = extractMolFragment(atomistic_mol, atom_idxs);
    std::unique_ptr<RDKit::RWMol> db_mol(
        RDKit::SmilesToMol(db_smiles, debug, sanitize));

    // Remove stereochemistry, atom map numbers, and neutralize the atoms
    // leaving groups are removed since they won't always be included in the
    // residue extracted from the atomistic molecule
    auto clean_mol = [](RDKit::RWMol& mol) {
        RDKit::MolOps::removeStereochemistry(mol);
        mol.beginBatchEdit();
        for (auto at : mol.atoms()) {
            unsigned int map_no;
            if (at->getPropIfPresent(RDKit::common_properties::molAtomMapNumber,
                                     map_no)) {
                // Label parent with this map number so we know which map number
                // the linkage uses
                for (const auto& nbr : mol.atomNeighbors(at)) {
                    nbr->setProp(RDKit::common_properties::molAtomMapNumber,
                                 map_no);
                }

                mol.removeAtom(at);
            }
        }
        mol.commitBatchEdit();

        MolOps::RemoveHsParameters params;
        MolOps::removeHs(mol, params, false);
        // set aromaticity
        RDKit::MolOps::setAromaticity(mol);
    };
    clean_mol(*monomer_frag);
    clean_mol(*db_mol);

    // The DB monomer has had the leaving groups removed, while the residue
    // extracted from the atomistic mol may still have them present if it is at
    // the beginning or end of a chain. As a result, we need to allow for the DB
    // monomer to have one less atom than the residue extracted from the
    // atomistic mol
    auto match = RDKit::SubstructMatch(*monomer_frag,
                                       *db_mol); // (queryAtomIdx, molAtomIdx)

    if (match.size()) {

        // Any unmapped atom in monomer_frag must map to a leaving group in the
        // db_mol create vector of atom indices in the monomer_frag that are not
        // in the match
        for (auto at : monomer_frag->atoms()) {
            int target_idx = at->getIdx();
            auto it = std::find_if(match[0].begin(), match[0].end(),
                                   [target_idx](const std::pair<int, int>& p) {
                                       return p.second == target_idx;
                                   });
            if (it == match[0].end()) {
                // This atom should be terminal, otherwise the match is
                // definitely wrong
                if (at->getDegree() > 1) {
                    return false;
                }

                // Unmatched atoms should be leaving groups, we can check this
                // by matching
                RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
                boost::tie(nbrIdx, endNbrs) =
                    monomer_frag->getAtomNeighbors(at);

                it = std::find_if(match[0].begin(), match[0].end(),
                                  [nbrIdx](const std::pair<int, int>& p) {
                                      return p.second ==
                                             static_cast<int>(*nbrIdx);
                                  });
                if (it == match[0].end()) {
                    return false;
                }

                auto db_atom = db_mol->getAtomWithIdx(it->first);
                if (!db_atom->hasProp(
                        RDKit::common_properties::molAtomMapNumber)) {
                    return false;
                }
            }
        }

        // we are mapping atomistic atom indices to attachment point numbers
        // (molAtomMapNumber)
        for (auto [db_idx, at_idx] : match[0]) {
            auto db_at = db_mol->getAtomWithIdx(db_idx);
            auto at = monomer_frag->getAtomWithIdx(at_idx);
            unsigned int map_no;
            if (db_at->getPropIfPresent(
                    RDKit::common_properties::molAtomMapNumber, map_no)) {
                auto ref_idx = at->getProp<unsigned int>(REFERENCE_IDX);
                atomistic_mol.getAtomWithIdx(ref_idx)->setProp<unsigned int>(
                    MONOMER_MAP_NUM, map_no);
            }
        }
        return true;
    }

    return false;
}

int getAttchpt(const RDKit::Atom& monomer, const RDKit::Bond& bond,
               const RDKit::RWMol& atomistic_mol)
{
    // This should be the index of the atom in this monomer
    auto target_atom_idx = bond.getBeginAtomIdx();
    if (atomistic_mol.getAtomWithIdx(target_atom_idx)
            ->getProp<unsigned int>(MONOMER_IDX) != monomer.getIdx()) {
        target_atom_idx = bond.getEndAtomIdx();
    }
    unsigned int attchpt;
    if (atomistic_mol.getAtomWithIdx(target_atom_idx)
            ->getPropIfPresent(MONOMER_MAP_NUM, attchpt)) {
        return attchpt;
    } else {
        return -1;
    }
}

void detectLinkages(MonomerMol& monomer_mol,
                    const RDKit::RWMol& atomistic_mol)
{
    // Find all linkages between monomers (used when PDB residue info is
    // available)
    for (const auto& bond : atomistic_mol.bonds()) {
        const auto begin_atom = bond->getBeginAtom();
        const auto end_atom = bond->getEndAtom();

        unsigned int begin_monomer_idx =
            std::numeric_limits<unsigned int>::max();
        unsigned int end_monomer_idx = std::numeric_limits<unsigned int>::max();
        if (!begin_atom->getPropIfPresent(MONOMER_IDX, begin_monomer_idx) ||
            !end_atom->getPropIfPresent(MONOMER_IDX, end_monomer_idx)) {
            // One of these is a hydrogen or unmapped for some reason
            continue;
        }

        // Check if this is already present as a backbone linkage
        if (begin_monomer_idx == end_monomer_idx ||
            monomer_mol.getBondBetweenAtoms(begin_monomer_idx,
                                            end_monomer_idx) != nullptr) {
            continue;
        }

        const auto begin_monomer =
            monomer_mol.getAtomWithIdx(begin_monomer_idx);
        const auto end_monomer = monomer_mol.getAtomWithIdx(end_monomer_idx);

        auto begin_attchpt = getAttchpt(*begin_monomer, *bond, atomistic_mol);
        auto end_attchpt = getAttchpt(*end_monomer, *bond, atomistic_mol);

        if (begin_attchpt == -1 || end_attchpt == -1) {
            // This happens when the input atomistic mol has bonds between
            // residues on atoms that are not marked as attachment points in the
            // monomer DB. This is important to capture so that we can determine
            // where we may need to add attachment points.
            auto begin_res_info =
                dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                    begin_atom->getMonomerInfo());
            auto end_res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                end_atom->getMonomerInfo());
            if (begin_res_info && end_res_info) {
                std::cerr << "Bond between residue " << begin_res_info->getResidueName()
                          << begin_res_info->getResidueNumber() << " in chain "
                          << begin_res_info->getChainId() << " and residue "
                          << end_res_info->getResidueName() << end_res_info->getResidueNumber()
                          << " in chain " << end_res_info->getChainId()
                          << " do not correspond to attachment points in their residues.\n";
            } else {
                std::cerr << "Bond between atoms " << begin_atom->getIdx()
                          << " and " << end_atom->getIdx() << " in monomers "
                          << begin_monomer->getProp<std::string>(ATOM_LABEL)
                          << " and " << end_monomer->getProp<std::string>(ATOM_LABEL)
                          << " do not correspond to attachment points in their residues.\n";
            }
            continue;
        }

        // Backbone connections (R2-R1) should be added in that order when
        // possible.
        if (begin_attchpt == 2 && end_attchpt == 1) {
            monomer_mol.addConnection(begin_monomer_idx, end_monomer_idx,
                          "R2-R1");
        } else if (begin_attchpt == 1 && end_attchpt == 2) {
            monomer_mol.addConnection(end_monomer_idx, begin_monomer_idx,
                          "R2-R1");
        } else if (begin_monomer_idx < end_monomer_idx) {
            monomer_mol.addConnection(begin_monomer_idx, end_monomer_idx,
                          "R" + std::to_string(begin_attchpt) + "-R" + std::to_string(end_attchpt));
        } else {
            monomer_mol.addConnection(end_monomer_idx, begin_monomer_idx,
                          "R" + std::to_string(end_attchpt) + "-R" + std::to_string(begin_attchpt));
        }
    }
}

MonomerLibrary::monomer_info_t
getMonomerInfoFromAtom(const MonomerLibrary& db, const RDKit::Atom* atom)
{
    // This comes from something like a PDB or MAE file
    auto res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
        atom->getMonomerInfo());
    auto res_name = res_info->getResidueName();
    res_name.erase(
        std::remove_if(res_name.begin(), res_name.end(), ::isspace),
        res_name.end());
    return db.getMonomerInfo(res_name);
}

std::unique_ptr<MonomerMol>
pdbInfoAtomisticToMM(const RDKit::ROMol& input_mol)
{
    // A couple preprocessing steps; remove waters and neutralize atoms
    RDKit::RWMol mol(input_mol);
    removeWaters(mol);
    neutralizeAtoms(mol);

    // Set reference index for SMILES fragments
    for (auto at : mol.atoms()) {
        at->setProp(REFERENCE_IDX, at->getIdx());
    }

    // Map chain_id -> {residue mols}
    ChainsAndResidues chains_and_residues;
    findChainsAndResidues(mol, chains_and_residues);

    std::map<std::string, unsigned int> chain_counts = {{"PEPTIDE", 0},
                                                        {"RNA", 0},
                                                        {"DNA", 0},
                                                        {"CHEM", 0}};
    auto monomer_mol = std::make_unique<MonomerMol>();

    // Access the monomer library via the MonomerMol instance
    const auto& db = monomer_mol->getMonomerLibrary();

    for (const auto& [chain_id, residues] : chains_and_residues) {
        // Use first residue to determine chain type. We assume that PDB data
        // is correct and there aren't multiple chain types in a single chain.
        // Default chain type is PEPTIDE if not specified.
        auto monomer_info = getMonomerInfoFromAtom(
            db, mol.getAtomWithIdx(residues.begin()->second[0]));
        auto monomer_class =
            monomer_info ? std::get<2>(*monomer_info) : std::string("PEPTIDE");
        std::string polymer_chain_id = monomer_class +
                                     std::to_string(++chain_counts[monomer_class]);
        // Assuming residues are ordered correctly
        unsigned int res_num = 1;
        for (const auto& [key, atom_idxs] : residues) {
            monomer_info = getMonomerInfoFromAtom(db, mol.getAtomWithIdx(atom_idxs[0]));
            bool end_of_chain = res_num == residues.size();
            unsigned int this_monomer;
            if (monomer_info &&
                sameMonomer(mol, atom_idxs, std::get<1>(*monomer_info))) {
                // Standard residue in monomer DB, Verify that the fragment
                // labeled as the residue matches what is in the monomer
                // database
                this_monomer = monomer_mol->addMonomer(std::get<0>(*monomer_info),
                                          res_num, monomer_class, polymer_chain_id);
            } else {
                auto smiles = getMonomerSmiles(mol, atom_idxs, chain_id, key,
                                               res_num, end_of_chain);
                this_monomer = monomer_mol->addMonomer(smiles, res_num,
                                          monomer_class, polymer_chain_id, MonomerType::SMILES);
            }

            // Track which atoms are in which monomer
            for (auto idx : atom_idxs) {
                mol.getAtomWithIdx(idx)->setProp<unsigned int>(MONOMER_IDX,
                                                               this_monomer);
            }
            ++res_num;
        }
    }
    detectLinkages(*monomer_mol, mol);
    return monomer_mol;
}

bool hasPdbResidueInfo(const RDKit::ROMol& mol)
{
    return std::any_of(mol.atoms().begin(), mol.atoms().end(),
                       [](const auto& atom) {
                           return atom->getMonomerInfo() &&
                                  atom->getMonomerInfo()->getMonomerType() ==
                                      RDKit::AtomMonomerInfo::PDBRESIDUE;
                       });
}
} // unnamed namespace

std::unique_ptr<MonomerMol> toMonomeric(const RDKit::ROMol& atomistic_mol)
{
    if (!hasPdbResidueInfo(atomistic_mol)) {
        // If there is no residue information, we cannot convert to monomeric
        // form, so throw an error
        throw std::runtime_error(
            "No residue information found in molecule, cannot convert to "
            "monomeric form");
    }
    auto monomer_mol = pdbInfoAtomisticToMM(atomistic_mol);
    monomer_mol->assignChains();
    return monomer_mol;
}

} // namespace RDKit